package tensornetwork

import (
	"errors"
	"fmt"
	"sort"
	"strings"
)

// Outer returns the outer (tensor) product of a and b: a tensor whose shape is
// the concatenation of the shapes of a and b, with
// out[i…, j…] = a[i…] * b[j…].
func Outer(a, b *Tensor) *Tensor {
	sh := append(append([]int(nil), a.shape...), b.shape...)
	out := New(sh...)
	for i, av := range a.data {
		base := i * len(b.data)
		for j, bv := range b.data {
			out.data[base+j] = av * bv
		}
	}
	return out
}

// Kronecker returns the Kronecker product of two tensors of equal rank. Along
// each axis the output dimension is the product of the input dimensions, with
// out[i₀*b₀+j₀, …] = a[i…] * b[j…]. It returns an error if the ranks differ.
func Kronecker(a, b *Tensor) (*Tensor, error) {
	if len(a.shape) != len(b.shape) {
		return nil, errors.New("tensornetwork: Kronecker requires equal rank")
	}
	n := len(a.shape)
	sh := make([]int, n)
	for i := 0; i < n; i++ {
		sh[i] = a.shape[i] * b.shape[i]
	}
	out := New(sh...)
	ai := make([]int, n)
	bi := make([]int, n)
	oi := make([]int, n)
	for af := 0; af < len(a.data); af++ {
		rem := af
		for k := 0; k < n; k++ {
			ai[k] = (rem / a.stride[k]) % a.shape[k]
		}
		for bf := 0; bf < len(b.data); bf++ {
			rem := bf
			for k := 0; k < n; k++ {
				bi[k] = (rem / b.stride[k]) % b.shape[k]
			}
			for k := 0; k < n; k++ {
				oi[k] = ai[k]*b.shape[k] + bi[k]
			}
			off, _ := out.flatIndex(oi)
			out.data[off] = a.data[af] * b.data[bf]
		}
	}
	return out, nil
}

// Unfold returns the mode-n matricization of t: a matrix whose rows are indexed
// by axis mode and whose columns are indexed by the remaining axes taken in
// increasing order, with the first remaining axis varying slowest. It returns an
// error if mode is out of range.
func (t *Tensor) Unfold(mode int) (*Matrix, error) {
	n := len(t.shape)
	if mode < 0 || mode >= n {
		return nil, errors.New("tensornetwork: mode out of range")
	}
	rows := t.shape[mode]
	cols := len(t.data) / rows
	m := NewMatrix(rows, cols)
	var rest []int
	for a := 0; a < n; a++ {
		if a != mode {
			rest = append(rest, a)
		}
	}
	idx := make([]int, n)
	for flat := 0; flat < len(t.data); flat++ {
		rem := flat
		for a := 0; a < n; a++ {
			idx[a] = (rem / t.stride[a]) % t.shape[a]
		}
		col := 0
		for _, a := range rest {
			col = col*t.shape[a] + idx[a]
		}
		m.data[idx[mode]*cols+col] = t.data[flat]
	}
	return m, nil
}

// Fold reconstructs a tensor of the given shape from its mode-n matricization,
// inverting [Tensor.Unfold]. It returns an error if the matrix dimensions are
// inconsistent with the shape.
func Fold(m *Matrix, mode int, shape []int) (*Tensor, error) {
	n := len(shape)
	if mode < 0 || mode >= n {
		return nil, errors.New("tensornetwork: mode out of range")
	}
	if m.rows != shape[mode] {
		return nil, fmt.Errorf("tensornetwork: fold row count %d != shape[%d]=%d", m.rows, mode, shape[mode])
	}
	if m.rows*m.cols != sizeOf(shape) {
		return nil, errors.New("tensornetwork: fold size mismatch")
	}
	out := New(shape...)
	var rest []int
	for a := 0; a < n; a++ {
		if a != mode {
			rest = append(rest, a)
		}
	}
	idx := make([]int, n)
	for flat := 0; flat < len(out.data); flat++ {
		rem := flat
		for a := 0; a < n; a++ {
			idx[a] = (rem / out.stride[a]) % shape[a]
		}
		col := 0
		for _, a := range rest {
			col = col*shape[a] + idx[a]
		}
		out.data[flat] = m.data[idx[mode]*m.cols+col]
	}
	return out, nil
}

// ModeProduct returns the n-mode product t ×ₙ m, contracting axis mode of t with
// the columns of m. If t has shape (…, Iₙ, …) and m is J x Iₙ the result has
// axis mode replaced by J. It returns an error on a dimension mismatch.
func ModeProduct(t *Tensor, m *Matrix, mode int) (*Tensor, error) {
	n := len(t.shape)
	if mode < 0 || mode >= n {
		return nil, errors.New("tensornetwork: mode out of range")
	}
	if m.cols != t.shape[mode] {
		return nil, fmt.Errorf("tensornetwork: mode product mismatch: matrix cols %d != dim %d", m.cols, t.shape[mode])
	}
	unf, err := t.Unfold(mode)
	if err != nil {
		return nil, err
	}
	prod, err := m.Mul(unf)
	if err != nil {
		return nil, err
	}
	newShape := append([]int(nil), t.shape...)
	newShape[mode] = m.rows
	return Fold(prod, mode, newShape)
}

// MultiModeProduct applies a sequence of n-mode products, contracting each given
// mode of t with the corresponding matrix. The modes and matrices slices must
// have equal length. It returns an error on any dimension mismatch.
func MultiModeProduct(t *Tensor, mats []*Matrix, modes []int) (*Tensor, error) {
	if len(mats) != len(modes) {
		return nil, errors.New("tensornetwork: mats and modes length mismatch")
	}
	cur := t
	var err error
	for i, mode := range modes {
		cur, err = ModeProduct(cur, mats[i], mode)
		if err != nil {
			return nil, err
		}
	}
	return cur, nil
}

// TensorDot contracts a and b over the axis lists axesA and axesB, summing over
// the paired axes. The result axes are the free axes of a followed by the free
// axes of b, each in their original order. It returns an error if the axis lists
// disagree in length or the paired dimensions do not match.
func TensorDot(a, b *Tensor, axesA, axesB []int) (*Tensor, error) {
	if len(axesA) != len(axesB) {
		return nil, errors.New("tensornetwork: axesA and axesB length mismatch")
	}
	na, nb := len(a.shape), len(b.shape)
	usedA := make([]bool, na)
	usedB := make([]bool, nb)
	for k := range axesA {
		ax, bx := axesA[k], axesB[k]
		if ax < 0 || ax >= na || bx < 0 || bx >= nb {
			return nil, errors.New("tensornetwork: contraction axis out of range")
		}
		if a.shape[ax] != b.shape[bx] {
			return nil, fmt.Errorf("tensornetwork: contracted dim mismatch %d != %d", a.shape[ax], b.shape[bx])
		}
		if usedA[ax] || usedB[bx] {
			return nil, errors.New("tensornetwork: repeated contraction axis")
		}
		usedA[ax] = true
		usedB[bx] = true
	}
	var freeA, freeB []int
	for i := 0; i < na; i++ {
		if !usedA[i] {
			freeA = append(freeA, i)
		}
	}
	for i := 0; i < nb; i++ {
		if !usedB[i] {
			freeB = append(freeB, i)
		}
	}
	// Permute a to (freeA..., contractedA...) and b to (contractedB..., freeB...).
	permA := append(append([]int(nil), freeA...), axesA...)
	permB := append(append([]int(nil), axesB...), freeB...)
	pa, err := a.Permute(permA...)
	if err != nil {
		return nil, err
	}
	pb, err := b.Permute(permB...)
	if err != nil {
		return nil, err
	}
	freeASize := 1
	for _, i := range freeA {
		freeASize *= a.shape[i]
	}
	freeBSize := 1
	for _, i := range freeB {
		freeBSize *= b.shape[i]
	}
	contractSize := 1
	for _, i := range axesA {
		contractSize *= a.shape[i]
	}
	ma, err := pa.Reshape(freeASize, contractSize)
	if err != nil {
		return nil, err
	}
	mb, err := pb.Reshape(contractSize, freeBSize)
	if err != nil {
		return nil, err
	}
	mam, _ := ma.ToMatrix()
	mbm, _ := mb.ToMatrix()
	prod, err := mam.Mul(mbm)
	if err != nil {
		return nil, err
	}
	var outShape []int
	for _, i := range freeA {
		outShape = append(outShape, a.shape[i])
	}
	for _, i := range freeB {
		outShape = append(outShape, b.shape[i])
	}
	if len(outShape) == 0 {
		return Scalar(prod.data[0]), nil
	}
	return NewWithData(prod.data, outShape...)
}

// Contract contracts a single axis of a with a single axis of b, a convenience
// wrapper over [TensorDot].
func Contract(a, b *Tensor, axisA, axisB int) (*Tensor, error) {
	return TensorDot(a, b, []int{axisA}, []int{axisB})
}

// TraceAxes contracts axes a and b of a single tensor (a partial trace),
// summing over the diagonal where the two indices are equal and removing both
// axes. It returns an error if the two axes have different sizes.
func (t *Tensor) TraceAxes(a, b int) (*Tensor, error) {
	n := len(t.shape)
	if a < 0 || a >= n || b < 0 || b >= n || a == b {
		return nil, errors.New("tensornetwork: invalid trace axes")
	}
	if t.shape[a] != t.shape[b] {
		return nil, errors.New("tensornetwork: trace axes have different sizes")
	}
	var sh []int
	for i, d := range t.shape {
		if i != a && i != b {
			sh = append(sh, d)
		}
	}
	out := New(sh...)
	idx := make([]int, n)
	for flat := 0; flat < len(t.data); flat++ {
		rem := flat
		for k := 0; k < n; k++ {
			idx[k] = (rem / t.stride[k]) % t.shape[k]
		}
		if idx[a] != idx[b] {
			continue
		}
		var outIdx []int
		for k := 0; k < n; k++ {
			if k != a && k != b {
				outIdx = append(outIdx, idx[k])
			}
		}
		off := 0
		if len(outIdx) > 0 {
			off, _ = out.flatIndex(outIdx)
		}
		out.data[off] += t.data[flat]
	}
	return out, nil
}

// Einsum evaluates an Einstein-summation expression over the given operands. The
// spec has the form "ij,jk->ik": comma-separated input index strings and an
// optional output index string after "->". Indices that appear in inputs but not
// in the output are summed; a repeated index within one operand contracts its
// diagonal. If "->" is omitted the output indices are those appearing exactly
// once, in alphabetical order. It returns an error on any malformed spec or
// dimension mismatch.
func Einsum(spec string, operands ...*Tensor) (*Tensor, error) {
	spec = strings.ReplaceAll(spec, " ", "")
	var lhs, rhs string
	explicit := false
	if i := strings.Index(spec, "->"); i >= 0 {
		lhs = spec[:i]
		rhs = spec[i+2:]
		explicit = true
	} else {
		lhs = spec
	}
	inSpecs := strings.Split(lhs, ",")
	if len(inSpecs) != len(operands) {
		return nil, fmt.Errorf("tensornetwork: einsum has %d input specs but %d operands", len(inSpecs), len(operands))
	}
	// Determine dimension of each label and validate.
	dims := make(map[rune]int)
	counts := make(map[rune]int)
	for oi, s := range inSpecs {
		if len(s) != len(operands[oi].shape) {
			return nil, fmt.Errorf("tensornetwork: spec %q rank != operand rank %d", s, len(operands[oi].shape))
		}
		for ai, r := range s {
			d := operands[oi].shape[ai]
			if prev, ok := dims[r]; ok {
				if prev != d {
					return nil, fmt.Errorf("tensornetwork: einsum label %c size mismatch %d != %d", r, prev, d)
				}
			} else {
				dims[r] = d
			}
			counts[r]++
		}
	}
	if !explicit {
		var once []rune
		for r, c := range counts {
			if c == 1 {
				once = append(once, r)
			}
		}
		sort.Slice(once, func(i, j int) bool { return once[i] < once[j] })
		rhs = string(once)
	}
	outLabels := []rune(rhs)
	seenOut := make(map[rune]bool)
	for _, r := range outLabels {
		if _, ok := dims[r]; !ok {
			return nil, fmt.Errorf("tensornetwork: output label %c not in inputs", r)
		}
		if seenOut[r] {
			return nil, fmt.Errorf("tensornetwork: repeated output label %c", r)
		}
		seenOut[r] = true
	}
	// Summed labels: all labels not in output.
	var sumLabels []rune
	for r := range dims {
		if !seenOut[r] {
			sumLabels = append(sumLabels, r)
		}
	}
	sort.Slice(sumLabels, func(i, j int) bool { return sumLabels[i] < sumLabels[j] })
	allLabels := append(append([]rune(nil), outLabels...), sumLabels...)
	labelPos := make(map[rune]int)
	for i, r := range allLabels {
		labelPos[r] = i
	}
	// Output tensor.
	var outShape []int
	for _, r := range outLabels {
		outShape = append(outShape, dims[r])
	}
	out := New(outShape...)
	outSize := len(out.data)
	if outSize == 0 {
		outSize = 1
	}
	// Iterate over the full cartesian product of all labels.
	counters := make([]int, len(allLabels))
	total := 1
	for _, r := range allLabels {
		total *= dims[r]
	}
	inSpecRunes := make([][]rune, len(inSpecs))
	for i, s := range inSpecs {
		inSpecRunes[i] = []rune(s)
	}
	for it := 0; it < total; it++ {
		// Compute product of operand entries for the current assignment.
		prod := 1.0
		valid := true
		for oi := range operands {
			op := operands[oi]
			idx := make([]int, len(op.shape))
			for ai, r := range inSpecRunes[oi] {
				idx[ai] = counters[labelPos[r]]
			}
			off, err := op.flatIndex(idx)
			if err != nil {
				valid = false
				break
			}
			prod *= op.data[off]
		}
		if valid && prod != 0 {
			outIdx := make([]int, len(outLabels))
			for i, r := range outLabels {
				outIdx[i] = counters[labelPos[r]]
			}
			off := 0
			if len(outIdx) > 0 {
				off, _ = out.flatIndex(outIdx)
			}
			out.data[off] += prod
		}
		// Increment mixed-radix counters.
		for k := len(allLabels) - 1; k >= 0; k-- {
			counters[k]++
			if counters[k] < dims[allLabels[k]] {
				break
			}
			counters[k] = 0
		}
	}
	if len(outShape) == 0 {
		return Scalar(out.data[0]), nil
	}
	return out, nil
}

// KhatriRaoTensors returns the column-wise Khatri-Rao product of a list of
// rank-2 tensors interpreted as matrices, a convenience wrapper over
// [KhatriRao]. It returns an error if any tensor is not rank 2.
func KhatriRaoTensors(tensors ...*Tensor) (*Tensor, error) {
	mats := make([]*Matrix, len(tensors))
	for i, t := range tensors {
		m, err := t.ToMatrix()
		if err != nil {
			return nil, err
		}
		mats[i] = m
	}
	kr, err := KhatriRao(mats...)
	if err != nil {
		return nil, err
	}
	return FromMatrix(kr), nil
}
