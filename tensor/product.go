package tensor

import "fmt"

// Outer returns the tensor (outer) product a⊗b. The result has rank
// a.Rank()+b.Rank() and shape a.Shape() followed by b.Shape(); its element at
// the concatenated multi-index (i,j) equals a[i]*b[j]. The outer product of two
// vectors is the rank-2 matrix a_i b_j.
func Outer(a, b *Tensor) *Tensor {
	outShape := make([]int, 0, a.Rank()+b.Rank())
	outShape = append(outShape, a.shape...)
	outShape = append(outShape, b.shape...)
	var out *Tensor
	if len(outShape) == 0 {
		out = FromScalar(a.data[0] * b.data[0])
		return out
	}
	out = New(outShape...)
	nb := b.Size()
	for i, av := range a.data {
		for j, bv := range b.data {
			out.data[i*nb+j] = av * bv
		}
	}
	return out
}

// TensorProduct is a synonym for [Outer], using the name common in
// multilinear-algebra texts for the product a⊗b.
func TensorProduct(a, b *Tensor) *Tensor { return Outer(a, b) }

// Kronecker returns the Kronecker product of two tensors of equal rank. The
// result has, on each axis k, dimension a.shape[k]*b.shape[k], and its element
// at index (i_k) with i_k = ia_k*b.shape[k]+ib_k equals a[ia]*b[ib]. For two
// matrices this is the usual block Kronecker product. It returns [ErrRank] if
// the ranks differ.
func Kronecker(a, b *Tensor) (*Tensor, error) {
	if a.Rank() != b.Rank() {
		return nil, fmt.Errorf("%w: Kronecker needs equal ranks, got %d and %d", ErrRank, a.Rank(), b.Rank())
	}
	r := a.Rank()
	if r == 0 {
		return FromScalar(a.data[0] * b.data[0]), nil
	}
	outShape := make([]int, r)
	for k := 0; k < r; k++ {
		outShape[k] = a.shape[k] * b.shape[k]
	}
	out := New(outShape...)
	idx := make([]int, r)
	for fa := 0; fa < a.Size(); fa++ {
		ia := tensorUnravel(fa, a.shape)
		for fb := 0; fb < b.Size(); fb++ {
			ib := tensorUnravel(fb, b.shape)
			for k := 0; k < r; k++ {
				idx[k] = ia[k]*b.shape[k] + ib[k]
			}
			out.Set(a.data[fa]*b.data[fb], idx...)
		}
	}
	return out, nil
}

// TensorDot contracts a and b over the paired axes in axesA and axesB. For each
// i the axis axesA[i] of a is summed against axis axesB[i] of b, which must have
// matching lengths. The result's axes are the free axes of a (in order) followed
// by the free axes of b. With empty axis lists TensorDot reduces to [Outer]; the
// classic matrix product is TensorDot(a,b,[]int{1},[]int{0}). It returns
// [ErrAxis] for malformed axis lists and [ErrShape] when paired axes differ in
// length.
func TensorDot(a, b *Tensor, axesA, axesB []int) (*Tensor, error) {
	if len(axesA) != len(axesB) {
		return nil, fmt.Errorf("%w: axesA and axesB must have equal length", ErrAxis)
	}
	if !tensorDistinctInRange(axesA, a.Rank()) || !tensorDistinctInRange(axesB, b.Rank()) {
		return nil, fmt.Errorf("%w: contraction axes out of range or repeated", ErrAxis)
	}
	for i := range axesA {
		if a.shape[axesA[i]] != b.shape[axesB[i]] {
			return nil, fmt.Errorf("%w: contracted axes %d and %d differ in length", ErrShape, axesA[i], axesB[i])
		}
	}
	freeA := tensorComplement(axesA, a.Rank())
	freeB := tensorComplement(axesB, b.Rank())

	outShape := make([]int, 0, len(freeA)+len(freeB))
	for _, ax := range freeA {
		outShape = append(outShape, a.shape[ax])
	}
	for _, ax := range freeB {
		outShape = append(outShape, b.shape[ax])
	}

	// Shape of the summed (contracted) index space.
	conShape := make([]int, len(axesA))
	for i, ax := range axesA {
		conShape[i] = a.shape[ax]
	}
	conSize := tensorProduct(conShape)

	var out *Tensor
	if len(outShape) == 0 {
		out = FromScalar(0)
	} else {
		out = New(outShape...)
	}

	aIdx := make([]int, a.Rank())
	bIdx := make([]int, b.Rank())
	oIdx := make([]int, len(outShape))

	for of := 0; of < out.Size(); of++ {
		var ofree []int
		if len(outShape) == 0 {
			ofree = nil
		} else {
			ofree = tensorUnravel(of, outShape)
		}
		copy(oIdx, ofree)
		// Place free-axis values into aIdx and bIdx.
		for i, ax := range freeA {
			aIdx[ax] = oIdx[i]
		}
		for i, ax := range freeB {
			bIdx[ax] = oIdx[len(freeA)+i]
		}
		sum := 0.0
		for cf := 0; cf < conSize; cf++ {
			cidx := tensorUnravel(cf, conShape)
			for i := range axesA {
				aIdx[axesA[i]] = cidx[i]
				bIdx[axesB[i]] = cidx[i]
			}
			sum += a.At(aIdx...) * b.At(bIdx...)
		}
		out.data[of] = sum
	}
	return out, nil
}

// tensorDistinctInRange reports whether every element of axes is in [0,n) and
// no element repeats.
func tensorDistinctInRange(axes []int, n int) bool {
	seen := make(map[int]bool, len(axes))
	for _, ax := range axes {
		if ax < 0 || ax >= n || seen[ax] {
			return false
		}
		seen[ax] = true
	}
	return true
}

// tensorComplement returns, in increasing order, the axes in [0,n) that are not
// present in axes.
func tensorComplement(axes []int, n int) []int {
	in := make(map[int]bool, len(axes))
	for _, ax := range axes {
		in[ax] = true
	}
	out := make([]int, 0, n-len(axes))
	for i := 0; i < n; i++ {
		if !in[i] {
			out = append(out, i)
		}
	}
	return out
}

// MatMul returns the matrix product of two rank-2 tensors a (m×k) and b (k×n),
// giving an m×n tensor. It returns [ErrRank] if either operand is not a matrix
// and [ErrShape] if their inner dimensions disagree.
func MatMul(a, b *Tensor) (*Tensor, error) {
	if a.Rank() != 2 || b.Rank() != 2 {
		return nil, fmt.Errorf("%w: MatMul requires rank-2 operands", ErrRank)
	}
	if a.shape[1] != b.shape[0] {
		return nil, fmt.Errorf("%w: inner dimensions %d and %d disagree", ErrShape, a.shape[1], b.shape[0])
	}
	return TensorDot(a, b, []int{1}, []int{0})
}

// Dot returns the inner product of two rank-1 tensors of equal length. It
// returns [ErrRank] if either operand is not a vector and [ErrShape] if their
// lengths differ.
func Dot(a, b *Tensor) (float64, error) {
	if a.Rank() != 1 || b.Rank() != 1 {
		return 0, fmt.Errorf("%w: Dot requires rank-1 operands", ErrRank)
	}
	if a.shape[0] != b.shape[0] {
		return 0, fmt.Errorf("%w: vectors of length %d and %d", ErrShape, a.shape[0], b.shape[0])
	}
	s := 0.0
	for i := range a.data {
		s += a.data[i] * b.data[i]
	}
	return s, nil
}

// Contract returns the tensor obtained by summing t over a pair of its own axes
// axisA and axisB, which must be distinct and have equal length. The two axes
// are removed and the rank drops by two; for a matrix this is the ordinary
// trace. It returns [ErrAxis] if the axes coincide or are out of range and
// [ErrShape] if they differ in length. See also [Tensor.Trace].
func (t *Tensor) Contract(axisA, axisB int) (*Tensor, error) {
	r := t.Rank()
	if axisA < 0 || axisA >= r || axisB < 0 || axisB >= r || axisA == axisB {
		return nil, fmt.Errorf("%w: contraction axes %d and %d for rank-%d tensor", ErrAxis, axisA, axisB, r)
	}
	if t.shape[axisA] != t.shape[axisB] {
		return nil, fmt.Errorf("%w: contracted axes of length %d and %d", ErrShape, t.shape[axisA], t.shape[axisB])
	}
	keep := make([]int, 0, r-2)
	for i := 0; i < r; i++ {
		if i != axisA && i != axisB {
			keep = append(keep, i)
		}
	}
	newShape := make([]int, len(keep))
	for i, ax := range keep {
		newShape[i] = t.shape[ax]
	}
	var out *Tensor
	if len(newShape) == 0 {
		out = FromScalar(0)
	} else {
		out = New(newShape...)
	}
	oidx := make([]int, len(keep))
	for f := 0; f < t.Size(); f++ {
		idx := tensorUnravel(f, t.shape)
		if idx[axisA] != idx[axisB] {
			continue
		}
		for i, ax := range keep {
			oidx[i] = idx[ax]
		}
		out.addAt(t.data[f], oidx)
	}
	return out, nil
}

// Trace returns the sum of the diagonal of a rank-2 tensor. It returns [ErrRank]
// if t is not a matrix and [ErrShape] if the matrix is not square. For contraction
// over an arbitrary axis pair of a higher-rank tensor use [Tensor.Contract].
func (t *Tensor) Trace() (float64, error) {
	if t.Rank() != 2 {
		return 0, fmt.Errorf("%w: Trace requires a rank-2 tensor", ErrRank)
	}
	if t.shape[0] != t.shape[1] {
		return 0, fmt.Errorf("%w: Trace requires a square matrix, got %v", ErrShape, t.shape)
	}
	c, err := t.Contract(0, 1)
	if err != nil {
		return 0, err
	}
	return c.data[0], nil
}
