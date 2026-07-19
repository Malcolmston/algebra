package tensornetwork

import (
	"errors"
	"fmt"
	"math"
	"math/rand"
	"strings"
)

// Tensor is a dense, row-major (C-ordered) N-dimensional array of float64
// values. A rank-0 tensor is a scalar, rank-1 a vector and rank-2 a matrix; the
// same code paths handle arbitrary rank.
type Tensor struct {
	shape  []int
	stride []int
	data   []float64
}

// rowMajorStrides returns the row-major strides for the given shape.
func rowMajorStrides(shape []int) []int {
	st := make([]int, len(shape))
	acc := 1
	for i := len(shape) - 1; i >= 0; i-- {
		st[i] = acc
		acc *= shape[i]
	}
	return st
}

// sizeOf returns the product of the dimensions in shape (1 for an empty shape).
func sizeOf(shape []int) int {
	n := 1
	for _, d := range shape {
		n *= d
	}
	return n
}

// New returns a new tensor of the given shape with all entries zero. A call with
// no arguments returns a scalar tensor. It panics if any dimension is negative.
func New(shape ...int) *Tensor {
	for _, d := range shape {
		if d < 0 {
			panic("tensornetwork: negative dimension")
		}
	}
	sh := append([]int(nil), shape...)
	return &Tensor{shape: sh, stride: rowMajorStrides(sh), data: make([]float64, sizeOf(sh))}
}

// Zeros is an alias for [New]: it returns a zero-filled tensor of the given
// shape.
func Zeros(shape ...int) *Tensor { return New(shape...) }

// Ones returns a tensor of the given shape with every entry equal to 1.
func Ones(shape ...int) *Tensor { return Full(1, shape...) }

// Full returns a tensor of the given shape with every entry equal to v.
func Full(v float64, shape ...int) *Tensor {
	t := New(shape...)
	for i := range t.data {
		t.data[i] = v
	}
	return t
}

// NewWithData returns a tensor of the given shape wrapping a copy of data, which
// must be in row-major order. It returns an error if len(data) does not equal
// the product of the shape.
func NewWithData(data []float64, shape ...int) (*Tensor, error) {
	sh := append([]int(nil), shape...)
	if len(data) != sizeOf(sh) {
		return nil, fmt.Errorf("tensornetwork: data length %d != shape product %d", len(data), sizeOf(sh))
	}
	cp := make([]float64, len(data))
	copy(cp, data)
	return &Tensor{shape: sh, stride: rowMajorStrides(sh), data: cp}, nil
}

// Scalar returns a rank-0 tensor holding the value v.
func Scalar(v float64) *Tensor {
	return &Tensor{shape: []int{}, stride: []int{}, data: []float64{v}}
}

// FromVector returns a rank-1 tensor holding a copy of v.
func FromVector(v []float64) *Tensor {
	t, _ := NewWithData(v, len(v))
	return t
}

// FromMatrix returns a rank-2 tensor holding a copy of the entries of m.
func FromMatrix(m *Matrix) *Tensor {
	t, _ := NewWithData(m.data, m.rows, m.cols)
	return t
}

// ARange returns a rank-1 tensor with values start, start+step, … up to but not
// including stop. It returns an error if step is zero or points the wrong way.
func ARange(start, stop, step float64) (*Tensor, error) {
	if step == 0 {
		return nil, errors.New("tensornetwork: zero step")
	}
	var vals []float64
	if step > 0 {
		for x := start; x < stop; x += step {
			vals = append(vals, x)
		}
	} else {
		for x := start; x > stop; x += step {
			vals = append(vals, x)
		}
	}
	return FromVector(vals), nil
}

// LinSpace returns a rank-1 tensor of n evenly spaced values from start to stop
// inclusive. It returns an error if n < 1.
func LinSpace(start, stop float64, n int) (*Tensor, error) {
	if n < 1 {
		return nil, errors.New("tensornetwork: LinSpace needs n >= 1")
	}
	vals := make([]float64, n)
	if n == 1 {
		vals[0] = start
		return FromVector(vals), nil
	}
	step := (stop - start) / float64(n-1)
	for i := 0; i < n; i++ {
		vals[i] = start + float64(i)*step
	}
	return FromVector(vals), nil
}

// RandTensor returns a tensor of the given shape filled with independent
// standard-normal samples from a deterministic source seeded by seed.
func RandTensor(seed int64, shape ...int) *Tensor {
	r := rand.New(rand.NewSource(seed))
	t := New(shape...)
	for i := range t.data {
		t.data[i] = r.NormFloat64()
	}
	return t
}

// Shape returns a copy of the tensor's shape.
func (t *Tensor) Shape() []int { return append([]int(nil), t.shape...) }

// Rank returns the number of dimensions (axes) of the tensor.
func (t *Tensor) Rank() int { return len(t.shape) }

// Size returns the total number of elements in the tensor.
func (t *Tensor) Size() int { return len(t.data) }

// Dim returns the size of axis i.
func (t *Tensor) Dim(i int) int { return t.shape[i] }

// Strides returns a copy of the tensor's row-major strides.
func (t *Tensor) Strides() []int { return append([]int(nil), t.stride...) }

// Data returns the underlying row-major slice. Use [Tensor.Clone] for an
// independent copy.
func (t *Tensor) Data() []float64 { return t.data }

// Clone returns an independent deep copy of t.
func (t *Tensor) Clone() *Tensor {
	cp := make([]float64, len(t.data))
	copy(cp, t.data)
	return &Tensor{shape: append([]int(nil), t.shape...), stride: append([]int(nil), t.stride...), data: cp}
}

// flatIndex converts a multi-index into a flat offset using the tensor strides.
func (t *Tensor) flatIndex(idx []int) (int, error) {
	if len(idx) != len(t.shape) {
		return 0, fmt.Errorf("tensornetwork: index rank %d != tensor rank %d", len(idx), len(t.shape))
	}
	off := 0
	for a, i := range idx {
		if i < 0 || i >= t.shape[a] {
			return 0, fmt.Errorf("tensornetwork: index %d out of range for axis %d (size %d)", i, a, t.shape[a])
		}
		off += i * t.stride[a]
	}
	return off, nil
}

// At returns the element at the given multi-index. It panics if the index is
// out of range; use [Tensor.AtErr] for a checked version.
func (t *Tensor) At(idx ...int) float64 {
	off, err := t.flatIndex(idx)
	if err != nil {
		panic(err)
	}
	return t.data[off]
}

// AtErr returns the element at the given multi-index or an error if the index is
// invalid.
func (t *Tensor) AtErr(idx ...int) (float64, error) {
	off, err := t.flatIndex(idx)
	if err != nil {
		return 0, err
	}
	return t.data[off], nil
}

// Set stores v at the given multi-index. It panics if the index is out of range.
func (t *Tensor) Set(v float64, idx ...int) {
	off, err := t.flatIndex(idx)
	if err != nil {
		panic(err)
	}
	t.data[off] = v
}

// AtFlat returns the element at flat offset i in row-major order.
func (t *Tensor) AtFlat(i int) float64 { return t.data[i] }

// SetFlat stores v at flat offset i in row-major order.
func (t *Tensor) SetFlat(i int, v float64) { t.data[i] = v }

// MultiIndex converts a flat row-major offset into the corresponding
// multi-index.
func (t *Tensor) MultiIndex(flat int) []int {
	idx := make([]int, len(t.shape))
	for a := 0; a < len(t.shape); a++ {
		idx[a] = (flat / t.stride[a]) % t.shape[a]
	}
	return idx
}

// Equal reports whether t and o have the same shape and all corresponding
// entries differ by at most tol in absolute value.
func (t *Tensor) Equal(o *Tensor, tol float64) bool {
	if !shapeEqual(t.shape, o.shape) {
		return false
	}
	for i := range t.data {
		if math.Abs(t.data[i]-o.data[i]) > tol {
			return false
		}
	}
	return true
}

// shapeEqual reports whether two shapes are identical.
func shapeEqual(a, b []int) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

// String returns a short human-readable description of the tensor's shape and,
// for small tensors, its entries.
func (t *Tensor) String() string {
	var b strings.Builder
	fmt.Fprintf(&b, "Tensor(shape=%v", t.shape)
	if len(t.data) <= 16 {
		fmt.Fprintf(&b, ", data=%v", t.data)
	}
	b.WriteByte(')')
	return b.String()
}

// Reshape returns a tensor sharing t's data reinterpreted with the new shape,
// which must have the same total number of elements. Exactly one dimension may
// be -1, in which case it is inferred. It returns an error if the sizes are
// incompatible.
func (t *Tensor) Reshape(shape ...int) (*Tensor, error) {
	sh := append([]int(nil), shape...)
	neg := -1
	prod := 1
	for i, d := range sh {
		if d == -1 {
			if neg >= 0 {
				return nil, errors.New("tensornetwork: multiple -1 in reshape")
			}
			neg = i
			continue
		}
		if d < 0 {
			return nil, errors.New("tensornetwork: negative dimension in reshape")
		}
		prod *= d
	}
	if neg >= 0 {
		if prod == 0 || len(t.data)%prod != 0 {
			return nil, errors.New("tensornetwork: cannot infer reshape dimension")
		}
		sh[neg] = len(t.data) / prod
		prod *= sh[neg]
	}
	if prod != len(t.data) {
		return nil, fmt.Errorf("tensornetwork: reshape size %d != %d", prod, len(t.data))
	}
	cp := make([]float64, len(t.data))
	copy(cp, t.data)
	return &Tensor{shape: sh, stride: rowMajorStrides(sh), data: cp}, nil
}

// Ravel returns a rank-1 tensor with the same elements in row-major order.
func (t *Tensor) Ravel() *Tensor {
	r, _ := t.Reshape(len(t.data))
	return r
}

// Flatten is an alias for [Tensor.Ravel].
func (t *Tensor) Flatten() *Tensor { return t.Ravel() }

// Permute returns a tensor whose axes are reordered according to perm, a
// permutation of 0..rank-1. It returns an error if perm is not a valid
// permutation.
func (t *Tensor) Permute(perm ...int) (*Tensor, error) {
	n := len(t.shape)
	if len(perm) != n {
		return nil, fmt.Errorf("tensornetwork: permutation length %d != rank %d", len(perm), n)
	}
	seen := make([]bool, n)
	for _, p := range perm {
		if p < 0 || p >= n || seen[p] {
			return nil, errors.New("tensornetwork: invalid permutation")
		}
		seen[p] = true
	}
	newShape := make([]int, n)
	for i, p := range perm {
		newShape[i] = t.shape[p]
	}
	out := New(newShape...)
	idx := make([]int, n)
	src := make([]int, n)
	for flat := 0; flat < len(out.data); flat++ {
		rem := flat
		for a := 0; a < n; a++ {
			idx[a] = (rem / out.stride[a]) % newShape[a]
		}
		for a := 0; a < n; a++ {
			src[perm[a]] = idx[a]
		}
		off, _ := t.flatIndex(src)
		out.data[flat] = t.data[off]
	}
	return out, nil
}

// Transpose returns the tensor with its axes reversed. For a rank-2 tensor this
// is the ordinary matrix transpose.
func (t *Tensor) Transpose() *Tensor {
	n := len(t.shape)
	perm := make([]int, n)
	for i := 0; i < n; i++ {
		perm[i] = n - 1 - i
	}
	out, _ := t.Permute(perm...)
	return out
}

// SwapAxes returns the tensor with axes a and b exchanged.
func (t *Tensor) SwapAxes(a, b int) (*Tensor, error) {
	n := len(t.shape)
	if a < 0 || a >= n || b < 0 || b >= n {
		return nil, errors.New("tensornetwork: axis out of range")
	}
	perm := make([]int, n)
	for i := range perm {
		perm[i] = i
	}
	perm[a], perm[b] = perm[b], perm[a]
	return t.Permute(perm...)
}

// MoveAxis returns the tensor with the axis at position src moved to position
// dst, other axes keeping their relative order.
func (t *Tensor) MoveAxis(src, dst int) (*Tensor, error) {
	n := len(t.shape)
	if src < 0 || src >= n || dst < 0 || dst >= n {
		return nil, errors.New("tensornetwork: axis out of range")
	}
	order := make([]int, 0, n)
	for i := 0; i < n; i++ {
		if i != src {
			order = append(order, i)
		}
	}
	perm := make([]int, 0, n)
	perm = append(perm, order[:dst]...)
	perm = append(perm, src)
	perm = append(perm, order[dst:]...)
	return t.Permute(perm...)
}

// ExpandDims returns a tensor with a new axis of length 1 inserted at position
// axis.
func (t *Tensor) ExpandDims(axis int) (*Tensor, error) {
	n := len(t.shape)
	if axis < 0 || axis > n {
		return nil, errors.New("tensornetwork: axis out of range")
	}
	sh := make([]int, 0, n+1)
	sh = append(sh, t.shape[:axis]...)
	sh = append(sh, 1)
	sh = append(sh, t.shape[axis:]...)
	return t.Reshape(sh...)
}

// Squeeze returns a tensor with all length-1 axes removed. If axes are given,
// only those axes are removed and it is an error if any is not length 1.
func (t *Tensor) Squeeze(axes ...int) (*Tensor, error) {
	remove := make(map[int]bool)
	if len(axes) == 0 {
		for i, d := range t.shape {
			if d == 1 {
				remove[i] = true
			}
		}
	} else {
		for _, a := range axes {
			if a < 0 || a >= len(t.shape) {
				return nil, errors.New("tensornetwork: axis out of range")
			}
			if t.shape[a] != 1 {
				return nil, fmt.Errorf("tensornetwork: axis %d has size %d, cannot squeeze", a, t.shape[a])
			}
			remove[a] = true
		}
	}
	var sh []int
	for i, d := range t.shape {
		if !remove[i] {
			sh = append(sh, d)
		}
	}
	return t.Reshape(sh...)
}

// Concatenate joins tensors along the given axis. All tensors must have the same
// shape except along that axis. It returns an error on a shape mismatch or an
// empty list.
func Concatenate(axis int, tensors ...*Tensor) (*Tensor, error) {
	if len(tensors) == 0 {
		return nil, errors.New("tensornetwork: Concatenate needs at least one tensor")
	}
	first := tensors[0]
	n := len(first.shape)
	if axis < 0 || axis >= n {
		return nil, errors.New("tensornetwork: axis out of range")
	}
	outShape := append([]int(nil), first.shape...)
	total := 0
	for _, t := range tensors {
		if len(t.shape) != n {
			return nil, errors.New("tensornetwork: rank mismatch in Concatenate")
		}
		for a := 0; a < n; a++ {
			if a != axis && t.shape[a] != first.shape[a] {
				return nil, errors.New("tensornetwork: shape mismatch in Concatenate")
			}
		}
		total += t.shape[axis]
	}
	outShape[axis] = total
	out := New(outShape...)
	offset := 0
	idx := make([]int, n)
	for _, t := range tensors {
		for flat := 0; flat < len(t.data); flat++ {
			rem := flat
			for a := 0; a < n; a++ {
				idx[a] = (rem / t.stride[a]) % t.shape[a]
			}
			outIdx := append([]int(nil), idx...)
			outIdx[axis] += offset
			off, _ := out.flatIndex(outIdx)
			out.data[off] = t.data[flat]
		}
		offset += t.shape[axis]
	}
	return out, nil
}

// Stack joins tensors of identical shape along a new axis inserted at position
// axis. It returns an error on a shape mismatch or an empty list.
func Stack(axis int, tensors ...*Tensor) (*Tensor, error) {
	if len(tensors) == 0 {
		return nil, errors.New("tensornetwork: Stack needs at least one tensor")
	}
	expanded := make([]*Tensor, len(tensors))
	for i, t := range tensors {
		if !shapeEqual(t.shape, tensors[0].shape) {
			return nil, errors.New("tensornetwork: shape mismatch in Stack")
		}
		e, err := t.ExpandDims(axis)
		if err != nil {
			return nil, err
		}
		expanded[i] = e
	}
	return Concatenate(axis, expanded...)
}

// ToMatrix converts a rank-2 tensor into a [Matrix]. It returns an error if the
// tensor is not rank 2.
func (t *Tensor) ToMatrix() (*Matrix, error) {
	if len(t.shape) != 2 {
		return nil, fmt.Errorf("tensornetwork: ToMatrix requires rank 2, got %d", len(t.shape))
	}
	return NewMatrixData(t.shape[0], t.shape[1], t.data)
}
