package tensor

import "fmt"

// Reshape returns a new tensor with the given shape and the same elements in
// row-major order. The product of the new dimensions must equal t.Size(). It
// returns [ErrShape] if the shape is invalid or the sizes do not match. The
// result is independent of t.
func (t *Tensor) Reshape(shape ...int) (*Tensor, error) {
	if !tensorValidShape(shape) {
		return nil, fmt.Errorf("%w: %v", ErrShape, shape)
	}
	if tensorProduct(shape) != t.Size() {
		return nil, fmt.Errorf("%w: cannot reshape size %d into %v", ErrShape, t.Size(), shape)
	}
	s := tensorCopyInts(shape)
	return &Tensor{shape: s, strides: tensorStrides(s), data: t.Data()}, nil
}

// Ravel returns a rank-1 tensor holding t's elements in row-major order. It is
// [Tensor.Reshape] to a single axis and never fails.
func (t *Tensor) Ravel() *Tensor {
	out, _ := t.Reshape(t.Size())
	return out
}

// Flatten returns a copy of t's elements in row-major order as a plain slice. It
// is a synonym for [Tensor.Data].
func (t *Tensor) Flatten() []float64 { return t.Data() }

// tensorIsPerm reports whether perm is a permutation of 0..n-1.
func tensorIsPerm(perm []int, n int) bool {
	if len(perm) != n {
		return false
	}
	seen := make([]bool, n)
	for _, p := range perm {
		if p < 0 || p >= n || seen[p] {
			return false
		}
		seen[p] = true
	}
	return true
}

// Permute returns a new tensor whose axes are a reordering of t's axes: axis i
// of the result is axis perm[i] of t. The perm argument must be a permutation of
// 0..t.Rank()-1. It returns [ErrAxis] otherwise. The result is materialised as a
// contiguous row-major tensor.
func (t *Tensor) Permute(perm ...int) (*Tensor, error) {
	r := t.Rank()
	if !tensorIsPerm(perm, r) {
		return nil, fmt.Errorf("%w: %v is not a permutation of %d axes", ErrAxis, perm, r)
	}
	newShape := make([]int, r)
	for i, p := range perm {
		newShape[i] = t.shape[p]
	}
	if r == 0 {
		return t.Clone(), nil
	}
	out := New(newShape...)
	src := make([]int, r)
	for f := 0; f < out.Size(); f++ {
		oidx := tensorUnravel(f, newShape)
		for i, p := range perm {
			src[p] = oidx[i]
		}
		out.data[f] = t.At(src...)
	}
	return out, nil
}

// Transpose returns a new tensor with all axes reversed. For a matrix this is
// the ordinary matrix transpose; for a scalar or vector it returns a copy.
func (t *Tensor) Transpose() *Tensor {
	r := t.Rank()
	perm := make([]int, r)
	for i := range perm {
		perm[i] = r - 1 - i
	}
	out, _ := t.Permute(perm...)
	return out
}

// SwapAxes returns a new tensor identical to t but with axes a and b exchanged.
// It returns [ErrAxis] if either axis is out of range.
func (t *Tensor) SwapAxes(a, b int) (*Tensor, error) {
	r := t.Rank()
	if a < 0 || a >= r || b < 0 || b >= r {
		return nil, fmt.Errorf("%w: cannot swap axes %d and %d of rank-%d tensor", ErrAxis, a, b, r)
	}
	perm := make([]int, r)
	for i := range perm {
		perm[i] = i
	}
	perm[a], perm[b] = b, a
	return t.Permute(perm...)
}

// MoveAxis returns a new tensor with the axis at position src moved to position
// dst, the remaining axes keeping their relative order. It returns [ErrAxis] if
// either position is out of range.
func (t *Tensor) MoveAxis(src, dst int) (*Tensor, error) {
	r := t.Rank()
	if src < 0 || src >= r || dst < 0 || dst >= r {
		return nil, fmt.Errorf("%w: cannot move axis %d to %d of rank-%d tensor", ErrAxis, src, dst, r)
	}
	order := make([]int, 0, r)
	for i := 0; i < r; i++ {
		if i != src {
			order = append(order, i)
		}
	}
	perm := make([]int, 0, r)
	perm = append(perm, order[:dst]...)
	perm = append(perm, src)
	perm = append(perm, order[dst:]...)
	return t.Permute(perm...)
}

// ExpandDims returns a new tensor with a length-1 axis inserted at the given
// position, increasing the rank by one. axis may equal t.Rank() to append the
// new axis. It returns [ErrAxis] if axis is out of range.
func (t *Tensor) ExpandDims(axis int) (*Tensor, error) {
	r := t.Rank()
	if axis < 0 || axis > r {
		return nil, fmt.Errorf("%w: cannot insert axis at %d of rank-%d tensor", ErrAxis, axis, r)
	}
	newShape := make([]int, 0, r+1)
	newShape = append(newShape, t.shape[:axis]...)
	newShape = append(newShape, 1)
	newShape = append(newShape, t.shape[axis:]...)
	return t.Reshape(newShape...)
}

// Squeeze returns a new tensor with every length-1 axis removed. A tensor whose
// axes are all length 1 (or a scalar) squeezes to a rank-0 scalar. The result
// shares no storage with t.
func (t *Tensor) Squeeze() *Tensor {
	newShape := make([]int, 0, t.Rank())
	for _, d := range t.shape {
		if d != 1 {
			newShape = append(newShape, d)
		}
	}
	out, _ := t.Reshape(newShape...)
	return out
}

// tensorSameExcept reports whether shapes a and b are equal on every axis other
// than axis.
func tensorSameExcept(a, b []int, axis int) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if i == axis {
			continue
		}
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

// Concatenate joins tensors along an existing axis. Every input must have the
// same rank and identical dimensions on all axes other than axis; the result's
// length along axis is the sum of the inputs'. It returns [ErrShape] on a
// mismatch, [ErrAxis] on a bad axis, and [ErrRank] if no tensors are given.
func Concatenate(axis int, tensors ...*Tensor) (*Tensor, error) {
	if len(tensors) == 0 {
		return nil, fmt.Errorf("%w: Concatenate needs at least one tensor", ErrRank)
	}
	first := tensors[0]
	r := first.Rank()
	if axis < 0 || axis >= r {
		return nil, fmt.Errorf("%w: axis %d for rank-%d tensor", ErrAxis, axis, r)
	}
	total := 0
	for _, t := range tensors {
		if !tensorSameExcept(first.shape, t.shape, axis) {
			return nil, fmt.Errorf("%w: shapes %v and %v cannot be concatenated on axis %d", ErrShape, first.shape, t.shape, axis)
		}
		total += t.shape[axis]
	}
	outShape := tensorCopyInts(first.shape)
	outShape[axis] = total
	out := New(outShape...)
	idx := make([]int, r)
	offset := 0
	for _, t := range tensors {
		for f := 0; f < t.Size(); f++ {
			src := tensorUnravel(f, t.shape)
			copy(idx, src)
			idx[axis] = src[axis] + offset
			out.Set(t.data[f], idx...)
		}
		offset += t.shape[axis]
	}
	return out, nil
}

// Stack joins tensors along a new axis inserted at position axis, increasing the
// rank by one. Every input must have an identical shape. axis may range from 0
// to the common rank inclusive. It returns [ErrShape] on a shape mismatch,
// [ErrAxis] on a bad axis, and [ErrRank] if no tensors are given.
func Stack(axis int, tensors ...*Tensor) (*Tensor, error) {
	if len(tensors) == 0 {
		return nil, fmt.Errorf("%w: Stack needs at least one tensor", ErrRank)
	}
	first := tensors[0]
	for _, t := range tensors {
		if !first.ShapeEqual(t) {
			return nil, fmt.Errorf("%w: Stack requires identical shapes, got %v and %v", ErrShape, first.shape, t.shape)
		}
	}
	expanded := make([]*Tensor, len(tensors))
	for i, t := range tensors {
		e, err := t.ExpandDims(axis)
		if err != nil {
			return nil, err
		}
		expanded[i] = e
	}
	return Concatenate(axis, expanded...)
}
