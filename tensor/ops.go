package tensor

import (
	"fmt"
	"math"
)

// tensorBinary applies f elementwise to two tensors of identical shape and
// returns the result. It returns [ErrShape] if the shapes differ.
func tensorBinary(a, b *Tensor, f func(x, y float64) float64) (*Tensor, error) {
	if !a.ShapeEqual(b) {
		return nil, fmt.Errorf("%w: elementwise op needs equal shapes, got %v and %v", ErrShape, a.shape, b.shape)
	}
	out := New(a.shape...)
	for i := range a.data {
		out.data[i] = f(a.data[i], b.data[i])
	}
	return out, nil
}

// Add returns the elementwise sum t+other. It returns [ErrShape] if the shapes
// differ.
func (t *Tensor) Add(other *Tensor) (*Tensor, error) {
	return tensorBinary(t, other, func(x, y float64) float64 { return x + y })
}

// Sub returns the elementwise difference t-other. It returns [ErrShape] if the
// shapes differ.
func (t *Tensor) Sub(other *Tensor) (*Tensor, error) {
	return tensorBinary(t, other, func(x, y float64) float64 { return x - y })
}

// Mul returns the elementwise (Hadamard) product t*other. It returns [ErrShape]
// if the shapes differ. For contraction-style products see [MatMul], [Outer] and
// [TensorDot].
func (t *Tensor) Mul(other *Tensor) (*Tensor, error) {
	return tensorBinary(t, other, func(x, y float64) float64 { return x * y })
}

// Div returns the elementwise quotient t/other. Division by zero follows IEEE
// 754 semantics (producing infinities or NaN). It returns [ErrShape] if the
// shapes differ.
func (t *Tensor) Div(other *Tensor) (*Tensor, error) {
	return tensorBinary(t, other, func(x, y float64) float64 { return x / y })
}

// Neg returns a tensor with every element of t negated.
func (t *Tensor) Neg() *Tensor { return t.Scale(-1) }

// Scale returns a tensor with every element of t multiplied by s.
func (t *Tensor) Scale(s float64) *Tensor {
	out := New(t.shape...)
	for i := range t.data {
		out.data[i] = t.data[i] * s
	}
	return out
}

// AddScalar returns a tensor with s added to every element of t.
func (t *Tensor) AddScalar(s float64) *Tensor {
	out := New(t.shape...)
	for i := range t.data {
		out.data[i] = t.data[i] + s
	}
	return out
}

// Apply returns a tensor obtained by applying f to every element of t.
func (t *Tensor) Apply(f func(float64) float64) *Tensor {
	out := New(t.shape...)
	for i := range t.data {
		out.data[i] = f(t.data[i])
	}
	return out
}

// Abs returns a tensor holding the absolute value of every element of t.
func (t *Tensor) Abs() *Tensor { return t.Apply(math.Abs) }

// Sum returns the sum of all elements of t.
func (t *Tensor) Sum() float64 {
	s := 0.0
	for _, v := range t.data {
		s += v
	}
	return s
}

// Mean returns the arithmetic mean of all elements of t.
func (t *Tensor) Mean() float64 { return t.Sum() / float64(t.Size()) }

// Product returns the product of all elements of t.
func (t *Tensor) Product() float64 {
	p := 1.0
	for _, v := range t.data {
		p *= v
	}
	return p
}

// Max returns the largest element of t.
func (t *Tensor) Max() float64 {
	m := t.data[0]
	for _, v := range t.data[1:] {
		if v > m {
			m = v
		}
	}
	return m
}

// Min returns the smallest element of t.
func (t *Tensor) Min() float64 {
	m := t.data[0]
	for _, v := range t.data[1:] {
		if v < m {
			m = v
		}
	}
	return m
}

// Norm returns the Frobenius (Euclidean) norm of t, the square root of the sum
// of squared elements.
func (t *Tensor) Norm() float64 {
	s := 0.0
	for _, v := range t.data {
		s += v * v
	}
	return math.Sqrt(s)
}

// SumAxis returns a tensor with the given axis summed out, reducing the rank by
// one. Summing the sole axis of a vector yields a scalar. It returns [ErrAxis]
// if axis is out of range.
func (t *Tensor) SumAxis(axis int) (*Tensor, error) {
	r := t.Rank()
	if axis < 0 || axis >= r {
		return nil, fmt.Errorf("%w: axis %d for rank-%d tensor", ErrAxis, axis, r)
	}
	newShape := make([]int, 0, r-1)
	newShape = append(newShape, t.shape[:axis]...)
	newShape = append(newShape, t.shape[axis+1:]...)
	out := New(newShape...)
	oidx := make([]int, len(newShape))
	for f := 0; f < t.Size(); f++ {
		idx := tensorUnravel(f, t.shape)
		copy(oidx, idx[:axis])
		copy(oidx[axis:], idx[axis+1:])
		out.addAt(t.data[f], oidx)
	}
	return out, nil
}
