package tensor

import (
	"errors"
	"fmt"
	"math"
	"strconv"
	"strings"
)

// Sentinel errors returned by operations in this package. Callers may test for
// them with [errors.Is].
var (
	// ErrShape reports that a shape is invalid (contains a non-positive
	// dimension) or that two shapes are incompatible for the requested
	// operation.
	ErrShape = errors.New("tensor: invalid or incompatible shape")
	// ErrIndex reports that a multi-index or flat index is out of range for a
	// tensor's shape.
	ErrIndex = errors.New("tensor: index out of range")
	// ErrAxis reports that an axis argument is out of range or that a set of
	// axes is not a valid permutation.
	ErrAxis = errors.New("tensor: invalid axis")
	// ErrRank reports that a tensor does not have the rank required by an
	// operation (for example a non-matrix passed to [MatMul]).
	ErrRank = errors.New("tensor: invalid rank for operation")
	// ErrDataLength reports that the length of a data slice does not match the
	// number of elements implied by a shape.
	ErrDataLength = errors.New("tensor: data length does not match shape")
	// ErrSpec reports a malformed Einstein-summation specification passed to
	// [Einsum].
	ErrSpec = errors.New("tensor: invalid einsum specification")
)

// Tensor is a dense, row-major, rank-r array of float64 values.
//
// The zero value is not usable; construct tensors with [New], [Zeros],
// [NewWithData] and the other constructors. A Tensor with an empty shape is a
// scalar holding a single value. Tensors are mutable through [Tensor.Set]; use
// [Tensor.Clone] to obtain an independent copy.
type Tensor struct {
	shape   []int
	strides []int
	data    []float64
}

// tensorProduct returns the product of the dimensions in shape, i.e. the number
// of elements a tensor of that shape holds. The product of the empty shape is 1
// (a scalar).
func tensorProduct(shape []int) int {
	n := 1
	for _, d := range shape {
		n *= d
	}
	return n
}

// tensorStrides returns the row-major strides for the given shape. The stride of
// the last axis is 1 and each earlier stride is the product of all later
// dimensions.
func tensorStrides(shape []int) []int {
	strides := make([]int, len(shape))
	acc := 1
	for i := len(shape) - 1; i >= 0; i-- {
		strides[i] = acc
		acc *= shape[i]
	}
	return strides
}

// tensorValidShape reports whether every dimension in shape is strictly
// positive. The empty shape (a scalar) is valid.
func tensorValidShape(shape []int) bool {
	for _, d := range shape {
		if d <= 0 {
			return false
		}
	}
	return true
}

// tensorCopyInts returns a fresh copy of s so callers cannot mutate a tensor's
// internal shape or strides through an aliased slice.
func tensorCopyInts(s []int) []int {
	out := make([]int, len(s))
	copy(out, s)
	return out
}

// New returns a new zero-filled tensor with the given shape. Passing no
// dimensions yields a scalar tensor. New panics if any dimension is
// non-positive; use [NewWithData] when a fallible constructor is required.
func New(shape ...int) *Tensor {
	if !tensorValidShape(shape) {
		panic(fmt.Errorf("%w: %v", ErrShape, shape))
	}
	s := tensorCopyInts(shape)
	return &Tensor{
		shape:   s,
		strides: tensorStrides(s),
		data:    make([]float64, tensorProduct(s)),
	}
}

// NewWithData returns a tensor with the given shape whose elements are the
// row-major contents of data. The data slice is copied. It returns [ErrShape]
// if the shape is invalid and [ErrDataLength] if len(data) does not equal the
// number of elements implied by shape.
func NewWithData(shape []int, data []float64) (*Tensor, error) {
	if !tensorValidShape(shape) {
		return nil, fmt.Errorf("%w: %v", ErrShape, shape)
	}
	if len(data) != tensorProduct(shape) {
		return nil, fmt.Errorf("%w: have %d, want %d", ErrDataLength, len(data), tensorProduct(shape))
	}
	s := tensorCopyInts(shape)
	d := make([]float64, len(data))
	copy(d, data)
	return &Tensor{shape: s, strides: tensorStrides(s), data: d}, nil
}

// Zeros returns a tensor of the given shape filled with zeros. It is a synonym
// for [New] and panics on an invalid shape.
func Zeros(shape ...int) *Tensor { return New(shape...) }

// Ones returns a tensor of the given shape with every element equal to 1. It
// panics on an invalid shape.
func Ones(shape ...int) *Tensor { return Full(1, shape...) }

// Full returns a tensor of the given shape with every element equal to value.
// It panics on an invalid shape.
func Full(value float64, shape ...int) *Tensor {
	t := New(shape...)
	for i := range t.data {
		t.data[i] = value
	}
	return t
}

// FromScalar returns a rank-0 tensor holding the single value v.
func FromScalar(v float64) *Tensor {
	return &Tensor{shape: []int{}, strides: []int{}, data: []float64{v}}
}

// FromVector returns a rank-1 tensor holding a copy of v. It panics if v is
// empty, since a length-0 axis is not a valid dimension.
func FromVector(v []float64) *Tensor {
	t := New(len(v))
	copy(t.data, v)
	return t
}

// FromMatrix returns a rank-2 tensor holding a copy of the row-major matrix m.
// It returns [ErrShape] if m is empty or if its rows have differing lengths.
func FromMatrix(m [][]float64) (*Tensor, error) {
	if len(m) == 0 || len(m[0]) == 0 {
		return nil, fmt.Errorf("%w: empty matrix", ErrShape)
	}
	cols := len(m[0])
	for _, row := range m {
		if len(row) != cols {
			return nil, fmt.Errorf("%w: ragged matrix", ErrShape)
		}
	}
	t := New(len(m), cols)
	for i, row := range m {
		copy(t.data[i*cols:(i+1)*cols], row)
	}
	return t, nil
}

// Rank returns the number of axes of t (0 for a scalar).
func (t *Tensor) Rank() int { return len(t.shape) }

// Shape returns a copy of t's dimensions, one entry per axis.
func (t *Tensor) Shape() []int { return tensorCopyInts(t.shape) }

// Strides returns a copy of t's row-major strides, one entry per axis.
func (t *Tensor) Strides() []int { return tensorCopyInts(t.strides) }

// Size returns the total number of elements in t (1 for a scalar).
func (t *Tensor) Size() int { return len(t.data) }

// Dim returns the length of the given axis. It panics if axis is out of range.
func (t *Tensor) Dim(axis int) int {
	if axis < 0 || axis >= len(t.shape) {
		panic(fmt.Errorf("%w: axis %d of rank-%d tensor", ErrAxis, axis, len(t.shape)))
	}
	return t.shape[axis]
}

// IsScalar reports whether t is a rank-0 tensor.
func (t *Tensor) IsScalar() bool { return len(t.shape) == 0 }

// Data returns a copy of t's underlying elements in row-major order. Mutating
// the result does not affect t.
func (t *Tensor) Data() []float64 {
	out := make([]float64, len(t.data))
	copy(out, t.data)
	return out
}

// flatIndex converts a multi-index to the position of that element in the
// backing slice. It panics with [ErrIndex] if idx has the wrong length or any
// component is out of range.
func (t *Tensor) flatIndex(idx []int) int {
	if len(idx) != len(t.shape) {
		panic(fmt.Errorf("%w: got %d indices for rank-%d tensor", ErrIndex, len(idx), len(t.shape)))
	}
	f := 0
	for i, v := range idx {
		if v < 0 || v >= t.shape[i] {
			panic(fmt.Errorf("%w: index %d on axis %d (dim %d)", ErrIndex, v, i, t.shape[i]))
		}
		f += v * t.strides[i]
	}
	return f
}

// At returns the element of t at the given multi-index. The number of indices
// must equal t.Rank(); for a scalar, call At with no arguments. At panics with
// [ErrIndex] on a wrong number of indices or an out-of-range component.
func (t *Tensor) At(idx ...int) float64 { return t.data[t.flatIndex(idx)] }

// Set stores value at the given multi-index of t, following the same index
// rules as [Tensor.At]. It panics with [ErrIndex] on invalid indices.
func (t *Tensor) Set(value float64, idx ...int) { t.data[t.flatIndex(idx)] = value }

// addAt accumulates value into the element at the given multi-index. It is an
// internal helper used by contraction and summation routines.
func (t *Tensor) addAt(value float64, idx []int) { t.data[t.flatIndex(idx)] += value }

// Clone returns an independent deep copy of t.
func (t *Tensor) Clone() *Tensor {
	return &Tensor{
		shape:   tensorCopyInts(t.shape),
		strides: tensorCopyInts(t.strides),
		data:    t.Data(),
	}
}

// ShapeEqual reports whether t and other have identical shapes.
func (t *Tensor) ShapeEqual(other *Tensor) bool {
	if len(t.shape) != len(other.shape) {
		return false
	}
	for i := range t.shape {
		if t.shape[i] != other.shape[i] {
			return false
		}
	}
	return true
}

// Equal reports whether t and other have the same shape and bit-for-bit equal
// elements. Use [Tensor.AlmostEqual] for floating-point comparisons with a
// tolerance.
func (t *Tensor) Equal(other *Tensor) bool {
	if !t.ShapeEqual(other) {
		return false
	}
	for i := range t.data {
		if t.data[i] != other.data[i] {
			return false
		}
	}
	return true
}

// AlmostEqual reports whether t and other have the same shape and every pair of
// corresponding elements differs by at most tol in absolute value. NaN elements
// never compare equal.
func (t *Tensor) AlmostEqual(other *Tensor, tol float64) bool {
	if !t.ShapeEqual(other) {
		return false
	}
	for i := range t.data {
		if math.Abs(t.data[i]-other.data[i]) > tol {
			return false
		}
	}
	return true
}

// ScalarValue returns the single element of a rank-0 tensor. It returns
// [ErrRank] if t is not a scalar.
func (t *Tensor) ScalarValue() (float64, error) {
	if !t.IsScalar() {
		return 0, fmt.Errorf("%w: ScalarValue requires rank 0, got rank %d", ErrRank, t.Rank())
	}
	return t.data[0], nil
}

// ToVector returns a copy of a rank-1 tensor's elements. It returns [ErrRank] if
// t is not a vector.
func (t *Tensor) ToVector() ([]float64, error) {
	if t.Rank() != 1 {
		return nil, fmt.Errorf("%w: ToVector requires rank 1, got rank %d", ErrRank, t.Rank())
	}
	return t.Data(), nil
}

// ToMatrix returns the elements of a rank-2 tensor as a freshly allocated
// row-major [][]float64. It returns [ErrRank] if t is not a matrix.
func (t *Tensor) ToMatrix() ([][]float64, error) {
	if t.Rank() != 2 {
		return nil, fmt.Errorf("%w: ToMatrix requires rank 2, got rank %d", ErrRank, t.Rank())
	}
	rows, cols := t.shape[0], t.shape[1]
	out := make([][]float64, rows)
	for i := 0; i < rows; i++ {
		out[i] = make([]float64, cols)
		copy(out[i], t.data[i*cols:(i+1)*cols])
	}
	return out, nil
}

// RavelIndex converts a multi-index into a flat row-major offset for the given
// shape. It returns [ErrIndex] if idx has the wrong length or is out of range,
// and [ErrShape] if shape is invalid.
func RavelIndex(shape, idx []int) (int, error) {
	if !tensorValidShape(shape) {
		return 0, fmt.Errorf("%w: %v", ErrShape, shape)
	}
	if len(idx) != len(shape) {
		return 0, fmt.Errorf("%w: got %d indices for rank-%d shape", ErrIndex, len(idx), len(shape))
	}
	strides := tensorStrides(shape)
	f := 0
	for i, v := range idx {
		if v < 0 || v >= shape[i] {
			return 0, fmt.Errorf("%w: index %d on axis %d", ErrIndex, v, i)
		}
		f += v * strides[i]
	}
	return f, nil
}

// UnravelIndex converts a flat row-major offset into the corresponding
// multi-index for the given shape. It returns [ErrShape] if shape is invalid and
// [ErrIndex] if flat is out of range.
func UnravelIndex(shape []int, flat int) ([]int, error) {
	if !tensorValidShape(shape) {
		return nil, fmt.Errorf("%w: %v", ErrShape, shape)
	}
	total := tensorProduct(shape)
	if flat < 0 || flat >= total {
		return nil, fmt.Errorf("%w: flat index %d for size %d", ErrIndex, flat, total)
	}
	return tensorUnravel(flat, shape), nil
}

// tensorUnravel converts a flat row-major offset into a multi-index for shape.
// It assumes flat is in range and shape is valid.
func tensorUnravel(flat int, shape []int) []int {
	idx := make([]int, len(shape))
	for i := len(shape) - 1; i >= 0; i-- {
		idx[i] = flat % shape[i]
		flat /= shape[i]
	}
	return idx
}

// String returns a compact, deterministic textual description of t giving its
// shape and, for small tensors, its elements.
func (t *Tensor) String() string {
	var b strings.Builder
	b.WriteString("Tensor(shape=[")
	for i, d := range t.shape {
		if i > 0 {
			b.WriteByte(',')
		}
		b.WriteString(strconv.Itoa(d))
	}
	b.WriteString("], data=[")
	for i, v := range t.data {
		if i > 0 {
			b.WriteByte(' ')
		}
		if i >= 12 {
			b.WriteString("...")
			break
		}
		b.WriteString(strconv.FormatFloat(v, 'g', -1, 64))
	}
	b.WriteString("])")
	return b.String()
}
