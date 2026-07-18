package matrix

import (
	"math"
	"strings"

	"github.com/malcolmston/algebra"
)

// Vector is a one-dimensional slice of algebra.Expr components. In matrix
// products a Vector is treated as a column vector. Construct vectors with
// [NewVector] or [VectorFromInts].
type Vector struct {
	data []algebra.Expr
}

// NewVector returns a vector holding the simplified given components. A nil
// component is stored as 0.
func NewVector(components ...algebra.Expr) *Vector {
	out := make([]algebra.Expr, len(components))
	for i, c := range components {
		if c == nil {
			out[i] = zero()
		} else {
			out[i] = simp(c)
		}
	}
	return &Vector{data: out}
}

// VectorFromInts returns a vector built from int64 components.
func VectorFromInts(vals ...int64) *Vector {
	out := make([]algebra.Expr, len(vals))
	for i, v := range vals {
		out[i] = algebra.Int(v)
	}
	return &Vector{data: out}
}

// Len returns the number of components.
func (v *Vector) Len() int { return len(v.data) }

// At returns the i-th component (0-based). It panics if i is out of range.
func (v *Vector) At(i int) algebra.Expr {
	if i < 0 || i >= len(v.data) {
		panic("matrix: vector index out of range")
	}
	return v.data[i]
}

// Set stores the simplified value e at index i (0-based). A nil value stores 0.
func (v *Vector) Set(i int, e algebra.Expr) {
	if i < 0 || i >= len(v.data) {
		panic("matrix: vector index out of range")
	}
	if e == nil {
		v.data[i] = zero()
		return
	}
	v.data[i] = simp(e)
}

// Clone returns a copy of the vector.
func (v *Vector) Clone() *Vector {
	out := make([]algebra.Expr, len(v.data))
	copy(out, v.data)
	return &Vector{data: out}
}

// Equal reports whether v and w have the same length and structurally equal
// components (compared after simplification).
func (v *Vector) Equal(w *Vector) bool {
	if w == nil || len(v.data) != len(w.data) {
		return false
	}
	for i := range v.data {
		if !simp(v.data[i]).Equal(simp(w.data[i])) {
			return false
		}
	}
	return true
}

// String renders the vector as a bracketed, comma-separated list.
func (v *Vector) String() string {
	parts := make([]string, len(v.data))
	for i, e := range v.data {
		parts[i] = e.String()
	}
	return "[" + strings.Join(parts, ", ") + "]"
}

// ColMatrix returns the vector as an n×1 column matrix.
func (v *Vector) ColMatrix() *Matrix {
	m := New(len(v.data), 1)
	for i, e := range v.data {
		m.data[i][0] = e
	}
	return m
}

// RowMatrix returns the vector as a 1×n row matrix.
func (v *Vector) RowMatrix() *Matrix {
	m := New(1, len(v.data))
	copy(m.data[0], v.data)
	return m
}

// Add returns the component-wise sum v+w. It returns [ErrDimension] if the
// lengths differ.
func (v *Vector) Add(w *Vector) (*Vector, error) {
	if len(v.data) != len(w.data) {
		return nil, ErrDimension
	}
	out := make([]algebra.Expr, len(v.data))
	for i := range v.data {
		out[i] = simp(algebra.Add(v.data[i], w.data[i]))
	}
	return &Vector{data: out}, nil
}

// Sub returns the component-wise difference v-w. It returns [ErrDimension] if
// the lengths differ.
func (v *Vector) Sub(w *Vector) (*Vector, error) {
	if len(v.data) != len(w.data) {
		return nil, ErrDimension
	}
	out := make([]algebra.Expr, len(v.data))
	for i := range v.data {
		out[i] = simp(algebra.Add(v.data[i], algebra.Mul(algebra.Int(-1), w.data[i])))
	}
	return &Vector{data: out}, nil
}

// ScalarMul returns the vector scaled by the expression s.
func (v *Vector) ScalarMul(s algebra.Expr) *Vector {
	out := make([]algebra.Expr, len(v.data))
	for i := range v.data {
		out[i] = simp(algebra.Mul(s, v.data[i]))
	}
	return &Vector{data: out}
}

// Neg returns -v.
func (v *Vector) Neg() *Vector { return v.ScalarMul(algebra.Int(-1)) }

// Dot returns the dot product v·w as a single simplified expression. It returns
// [ErrDimension] if the lengths differ.
func (v *Vector) Dot(w *Vector) (algebra.Expr, error) {
	if len(v.data) != len(w.data) {
		return nil, ErrDimension
	}
	terms := make([]algebra.Expr, len(v.data))
	for i := range v.data {
		terms[i] = algebra.Mul(v.data[i], w.data[i])
	}
	return simp(algebra.Add(terms...)), nil
}

// Cross returns the 3-dimensional cross product v×w. Both vectors must have
// length 3; otherwise it returns [ErrDimension].
func (v *Vector) Cross(w *Vector) (*Vector, error) {
	if len(v.data) != 3 || len(w.data) != 3 {
		return nil, ErrDimension
	}
	a, b := v.data, w.data
	x := algebra.Add(algebra.Mul(a[1], b[2]), algebra.Mul(algebra.Int(-1), a[2], b[1]))
	y := algebra.Add(algebra.Mul(a[2], b[0]), algebra.Mul(algebra.Int(-1), a[0], b[2]))
	z := algebra.Add(algebra.Mul(a[0], b[1]), algebra.Mul(algebra.Int(-1), a[1], b[0]))
	return &Vector{data: []algebra.Expr{simp(x), simp(y), simp(z)}}, nil
}

// Norm returns the Euclidean norm sqrt(v·v) as an exact symbolic expression.
// The parent package's Sqrt pulls out perfect-square factors, so the norm of an
// integer vector is returned in simplest radical form.
func (v *Vector) Norm() algebra.Expr {
	d, _ := v.Dot(v)
	return simp(algebra.Sqrt(d))
}

// Normalize returns the unit vector v/‖v‖ with each component simplified. It
// returns [ErrSingular] if v is the zero vector.
func (v *Vector) Normalize() (*Vector, error) {
	n := v.Norm()
	if simp(n).Equal(zero()) {
		return nil, ErrSingular
	}
	inv := algebra.Pow(n, algebra.Int(-1))
	return v.ScalarMul(inv), nil
}

// Angle returns the angle between v and w in radians as a numeric float64,
// computed as arccos(v·w / (‖v‖‖w‖)). Both vectors must be numeric (contain no
// free symbols) and non-zero; otherwise it returns an error. The result is
// numeric because the parent package has no inverse trigonometric functions.
func (v *Vector) Angle(w *Vector) (float64, error) {
	if len(v.data) != len(w.data) {
		return 0, ErrDimension
	}
	dot, _ := v.Dot(w)
	dv, err := algebra.Evalf(dot)
	if err != nil {
		return 0, err
	}
	nv, err := algebra.Evalf(v.Norm())
	if err != nil {
		return 0, err
	}
	nw, err := algebra.Evalf(w.Norm())
	if err != nil {
		return 0, err
	}
	if nv == 0 || nw == 0 {
		return 0, ErrSingular
	}
	c := dv / (nv * nw)
	// Clamp to guard against tiny floating-point overshoot outside [-1, 1].
	if c > 1 {
		c = 1
	} else if c < -1 {
		c = -1
	}
	return math.Acos(c), nil
}
