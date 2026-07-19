package lattices

import (
	"fmt"
	"math"
	"strings"
)

// Vec is a real vector represented as a slice of float64 coordinates. It is the
// numeric vector type used by the reduction and enumeration routines.
type Vec []float64

// NewVec returns a Vec with the given coordinates.
func NewVec(xs ...float64) Vec {
	v := make(Vec, len(xs))
	copy(v, xs)
	return v
}

// ZeroVec returns the n-dimensional zero vector.
func ZeroVec(n int) Vec {
	return make(Vec, n)
}

// UnitVec returns the n-dimensional standard basis vector e_i (a 1 in position
// i and zeros elsewhere). It panics if i is out of range.
func UnitVec(n, i int) Vec {
	v := make(Vec, n)
	v[i] = 1
	return v
}

// VecFromInts builds a Vec from integer coordinates.
func VecFromInts(xs ...int64) Vec {
	v := make(Vec, len(xs))
	for i, x := range xs {
		v[i] = float64(x)
	}
	return v
}

// Dim returns the number of coordinates in v.
func (v Vec) Dim() int { return len(v) }

// Clone returns an independent copy of v.
func (v Vec) Clone() Vec {
	w := make(Vec, len(v))
	copy(w, v)
	return w
}

// Add returns the sum v+w. It panics if the dimensions differ.
func (v Vec) Add(w Vec) Vec {
	v.mustMatch(w)
	r := make(Vec, len(v))
	for i := range v {
		r[i] = v[i] + w[i]
	}
	return r
}

// Sub returns the difference v-w. It panics if the dimensions differ.
func (v Vec) Sub(w Vec) Vec {
	v.mustMatch(w)
	r := make(Vec, len(v))
	for i := range v {
		r[i] = v[i] - w[i]
	}
	return r
}

// Scale returns the scalar multiple s*v.
func (v Vec) Scale(s float64) Vec {
	r := make(Vec, len(v))
	for i := range v {
		r[i] = s * v[i]
	}
	return r
}

// Neg returns -v.
func (v Vec) Neg() Vec { return v.Scale(-1) }

// AddScaled returns v + s*w. It panics if the dimensions differ.
func (v Vec) AddScaled(w Vec, s float64) Vec {
	v.mustMatch(w)
	r := make(Vec, len(v))
	for i := range v {
		r[i] = v[i] + s*w[i]
	}
	return r
}

// Dot returns the Euclidean inner product <v, w>. It panics if the dimensions
// differ.
func (v Vec) Dot(w Vec) float64 {
	v.mustMatch(w)
	var s float64
	for i := range v {
		s += v[i] * w[i]
	}
	return s
}

// Norm2 returns the squared Euclidean norm <v, v>.
func (v Vec) Norm2() float64 { return v.Dot(v) }

// Norm returns the Euclidean (L2) norm of v.
func (v Vec) Norm() float64 { return math.Sqrt(v.Norm2()) }

// L1 returns the L1 (taxicab) norm of v, the sum of absolute values.
func (v Vec) L1() float64 {
	var s float64
	for _, x := range v {
		s += math.Abs(x)
	}
	return s
}

// LInf returns the L-infinity (maximum absolute value) norm of v.
func (v Vec) LInf() float64 {
	var m float64
	for _, x := range v {
		if a := math.Abs(x); a > m {
			m = a
		}
	}
	return m
}

// Sum returns the sum of the coordinates of v.
func (v Vec) Sum() float64 {
	var s float64
	for _, x := range v {
		s += x
	}
	return s
}

// Dist2 returns the squared Euclidean distance between v and w.
func (v Vec) Dist2(w Vec) float64 {
	v.mustMatch(w)
	var s float64
	for i := range v {
		d := v[i] - w[i]
		s += d * d
	}
	return s
}

// Dist returns the Euclidean distance between v and w.
func (v Vec) Dist(w Vec) float64 { return math.Sqrt(v.Dist2(w)) }

// IsZero reports whether every coordinate of v is exactly zero.
func (v Vec) IsZero() bool {
	for _, x := range v {
		if x != 0 {
			return false
		}
	}
	return true
}

// Equal reports whether v and w have identical dimensions and coordinates.
func (v Vec) Equal(w Vec) bool {
	if len(v) != len(w) {
		return false
	}
	for i := range v {
		if v[i] != w[i] {
			return false
		}
	}
	return true
}

// ApproxEqual reports whether v and w have the same dimension and every
// coordinate agrees within the absolute tolerance tol.
func (v Vec) ApproxEqual(w Vec, tol float64) bool {
	if len(v) != len(w) {
		return false
	}
	for i := range v {
		if math.Abs(v[i]-w[i]) > tol {
			return false
		}
	}
	return true
}

// ProjCoeff returns the scalar <v, w>/<w, w>, the coefficient of the orthogonal
// projection of v onto w. It returns 0 when w is the zero vector.
func (v Vec) ProjCoeff(w Vec) float64 {
	d := w.Norm2()
	if d == 0 {
		return 0
	}
	return v.Dot(w) / d
}

// Proj returns the orthogonal projection of v onto the line spanned by w. It
// returns the zero vector when w is the zero vector.
func (v Vec) Proj(w Vec) Vec {
	return w.Scale(v.ProjCoeff(w))
}

// Reject returns the component of v orthogonal to w, that is v minus its
// projection onto w.
func (v Vec) Reject(w Vec) Vec {
	return v.Sub(v.Proj(w))
}

// Normalize returns the unit vector in the direction of v. It returns a copy of
// v (the zero vector) when v has zero norm.
func (v Vec) Normalize() Vec {
	n := v.Norm()
	if n == 0 {
		return v.Clone()
	}
	return v.Scale(1 / n)
}

// Cosine returns the cosine of the angle between v and w. It returns 0 when
// either vector is zero.
func (v Vec) Cosine(w Vec) float64 {
	d := v.Norm() * w.Norm()
	if d == 0 {
		return 0
	}
	c := v.Dot(w) / d
	if c > 1 {
		c = 1
	} else if c < -1 {
		c = -1
	}
	return c
}

// Angle returns the angle in radians between v and w in the range [0, pi].
func (v Vec) Angle(w Vec) float64 {
	return math.Acos(v.Cosine(w))
}

// IsOrthogonal reports whether <v, w> is zero within the tolerance tol.
func (v Vec) IsOrthogonal(w Vec, tol float64) bool {
	return math.Abs(v.Dot(w)) <= tol
}

// Round returns the vector obtained by rounding each coordinate of v to the
// nearest integer (halves away from zero).
func (v Vec) Round() Vec {
	r := make(Vec, len(v))
	for i := range v {
		r[i] = math.Round(v[i])
	}
	return r
}

// String renders v as a bracketed, space-separated list of coordinates.
func (v Vec) String() string {
	parts := make([]string, len(v))
	for i, x := range v {
		parts[i] = fmt.Sprintf("%g", x)
	}
	return "[" + strings.Join(parts, " ") + "]"
}

// Combine returns the linear combination sum_i coeffs[i]*vecs[i]. It panics if
// coeffs and vecs have different lengths, and returns nil for an empty input.
func Combine(coeffs []float64, vecs []Vec) Vec {
	if len(coeffs) != len(vecs) {
		panic("lattices: Combine length mismatch")
	}
	if len(vecs) == 0 {
		return nil
	}
	r := ZeroVec(len(vecs[0]))
	for i, c := range coeffs {
		r = r.AddScaled(vecs[i], c)
	}
	return r
}

// IntCombine returns the integer linear combination sum_i coeffs[i]*vecs[i]. It
// panics if the lengths differ and returns nil for an empty input.
func IntCombine(coeffs []int64, vecs []Vec) Vec {
	if len(coeffs) != len(vecs) {
		panic("lattices: IntCombine length mismatch")
	}
	if len(vecs) == 0 {
		return nil
	}
	r := ZeroVec(len(vecs[0]))
	for i, c := range coeffs {
		r = r.AddScaled(vecs[i], float64(c))
	}
	return r
}

func (v Vec) mustMatch(w Vec) {
	if len(v) != len(w) {
		panic(ErrDimMismatch)
	}
}
