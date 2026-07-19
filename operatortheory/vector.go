package operatortheory

import (
	"math"
	"math/cmplx"
)

// Vector is a finite-dimensional complex vector, an element of the Hilbert
// space C^n with inner product <x,y> = sum conj(x_i) y_i.
type Vector []complex128

// NewVector returns a zero vector of length n. It panics if n is negative.
func NewVector(n int) Vector {
	if n < 0 {
		panic("operatortheory: negative vector length")
	}
	return make(Vector, n)
}

// VectorFromReal returns a complex vector whose entries are the given real
// numbers with zero imaginary part.
func VectorFromReal(data []float64) Vector {
	v := make(Vector, len(data))
	for i, x := range data {
		v[i] = complex(x, 0)
	}
	return v
}

// VectorOf returns a Vector containing the supplied entries.
func VectorOf(entries ...complex128) Vector {
	v := make(Vector, len(entries))
	copy(v, entries)
	return v
}

// BasisVector returns the i-th standard basis vector of C^n. It panics if i is
// out of range.
func BasisVector(n, i int) Vector {
	if i < 0 || i >= n {
		panic("operatortheory: basis index out of range")
	}
	v := make(Vector, n)
	v[i] = 1
	return v
}

// Len returns the dimension of the vector.
func (v Vector) Len() int { return len(v) }

// Clone returns an independent copy of v.
func (v Vector) Clone() Vector {
	w := make(Vector, len(v))
	copy(w, v)
	return w
}

// Conjugate returns the entrywise complex conjugate of v.
func (v Vector) Conjugate() Vector {
	w := make(Vector, len(v))
	for i, z := range v {
		w[i] = cmplx.Conj(z)
	}
	return w
}

// Add returns v + w. It panics if the lengths differ.
func (v Vector) Add(w Vector) Vector {
	requireSameLen(v, w)
	r := make(Vector, len(v))
	for i := range v {
		r[i] = v[i] + w[i]
	}
	return r
}

// Sub returns v - w. It panics if the lengths differ.
func (v Vector) Sub(w Vector) Vector {
	requireSameLen(v, w)
	r := make(Vector, len(v))
	for i := range v {
		r[i] = v[i] - w[i]
	}
	return r
}

// Scale returns the vector s*v.
func (v Vector) Scale(s complex128) Vector {
	r := make(Vector, len(v))
	for i := range v {
		r[i] = s * v[i]
	}
	return r
}

// Neg returns -v.
func (v Vector) Neg() Vector { return v.Scale(-1) }

// Dot returns the Hermitian inner product <v,w> = sum conj(v_i) w_i. It is
// conjugate-linear in its first argument and linear in the second. It panics if
// the lengths differ.
func (v Vector) Dot(w Vector) complex128 {
	requireSameLen(v, w)
	var s complex128
	for i := range v {
		s += cmplx.Conj(v[i]) * w[i]
	}
	return s
}

// Norm returns the Euclidean (l2) norm of v.
func (v Vector) Norm() float64 {
	var s float64
	for _, z := range v {
		s += real(z)*real(z) + imag(z)*imag(z)
	}
	return math.Sqrt(s)
}

// NormSquared returns the squared Euclidean norm <v,v>.
func (v Vector) NormSquared() float64 {
	var s float64
	for _, z := range v {
		s += real(z)*real(z) + imag(z)*imag(z)
	}
	return s
}

// OneNorm returns the l1 norm, the sum of the moduli of the entries.
func (v Vector) OneNorm() float64 {
	var s float64
	for _, z := range v {
		s += cmplx.Abs(z)
	}
	return s
}

// InfNorm returns the l-infinity norm, the largest modulus among the entries.
func (v Vector) InfNorm() float64 {
	var m float64
	for _, z := range v {
		if a := cmplx.Abs(z); a > m {
			m = a
		}
	}
	return m
}

// Normalize returns a unit vector in the direction of v together with the
// original norm. If v is the zero vector it returns a copy of v and a zero
// norm.
func (v Vector) Normalize() (Vector, float64) {
	n := v.Norm()
	if n == 0 {
		return v.Clone(), 0
	}
	return v.Scale(complex(1/n, 0)), n
}

// IsZero reports whether every entry of v has modulus at most tol.
func (v Vector) IsZero(tol float64) bool {
	for _, z := range v {
		if cmplx.Abs(z) > tol {
			return false
		}
	}
	return true
}

// Equal reports whether v and w have the same length and agree entrywise to
// within tol.
func (v Vector) Equal(w Vector, tol float64) bool {
	if len(v) != len(w) {
		return false
	}
	for i := range v {
		if cmplx.Abs(v[i]-w[i]) > tol {
			return false
		}
	}
	return true
}

// IsOrthogonal reports whether v and w are orthogonal to within tol, i.e.
// |<v,w>| <= tol.
func (v Vector) IsOrthogonal(w Vector, tol float64) bool {
	return cmplx.Abs(v.Dot(w)) <= tol
}

// Angle returns the angle in radians between the real directions of v and w,
// computed from |<v,w>| / (||v|| ||w||). It lies in [0, pi/2].
func (v Vector) Angle(w Vector) float64 {
	nv, nw := v.Norm(), w.Norm()
	if nv == 0 || nw == 0 {
		return 0
	}
	c := cmplx.Abs(v.Dot(w)) / (nv * nw)
	if c > 1 {
		c = 1
	}
	return math.Acos(c)
}

// ProjectOnto returns the orthogonal projection of v onto the line spanned by
// the nonzero vector u, namely (<u,v>/<u,u>) u.
func (v Vector) ProjectOnto(u Vector) Vector {
	d := u.Dot(u)
	if d == 0 {
		return NewVector(len(v))
	}
	return u.Scale(u.Dot(v) / d)
}

// OuterProduct returns the rank-one operator v w^H, whose (i,j) entry is
// v_i conj(w_j).
func (v Vector) OuterProduct(w Vector) *Matrix {
	m := NewMatrix(len(v), len(w))
	for i := range v {
		for j := range w {
			m.data[i*m.cols+j] = v[i] * cmplx.Conj(w[j])
		}
	}
	return m
}

// GramSchmidt orthonormalises the supplied vectors using the modified
// Gram-Schmidt process with the Hermitian inner product. Vectors that are
// (numerically) linearly dependent on the earlier ones are dropped, so the
// result is an orthonormal basis for the span of the inputs.
func GramSchmidt(vectors []Vector, tol float64) []Vector {
	var out []Vector
	for _, v := range vectors {
		u := v.Clone()
		for _, q := range out {
			u = u.Sub(q.Scale(q.Dot(u)))
		}
		n := u.Norm()
		if n > tol {
			out = append(out, u.Scale(complex(1/n, 0)))
		}
	}
	return out
}

func requireSameLen(v, w Vector) {
	if len(v) != len(w) {
		panic("operatortheory: vector length mismatch")
	}
}
