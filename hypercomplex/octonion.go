package hypercomplex

import "math"

// Octonion is an element of the eight-dimensional real algebra of octonions,
// stored as the coefficients E0..E7 of the basis units e0=1, e1..e7. The
// octonions form a normed division algebra that is neither commutative nor
// associative, but is alternative (any two elements generate an associative
// subalgebra).
type Octonion struct {
	E0, E1, E2, E3, E4, E5, E6, E7 float64
}

// Oct constructs an octonion from its eight real coefficients.
func Oct(e0, e1, e2, e3, e4, e5, e6, e7 float64) Octonion {
	return Octonion{e0, e1, e2, e3, e4, e5, e6, e7}
}

// IdentityOct returns the multiplicative identity octonion e0 = 1.
func IdentityOct() Octonion {
	return Octonion{E0: 1}
}

// OctFromScalar returns the octonion whose real part is s and whose other seven
// coefficients are zero.
func OctFromScalar(s float64) Octonion {
	return Octonion{E0: s}
}

// hypercomplexOctToPair splits an octonion into the pair (a, b) of quaternions
// used by the Cayley-Dickson construction, with a = (E0,E1,E2,E3) and
// b = (E4,E5,E6,E7).
func hypercomplexOctToPair(o Octonion) (a, b Quaternion) {
	a = Quaternion{o.E0, o.E1, o.E2, o.E3}
	b = Quaternion{o.E4, o.E5, o.E6, o.E7}
	return a, b
}

// hypercomplexPairToOct reassembles an octonion from a Cayley-Dickson pair of
// quaternions.
func hypercomplexPairToOct(a, b Quaternion) Octonion {
	return Octonion{a.W, a.X, a.Y, a.Z, b.W, b.X, b.Y, b.Z}
}

// Add returns the component-wise sum o + p.
func (o Octonion) Add(p Octonion) Octonion {
	return Octonion{
		o.E0 + p.E0, o.E1 + p.E1, o.E2 + p.E2, o.E3 + p.E3,
		o.E4 + p.E4, o.E5 + p.E5, o.E6 + p.E6, o.E7 + p.E7,
	}
}

// Sub returns the component-wise difference o - p.
func (o Octonion) Sub(p Octonion) Octonion {
	return Octonion{
		o.E0 - p.E0, o.E1 - p.E1, o.E2 - p.E2, o.E3 - p.E3,
		o.E4 - p.E4, o.E5 - p.E5, o.E6 - p.E6, o.E7 - p.E7,
	}
}

// Scale returns o with every coefficient multiplied by the real factor s.
func (o Octonion) Scale(s float64) Octonion {
	return Octonion{
		o.E0 * s, o.E1 * s, o.E2 * s, o.E3 * s,
		o.E4 * s, o.E5 * s, o.E6 * s, o.E7 * s,
	}
}

// Neg returns the additive inverse -o.
func (o Octonion) Neg() Octonion {
	return o.Scale(-1)
}

// Mul returns the octonion product o*p computed by the Cayley-Dickson formula
// (a,b)(c,d) = (a*c - conj(d)*b, d*a + b*conj(c)) over the quaternions.
// Multiplication is neither commutative nor associative.
func (o Octonion) Mul(p Octonion) Octonion {
	a, b := hypercomplexOctToPair(o)
	c, d := hypercomplexOctToPair(p)
	first := a.Mul(c).Sub(d.Conj().Mul(b))
	second := d.Mul(a).Add(b.Mul(c.Conj()))
	return hypercomplexPairToOct(first, second)
}

// Conj returns the octonion conjugate, which negates all seven imaginary
// coefficients and leaves the real part unchanged.
func (o Octonion) Conj() Octonion {
	return Octonion{o.E0, -o.E1, -o.E2, -o.E3, -o.E4, -o.E5, -o.E6, -o.E7}
}

// Dot returns the Euclidean dot product of o and p treated as 8-vectors.
func (o Octonion) Dot(p Octonion) float64 {
	return o.E0*p.E0 + o.E1*p.E1 + o.E2*p.E2 + o.E3*p.E3 +
		o.E4*p.E4 + o.E5*p.E5 + o.E6*p.E6 + o.E7*p.E7
}

// NormSq returns the squared norm, the sum of the squares of the eight
// coefficients, equal to the real part of o * conj(o).
func (o Octonion) NormSq() float64 {
	return o.Dot(o)
}

// Norm returns the Euclidean norm |o| = sqrt(NormSq).
func (o Octonion) Norm() float64 {
	return math.Sqrt(o.NormSq())
}

// Inverse returns the multiplicative inverse conj(o)/|o|^2. A zero octonion is
// returned unchanged.
func (o Octonion) Inverse() Octonion {
	n2 := o.NormSq()
	if n2 == 0 {
		return o
	}
	return o.Conj().Scale(1 / n2)
}

// Real returns the real (scalar) part E0 of o.
func (o Octonion) Real() float64 {
	return o.E0
}

// IsUnit reports whether o has unit norm to within the absolute tolerance tol.
func (o Octonion) IsUnit(tol float64) bool {
	return math.Abs(o.Norm()-1) <= tol
}

// Equal reports whether o and p agree in every coefficient to within the
// absolute tolerance tol.
func (o Octonion) Equal(p Octonion, tol float64) bool {
	return math.Abs(o.E0-p.E0) <= tol && math.Abs(o.E1-p.E1) <= tol &&
		math.Abs(o.E2-p.E2) <= tol && math.Abs(o.E3-p.E3) <= tol &&
		math.Abs(o.E4-p.E4) <= tol && math.Abs(o.E5-p.E5) <= tol &&
		math.Abs(o.E6-p.E6) <= tol && math.Abs(o.E7-p.E7) <= tol
}

// Commutator returns the commutator o*p - p*o, which vanishes exactly when o
// and p commute. It is generally non-zero because octonion multiplication is
// non-commutative.
func Commutator(o, p Octonion) Octonion {
	return o.Mul(p).Sub(p.Mul(o))
}

// Associator returns the associator (o*p)*q - o*(p*q), a measure of the failure
// of associativity. It vanishes whenever any two of the three arguments are
// equal (the alternative-algebra property) but is generally non-zero otherwise.
func Associator(o, p, q Octonion) Octonion {
	return o.Mul(p).Mul(q).Sub(o.Mul(p.Mul(q)))
}
