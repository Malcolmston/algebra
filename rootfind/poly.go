package rootfind

import (
	"math"
	"math/cmplx"
	"strconv"
	"strings"
)

// Poly is a real polynomial stored as coefficients in ascending order of power:
// p[i] is the coefficient of x^i, so p[len(p)-1] is the leading coefficient.
// The zero polynomial is represented by an empty or all-zero slice.
type Poly []float64

// NewPoly returns a Poly built from the given ascending-order coefficients. The
// arguments are copied, so the caller may reuse the slice afterwards.
func NewPoly(coeffs ...float64) Poly {
	p := make(Poly, len(coeffs))
	copy(p, coeffs)
	return p
}

// PolyFromDesc returns a Poly from coefficients given in descending order of
// power, i.e. highest-degree coefficient first, as polynomials are usually
// written. It is the inverse of [Poly.CoeffsDesc].
func PolyFromDesc(coeffs ...float64) Poly {
	p := make(Poly, len(coeffs))
	for i, c := range coeffs {
		p[len(coeffs)-1-i] = c
	}
	return p
}

// Clone returns an independent copy of p.
func (p Poly) Clone() Poly {
	q := make(Poly, len(p))
	copy(q, p)
	return q
}

// Degree returns the degree of p, that is the largest index whose coefficient
// is nonzero. The zero polynomial has degree -1 by convention.
func (p Poly) Degree() int {
	for i := len(p) - 1; i >= 0; i-- {
		if p[i] != 0 {
			return i
		}
	}
	return -1
}

// IsZero reports whether p is the zero polynomial, i.e. every coefficient is 0.
func (p Poly) IsZero() bool {
	return p.Degree() < 0
}

// Trim returns p with trailing (high-order) zero coefficients removed. The
// result always has length equal to Degree()+1, or length zero for the zero
// polynomial. The underlying array is shared with p.
func (p Poly) Trim() Poly {
	d := p.Degree()
	if d < 0 {
		return Poly{}
	}
	return p[:d+1]
}

// Coeff returns the coefficient of x^i, or 0 when i is out of range.
func (p Poly) Coeff(i int) float64 {
	if i < 0 || i >= len(p) {
		return 0
	}
	return p[i]
}

// LeadingCoeff returns the coefficient of the highest-degree term, or 0 for the
// zero polynomial.
func (p Poly) LeadingCoeff() float64 {
	d := p.Degree()
	if d < 0 {
		return 0
	}
	return p[d]
}

// CoeffsDesc returns the coefficients of p in descending order of power,
// highest-degree first, trimmed of leading zeros. The zero polynomial yields a
// single 0.
func (p Poly) CoeffsDesc() []float64 {
	d := p.Degree()
	if d < 0 {
		return []float64{0}
	}
	out := make([]float64, d+1)
	for i := 0; i <= d; i++ {
		out[i] = p[d-i]
	}
	return out
}

// Eval evaluates p(x) using Horner's method, which is both fast and numerically
// stable.
func (p Poly) Eval(x float64) float64 {
	d := p.Degree()
	if d < 0 {
		return 0
	}
	y := p[d]
	for i := d - 1; i >= 0; i-- {
		y = y*x + p[i]
	}
	return y
}

// EvalComplex evaluates p at a complex argument using Horner's method.
func (p Poly) EvalComplex(x complex128) complex128 {
	d := p.Degree()
	if d < 0 {
		return 0
	}
	y := complex(p[d], 0)
	for i := d - 1; i >= 0; i-- {
		y = y*x + complex(p[i], 0)
	}
	return y
}

// EvalDeriv evaluates p(x) and p'(x) simultaneously with a single Horner sweep,
// returning the value and the first derivative.
func (p Poly) EvalDeriv(x float64) (val, deriv float64) {
	d := p.Degree()
	if d < 0 {
		return 0, 0
	}
	val = p[d]
	deriv = 0
	for i := d - 1; i >= 0; i-- {
		deriv = deriv*x + val
		val = val*x + p[i]
	}
	return val, deriv
}

// EvalDeriv2 evaluates p, p', and p” at x in a single Horner sweep.
func (p Poly) EvalDeriv2(x float64) (val, d1, d2 float64) {
	d := p.Degree()
	if d < 0 {
		return 0, 0, 0
	}
	val = p[d]
	for i := d - 1; i >= 0; i-- {
		d2 = d2*x + d1
		d1 = d1*x + val
		val = val*x + p[i]
	}
	d2 *= 2
	return val, d1, d2
}

// Derivative returns the formal derivative p'(x).
func (p Poly) Derivative() Poly {
	d := p.Degree()
	if d <= 0 {
		return Poly{}
	}
	q := make(Poly, d)
	for i := 1; i <= d; i++ {
		q[i-1] = p[i] * float64(i)
	}
	return q
}

// Integral returns an antiderivative of p whose constant term equals c.
func (p Poly) Integral(c float64) Poly {
	d := p.Degree()
	if d < 0 {
		return Poly{c}
	}
	q := make(Poly, d+2)
	q[0] = c
	for i := 0; i <= d; i++ {
		q[i+1] = p[i] / float64(i+1)
	}
	return q
}

// Add returns the sum p+q.
func (p Poly) Add(q Poly) Poly {
	n := len(p)
	if len(q) > n {
		n = len(q)
	}
	r := make(Poly, n)
	for i := 0; i < n; i++ {
		r[i] = p.Coeff(i) + q.Coeff(i)
	}
	return r.Trim().Clone()
}

// Sub returns the difference p-q.
func (p Poly) Sub(q Poly) Poly {
	n := len(p)
	if len(q) > n {
		n = len(q)
	}
	r := make(Poly, n)
	for i := 0; i < n; i++ {
		r[i] = p.Coeff(i) - q.Coeff(i)
	}
	return r.Trim().Clone()
}

// Neg returns -p.
func (p Poly) Neg() Poly {
	r := make(Poly, len(p))
	for i, c := range p {
		r[i] = -c
	}
	return r
}

// Scale returns the polynomial s*p.
func (p Poly) Scale(s float64) Poly {
	r := make(Poly, len(p))
	for i, c := range p {
		r[i] = s * c
	}
	return r
}

// Mul returns the product p*q computed by convolution of the coefficient
// sequences.
func (p Poly) Mul(q Poly) Poly {
	dp, dq := p.Degree(), q.Degree()
	if dp < 0 || dq < 0 {
		return Poly{}
	}
	r := make(Poly, dp+dq+1)
	for i := 0; i <= dp; i++ {
		if p[i] == 0 {
			continue
		}
		for j := 0; j <= dq; j++ {
			r[i+j] += p[i] * q[j]
		}
	}
	return r
}

// Monic returns p divided by its leading coefficient, so the result has leading
// coefficient 1. It returns ErrZeroPolynomial when p is zero.
func (p Poly) Monic() (Poly, error) {
	d := p.Degree()
	if d < 0 {
		return nil, ErrZeroPolynomial
	}
	lc := p[d]
	r := make(Poly, d+1)
	for i := 0; i <= d; i++ {
		r[i] = p[i] / lc
	}
	return r, nil
}

// Equal reports whether p and q are equal as polynomials, that is every
// coefficient agrees within tol after ignoring trailing zeros.
func (p Poly) Equal(q Poly, tol float64) bool {
	dp, dq := p.Degree(), q.Degree()
	if dp != dq {
		return false
	}
	for i := 0; i <= dp; i++ {
		if math.Abs(p[i]-q[i]) > tol {
			return false
		}
	}
	return true
}

// DivMod divides p by d and returns the quotient q and remainder r satisfying
// p = q*d + r with deg(r) < deg(d). It returns ErrZeroPolynomial when d is zero.
func (p Poly) DivMod(d Poly) (q, r Poly, err error) {
	dd := d.Degree()
	if dd < 0 {
		return nil, nil, ErrZeroPolynomial
	}
	dp := p.Degree()
	if dp < dd {
		return Poly{}, p.Trim().Clone(), nil
	}
	rem := p.Trim().Clone()
	quo := make(Poly, dp-dd+1)
	lc := d[dd]
	for deg := dp; deg >= dd; deg-- {
		if len(rem) <= deg || rem[deg] == 0 {
			continue
		}
		factor := rem[deg] / lc
		quo[deg-dd] = factor
		for i := 0; i <= dd; i++ {
			rem[deg-dd+i] -= factor * d[i]
		}
	}
	return quo.Trim().Clone(), rem.Trim().Clone(), nil
}

// Quo returns the quotient of p divided by d, discarding the remainder.
func (p Poly) Quo(d Poly) (Poly, error) {
	q, _, err := p.DivMod(d)
	return q, err
}

// Rem returns the remainder of p divided by d.
func (p Poly) Rem(d Poly) (Poly, error) {
	_, r, err := p.DivMod(d)
	return r, err
}

// GCD returns a greatest common divisor of p and q using the Euclidean
// algorithm, normalized to be monic. The GCD of two zero polynomials is zero.
func (p Poly) GCD(q Poly) Poly {
	a := p.Trim().Clone()
	b := q.Trim().Clone()
	if a.Degree() < 0 && b.Degree() < 0 {
		return Poly{}
	}
	const relTol = 1e-9
	for b.Degree() >= 0 {
		_, r, err := a.DivMod(b)
		if err != nil {
			break
		}
		// Recognize a remainder that is negligible relative to the current
		// divisor as an exact zero. This makes the floating-point Euclidean
		// algorithm robust to the rounding noise that would otherwise mask a
		// genuine common factor (or fabricate a spurious constant one).
		ref := polyMaxAbs(b) + polyMaxAbs(r)
		r = polyChop(r, relTol*ref)
		a, b = b, r
	}
	if a.Degree() < 0 {
		return Poly{}
	}
	m, _ := a.Monic()
	return m
}

// Compose returns the polynomial p(q(x)) formed by substituting q into p, using
// Horner's method over polynomial arithmetic.
func (p Poly) Compose(q Poly) Poly {
	d := p.Degree()
	if d < 0 {
		return Poly{}
	}
	r := Poly{p[d]}
	for i := d - 1; i >= 0; i-- {
		r = r.Mul(q).Add(Poly{p[i]})
	}
	return r
}

// Reverse returns the reversal (reciprocal) polynomial x^n * p(1/x), whose
// coefficient sequence is that of p reversed. Its nonzero roots are the
// reciprocals of the nonzero roots of p.
func (p Poly) Reverse() Poly {
	d := p.Degree()
	if d < 0 {
		return Poly{}
	}
	r := make(Poly, d+1)
	for i := 0; i <= d; i++ {
		r[i] = p[d-i]
	}
	return r
}

// ShiftScale returns the polynomial p(a*x + b), the composition of p with the
// affine map x |-> a*x + b.
func (p Poly) ShiftScale(a, b float64) Poly {
	return p.Compose(Poly{b, a})
}

// DeflateReal divides p by the linear factor (x - r), returning the quotient and
// the remainder p(r). When r is an exact root the remainder is (near) zero. The
// division is performed by synthetic (Horner) division.
func (p Poly) DeflateReal(r float64) (quo Poly, remainder float64) {
	d := p.Degree()
	if d < 0 {
		return Poly{}, 0
	}
	if d == 0 {
		return Poly{}, p[0]
	}
	q := make(Poly, d)
	acc := p[d]
	for i := d - 1; i >= 0; i-- {
		q[i] = acc
		acc = acc*r + p[i]
	}
	return q, acc
}

// FromRoots builds the monic real polynomial whose roots are exactly the given
// values, i.e. the product (x - r0)(x - r1)... .
func FromRoots(roots ...float64) Poly {
	p := Poly{1}
	for _, r := range roots {
		p = p.Mul(Poly{-r, 1})
	}
	return p
}

// FromRootsWithLead builds the real polynomial lead*(x-r0)(x-r1)... with the
// given leading coefficient.
func FromRootsWithLead(lead float64, roots ...float64) Poly {
	return FromRoots(roots...).Scale(lead)
}

// String renders p in conventional descending-power notation, for example
// "2x^3 - x + 5". The zero polynomial renders as "0".
func (p Poly) String() string {
	d := p.Degree()
	if d < 0 {
		return "0"
	}
	var b strings.Builder
	first := true
	for i := d; i >= 0; i-- {
		c := p[i]
		if c == 0 {
			continue
		}
		mag := math.Abs(c)
		switch {
		case first:
			if c < 0 {
				b.WriteByte('-')
			}
		case c < 0:
			b.WriteString(" - ")
		default:
			b.WriteString(" + ")
		}
		first = false
		if mag != 1 || i == 0 {
			b.WriteString(strconv.FormatFloat(mag, 'g', -1, 64))
		}
		switch i {
		case 0:
		case 1:
			b.WriteString("x")
		default:
			b.WriteString("x^")
			b.WriteString(strconv.Itoa(i))
		}
	}
	return b.String()
}

// polyMaxAbs returns the largest absolute coefficient value of p.
func polyMaxAbs(p Poly) float64 {
	m := 0.0
	for _, c := range p {
		if a := math.Abs(c); a > m {
			m = a
		}
	}
	return m
}

// polyChop zeros out coefficients of p whose magnitude is below tol and returns
// the trimmed result. It is used to control coefficient growth of rounding
// noise in the Euclidean algorithm.
func polyChop(p Poly, tol float64) Poly {
	q := p.Clone()
	for i := range q {
		if math.Abs(q[i]) < tol {
			q[i] = 0
		}
	}
	return q.Trim().Clone()
}

// ToComplex converts a real polynomial to its complex counterpart.
func (p Poly) ToComplex() CPoly {
	c := make(CPoly, len(p))
	for i, v := range p {
		c[i] = complex(v, 0)
	}
	return c
}

// polyIsFinite reports whether every coefficient of p is a finite number.
func polyIsFinite(p Poly) bool {
	for _, c := range p {
		if math.IsNaN(c) || math.IsInf(c, 0) {
			return false
		}
	}
	return true
}

// cabs is a small helper returning the modulus of a complex number.
func cabs(z complex128) float64 { return cmplx.Abs(z) }
