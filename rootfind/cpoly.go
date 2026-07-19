package rootfind

import (
	"math/cmplx"
)

// CPoly is a complex polynomial stored as coefficients in ascending order of
// power: c[i] is the coefficient of x^i and c[len(c)-1] is the leading
// coefficient.
type CPoly []complex128

// NewCPoly returns a CPoly built from the given ascending-order coefficients.
// The arguments are copied.
func NewCPoly(coeffs ...complex128) CPoly {
	c := make(CPoly, len(coeffs))
	copy(c, coeffs)
	return c
}

// Clone returns an independent copy of c.
func (c CPoly) Clone() CPoly {
	q := make(CPoly, len(c))
	copy(q, c)
	return q
}

// Degree returns the degree of c, the largest index with a nonzero coefficient,
// or -1 for the zero polynomial.
func (c CPoly) Degree() int {
	for i := len(c) - 1; i >= 0; i-- {
		if c[i] != 0 {
			return i
		}
	}
	return -1
}

// IsZero reports whether c is the zero polynomial.
func (c CPoly) IsZero() bool { return c.Degree() < 0 }

// Trim returns c with trailing zero coefficients removed; the underlying array
// is shared with c.
func (c CPoly) Trim() CPoly {
	d := c.Degree()
	if d < 0 {
		return CPoly{}
	}
	return c[:d+1]
}

// Coeff returns the coefficient of x^i, or 0 when i is out of range.
func (c CPoly) Coeff(i int) complex128 {
	if i < 0 || i >= len(c) {
		return 0
	}
	return c[i]
}

// LeadingCoeff returns the leading coefficient, or 0 for the zero polynomial.
func (c CPoly) LeadingCoeff() complex128 {
	d := c.Degree()
	if d < 0 {
		return 0
	}
	return c[d]
}

// Eval evaluates c(x) using Horner's method.
func (c CPoly) Eval(x complex128) complex128 {
	d := c.Degree()
	if d < 0 {
		return 0
	}
	y := c[d]
	for i := d - 1; i >= 0; i-- {
		y = y*x + c[i]
	}
	return y
}

// EvalDeriv evaluates c(x) and c'(x) together in a single Horner sweep.
func (c CPoly) EvalDeriv(x complex128) (val, deriv complex128) {
	d := c.Degree()
	if d < 0 {
		return 0, 0
	}
	val = c[d]
	for i := d - 1; i >= 0; i-- {
		deriv = deriv*x + val
		val = val*x + c[i]
	}
	return val, deriv
}

// EvalDeriv2 evaluates c, c', and c” at x in one Horner sweep.
func (c CPoly) EvalDeriv2(x complex128) (val, d1, d2 complex128) {
	d := c.Degree()
	if d < 0 {
		return 0, 0, 0
	}
	val = c[d]
	for i := d - 1; i >= 0; i-- {
		d2 = d2*x + d1
		d1 = d1*x + val
		val = val*x + c[i]
	}
	d2 *= 2
	return val, d1, d2
}

// Derivative returns the formal derivative c'(x).
func (c CPoly) Derivative() CPoly {
	d := c.Degree()
	if d <= 0 {
		return CPoly{}
	}
	q := make(CPoly, d)
	for i := 1; i <= d; i++ {
		q[i-1] = c[i] * complex(float64(i), 0)
	}
	return q
}

// Add returns the sum c+q.
func (c CPoly) Add(q CPoly) CPoly {
	n := len(c)
	if len(q) > n {
		n = len(q)
	}
	r := make(CPoly, n)
	for i := 0; i < n; i++ {
		r[i] = c.Coeff(i) + q.Coeff(i)
	}
	return r
}

// Sub returns the difference c-q.
func (c CPoly) Sub(q CPoly) CPoly {
	n := len(c)
	if len(q) > n {
		n = len(q)
	}
	r := make(CPoly, n)
	for i := 0; i < n; i++ {
		r[i] = c.Coeff(i) - q.Coeff(i)
	}
	return r
}

// Scale returns the polynomial s*c.
func (c CPoly) Scale(s complex128) CPoly {
	r := make(CPoly, len(c))
	for i, v := range c {
		r[i] = s * v
	}
	return r
}

// Neg returns -c.
func (c CPoly) Neg() CPoly { return c.Scale(-1) }

// Mul returns the product c*q by convolution.
func (c CPoly) Mul(q CPoly) CPoly {
	dc, dq := c.Degree(), q.Degree()
	if dc < 0 || dq < 0 {
		return CPoly{}
	}
	r := make(CPoly, dc+dq+1)
	for i := 0; i <= dc; i++ {
		for j := 0; j <= dq; j++ {
			r[i+j] += c[i] * q[j]
		}
	}
	return r
}

// Monic returns c divided by its leading coefficient. It returns
// ErrZeroPolynomial when c is zero.
func (c CPoly) Monic() (CPoly, error) {
	d := c.Degree()
	if d < 0 {
		return nil, ErrZeroPolynomial
	}
	lc := c[d]
	r := make(CPoly, d+1)
	for i := 0; i <= d; i++ {
		r[i] = c[i] / lc
	}
	return r, nil
}

// Deflate divides c by the linear factor (x - r), returning the quotient and the
// remainder c(r) via synthetic division. When r is a root the remainder is
// (near) zero.
func (c CPoly) Deflate(r complex128) (quo CPoly, remainder complex128) {
	d := c.Degree()
	if d < 0 {
		return CPoly{}, 0
	}
	if d == 0 {
		return CPoly{}, c[0]
	}
	q := make(CPoly, d)
	acc := c[d]
	for i := d - 1; i >= 0; i-- {
		q[i] = acc
		acc = acc*r + c[i]
	}
	return q, acc
}

// CFromRoots builds the monic complex polynomial with the given roots, the
// product (x - r0)(x - r1)... .
func CFromRoots(roots ...complex128) CPoly {
	p := CPoly{1}
	for _, r := range roots {
		p = p.Mul(CPoly{-r, 1})
	}
	return p
}

// ToReal returns the real polynomial formed by taking the real part of each
// coefficient of c, discarding imaginary parts.
func (c CPoly) ToReal() Poly {
	p := make(Poly, len(c))
	for i, v := range c {
		p[i] = real(v)
	}
	return p
}

// cpolyMaxAbs returns the largest coefficient modulus of c.
func cpolyMaxAbs(c CPoly) float64 {
	m := 0.0
	for _, v := range c {
		if a := cmplx.Abs(v); a > m {
			m = a
		}
	}
	return m
}
