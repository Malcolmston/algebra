package controltheory

import (
	"math"
	"math/cmplx"
)

// Poly represents a real polynomial in the Laplace variable s using the
// ascending-power convention: element i is the coefficient of s^i. For example
// Poly{6, 5, 1} is the polynomial 6 + 5s + s^2. The zero-length polynomial and
// a polynomial of all zeros both represent the constant 0.
type Poly []float64

// NewPoly returns a Poly built from the given ascending-power coefficients.
// coeffs[i] becomes the coefficient of s^i. The slice is copied so later
// mutation of the argument does not affect the result.
func NewPoly(coeffs ...float64) Poly {
	p := make(Poly, len(coeffs))
	copy(p, coeffs)
	return p
}

// controltheoryTrim returns the polynomial with trailing (highest-power)
// coefficients that are exactly representable as zero removed, so the length
// reflects the true degree plus one.
func controltheoryTrim(p Poly) Poly {
	n := len(p)
	for n > 0 && p[n-1] == 0 {
		n--
	}
	return p[:n:n]
}

// Degree returns the degree of the polynomial, i.e. the highest power of s
// with a nonzero coefficient. The zero polynomial has degree 0 by convention.
func (p Poly) Degree() int {
	t := controltheoryTrim(p)
	if len(t) == 0 {
		return 0
	}
	return len(t) - 1
}

// LeadingCoeff returns the coefficient of the highest nonzero power of s.
// For the zero polynomial it returns 0.
func (p Poly) LeadingCoeff() float64 {
	t := controltheoryTrim(p)
	if len(t) == 0 {
		return 0
	}
	return t[len(t)-1]
}

// Eval evaluates the polynomial at the real point x using Horner's method.
func (p Poly) Eval(x float64) float64 {
	var acc float64
	for i := len(p) - 1; i >= 0; i-- {
		acc = acc*x + p[i]
	}
	return acc
}

// EvalComplex evaluates the polynomial at the complex point z using Horner's
// method, which is the form needed for frequency-response calculations.
func (p Poly) EvalComplex(z complex128) complex128 {
	var acc complex128
	for i := len(p) - 1; i >= 0; i-- {
		acc = acc*z + complex(p[i], 0)
	}
	return acc
}

// Add returns the sum p + q as a new polynomial.
func (p Poly) Add(q Poly) Poly {
	n := len(p)
	if len(q) > n {
		n = len(q)
	}
	r := make(Poly, n)
	for i := 0; i < n; i++ {
		var a, b float64
		if i < len(p) {
			a = p[i]
		}
		if i < len(q) {
			b = q[i]
		}
		r[i] = a + b
	}
	return controltheoryTrim(r)
}

// Sub returns the difference p - q as a new polynomial.
func (p Poly) Sub(q Poly) Poly {
	return p.Add(q.Scale(-1))
}

// Scale returns the polynomial with every coefficient multiplied by k.
func (p Poly) Scale(k float64) Poly {
	r := make(Poly, len(p))
	for i, c := range p {
		r[i] = c * k
	}
	return r
}

// Mul returns the product p*q. Polynomial multiplication is the convolution of
// the two coefficient sequences.
func (p Poly) Mul(q Poly) Poly {
	if len(p) == 0 || len(q) == 0 {
		return Poly{}
	}
	r := make(Poly, len(p)+len(q)-1)
	for i, a := range p {
		for j, b := range q {
			r[i+j] += a * b
		}
	}
	return controltheoryTrim(r)
}

// Derivative returns the derivative dp/ds as a new polynomial.
func (p Poly) Derivative() Poly {
	if len(p) <= 1 {
		return Poly{}
	}
	r := make(Poly, len(p)-1)
	for i := 1; i < len(p); i++ {
		r[i-1] = p[i] * float64(i)
	}
	return controltheoryTrim(r)
}

// DivMod divides p by q and returns the quotient and remainder polynomials
// such that p = quotient*q + remainder with degree(remainder) < degree(q).
// It panics if q is the zero polynomial.
func (p Poly) DivMod(q Poly) (quotient, remainder Poly) {
	dq := controltheoryTrim(q)
	if len(dq) == 0 {
		panic("controltheory: division by zero polynomial")
	}
	rem := append(Poly{}, controltheoryTrim(p)...)
	dLead := dq[len(dq)-1]
	quo := make(Poly, 0)
	if len(rem) < len(dq) {
		return controltheoryTrim(Poly{0}), controltheoryTrim(rem)
	}
	quo = make(Poly, len(rem)-len(dq)+1)
	for len(rem) >= len(dq) && len(rem) > 0 {
		shift := len(rem) - len(dq)
		factor := rem[len(rem)-1] / dLead
		quo[shift] = factor
		for i := 0; i < len(dq); i++ {
			rem[shift+i] -= factor * dq[i]
		}
		rem = controltheoryTrim(rem)
	}
	return controltheoryTrim(quo), controltheoryTrim(rem)
}

// PolyFromRoots returns the monic real polynomial whose roots are the given
// complex values. Roots are expected to appear in conjugate pairs so the
// resulting coefficients are real; any residual imaginary parts from rounding
// are discarded.
func PolyFromRoots(roots ...complex128) Poly {
	acc := []complex128{complex(1, 0)}
	for _, r := range roots {
		next := make([]complex128, len(acc)+1)
		for i, c := range acc {
			next[i] += c * (-r)
			next[i+1] += c
		}
		acc = next
	}
	p := make(Poly, len(acc))
	for i, c := range acc {
		p[i] = real(c)
	}
	return controltheoryTrim(p)
}

// Roots returns all complex roots of the polynomial using the Durand-Kerner
// (Weierstrass) iteration with deterministic starting points. Real roots are
// returned with a zero (or negligible) imaginary part. The zero polynomial and
// nonzero constants have no roots and return an empty slice.
func (p Poly) Roots() []complex128 {
	t := controltheoryTrim(p)
	n := len(t) - 1
	if n <= 0 {
		return nil
	}
	// Monic complex coefficients, ascending.
	lead := t[n]
	mon := make([]complex128, n+1)
	for i := 0; i <= n; i++ {
		mon[i] = complex(t[i]/lead, 0)
	}
	if n == 1 {
		return []complex128{-mon[0]}
	}
	evalMonic := func(z complex128) complex128 {
		var acc complex128
		for i := n; i >= 0; i-- {
			acc = acc*z + mon[i]
		}
		return acc
	}
	// Deterministic initial guesses spread around the complex plane.
	z := make([]complex128, n)
	seed := complex(0.4, 0.9)
	cur := complex(1, 0)
	for i := 0; i < n; i++ {
		cur *= seed
		z[i] = cur
	}
	const maxIter = 1000
	const tol = 1e-14
	for iter := 0; iter < maxIter; iter++ {
		maxDelta := 0.0
		for i := 0; i < n; i++ {
			denom := complex(1, 0)
			for j := 0; j < n; j++ {
				if j != i {
					denom *= z[i] - z[j]
				}
			}
			if denom == 0 {
				denom = complex(1e-30, 0)
			}
			delta := evalMonic(z[i]) / denom
			z[i] -= delta
			if d := cmplx.Abs(delta); d > maxDelta {
				maxDelta = d
			}
		}
		if maxDelta < tol {
			break
		}
	}
	// Clean negligible imaginary parts.
	for i := range z {
		if math.Abs(imag(z[i])) < 1e-10*(1+math.Abs(real(z[i]))) {
			z[i] = complex(real(z[i]), 0)
		}
	}
	return z
}
