package rootfind

import "math"

// CauchyBound returns Cauchy's bound: every (real or complex) root z of p
// satisfies |z| <= 1 + max_{i<n} |a_i / a_n|, where a_n is the leading
// coefficient. It returns 0 for constant polynomials.
func CauchyBound(p Poly) float64 {
	d := p.Degree()
	if d <= 0 {
		return 0
	}
	lc := p[d]
	m := 0.0
	for i := 0; i < d; i++ {
		if v := math.Abs(p[i] / lc); v > m {
			m = v
		}
	}
	return 1 + m
}

// LagrangeBound returns Lagrange's bound: |z| <= max(1, sum_{i<n} |a_i/a_n|)
// for every root z of p. It is often tighter than [CauchyBound].
func LagrangeBound(p Poly) float64 {
	d := p.Degree()
	if d <= 0 {
		return 0
	}
	lc := p[d]
	s := 0.0
	for i := 0; i < d; i++ {
		s += math.Abs(p[i] / lc)
	}
	return math.Max(1, s)
}

// FujiwaraBound returns Fujiwara's bound on the root moduli:
// |z| <= 2 * max_i ( |a_{n-i}/a_n| ^ (1/i) ), which is one of the sharpest
// simple upper bounds available.
func FujiwaraBound(p Poly) float64 {
	d := p.Degree()
	if d <= 0 {
		return 0
	}
	lc := p[d]
	m := 0.0
	for i := 1; i <= d; i++ {
		c := math.Abs(p[d-i] / lc)
		if i == d {
			c /= 2
		}
		v := math.Pow(c, 1/float64(i))
		if v > m {
			m = v
		}
	}
	return 2 * m
}

// KojimaBound returns Kojima's bound on the root moduli, formed from the ratios
// of consecutive coefficients: |z| <= 2 * max_i q_i, where q_1 = |a_{n-1}/a_n|
// and q_i = |a_{n-i}/a_{n-i+1}| for i>1 (terms with a zero denominator are
// skipped).
func KojimaBound(p Poly) float64 {
	d := p.Degree()
	if d <= 0 {
		return 0
	}
	m := 0.0
	for i := 1; i <= d; i++ {
		num := math.Abs(p[d-i])
		den := math.Abs(p[d-i+1])
		if den == 0 {
			continue
		}
		q := num / den
		if i == d {
			q /= 2
		}
		if q > m {
			m = q
		}
	}
	return 2 * m
}

// LowerRootBound returns a positive lower bound on the moduli of the nonzero
// roots of p: no nonzero root satisfies |z| < LowerRootBound(p). It is obtained
// by applying [CauchyBound] to the reversed polynomial and taking the
// reciprocal. It returns 0 when 0 is a root or p is constant.
func LowerRootBound(p Poly) float64 {
	t := p.Trim()
	if t.Degree() <= 0 {
		return 0
	}
	if t[0] == 0 {
		return 0
	}
	ub := CauchyBound(t.Reverse())
	if ub == 0 {
		return 0
	}
	return 1 / ub
}

// AnnulusBounds returns a lower and an upper bound r0 <= |z| <= r1 that enclose
// every nonzero root z of p in a complex annulus. When 0 is a root the lower
// bound is 0.
func AnnulusBounds(p Poly) (lo, hi float64) {
	return LowerRootBound(p), CauchyBound(p)
}

// RealRootInterval returns a symmetric interval [-b, b] that contains all real
// roots of p, where b is the Lagrange bound. Every real root x satisfies
// -b <= x <= b.
func RealRootInterval(p Poly) (lo, hi float64) {
	b := LagrangeBound(p)
	return -b, b
}
