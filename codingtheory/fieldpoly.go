package codingtheory

// This file provides polynomial arithmetic with coefficients in GF(2^m).
// Polynomials are []int slices in which index i holds the coefficient of x^i.
// The empty or all-zero slice represents the zero polynomial. Helper results
// are trimmed of high-order zero coefficients except where noted.

// polyTrim removes trailing (high-degree) zero coefficients, returning at least
// a length-one slice for the zero polynomial.
func polyTrim(p []int) []int {
	n := len(p)
	for n > 1 && p[n-1] == 0 {
		n--
	}
	return p[:n]
}

// PolyDegree returns the degree of a polynomial over GF(2^m). The zero
// polynomial has degree -1.
func (f *Field) PolyDegree(p []int) int {
	for i := len(p) - 1; i >= 0; i-- {
		if p[i] != 0 {
			return i
		}
	}
	return -1
}

// PolyAdd returns the sum of two polynomials over GF(2^m).
func (f *Field) PolyAdd(a, b []int) []int {
	n := len(a)
	if len(b) > n {
		n = len(b)
	}
	out := make([]int, n)
	for i := range a {
		out[i] ^= a[i]
	}
	for i := range b {
		out[i] ^= b[i]
	}
	return polyTrim(out)
}

// PolySub returns the difference of two polynomials; over GF(2^m) this equals
// PolyAdd.
func (f *Field) PolySub(a, b []int) []int { return f.PolyAdd(a, b) }

// PolyScale multiplies every coefficient of p by the scalar s in GF(2^m).
func (f *Field) PolyScale(p []int, s int) []int {
	out := make([]int, len(p))
	for i := range p {
		out[i] = f.Mul(p[i], s)
	}
	return polyTrim(out)
}

// PolyMul returns the product of two polynomials over GF(2^m).
func (f *Field) PolyMul(a, b []int) []int {
	if len(a) == 0 || len(b) == 0 {
		return []int{0}
	}
	out := make([]int, len(a)+len(b)-1)
	for i := range a {
		if a[i] == 0 {
			continue
		}
		la := f.log[a[i]]
		for j := range b {
			if b[j] == 0 {
				continue
			}
			out[i+j] ^= f.exp[la+f.log[b[j]]]
		}
	}
	return polyTrim(out)
}

// PolyEval evaluates the polynomial p at the field element x using Horner's
// rule.
func (f *Field) PolyEval(p []int, x int) int {
	var y int
	for i := len(p) - 1; i >= 0; i-- {
		y = f.Mul(y, x) ^ p[i]
	}
	return y
}

// PolyDivMod divides a by the non-zero polynomial b over GF(2^m), returning the
// quotient and remainder with deg(remainder) < deg(b). It panics if b is zero.
func (f *Field) PolyDivMod(a, b []int) (quotient, remainder []int) {
	b = polyTrim(b)
	if len(b) == 1 && b[0] == 0 {
		panic("codingtheory: polynomial division by zero over GF(2^m)")
	}
	r := append([]int(nil), a...)
	db := f.PolyDegree(b)
	lead := b[db]
	q := make([]int, len(a))
	for {
		dr := f.PolyDegree(r)
		if dr < db {
			break
		}
		coef := f.Div(r[dr], lead)
		shift := dr - db
		q[shift] = coef
		for i := 0; i <= db; i++ {
			r[i+shift] ^= f.Mul(coef, b[i])
		}
	}
	return polyTrim(q), polyTrim(r)
}

// PolyMod returns a modulo the non-zero polynomial b over GF(2^m).
func (f *Field) PolyMod(a, b []int) []int {
	_, r := f.PolyDivMod(a, b)
	return r
}

// PolyDeriv returns the formal derivative of p over GF(2^m). Because the field
// has characteristic two, coefficients of even-power terms vanish.
func (f *Field) PolyDeriv(p []int) []int {
	if len(p) <= 1 {
		return []int{0}
	}
	out := make([]int, len(p)-1)
	for i := 1; i < len(p); i++ {
		if i&1 == 1 {
			out[i-1] = p[i]
		}
	}
	return polyTrim(out)
}

// PolyMonic returns p scaled so its leading coefficient is one. The zero
// polynomial is returned unchanged.
func (f *Field) PolyMonic(p []int) []int {
	d := f.PolyDegree(p)
	if d < 0 {
		return []int{0}
	}
	return f.PolyScale(p, f.Inv(p[d]))
}

// BerlekampMassey returns the error-locator polynomial (little-endian, with
// constant term one) implied by the syndrome sequence s over GF(2^m) using the
// Berlekamp-Massey algorithm. The degree of the result is the number of errors
// the syndromes describe.
func (f *Field) BerlekampMassey(s []int) []int {
	lambda := []int{1}
	b := []int{1}
	L := 0
	m := 1
	bb := 1
	for n := 0; n < len(s); n++ {
		delta := s[n]
		for i := 1; i <= L; i++ {
			if i < len(lambda) {
				delta ^= f.Mul(lambda[i], s[n-i])
			}
		}
		if delta == 0 {
			m++
			continue
		}
		scale := f.Div(delta, bb)
		shifted := polyShiftLE(b, m)
		T := f.PolyAdd(lambda, f.PolyScale(shifted, scale))
		if 2*L <= n {
			b = lambda
			lambda = T
			L = n + 1 - L
			bb = delta
			m = 1
		} else {
			lambda = T
			m++
		}
	}
	return lambda
}

// ChienSearch returns the little-endian error positions (powers of x in [0,n))
// for which the locator polynomial lambda evaluates to zero at alpha^{-i}.
func (f *Field) ChienSearch(lambda []int, n int) []int {
	var positions []int
	for i := 0; i < n; i++ {
		if f.PolyEval(lambda, f.Exp(-i)) == 0 {
			positions = append(positions, i)
		}
	}
	return positions
}

// MinimalPoly returns the minimal polynomial over GF(2) of the element
// alpha^e, packed as a GF(2) polynomial int (bit i is the coefficient of x^i).
// It is the product over the cyclotomic coset of e of (x - alpha^s).
func (f *Field) MinimalPoly(e int) int {
	order := f.Order()
	e = ((e % order) + order) % order
	// cyclotomic coset of e under multiplication by 2 mod order
	coset := map[int]bool{}
	s := e
	for !coset[s] {
		coset[s] = true
		s = (s * 2) % order
	}
	// build product of (x - alpha^s) with coefficients in GF(2^m)
	poly := []int{1}
	for s := range coset {
		poly = f.PolyMul(poly, []int{f.Exp(s), 1})
	}
	// coefficients must lie in GF(2); pack to int
	var out int
	for i, c := range poly {
		if c == 1 {
			out |= 1 << uint(i)
		} else if c != 0 {
			// should not happen for a genuine minimal polynomial
			return 0
		}
	}
	return out
}

// PolyEqual reports whether a and b are the same polynomial, ignoring trailing
// zero coefficients.
func (f *Field) PolyEqual(a, b []int) bool {
	a = polyTrim(append([]int(nil), a...))
	b = polyTrim(append([]int(nil), b...))
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}
