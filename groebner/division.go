package groebner

// DivisionResult holds the outcome of multivariate division of a dividend by an
// ordered list of divisors: the quotient polynomials (one per divisor) and the
// remainder.
type DivisionResult struct {
	Quotients []Poly
	Remainder Poly
}

// MultivariateDivide performs the multivariate division algorithm, dividing f
// by the ordered list of divisors with respect to the monomial order o. It
// returns quotients q_i and a remainder r such that
//
//	f = q_1*g_1 + ... + q_s*g_s + r,
//
// where no monomial of r is divisible by the leading monomial of any divisor.
// The result depends on the order of the divisors. Zero divisors are skipped.
func MultivariateDivide(f Poly, divisors []Poly, o Order) DivisionResult {
	n := f.nvars
	q := make([]Poly, len(divisors))
	for i := range q {
		q[i] = Zero(n)
	}
	r := Zero(n)
	p := f.Clone()

	for !p.IsZero() {
		lp := p.LeadingTerm(o)
		divided := false
		for i, g := range divisors {
			if g.IsZero() {
				continue
			}
			lg := g.LeadingTerm(o)
			if quo, ok := lg.Div(lp); ok {
				q[i] = q[i].Add(FromTerm(n, quo))
				p = p.Sub(g.MulTerm(quo))
				divided = true
				break
			}
		}
		if !divided {
			r = r.Add(FromTerm(n, lp))
			p = p.Sub(FromTerm(n, lp))
		}
	}
	return DivisionResult{Quotients: q, Remainder: r}
}

// Remainder returns just the remainder of dividing f by the divisors with
// respect to o. It is the normal form of f modulo the divisor list.
func Remainder(f Poly, divisors []Poly, o Order) Poly {
	return MultivariateDivide(f, divisors, o).Remainder
}

// DivideOne divides f by a single polynomial g and returns the quotient and
// remainder with respect to o.
func DivideOne(f, g Poly, o Order) (quotient, remainder Poly) {
	res := MultivariateDivide(f, []Poly{g}, o)
	return res.Quotients[0], res.Remainder
}

// Divides reports whether g divides f exactly (the remainder of f divided by g
// is zero) with respect to the order o.
func Divides(g, f Poly, o Order) bool {
	_, rem := DivideOne(f, g, o)
	return rem.IsZero()
}

// ExactQuotient returns f/g when g divides f exactly. The boolean result is
// false when the division leaves a nonzero remainder.
func ExactQuotient(f, g Poly, o Order) (Poly, bool) {
	quo, rem := DivideOne(f, g, o)
	if !rem.IsZero() {
		return Zero(f.nvars), false
	}
	return quo, true
}

// IsReducible reports whether some monomial of f is divisible by the leading
// monomial of one of the divisors, i.e. whether a division step is possible.
func IsReducible(f Poly, divisors []Poly, o Order) bool {
	for _, t := range f.terms {
		for _, g := range divisors {
			if g.IsZero() {
				continue
			}
			if g.LeadingMonomial(o).Divides(t.Mono) {
				return true
			}
		}
	}
	return false
}

// SPolynomial returns the S-polynomial of f and g with respect to the order o:
//
//	S(f,g) = (L/LT(f))*f - (L/LT(g))*g,
//
// where L is the least common multiple of the leading monomials of f and g.
// The S-polynomial cancels the leading terms and is the central object of
// Buchberger's algorithm. It returns the zero polynomial if either input is
// zero.
func SPolynomial(f, g Poly, o Order) Poly {
	if f.IsZero() || g.IsZero() {
		return Zero(f.nvars)
	}
	lf := f.LeadingTerm(o)
	lg := g.LeadingTerm(o)
	lcm := lf.Mono.LCM(lg.Mono)

	af, _ := lf.Mono.Div(lcm)
	ag, _ := lg.Mono.Div(lcm)

	tf := Term{Coeff: ratDiv(bigOne, lf.Coeff), Mono: af}
	tg := Term{Coeff: ratDiv(bigOne, lg.Coeff), Mono: ag}

	return f.MulTerm(tf).Sub(g.MulTerm(tg))
}
