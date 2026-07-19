package rootfind

// SignVariations returns the number of sign changes in the sequence of values,
// skipping zeros. It is the basic count underlying Descartes' rule of signs and
// the Budan-Fourier theorem.
func SignVariations(vals []float64) int {
	prev := 0
	count := 0
	for _, v := range vals {
		s := 0
		if v > 0 {
			s = 1
		} else if v < 0 {
			s = -1
		}
		if s == 0 {
			continue
		}
		if prev != 0 && s != prev {
			count++
		}
		prev = s
	}
	return count
}

// CoeffSignVariations returns the number of sign changes among the coefficients
// of p taken in order of ascending power (zeros skipped). By Descartes' rule of
// signs this is an upper bound, of the correct parity, on the number of positive
// real roots of p counted with multiplicity.
func CoeffSignVariations(p Poly) int {
	return SignVariations(p)
}

// DescartesPositiveBound returns the Descartes upper bound on the number of
// positive real roots of p (counted with multiplicity): the sign-variation count
// of its coefficient sequence. The true number of positive roots equals this
// value minus a nonnegative even integer.
func DescartesPositiveBound(p Poly) int {
	return SignVariations(p)
}

// DescartesNegativeBound returns the Descartes upper bound on the number of
// negative real roots of p, obtained by applying the rule of signs to p(-x).
func DescartesNegativeBound(p Poly) int {
	return SignVariations(reflectCoeffs(p))
}

// DescartesRuleOfSigns returns both Descartes bounds at once: an upper bound on
// the number of positive real roots and an upper bound on the number of negative
// real roots, each of the correct parity and counted with multiplicity.
func DescartesRuleOfSigns(p Poly) (positive, negative int) {
	return DescartesPositiveBound(p), DescartesNegativeBound(p)
}

// reflectCoeffs returns the coefficient slice of p(-x): every odd-power
// coefficient is negated.
func reflectCoeffs(p Poly) []float64 {
	out := make([]float64, len(p))
	for i, c := range p {
		if i%2 == 1 {
			out[i] = -c
		} else {
			out[i] = c
		}
	}
	return out
}

// ReflectX returns the polynomial p(-x), whose positive roots are the negatives
// of the negative roots of p.
func (p Poly) ReflectX() Poly {
	return Poly(reflectCoeffs(p))
}

// FourierSequence returns the Fourier (derivative) sequence of p:
// p, p', p”, ... , p^(n), a chain of length deg(p)+1 whose sign variations feed
// the Budan-Fourier theorem.
func FourierSequence(p Poly) []Poly {
	d := p.Degree()
	if d < 0 {
		return nil
	}
	seq := make([]Poly, 0, d+1)
	cur := p.Trim().Clone()
	for cur.Degree() >= 0 {
		seq = append(seq, cur)
		if cur.Degree() == 0 {
			break
		}
		cur = cur.Derivative()
	}
	return seq
}

// fourierVariationsAt evaluates the Fourier sequence at x and returns its sign
// variation count.
func fourierVariationsAt(seq []Poly, x float64) int {
	vals := make([]float64, len(seq))
	for i, q := range seq {
		vals[i] = q.Eval(x)
	}
	return SignVariations(vals)
}

// BudanFourierCount returns V(a) - V(b), where V is the sign-variation count of
// the Fourier sequence of p. By the Budan-Fourier theorem this value is an upper
// bound, differing from the true count of real roots of p in (a, b] (with
// multiplicity) by a nonnegative even integer, and has the same parity as that
// count. The endpoints must satisfy a < b.
func BudanFourierCount(p Poly, a, b float64) int {
	seq := FourierSequence(p)
	if len(seq) == 0 {
		return 0
	}
	return fourierVariationsAt(seq, a) - fourierVariationsAt(seq, b)
}

// BudanFourierUpperBound returns an upper bound on the number of real roots of p
// (with multiplicity) in the open positive axis (0, +inf), via the Fourier
// sequence evaluated at 0 and at a value beyond all roots.
func BudanFourierUpperBound(p Poly) int {
	b := LagrangeBound(p) * 1.001
	if b <= 0 {
		b = 1
	}
	return BudanFourierCount(p, 0, b)
}
