package padic

import "math/big"

// StrassmannBound returns the number of zeros, counted with multiplicity, that
// a convergent power series sum a_n x^n has in the closed unit ball of Q_p (the
// ring Z_p). By Strassmann's theorem this equals the largest index N with
// val(a_N) minimal among all coefficients (equivalently |a_N| maximal). The
// coefficients are given as big.Rat values a_0, a_1, ...; they must not all be
// zero. For a finite slice the implicit tail is zero, so convergence holds.
func StrassmannBound(p *big.Int, coeffs []*big.Rat) (int, error) {
	minVal := 0
	found := false
	best := -1
	for n, c := range coeffs {
		if c.Sign() == 0 {
			continue
		}
		v := ValuationBigRat(p, c)
		if !found || v < minVal {
			minVal = v
			best = n
			found = true
		} else if v == minVal {
			best = n
		}
	}
	if !found {
		return 0, ErrDomain
	}
	return best, nil
}

// StrassmannBoundInts is StrassmannBound for integer coefficients.
func StrassmannBoundInts(p *big.Int, coeffs []*big.Int) (int, error) {
	rats := make([]*big.Rat, len(coeffs))
	for i, c := range coeffs {
		rats[i] = new(big.Rat).SetInt(c)
	}
	return StrassmannBound(p, rats)
}

// SeriesValuationAt returns a lower bound on the p-adic valuation of the value
// f(x) = sum a_n x^n, namely min_n (val(a_n) + n*xval), where xval is the
// valuation of the argument x. It returns the minimising bound and whether any
// non-zero coefficient exists.
func SeriesValuationAt(p *big.Int, coeffs []*big.Rat, xval int) (int, bool) {
	best := 0
	found := false
	for n, c := range coeffs {
		if c.Sign() == 0 {
			continue
		}
		v := ValuationBigRat(p, c) + n*xval
		if !found || v < best {
			best = v
			found = true
		}
	}
	return best, found
}

// ConvergesInUnitBall reports whether the power series with the given
// coefficients converges everywhere on the closed unit ball, i.e. the
// coefficient valuations tend to +infinity. For a finite slice this is always
// true (the tail is zero) provided the slice is not empty.
func ConvergesInUnitBall(p *big.Int, coeffs []*big.Rat) bool {
	return len(coeffs) > 0
}

// ConvergesAt reports whether sum a_n x^n converges for an argument of
// valuation xval, i.e. val(a_n) + n*xval -> +infinity. For a finite slice this
// holds whenever every retained term is finite; the practical test is that the
// last non-zero coefficient's contribution grows, which for xval > radius-limit
// is guaranteed. Here, for a finite polynomial, it always converges.
func ConvergesAt(p *big.Int, coeffs []*big.Rat, xval int) bool {
	return len(coeffs) > 0
}

// ConvergenceRadiusSlope returns the "slope" that bounds the region of
// convergence of the power series: the series sum a_n x^n converges for
// arguments x with val(x) > -s, where s = limsup val(a_n)/n. For a finite,
// non-zero slice it returns the maximum of -val(a_n)/n over the non-zero terms
// (n >= 1) as an exact big.Rat, plus whether such a term exists.
func ConvergenceRadiusSlope(p *big.Int, coeffs []*big.Rat) (*big.Rat, bool) {
	var best *big.Rat
	for n := 1; n < len(coeffs); n++ {
		if coeffs[n].Sign() == 0 {
			continue
		}
		v := ValuationBigRat(p, coeffs[n])
		s := new(big.Rat).SetFrac(big.NewInt(int64(-v)), big.NewInt(int64(n)))
		if best == nil || s.Cmp(best) > 0 {
			best = s
		}
	}
	return best, best != nil
}

// StrassmannBoundPadics is StrassmannBound for p-adic coefficients, using their
// exact valuations. All coefficients must share the prime.
func StrassmannBoundPadics(coeffs []*Padic) (int, error) {
	if len(coeffs) == 0 {
		return 0, ErrDomain
	}
	p := coeffs[0].p
	minVal := 0
	found := false
	best := -1
	for n, c := range coeffs {
		if c.p.Cmp(p) != 0 {
			return 0, ErrPrimeMismatch
		}
		if c.IsZero() {
			continue
		}
		v := c.Valuation()
		if !found || v < minVal {
			minVal = v
			best = n
			found = true
		} else if v == minVal {
			best = n
		}
	}
	if !found {
		return 0, ErrDomain
	}
	return best, nil
}

// WeierstrassDegree returns the Weierstrass degree of a convergent power
// series: the number of zeros in the open unit ball, equal to the largest index
// achieving the minimum coefficient valuation among the terms of positive
// index and index 0. It coincides with the Strassmann bound here and is
// provided as a domain-specific synonym.
func WeierstrassDegree(p *big.Int, coeffs []*big.Rat) (int, error) {
	return StrassmannBound(p, coeffs)
}
