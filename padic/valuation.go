package padic

import "math/big"

// ValuationInt returns the p-adic valuation of the integer n: the exponent of
// the highest power of p dividing n. By convention it returns -1 when n is
// zero (whose valuation is +infinity); callers that care must check n first.
func ValuationInt(p, n *big.Int) int {
	if n.Sign() == 0 {
		return -1
	}
	v := 0
	t := new(big.Int).Abs(n)
	q := new(big.Int)
	r := new(big.Int)
	for {
		q.QuoRem(t, p, r)
		if r.Sign() != 0 {
			break
		}
		t.Set(q)
		v++
	}
	return v
}

// UnitPartInt returns n with all factors of p removed, i.e. n / p^ValuationInt.
// For n = 0 it returns 0.
func UnitPartInt(p, n *big.Int) *big.Int {
	if n.Sign() == 0 {
		return big.NewInt(0)
	}
	neg := n.Sign() < 0
	t := new(big.Int).Abs(n)
	q := new(big.Int)
	r := new(big.Int)
	for {
		q.QuoRem(t, p, r)
		if r.Sign() != 0 {
			break
		}
		t.Set(q)
	}
	if neg {
		t.Neg(t)
	}
	return t
}

// SplitInt returns the valuation v and unit u such that n = p^v * u with u
// coprime to p. For n = 0 it returns (-1, 0).
func SplitInt(p, n *big.Int) (int, *big.Int) {
	return ValuationInt(p, n), UnitPartInt(p, n)
}

// ValuationRat returns the p-adic valuation of the rational a/b, defined as
// ValuationInt(a) - ValuationInt(b). The rational is given by numerator a and
// denominator b, which must be non-zero. If a is zero it returns -1.
func ValuationRat(p, a, b *big.Int) int {
	if a.Sign() == 0 {
		return -1
	}
	return ValuationInt(p, a) - ValuationInt(p, b)
}

// ValuationBigRat returns the p-adic valuation of a big.Rat. A zero rational
// yields -1 by convention.
func ValuationBigRat(p *big.Int, r *big.Rat) int {
	if r.Sign() == 0 {
		return -1
	}
	return ValuationInt(p, r.Num()) - ValuationInt(p, r.Denom())
}

// AbsValueRat returns the p-adic absolute value |a/b|_p = p^(-v) as an exact
// big.Rat, where v is the valuation of the rational a/b. A zero numerator
// yields 0.
func AbsValueRat(p, a, b *big.Int) *big.Rat {
	if a.Sign() == 0 {
		return new(big.Rat)
	}
	v := ValuationRat(p, a, b)
	return absFromVal(p, v)
}

// AbsValueBigRat returns the p-adic absolute value of a big.Rat as an exact
// big.Rat.
func AbsValueBigRat(p *big.Int, r *big.Rat) *big.Rat {
	if r.Sign() == 0 {
		return new(big.Rat)
	}
	return absFromVal(p, ValuationBigRat(p, r))
}

// absFromVal returns p^(-v) as a big.Rat.
func absFromVal(p *big.Int, v int) *big.Rat {
	if v >= 0 {
		return new(big.Rat).SetFrac(bigOne, PPow(p, v))
	}
	return new(big.Rat).SetInt(PPow(p, -v))
}

// AbsValueFloat returns the p-adic absolute value p^(-v) of the rational a/b as
// a float64.
func AbsValueFloat(p, a, b *big.Int) float64 {
	r := AbsValueRat(p, a, b)
	f, _ := r.Float64()
	return f
}

// SameValuation reports whether two integers have the same p-adic valuation.
func SameValuation(p, m, n *big.Int) bool {
	return ValuationInt(p, m) == ValuationInt(p, n)
}
