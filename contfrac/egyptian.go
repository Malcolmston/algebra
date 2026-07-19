package contfrac

import "math/big"

// EgyptianFraction returns the denominators of a sum of distinct unit fractions
// equal to the proper fraction p/q (with 0 < p/q < 1) using Fibonacci's greedy
// algorithm. Each returned denominator d contributes 1/d, and the denominators
// are strictly increasing. It returns an error if p/q is not a proper fraction.
func EgyptianFraction(p, q int64) ([]int64, error) {
	p, q = ReduceFraction(p, q)
	if p <= 0 || p >= q {
		return nil, ErrNotProper
	}
	var out []int64
	for p != 0 {
		d := ceilDiv(q, p) // smallest unit fraction not exceeding p/q
		out = append(out, d)
		// p/q - 1/d = (p*d - q)/(q*d)
		p, q = ReduceFraction(p*d-q, q*d)
	}
	return out, nil
}

// EgyptianFractionFrac is the [Frac] form of [EgyptianFraction].
func EgyptianFractionFrac(f Frac) ([]int64, error) {
	return EgyptianFraction(f.Num, f.Den)
}

// EgyptianFractionBig returns the denominators of a distinct-unit-fraction
// decomposition of the proper fraction r using the greedy algorithm, with
// arbitrary-precision denominators so it never overflows. It returns an error
// if r is not in the open interval (0, 1).
func EgyptianFractionBig(r *big.Rat) ([]*big.Int, error) {
	num := new(big.Int).Set(r.Num())
	den := new(big.Int).Set(r.Denom())
	zero := big.NewInt(0)
	if num.Sign() <= 0 || num.Cmp(den) >= 0 {
		return nil, ErrNotProper
	}
	var out []*big.Int
	one := big.NewInt(1)
	for num.Sign() != 0 {
		// d = ceil(den/num)
		d := new(big.Int)
		m := new(big.Int)
		d.DivMod(den, num, m)
		if m.Sign() != 0 {
			d.Add(d, one)
		}
		out = append(out, new(big.Int).Set(d))
		// num/den - 1/d = (num*d - den)/(den*d)
		newNum := new(big.Int).Sub(new(big.Int).Mul(num, d), den)
		newDen := new(big.Int).Mul(den, d)
		g := new(big.Int).GCD(nil, nil, newNum, newDen)
		if g.Sign() != 0 {
			newNum.Div(newNum, g)
			newDen.Div(newDen, g)
		}
		num, den = newNum, newDen
	}
	_ = zero
	return out, nil
}

// SumUnitFractions returns the exact sum of the unit fractions 1/d for each d in
// dens as a *big.Rat. A denominator of zero is skipped.
func SumUnitFractions(dens []int64) *big.Rat {
	sum := new(big.Rat)
	for _, d := range dens {
		if d == 0 {
			continue
		}
		sum.Add(sum, big.NewRat(1, d))
	}
	return sum
}

// IsUnitFraction reports whether p/q reduces to 1/n for some positive integer n.
func IsUnitFraction(p, q int64) bool {
	p, q = ReduceFraction(p, q)
	return p == 1 && q >= 1
}

// UnitFraction returns the unit fraction 1/d as a [Frac]. It panics if d == 0.
func UnitFraction(d int64) Frac {
	return NewFrac(1, d)
}

// ceilDiv returns ceil(a/b) for a >= 0 and b > 0.
func ceilDiv(a, b int64) int64 {
	return (a + b - 1) / b
}
