package ntheory

import "math/big"

// ntheoryFloorDivMod returns the floored quotient a and non-negative remainder r
// of p divided by q, where q must be positive. It satisfies p == a*q + r with
// 0 <= r < q, matching the flooring convention used by the Euclidean continued
// fraction expansion (Go's built-in / truncates toward zero, which differs for
// negative p).
func ntheoryFloorDivMod(p, q int64) (a, r int64) {
	a = p / q
	r = p - a*q
	if r < 0 {
		a--
		r += q
	}
	return a, r
}

// ContinuedFraction returns the finite continued-fraction expansion
// [a0; a1, a2, ...] of the rational p/q (q != 0) as a slice of int64 partial
// quotients via the Euclidean algorithm. The expansion is the canonical one:
// a0 = floor(p/q) may be negative or zero, while every later partial quotient is
// positive. It panics if q == 0.
func ContinuedFraction(p, q int64) []int64 {
	if q == 0 {
		panic("ntheory: ContinuedFraction requires q != 0")
	}
	if q < 0 {
		p, q = -p, -q
	}
	var cf []int64
	for q != 0 {
		a, r := ntheoryFloorDivMod(p, q)
		cf = append(cf, a)
		p, q = q, r
	}
	return cf
}

// ContinuedFractionRat is the big.Rat form of ContinuedFraction for exact large
// rationals. It returns the canonical finite expansion [a0; a1, a2, ...] of r as
// a slice of int64 partial quotients, computed with the Euclidean algorithm on
// arbitrary-precision integers. r must be non-nil.
func ContinuedFractionRat(r *big.Rat) []int64 {
	// big.Rat keeps the denominator positive, so Euclidean (floored) division of
	// the numerator by the denominator yields the canonical expansion.
	num := new(big.Int).Set(r.Num())
	den := new(big.Int).Set(r.Denom())
	var cf []int64
	q := new(big.Int)
	m := new(big.Int)
	for den.Sign() != 0 {
		q.DivMod(num, den, m) // Euclidean: 0 <= m < den, den > 0.
		cf = append(cf, q.Int64())
		num.Set(den)
		den.Set(m)
	}
	return cf
}

// Convergents returns the successive convergents of the partial-quotient slice
// cf as reduced big.Rat values. The i-th element is the value of the truncated
// continued fraction [cf[0]; cf[1], ..., cf[i]]. They are computed with the
// classic h/k recurrence
//
//	h_i = a_i*h_{i-1} + h_{i-2},  k_i = a_i*k_{i-1} + k_{i-2},
//
// seeded by h_{-1}=1, h_{-2}=0, k_{-1}=0, k_{-2}=1, so each convergent is
// obtained from the previous two without re-multiplying whole fractions. The
// result has len(cf) entries; an empty cf yields an empty slice.
func Convergents(cf []int64) []*big.Rat {
	out := make([]*big.Rat, 0, len(cf))
	hPrev, hPrev2 := big.NewInt(1), big.NewInt(0)
	kPrev, kPrev2 := big.NewInt(0), big.NewInt(1)
	a := new(big.Int)
	for _, ai := range cf {
		a.SetInt64(ai)
		h := new(big.Int).Mul(a, hPrev)
		h.Add(h, hPrev2)
		k := new(big.Int).Mul(a, kPrev)
		k.Add(k, kPrev2)
		// big.Rat.SetFrac reduces the fraction to lowest terms.
		out = append(out, new(big.Rat).SetFrac(h, k))
		hPrev2, hPrev = hPrev, h
		kPrev2, kPrev = kPrev, k
	}
	return out
}

// RatFromContinuedFraction reconstructs the exact big.Rat value of the
// partial-quotient slice cf, i.e. the value of [cf[0]; cf[1], ..., cf[n-1]]. It
// evaluates the expansion from the innermost quotient outward using exact
// rational arithmetic. An empty cf yields 0.
func RatFromContinuedFraction(cf []int64) *big.Rat {
	if len(cf) == 0 {
		return new(big.Rat)
	}
	// Fold from the back: val = a_last, then val = a_i + 1/val.
	val := new(big.Rat).SetInt64(cf[len(cf)-1])
	for i := len(cf) - 2; i >= 0; i-- {
		val.Inv(val)
		val.Add(val, new(big.Rat).SetInt64(cf[i]))
	}
	return val
}

// SqrtContinuedFraction returns the continued-fraction expansion of the square
// root of n for a positive n. The result is the pair (a0, period), where
// a0 = floor(sqrt(n)) is the leading term and period is the repeating block, so
// that sqrt(n) = [a0; period[0], period[1], ..., a0+period repeated]. It uses
// the standard integer (m, d, a) recurrence
//
//	m_{i+1} = d_i*a_i - m_i,
//	d_{i+1} = (n - m_{i+1}^2) / d_i,
//	a_{i+1} = floor((a0 + m_{i+1}) / d_{i+1}),
//
// which is exact with no floating-point rounding; the period terminates when a
// term equal to 2*a0 is reached. If n is a perfect square, sqrt(n) is the
// integer a0 and period is nil.
func SqrtContinuedFraction(n uint64) (a0 int64, period []int64) {
	N := new(big.Int).SetUint64(n)
	a0big := IsqrtBig(N)
	// Perfect square: a0big^2 == n means the expansion is just the integer root.
	if new(big.Int).Mul(a0big, a0big).Cmp(N) == 0 {
		return a0big.Int64(), nil
	}
	a0 = a0big.Int64()

	// Integer recurrence in big.Int to stay exact even for n near 2^64, where
	// m*m would overflow a fixed-width integer.
	m := big.NewInt(0)
	d := big.NewInt(1)
	a := new(big.Int).Set(a0big)
	twoA0 := int64(2) * a0

	tmp := new(big.Int)
	for {
		// m = d*a - m
		tmp.Mul(d, a)
		m.Sub(tmp, m)
		// d = (N - m*m) / d
		tmp.Mul(m, m)
		tmp.Sub(N, tmp)
		d.Quo(tmp, d) // exact division; remainder is always zero here.
		// a = floor((a0 + m) / d)
		tmp.Add(a0big, m)
		a.Div(tmp, d) // d > 0, so Div is a floor.

		av := a.Int64()
		period = append(period, av)
		if av == twoA0 {
			break
		}
	}
	return a0, period
}

// PellFundamental returns the fundamental solution (x, y) of the Pell equation
// x^2 - n*y^2 = 1 as arbitrary-precision integers, derived from the convergents
// of the periodic continued-fraction expansion of sqrt(n). If the period length
// L is even, the convergent with index L-1 is the fundamental solution;
// otherwise the convergent with index 2L-1 is used. The convergents are built
// with the exact h/k recurrence, so the result never overflows. It panics if n
// is a perfect square, for which the equation has no nontrivial solution.
func PellFundamental(n uint64) (x, y *big.Int) {
	a0, period := SqrtContinuedFraction(n)
	if period == nil {
		panic("ntheory: PellFundamental requires a non-square n")
	}
	L := len(period)
	// idx is the 0-based convergent index yielding the fundamental solution.
	idx := L - 1
	if L%2 != 0 {
		idx = 2*L - 1
	}

	// h/k recurrence, seeded by h_{-1}=1, h_{-2}=0, k_{-1}=0, k_{-2}=1, walking
	// the partial quotients a0, then the period cycled as needed.
	hPrev, hPrev2 := big.NewInt(1), big.NewInt(0)
	kPrev, kPrev2 := big.NewInt(0), big.NewInt(1)
	a := new(big.Int)
	for i := 0; i <= idx; i++ {
		if i == 0 {
			a.SetInt64(a0)
		} else {
			a.SetInt64(period[(i-1)%L])
		}
		h := new(big.Int).Mul(a, hPrev)
		h.Add(h, hPrev2)
		k := new(big.Int).Mul(a, kPrev)
		k.Add(k, kPrev2)
		hPrev2, hPrev = hPrev, h
		kPrev2, kPrev = kPrev, k
	}
	return hPrev, kPrev
}
