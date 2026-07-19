package contfrac

import (
	"math"
	"math/big"
)

// BestApproximation returns the rational p/q closest to x whose denominator q
// satisfies 1 <= q <= maxDen. Among equally close fractions the one with the
// smaller denominator is returned. It uses the continued-fraction algorithm, so
// the answer is always a convergent or a best semiconvergent of x. maxDen must
// be at least 1.
func BestApproximation(x float64, maxDen int64) (int64, int64) {
	if maxDen < 1 {
		maxDen = 1
	}
	neg := false
	if x < 0 {
		neg = true
		x = -x
	}
	hPrev, hPrev2 := int64(1), int64(0)
	kPrev, kPrev2 := int64(0), int64(1)
	bestP, bestQ := int64(0), int64(1)
	f := x
	for {
		a := int64(math.Floor(f))
		h := a*hPrev + hPrev2
		k := a*kPrev + kPrev2
		if k > maxDen || k <= 0 {
			// The next convergent overshoots maxDen. The best candidate is
			// either the previous convergent (hPrev/kPrev) or the largest
			// admissible semiconvergent.
			cp, cq := hPrev, kPrev
			if kPrev > 0 {
				aPrime := (maxDen - kPrev2) / kPrev
				sp := aPrime*hPrev + hPrev2
				sq := aPrime*kPrev + kPrev2
				if sq >= 1 && sq <= maxDen && closer(sp, sq, cp, cq, x) {
					cp, cq = sp, sq
				}
			}
			if cq >= 1 {
				bestP, bestQ = cp, cq
			}
			break
		}
		hPrev2, hPrev = hPrev, h
		kPrev2, kPrev = kPrev, k
		bestP, bestQ = h, k
		frac := f - float64(a)
		if frac < 1e-15 {
			break
		}
		f = 1 / frac
	}
	if neg {
		bestP = -bestP
	}
	return bestP, bestQ
}

// closer reports whether p1/q1 is strictly closer to x than p2/q2, breaking
// exact ties in favour of the smaller denominator.
func closer(p1, q1, p2, q2 int64, x float64) bool {
	d1 := math.Abs(float64(p1)/float64(q1) - x)
	d2 := math.Abs(float64(p2)/float64(q2) - x)
	if d1 == d2 {
		return q1 < q2
	}
	return d1 < d2
}

// BestApproximationFrac returns [BestApproximation] as a reduced [Frac].
func BestApproximationFrac(x float64, maxDen int64) Frac {
	p, q := BestApproximation(x, maxDen)
	return NewFrac(p, q)
}

// Rationalize is an alias for [BestApproximationFrac]: it returns the simplest
// good rational approximation of x with denominator at most maxDen.
func Rationalize(x float64, maxDen int64) Frac {
	return BestApproximationFrac(x, maxDen)
}

// BestApproximationRat returns the rational closest to the target *big.Rat whose
// denominator does not exceed maxDen, computed exactly from the continued
// fraction of the target.
func BestApproximationRat(target *big.Rat, maxDen int64) *big.Rat {
	cf := FromRationalBig(target)
	hPrev, hPrev2 := big.NewInt(1), big.NewInt(0)
	kPrev, kPrev2 := big.NewInt(0), big.NewInt(1)
	maxD := big.NewInt(maxDen)
	best := new(big.Rat).SetFrac(big.NewInt(0), big.NewInt(1))
	haveBest := false
	for _, a := range cf {
		ab := big.NewInt(a)
		h := new(big.Int).Add(new(big.Int).Mul(ab, hPrev), hPrev2)
		k := new(big.Int).Add(new(big.Int).Mul(ab, kPrev), kPrev2)
		if k.Cmp(maxD) > 0 {
			// semiconvergent with the largest admissible partial quotient
			if kPrev.Sign() > 0 {
				aPrime := new(big.Int).Div(new(big.Int).Sub(maxD, kPrev2), kPrev)
				sp := new(big.Int).Add(new(big.Int).Mul(aPrime, hPrev), hPrev2)
				sq := new(big.Int).Add(new(big.Int).Mul(aPrime, kPrev), kPrev2)
				cand := new(big.Rat).SetFrac(sp, sq)
				prev := new(big.Rat).SetFrac(new(big.Int).Set(hPrev), new(big.Int).Set(kPrev))
				if sq.Sign() > 0 && ratCloser(cand, prev, target) {
					best, haveBest = cand, true
				} else {
					best, haveBest = prev, true
				}
			}
			break
		}
		best = new(big.Rat).SetFrac(new(big.Int).Set(h), new(big.Int).Set(k))
		haveBest = true
		hPrev, hPrev2 = h, hPrev
		kPrev, kPrev2 = k, kPrev
	}
	if !haveBest {
		return new(big.Rat)
	}
	return best
}

// ratCloser reports whether a is at least as close to target as b (strictly
// closer, or equally close with a smaller denominator).
func ratCloser(a, b, target *big.Rat) bool {
	da := new(big.Rat).Sub(a, target)
	da.Abs(da)
	db := new(big.Rat).Sub(b, target)
	db.Abs(db)
	switch da.Cmp(db) {
	case -1:
		return true
	case 1:
		return false
	default:
		return a.Denom().Cmp(b.Denom()) < 0
	}
}

// BestApproximations returns the sequence of convergents of x whose
// denominators do not exceed maxDen. Each convergent is a best rational
// approximation of the second kind, so the list gives successively better
// approximations with growing denominators.
func BestApproximations(x float64, maxDen int64) []Frac {
	cf := FromFloat(x, 64)
	var out []Frac
	for _, c := range cf.Convergents() {
		if c.Den >= 1 && c.Den <= maxDen {
			out = append(out, c)
		} else if c.Den > maxDen {
			break
		}
	}
	return out
}

// SmallestDenominatorWithin returns the fraction with the smallest possible
// denominator whose value lies within eps of x. It searches increasing
// denominator bounds using the continued-fraction best-approximation routine.
// eps must be positive.
func SmallestDenominatorWithin(x, eps float64) Frac {
	if eps <= 0 {
		eps = 1e-12
	}
	for q := int64(1); ; q++ {
		p, qq := BestApproximation(x, q)
		if math.Abs(float64(p)/float64(qq)-x) <= eps {
			return NewFrac(p, qq)
		}
		if q > 1<<40 { // safety valve; unreachable for sane eps
			return NewFrac(p, qq)
		}
	}
}
