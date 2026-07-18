package ecc

import "math/big"

// CountPointsNaive returns the order #E(GF(p)) of the curve group, that is the
// number of affine points plus one for the point at infinity, by enumerating
// every x in GF(p) and counting the corresponding square roots. Its running
// time is proportional to p, so it is intended for small primes and for
// validating faster routines.
func (c *CurveFp) CountPointsNaive() *big.Int {
	count := big.NewInt(1) // point at infinity
	x := new(big.Int)
	for x.Sign() >= 0 && x.Cmp(c.P) < 0 {
		rhs := c.RHS(x)
		switch Legendre(rhs, c.P) {
		case 0:
			count.Add(count, big.NewInt(1)) // single root y = 0
		case 1:
			count.Add(count, big.NewInt(2)) // two roots
		}
		x.Add(x, big.NewInt(1))
	}
	return count
}

// TraceOfFrobenius returns the trace of Frobenius t = p + 1 - #E(GF(p)), where
// #E is the group order. By the Hasse bound |t| <= 2*sqrt(p).
func (c *CurveFp) TraceOfFrobenius() *big.Int {
	order := c.CountPointsNaive()
	t := new(big.Int).Add(c.P, big.NewInt(1))
	return t.Sub(t, order)
}

// HasseInterval returns the closed interval [lo, hi] of admissible group orders
// #E(GF(p)) guaranteed by the Hasse theorem: lo = p + 1 - 2*floor(sqrt(p)) - 1
// and hi = p + 1 + 2*floor(sqrt(p)) + 1, a superset of the true bound that is
// safe to use as an integer search window.
func (c *CurveFp) HasseInterval() (lo, hi *big.Int) {
	s := eccISqrt(c.P)
	twoS := new(big.Int).Mul(big.NewInt(2), s)
	twoS.Add(twoS, big.NewInt(1)) // pad to cover the fractional part of sqrt(p)
	pp1 := new(big.Int).Add(c.P, big.NewInt(1))
	lo = new(big.Int).Sub(pp1, twoS)
	if lo.Sign() < 0 {
		lo = big.NewInt(0)
	}
	hi = new(big.Int).Add(pp1, twoS)
	return lo, hi
}

// PointOrderNaive returns the order of the point pt in the curve group: the
// least positive n with n*pt equal to the point at infinity. It repeatedly adds
// pt, bounded by the upper end of the Hasse interval, and returns nil if no such
// n is found within that bound (which cannot happen for a point actually on the
// curve). The point at infinity has order 1.
func (c *CurveFp) PointOrderNaive(pt PointFp) *big.Int {
	if pt.Infinity {
		return big.NewInt(1)
	}
	_, hi := c.HasseInterval()
	acc := c.eccCopy(pt)
	n := big.NewInt(1)
	for n.Cmp(hi) <= 0 {
		if acc.Infinity {
			return new(big.Int).Set(n)
		}
		acc = c.Add(acc, pt)
		n.Add(n, big.NewInt(1))
	}
	return nil
}

// PointOrderBSGS returns the order of the point pt using the baby-step
// giant-step algorithm, which runs in time proportional to sqrt(#E) rather than
// #E. It first locates a multiple of the order inside the Hasse window and then
// reduces that multiple to the exact order by dividing out prime factors. The
// point at infinity has order 1.
func (c *CurveFp) PointOrderBSGS(pt PointFp) *big.Int {
	if pt.Infinity {
		return big.NewInt(1)
	}
	_, hi := c.HasseInterval()
	m := new(big.Int).Add(eccISqrt(hi), big.NewInt(1)) // m = ceil-ish sqrt(bound)

	// Baby steps: table[key(j*pt)] = j for j in [0, m).
	table := make(map[string]*big.Int)
	acc := PointFp{Infinity: true}
	for j := big.NewInt(0); j.Cmp(m) < 0; j.Add(j, big.NewInt(1)) {
		table[eccPointKey(c, acc)] = new(big.Int).Set(j)
		acc = c.Add(acc, pt)
	}

	// Giant steps: search for i,j with (i*m + j)*pt == O, i.e.
	// -(i*m*pt) == j*pt.
	mPt := c.ScalarMul(m, pt)
	giant := PointFp{Infinity: true} // i*m*pt for i = 0,1,...
	var multiple *big.Int
	for i := big.NewInt(0); i.Cmp(m) <= 0; i.Add(i, big.NewInt(1)) {
		neg := c.Negate(giant)
		if j, ok := table[eccPointKey(c, neg)]; ok {
			t := new(big.Int).Mul(i, m)
			t.Add(t, j)
			if t.Sign() > 0 {
				multiple = t
				break
			}
		}
		giant = c.Add(giant, mPt)
	}
	if multiple == nil {
		return nil
	}
	return c.eccReduceOrder(multiple, pt)
}

// eccReduceOrder reduces a known multiple of the point order to the exact order
// by repeatedly dividing out prime factors for which the reduced multiple still
// annihilates pt.
func (c *CurveFp) eccReduceOrder(multiple *big.Int, pt PointFp) *big.Int {
	n := new(big.Int).Set(multiple)
	for _, q := range eccPrimeFactors(n) {
		for {
			quo := new(big.Int).Div(n, q)
			if new(big.Int).Mod(n, q).Sign() != 0 {
				break
			}
			if c.ScalarMul(quo, pt).Infinity {
				n = quo
			} else {
				break
			}
		}
	}
	return n
}

// eccPointKey returns a canonical string key uniquely identifying a point,
// usable as a map key.
func eccPointKey(c *CurveFp, pt PointFp) string {
	if pt.Infinity {
		return "inf"
	}
	return eccMod(pt.X, c.P).String() + "," + eccMod(pt.Y, c.P).String()
}

// eccPrimeFactors returns the distinct prime factors of n (n > 0) by trial
// division. It is used only to reduce order multiples that are bounded by the
// Hasse window, so trial division is adequate.
func eccPrimeFactors(n *big.Int) []*big.Int {
	var factors []*big.Int
	m := new(big.Int).Set(n)
	d := big.NewInt(2)
	for new(big.Int).Mul(d, d).Cmp(m) <= 0 {
		if new(big.Int).Mod(m, d).Sign() == 0 {
			factors = append(factors, new(big.Int).Set(d))
			for new(big.Int).Mod(m, d).Sign() == 0 {
				m.Div(m, d)
			}
		}
		d.Add(d, big.NewInt(1))
	}
	if m.Cmp(big.NewInt(1)) > 0 {
		factors = append(factors, new(big.Int).Set(m))
	}
	return factors
}
