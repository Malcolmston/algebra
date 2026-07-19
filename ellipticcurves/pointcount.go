package ellipticcurves

import (
	"errors"
	"math/big"
	"math/rand"
)

// ErrOrderNotFound indicates that a baby-step/giant-step search failed to find
// an order within the expected bound, typically because the inputs violate a
// precondition.
var ErrOrderNotFound = errors.New("ellipticcurves: order not found within bound")

// ErrAmbiguousOrder indicates that a baby-step/giant-step curve-order search
// could not isolate a unique candidate in the Hasse interval.
var ErrAmbiguousOrder = errors.New("ellipticcurves: curve order is ambiguous")

// CountAffinePoints returns the number of affine points on the curve by direct
// enumeration: for each x in F_p it adds 1 + (rhs/p) using the Legendre symbol.
// It runs in O(p) field operations and is exact.
func (c *CurveFp) CountAffinePoints() *big.Int {
	count := big.NewInt(0)
	x := big.NewInt(0)
	for x.Cmp(c.P) < 0 {
		rhs := c.RightHandSide(x)
		switch LegendreSymbol(rhs, c.P) {
		case 0:
			count.Add(count, bigOne)
		case 1:
			count.Add(count, bigTwo)
		}
		x.Add(x, bigOne)
	}
	return count
}

// CurveOrderNaive returns the exact order of the group E(F_p), namely the number
// of affine points plus one for the point at infinity. It is the reference
// implementation used to validate the baby-step/giant-step routines.
func (c *CurveFp) CurveOrderNaive() *big.Int {
	return new(big.Int).Add(c.CountAffinePoints(), bigOne)
}

// TraceOfFrobenius returns a_p = p + 1 - #E(F_p) computed from the naive order.
// By the Hasse bound |a_p| <= 2*sqrt(p).
func (c *CurveFp) TraceOfFrobenius() *big.Int {
	n := c.CurveOrderNaive()
	return new(big.Int).Sub(new(big.Int).Add(c.P, bigOne), n)
}

// HasseInterval returns the closed interval [p+1-2*sqrt(p), p+1+2*sqrt(p)] in
// which #E(F_p) must lie, as (lo, hi).
func (c *CurveFp) HasseInterval() (*big.Int, *big.Int) {
	twoRootP := new(big.Int).Mul(bigTwo, IntSqrt(c.P))
	// Widen by one to be safe against the floor in the integer square root.
	twoRootP.Add(twoRootP, bigOne)
	center := new(big.Int).Add(c.P, bigOne)
	lo := new(big.Int).Sub(center, twoRootP)
	hi := new(big.Int).Add(center, twoRootP)
	if lo.Sign() < 1 {
		lo = big.NewInt(1)
	}
	return lo, hi
}

// EnumeratePoints returns every point of E(F_p) including the point at infinity.
// It is intended for small p; the slice has length #E(F_p).
func (c *CurveFp) EnumeratePoints() []PointFp {
	pts := []PointFp{PointAtInfinityFp()}
	x := big.NewInt(0)
	for x.Cmp(c.P) < 0 {
		rhs := c.RightHandSide(x)
		switch LegendreSymbol(rhs, c.P) {
		case 0:
			pts = append(pts, PointFp{X: new(big.Int).Set(x), Y: big.NewInt(0)})
		case 1:
			y, _ := ModSqrt(rhs, c.P)
			pts = append(pts, PointFp{X: new(big.Int).Set(x), Y: y})
			pts = append(pts, PointFp{X: new(big.Int).Set(x), Y: ModNeg(y, c.P)})
		}
		x.Add(x, bigOne)
	}
	return pts
}

// PointOrderBSGS returns the order of pt, the least m > 0 with m*pt = O, using
// baby-step/giant-step over the Hasse interval followed by prime-factor
// reduction. It is exact for any point on the curve.
func (c *CurveFp) PointOrderBSGS(pt PointFp) (*big.Int, error) {
	if pt.Infinity {
		return big.NewInt(1), nil
	}
	_, hi := c.HasseInterval()
	e, err := c.bsgsMultiple(pt, hi)
	if err != nil {
		return nil, err
	}
	return c.reduceOrder(pt, e), nil
}

// bsgsMultiple finds some positive e <= bound with e*pt = O via
// baby-step/giant-step. It returns ErrOrderNotFound if none exists in range.
func (c *CurveFp) bsgsMultiple(pt PointFp, bound *big.Int) (*big.Int, error) {
	m := new(big.Int).Add(IntSqrt(bound), bigOne)
	// Baby steps: table of j*pt for j in [0, m).
	table := make(map[string]int64)
	cur := PointAtInfinityFp()
	mInt := m.Int64()
	for j := int64(0); j < mInt; j++ {
		key := c.pointKey(cur)
		if _, ok := table[key]; !ok {
			table[key] = j
		}
		cur = c.Add(cur, pt)
	}
	// Giant steps: for i in [1, m], test i*m*pt against baby table.
	mPt := c.ScalarMul(m, pt)
	giant := clonePointFp(mPt)
	for i := int64(1); i <= mInt; i++ {
		key := c.pointKey(giant)
		if j, ok := table[key]; ok {
			e := i*mInt - j
			if e > 0 {
				return big.NewInt(e), nil
			}
		}
		giant = c.Add(giant, mPt)
	}
	return nil, ErrOrderNotFound
}

// reduceOrder returns the exact order of pt given a known multiple e with
// e*pt = O, by dividing out prime factors while the reduced value still
// annihilates pt.
func (c *CurveFp) reduceOrder(pt PointFp, e *big.Int) *big.Int {
	order := new(big.Int).Set(e)
	for prime := range Factorize(e) {
		for new(big.Int).Mod(order, prime).Sign() == 0 {
			q := new(big.Int).Div(order, prime)
			if c.ScalarMul(q, pt).Infinity {
				order = q
			} else {
				break
			}
		}
	}
	return order
}

// pointKey returns a stable map key for a point.
func (c *CurveFp) pointKey(pt PointFp) string {
	if pt.Infinity {
		return "O"
	}
	return pt.X.String() + ":" + pt.Y.String()
}

// PointOrder returns the order of pt computed from the group order: it factors
// #E(F_p) and returns the least divisor n with n*pt = O. This is exact and used
// as a cross-check of PointOrderBSGS.
func (c *CurveFp) PointOrder(pt PointFp, groupOrder *big.Int) *big.Int {
	if pt.Infinity {
		return big.NewInt(1)
	}
	order := new(big.Int).Set(groupOrder)
	for prime, exp := range Factorize(groupOrder) {
		for e := 0; e < exp; e++ {
			cand := new(big.Int).Div(order, prime)
			if c.ScalarMul(cand, pt).Infinity {
				order = cand
			} else {
				break
			}
		}
	}
	return order
}

// DiscreteLogBSGS returns x in [0, order) with x*p = q, using
// baby-step/giant-step, given the order of the base point p. It returns
// ErrOrderNotFound when q is not in the subgroup generated by p.
func (c *CurveFp) DiscreteLogBSGS(p, q PointFp, order *big.Int) (*big.Int, error) {
	m := new(big.Int).Add(IntSqrt(order), bigOne)
	mInt := m.Int64()
	// Baby steps: q - j*p for j in [0, m).
	table := make(map[string]int64)
	cur := clonePointFp(q)
	negP := c.Neg(p)
	for j := int64(0); j < mInt; j++ {
		key := c.pointKey(cur)
		if _, ok := table[key]; !ok {
			table[key] = j
		}
		cur = c.Add(cur, negP)
	}
	// Giant steps: i*m*p for i in [0, m).
	mP := c.ScalarMul(m, p)
	giant := PointAtInfinityFp()
	for i := int64(0); i <= mInt; i++ {
		key := c.pointKey(giant)
		if j, ok := table[key]; ok {
			x := new(big.Int).Add(new(big.Int).Mul(big.NewInt(i), m), big.NewInt(j))
			return x.Mod(x, order), nil
		}
		giant = c.Add(giant, mP)
	}
	return nil, ErrOrderNotFound
}

// CurveOrderBSGS returns #E(F_p) using a Mestre-style baby-step/giant-step
// method: it takes random points, uses their point orders to constrain the
// order to a unique multiple in the Hasse interval, and returns it. The seed
// makes the search reproducible. It agrees with CurveOrderNaive but runs in
// roughly O(p^{1/4}) group operations per point.
func (c *CurveFp) CurveOrderBSGS(rng *rand.Rand) (*big.Int, error) {
	lo, hi := c.HasseInterval()
	lcm := big.NewInt(1)
	for attempt := 0; attempt < 100; attempt++ {
		pt := c.RandomPointFp(rng)
		ord, err := c.PointOrderBSGS(pt)
		if err != nil {
			continue
		}
		lcm = Lcm(lcm, ord)
		// Count multiples of lcm within [lo, hi].
		cand, unique := uniqueMultipleInInterval(lcm, lo, hi)
		if unique {
			return cand, nil
		}
	}
	return nil, ErrAmbiguousOrder
}

// uniqueMultipleInInterval returns the unique multiple of step in [lo, hi] and
// true, or (nil, false) if there are zero or several such multiples.
func uniqueMultipleInInterval(step, lo, hi *big.Int) (*big.Int, bool) {
	if step.Sign() == 0 {
		return nil, false
	}
	// First multiple >= lo.
	q := new(big.Int).Add(lo, new(big.Int).Sub(step, bigOne))
	q.Div(q, step)
	first := new(big.Int).Mul(q, step)
	if first.Cmp(lo) < 0 {
		first.Add(first, step)
	}
	if first.Cmp(hi) > 0 {
		return nil, false
	}
	second := new(big.Int).Add(first, step)
	if second.Cmp(hi) <= 0 {
		return nil, false
	}
	return first, true
}

// GroupExponent returns the exponent of E(F_p), the least common multiple of
// the orders of all points, computed from the group order and structure. For an
// abelian group it equals the largest invariant factor.
func (c *CurveFp) GroupExponent(rng *rand.Rand) (*big.Int, error) {
	n, err := c.CurveOrderBSGS(rng)
	if err != nil {
		n = c.CurveOrderNaive()
	}
	lcm := big.NewInt(1)
	for attempt := 0; attempt < 40; attempt++ {
		pt := c.RandomPointFp(rng)
		lcm = Lcm(lcm, c.PointOrder(pt, n))
	}
	return lcm, nil
}

// IsCyclicFp reports whether E(F_p) is cyclic, tested by comparing the group
// exponent against the group order. It uses randomness from rng.
func (c *CurveFp) IsCyclicFp(rng *rand.Rand) bool {
	n := c.CurveOrderNaive()
	exp, err := c.GroupExponent(rng)
	if err != nil {
		return false
	}
	return exp.Cmp(n) == 0
}

// GroupStructureFp returns the invariant factors (n1, n2) with n1 dividing n2
// and n1*n2 = #E(F_p), so that E(F_p) is isomorphic to Z/n1 x Z/n2. For a
// cyclic group n1 = 1. It uses the naive group order and randomness from rng to
// estimate the exponent n2.
func (c *CurveFp) GroupStructureFp(rng *rand.Rand) (*big.Int, *big.Int) {
	n := c.CurveOrderNaive()
	exp := big.NewInt(1)
	for attempt := 0; attempt < 60; attempt++ {
		pt := c.RandomPointFp(rng)
		exp = Lcm(exp, c.PointOrder(pt, n))
	}
	n1 := new(big.Int).Div(n, exp)
	return n1, exp
}

// FindPointOfOrder returns a point of exact order n, or ErrOrderNotFound if no
// such point exists (n must divide the group order and divide the exponent). It
// searches random points using rng.
func (c *CurveFp) FindPointOfOrder(n *big.Int, rng *rand.Rand) (PointFp, error) {
	order := c.CurveOrderNaive()
	if new(big.Int).Mod(order, n).Sign() != 0 {
		return PointFp{}, ErrOrderNotFound
	}
	cofactor := new(big.Int).Div(order, n)
	for attempt := 0; attempt < 400; attempt++ {
		pt := c.RandomPointFp(rng)
		q := c.ScalarMul(cofactor, pt)
		if q.Infinity {
			continue
		}
		if c.PointOrder(q, order).Cmp(n) == 0 {
			return q, nil
		}
	}
	return PointFp{}, ErrOrderNotFound
}
