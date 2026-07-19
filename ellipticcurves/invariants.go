package ellipticcurves

import "math/big"

// IsSupersingular reports whether the curve is supersingular over F_p, i.e. the
// trace of Frobenius a_p is congruent to 0 modulo p. For p > 3 this is
// equivalent to a_p = 0. It uses the naive point count.
func (c *CurveFp) IsSupersingular() bool {
	ap := c.TraceOfFrobenius()
	return new(big.Int).Mod(ap, c.P).Sign() == 0
}

// IsOrdinary reports whether the curve is ordinary over F_p, the negation of
// IsSupersingular.
func (c *CurveFp) IsOrdinary() bool {
	return !c.IsSupersingular()
}

// IsAnomalous reports whether #E(F_p) = p, equivalently a_p = 1. Such curves are
// called anomalous and have a trivial (additive) discrete logarithm.
func (c *CurveFp) IsAnomalous() bool {
	return c.CurveOrderNaive().Cmp(c.P) == 0
}

// HasFullTwoTorsion reports whether all three non-trivial 2-torsion points are
// rational over F_p, equivalently the cubic x^3 + A*x + B splits completely.
// This holds precisely when the group order is divisible by 4 with the 2-part
// non-cyclic; here it is tested directly by counting roots of the cubic.
func (c *CurveFp) HasFullTwoTorsion() bool {
	return len(c.TwoTorsionPoints()) == 3
}

// TwoTorsionPoints returns the non-identity 2-torsion points, the affine points
// with y = 0, i.e. the roots of x^3 + A*x + B. There are zero, one or three
// such points.
func (c *CurveFp) TwoTorsionPoints() []PointFp {
	var pts []PointFp
	x := big.NewInt(0)
	for x.Cmp(c.P) < 0 {
		if c.RightHandSide(x).Sign() == 0 {
			pts = append(pts, PointFp{X: new(big.Int).Set(x), Y: big.NewInt(0)})
		}
		x.Add(x, bigOne)
	}
	return pts
}

// FrobeniusCharacteristicPoly returns the coefficients (1, -a_p, p) of the
// characteristic polynomial T^2 - a_p*T + p of the Frobenius endomorphism, as a
// three-element slice.
func (c *CurveFp) FrobeniusCharacteristicPoly() []*big.Int {
	ap := c.TraceOfFrobenius()
	return []*big.Int{big.NewInt(1), new(big.Int).Neg(ap), new(big.Int).Set(c.P)}
}

// QuadraticTwistFp returns the quadratic twist of the curve by a non-residue d:
// the curve y^2 = x^3 + A*d^2*x + B*d^3. When d is a square the twist is
// isomorphic to the original curve; the interesting case is d a non-residue.
func (c *CurveFp) QuadraticTwistFp(d *big.Int) (*CurveFp, error) {
	d2 := ModSquare(d, c.P)
	d3 := ModMul(d2, d, c.P)
	a := ModMul(c.A, d2, c.P)
	b := ModMul(c.B, d3, c.P)
	return NewCurveFp(a, b, c.P)
}

// TwistOrder returns the order of the quadratic twist of the curve, which is
// 2*(p+1) - #E(F_p) = p + 1 + a_p. It uses the naive point count.
func (c *CurveFp) TwistOrder() *big.Int {
	ap := c.TraceOfFrobenius()
	return new(big.Int).Add(new(big.Int).Add(c.P, bigOne), ap)
}

// SameJInvariantFp reports whether two curves over the same field share a
// j-invariant, a necessary condition for being isomorphic or twists.
func (c *CurveFp) SameJInvariantFp(other *CurveFp) bool {
	if c.P.Cmp(other.P) != 0 {
		return false
	}
	j1, err1 := c.JInvariant()
	j2, err2 := other.JInvariant()
	if err1 != nil || err2 != nil {
		return false
	}
	return j1.Cmp(j2) == 0
}

// IsIsomorphicFp reports whether two curves over the same field are isomorphic
// over F_p, i.e. there exists u != 0 with A' = u^4*A and B' = u^6*B. When A and
// B are both zero the test degenerates to comparing the corresponding
// coefficient; the general special cases (A = 0 or B = 0) are handled directly.
func (c *CurveFp) IsIsomorphicFp(other *CurveFp) bool {
	if c.P.Cmp(other.P) != 0 {
		return false
	}
	u, ok := c.IsomorphismScaleFp(other)
	return ok && u != nil
}

// IsomorphismScaleFp returns a field element u with A' = u^4*A and B' = u^6*B
// witnessing an isomorphism from the receiver to other, and true when such u
// exists. When no isomorphism over F_p exists it returns (nil, false).
func (c *CurveFp) IsomorphismScaleFp(other *CurveFp) (*big.Int, bool) {
	if c.P.Cmp(other.P) != 0 {
		return nil, false
	}
	p := c.P
	aZero := c.A.Sign() == 0
	bZero := c.B.Sign() == 0
	aZero2 := other.A.Sign() == 0
	bZero2 := other.B.Sign() == 0
	if aZero != aZero2 || bZero != bZero2 {
		return nil, false
	}
	switch {
	case aZero && bZero:
		return nil, false // singular, excluded by construction
	case aZero:
		// Need u^6 = B'/B; u exists iff B'/B is a sixth power.
		ratio, err := ModDiv(other.B, c.B, p)
		if err != nil {
			return nil, false
		}
		u := nthRootFp(ratio, big.NewInt(6), p)
		if u == nil {
			return nil, false
		}
		return u, true
	case bZero:
		ratio, err := ModDiv(other.A, c.A, p)
		if err != nil {
			return nil, false
		}
		u := nthRootFp(ratio, big.NewInt(4), p)
		if u == nil {
			return nil, false
		}
		return u, true
	default:
		// u^4 = A'/A and u^6 = B'/B, so u^2 = (u^6)/(u^4) = (B'/B)*(A/A').
		ra, err := ModDiv(other.A, c.A, p)
		if err != nil {
			return nil, false
		}
		rb, err := ModDiv(other.B, c.B, p)
		if err != nil {
			return nil, false
		}
		u2, err := ModDiv(rb, ra, p)
		if err != nil {
			return nil, false
		}
		// Require (u^2)^2 = A'/A.
		if ModSquare(u2, p).Cmp(ra) != 0 {
			return nil, false
		}
		u, err := ModSqrt(u2, p)
		if err != nil {
			return nil, false
		}
		return u, true
	}
}

// IsQuadraticTwistFp reports whether other is a quadratic twist of the receiver
// over F_p: the curves share a j-invariant but are not isomorphic over F_p. For
// j not in {0, 1728} the quadratic twist is the only non-trivial twist.
func (c *CurveFp) IsQuadraticTwistFp(other *CurveFp) bool {
	if !c.SameJInvariantFp(other) {
		return false
	}
	return !c.IsIsomorphicFp(other)
}

// nthRootFp returns an n-th root of a in F_p, or nil when none exists. It works
// by brute force over the cyclic multiplicative group when n is small relative
// to p-1, using the existence criterion a^((p-1)/g) = 1 with g = gcd(n, p-1).
func nthRootFp(a, n, p *big.Int) *big.Int {
	if a.Sign() == 0 {
		return big.NewInt(0)
	}
	pm1 := new(big.Int).Sub(p, bigOne)
	g := Gcd(n, pm1)
	// Existence: a is an n-th power iff a^((p-1)/g) = 1.
	exp := new(big.Int).Div(pm1, g)
	if new(big.Int).Exp(a, exp, p).Cmp(bigOne) != 0 {
		return nil
	}
	// Construct a root via a primitive root.
	gen := PrimitiveRoot(p)
	if gen == nil {
		return nil
	}
	// Find k with gen^k = a using baby-step/giant-step in F_p^*.
	k := discreteLogFp(gen, a, p)
	if k == nil {
		return nil
	}
	// We need u = gen^t with n*t = k mod (p-1). Solve linear congruence.
	if new(big.Int).Mod(k, g).Sign() != 0 {
		return nil
	}
	nn := new(big.Int).Div(n, g)
	kk := new(big.Int).Div(k, g)
	mm := new(big.Int).Div(pm1, g)
	inv := new(big.Int).ModInverse(new(big.Int).Mod(nn, mm), mm)
	if inv == nil {
		return nil
	}
	t := new(big.Int).Mul(kk, inv)
	t.Mod(t, mm)
	return new(big.Int).Exp(gen, t, p)
}

// discreteLogFp returns k with base^k = a in F_p^*, using baby-step/giant-step,
// or nil when a is not in the group (which cannot happen for a primitive base).
func discreteLogFp(base, a, p *big.Int) *big.Int {
	order := new(big.Int).Sub(p, bigOne)
	m := new(big.Int).Add(IntSqrt(order), bigOne)
	table := make(map[string]int64)
	cur := big.NewInt(1)
	mInt := m.Int64()
	for j := int64(0); j < mInt; j++ {
		key := cur.String()
		if _, ok := table[key]; !ok {
			table[key] = j
		}
		cur = ModMul(cur, base, p)
	}
	// factor = base^{-m}
	baseM := new(big.Int).Exp(base, m, p)
	invBaseM := new(big.Int).ModInverse(baseM, p)
	giant := new(big.Int).Set(a)
	for i := int64(0); i <= mInt; i++ {
		key := giant.String()
		if j, ok := table[key]; ok {
			return big.NewInt(i*mInt + j)
		}
		giant = ModMul(giant, invBaseM, p)
	}
	return nil
}
