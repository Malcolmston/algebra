package ellipticcurves

import "math/big"

// QuadraticTwistQ returns the quadratic twist of the curve by the rational d:
// the curve y^2 = x^3 + A*d^2*x + B*d^3. Twisting by a square gives an
// isomorphic curve; twisting by a non-square gives a genuinely different curve
// with the same j-invariant.
func (c *CurveQ) QuadraticTwistQ(d *big.Rat) (*CurveQ, error) {
	d2 := new(big.Rat).Mul(d, d)
	d3 := new(big.Rat).Mul(d2, d)
	a := new(big.Rat).Mul(c.A, d2)
	b := new(big.Rat).Mul(c.B, d3)
	return NewCurveQ(a, b)
}

// SameJInvariantQ reports whether two rational curves share a j-invariant.
func (c *CurveQ) SameJInvariantQ(other *CurveQ) bool {
	j1, err1 := c.JInvariant()
	j2, err2 := other.JInvariant()
	if err1 != nil || err2 != nil {
		return false
	}
	return j1.Cmp(j2) == 0
}

// IsomorphismScaleQ returns a rational u with A' = u^4*A and B' = u^6*B
// witnessing an isomorphism from the receiver to other over Q, and true when
// such u exists. The special cases A = 0 (j = 0) and B = 0 (j = 1728) are
// handled by requiring the corresponding sixth or fourth power to be rational.
func (c *CurveQ) IsomorphismScaleQ(other *CurveQ) (*big.Rat, bool) {
	aZero := c.A.Sign() == 0
	bZero := c.B.Sign() == 0
	if aZero != (other.A.Sign() == 0) || bZero != (other.B.Sign() == 0) {
		return nil, false
	}
	switch {
	case aZero && bZero:
		return nil, false
	case aZero:
		ratio := new(big.Rat).Quo(other.B, c.B)
		u, ok := rationalNthRoot(ratio, 6)
		return u, ok
	case bZero:
		ratio := new(big.Rat).Quo(other.A, c.A)
		u, ok := rationalNthRoot(ratio, 4)
		return u, ok
	default:
		ra := new(big.Rat).Quo(other.A, c.A)
		rb := new(big.Rat).Quo(other.B, c.B)
		u2 := new(big.Rat).Quo(rb, ra)
		// Require (u^2)^2 = A'/A.
		if new(big.Rat).Mul(u2, u2).Cmp(ra) != 0 {
			return nil, false
		}
		u, ok := rationalSqrt(u2)
		if !ok {
			return nil, false
		}
		return u, true
	}
}

// IsIsomorphicQ reports whether two rational curves are isomorphic over Q.
func (c *CurveQ) IsIsomorphicQ(other *CurveQ) bool {
	_, ok := c.IsomorphismScaleQ(other)
	return ok
}

// IsQuadraticTwistQ reports whether other is a quadratic twist of the receiver
// over Q, i.e. they share a j-invariant but are not isomorphic over Q.
func (c *CurveQ) IsQuadraticTwistQ(other *CurveQ) bool {
	return c.SameJInvariantQ(other) && !c.IsIsomorphicQ(other)
}

// rationalSqrt returns the rational square root of r and true when r is a
// non-negative rational whose numerator and denominator are both perfect
// squares.
func rationalSqrt(r *big.Rat) (*big.Rat, bool) {
	if r.Sign() < 0 {
		return nil, false
	}
	if r.Sign() == 0 {
		return new(big.Rat), true
	}
	ns, ok := IsPerfectSquare(new(big.Int).Abs(r.Num()))
	if !ok {
		return nil, false
	}
	ds, ok := IsPerfectSquare(r.Denom())
	if !ok {
		return nil, false
	}
	return new(big.Rat).SetFrac(ns, ds), true
}

// rationalNthRoot returns the rational n-th root of r and true when it exists as
// a rational; n must be positive. For even n a negative r has no real rational
// root.
func rationalNthRoot(r *big.Rat, n int) (*big.Rat, bool) {
	if n <= 0 {
		return nil, false
	}
	if r.Sign() == 0 {
		return new(big.Rat), true
	}
	neg := r.Sign() < 0
	if neg && n%2 == 0 {
		return nil, false
	}
	numRoot, ok := integerNthRoot(new(big.Int).Abs(r.Num()), n)
	if !ok {
		return nil, false
	}
	denRoot, ok := integerNthRoot(r.Denom(), n)
	if !ok {
		return nil, false
	}
	res := new(big.Rat).SetFrac(numRoot, denRoot)
	if neg {
		res.Neg(res)
	}
	return res, true
}

// integerNthRoot returns the exact integer n-th root of a non-negative integer
// and true when it exists.
func integerNthRoot(v *big.Int, n int) (*big.Int, bool) {
	if v.Sign() < 0 {
		return nil, false
	}
	if v.Sign() == 0 {
		return big.NewInt(0), true
	}
	if n == 1 {
		return new(big.Int).Set(v), true
	}
	nn := big.NewInt(int64(n))
	// Binary search for the root.
	lo := big.NewInt(1)
	hi := new(big.Int).Add(v, bigOne)
	for lo.Cmp(hi) < 0 {
		mid := new(big.Int).Add(lo, hi)
		mid.Rsh(mid, 1)
		p := new(big.Int).Exp(mid, nn, nil)
		switch p.Cmp(v) {
		case 0:
			return mid, true
		case -1:
			lo.Add(mid, bigOne)
		default:
			hi.Set(mid)
		}
	}
	return nil, false
}

// TransformPointFp maps a point under the isomorphism (x, y) -> (u^2*x, u^3*y)
// that carries the curve y^2 = x^3 + A*x + B to y^2 = x^3 + u^4*A*x + u^6*B over
// F_p. It is the inverse-model change of variables used by IsomorphismScaleFp.
func (c *CurveFp) TransformPointFp(u *big.Int, pt PointFp) PointFp {
	if pt.Infinity {
		return PointAtInfinityFp()
	}
	u2 := ModSquare(u, c.P)
	u3 := ModMul(u2, u, c.P)
	return PointFp{X: ModMul(u2, pt.X, c.P), Y: ModMul(u3, pt.Y, c.P)}
}

// ScaleCurveFp returns the curve y^2 = x^3 + u^4*A*x + u^6*B over the same field,
// the image of the receiver under the scaling isomorphism with parameter u.
func (c *CurveFp) ScaleCurveFp(u *big.Int) (*CurveFp, error) {
	u2 := ModSquare(u, c.P)
	u4 := ModSquare(u2, c.P)
	u6 := ModMul(u4, u2, c.P)
	return NewCurveFp(ModMul(u4, c.A, c.P), ModMul(u6, c.B, c.P), c.P)
}
