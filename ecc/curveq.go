package ecc

import (
	"fmt"
	"math/big"
)

// PointQ is an affine point on a curve over the rational field Q. The zero of
// the group, the point at infinity, is represented by Infinity set to true, in
// which case X and Y are ignored.
type PointQ struct {
	// X is the affine x-coordinate as an exact rational for finite points.
	X *big.Rat
	// Y is the affine y-coordinate as an exact rational for finite points.
	Y *big.Rat
	// Infinity marks the identity element (the point at infinity).
	Infinity bool
}

// CurveQ is a short Weierstrass curve y^2 = x^3 + A*x + B defined over the
// rational numbers Q. The parameters must satisfy 4*A^3 + 27*B^2 != 0, which
// NewCurveQ enforces.
type CurveQ struct {
	// A is the linear coefficient of the curve equation.
	A *big.Rat
	// B is the constant coefficient of the curve equation.
	B *big.Rat
}

// NewCurveQ constructs a curve y^2 = x^3 + a*x + b over Q. The curve must be
// non-singular (its discriminant must be non-zero); otherwise a non-nil error
// is returned. The coefficients are copied, so the caller may reuse them.
func NewCurveQ(a, b *big.Rat) (*CurveQ, error) {
	c := &CurveQ{A: new(big.Rat).Set(a), B: new(big.Rat).Set(b)}
	if c.Discriminant().Sign() == 0 {
		return nil, fmt.Errorf("ecc: singular curve over Q, discriminant is zero")
	}
	return c, nil
}

// Discriminant returns the exact rational discriminant -16*(4*A^3 + 27*B^2) of
// the curve over Q. A curve is non-singular exactly when this value is non-zero.
func (c *CurveQ) Discriminant() *big.Rat {
	a3 := new(big.Rat).Mul(new(big.Rat).Mul(c.A, c.A), c.A)
	b2 := new(big.Rat).Mul(c.B, c.B)
	fourA3 := new(big.Rat).Mul(big.NewRat(4, 1), a3)
	term := new(big.Rat).Add(fourA3, new(big.Rat).Mul(big.NewRat(27, 1), b2))
	return new(big.Rat).Mul(big.NewRat(-16, 1), term)
}

// JInvariant returns the exact rational j-invariant
// 1728 * 4*A^3 / (4*A^3 + 27*B^2) of the curve over Q. It panics only if the
// curve is singular, which NewCurveQ precludes.
func (c *CurveQ) JInvariant() *big.Rat {
	a3 := new(big.Rat).Mul(new(big.Rat).Mul(c.A, c.A), c.A)
	b2 := new(big.Rat).Mul(c.B, c.B)
	fourA3 := new(big.Rat).Mul(big.NewRat(4, 1), a3)
	denom := new(big.Rat).Add(fourA3, new(big.Rat).Mul(big.NewRat(27, 1), b2))
	if denom.Sign() == 0 {
		panic("ecc: JInvariant of a singular curve over Q")
	}
	num := new(big.Rat).Mul(big.NewRat(1728, 1), fourA3)
	return new(big.Rat).Quo(num, denom)
}

// Identity returns the point at infinity, the neutral element of the curve
// group over Q.
func (c *CurveQ) Identity() PointQ {
	return PointQ{Infinity: true}
}

// RHS returns the exact value x^3 + A*x + B of the curve equation's right-hand
// side evaluated at the rational x.
func (c *CurveQ) RHS(x *big.Rat) *big.Rat {
	x2 := new(big.Rat).Mul(x, x)
	x3 := new(big.Rat).Mul(x2, x)
	ax := new(big.Rat).Mul(c.A, x)
	return new(big.Rat).Add(x3, new(big.Rat).Add(ax, c.B))
}

// IsOnCurve reports whether the rational point pt lies on the curve. The point
// at infinity is always considered to be on the curve.
func (c *CurveQ) IsOnCurve(pt PointQ) bool {
	if pt.Infinity {
		return true
	}
	lhs := new(big.Rat).Mul(pt.Y, pt.Y)
	return lhs.Cmp(c.RHS(pt.X)) == 0
}

// NewPoint constructs the affine rational point (x, y) on the curve. It returns
// an error if the point does not satisfy the curve equation.
func (c *CurveQ) NewPoint(x, y *big.Rat) (PointQ, error) {
	pt := PointQ{X: new(big.Rat).Set(x), Y: new(big.Rat).Set(y)}
	if !c.IsOnCurve(pt) {
		return PointQ{}, fmt.Errorf("ecc: rational point (%s, %s) is not on the curve", x.RatString(), y.RatString())
	}
	return pt, nil
}

// Equal reports whether two rational points are equal, treating all
// point-at-infinity values as a single element.
func (c *CurveQ) Equal(p, q PointQ) bool {
	if p.Infinity || q.Infinity {
		return p.Infinity == q.Infinity
	}
	return p.X.Cmp(q.X) == 0 && p.Y.Cmp(q.Y) == 0
}

// Negate returns the additive inverse -pt, the reflection of pt across the
// x-axis. The inverse of the point at infinity is itself.
func (c *CurveQ) Negate(pt PointQ) PointQ {
	if pt.Infinity {
		return PointQ{Infinity: true}
	}
	return PointQ{X: new(big.Rat).Set(pt.X), Y: new(big.Rat).Neg(pt.Y)}
}

// Double returns 2*pt over Q using the tangent-line group law. Doubling a point
// with zero y-coordinate yields the point at infinity.
func (c *CurveQ) Double(pt PointQ) PointQ {
	if pt.Infinity || pt.Y.Sign() == 0 {
		return PointQ{Infinity: true}
	}
	// lambda = (3x^2 + A) / (2y)
	x2 := new(big.Rat).Mul(pt.X, pt.X)
	num := new(big.Rat).Add(new(big.Rat).Mul(big.NewRat(3, 1), x2), c.A)
	den := new(big.Rat).Mul(big.NewRat(2, 1), pt.Y)
	lam := new(big.Rat).Quo(num, den)
	return c.eccAffineFromLambda(pt.X, pt.Y, pt.X, lam)
}

// Add returns the sum p + q over Q under the elliptic-curve group law. It
// handles the identity, the doubling case, and inverse pairs summing to
// infinity.
func (c *CurveQ) Add(p, q PointQ) PointQ {
	if p.Infinity {
		return c.eccCopy(q)
	}
	if q.Infinity {
		return c.eccCopy(p)
	}
	if p.X.Cmp(q.X) == 0 {
		if new(big.Rat).Add(p.Y, q.Y).Sign() == 0 {
			return PointQ{Infinity: true}
		}
		return c.Double(p)
	}
	// lambda = (qy - py) / (qx - px)
	num := new(big.Rat).Sub(q.Y, p.Y)
	den := new(big.Rat).Sub(q.X, p.X)
	lam := new(big.Rat).Quo(num, den)
	return c.eccAffineFromLambda(p.X, p.Y, q.X, lam)
}

// eccAffineFromLambda computes the resulting affine rational point from the
// source abscissas x1, x2, the ordinate y1, and the slope lambda.
func (c *CurveQ) eccAffineFromLambda(x1, y1, x2, lam *big.Rat) PointQ {
	lam2 := new(big.Rat).Mul(lam, lam)
	x3 := new(big.Rat).Sub(new(big.Rat).Sub(lam2, x1), x2)
	y3 := new(big.Rat).Sub(new(big.Rat).Mul(lam, new(big.Rat).Sub(x1, x3)), y1)
	return PointQ{X: x3, Y: y3}
}

// eccCopy returns a deep copy of a rational point so callers cannot alias
// internal coordinate storage.
func (c *CurveQ) eccCopy(pt PointQ) PointQ {
	if pt.Infinity {
		return PointQ{Infinity: true}
	}
	return PointQ{X: new(big.Rat).Set(pt.X), Y: new(big.Rat).Set(pt.Y)}
}

// ScalarMul returns k*pt over Q using double-and-add. A negative scalar negates
// the point, and a zero scalar yields the point at infinity. Because rational
// coordinates grow rapidly under repeated addition, this is intended for small
// scalars.
func (c *CurveQ) ScalarMul(k *big.Int, pt PointQ) PointQ {
	if pt.Infinity || k.Sign() == 0 {
		return PointQ{Infinity: true}
	}
	base := pt
	n := new(big.Int).Set(k)
	if n.Sign() < 0 {
		n.Neg(n)
		base = c.Negate(pt)
	}
	result := PointQ{Infinity: true}
	addend := c.eccCopy(base)
	for i := 0; i < n.BitLen(); i++ {
		if n.Bit(i) == 1 {
			result = c.Add(result, addend)
		}
		addend = c.Double(addend)
	}
	return result
}
