package ellipticcurves

import (
	"errors"
	"fmt"
	"math/big"
)

// ErrSingularCurveQ indicates that a rational curve has vanishing discriminant.
var ErrSingularCurveQ = errors.New("ellipticcurves: singular rational curve")

// ErrPointNotOnCurveQ indicates that a rational point does not satisfy the
// curve equation.
var ErrPointNotOnCurveQ = errors.New("ellipticcurves: rational point is not on the curve")

// CurveQ is an elliptic curve y^2 = x^3 + A*x + B over the rationals Q, with A
// and B stored as reduced big.Rat values.
type CurveQ struct {
	A *big.Rat
	B *big.Rat
}

// PointQ is an affine rational point on a CurveQ, or the point at infinity when
// Infinity is set.
type PointQ struct {
	X        *big.Rat
	Y        *big.Rat
	Infinity bool
}

// NewCurveQ constructs y^2 = x^3 + A*x + B over Q. It returns ErrSingularCurveQ
// when the discriminant vanishes.
func NewCurveQ(a, b *big.Rat) (*CurveQ, error) {
	c := &CurveQ{A: new(big.Rat).Set(a), B: new(big.Rat).Set(b)}
	if c.Discriminant().Sign() == 0 {
		return nil, ErrSingularCurveQ
	}
	return c, nil
}

// NewCurveQInt constructs y^2 = x^3 + A*x + B over Q from integer coefficients.
func NewCurveQInt(a, b *big.Int) (*CurveQ, error) {
	return NewCurveQ(new(big.Rat).SetInt(a), new(big.Rat).SetInt(b))
}

// PointAtInfinityQ returns the identity element of E(Q).
func PointAtInfinityQ() PointQ {
	return PointQ{Infinity: true}
}

// NewPointQ constructs the affine point (x, y) and checks that it lies on the
// curve, returning ErrPointNotOnCurveQ otherwise.
func (c *CurveQ) NewPointQ(x, y *big.Rat) (PointQ, error) {
	p := PointQ{X: new(big.Rat).Set(x), Y: new(big.Rat).Set(y)}
	if !c.IsOnCurve(p) {
		return PointQ{}, ErrPointNotOnCurveQ
	}
	return p, nil
}

// RightHandSide returns x^3 + A*x + B evaluated at x as a rational.
func (c *CurveQ) RightHandSide(x *big.Rat) *big.Rat {
	x2 := new(big.Rat).Mul(x, x)
	x3 := new(big.Rat).Mul(x2, x)
	res := new(big.Rat).Set(x3)
	res.Add(res, new(big.Rat).Mul(c.A, x))
	res.Add(res, c.B)
	return res
}

// IsOnCurve reports whether pt satisfies the curve equation. Infinity is always
// on the curve.
func (c *CurveQ) IsOnCurve(pt PointQ) bool {
	if pt.Infinity {
		return true
	}
	lhs := new(big.Rat).Mul(pt.Y, pt.Y)
	return lhs.Cmp(c.RightHandSide(pt.X)) == 0
}

// Discriminant returns the discriminant -16*(4*A^3 + 27*B^2) as a rational.
func (c *CurveQ) Discriminant() *big.Rat {
	a3 := new(big.Rat).Mul(new(big.Rat).Mul(c.A, c.A), c.A)
	b2 := new(big.Rat).Mul(c.B, c.B)
	t := new(big.Rat).Add(new(big.Rat).Mul(ratFromInt(4), a3), new(big.Rat).Mul(ratFromInt(27), b2))
	return new(big.Rat).Mul(ratFromInt(-16), t)
}

// JInvariant returns the j-invariant 1728 * 4*A^3 / (4*A^3 + 27*B^2) as a
// rational. It returns ErrSingularCurveQ for a singular curve.
func (c *CurveQ) JInvariant() (*big.Rat, error) {
	a3 := new(big.Rat).Mul(new(big.Rat).Mul(c.A, c.A), c.A)
	fourA3 := new(big.Rat).Mul(ratFromInt(4), a3)
	denom := new(big.Rat).Add(fourA3, new(big.Rat).Mul(ratFromInt(27), new(big.Rat).Mul(c.B, c.B)))
	if denom.Sign() == 0 {
		return nil, ErrSingularCurveQ
	}
	num := new(big.Rat).Mul(ratFromInt(1728), fourA3)
	return new(big.Rat).Quo(num, denom), nil
}

// Neg returns the additive inverse (x, -y) of pt.
func (c *CurveQ) Neg(pt PointQ) PointQ {
	if pt.Infinity {
		return PointAtInfinityQ()
	}
	return PointQ{X: new(big.Rat).Set(pt.X), Y: new(big.Rat).Neg(pt.Y)}
}

// Add returns p + q using the rational chord-and-tangent law.
func (c *CurveQ) Add(p, q PointQ) PointQ {
	if p.Infinity {
		return clonePointQ(q)
	}
	if q.Infinity {
		return clonePointQ(p)
	}
	if p.X.Cmp(q.X) == 0 {
		if new(big.Rat).Add(p.Y, q.Y).Sign() == 0 {
			return PointAtInfinityQ()
		}
		return c.Double(p)
	}
	lambda := new(big.Rat).Quo(new(big.Rat).Sub(q.Y, p.Y), new(big.Rat).Sub(q.X, p.X))
	return c.lineResultQ(lambda, p, q)
}

// Double returns 2*pt using the tangent slope.
func (c *CurveQ) Double(pt PointQ) PointQ {
	if pt.Infinity {
		return PointAtInfinityQ()
	}
	if pt.Y.Sign() == 0 {
		return PointAtInfinityQ()
	}
	num := new(big.Rat).Add(new(big.Rat).Mul(ratFromInt(3), new(big.Rat).Mul(pt.X, pt.X)), c.A)
	den := new(big.Rat).Mul(ratFromInt(2), pt.Y)
	lambda := new(big.Rat).Quo(num, den)
	return c.lineResultQ(lambda, pt, pt)
}

func (c *CurveQ) lineResultQ(lambda *big.Rat, p, q PointQ) PointQ {
	x3 := new(big.Rat).Sub(new(big.Rat).Sub(new(big.Rat).Mul(lambda, lambda), p.X), q.X)
	y3 := new(big.Rat).Sub(new(big.Rat).Mul(lambda, new(big.Rat).Sub(p.X, x3)), p.Y)
	return PointQ{X: x3, Y: y3}
}

// Sub returns p - q.
func (c *CurveQ) Sub(p, q PointQ) PointQ {
	return c.Add(p, c.Neg(q))
}

// ScalarMul returns k*pt using double-and-add. Negative k negates the point.
func (c *CurveQ) ScalarMul(k *big.Int, pt PointQ) PointQ {
	if pt.Infinity || k.Sign() == 0 {
		return PointAtInfinityQ()
	}
	n := new(big.Int).Abs(k)
	base := pt
	if k.Sign() < 0 {
		base = c.Neg(pt)
	}
	result := PointAtInfinityQ()
	for i := n.BitLen() - 1; i >= 0; i-- {
		result = c.Double(result)
		if n.Bit(i) == 1 {
			result = c.Add(result, base)
		}
	}
	return result
}

// PointEqual reports whether two rational points are equal, treating all points
// at infinity as equal.
func (c *CurveQ) PointEqual(p, q PointQ) bool {
	if p.Infinity || q.Infinity {
		return p.Infinity == q.Infinity
	}
	return p.X.Cmp(q.X) == 0 && p.Y.Cmp(q.Y) == 0
}

// IsIntegralPoint reports whether both coordinates of pt are integers. The point
// at infinity is treated as integral.
func (c *CurveQ) IsIntegralPoint(pt PointQ) bool {
	if pt.Infinity {
		return true
	}
	return pt.X.IsInt() && pt.Y.IsInt()
}

// IsIntegralCurve reports whether A and B are integers, a precondition for the
// Nagell-Lutz theorem.
func (c *CurveQ) IsIntegralCurve() bool {
	return c.A.IsInt() && c.B.IsInt()
}

// String renders the rational curve equation.
func (c *CurveQ) String() string {
	return fmt.Sprintf("y^2 = x^3 + %s*x + %s over Q", c.A.RatString(), c.B.RatString())
}

// String renders a rational point, or "O" for infinity.
func (pt PointQ) String() string {
	if pt.Infinity {
		return "O"
	}
	return fmt.Sprintf("(%s, %s)", pt.X.RatString(), pt.Y.RatString())
}

func ratFromInt(n int64) *big.Rat {
	return new(big.Rat).SetInt64(n)
}

func clonePointQ(p PointQ) PointQ {
	if p.Infinity {
		return PointAtInfinityQ()
	}
	return PointQ{X: new(big.Rat).Set(p.X), Y: new(big.Rat).Set(p.Y)}
}
