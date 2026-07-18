package ecc

import (
	"fmt"
	"math/big"
)

// PointFp is an affine point on a curve over GF(p). The zero of the group, the
// point at infinity, is represented by Infinity set to true, in which case X
// and Y are ignored.
type PointFp struct {
	// X is the affine x-coordinate, reduced into [0, p) for finite points.
	X *big.Int
	// Y is the affine y-coordinate, reduced into [0, p) for finite points.
	Y *big.Int
	// Infinity marks the identity element (the point at infinity).
	Infinity bool
}

// CurveFp is a short Weierstrass curve y^2 = x^3 + A*x + B defined over the
// prime field GF(P). The parameters must satisfy the non-singularity condition
// 4*A^3 + 27*B^2 != 0 (mod P), which NewCurveFp enforces.
type CurveFp struct {
	// A is the linear coefficient of the curve equation, reduced mod P.
	A *big.Int
	// B is the constant coefficient of the curve equation, reduced mod P.
	B *big.Int
	// P is the prime characterizing the base field GF(P).
	P *big.Int
}

// NewCurveFp constructs a curve y^2 = x^3 + a*x + b over GF(p). The prime p must
// be greater than 3 and the curve must be non-singular (its discriminant must
// be non-zero modulo p); otherwise a non-nil error is returned. The coefficients
// a and b are copied and reduced modulo p, so the caller may reuse them freely.
func NewCurveFp(a, b, p *big.Int) (*CurveFp, error) {
	if p.Sign() <= 0 || p.Cmp(big.NewInt(3)) <= 0 {
		return nil, fmt.Errorf("ecc: NewCurveFp requires a prime p > 3, got %s", p.String())
	}
	c := &CurveFp{
		A: eccMod(a, p),
		B: eccMod(b, p),
		P: new(big.Int).Set(p),
	}
	if c.Discriminant().Sign() == 0 {
		return nil, fmt.Errorf("ecc: singular curve, discriminant is zero mod %s", p.String())
	}
	return c, nil
}

// Discriminant returns the curve discriminant -16*(4*A^3 + 27*B^2) reduced
// modulo P. A curve is non-singular exactly when this value is non-zero.
func (c *CurveFp) Discriminant() *big.Int {
	a3 := ModMul(ModMul(c.A, c.A, c.P), c.A, c.P)
	b2 := ModMul(c.B, c.B, c.P)
	fourA3 := ModMul(big.NewInt(4), a3, c.P)
	term := ModAdd(fourA3, ModMul(big.NewInt(27), b2, c.P), c.P)
	return ModMul(big.NewInt(-16), term, c.P)
}

// JInvariant returns the j-invariant 1728 * 4*A^3 / (4*A^3 + 27*B^2) of the
// curve, reduced modulo P. Curves with equal j-invariant are isomorphic over
// the algebraic closure of GF(P). It panics only if the curve is singular,
// which NewCurveFp precludes.
func (c *CurveFp) JInvariant() *big.Int {
	a3 := ModMul(ModMul(c.A, c.A, c.P), c.A, c.P)
	b2 := ModMul(c.B, c.B, c.P)
	fourA3 := ModMul(big.NewInt(4), a3, c.P)
	denom := ModAdd(fourA3, ModMul(big.NewInt(27), b2, c.P), c.P)
	num := ModMul(big.NewInt(1728), fourA3, c.P)
	j, ok := ModDiv(num, denom, c.P)
	if !ok {
		panic("ecc: JInvariant of a singular curve")
	}
	return j
}

// IsSmooth reports whether the curve is non-singular over GF(P), i.e. whether
// its discriminant is non-zero. Curves built with NewCurveFp are always smooth.
func (c *CurveFp) IsSmooth() bool {
	return c.Discriminant().Sign() != 0
}

// Identity returns the point at infinity, the neutral element of the curve
// group.
func (c *CurveFp) Identity() PointFp {
	return PointFp{Infinity: true}
}

// NewPoint constructs the affine point (x, y) on the curve, reducing the
// coordinates modulo P. It returns an error if the point does not satisfy the
// curve equation.
func (c *CurveFp) NewPoint(x, y *big.Int) (PointFp, error) {
	pt := PointFp{X: eccMod(x, c.P), Y: eccMod(y, c.P)}
	if !c.IsOnCurve(pt) {
		return PointFp{}, fmt.Errorf("ecc: point (%s, %s) is not on the curve", pt.X, pt.Y)
	}
	return pt, nil
}

// RHS returns the right-hand side x^3 + A*x + B of the curve equation evaluated
// at x, reduced modulo P.
func (c *CurveFp) RHS(x *big.Int) *big.Int {
	xm := eccMod(x, c.P)
	x2 := ModMul(xm, xm, c.P)
	x3 := ModMul(x2, xm, c.P)
	return ModAdd(x3, ModAdd(ModMul(c.A, xm, c.P), c.B, c.P), c.P)
}

// IsOnCurve reports whether pt lies on the curve. The point at infinity is
// always considered to be on the curve.
func (c *CurveFp) IsOnCurve(pt PointFp) bool {
	if pt.Infinity {
		return true
	}
	lhs := ModMul(eccMod(pt.Y, c.P), eccMod(pt.Y, c.P), c.P)
	return lhs.Cmp(c.RHS(pt.X)) == 0
}

// Equal reports whether two points on the curve are equal, treating all
// point-at-infinity values as a single element.
func (c *CurveFp) Equal(p, q PointFp) bool {
	if p.Infinity || q.Infinity {
		return p.Infinity == q.Infinity
	}
	return eccMod(p.X, c.P).Cmp(eccMod(q.X, c.P)) == 0 &&
		eccMod(p.Y, c.P).Cmp(eccMod(q.Y, c.P)) == 0
}

// Negate returns the additive inverse -pt, the reflection of pt across the
// x-axis. The inverse of the point at infinity is itself.
func (c *CurveFp) Negate(pt PointFp) PointFp {
	if pt.Infinity {
		return PointFp{Infinity: true}
	}
	return PointFp{X: new(big.Int).Set(eccMod(pt.X, c.P)), Y: ModNeg(pt.Y, c.P)}
}

// Double returns 2*pt using the tangent-line group law. Doubling a point whose
// y-coordinate is zero yields the point at infinity.
func (c *CurveFp) Double(pt PointFp) PointFp {
	if pt.Infinity {
		return PointFp{Infinity: true}
	}
	if eccMod(pt.Y, c.P).Sign() == 0 {
		return PointFp{Infinity: true}
	}
	x := eccMod(pt.X, c.P)
	y := eccMod(pt.Y, c.P)
	// lambda = (3x^2 + A) / (2y)
	num := ModAdd(ModMul(big.NewInt(3), ModMul(x, x, c.P), c.P), c.A, c.P)
	den := ModMul(big.NewInt(2), y, c.P)
	lam, ok := ModDiv(num, den, c.P)
	if !ok {
		return PointFp{Infinity: true}
	}
	return c.eccAffineFromLambda(x, y, x, lam)
}

// Add returns the sum p + q under the elliptic-curve group law. It handles the
// identity, the doubling case, and inverse pairs summing to infinity.
func (c *CurveFp) Add(p, q PointFp) PointFp {
	if p.Infinity {
		return c.eccCopy(q)
	}
	if q.Infinity {
		return c.eccCopy(p)
	}
	px, py := eccMod(p.X, c.P), eccMod(p.Y, c.P)
	qx, qy := eccMod(q.X, c.P), eccMod(q.Y, c.P)
	if px.Cmp(qx) == 0 {
		if ModAdd(py, qy, c.P).Sign() == 0 {
			return PointFp{Infinity: true}
		}
		return c.Double(p)
	}
	// lambda = (qy - py) / (qx - px)
	lam, ok := ModDiv(ModSub(qy, py, c.P), ModSub(qx, px, c.P), c.P)
	if !ok {
		return PointFp{Infinity: true}
	}
	return c.eccAffineFromLambda(px, py, qx, lam)
}

// eccAffineFromLambda computes the resulting affine point given the source
// abscissas x1, x2, the ordinate y1 of the first point, and the chord/tangent
// slope lambda.
func (c *CurveFp) eccAffineFromLambda(x1, y1, x2, lam *big.Int) PointFp {
	x3 := ModSub(ModSub(ModMul(lam, lam, c.P), x1, c.P), x2, c.P)
	y3 := ModSub(ModMul(lam, ModSub(x1, x3, c.P), c.P), y1, c.P)
	return PointFp{X: x3, Y: y3}
}

// eccCopy returns a deep copy of a point so that callers cannot alias internal
// coordinate storage.
func (c *CurveFp) eccCopy(pt PointFp) PointFp {
	if pt.Infinity {
		return PointFp{Infinity: true}
	}
	return PointFp{X: new(big.Int).Set(pt.X), Y: new(big.Int).Set(pt.Y)}
}

// ScalarMul returns k*pt using the constant-structure double-and-add algorithm.
// A negative scalar is handled by negating the point, so k*pt = (-k)*(-pt), and
// a zero scalar yields the point at infinity.
func (c *CurveFp) ScalarMul(k *big.Int, pt PointFp) PointFp {
	if pt.Infinity || k.Sign() == 0 {
		return PointFp{Infinity: true}
	}
	base := pt
	n := new(big.Int).Set(k)
	if n.Sign() < 0 {
		n.Neg(n)
		base = c.Negate(pt)
	}
	result := PointFp{Infinity: true}
	addend := c.eccCopy(base)
	for i := 0; i < n.BitLen(); i++ {
		if n.Bit(i) == 1 {
			result = c.Add(result, addend)
		}
		addend = c.Double(addend)
	}
	return result
}
