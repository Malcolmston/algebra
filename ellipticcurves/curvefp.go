package ellipticcurves

import (
	"errors"
	"fmt"
	"math/big"
	"math/rand"
)

// ErrSingularCurve indicates that the discriminant of a curve vanishes, so the
// cubic is singular and does not define an elliptic curve.
var ErrSingularCurve = errors.New("ellipticcurves: singular curve (zero discriminant)")

// ErrPointNotOnCurve indicates that an affine point does not satisfy the curve
// equation.
var ErrPointNotOnCurve = errors.New("ellipticcurves: point is not on the curve")

// CurveFp is an elliptic curve in short Weierstrass form y^2 = x^3 + A*x + B
// over the prime field F_p, with p > 3. The coefficients are stored reduced
// into [0, p).
type CurveFp struct {
	A *big.Int
	B *big.Int
	P *big.Int
}

// PointFp is an affine point on a CurveFp, or the point at infinity when
// Infinity is set. The coordinates are stored reduced into [0, p).
type PointFp struct {
	X        *big.Int
	Y        *big.Int
	Infinity bool
}

// NewCurveFp constructs the curve y^2 = x^3 + A*x + B over F_p. It returns an
// error when p is not an odd prime greater than 3 or when the curve is
// singular.
func NewCurveFp(a, b, p *big.Int) (*CurveFp, error) {
	if p.Sign() <= 0 || !p.ProbablyPrime(20) {
		return nil, ErrNotPrime
	}
	if p.Cmp(bigThree) <= 0 {
		return nil, ErrNotPrime
	}
	c := &CurveFp{
		A: Mod(a, p),
		B: Mod(b, p),
		P: new(big.Int).Set(p),
	}
	if c.IsSingular() {
		return nil, ErrSingularCurve
	}
	return c, nil
}

// PointAtInfinityFp returns the identity element of the group, the point at
// infinity.
func PointAtInfinityFp() PointFp {
	return PointFp{Infinity: true}
}

// NewPointFp constructs an affine point (x, y) on the curve, reducing the
// coordinates modulo p. It returns ErrPointNotOnCurve when the point does not
// satisfy the curve equation.
func (c *CurveFp) NewPointFp(x, y *big.Int) (PointFp, error) {
	p := PointFp{X: Mod(x, c.P), Y: Mod(y, c.P)}
	if !c.IsOnCurve(p) {
		return PointFp{}, ErrPointNotOnCurve
	}
	return p, nil
}

// RightHandSide returns x^3 + A*x + B mod p, the right-hand side of the curve
// equation evaluated at x.
func (c *CurveFp) RightHandSide(x *big.Int) *big.Int {
	xm := Mod(x, c.P)
	res := ModSquare(xm, c.P)
	res = ModMul(res, xm, c.P)                   // x^3
	res = ModAdd(res, ModMul(c.A, xm, c.P), c.P) // + A x
	res = ModAdd(res, c.B, c.P)                  // + B
	return res
}

// IsOnCurve reports whether pt satisfies y^2 = x^3 + A*x + B mod p. The point
// at infinity is always on the curve.
func (c *CurveFp) IsOnCurve(pt PointFp) bool {
	if pt.Infinity {
		return true
	}
	lhs := ModSquare(pt.Y, c.P)
	rhs := c.RightHandSide(pt.X)
	return lhs.Cmp(rhs) == 0
}

// IsSingular reports whether the curve discriminant vanishes modulo p.
func (c *CurveFp) IsSingular() bool {
	return c.Discriminant().Sign() == 0
}

// Discriminant returns the discriminant -16*(4*A^3 + 27*B^2) reduced modulo p.
func (c *CurveFp) Discriminant() *big.Int {
	a3 := ModMul(ModSquare(c.A, c.P), c.A, c.P)
	b2 := ModSquare(c.B, c.P)
	t := ModAdd(ModMul(bigFour, a3, c.P), ModMul(big.NewInt(27), b2, c.P), c.P)
	return ModMul(big.NewInt(-16), t, c.P)
}

// JInvariant returns the j-invariant 1728 * 4*A^3 / (4*A^3 + 27*B^2) reduced
// modulo p. It returns an error only for a singular curve, whose j-invariant is
// undefined.
func (c *CurveFp) JInvariant() (*big.Int, error) {
	a3 := ModMul(ModSquare(c.A, c.P), c.A, c.P)
	fourA3 := ModMul(bigFour, a3, c.P)
	denom := ModAdd(fourA3, ModMul(big.NewInt(27), ModSquare(c.B, c.P), c.P), c.P)
	if denom.Sign() == 0 {
		return nil, ErrSingularCurve
	}
	num := ModMul(big.NewInt(1728), fourA3, c.P)
	return ModDiv(num, denom, c.P)
}

// C4 returns the invariant c4 = -48*A reduced modulo p.
func (c *CurveFp) C4() *big.Int {
	return ModMul(big.NewInt(-48), c.A, c.P)
}

// C6 returns the invariant c6 = -864*B reduced modulo p.
func (c *CurveFp) C6() *big.Int {
	return ModMul(big.NewInt(-864), c.B, c.P)
}

// Equal reports whether two curves have identical coefficients and modulus.
func (c *CurveFp) Equal(other *CurveFp) bool {
	return c.P.Cmp(other.P) == 0 && c.A.Cmp(other.A) == 0 && c.B.Cmp(other.B) == 0
}

// PointEqual reports whether two points are equal, treating all points at
// infinity as equal.
func (c *CurveFp) PointEqual(p, q PointFp) bool {
	if p.Infinity || q.Infinity {
		return p.Infinity == q.Infinity
	}
	return p.X.Cmp(q.X) == 0 && p.Y.Cmp(q.Y) == 0
}

// Neg returns the additive inverse (x, -y) of pt. The inverse of infinity is
// infinity.
func (c *CurveFp) Neg(pt PointFp) PointFp {
	if pt.Infinity {
		return PointAtInfinityFp()
	}
	return PointFp{X: new(big.Int).Set(pt.X), Y: ModNeg(pt.Y, c.P)}
}

// Add returns the group sum p + q using the affine chord-and-tangent law.
func (c *CurveFp) Add(p, q PointFp) PointFp {
	if p.Infinity {
		return clonePointFp(q)
	}
	if q.Infinity {
		return clonePointFp(p)
	}
	if p.X.Cmp(q.X) == 0 {
		// Either doubling or inverse points.
		if ModAdd(p.Y, q.Y, c.P).Sign() == 0 {
			return PointAtInfinityFp()
		}
		return c.Double(p)
	}
	// lambda = (q.Y - p.Y) / (q.X - p.X)
	num := ModSub(q.Y, p.Y, c.P)
	den := ModSub(q.X, p.X, c.P)
	lambda, err := ModDiv(num, den, c.P)
	if err != nil {
		return PointAtInfinityFp()
	}
	return c.lineResult(lambda, p, q)
}

// Double returns 2*pt using the tangent-line slope.
func (c *CurveFp) Double(pt PointFp) PointFp {
	if pt.Infinity {
		return PointAtInfinityFp()
	}
	if pt.Y.Sign() == 0 {
		// Vertical tangent: 2-torsion point doubles to infinity.
		return PointAtInfinityFp()
	}
	num := ModAdd(ModMul(bigThree, ModSquare(pt.X, c.P), c.P), c.A, c.P)
	den := ModDouble(pt.Y, c.P)
	lambda, err := ModDiv(num, den, c.P)
	if err != nil {
		return PointAtInfinityFp()
	}
	return c.lineResult(lambda, pt, pt)
}

// lineResult applies x3 = lambda^2 - x1 - x2, y3 = lambda*(x1 - x3) - y1.
func (c *CurveFp) lineResult(lambda *big.Int, p, q PointFp) PointFp {
	x3 := ModSub(ModSub(ModSquare(lambda, c.P), p.X, c.P), q.X, c.P)
	y3 := ModSub(ModMul(lambda, ModSub(p.X, x3, c.P), c.P), p.Y, c.P)
	return PointFp{X: x3, Y: y3}
}

// Sub returns p - q.
func (c *CurveFp) Sub(p, q PointFp) PointFp {
	return c.Add(p, c.Neg(q))
}

// ScalarMul returns k*pt using left-to-right double-and-add. Negative k is
// handled by negating the point.
func (c *CurveFp) ScalarMul(k *big.Int, pt PointFp) PointFp {
	if pt.Infinity || k.Sign() == 0 {
		return PointAtInfinityFp()
	}
	n := new(big.Int).Abs(k)
	base := pt
	if k.Sign() < 0 {
		base = c.Neg(pt)
	}
	result := PointAtInfinityFp()
	for i := n.BitLen() - 1; i >= 0; i-- {
		result = c.Double(result)
		if n.Bit(i) == 1 {
			result = c.Add(result, base)
		}
	}
	return result
}

// Multiply is an alias for ScalarMul.
func (c *CurveFp) Multiply(k *big.Int, pt PointFp) PointFp {
	return c.ScalarMul(k, pt)
}

// LiftX returns the up-to-two affine points whose x-coordinate equals x. It
// returns an error when x is not the abscissa of any point (the right-hand side
// is a non-residue). When the right-hand side is zero a single 2-torsion point
// is returned twice.
func (c *CurveFp) LiftX(x *big.Int) (PointFp, PointFp, error) {
	rhs := c.RightHandSide(x)
	y, err := ModSqrt(rhs, c.P)
	if err != nil {
		return PointFp{}, PointFp{}, err
	}
	xm := Mod(x, c.P)
	p1 := PointFp{X: xm, Y: y}
	p2 := PointFp{X: new(big.Int).Set(xm), Y: ModNeg(y, c.P)}
	return p1, p2, nil
}

// IsTwoTorsion reports whether pt is a non-identity point of order two, i.e. a
// point with y = 0.
func (c *CurveFp) IsTwoTorsion(pt PointFp) bool {
	return !pt.Infinity && pt.Y.Sign() == 0
}

// RandomPointFp returns a uniformly random affine point on the curve using the
// supplied source. It never returns the point at infinity.
func (c *CurveFp) RandomPointFp(rng *rand.Rand) PointFp {
	for {
		x := new(big.Int).Rand(rng, c.P)
		rhs := c.RightHandSide(x)
		if rhs.Sign() == 0 {
			return PointFp{X: x, Y: big.NewInt(0)}
		}
		if y, err := ModSqrt(rhs, c.P); err == nil {
			if rng.Intn(2) == 1 {
				y = ModNeg(y, c.P)
			}
			return PointFp{X: x, Y: y}
		}
	}
}

// String renders the curve equation and field.
func (c *CurveFp) String() string {
	return fmt.Sprintf("y^2 = x^3 + %s*x + %s over F_%s", c.A, c.B, c.P)
}

// String renders a point in affine coordinates, or "O" for infinity.
func (pt PointFp) String() string {
	if pt.Infinity {
		return "O"
	}
	return fmt.Sprintf("(%s, %s)", pt.X, pt.Y)
}

func clonePointFp(p PointFp) PointFp {
	if p.Infinity {
		return PointAtInfinityFp()
	}
	return PointFp{X: new(big.Int).Set(p.X), Y: new(big.Int).Set(p.Y)}
}
