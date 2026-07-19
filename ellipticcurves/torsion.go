package ellipticcurves

import (
	"errors"
	"math/big"
	"sort"
)

// ErrNonIntegralCurve indicates that a routine requiring integer coefficients
// (such as the Nagell-Lutz torsion search) was called on a curve with
// non-integral A or B.
var ErrNonIntegralCurve = errors.New("ellipticcurves: curve coefficients must be integers")

// maxTorsionOrder is the largest possible order of a rational torsion point by
// Mazur's theorem; the torsion subgroup has order dividing 12 and no element of
// order 11 or greater than 12.
const maxTorsionOrder = 12

// IntCoefficients returns the integer coefficients (A, B) of an integral curve
// and true, or (nil, nil, false) when the coefficients are not integers.
func (c *CurveQ) IntCoefficients() (*big.Int, *big.Int, bool) {
	if !c.IsIntegralCurve() {
		return nil, nil, false
	}
	return new(big.Int).Set(c.A.Num()), new(big.Int).Set(c.B.Num()), true
}

// PointOrderQ returns the order of pt in E(Q), the least k > 0 with k*pt = O, or
// 0 when pt has infinite order (no such k up to Mazur's bound of 12). The point
// at infinity has order 1.
func (c *CurveQ) PointOrderQ(pt PointQ) int {
	if pt.Infinity {
		return 1
	}
	acc := clonePointQ(pt)
	for k := 1; k <= maxTorsionOrder; k++ {
		if acc.Infinity {
			return k
		}
		acc = c.Add(acc, pt)
	}
	return 0
}

// IsTorsionPoint reports whether pt has finite order in E(Q).
func (c *CurveQ) IsTorsionPoint(pt PointQ) bool {
	return c.PointOrderQ(pt) != 0
}

// NagellLutzCandidates returns the integral points (x, y) satisfying the
// Nagell-Lutz necessary conditions for a torsion point: x and y integers with
// either y = 0 or y^2 dividing the discriminant. Not every candidate is
// necessarily torsion; filter with IsTorsionPoint. It returns
// ErrNonIntegralCurve for a non-integral curve.
func (c *CurveQ) NagellLutzCandidates() ([]PointQ, error) {
	a, b, ok := c.IntCoefficients()
	if !ok {
		return nil, ErrNonIntegralCurve
	}
	disc := c.Discriminant()
	// disc is an integer for an integral curve.
	dInt := new(big.Int).Abs(disc.Num())

	var candidates []PointQ
	seen := make(map[string]bool)
	addPoint := func(x, y *big.Int) {
		pt, err := c.NewPointQ(new(big.Rat).SetInt(x), new(big.Rat).SetInt(y))
		if err != nil {
			return
		}
		key := x.String() + ":" + y.String()
		if !seen[key] {
			seen[key] = true
			candidates = append(candidates, pt)
		}
	}

	// y = 0: roots of x^3 + A x + B.
	for _, x := range integerRootsOfCubic(a, b, big.NewInt(0)) {
		addPoint(x, big.NewInt(0))
	}

	// y != 0 with y^2 | disc.
	root := IntSqrt(dInt)
	y := big.NewInt(1)
	for y.Cmp(root) <= 0 {
		y2 := new(big.Int).Mul(y, y)
		if new(big.Int).Mod(dInt, y2).Sign() == 0 {
			for _, x := range integerRootsOfCubic(a, b, y2) {
				addPoint(x, new(big.Int).Set(y))
				addPoint(x, new(big.Int).Neg(y))
			}
		}
		y.Add(y, bigOne)
	}
	return candidates, nil
}

// integerRootsOfCubic returns the integer solutions x of x^3 + A*x + B = target.
// It uses the rational-root theorem: any integer root divides the constant term
// B - target (or is 0 when that constant is 0).
func integerRootsOfCubic(a, b, target *big.Int) []*big.Int {
	c0 := new(big.Int).Sub(b, target)
	eval := func(x *big.Int) *big.Int {
		r := new(big.Int).Mul(x, x)
		r.Mul(r, x)
		r.Add(r, new(big.Int).Mul(a, x))
		r.Add(r, b)
		return r
	}
	var roots []*big.Int
	seen := make(map[string]bool)
	consider := func(x *big.Int) {
		if eval(x).Cmp(target) == 0 {
			key := x.String()
			if !seen[key] {
				seen[key] = true
				roots = append(roots, new(big.Int).Set(x))
			}
		}
	}
	if c0.Sign() == 0 {
		// x^3 + A*x = x*(x^2 + A); roots are 0 and the integer roots of x^2 = -A.
		consider(big.NewInt(0))
		if r, ok := IsPerfectSquare(new(big.Int).Neg(a)); ok {
			consider(r)
			consider(new(big.Int).Neg(r))
		}
	} else {
		for _, d := range Divisors(c0) {
			consider(d)
			consider(new(big.Int).Neg(d))
		}
	}
	sort.Slice(roots, func(i, j int) bool { return roots[i].Cmp(roots[j]) < 0 })
	return roots
}

// TorsionSubgroup returns all torsion points of E(Q), including the point at
// infinity, obtained by filtering the Nagell-Lutz candidates for finite order.
// It returns ErrNonIntegralCurve for a non-integral curve.
func (c *CurveQ) TorsionSubgroup() ([]PointQ, error) {
	candidates, err := c.NagellLutzCandidates()
	if err != nil {
		return nil, err
	}
	points := []PointQ{PointAtInfinityQ()}
	for _, pt := range candidates {
		if c.IsTorsionPoint(pt) {
			points = append(points, pt)
		}
	}
	// Deduplicate (candidates already unique, but be safe) and sort by string.
	sort.Slice(points, func(i, j int) bool { return points[i].String() < points[j].String() })
	return points, nil
}

// TorsionOrder returns the order of the torsion subgroup of E(Q). It returns
// ErrNonIntegralCurve for a non-integral curve.
func (c *CurveQ) TorsionOrder() (int, error) {
	pts, err := c.TorsionSubgroup()
	if err != nil {
		return 0, err
	}
	return len(pts), nil
}

// TorsionExponent returns the exponent of the torsion subgroup, the least
// common multiple of the orders of its elements. It returns ErrNonIntegralCurve
// for a non-integral curve.
func (c *CurveQ) TorsionExponent() (int, error) {
	pts, err := c.TorsionSubgroup()
	if err != nil {
		return 0, err
	}
	exp := big.NewInt(1)
	for _, pt := range pts {
		exp = Lcm(exp, big.NewInt(int64(c.PointOrderQ(pt))))
	}
	return int(exp.Int64()), nil
}

// TwoTorsionPointsQ returns the rational 2-torsion points (the affine points
// with y = 0, i.e. the rational roots of x^3 + A*x + B). It requires an
// integral curve and returns ErrNonIntegralCurve otherwise.
func (c *CurveQ) TwoTorsionPointsQ() ([]PointQ, error) {
	a, b, ok := c.IntCoefficients()
	if !ok {
		return nil, ErrNonIntegralCurve
	}
	var pts []PointQ
	for _, x := range integerRootsOfCubic(a, b, big.NewInt(0)) {
		pt, err := c.NewPointQ(new(big.Rat).SetInt(x), new(big.Rat).SetInt64(0))
		if err == nil {
			pts = append(pts, pt)
		}
	}
	return pts, nil
}
