package projectivegeom

import (
	"errors"
	"math"
)

// ErrNotCollinear is returned by cross-ratio routines when their four points
// are not collinear to the requested tolerance.
var ErrNotCollinear = errors.New("projectivegeom: points are not collinear")

// CrossRatio1D returns the cross ratio (a, b; c, d) of four scalars on the
// affine line, defined as ((c-a)(d-b)) / ((c-b)(d-a)). Infinite arguments are
// not supported; use CrossRatio for the projective form.
func CrossRatio1D(a, b, c, d float64) float64 {
	return ((c - a) * (d - b)) / ((c - b) * (d - a))
}

// combine solves c = c0*a + c1*b in the least-squares sense for three
// homogeneous 3-vectors known to be linearly dependent, returning (c0, c1). It
// is used to express a collinear point in the basis of two others.
func combine(a, b, c Vec3) (c0, c1 float64) {
	aa := a.Dot(a)
	ab := a.Dot(b)
	bb := b.Dot(b)
	ac := a.Dot(c)
	bc := b.Dot(c)
	det := aa*bb - ab*ab
	if math.Abs(det) < 1e-300 {
		return 0, 0
	}
	c0 = (ac*bb - bc*ab) / det
	c1 = (aa*bc - ab*ac) / det
	return
}

// CrossRatio returns the projective cross ratio (a, b; c, d) of four collinear
// points. For the harmonic range the value is -1. It returns ErrNotCollinear
// when the four points are not collinear within tol. The value may be infinite
// when d coincides with a projectively.
func CrossRatio(a, b, c, d Point, tol float64) (float64, error) {
	if !Collinear(a, b, c, tol) || !Collinear(a, b, d, tol) {
		return 0, ErrNotCollinear
	}
	// Express c and d in the basis {a, b}: c = c0 a + c1 b, d = d0 a + d1 b.
	c0, c1 := combine(a.V, b.V, c.V)
	d0, d1 := combine(a.V, b.V, d.V)
	// (a,b;c,d) = (c1/c0) / (d1/d0).
	num := c1 * d0
	den := c0 * d1
	if den == 0 {
		return math.Inf(1), nil
	}
	return num / den, nil
}

// CrossRatioLines returns the cross ratio of four concurrent lines, equal to
// the cross ratio of their four intersection points with any transversal. It
// returns ErrNotCollinear (interpreted dually as "not concurrent") when the
// lines are not concurrent within tol.
func CrossRatioLines(a, b, c, d Line, tol float64) (float64, error) {
	return CrossRatio(Point{a.V}, Point{b.V}, Point{c.V}, Point{d.V}, tol)
}

// HarmonicConjugate returns the fourth point d that forms a harmonic range with
// a, b, c, i.e. the unique collinear point with (a, b; c, d) = -1. It returns
// ErrNotCollinear when a, b, c are not collinear within tol, and ErrDegenerate
// when two of them coincide.
func HarmonicConjugate(a, b, c Point, tol float64) (Point, error) {
	if !Collinear(a, b, c, tol) {
		return Point{}, ErrNotCollinear
	}
	c0, c1 := combine(a.V, b.V, c.V)
	if math.Abs(c0) < 1e-300 && math.Abs(c1) < 1e-300 {
		return Point{}, ErrDegenerate
	}
	// With c = c0 a + c1 b, the harmonic conjugate is c0 a - c1 b.
	d := a.V.Scale(c0).Sub(b.V.Scale(c1))
	if d.IsZero(Eps) {
		return Point{}, ErrDegenerate
	}
	return Point{d}, nil
}

// AreHarmonic reports whether the four collinear points form a harmonic range,
// i.e. their cross ratio equals -1 within tol. It returns false (with no error)
// when the points are not collinear.
func AreHarmonic(a, b, c, d Point, tol float64) bool {
	cr, err := CrossRatio(a, b, c, d, tol)
	if err != nil {
		return false
	}
	return math.Abs(cr+1) <= tol*10
}

// CrossRatioPermuteSwapFirstPair returns the cross ratio value that results
// from swapping the first pair, given the original value x = (a,b;c,d). By the
// symmetry (b,a;c,d) = 1/x.
func CrossRatioPermuteSwapFirstPair(x float64) float64 { return 1 / x }

// CrossRatioPermuteSwapPairs returns the cross ratio (c,d;a,b) given x =
// (a,b;c,d); the two are equal, so this returns x unchanged. It exists to
// document the invariance of the cross ratio under swapping the two pairs.
func CrossRatioPermuteSwapPairs(x float64) float64 { return x }

// CrossRatioPermuteSwapInner returns (a,c;b,d) given x = (a,b;c,d), equal to
// 1 - x, one of the six values of the anharmonic group.
func CrossRatioPermuteSwapInner(x float64) float64 { return 1 - x }

// CrossRatioSixValues returns the six values of the cross ratio taken by the
// symmetric group acting on the four points, namely {x, 1/x, 1-x, 1/(1-x),
// x/(x-1), (x-1)/x}.
func CrossRatioSixValues(x float64) [6]float64 {
	return [6]float64{
		x,
		1 / x,
		1 - x,
		1 / (1 - x),
		x / (x - 1),
		(x - 1) / x,
	}
}
