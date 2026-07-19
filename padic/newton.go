package padic

import (
	"math/big"
	"sort"
)

// NPPoint is a lattice point (X, Y) of a Newton polygon, where X is a
// coefficient index and Y is the p-adic valuation of that coefficient.
type NPPoint struct {
	X int
	Y int
}

// NewtonPolygon is the lower convex hull of the points (i, val_p(a_i)) of a
// polynomial sum a_i x^i, ignoring indices with a_i == 0. Its slopes are the
// negatives of the valuations of the polynomial's roots (Newton's theorem).
type NewtonPolygon struct {
	P        *big.Int
	Points   []NPPoint // finite (non-zero-coefficient) points, sorted by X
	Vertices []NPPoint // vertices of the lower hull, sorted by X
}

// NewtonPolygonFromInts builds the Newton polygon of the integer polynomial
// coeffs (low-to-high) for the prime p.
func NewtonPolygonFromInts(p *big.Int, coeffs []*big.Int) *NewtonPolygon {
	pts := make([]NPPoint, 0, len(coeffs))
	for i, c := range coeffs {
		if c.Sign() == 0 {
			continue
		}
		pts = append(pts, NPPoint{X: i, Y: ValuationInt(p, c)})
	}
	return newNewtonPolygon(p, pts)
}

// NewtonPolygonFromRats builds the Newton polygon of the rational-coefficient
// polynomial coeffs (low-to-high) for the prime p.
func NewtonPolygonFromRats(p *big.Int, coeffs []*big.Rat) *NewtonPolygon {
	pts := make([]NPPoint, 0, len(coeffs))
	for i, c := range coeffs {
		if c.Sign() == 0 {
			continue
		}
		pts = append(pts, NPPoint{X: i, Y: ValuationBigRat(p, c)})
	}
	return newNewtonPolygon(p, pts)
}

// NewtonPolygonFromPadics builds the Newton polygon from p-adic coefficients
// (low-to-high), using their valuations. All coefficients must share the prime.
func NewtonPolygonFromPadics(coeffs []*Padic) (*NewtonPolygon, error) {
	if len(coeffs) == 0 {
		return nil, ErrDomain
	}
	p := coeffs[0].p
	pts := make([]NPPoint, 0, len(coeffs))
	for i, c := range coeffs {
		if c.p.Cmp(p) != 0 {
			return nil, ErrPrimeMismatch
		}
		if c.IsZero() {
			continue
		}
		pts = append(pts, NPPoint{X: i, Y: c.val})
	}
	return newNewtonPolygon(new(big.Int).Set(p), pts), nil
}

// newNewtonPolygon computes the lower convex hull of pts (sorted by X).
func newNewtonPolygon(p *big.Int, pts []NPPoint) *NewtonPolygon {
	sort.Slice(pts, func(i, j int) bool { return pts[i].X < pts[j].X })
	np := &NewtonPolygon{P: p, Points: pts}
	np.Vertices = lowerHull(pts)
	return np
}

// lowerHull returns the vertices of the lower convex hull of points sorted by
// increasing X, using an Andrew-style monotonic chain with rational slope
// comparison.
func lowerHull(pts []NPPoint) []NPPoint {
	var hull []NPPoint
	for _, pt := range pts {
		for len(hull) >= 2 {
			a := hull[len(hull)-2]
			b := hull[len(hull)-1]
			// cross product of (b-a) x (pt-a); <=0 means b not a lower vertex.
			cross := (b.X-a.X)*(pt.Y-a.Y) - (b.Y-a.Y)*(pt.X-a.X)
			if cross <= 0 {
				hull = hull[:len(hull)-1]
			} else {
				break
			}
		}
		hull = append(hull, pt)
	}
	return hull
}

// Slopes returns the slopes of the segments of the lower hull, in increasing
// order, as exact big.Rat values. There is one slope per segment (there are
// len(Vertices)-1 segments).
func (np *NewtonPolygon) Slopes() []*big.Rat {
	if len(np.Vertices) < 2 {
		return nil
	}
	out := make([]*big.Rat, 0, len(np.Vertices)-1)
	for i := 1; i < len(np.Vertices); i++ {
		a := np.Vertices[i-1]
		b := np.Vertices[i]
		out = append(out, new(big.Rat).SetFrac(
			big.NewInt(int64(b.Y-a.Y)), big.NewInt(int64(b.X-a.X))))
	}
	return out
}

// RootValuations returns the p-adic valuations of the roots of the polynomial,
// each repeated by the horizontal length of its segment. By Newton's theorem a
// segment of slope m and horizontal length L contributes L roots of valuation
// -m. Roots are returned in increasing order of valuation.
func (np *NewtonPolygon) RootValuations() []*big.Rat {
	var out []*big.Rat
	for i := 1; i < len(np.Vertices); i++ {
		a := np.Vertices[i-1]
		b := np.Vertices[i]
		slope := new(big.Rat).SetFrac(
			big.NewInt(int64(b.Y-a.Y)), big.NewInt(int64(b.X-a.X)))
		rootVal := new(big.Rat).Neg(slope)
		for k := 0; k < b.X-a.X; k++ {
			out = append(out, new(big.Rat).Set(rootVal))
		}
	}
	return out
}

// SegmentLengths returns the horizontal lengths of the hull's segments, i.e.
// the number of roots of each distinct valuation.
func (np *NewtonPolygon) SegmentLengths() []int {
	if len(np.Vertices) < 2 {
		return nil
	}
	out := make([]int, 0, len(np.Vertices)-1)
	for i := 1; i < len(np.Vertices); i++ {
		out = append(out, np.Vertices[i].X-np.Vertices[i-1].X)
	}
	return out
}

// IsPure reports whether the Newton polygon consists of a single segment, i.e.
// all roots share one valuation (the polynomial is "pure" or isoclinic).
func (np *NewtonPolygon) IsPure() bool {
	return len(np.Vertices) == 2
}

// Valuation returns the polygon's value (minimum height of the hull) at a given
// integer abscissa x, interpolating linearly along the hull segments. It is the
// lower bound val(a_x) implied by convexity.
func (np *NewtonPolygon) Valuation(x int) *big.Rat {
	v := np.Vertices
	if len(v) == 0 {
		return nil
	}
	if x <= v[0].X {
		return new(big.Rat).SetInt64(int64(v[0].Y))
	}
	if x >= v[len(v)-1].X {
		return new(big.Rat).SetInt64(int64(v[len(v)-1].Y))
	}
	for i := 1; i < len(v); i++ {
		if x <= v[i].X {
			a := v[i-1]
			b := v[i]
			// y = a.Y + (x-a.X)*(b.Y-a.Y)/(b.X-a.X)
			slope := new(big.Rat).SetFrac(
				big.NewInt(int64(b.Y-a.Y)), big.NewInt(int64(b.X-a.X)))
			dx := new(big.Rat).SetInt64(int64(x - a.X))
			y := new(big.Rat).Mul(slope, dx)
			y.Add(y, new(big.Rat).SetInt64(int64(a.Y)))
			return y
		}
	}
	return new(big.Rat).SetInt64(int64(v[len(v)-1].Y))
}

// NumRoots returns the total number of roots (with multiplicity) accounted for
// by the Newton polygon, equal to the horizontal width of the lower hull.
func (np *NewtonPolygon) NumRoots() int {
	if len(np.Vertices) < 2 {
		return 0
	}
	return np.Vertices[len(np.Vertices)-1].X - np.Vertices[0].X
}

// Width returns the horizontal extent of the Newton polygon.
func (np *NewtonPolygon) Width() int {
	if len(np.Vertices) < 2 {
		return 0
	}
	return np.Vertices[len(np.Vertices)-1].X - np.Vertices[0].X
}

// LowerVertices returns a copy of the vertices of the lower hull.
func (np *NewtonPolygon) LowerVertices() []NPPoint {
	out := make([]NPPoint, len(np.Vertices))
	copy(out, np.Vertices)
	return out
}

// DistinctSlopes returns the distinct slopes of the Newton polygon in
// increasing order (each root valuation is the negative of a slope).
func (np *NewtonPolygon) DistinctSlopes() []*big.Rat {
	return np.Slopes()
}

// PrincipalSlope returns the smallest (most negative) slope of the polygon,
// corresponding to the roots of largest valuation, or nil if undefined.
func (np *NewtonPolygon) PrincipalSlope() *big.Rat {
	s := np.Slopes()
	if len(s) == 0 {
		return nil
	}
	return s[0]
}
