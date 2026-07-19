package auctions

import (
	"errors"
	"math"
	"sort"
)

// Point is a two-dimensional payoff pair (utilities of players 0 and 1).
type Point struct {
	X, Y float64
}

// Dominates reports whether a weakly dominates b (a.X >= b.X and a.Y >= b.Y)
// with at least one strict inequality.
func Dominates(a, b Point) bool {
	return a.X >= b.X && a.Y >= b.Y && (a.X > b.X || a.Y > b.Y)
}

// cross returns the 2D cross product (b-a) × (c-a).
func cross(a, b, c Point) float64 {
	return (b.X-a.X)*(c.Y-a.Y) - (b.Y-a.Y)*(c.X-a.X)
}

// ConvexHull2D returns the vertices of the convex hull of the input points in
// counter-clockwise order using Andrew's monotone-chain algorithm. Duplicate
// and interior points are dropped. Fewer than three distinct points are
// returned as-is (deduplicated and sorted).
func ConvexHull2D(pts []Point) []Point {
	n := len(pts)
	if n == 0 {
		return nil
	}
	ps := make([]Point, n)
	copy(ps, pts)
	sort.Slice(ps, func(i, j int) bool {
		if ps[i].X != ps[j].X {
			return ps[i].X < ps[j].X
		}
		return ps[i].Y < ps[j].Y
	})
	// deduplicate
	uniq := ps[:0]
	for i, p := range ps {
		if i == 0 || p != ps[i-1] {
			uniq = append(uniq, p)
		}
	}
	ps = uniq
	if len(ps) < 3 {
		return ps
	}
	var hull []Point
	// lower
	for _, p := range ps {
		for len(hull) >= 2 && cross(hull[len(hull)-2], hull[len(hull)-1], p) <= 0 {
			hull = hull[:len(hull)-1]
		}
		hull = append(hull, p)
	}
	// upper
	lower := len(hull) + 1
	for i := len(ps) - 2; i >= 0; i-- {
		p := ps[i]
		for len(hull) >= lower && cross(hull[len(hull)-2], hull[len(hull)-1], p) <= 0 {
			hull = hull[:len(hull)-1]
		}
		hull = append(hull, p)
	}
	return hull[:len(hull)-1]
}

// ParetoFrontier2D returns the Pareto-efficient points of the input set — those
// not dominated by any other point — sorted by increasing X.
func ParetoFrontier2D(pts []Point) []Point {
	var out []Point
	for i, p := range pts {
		dominated := false
		for j, q := range pts {
			if i == j {
				continue
			}
			if Dominates(q, p) {
				dominated = true
				break
			}
		}
		if !dominated {
			out = append(out, p)
		}
	}
	sort.Slice(out, func(i, j int) bool {
		if out[i].X != out[j].X {
			return out[i].X < out[j].X
		}
		return out[i].Y < out[j].Y
	})
	return out
}

// IdealPoint returns the utopia point of the feasible set: the maximum
// achievable X paired with the maximum achievable Y (generally infeasible
// together).
func IdealPoint(pts []Point) (Point, error) {
	if len(pts) == 0 {
		return Point{}, errors.New("auctions: empty feasible set")
	}
	m := Point{X: math.Inf(-1), Y: math.Inf(-1)}
	for _, p := range pts {
		if p.X > m.X {
			m.X = p.X
		}
		if p.Y > m.Y {
			m.Y = p.Y
		}
	}
	return m, nil
}

// NashProduct returns the Nash product (p.X - d.X)(p.Y - d.Y) relative to the
// disagreement point d.
func NashProduct(p, d Point) float64 {
	return (p.X - d.X) * (p.Y - d.Y)
}

// maximizeNashOnSegment maximizes the Nash product over the segment [a,b],
// restricted to points weakly exceeding the disagreement point d. It returns
// the best point and product, and whether any feasible point was found.
func maximizeNashOnSegment(a, b, d Point) (Point, float64, bool) {
	best := Point{}
	bestVal := math.Inf(-1)
	found := false
	consider := func(t float64) {
		if t < 0 || t > 1 {
			return
		}
		p := Point{X: a.X + t*(b.X-a.X), Y: a.Y + t*(b.Y-a.Y)}
		if p.X < d.X-1e-12 || p.Y < d.Y-1e-12 {
			return
		}
		v := NashProduct(p, d)
		if v > bestVal {
			bestVal = v
			best = p
			found = true
		}
	}
	// f(t) = (a.X + t dx - d.X)(a.Y + t dy - d.Y); f'(t)=0 -> linear in t.
	dx := b.X - a.X
	dy := b.Y - a.Y
	ax := a.X - d.X
	ay := a.Y - d.Y
	// f'(t) = dx*(ay + t dy) + dy*(ax + t dx) = (dx*ay + dy*ax) + 2 dx dy t
	denom := 2 * dx * dy
	if math.Abs(denom) > 1e-15 {
		t := -(dx*ay + dy*ax) / denom
		consider(t)
	}
	consider(0)
	consider(1)
	return best, bestVal, found
}

// NashBargainingSolution returns the Nash bargaining solution over the convex
// hull of the feasible payoff set: the point maximizing the product of gains
// (p.X - d.X)(p.Y - d.Y) subject to p.X >= d.X and p.Y >= d.Y. It requires at
// least one feasible point that strictly dominates the disagreement point d.
func NashBargainingSolution(pts []Point, d Point) (Point, error) {
	hull := ConvexHull2D(pts)
	if len(hull) == 0 {
		return Point{}, errors.New("auctions: empty feasible set")
	}
	best := Point{}
	bestVal := math.Inf(-1)
	found := false
	m := len(hull)
	if m == 1 {
		p := hull[0]
		if p.X >= d.X && p.Y >= d.Y {
			return p, nil
		}
		return Point{}, errors.New("auctions: no feasible point dominates the disagreement point")
	}
	for i := 0; i < m; i++ {
		a := hull[i]
		b := hull[(i+1)%m]
		p, v, ok := maximizeNashOnSegment(a, b, d)
		if ok && v > bestVal {
			bestVal = v
			best = p
			found = true
		}
	}
	if !found || bestVal <= 1e-12 {
		return Point{}, errors.New("auctions: no feasible point strictly dominates the disagreement point")
	}
	return best, nil
}

// NashBargainingDiscrete returns the Nash bargaining solution restricted to the
// finite set of feasible points (no convexification): the point maximizing the
// Nash product among those dominating the disagreement point d.
func NashBargainingDiscrete(pts []Point, d Point) (Point, error) {
	best := Point{}
	bestVal := math.Inf(-1)
	found := false
	for _, p := range pts {
		if p.X < d.X || p.Y < d.Y {
			continue
		}
		v := NashProduct(p, d)
		if v > bestVal {
			bestVal = v
			best = p
			found = true
		}
	}
	if !found {
		return Point{}, errors.New("auctions: no feasible point dominates the disagreement point")
	}
	return best, nil
}

// maxLambda returns the largest λ >= 0 such that d + λ·dir lies inside the
// convex hull, using the hull's outward half-plane bounds. It returns false if
// the direction is degenerate or the hull is lower-dimensional.
func maxLambda(d, dir Point, hull []Point) (float64, bool) {
	m := len(hull)
	if m < 3 {
		return 0, false
	}
	best := math.Inf(1)
	for i := 0; i < m; i++ {
		a := hull[i]
		b := hull[(i+1)%m]
		// outward normal for CCW polygon
		nx := b.Y - a.Y
		ny := -(b.X - a.X)
		nd := nx*dir.X + ny*dir.Y
		if nd <= 1e-15 {
			continue
		}
		// n·(d + λ dir - a) <= 0  =>  λ <= -n·(d-a)/(n·dir)
		val := -(nx*(d.X-a.X) + ny*(d.Y-a.Y)) / nd
		if val < best {
			best = val
		}
	}
	if math.IsInf(best, 1) || best < 0 {
		return 0, false
	}
	return best, true
}

// KalaiSmorodinskySolution returns the Kalai-Smorodinsky bargaining solution:
// the Pareto-efficient point on the segment from the disagreement point d to
// the ideal point, equalizing each player's fraction of their maximum possible
// gain.
func KalaiSmorodinskySolution(pts []Point, d Point) (Point, error) {
	hull := ConvexHull2D(pts)
	ideal, err := IdealPoint(pts)
	if err != nil {
		return Point{}, err
	}
	dir := Point{X: ideal.X - d.X, Y: ideal.Y - d.Y}
	if dir.X <= 0 && dir.Y <= 0 {
		return Point{}, errors.New("auctions: ideal point does not dominate disagreement point")
	}
	lam, ok := maxLambda(d, dir, hull)
	if !ok {
		return Point{}, errors.New("auctions: could not locate the frontier (degenerate feasible set)")
	}
	return Point{X: d.X + lam*dir.X, Y: d.Y + lam*dir.Y}, nil
}

// EgalitarianSolution returns the egalitarian bargaining solution: the
// Pareto-efficient point giving both players an equal gain over the
// disagreement point d (moving along the (1,1) direction).
func EgalitarianSolution(pts []Point, d Point) (Point, error) {
	hull := ConvexHull2D(pts)
	dir := Point{X: 1, Y: 1}
	lam, ok := maxLambda(d, dir, hull)
	if !ok {
		return Point{}, errors.New("auctions: could not locate the frontier (degenerate feasible set)")
	}
	return Point{X: d.X + lam, Y: d.Y + lam}, nil
}

// UtilitarianSolution returns a feasible point maximizing the sum of utilities
// X + Y. Because the objective is linear, an optimum is attained at a hull
// vertex.
func UtilitarianSolution(pts []Point) (Point, error) {
	if len(pts) == 0 {
		return Point{}, errors.New("auctions: empty feasible set")
	}
	hull := ConvexHull2D(pts)
	cands := hull
	if len(cands) == 0 {
		cands = pts
	}
	best := cands[0]
	for _, p := range cands[1:] {
		if p.X+p.Y > best.X+best.Y {
			best = p
		}
	}
	return best, nil
}

// IsParetoEfficient reports whether p is not strictly dominated by any point in
// the feasible set.
func IsParetoEfficient(pts []Point, p Point) bool {
	for _, q := range pts {
		if Dominates(q, p) {
			return false
		}
	}
	return true
}
