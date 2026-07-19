package voronoi

import (
	"math"
	"sort"
)

// ConvexHull returns the vertices of the convex hull of the given points in
// counterclockwise order, using Andrew's monotone-chain algorithm. The hull
// starts at the lexicographically smallest point and does not repeat it.
// Collinear points on hull edges are omitted. Fewer than three distinct points
// yield the deduplicated, sorted input.
func ConvexHull(pts []Point) []Point {
	p := SortPointsXY(DedupePoints(pts))
	n := len(p)
	if n < 3 {
		return p
	}

	hull := make([]Point, 0, 2*n)

	// Lower hull.
	for _, pt := range p {
		for len(hull) >= 2 && Orient2D(hull[len(hull)-2], hull[len(hull)-1], pt) <= 0 {
			hull = hull[:len(hull)-1]
		}
		hull = append(hull, pt)
	}

	// Upper hull.
	lower := len(hull) + 1
	for i := n - 2; i >= 0; i-- {
		pt := p[i]
		for len(hull) >= lower && Orient2D(hull[len(hull)-2], hull[len(hull)-1], pt) <= 0 {
			hull = hull[:len(hull)-1]
		}
		hull = append(hull, pt)
	}

	return hull[:len(hull)-1]
}

// ConvexHullIndices returns the indices (into pts) of the convex-hull vertices
// in counterclockwise order. Duplicate coordinates resolve to the first index
// with those coordinates.
func ConvexHullIndices(pts []Point) []int {
	hull := ConvexHull(pts)
	index := make(map[Point]int, len(pts))
	for i, p := range pts {
		if _, ok := index[p]; !ok {
			index[p] = i
		}
	}
	out := make([]int, 0, len(hull))
	for _, h := range hull {
		out = append(out, index[h])
	}
	return out
}

// ConvexHullArea returns the area enclosed by the convex hull of the points.
func ConvexHullArea(pts []Point) float64 {
	return PolygonArea(ConvexHull(pts))
}

// ConvexHullPerimeter returns the perimeter of the convex hull of the points.
func ConvexHullPerimeter(pts []Point) float64 {
	return PolygonPerimeter(ConvexHull(pts))
}

// IsConvexPolygon reports whether the polygon given by its vertices in order is
// convex. It tolerates collinear vertices and either winding direction.
func IsConvexPolygon(poly []Point) bool {
	n := len(poly)
	if n < 3 {
		return false
	}
	var sign int
	for i := 0; i < n; i++ {
		a := poly[i]
		b := poly[(i+1)%n]
		c := poly[(i+2)%n]
		o := Orient2D(a, b, c)
		if o > Eps {
			if sign < 0 {
				return false
			}
			sign = 1
		} else if o < -Eps {
			if sign > 0 {
				return false
			}
			sign = -1
		}
	}
	return true
}

// Diameter returns the maximum distance between any two of the given points
// (the diameter of the set) together with the indices of the realizing pair.
// It evaluates candidate pairs on the convex hull using rotating calipers. It
// returns 0 and (-1, -1) for fewer than two points.
func Diameter(pts []Point) (dist float64, i, j int) {
	if len(pts) < 2 {
		return 0, -1, -1
	}
	hull := ConvexHull(pts)
	m := len(hull)
	if m == 1 {
		return 0, -1, -1
	}
	// Rotating calipers over the hull vertices.
	best := -1.0
	var bi, bj Point
	if m == 2 {
		bi, bj = hull[0], hull[1]
		best = bi.Distance(bj)
	} else {
		k := 1
		for i := 0; i < m; i++ {
			ni := (i + 1) % m
			for {
				nk := (k + 1) % m
				area := math.Abs(Orient2D(hull[i], hull[ni], hull[nk]))
				cur := math.Abs(Orient2D(hull[i], hull[ni], hull[k]))
				if area > cur {
					k = nk
				} else {
					break
				}
			}
			if d := hull[i].Distance(hull[k]); d > best {
				best, bi, bj = d, hull[i], hull[k]
			}
			if d := hull[ni].Distance(hull[k]); d > best {
				best, bi, bj = d, hull[ni], hull[k]
			}
		}
	}
	// Map the realizing hull points back to original indices.
	i, j = indexOfPoint(pts, bi), indexOfPoint(pts, bj)
	return best, i, j
}

// indexOfPoint returns the index of the first exact match of p in pts, or -1.
func indexOfPoint(pts []Point, p Point) int {
	for i, q := range pts {
		if q.Equal(p) {
			return i
		}
	}
	return -1
}

// GrahamHullByAngle returns the convex hull computed by sorting points by polar
// angle about the lowest point. The result matches ConvexHull; it is provided
// as an alternative construction. Vertices are returned counterclockwise.
func GrahamHullByAngle(pts []Point) []Point {
	p := DedupePoints(pts)
	n := len(p)
	if n < 3 {
		return SortPointsXY(p)
	}
	// Pivot: lowest Y, then lowest X.
	piv := 0
	for i := 1; i < n; i++ {
		if p[i].Y < p[piv].Y || (p[i].Y == p[piv].Y && p[i].X < p[piv].X) {
			piv = i
		}
	}
	p[0], p[piv] = p[piv], p[0]
	pivot := p[0]
	rest := p[1:]
	sort.Slice(rest, func(i, j int) bool {
		o := Orient2D(pivot, rest[i], rest[j])
		if o != 0 {
			return o > 0
		}
		return pivot.DistanceSq(rest[i]) < pivot.DistanceSq(rest[j])
	})
	stack := []Point{pivot}
	for _, pt := range rest {
		for len(stack) >= 2 && Orient2D(stack[len(stack)-2], stack[len(stack)-1], pt) <= 0 {
			stack = stack[:len(stack)-1]
		}
		stack = append(stack, pt)
	}
	return stack
}
