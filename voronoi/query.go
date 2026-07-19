package voronoi

import (
	"math"
	"sort"
)

// NearestNeighbor returns the index of the point in pts closest to the query
// point q, together with the distance. The index at skip is ignored (pass a
// negative value to consider all points); this makes it convenient to find a
// point's nearest neighbour within its own set. It returns -1 for an empty set.
func NearestNeighbor(pts []Point, q Point, skip int) (int, float64) {
	best := -1
	bestD := math.Inf(1)
	for i, p := range pts {
		if i == skip {
			continue
		}
		if d := p.DistanceSq(q); d < bestD {
			bestD = d
			best = i
		}
	}
	if best < 0 {
		return -1, 0
	}
	return best, math.Sqrt(bestD)
}

// neighborDist pairs an index with a distance for sorting.
type neighborDist struct {
	index int
	dist  float64
}

// KNearest returns the indices of the k points in pts closest to q, ordered
// from nearest to farthest. The index at skip is ignored (pass a negative value
// to include all points). If k exceeds the number of eligible points, all are
// returned.
func KNearest(pts []Point, q Point, k, skip int) []int {
	if k <= 0 {
		return nil
	}
	cand := make([]neighborDist, 0, len(pts))
	for i, p := range pts {
		if i == skip {
			continue
		}
		cand = append(cand, neighborDist{i, p.DistanceSq(q)})
	}
	sort.Slice(cand, func(i, j int) bool {
		if cand[i].dist != cand[j].dist {
			return cand[i].dist < cand[j].dist
		}
		return cand[i].index < cand[j].index
	})
	if k > len(cand) {
		k = len(cand)
	}
	out := make([]int, k)
	for i := 0; i < k; i++ {
		out[i] = cand[i].index
	}
	return out
}

// AllNearestNeighbors returns, for each point in pts, the index of its nearest
// other point. Entry i is -1 only when pts has a single element.
func AllNearestNeighbors(pts []Point) []int {
	out := make([]int, len(pts))
	for i := range pts {
		out[i], _ = NearestNeighbor(pts, pts[i], i)
	}
	return out
}

// NearestNeighborEdges returns the directed nearest-neighbour edges of the set:
// entry i is the edge from point i to its nearest other point. It returns nil
// when there are fewer than two points.
func NearestNeighborEdges(pts []Point) []Edge {
	if len(pts) < 2 {
		return nil
	}
	nn := AllNearestNeighbors(pts)
	out := make([]Edge, 0, len(pts))
	for i, j := range nn {
		if j >= 0 {
			out = append(out, Edge{i, j})
		}
	}
	return out
}

// ClosestPair returns the two indices of the closest pair of points and their
// distance. It uses an O(n log n) divide-and-conquer sweep. It returns
// (-1, -1, +Inf) for fewer than two points.
func ClosestPair(pts []Point) (int, int, float64) {
	n := len(pts)
	if n < 2 {
		return -1, -1, math.Inf(1)
	}
	idx := make([]int, n)
	for i := range idx {
		idx[i] = i
	}
	byX := make([]int, n)
	copy(byX, idx)
	sort.Slice(byX, func(a, b int) bool {
		if pts[byX[a]].X != pts[byX[b]].X {
			return pts[byX[a]].X < pts[byX[b]].X
		}
		return pts[byX[a]].Y < pts[byX[b]].Y
	})
	a, b, d := closestRec(pts, byX)
	return a, b, d
}

// closestRec is the recursive core of ClosestPair over indices sorted by X.
func closestRec(pts []Point, byX []int) (int, int, float64) {
	n := len(byX)
	if n <= 3 {
		bi, bj := -1, -1
		best := math.Inf(1)
		for i := 0; i < n; i++ {
			for j := i + 1; j < n; j++ {
				if d := pts[byX[i]].Distance(pts[byX[j]]); d < best {
					best, bi, bj = d, byX[i], byX[j]
				}
			}
		}
		return bi, bj, best
	}
	mid := n / 2
	midX := pts[byX[mid]].X
	li, lj, ld := closestRec(pts, byX[:mid])
	ri, rj, rd := closestRec(pts, byX[mid:])
	bi, bj, best := li, lj, ld
	if rd < best {
		bi, bj, best = ri, rj, rd
	}
	// Build a strip of points within best of the dividing line, sorted by Y.
	var strip []int
	for _, id := range byX {
		if math.Abs(pts[id].X-midX) < best {
			strip = append(strip, id)
		}
	}
	sort.Slice(strip, func(a, b int) bool { return pts[strip[a]].Y < pts[strip[b]].Y })
	for i := 0; i < len(strip); i++ {
		for j := i + 1; j < len(strip) && pts[strip[j]].Y-pts[strip[i]].Y < best; j++ {
			if d := pts[strip[i]].Distance(pts[strip[j]]); d < best {
				best, bi, bj = d, strip[i], strip[j]
			}
		}
	}
	return bi, bj, best
}

// LargestEmptyCircle returns the largest circle whose centre lies inside the
// convex hull of the sites and whose interior contains no site. The centre is
// the Delaunay circumcentre (Voronoi vertex) that maximises the distance to its
// nearest site, restricted to circumcentres lying within the hull. It returns
// ErrTooFewPoints or ErrDegenerate as Triangulate would.
func LargestEmptyCircle(sites []Point) (Circle, error) {
	tri, err := Triangulate(sites)
	if err != nil {
		return Circle{}, err
	}
	hull := tri.ConvexHull()
	best := Circle{}
	found := false
	for i := range tri.triangles {
		cc, err := tri.Circumcenter(i)
		if err != nil {
			continue
		}
		if !PointInConvexPolygon(cc, hull, Eps) {
			continue
		}
		_, d := NearestNeighbor(tri.points, cc, -1)
		if !found || d > best.Radius {
			best = Circle{Center: cc, Radius: d}
			found = true
		}
	}
	if !found {
		return Circle{}, ErrDegenerate
	}
	return best, nil
}

// EmptyCircleAt returns the largest circle centred at c that contains no site
// in its interior; its radius is the distance from c to the nearest site.
func EmptyCircleAt(sites []Point, c Point) Circle {
	_, d := NearestNeighbor(sites, c, -1)
	return Circle{Center: c, Radius: d}
}

// RelativeNeighborhoodEdges returns the edges (i, j) such that no third point k
// is simultaneously closer to both i and j than they are to each other. This is
// the relative neighbourhood graph, a subgraph of the Delaunay triangulation.
func RelativeNeighborhoodEdges(pts []Point) []Edge {
	n := len(pts)
	var out []Edge
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			dij := pts[i].Distance(pts[j])
			ok := true
			for k := 0; k < n; k++ {
				if k == i || k == j {
					continue
				}
				if pts[i].Distance(pts[k]) < dij && pts[j].Distance(pts[k]) < dij {
					ok = false
					break
				}
			}
			if ok {
				out = append(out, Edge{i, j})
			}
		}
	}
	return out
}

// GabrielEdges returns the edges (i, j) whose diametral circle (the circle with
// segment ij as diameter) contains no other point. This is the Gabriel graph, a
// subgraph of the Delaunay triangulation.
func GabrielEdges(pts []Point) []Edge {
	n := len(pts)
	var out []Edge
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			c := CircleFromDiameter(pts[i], pts[j])
			ok := true
			for k := 0; k < n; k++ {
				if k == i || k == j {
					continue
				}
				if c.ContainsStrict(pts[k], Eps) {
					ok = false
					break
				}
			}
			if ok {
				out = append(out, Edge{i, j})
			}
		}
	}
	return out
}

// EuclideanMST returns the edges of a Euclidean minimum spanning tree of the
// points. The tree's edges are a subset of the Delaunay edges; the algorithm
// builds it from the Delaunay triangulation with Prim's method, falling back to
// the complete graph for degenerate inputs. Edges are returned in the order
// they are added to the tree.
func EuclideanMST(pts []Point) []Edge {
	n := len(pts)
	if n < 2 {
		return nil
	}
	// Candidate edges: Delaunay edges when available, else the complete graph.
	var cand []Edge
	if tri, err := Triangulate(pts); err == nil {
		cand = tri.Edges()
	} else {
		for i := 0; i < n; i++ {
			for j := i + 1; j < n; j++ {
				cand = append(cand, Edge{i, j})
			}
		}
	}
	// Kruskal with union-find.
	sort.Slice(cand, func(a, b int) bool {
		return pts[cand[a].A].DistanceSq(pts[cand[a].B]) <
			pts[cand[b].A].DistanceSq(pts[cand[b].B])
	})
	parent := make([]int, n)
	for i := range parent {
		parent[i] = i
	}
	var find func(int) int
	find = func(x int) int {
		for parent[x] != x {
			parent[x] = parent[parent[x]]
			x = parent[x]
		}
		return x
	}
	var out []Edge
	for _, e := range cand {
		ra, rb := find(e.A), find(e.B)
		if ra != rb {
			parent[ra] = rb
			out = append(out, e.Canonical())
			if len(out) == n-1 {
				break
			}
		}
	}
	return out
}

// TotalEdgeLength returns the sum of the lengths of the given edges over the
// point set.
func TotalEdgeLength(pts []Point, edges []Edge) float64 {
	var s float64
	for _, e := range edges {
		s += pts[e.A].Distance(pts[e.B])
	}
	return s
}
