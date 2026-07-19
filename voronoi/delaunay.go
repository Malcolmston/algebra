package voronoi

import (
	"errors"
	"math"
	"sort"
)

// ErrTooFewPoints is returned when a triangulation is requested for fewer than
// three points.
var ErrTooFewPoints = errors.New("voronoi: need at least 3 points")

// ErrDegenerate is returned when all input points are collinear, so no triangle
// can be formed.
var ErrDegenerate = errors.New("voronoi: all points are collinear")

// Triangulation is a Delaunay triangulation of a planar point set. Its
// triangles reference the input points by index. The zero value is not usable;
// obtain one from Triangulate.
type Triangulation struct {
	points    []Point
	triangles []Triangle
}

// Triangulate computes the Delaunay triangulation of the given points using the
// incremental Bowyer-Watson algorithm. Duplicate points are removed before
// triangulating. It returns ErrTooFewPoints if fewer than three distinct
// points remain and ErrDegenerate if all points are collinear.
func Triangulate(pts []Point) (*Triangulation, error) {
	points := DedupePoints(pts)
	n := len(points)
	if n < 3 {
		return nil, ErrTooFewPoints
	}

	// Determine a super-triangle that strictly contains all points.
	box, err := BoundingBox(points)
	if err != nil {
		return nil, err
	}
	dx := box.Width()
	dy := box.Height()
	dmax := math.Max(dx, dy)
	if dmax == 0 {
		return nil, ErrDegenerate
	}
	mid := box.Center()
	// Make the super-triangle comfortably large.
	m := dmax * 20
	st0 := Point{mid.X - 2*m, mid.Y - m}
	st1 := Point{mid.X + 2*m, mid.Y - m}
	st2 := Point{mid.X, mid.Y + 2*m}

	// Working point list: original points followed by the three super vertices.
	work := make([]Point, n, n+3)
	copy(work, points)
	i0 := n
	i1 := n + 1
	i2 := n + 2
	work = append(work, st0, st1, st2)

	tris := []Triangle{orientCCW(work, Triangle{i0, i1, i2})}

	for pi := 0; pi < n; pi++ {
		p := work[pi]

		// Find all triangles whose circumcircle contains p.
		var bad []int
		for ti, t := range tris {
			if inCircumcircle(work, t, p) {
				bad = append(bad, ti)
			}
		}

		// Collect the boundary of the polygonal cavity: edges that belong to
		// exactly one bad triangle.
		type dedge struct{ a, b int }
		count := make(map[Edge]int, len(bad)*3)
		for _, ti := range bad {
			for _, e := range tris[ti].Edges() {
				count[e.Canonical()]++
			}
		}
		var boundary []dedge
		for _, ti := range bad {
			for _, e := range tris[ti].Edges() {
				if count[e.Canonical()] == 1 {
					boundary = append(boundary, dedge{e.A, e.B})
				}
			}
		}

		// Remove bad triangles (delete by index, largest first).
		sort.Sort(sort.Reverse(sort.IntSlice(bad)))
		for _, ti := range bad {
			last := len(tris) - 1
			tris[ti] = tris[last]
			tris = tris[:last]
		}

		// Re-triangulate the cavity by joining p to each boundary edge.
		for _, e := range boundary {
			tris = append(tris, orientCCW(work, Triangle{e.a, e.b, pi}))
		}
	}

	// Drop triangles that touch any super-triangle vertex.
	final := make([]Triangle, 0, len(tris))
	for _, t := range tris {
		if t.A >= n || t.B >= n || t.C >= n {
			continue
		}
		final = append(final, orientCCW(points, t))
	}

	if len(final) == 0 {
		return nil, ErrDegenerate
	}

	// Sort for deterministic output.
	sort.Slice(final, func(i, j int) bool {
		if final[i].A != final[j].A {
			return final[i].A < final[j].A
		}
		if final[i].B != final[j].B {
			return final[i].B < final[j].B
		}
		return final[i].C < final[j].C
	})

	return &Triangulation{points: points, triangles: final}, nil
}

// DelaunayTriangulation is an alias for Triangulate provided for readability at
// call sites that want the fully qualified name.
func DelaunayTriangulation(pts []Point) (*Triangulation, error) {
	return Triangulate(pts)
}

// orientCCW returns t with its vertices reordered so that they wind
// counterclockwise with respect to the given point slice.
func orientCCW(pts []Point, t Triangle) Triangle {
	if Orient2D(pts[t.A], pts[t.B], pts[t.C]) < 0 {
		t.B, t.C = t.C, t.B
	}
	return t
}

// inCircumcircle reports whether p lies inside the circumcircle of triangle t
// (whose vertices are read from pts).
func inCircumcircle(pts []Point, t Triangle, p Point) bool {
	return InCircleTest(pts[t.A], pts[t.B], pts[t.C], p)
}

// Points returns the (deduplicated) points the triangulation was built on. The
// returned slice must not be modified.
func (tr *Triangulation) Points() []Point { return tr.points }

// Point returns the i-th point of the triangulation.
func (tr *Triangulation) Point(i int) Point { return tr.points[i] }

// NumPoints returns the number of points in the triangulation.
func (tr *Triangulation) NumPoints() int { return len(tr.points) }

// Triangles returns the triangulation's triangles as index triples. The
// returned slice must not be modified.
func (tr *Triangulation) Triangles() []Triangle { return tr.triangles }

// NumTriangles returns the number of triangles in the triangulation.
func (tr *Triangulation) NumTriangles() int { return len(tr.triangles) }

// Triangle returns the i-th triangle.
func (tr *Triangulation) Triangle(i int) Triangle { return tr.triangles[i] }

// TrianglePoints returns the three corner points of triangle i.
func (tr *Triangulation) TrianglePoints(i int) (a, b, c Point) {
	t := tr.triangles[i]
	return tr.points[t.A], tr.points[t.B], tr.points[t.C]
}

// Circumcenter returns the circumcentre of triangle i. This is the
// corresponding Voronoi vertex.
func (tr *Triangulation) Circumcenter(i int) (Point, error) {
	a, b, c := tr.TrianglePoints(i)
	return Circumcenter(a, b, c)
}

// Circumcircle returns the circumcircle of triangle i.
func (tr *Triangulation) Circumcircle(i int) (Circle, error) {
	a, b, c := tr.TrianglePoints(i)
	return Circumcircle(a, b, c)
}

// TriangleArea returns the area of triangle i.
func (tr *Triangulation) TriangleArea(i int) float64 {
	a, b, c := tr.TrianglePoints(i)
	return TriangleArea(a, b, c)
}

// TotalArea returns the sum of the areas of all triangles, i.e. the area of the
// convex hull of the point set.
func (tr *Triangulation) TotalArea() float64 {
	var s float64
	for i := range tr.triangles {
		s += tr.TriangleArea(i)
	}
	return s
}

// Circumcenters returns the circumcentres of all triangles, in triangle order.
// Degenerate triangles (which should not occur in a valid Delaunay
// triangulation) contribute the zero point.
func (tr *Triangulation) Circumcenters() []Point {
	out := make([]Point, len(tr.triangles))
	for i := range tr.triangles {
		if c, err := tr.Circumcenter(i); err == nil {
			out[i] = c
		}
	}
	return out
}

// Edges returns the unique undirected edges of the triangulation, each in
// canonical (A <= B) form, sorted.
func (tr *Triangulation) Edges() []Edge {
	set := make(map[Edge]struct{})
	for _, t := range tr.triangles {
		for _, e := range t.Edges() {
			set[e.Canonical()] = struct{}{}
		}
	}
	out := make([]Edge, 0, len(set))
	for e := range set {
		out = append(out, e)
	}
	sortEdges(out)
	return out
}

// NumEdges returns the number of unique edges in the triangulation.
func (tr *Triangulation) NumEdges() int { return len(tr.Edges()) }

// BoundaryEdges returns the edges that belong to exactly one triangle; these
// form the boundary of the triangulated region (the convex hull).
func (tr *Triangulation) BoundaryEdges() []Edge {
	count := make(map[Edge]int)
	for _, t := range tr.triangles {
		for _, e := range t.Edges() {
			count[e.Canonical()]++
		}
	}
	out := make([]Edge, 0)
	for e, c := range count {
		if c == 1 {
			out = append(out, e)
		}
	}
	sortEdges(out)
	return out
}

// ConvexHull returns the convex hull of the triangulation's points in
// counterclockwise order.
func (tr *Triangulation) ConvexHull() []Point {
	return ConvexHull(tr.points)
}

// Neighbors returns, for each triangle index, the indices of triangles sharing
// an edge with it. The result is indexed by triangle.
func (tr *Triangulation) Neighbors() [][]int {
	edgeToTri := make(map[Edge][]int)
	for ti, t := range tr.triangles {
		for _, e := range t.Edges() {
			c := e.Canonical()
			edgeToTri[c] = append(edgeToTri[c], ti)
		}
	}
	adj := make([][]int, len(tr.triangles))
	seen := make([]map[int]struct{}, len(tr.triangles))
	for i := range seen {
		seen[i] = make(map[int]struct{})
	}
	for _, tris := range edgeToTri {
		for a := 0; a < len(tris); a++ {
			for b := a + 1; b < len(tris); b++ {
				x, y := tris[a], tris[b]
				if _, ok := seen[x][y]; !ok {
					seen[x][y] = struct{}{}
					adj[x] = append(adj[x], y)
				}
				if _, ok := seen[y][x]; !ok {
					seen[y][x] = struct{}{}
					adj[y] = append(adj[y], x)
				}
			}
		}
	}
	for i := range adj {
		sort.Ints(adj[i])
	}
	return adj
}

// VertexTriangles returns, for each point index, the indices of the triangles
// incident to it.
func (tr *Triangulation) VertexTriangles() [][]int {
	out := make([][]int, len(tr.points))
	for ti, t := range tr.triangles {
		out[t.A] = append(out[t.A], ti)
		out[t.B] = append(out[t.B], ti)
		out[t.C] = append(out[t.C], ti)
	}
	return out
}

// VertexDegree returns the number of distinct neighbours of point index v in
// the triangulation graph.
func (tr *Triangulation) VertexDegree(v int) int {
	nb := make(map[int]struct{})
	for _, t := range tr.triangles {
		if !t.HasVertex(v) {
			continue
		}
		for _, u := range t.Vertices() {
			if u != v {
				nb[u] = struct{}{}
			}
		}
	}
	return len(nb)
}

// Adjacency returns the undirected adjacency lists of the point graph induced
// by the triangulation edges. Entry i lists the neighbours of point i, sorted.
func (tr *Triangulation) Adjacency() [][]int {
	adj := make([][]int, len(tr.points))
	seen := make([]map[int]struct{}, len(tr.points))
	for i := range seen {
		seen[i] = make(map[int]struct{})
	}
	add := func(a, b int) {
		if _, ok := seen[a][b]; !ok {
			seen[a][b] = struct{}{}
			adj[a] = append(adj[a], b)
		}
	}
	for _, e := range tr.Edges() {
		add(e.A, e.B)
		add(e.B, e.A)
	}
	for i := range adj {
		sort.Ints(adj[i])
	}
	return adj
}

// Locate returns the index of a triangle that contains the query point p, or
// -1 if p lies outside the triangulated region. Boundary points are reported as
// contained.
func (tr *Triangulation) Locate(p Point) int {
	for i := range tr.triangles {
		a, b, c := tr.TrianglePoints(i)
		if PointInTriangle(p, a, b, c, Eps) {
			return i
		}
	}
	return -1
}

// Contains reports whether p lies within the triangulated region (its convex
// hull).
func (tr *Triangulation) Contains(p Point) bool {
	return tr.Locate(p) >= 0
}

// IsDelaunay reports whether the triangulation satisfies the empty-circumcircle
// property: no point lies strictly inside any triangle's circumcircle. This is
// a validation helper.
func (tr *Triangulation) IsDelaunay() bool {
	for i := range tr.triangles {
		t := tr.triangles[i]
		a, b, c := tr.points[t.A], tr.points[t.B], tr.points[t.C]
		cc, err := Circumcircle(a, b, c)
		if err != nil {
			return false
		}
		for pi, p := range tr.points {
			if pi == t.A || pi == t.B || pi == t.C {
				continue
			}
			if cc.ContainsStrict(p, Eps) {
				return false
			}
		}
	}
	return true
}

// EdgeLengths returns the length of every unique edge, in the order returned by
// Edges.
func (tr *Triangulation) EdgeLengths() []float64 {
	edges := tr.Edges()
	out := make([]float64, len(edges))
	for i, e := range edges {
		out[i] = tr.points[e.A].Distance(tr.points[e.B])
	}
	return out
}

// LongestEdge returns the length of the longest edge and its endpoints as an
// edge. The boolean is false when the triangulation has no edges.
func (tr *Triangulation) LongestEdge() (Edge, float64, bool) {
	edges := tr.Edges()
	if len(edges) == 0 {
		return Edge{}, 0, false
	}
	best := edges[0]
	bl := tr.points[best.A].Distance(tr.points[best.B])
	for _, e := range edges[1:] {
		if l := tr.points[e.A].Distance(tr.points[e.B]); l > bl {
			best, bl = e, l
		}
	}
	return best, bl, true
}

// EulerCharacteristic returns V - E + F for the triangulation, counting the
// unbounded face. For a triangulation of a simply connected region this equals
// 2.
func (tr *Triangulation) EulerCharacteristic() int {
	v := len(tr.points)
	e := len(tr.Edges())
	f := len(tr.triangles) + 1 // +1 for the outer face
	return v - e + f
}

// sortEdges sorts a slice of canonical edges lexicographically in place.
func sortEdges(edges []Edge) {
	sort.Slice(edges, func(i, j int) bool {
		if edges[i].A != edges[j].A {
			return edges[i].A < edges[j].A
		}
		return edges[i].B < edges[j].B
	})
}
