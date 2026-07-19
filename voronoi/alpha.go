package voronoi

import "sort"

// AlphaShape is the boundary of an alpha complex: the family of Delaunay
// triangles whose circumradius does not exceed the parameter Alpha, together
// with the edges bounding their union. As Alpha grows the shape ranges from the
// point set (small Alpha) up to the full convex hull (large Alpha).
type AlphaShape struct {
	Alpha     float64
	tri       *Triangulation
	triangles []int  // indices of triangles with circumradius <= Alpha
	edges     []Edge // boundary edges of the alpha complex
}

// AlphaShapeOf computes the alpha shape of the Delaunay triangulation for the
// given alpha radius. A triangle is included when its circumradius is at most
// alpha; boundary edges are the edges incident to exactly one included
// triangle. It returns the same errors as Triangulate.
func AlphaShapeOf(pts []Point, alpha float64) (*AlphaShape, error) {
	tri, err := Triangulate(pts)
	if err != nil {
		return nil, err
	}
	return AlphaShapeFromTriangulation(tri, alpha), nil
}

// AlphaShapeFromTriangulation computes the alpha shape from an existing
// triangulation without retriangulating.
func AlphaShapeFromTriangulation(tri *Triangulation, alpha float64) *AlphaShape {
	as := &AlphaShape{Alpha: alpha, tri: tri}
	for i := range tri.triangles {
		r, err := Circumradius(tri.points[tri.triangles[i].A],
			tri.points[tri.triangles[i].B], tri.points[tri.triangles[i].C])
		if err != nil {
			continue
		}
		if r <= alpha+Eps {
			as.triangles = append(as.triangles, i)
		}
	}
	// Boundary edges: incident to exactly one included triangle.
	count := make(map[Edge]int)
	for _, ti := range as.triangles {
		for _, e := range tri.triangles[ti].Edges() {
			count[e.Canonical()]++
		}
	}
	for e, c := range count {
		if c == 1 {
			as.edges = append(as.edges, e)
		}
	}
	sortEdges(as.edges)
	return as
}

// Triangles returns the indices of the Delaunay triangles that make up the
// alpha complex (those with circumradius <= Alpha).
func (as *AlphaShape) Triangles() []int { return as.triangles }

// NumTriangles returns the number of triangles in the alpha complex.
func (as *AlphaShape) NumTriangles() int { return len(as.triangles) }

// Edges returns the boundary edges of the alpha shape.
func (as *AlphaShape) Edges() []Edge { return as.edges }

// NumEdges returns the number of boundary edges of the alpha shape.
func (as *AlphaShape) NumEdges() int { return len(as.edges) }

// Area returns the total area covered by the alpha complex (the sum of its
// triangle areas).
func (as *AlphaShape) Area() float64 {
	var s float64
	for _, ti := range as.triangles {
		s += as.tri.TriangleArea(ti)
	}
	return s
}

// Perimeter returns the total length of the alpha shape's boundary edges.
func (as *AlphaShape) Perimeter() float64 {
	var s float64
	for _, e := range as.edges {
		s += as.tri.points[e.A].Distance(as.tri.points[e.B])
	}
	return s
}

// BoundaryPoints returns the indices of points that lie on the alpha-shape
// boundary, sorted.
func (as *AlphaShape) BoundaryPoints() []int {
	set := make(map[int]struct{})
	for _, e := range as.edges {
		set[e.A] = struct{}{}
		set[e.B] = struct{}{}
	}
	out := make([]int, 0, len(set))
	for v := range set {
		out = append(out, v)
	}
	sort.Ints(out)
	return out
}

// AlphaShapeEdges is a convenience wrapper returning just the boundary edges of
// the alpha shape of pts for the given alpha.
func AlphaShapeEdges(pts []Point, alpha float64) ([]Edge, error) {
	as, err := AlphaShapeOf(pts, alpha)
	if err != nil {
		return nil, err
	}
	return as.Edges(), nil
}

// CircumradiusSpectrum returns the sorted list of circumradii of all Delaunay
// triangles. The values are useful thresholds for choosing an alpha parameter:
// picking alpha between two consecutive radii changes which triangles are
// included.
func CircumradiusSpectrum(tri *Triangulation) []float64 {
	out := make([]float64, 0, tri.NumTriangles())
	for i := range tri.triangles {
		if r, err := Circumradius(tri.points[tri.triangles[i].A],
			tri.points[tri.triangles[i].B], tri.points[tri.triangles[i].C]); err == nil {
			out = append(out, r)
		}
	}
	sort.Float64s(out)
	return out
}
