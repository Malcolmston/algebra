package voronoi

// Edge is an undirected edge between two vertices identified by their indices
// into a point slice.
type Edge struct {
	A, B int
}

// NewEdge returns the edge joining vertex indices a and b.
func NewEdge(a, b int) Edge { return Edge{A: a, B: b} }

// Canonical returns the edge with its endpoints ordered so that A <= B. This
// gives a stable key for undirected edges.
func (e Edge) Canonical() Edge {
	if e.A > e.B {
		return Edge{A: e.B, B: e.A}
	}
	return e
}

// Reversed returns the edge with its endpoints swapped.
func (e Edge) Reversed() Edge { return Edge{A: e.B, B: e.A} }

// Equal reports whether e and f join the same pair of vertices, ignoring
// direction.
func (e Edge) Equal(f Edge) bool {
	return (e.A == f.A && e.B == f.B) || (e.A == f.B && e.B == f.A)
}

// Contains reports whether vertex index v is an endpoint of e.
func (e Edge) Contains(v int) bool { return e.A == v || e.B == v }

// Other returns the endpoint of e that is not v. It returns -1 if v is not an
// endpoint of e.
func (e Edge) Other(v int) int {
	switch v {
	case e.A:
		return e.B
	case e.B:
		return e.A
	default:
		return -1
	}
}

// Triangle references three vertices by their indices into a point slice.
type Triangle struct {
	A, B, C int
}

// NewTriangle returns the triangle on vertex indices a, b and c.
func NewTriangle(a, b, c int) Triangle { return Triangle{A: a, B: b, C: c} }

// Vertices returns the triangle's three vertex indices as a slice.
func (t Triangle) Vertices() []int { return []int{t.A, t.B, t.C} }

// Edges returns the triangle's three directed edges in order A->B, B->C, C->A.
func (t Triangle) Edges() [3]Edge {
	return [3]Edge{{t.A, t.B}, {t.B, t.C}, {t.C, t.A}}
}

// HasVertex reports whether v is one of the triangle's vertices.
func (t Triangle) HasVertex(v int) bool {
	return t.A == v || t.B == v || t.C == v
}

// HasEdge reports whether e is one of the triangle's edges (ignoring
// direction).
func (t Triangle) HasEdge(e Edge) bool {
	for _, te := range t.Edges() {
		if te.Equal(e) {
			return true
		}
	}
	return false
}

// OppositeVertex returns the triangle vertex not incident to edge e. It returns
// -1 if e is not an edge of the triangle.
func (t Triangle) OppositeVertex(e Edge) int {
	for _, v := range t.Vertices() {
		if !e.Contains(v) {
			return v
		}
	}
	return -1
}

// SharesEdge reports whether triangles t and u share an edge, and returns that
// shared edge. The boolean is false when they do not.
func (t Triangle) SharesEdge(u Triangle) (Edge, bool) {
	for _, te := range t.Edges() {
		for _, ue := range u.Edges() {
			if te.Equal(ue) {
				return te.Canonical(), true
			}
		}
	}
	return Edge{}, false
}
