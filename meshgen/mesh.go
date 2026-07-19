package meshgen

import (
	"fmt"
	"math"
	"sort"
)

// Tri is a triangle described by three indices into a vertex slice.
type Tri struct {
	A, B, C int
}

// NewTri returns the triangle with the given vertex indices.
func NewTri(a, b, c int) Tri { return Tri{a, b, c} }

// Indices returns the three vertex indices as an array.
func (t Tri) Indices() [3]int { return [3]int{t.A, t.B, t.C} }

// Has reports whether the triangle references vertex index v.
func (t Tri) Has(v int) bool { return t.A == v || t.B == v || t.C == v }

// Opposite returns the vertex of the triangle that is neither u nor v, and a
// boolean reporting whether exactly one such vertex exists.
func (t Tri) Opposite(u, v int) (int, bool) {
	switch {
	case !t.Has(u) || !t.Has(v):
		return 0, false
	case t.A != u && t.A != v:
		return t.A, true
	case t.B != u && t.B != v:
		return t.B, true
	default:
		return t.C, true
	}
}

// EdgesOf returns the three (unordered, sorted) edges of the triangle.
func (t Tri) EdgesOf() [3]Edge {
	return [3]Edge{
		NewEdge(t.A, t.B),
		NewEdge(t.B, t.C),
		NewEdge(t.C, t.A),
	}
}

// SharesEdge reports whether triangles t and u share an edge, returning that
// edge and true when they do.
func (t Tri) SharesEdge(u Tri) (Edge, bool) {
	for _, e := range t.EdgesOf() {
		for _, f := range u.EdgesOf() {
			if e == f {
				return e, true
			}
		}
	}
	return Edge{}, false
}

// Edge is an unordered pair of vertex indices with U <= V.
type Edge struct {
	U, V int
}

// NewEdge returns the canonical edge for the pair (a, b) with U <= V.
func NewEdge(a, b int) Edge {
	if a > b {
		a, b = b, a
	}
	return Edge{a, b}
}

// Other returns the endpoint of the edge that is not v, and true when v is an
// endpoint of the edge.
func (e Edge) Other(v int) (int, bool) {
	switch v {
	case e.U:
		return e.V, true
	case e.V:
		return e.U, true
	default:
		return 0, false
	}
}

// String returns a human-readable representation of the edge.
func (e Edge) String() string { return fmt.Sprintf("(%d-%d)", e.U, e.V) }

// Segment is a straight line segment between two planar points.
type Segment struct {
	A, B Vec2
}

// NewSegment returns the segment from a to b.
func NewSegment(a, b Vec2) Segment { return Segment{a, b} }

// Length returns the Euclidean length of the segment.
func (s Segment) Length() float64 { return s.A.Distance(s.B) }

// Midpoint returns the midpoint of the segment.
func (s Segment) Midpoint() Vec2 { return s.A.Midpoint(s.B) }

// Direction returns the unit direction from A to B.
func (s Segment) Direction() Vec2 { return s.B.Sub(s.A).Normalize() }

// Mesh is a triangular mesh: a slice of vertices and a slice of triangles that
// index into it.
type Mesh struct {
	Vertices  []Vec2
	Triangles []Tri
}

// NewMesh returns a mesh with the given vertices and triangles.
func NewMesh(verts []Vec2, tris []Tri) *Mesh {
	return &Mesh{Vertices: verts, Triangles: tris}
}

// NumVertices returns the number of vertices in the mesh.
func (m *Mesh) NumVertices() int { return len(m.Vertices) }

// NumTriangles returns the number of triangles in the mesh.
func (m *Mesh) NumTriangles() int { return len(m.Triangles) }

// TriangleVertices returns the three vertex positions of triangle i.
func (m *Mesh) TriangleVertices(i int) (a, b, c Vec2) {
	t := m.Triangles[i]
	return m.Vertices[t.A], m.Vertices[t.B], m.Vertices[t.C]
}

// TriangleArea returns the area of triangle i.
func (m *Mesh) TriangleArea(i int) float64 {
	a, b, c := m.TriangleVertices(i)
	return TriangleArea(a, b, c)
}

// TotalArea returns the sum of the areas of all triangles.
func (m *Mesh) TotalArea() float64 {
	var s float64
	for i := range m.Triangles {
		s += m.TriangleArea(i)
	}
	return s
}

// TriangleCentroid returns the centroid of triangle i.
func (m *Mesh) TriangleCentroid(i int) Vec2 {
	a, b, c := m.TriangleVertices(i)
	return a.Add(b).Add(c).Div(3)
}

// Edges returns every distinct edge of the mesh, sorted.
func (m *Mesh) Edges() []Edge {
	set := make(map[Edge]struct{})
	for _, t := range m.Triangles {
		for _, e := range t.EdgesOf() {
			set[e] = struct{}{}
		}
	}
	out := make([]Edge, 0, len(set))
	for e := range set {
		out = append(out, e)
	}
	sortEdges(out)
	return out
}

// NumEdges returns the number of distinct edges in the mesh.
func (m *Mesh) NumEdges() int { return len(m.Edges()) }

// EdgeTriangleCount maps each edge to the number of triangles incident to it.
func (m *Mesh) EdgeTriangleCount() map[Edge]int {
	counts := make(map[Edge]int)
	for _, t := range m.Triangles {
		for _, e := range t.EdgesOf() {
			counts[e]++
		}
	}
	return counts
}

// BoundaryEdges returns the edges incident to exactly one triangle, sorted.
func (m *Mesh) BoundaryEdges() []Edge {
	counts := m.EdgeTriangleCount()
	var out []Edge
	for e, n := range counts {
		if n == 1 {
			out = append(out, e)
		}
	}
	sortEdges(out)
	return out
}

// InteriorEdges returns the edges shared by two or more triangles, sorted.
func (m *Mesh) InteriorEdges() []Edge {
	counts := m.EdgeTriangleCount()
	var out []Edge
	for e, n := range counts {
		if n >= 2 {
			out = append(out, e)
		}
	}
	sortEdges(out)
	return out
}

// IsBoundaryEdge reports whether the edge (a,b) lies on the mesh boundary.
func (m *Mesh) IsBoundaryEdge(a, b int) bool {
	return m.EdgeTriangleCount()[NewEdge(a, b)] == 1
}

// BoundaryVertices returns the sorted indices of vertices lying on the mesh
// boundary.
func (m *Mesh) BoundaryVertices() []int {
	set := make(map[int]struct{})
	for _, e := range m.BoundaryEdges() {
		set[e.U] = struct{}{}
		set[e.V] = struct{}{}
	}
	out := make([]int, 0, len(set))
	for v := range set {
		out = append(out, v)
	}
	sort.Ints(out)
	return out
}

// IsManifold reports whether every edge is shared by at most two triangles.
func (m *Mesh) IsManifold() bool {
	for _, n := range m.EdgeTriangleCount() {
		if n > 2 {
			return false
		}
	}
	return true
}

// VertexTriangles returns, for each vertex index, the sorted list of triangle
// indices incident to it.
func (m *Mesh) VertexTriangles() [][]int {
	adj := make([][]int, len(m.Vertices))
	for ti, t := range m.Triangles {
		for _, v := range t.Indices() {
			if v >= 0 && v < len(adj) {
				adj[v] = append(adj[v], ti)
			}
		}
	}
	return adj
}

// VertexNeighbors returns, for each vertex, the sorted set of adjacent vertex
// indices connected by an edge.
func (m *Mesh) VertexNeighbors() [][]int {
	sets := make([]map[int]struct{}, len(m.Vertices))
	for i := range sets {
		sets[i] = make(map[int]struct{})
	}
	for _, e := range m.Edges() {
		sets[e.U][e.V] = struct{}{}
		sets[e.V][e.U] = struct{}{}
	}
	out := make([][]int, len(m.Vertices))
	for i, s := range sets {
		lst := make([]int, 0, len(s))
		for v := range s {
			lst = append(lst, v)
		}
		sort.Ints(lst)
		out[i] = lst
	}
	return out
}

// VertexDegree returns the number of edges incident to vertex v.
func (m *Mesh) VertexDegree(v int) int {
	return len(m.VertexNeighbors()[v])
}

// TriangleAdjacency returns, for each triangle, the sorted indices of triangles
// sharing an edge with it.
func (m *Mesh) TriangleAdjacency() [][]int {
	edgeTris := make(map[Edge][]int)
	for ti, t := range m.Triangles {
		for _, e := range t.EdgesOf() {
			edgeTris[e] = append(edgeTris[e], ti)
		}
	}
	adj := make([][]int, len(m.Triangles))
	for _, tris := range edgeTris {
		for _, a := range tris {
			for _, b := range tris {
				if a != b {
					adj[a] = append(adj[a], b)
				}
			}
		}
	}
	for i := range adj {
		adj[i] = uniqueSortedInts(adj[i])
	}
	return adj
}

// EulerCharacteristic returns V - E + F for the mesh, where F counts the
// triangular faces.
func (m *Mesh) EulerCharacteristic() int {
	return len(m.Vertices) - m.NumEdges() + len(m.Triangles)
}

// ConnectedComponents returns the number of connected components of the mesh
// graph, counting only vertices referenced by at least one triangle.
func (m *Mesh) ConnectedComponents() int {
	n := len(m.Vertices)
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
	union := func(a, b int) {
		ra, rb := find(a), find(b)
		if ra != rb {
			parent[ra] = rb
		}
	}
	used := make([]bool, n)
	for _, t := range m.Triangles {
		used[t.A], used[t.B], used[t.C] = true, true, true
		union(t.A, t.B)
		union(t.B, t.C)
	}
	roots := make(map[int]struct{})
	for i := 0; i < n; i++ {
		if used[i] {
			roots[find(i)] = struct{}{}
		}
	}
	return len(roots)
}

// BoundaryLength returns the total length of all boundary edges.
func (m *Mesh) BoundaryLength() float64 {
	var s float64
	for _, e := range m.BoundaryEdges() {
		s += m.Vertices[e.U].Distance(m.Vertices[e.V])
	}
	return s
}

// EdgeLengths returns the length of every distinct edge, in edge order.
func (m *Mesh) EdgeLengths() []float64 {
	edges := m.Edges()
	out := make([]float64, len(edges))
	for i, e := range edges {
		out[i] = m.Vertices[e.U].Distance(m.Vertices[e.V])
	}
	return out
}

// OrientCCW returns a copy of the mesh with every triangle re-ordered to be
// counterclockwise.
func (m *Mesh) OrientCCW() *Mesh {
	tris := make([]Tri, len(m.Triangles))
	for i, t := range m.Triangles {
		a, b, c := m.Vertices[t.A], m.Vertices[t.B], m.Vertices[t.C]
		if Orient2D(a, b, c) < 0 {
			t = Tri{t.A, t.C, t.B}
		}
		tris[i] = t
	}
	return &Mesh{Vertices: m.Vertices, Triangles: tris}
}

// Clone returns a deep copy of the mesh.
func (m *Mesh) Clone() *Mesh {
	v := make([]Vec2, len(m.Vertices))
	copy(v, m.Vertices)
	t := make([]Tri, len(m.Triangles))
	copy(t, m.Triangles)
	return &Mesh{Vertices: v, Triangles: t}
}

// Bounds returns the axis-aligned bounding box of the mesh vertices.
func (m *Mesh) Bounds() (min, max Vec2) { return BoundingBox2(m.Vertices) }

// Validate reports the first structural error found in the mesh, or nil when
// every triangle references distinct, in-range vertices.
func (m *Mesh) Validate() error {
	n := len(m.Vertices)
	for i, t := range m.Triangles {
		for _, v := range t.Indices() {
			if v < 0 || v >= n {
				return fmt.Errorf("meshgen: triangle %d references vertex %d out of range", i, v)
			}
		}
		if t.A == t.B || t.B == t.C || t.A == t.C {
			return fmt.Errorf("meshgen: triangle %d has repeated vertex", i)
		}
	}
	return nil
}

func sortEdges(e []Edge) {
	sort.Slice(e, func(i, j int) bool {
		if e[i].U != e[j].U {
			return e[i].U < e[j].U
		}
		return e[i].V < e[j].V
	})
}

func uniqueSortedInts(in []int) []int {
	if len(in) == 0 {
		return in
	}
	sort.Ints(in)
	out := in[:1]
	for _, v := range in[1:] {
		if v != out[len(out)-1] {
			out = append(out, v)
		}
	}
	return out
}

// MinFloat returns the minimum of the slice, or +Inf for an empty slice.
func MinFloat(xs []float64) float64 {
	m := math.Inf(1)
	for _, x := range xs {
		if x < m {
			m = x
		}
	}
	return m
}

// MaxFloat returns the maximum of the slice, or -Inf for an empty slice.
func MaxFloat(xs []float64) float64 {
	m := math.Inf(-1)
	for _, x := range xs {
		if x > m {
			m = x
		}
	}
	return m
}

// MeanFloat returns the arithmetic mean of the slice, or 0 for an empty slice.
func MeanFloat(xs []float64) float64 {
	if len(xs) == 0 {
		return 0
	}
	var s float64
	for _, x := range xs {
		s += x
	}
	return s / float64(len(xs))
}
