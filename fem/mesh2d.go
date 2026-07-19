package fem

import (
	"errors"
	"math"
	"sort"
)

// Edge is an undirected mesh edge identified by its two node indices with the
// smaller index stored first.
type Edge struct {
	A, B int
}

// MakeEdge returns the canonical Edge for nodes a and b (smaller index first).
func MakeEdge(a, b int) Edge {
	if a <= b {
		return Edge{a, b}
	}
	return Edge{b, a}
}

// Mesh2D is a two-dimensional triangular mesh: a list of node coordinates and a
// list of triangles given by triples of node indices.
type Mesh2D struct {
	Nodes     [][2]float64
	Triangles [][3]int
}

// NewMesh2D builds and validates a triangular mesh from nodes and triangles.
func NewMesh2D(nodes [][2]float64, triangles [][3]int) (*Mesh2D, error) {
	if len(nodes) < 3 {
		return nil, errors.New("fem: a 2D mesh needs at least three nodes")
	}
	for _, t := range triangles {
		for _, idx := range t {
			if idx < 0 || idx >= len(nodes) {
				return nil, errors.New("fem: triangle references an out-of-range node")
			}
		}
	}
	nc := make([][2]float64, len(nodes))
	copy(nc, nodes)
	tc := make([][3]int, len(triangles))
	copy(tc, triangles)
	return &Mesh2D{Nodes: nc, Triangles: tc}, nil
}

// NumNodes returns the number of nodes.
func (m *Mesh2D) NumNodes() int { return len(m.Nodes) }

// NumTriangles returns the number of triangles.
func (m *Mesh2D) NumTriangles() int { return len(m.Triangles) }

// TriangleVertices returns the three vertex coordinates of triangle t.
func (m *Mesh2D) TriangleVertices(t int) (v1, v2, v3 [2]float64) {
	tri := m.Triangles[t]
	return m.Nodes[tri[0]], m.Nodes[tri[1]], m.Nodes[tri[2]]
}

// TriangleArea returns the unsigned area of triangle t.
func (m *Mesh2D) TriangleArea(t int) float64 {
	v1, v2, v3 := m.TriangleVertices(t)
	return TriangleArea(v1, v2, v3)
}

// TotalArea returns the sum of all triangle areas (the meshed domain area).
func (m *Mesh2D) TotalArea() float64 {
	var s float64
	for t := 0; t < m.NumTriangles(); t++ {
		s += m.TriangleArea(t)
	}
	return s
}

// Centroid returns the centroid coordinate of triangle t.
func (m *Mesh2D) Centroid(t int) [2]float64 {
	v1, v2, v3 := m.TriangleVertices(t)
	return [2]float64{(v1[0] + v2[0] + v3[0]) / 3, (v1[1] + v2[1] + v3[1]) / 3}
}

// edgeCount returns a map from each edge to the number of triangles containing
// it. Boundary edges appear exactly once.
func (m *Mesh2D) edgeCount() map[Edge]int {
	counts := make(map[Edge]int)
	for _, tri := range m.Triangles {
		counts[MakeEdge(tri[0], tri[1])]++
		counts[MakeEdge(tri[1], tri[2])]++
		counts[MakeEdge(tri[2], tri[0])]++
	}
	return counts
}

// BoundaryEdges returns the edges that belong to exactly one triangle, sorted
// for reproducibility.
func (m *Mesh2D) BoundaryEdges() []Edge {
	counts := m.edgeCount()
	var edges []Edge
	for e, c := range counts {
		if c == 1 {
			edges = append(edges, e)
		}
	}
	sort.Slice(edges, func(i, j int) bool {
		if edges[i].A != edges[j].A {
			return edges[i].A < edges[j].A
		}
		return edges[i].B < edges[j].B
	})
	return edges
}

// BoundaryNodes returns the sorted indices of nodes lying on the boundary.
func (m *Mesh2D) BoundaryNodes() []int {
	set := make(map[int]bool)
	for _, e := range m.BoundaryEdges() {
		set[e.A] = true
		set[e.B] = true
	}
	out := make([]int, 0, len(set))
	for n := range set {
		out = append(out, n)
	}
	sort.Ints(out)
	return out
}

// InteriorNodes returns the sorted indices of nodes not on the boundary.
func (m *Mesh2D) InteriorNodes() []int {
	onB := make(map[int]bool)
	for _, n := range m.BoundaryNodes() {
		onB[n] = true
	}
	out := make([]int, 0)
	for n := 0; n < m.NumNodes(); n++ {
		if !onB[n] {
			out = append(out, n)
		}
	}
	return out
}

// AllEdges returns every distinct edge of the mesh, sorted.
func (m *Mesh2D) AllEdges() []Edge {
	counts := m.edgeCount()
	edges := make([]Edge, 0, len(counts))
	for e := range counts {
		edges = append(edges, e)
	}
	sort.Slice(edges, func(i, j int) bool {
		if edges[i].A != edges[j].A {
			return edges[i].A < edges[j].A
		}
		return edges[i].B < edges[j].B
	})
	return edges
}

// EnsureCounterClockwise returns a copy of the mesh with every triangle
// re-ordered to have counter-clockwise (positive) orientation.
func (m *Mesh2D) EnsureCounterClockwise() *Mesh2D {
	tris := make([][3]int, len(m.Triangles))
	copy(tris, m.Triangles)
	for i, tri := range tris {
		v1, v2, v3 := m.Nodes[tri[0]], m.Nodes[tri[1]], m.Nodes[tri[2]]
		if TriangleSignedArea(v1, v2, v3) < 0 {
			tris[i] = [3]int{tri[0], tri[2], tri[1]}
		}
	}
	nc := make([][2]float64, len(m.Nodes))
	copy(nc, m.Nodes)
	return &Mesh2D{Nodes: nc, Triangles: tris}
}

// RectangleMesh returns a structured triangular mesh of the rectangle
// [x0,x1]×[y0,y1] using nx by ny cells, each split into two triangles. It
// panics if nx < 1, ny < 1, x1 <= x0 or y1 <= y0.
func RectangleMesh(x0, y0, x1, y1 float64, nx, ny int) *Mesh2D {
	if nx < 1 || ny < 1 {
		panic("fem: RectangleMesh requires nx, ny >= 1")
	}
	if x1 <= x0 || y1 <= y0 {
		panic("fem: RectangleMesh requires x1 > x0 and y1 > y0")
	}
	nodes := make([][2]float64, 0, (nx+1)*(ny+1))
	hx := (x1 - x0) / float64(nx)
	hy := (y1 - y0) / float64(ny)
	idx := func(i, j int) int { return j*(nx+1) + i }
	for j := 0; j <= ny; j++ {
		for i := 0; i <= nx; i++ {
			nodes = append(nodes, [2]float64{x0 + float64(i)*hx, y0 + float64(j)*hy})
		}
	}
	tris := make([][3]int, 0, 2*nx*ny)
	for j := 0; j < ny; j++ {
		for i := 0; i < nx; i++ {
			a := idx(i, j)
			b := idx(i+1, j)
			c := idx(i+1, j+1)
			d := idx(i, j+1)
			tris = append(tris, [3]int{a, b, c}, [3]int{a, c, d})
		}
	}
	return &Mesh2D{Nodes: nodes, Triangles: tris}
}

// UnitSquareMesh returns a structured mesh of the unit square with n cells per
// side.
func UnitSquareMesh(n int) *Mesh2D {
	return RectangleMesh(0, 0, 1, 1, n, n)
}

// Refine returns the uniformly (red) refined mesh: each triangle is split into
// four by connecting its edge midpoints. Shared midpoints are reused.
func (m *Mesh2D) Refine() *Mesh2D {
	nodes := make([][2]float64, len(m.Nodes))
	copy(nodes, m.Nodes)
	midIndex := make(map[Edge]int)
	getMid := func(a, b int) int {
		e := MakeEdge(a, b)
		if idx, ok := midIndex[e]; ok {
			return idx
		}
		pa, pb := m.Nodes[a], m.Nodes[b]
		mid := [2]float64{0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1])}
		idx := len(nodes)
		nodes = append(nodes, mid)
		midIndex[e] = idx
		return idx
	}
	tris := make([][3]int, 0, 4*len(m.Triangles))
	for _, tri := range m.Triangles {
		a, b, c := tri[0], tri[1], tri[2]
		ab := getMid(a, b)
		bc := getMid(b, c)
		ca := getMid(c, a)
		tris = append(tris,
			[3]int{a, ab, ca},
			[3]int{ab, b, bc},
			[3]int{ca, bc, c},
			[3]int{ab, bc, ca},
		)
	}
	return &Mesh2D{Nodes: nodes, Triangles: tris}
}

// RefineN applies Refine k times.
func (m *Mesh2D) RefineN(k int) *Mesh2D {
	out := m
	for i := 0; i < k; i++ {
		out = out.Refine()
	}
	return out
}

// MaxDiameter returns the largest triangle diameter (longest edge length), a
// common mesh-size measure.
func (m *Mesh2D) MaxDiameter() float64 {
	var max float64
	for t := 0; t < m.NumTriangles(); t++ {
		v1, v2, v3 := m.TriangleVertices(t)
		for _, d := range []float64{dist2(v1, v2), dist2(v2, v3), dist2(v3, v1)} {
			if d > max {
				max = d
			}
		}
	}
	return math.Sqrt(max)
}

// P2Mesh converts the P1 triangular mesh into P2 data: coordinates of all P2
// degrees of freedom (vertices followed by edge midpoints) and, for each
// triangle, the six global P2 node indices in the ordering used by
// ShapeP2Triangle.
func (m *Mesh2D) P2Mesh() (nodes [][2]float64, conn [][6]int) {
	nodes = make([][2]float64, len(m.Nodes))
	copy(nodes, m.Nodes)
	midIndex := make(map[Edge]int)
	getMid := func(a, b int) int {
		e := MakeEdge(a, b)
		if idx, ok := midIndex[e]; ok {
			return idx
		}
		pa, pb := m.Nodes[a], m.Nodes[b]
		mid := [2]float64{0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1])}
		idx := len(nodes)
		nodes = append(nodes, mid)
		midIndex[e] = idx
		return idx
	}
	conn = make([][6]int, len(m.Triangles))
	for t, tri := range m.Triangles {
		a, b, c := tri[0], tri[1], tri[2]
		m23 := getMid(b, c) // opposite v1
		m13 := getMid(a, c) // opposite v2
		m12 := getMid(a, b) // opposite v3
		conn[t] = [6]int{a, b, c, m23, m13, m12}
	}
	return nodes, conn
}

func dist2(a, b [2]float64) float64 {
	dx := a[0] - b[0]
	dy := a[1] - b[1]
	return dx*dx + dy*dy
}
