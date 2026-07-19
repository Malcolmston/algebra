package matroids

// GraphicMatroid is the cycle matroid (polygon matroid) of an undirected
// multigraph. The ground-set elements are the graph's edges, numbered in the
// order they were added. A set of edges is independent exactly when it contains
// no cycle, i.e. it forms a forest. The rank of an edge set S is
// (#vertices) - (#connected components of the spanning subgraph (V, S)).
type GraphicMatroid struct {
	verts int
	edges [][2]int
}

// NewGraphicMatroid returns an empty graphic matroid on the given number of
// vertices (labeled 0..verts-1) with no edges. It panics if verts < 0.
func NewGraphicMatroid(verts int) *GraphicMatroid {
	if verts < 0 {
		panic("matroids: negative vertex count")
	}
	return &GraphicMatroid{verts: verts}
}

// NewGraphicMatroidFromEdges returns the graphic matroid on verts vertices with
// the given edges. Each edge is a {u, v} pair with endpoints in [0, verts).
// Loops (u == v) and parallel edges are allowed. It panics on an out-of-range
// endpoint.
func NewGraphicMatroidFromEdges(verts int, edges [][2]int) *GraphicMatroid {
	m := NewGraphicMatroid(verts)
	for _, e := range edges {
		m.AddEdge(e[0], e[1])
	}
	return m
}

// AddEdge appends an edge {u, v} and returns its element index. It panics on an
// out-of-range endpoint.
func (m *GraphicMatroid) AddEdge(u, v int) int {
	if u < 0 || u >= m.verts || v < 0 || v >= m.verts {
		panic(ErrInvalidElement)
	}
	m.edges = append(m.edges, [2]int{u, v})
	return len(m.edges) - 1
}

// Size returns the number of edges (ground-set elements).
func (m *GraphicMatroid) Size() int { return len(m.edges) }

// NumVertices returns the number of vertices of the underlying graph.
func (m *GraphicMatroid) NumVertices() int { return m.verts }

// Edge returns the endpoints of edge e.
func (m *GraphicMatroid) Edge(e int) (int, int) { return m.edges[e][0], m.edges[e][1] }

// Rank returns verts - (number of connected components of the spanning
// subgraph induced by the edges in set). Equivalently it is the size of a
// spanning forest of that subgraph.
func (m *GraphicMatroid) Rank(set []int) int {
	uf := newUnionFind(m.verts)
	seen := make(map[int]bool, len(set))
	for _, e := range set {
		if e < 0 || e >= len(m.edges) || seen[e] {
			continue
		}
		seen[e] = true
		uf.union(m.edges[e][0], m.edges[e][1])
	}
	// rank = verts - components
	return m.verts - uf.count
}

// IsForest reports whether the given edge set is independent in m (contains no
// cycle).
func (m *GraphicMatroid) IsForest(set []int) bool { return Independent(m, set) }

// unionFind is a disjoint-set structure with union by rank and path
// compression, used to compute connectivity for graphic matroids.
type unionFind struct {
	parent []int
	rank   []int
	count  int
}

func newUnionFind(n int) *unionFind {
	p := make([]int, n)
	for i := range p {
		p[i] = i
	}
	return &unionFind{parent: p, rank: make([]int, n), count: n}
}

func (u *unionFind) find(x int) int {
	for u.parent[x] != x {
		u.parent[x] = u.parent[u.parent[x]]
		x = u.parent[x]
	}
	return x
}

func (u *unionFind) union(a, b int) bool {
	ra, rb := u.find(a), u.find(b)
	if ra == rb {
		return false
	}
	if u.rank[ra] < u.rank[rb] {
		ra, rb = rb, ra
	}
	u.parent[rb] = ra
	if u.rank[ra] == u.rank[rb] {
		u.rank[ra]++
	}
	u.count--
	return true
}

// CycleMatroidOfCompleteGraph returns the graphic matroid of the complete graph
// K_n. Its ground set has n(n-1)/2 edges and its rank is n-1 (for n ≥ 1).
func CycleMatroidOfCompleteGraph(n int) *GraphicMatroid {
	m := NewGraphicMatroid(n)
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			m.AddEdge(i, j)
		}
	}
	return m
}

// CycleMatroidOfCycleGraph returns the graphic matroid of the cycle graph C_n
// (n ≥ 1 vertices joined in a single ring). The whole edge set is its unique
// circuit for n ≥ 3.
func CycleMatroidOfCycleGraph(n int) *GraphicMatroid {
	m := NewGraphicMatroid(n)
	if n == 1 {
		m.AddEdge(0, 0)
		return m
	}
	for i := 0; i < n; i++ {
		m.AddEdge(i, (i+1)%n)
	}
	return m
}
