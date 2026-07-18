package graph

import (
	"errors"
	"sort"
)

// ErrVertexNotFound reports that an operation referenced a vertex that is not
// present in the graph.
var ErrVertexNotFound = errors.New("graph: vertex not found")

// ErrNegativeWeight reports that an algorithm requiring non-negative edge
// weights (such as [Graph.Dijkstra] or [Graph.AStar]) encountered a negative
// weight.
var ErrNegativeWeight = errors.New("graph: negative edge weight")

// ErrNegativeCycle reports that a shortest-path computation detected a cycle of
// negative total weight reachable from the source, so finite shortest paths do
// not exist.
var ErrNegativeCycle = errors.New("graph: negative-weight cycle")

// ErrNotDAG reports that an operation requiring a directed acyclic graph was
// applied to a graph that is undirected or contains a cycle.
var ErrNotDAG = errors.New("graph: not a directed acyclic graph")

// Edge is a weighted connection from vertex From to vertex To. For undirected
// graphs each edge is reported once with From <= To by [Graph.Edges].
type Edge struct {
	// From is the tail (source) vertex of the edge.
	From int
	// To is the head (destination) vertex of the edge.
	To int
	// Weight is the numeric weight of the edge; unweighted edges use 1.
	Weight float64
}

// Graph is an adjacency-list graph with integer vertices and float64 edge
// weights. It may be directed or undirected. The zero value is not usable;
// construct graphs with [New] or [NewDirected].
type Graph struct {
	directed bool
	adj      map[int]map[int]float64
}

// New returns an empty undirected graph.
func New() *Graph {
	return &Graph{directed: false, adj: make(map[int]map[int]float64)}
}

// NewDirected returns an empty directed graph.
func NewDirected() *Graph {
	return &Graph{directed: true, adj: make(map[int]map[int]float64)}
}

// Directed reports whether the graph is directed.
func (g *Graph) Directed() bool { return g.directed }

// AddVertex adds v to the graph if it is not already present. Adding an
// existing vertex is a no-op.
func (g *Graph) AddVertex(v int) {
	if _, ok := g.adj[v]; !ok {
		g.adj[v] = make(map[int]float64)
	}
}

// HasVertex reports whether v is present in the graph.
func (g *Graph) HasVertex(v int) bool {
	_, ok := g.adj[v]
	return ok
}

// AddEdge adds an edge between u and v with weight 1, creating either vertex if
// necessary. For an undirected graph the reverse edge is added as well. Adding
// an existing edge overwrites its weight with 1.
func (g *Graph) AddEdge(u, v int) {
	g.AddWeightedEdge(u, v, 1)
}

// AddWeightedEdge adds an edge from u to v with the given weight, creating
// either vertex if necessary. For an undirected graph the reverse edge is added
// with the same weight. Adding an existing edge overwrites its weight.
func (g *Graph) AddWeightedEdge(u, v int, w float64) {
	g.AddVertex(u)
	g.AddVertex(v)
	g.adj[u][v] = w
	if !g.directed {
		g.adj[v][u] = w
	}
}

// HasEdge reports whether an edge from u to v exists. For undirected graphs the
// orientation is irrelevant.
func (g *Graph) HasEdge(u, v int) bool {
	nbrs, ok := g.adj[u]
	if !ok {
		return false
	}
	_, ok = nbrs[v]
	return ok
}

// Weight returns the weight of the edge from u to v and whether it exists.
func (g *Graph) Weight(u, v int) (float64, bool) {
	nbrs, ok := g.adj[u]
	if !ok {
		return 0, false
	}
	w, ok := nbrs[v]
	return w, ok
}

// RemoveEdge removes the edge from u to v if present. For an undirected graph
// the reverse edge is removed as well. Removing an absent edge is a no-op.
func (g *Graph) RemoveEdge(u, v int) {
	if nbrs, ok := g.adj[u]; ok {
		delete(nbrs, v)
	}
	if !g.directed {
		if nbrs, ok := g.adj[v]; ok {
			delete(nbrs, u)
		}
	}
}

// RemoveVertex removes v and all edges incident to it. Removing an absent
// vertex is a no-op.
func (g *Graph) RemoveVertex(v int) {
	if _, ok := g.adj[v]; !ok {
		return
	}
	delete(g.adj, v)
	for _, nbrs := range g.adj {
		delete(nbrs, v)
	}
}

// Vertices returns all vertices in ascending order.
func (g *Graph) Vertices() []int {
	vs := make([]int, 0, len(g.adj))
	for v := range g.adj {
		vs = append(vs, v)
	}
	sort.Ints(vs)
	return vs
}

// Neighbors returns the out-neighbors of v in ascending order. It returns
// [ErrVertexNotFound] if v is absent.
func (g *Graph) Neighbors(v int) ([]int, error) {
	nbrs, ok := g.adj[v]
	if !ok {
		return nil, ErrVertexNotFound
	}
	out := make([]int, 0, len(nbrs))
	for n := range nbrs {
		out = append(out, n)
	}
	sort.Ints(out)
	return out, nil
}

// Edges returns all edges of the graph. For a directed graph every edge is
// returned. For an undirected graph each edge is returned once, oriented so
// that From <= To. Edges are sorted by (From, To).
func (g *Graph) Edges() []Edge {
	var edges []Edge
	for u, nbrs := range g.adj {
		for v, w := range nbrs {
			if g.directed || u <= v {
				edges = append(edges, Edge{From: u, To: v, Weight: w})
			}
		}
	}
	sort.Slice(edges, func(i, j int) bool {
		if edges[i].From != edges[j].From {
			return edges[i].From < edges[j].From
		}
		return edges[i].To < edges[j].To
	})
	return edges
}

// NumVertices returns the number of vertices (the order of the graph).
func (g *Graph) NumVertices() int { return len(g.adj) }

// NumEdges returns the number of edges (the size of the graph). Each undirected
// edge is counted once.
func (g *Graph) NumEdges() int {
	count := 0
	selfLoops := 0
	for u, nbrs := range g.adj {
		for v := range nbrs {
			count++
			if u == v {
				selfLoops++
			}
		}
	}
	if g.directed {
		return count
	}
	// Each non-loop undirected edge was counted twice; self-loops once.
	return (count-selfLoops)/2 + selfLoops
}

// OutDegree returns the number of edges leaving v. It returns
// [ErrVertexNotFound] if v is absent.
func (g *Graph) OutDegree(v int) (int, error) {
	nbrs, ok := g.adj[v]
	if !ok {
		return 0, ErrVertexNotFound
	}
	return len(nbrs), nil
}

// InDegree returns the number of edges entering v. For an undirected graph it
// equals the out-degree. It returns [ErrVertexNotFound] if v is absent.
func (g *Graph) InDegree(v int) (int, error) {
	if _, ok := g.adj[v]; !ok {
		return 0, ErrVertexNotFound
	}
	count := 0
	for _, nbrs := range g.adj {
		if _, ok := nbrs[v]; ok {
			count++
		}
	}
	return count, nil
}

// Degree returns the degree of v: for an undirected graph the number of
// incident edges (a self-loop counts twice), and for a directed graph the sum
// of in-degree and out-degree. It returns [ErrVertexNotFound] if v is absent.
func (g *Graph) Degree(v int) (int, error) {
	nbrs, ok := g.adj[v]
	if !ok {
		return 0, ErrVertexNotFound
	}
	if !g.directed {
		d := len(nbrs)
		if _, loop := nbrs[v]; loop {
			d++ // self-loop contributes 2
		}
		return d, nil
	}
	in, _ := g.InDegree(v)
	return in + len(nbrs), nil
}

// DegreeSequence returns the degrees of all vertices sorted in non-increasing
// order. The definition of degree matches [Graph.Degree].
func (g *Graph) DegreeSequence() []int {
	vs := g.Vertices()
	seq := make([]int, 0, len(vs))
	for _, v := range vs {
		d, _ := g.Degree(v)
		seq = append(seq, d)
	}
	sort.Sort(sort.Reverse(sort.IntSlice(seq)))
	return seq
}

// Copy returns a deep copy of the graph.
func (g *Graph) Copy() *Graph {
	c := &Graph{directed: g.directed, adj: make(map[int]map[int]float64, len(g.adj))}
	for u, nbrs := range g.adj {
		m := make(map[int]float64, len(nbrs))
		for v, w := range nbrs {
			m[v] = w
		}
		c.adj[u] = m
	}
	return c
}

// Reverse returns a new directed graph with every edge reversed. For an
// undirected graph it returns a copy, since undirected edges are symmetric.
func (g *Graph) Reverse() *Graph {
	if !g.directed {
		return g.Copy()
	}
	r := NewDirected()
	for _, v := range g.Vertices() {
		r.AddVertex(v)
	}
	for u, nbrs := range g.adj {
		for v, w := range nbrs {
			r.AddWeightedEdge(v, u, w)
		}
	}
	return r
}

// AdjacencyMatrix returns the weighted adjacency matrix of the graph together
// with the vertex ordering (ascending) that indexes its rows and columns. Entry
// [i][j] is the weight of the edge from vertex order[i] to order[j], or 0 if no
// such edge exists.
func (g *Graph) AdjacencyMatrix() (matrix [][]float64, order []int) {
	order = g.Vertices()
	index := make(map[int]int, len(order))
	for i, v := range order {
		index[v] = i
	}
	n := len(order)
	matrix = make([][]float64, n)
	for i := range matrix {
		matrix[i] = make([]float64, n)
	}
	for u, nbrs := range g.adj {
		for v, w := range nbrs {
			matrix[index[u]][index[v]] = w
		}
	}
	return matrix, order
}

// graphSortedKeys returns the sorted keys of an int-keyed map.
func graphSortedKeys[V any](m map[int]V) []int {
	ks := make([]int, 0, len(m))
	for k := range m {
		ks = append(ks, k)
	}
	sort.Ints(ks)
	return ks
}
