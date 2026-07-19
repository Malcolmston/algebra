package graphspectral

import "math"

// Graph is a finite, undirected, optionally weighted graph on the vertex set
// {0, 1, …, n-1}. It is stored as a symmetric adjacency matrix; an edge weight
// of zero means "no edge". Self-loops (nonzero diagonal entries) are permitted
// and are reflected in the degree.
type Graph struct {
	n   int
	adj *Matrix
}

// NewGraph returns an edgeless graph on n vertices. It panics if n is negative.
func NewGraph(n int) *Graph {
	if n < 0 {
		panic("graphspectral: negative vertex count")
	}
	return &Graph{n: n, adj: NewMatrix(n, n)}
}

// NewGraphFromAdjacency builds a graph from a square symmetric adjacency matrix.
// The matrix is copied and symmetrized defensively. It returns ErrNotSquare or
// ErrNotSymmetric on an invalid matrix.
func NewGraphFromAdjacency(a *Matrix) (*Graph, error) {
	if !a.IsSquare() {
		return nil, ErrNotSquare
	}
	if !a.IsSymmetric(1e-9) {
		return nil, ErrNotSymmetric
	}
	return &Graph{n: a.rows, adj: a.Symmetrize()}, nil
}

// Order returns the number of vertices.
func (g *Graph) Order() int { return g.n }

// AddEdge adds (or, for a weighted graph, sets to weight 1) an undirected edge
// between u and v. It panics on an out-of-range vertex.
func (g *Graph) AddEdge(u, v int) { g.SetWeight(u, v, 1) }

// SetWeight sets the weight of the undirected edge between u and v (and, by
// symmetry, between v and u). A weight of zero removes the edge. It panics on an
// out-of-range vertex.
func (g *Graph) SetWeight(u, v int, w float64) {
	g.check(u)
	g.check(v)
	g.adj.Set(u, v, w)
	if u != v {
		g.adj.Set(v, u, w)
	}
}

// Weight returns the weight of the edge between u and v (0 if there is none).
func (g *Graph) Weight(u, v int) float64 { return g.adj.At(u, v) }

// HasEdge reports whether there is a (nonzero-weight) edge between u and v.
func (g *Graph) HasEdge(u, v int) bool { return g.adj.At(u, v) != 0 }

// RemoveEdge deletes the edge between u and v.
func (g *Graph) RemoveEdge(u, v int) { g.SetWeight(u, v, 0) }

func (g *Graph) check(v int) {
	if v < 0 || v >= g.n {
		panic("graphspectral: vertex index out of range")
	}
}

// Adjacency returns a copy of the adjacency matrix A, where A[i][j] is the
// weight of the edge between i and j.
func (g *Graph) Adjacency() *Matrix { return g.adj.Clone() }

// WeightedDegree returns the weighted degree of vertex v: the sum of the weights
// of its incident edges. A self-loop of weight w contributes w.
func (g *Graph) WeightedDegree(v int) float64 {
	var s float64
	for j := 0; j < g.n; j++ {
		s += g.adj.At(v, j)
	}
	return s
}

// Degree returns the (integer) number of neighbours of v in an unweighted graph,
// counting a self-loop once. For weighted degree use [Graph.WeightedDegree].
func (g *Graph) Degree(v int) int {
	d := 0
	for j := 0; j < g.n; j++ {
		if g.adj.At(v, j) != 0 {
			d++
		}
	}
	return d
}

// Degrees returns the weighted degree of every vertex.
func (g *Graph) Degrees() []float64 {
	out := make([]float64, g.n)
	for v := 0; v < g.n; v++ {
		out[v] = g.WeightedDegree(v)
	}
	return out
}

// DegreeMatrix returns the diagonal degree matrix D whose (v, v) entry is the
// weighted degree of v.
func (g *Graph) DegreeMatrix() *Matrix { return DiagMatrix(g.Degrees()) }

// Neighbors returns the vertices adjacent to v, in ascending order.
func (g *Graph) Neighbors(v int) []int {
	var out []int
	for j := 0; j < g.n; j++ {
		if j != v && g.adj.At(v, j) != 0 {
			out = append(out, j)
		}
	}
	return out
}

// Size returns the number of edges. Each undirected edge is counted once; a
// self-loop counts as one edge.
func (g *Graph) Size() int {
	c := 0
	for i := 0; i < g.n; i++ {
		for j := i; j < g.n; j++ {
			if g.adj.At(i, j) != 0 {
				c++
			}
		}
	}
	return c
}

// TotalWeight returns the sum of all edge weights, counting each undirected edge
// once.
func (g *Graph) TotalWeight() float64 {
	var s float64
	for i := 0; i < g.n; i++ {
		for j := i; j < g.n; j++ {
			s += g.adj.At(i, j)
		}
	}
	return s
}

// Laplacian returns the combinatorial Laplacian L = D - A.
func (g *Graph) Laplacian() *Matrix {
	l := g.adj.Scale(-1)
	for i := 0; i < g.n; i++ {
		l.Set(i, i, g.WeightedDegree(i)+l.At(i, i))
	}
	return l
}

// SignlessLaplacian returns the signless Laplacian Q = D + A.
func (g *Graph) SignlessLaplacian() *Matrix {
	q := g.adj.Clone()
	for i := 0; i < g.n; i++ {
		q.Set(i, i, g.WeightedDegree(i)+q.At(i, i))
	}
	return q
}

// NormalizedLaplacian returns the symmetric normalized Laplacian
// Lsym = I - D^{-1/2} A D^{-1/2}. Isolated vertices (degree 0) contribute a unit
// diagonal entry and zero off-diagonal entries.
func (g *Graph) NormalizedLaplacian() *Matrix {
	deg := g.Degrees()
	inv := make([]float64, g.n)
	for i, d := range deg {
		if d > 0 {
			inv[i] = 1 / math.Sqrt(d)
		}
	}
	l := NewMatrix(g.n, g.n)
	for i := 0; i < g.n; i++ {
		for j := 0; j < g.n; j++ {
			v := -inv[i] * g.adj.At(i, j) * inv[j]
			if i == j {
				if deg[i] > 0 {
					v += 1
				}
			}
			l.Set(i, j, v)
		}
	}
	return l
}

// RandomWalkLaplacian returns the random-walk Laplacian Lrw = I - D^{-1} A.
// Isolated vertices contribute a unit diagonal entry. Note Lrw is generally not
// symmetric.
func (g *Graph) RandomWalkLaplacian() *Matrix {
	deg := g.Degrees()
	l := NewMatrix(g.n, g.n)
	for i := 0; i < g.n; i++ {
		for j := 0; j < g.n; j++ {
			var v float64
			if deg[i] > 0 {
				v = -g.adj.At(i, j) / deg[i]
			}
			if i == j {
				if deg[i] > 0 {
					v += 1
				}
			}
			l.Set(i, j, v)
		}
	}
	return l
}

// TransitionMatrix returns the random-walk transition matrix P = D^{-1} A, whose
// row i is the probability distribution of a single step from vertex i. Rows for
// isolated vertices are all zero.
func (g *Graph) TransitionMatrix() *Matrix {
	deg := g.Degrees()
	p := NewMatrix(g.n, g.n)
	for i := 0; i < g.n; i++ {
		if deg[i] == 0 {
			continue
		}
		for j := 0; j < g.n; j++ {
			p.Set(i, j, g.adj.At(i, j)/deg[i])
		}
	}
	return p
}

// IncidenceMatrix returns the (oriented) vertex-edge incidence matrix B of an
// unweighted view of the graph. B has one column per edge (in ascending
// (u, v) order with u < v); the column for edge (u, v) has -1 in row u and +1 in
// row v. Self-loops are skipped. The identity B·Bᵀ = L holds for the unweighted
// Laplacian.
func (g *Graph) IncidenceMatrix() *Matrix {
	type edge struct{ u, v int }
	var edges []edge
	for i := 0; i < g.n; i++ {
		for j := i + 1; j < g.n; j++ {
			if g.adj.At(i, j) != 0 {
				edges = append(edges, edge{i, j})
			}
		}
	}
	b := NewMatrix(g.n, len(edges))
	for c, e := range edges {
		b.Set(e.u, c, -1)
		b.Set(e.v, c, 1)
	}
	return b
}

// IsRegular reports whether every vertex has the same weighted degree (to within
// 1e-9), and returns that common degree.
func (g *Graph) IsRegular() (bool, float64) {
	if g.n == 0 {
		return true, 0
	}
	d0 := g.WeightedDegree(0)
	for v := 1; v < g.n; v++ {
		if math.Abs(g.WeightedDegree(v)-d0) > 1e-9 {
			return false, 0
		}
	}
	return true, d0
}

// MaxDegree returns the largest weighted degree, or 0 for an empty graph.
func (g *Graph) MaxDegree() float64 { return VecMax(g.Degrees()) }

// MinDegree returns the smallest weighted degree, or 0 for an empty graph.
func (g *Graph) MinDegree() float64 { return VecMin(g.Degrees()) }

// Density returns the edge density of an unweighted view: the number of edges
// (ignoring self-loops) divided by the number of possible edges n(n-1)/2. It
// returns 0 for graphs with fewer than two vertices.
func (g *Graph) Density() float64 {
	if g.n < 2 {
		return 0
	}
	edges := 0
	for i := 0; i < g.n; i++ {
		for j := i + 1; j < g.n; j++ {
			if g.adj.At(i, j) != 0 {
				edges++
			}
		}
	}
	possible := float64(g.n*(g.n-1)) / 2
	return float64(edges) / possible
}

// Complement returns the complement of an unweighted simple graph: vertices are
// adjacent in the result exactly when they are non-adjacent (and distinct) in g.
func (g *Graph) Complement() *Graph {
	c := NewGraph(g.n)
	for i := 0; i < g.n; i++ {
		for j := i + 1; j < g.n; j++ {
			if g.adj.At(i, j) == 0 {
				c.SetWeight(i, j, 1)
			}
		}
	}
	return c
}

// ConnectedComponents returns a labelling of the vertices into connected
// components, using a breadth-first flood fill over the (unweighted) edge set.
// Component labels are consecutive integers starting at 0, assigned in order of
// first discovery.
func (g *Graph) ConnectedComponents() []int {
	label := make([]int, g.n)
	for i := range label {
		label[i] = -1
	}
	next := 0
	for s := 0; s < g.n; s++ {
		if label[s] != -1 {
			continue
		}
		queue := []int{s}
		label[s] = next
		for len(queue) > 0 {
			v := queue[0]
			queue = queue[1:]
			for w := 0; w < g.n; w++ {
				if w != v && g.adj.At(v, w) != 0 && label[w] == -1 {
					label[w] = next
					queue = append(queue, w)
				}
			}
		}
		next++
	}
	return label
}

// NumConnectedComponents returns the number of connected components. An empty
// graph has zero components.
func (g *Graph) NumConnectedComponents() int {
	if g.n == 0 {
		return 0
	}
	label := g.ConnectedComponents()
	return VecMaxInt(label) + 1
}

// IsConnected reports whether the graph has exactly one connected component.
func (g *Graph) IsConnected() bool { return g.NumConnectedComponents() == 1 }

// IsBipartite reports whether the graph is bipartite, and if so returns a
// two-colouring (side 0 or 1) of its vertices. Isolated vertices are placed on
// side 0. When the graph is not bipartite the returned slice is nil.
func (g *Graph) IsBipartite() (bool, []int) {
	color := make([]int, g.n)
	for i := range color {
		color[i] = -1
	}
	for s := 0; s < g.n; s++ {
		if color[s] != -1 {
			continue
		}
		color[s] = 0
		queue := []int{s}
		for len(queue) > 0 {
			v := queue[0]
			queue = queue[1:]
			for w := 0; w < g.n; w++ {
				if w == v || g.adj.At(v, w) == 0 {
					continue
				}
				if color[w] == -1 {
					color[w] = 1 - color[v]
					queue = append(queue, w)
				} else if color[w] == color[v] {
					return false, nil
				}
			}
		}
	}
	return true, color
}

// VecMaxInt returns the maximum entry of an int slice, or -1 for an empty slice.
func VecMaxInt(v []int) int {
	if len(v) == 0 {
		return -1
	}
	m := v[0]
	for _, x := range v[1:] {
		if x > m {
			m = x
		}
	}
	return m
}
