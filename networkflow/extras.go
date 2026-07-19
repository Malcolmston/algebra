package networkflow

import "sort"

// FlowNetworkFromMatrix builds a directed [FlowNetwork] from a square capacity
// matrix: an edge i->j is added for every strictly positive cap[i][j]. It
// returns [ErrDimensionMismatch] if the matrix is not square.
func FlowNetworkFromMatrix(cap [][]int64) (*FlowNetwork, error) {
	n := len(cap)
	for _, row := range cap {
		if len(row) != n {
			return nil, ErrDimensionMismatch
		}
	}
	g := NewFlowNetwork(n)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if cap[i][j] > 0 {
				g.AddEdge(i, j, cap[i][j])
			}
		}
	}
	return g, nil
}

// CapacityMatrix returns the n-by-n matrix whose (i,j) entry is the total
// capacity of caller edges from i to j.
func (g *FlowNetwork) CapacityMatrix() [][]int64 {
	m := make([][]int64, g.n)
	for i := range m {
		m[i] = make([]int64, g.n)
	}
	for i := 0; i < len(g.edges); i += 2 {
		e := g.edges[i]
		m[e.from][e.to] += e.cap
	}
	return m
}

// HasEdge reports whether at least one caller edge runs from u to v.
func (g *FlowNetwork) HasEdge(u, v int) bool {
	for i := 0; i < len(g.edges); i += 2 {
		if g.edges[i].from == u && g.edges[i].to == v {
			return true
		}
	}
	return false
}

// Reachable returns the vertices reachable from src along positive-capacity
// caller edges (ignoring flow), sorted ascending.
func (g *FlowNetwork) Reachable(src int) []int {
	seen := make([]bool, g.n)
	if !g.validVertex(src) {
		return nil
	}
	seen[src] = true
	queue := []int{src}
	for len(queue) > 0 {
		v := queue[0]
		queue = queue[1:]
		for _, id := range g.adj[v] {
			e := g.edges[id]
			if id%2 == 0 && e.cap > 0 && !seen[e.to] {
				seen[e.to] = true
				queue = append(queue, e.to)
			}
		}
	}
	var out []int
	for v, ok := range seen {
		if ok {
			out = append(out, v)
		}
	}
	return out
}

// MinCut is an alias for [MinCutST].
func MinCut(g *FlowNetwork, s, t int) *MinCutResult { return MinCutST(g, s, t) }

// WeightedGraphFromMatrix builds an undirected [WeightedGraph] from a symmetric
// weight matrix (only the upper triangle is read). It returns
// [ErrDimensionMismatch] if the matrix is not square.
func WeightedGraphFromMatrix(w [][]float64) (*WeightedGraph, error) {
	n := len(w)
	for _, row := range w {
		if len(row) != n {
			return nil, ErrDimensionMismatch
		}
	}
	g := NewWeightedGraph(n)
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			if w[i][j] != 0 {
				g.SetWeight(i, j, w[i][j])
			}
		}
	}
	return g, nil
}

// Edges returns every undirected edge with non-zero weight as {u, v, weight},
// with u < v.
func (g *WeightedGraph) Edges() [][3]float64 {
	var out [][3]float64
	for i := 0; i < g.n; i++ {
		for j := i + 1; j < g.n; j++ {
			if g.w[i][j] != 0 {
				out = append(out, [3]float64{float64(i), float64(j), g.w[i][j]})
			}
		}
	}
	return out
}

// EdgeConnectivity returns the weighted edge connectivity of an undirected
// graph, i.e. the value of its global minimum cut ([StoerWagner]). For a graph
// with unit weights this is the classical edge connectivity.
func EdgeConnectivity(g *WeightedGraph) float64 { return StoerWagner(g).Value }

// SVertices returns the S side of the cut.
func (c *MinCutResult) SVertices() []int { return c.S }

// TVertices returns the T side of the cut.
func (c *MinCutResult) TVertices() []int { return c.T }

// BipartiteFromMatrix builds a [BipartiteGraph] from a boolean adjacency matrix
// where adj[i][j] true means left i is joined to right j. It returns
// [ErrDimensionMismatch] if the rows are ragged.
func BipartiteFromMatrix(adj [][]bool) (*BipartiteGraph, error) {
	nl := len(adj)
	nr := 0
	if nl > 0 {
		nr = len(adj[0])
	}
	for _, row := range adj {
		if len(row) != nr {
			return nil, ErrDimensionMismatch
		}
	}
	g := NewBipartiteGraph(nl, nr)
	for i := 0; i < nl; i++ {
		for j := 0; j < nr; j++ {
			if adj[i][j] {
				g.AddEdge(i, j)
			}
		}
	}
	return g, nil
}

// Edges returns every edge of the bipartite graph as {left, right} pairs,
// sorted lexicographically.
func (g *BipartiteGraph) Edges() [][2]int {
	var out [][2]int
	for u := 0; u < g.nl; u++ {
		for _, v := range g.adj[u] {
			out = append(out, [2]int{u, v})
		}
	}
	sort.Slice(out, func(a, b int) bool {
		if out[a][0] != out[b][0] {
			return out[a][0] < out[b][0]
		}
		return out[a][1] < out[b][1]
	})
	return out
}

// RightNeighbors returns the left vertices adjacent to right vertex v, sorted.
func (g *BipartiteGraph) RightNeighbors(v int) []int {
	var out []int
	for u := 0; u < g.nl; u++ {
		if g.set[u][v] {
			out = append(out, u)
		}
	}
	return out
}

// UnmatchedLeft returns the left vertices left unmatched.
func (m *MatchingResult) UnmatchedLeft() []int {
	var out []int
	for u, v := range m.MatchL {
		if v < 0 {
			out = append(out, u)
		}
	}
	return out
}

// UnmatchedRight returns the right vertices left unmatched.
func (m *MatchingResult) UnmatchedRight() []int {
	var out []int
	for v, u := range m.MatchR {
		if u < 0 {
			out = append(out, v)
		}
	}
	return out
}

// Size returns the number of assigned rows.
func (a *AssignmentResult) Size() int {
	n := 0
	for _, j := range a.Assignment {
		if j >= 0 {
			n++
		}
	}
	return n
}

// OutDegree returns the number of caller edges leaving v.
func (g *MinCostNetwork) OutDegree(v int) int {
	c := 0
	for _, id := range g.adj[v] {
		if id%2 == 0 {
			c++
		}
	}
	return c
}

// InDegree returns the number of caller edges entering v.
func (g *MinCostNetwork) InDegree(v int) int {
	c := 0
	for _, id := range g.adj[v] {
		if id%2 == 1 {
			c++
		}
	}
	return c
}

// TotalCapacity returns the sum of all caller-edge capacities.
func (g *MinCostNetwork) TotalCapacity() int64 {
	var c int64
	for i := 0; i < len(g.edges); i += 2 {
		c += g.edges[i].cap
	}
	return c
}

// Root returns the root vertex of the Gomory-Hu tree, which is always 0.
func (t *GomoryHuTree) Root() int { return 0 }

// PathFlowTotal sums the flow over the given decomposition paths that are not
// cycles; for a valid decomposition this equals the flow value.
func PathFlowTotal(paths []FlowPath) int64 {
	var total int64
	for _, p := range paths {
		if !p.Cycle {
			total += p.Flow
		}
	}
	return total
}
