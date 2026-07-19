package graphspectral

// EmptyGraph returns the edgeless graph on n vertices.
func EmptyGraph(n int) *Graph { return NewGraph(n) }

// PathGraph returns the path P_n: vertices 0..n-1 joined in a line
// 0-1-2-…-(n-1).
func PathGraph(n int) *Graph {
	g := NewGraph(n)
	for i := 0; i+1 < n; i++ {
		g.AddEdge(i, i+1)
	}
	return g
}

// CycleGraph returns the cycle C_n on n vertices (n >= 3 for a genuine cycle;
// for n < 3 it degenerates to a path or single edge).
func CycleGraph(n int) *Graph {
	g := PathGraph(n)
	if n >= 3 {
		g.AddEdge(n-1, 0)
	}
	return g
}

// CompleteGraph returns the complete graph K_n in which every pair of distinct
// vertices is joined by an edge.
func CompleteGraph(n int) *Graph {
	g := NewGraph(n)
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			g.AddEdge(i, j)
		}
	}
	return g
}

// StarGraph returns the star K_{1,n-1}: a central vertex 0 joined to each of the
// remaining n-1 vertices.
func StarGraph(n int) *Graph {
	g := NewGraph(n)
	for i := 1; i < n; i++ {
		g.AddEdge(0, i)
	}
	return g
}

// WheelGraph returns the wheel W_n: a cycle on vertices 1..n-1 plus a hub vertex
// 0 joined to every cycle vertex. It requires n >= 4 to form a proper wheel.
func WheelGraph(n int) *Graph {
	g := NewGraph(n)
	for i := 1; i < n; i++ {
		g.AddEdge(0, i)
	}
	for i := 1; i+1 < n; i++ {
		g.AddEdge(i, i+1)
	}
	if n >= 4 {
		g.AddEdge(n-1, 1)
	}
	return g
}

// CompleteBipartiteGraph returns K_{m,n}: vertices 0..m-1 form one side, vertices
// m..m+n-1 form the other, and every cross pair is joined.
func CompleteBipartiteGraph(m, n int) *Graph {
	g := NewGraph(m + n)
	for i := 0; i < m; i++ {
		for j := 0; j < n; j++ {
			g.AddEdge(i, m+j)
		}
	}
	return g
}

// GridGraph returns the r-by-c grid (lattice) graph, the Cartesian product of
// P_r and P_c, with vertex (i, j) indexed as i*c + j.
func GridGraph(r, c int) *Graph {
	g := NewGraph(r * c)
	idx := func(i, j int) int { return i*c + j }
	for i := 0; i < r; i++ {
		for j := 0; j < c; j++ {
			if j+1 < c {
				g.AddEdge(idx(i, j), idx(i, j+1))
			}
			if i+1 < r {
				g.AddEdge(idx(i, j), idx(i+1, j))
			}
		}
	}
	return g
}
