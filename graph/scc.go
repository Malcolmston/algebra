package graph

import "sort"

// TarjanSCC returns the strongly connected components of a directed graph using
// Tarjan's algorithm. For an undirected graph the strongly connected components
// coincide with the connected components. Each component is returned as a slice
// of vertices in ascending order, and the components are ordered by their
// smallest vertex.
func (g *Graph) TarjanSCC() [][]int {
	index := make(map[int]int)
	low := make(map[int]int)
	onStack := make(map[int]bool)
	var stack []int
	counter := 0
	var comps [][]int

	var strongConnect func(v int)
	strongConnect = func(v int) {
		index[v] = counter
		low[v] = counter
		counter++
		stack = append(stack, v)
		onStack[v] = true
		for _, w := range graphSortedKeys(g.adj[v]) {
			if _, seen := index[w]; !seen {
				strongConnect(w)
				if low[w] < low[v] {
					low[v] = low[w]
				}
			} else if onStack[w] {
				if index[w] < low[v] {
					low[v] = index[w]
				}
			}
		}
		if low[v] == index[v] {
			var comp []int
			for {
				w := stack[len(stack)-1]
				stack = stack[:len(stack)-1]
				onStack[w] = false
				comp = append(comp, w)
				if w == v {
					break
				}
			}
			sort.Ints(comp)
			comps = append(comps, comp)
		}
	}

	for _, v := range g.Vertices() {
		if _, seen := index[v]; !seen {
			strongConnect(v)
		}
	}
	sort.Slice(comps, func(i, j int) bool { return comps[i][0] < comps[j][0] })
	return comps
}

// KosarajuSCC returns the strongly connected components of a directed graph
// using Kosaraju's two-pass algorithm. The result matches [Graph.TarjanSCC]:
// each component is sorted ascending and the components are ordered by their
// smallest vertex.
func (g *Graph) KosarajuSCC() [][]int {
	visited := make(map[int]bool)
	var order []int
	var dfs1 func(u int)
	dfs1 = func(u int) {
		visited[u] = true
		for _, v := range graphSortedKeys(g.adj[u]) {
			if !visited[v] {
				dfs1(v)
			}
		}
		order = append(order, u)
	}
	for _, v := range g.Vertices() {
		if !visited[v] {
			dfs1(v)
		}
	}

	rev := g.Reverse()
	assigned := make(map[int]bool)
	var comps [][]int
	var dfs2 func(u int, comp *[]int)
	dfs2 = func(u int, comp *[]int) {
		assigned[u] = true
		*comp = append(*comp, u)
		for _, v := range graphSortedKeys(rev.adj[u]) {
			if !assigned[v] {
				dfs2(v, comp)
			}
		}
	}
	for i := len(order) - 1; i >= 0; i-- {
		u := order[i]
		if !assigned[u] {
			var comp []int
			dfs2(u, &comp)
			sort.Ints(comp)
			comps = append(comps, comp)
		}
	}
	sort.Slice(comps, func(i, j int) bool { return comps[i][0] < comps[j][0] })
	return comps
}

// TransitiveClosure returns the reachability relation of the graph as a nested
// map: closure[u][v] is true when there is a directed path (of one or more
// edges) from u to v. A vertex reaches itself only if it lies on a cycle or has
// a self-loop.
func (g *Graph) TransitiveClosure() map[int]map[int]bool {
	closure := make(map[int]map[int]bool, g.NumVertices())
	for _, u := range g.Vertices() {
		closure[u] = make(map[int]bool)
		// BFS from u following directed edges; do not mark u unless reached.
		queue := graphSortedKeys(g.adj[u])
		seen := make(map[int]bool)
		for _, v := range queue {
			seen[v] = true
		}
		for len(queue) > 0 {
			x := queue[0]
			queue = queue[1:]
			closure[u][x] = true
			for _, y := range graphSortedKeys(g.adj[x]) {
				if !seen[y] {
					seen[y] = true
					queue = append(queue, y)
				}
			}
		}
	}
	return closure
}

// Condensation returns the condensation of a directed graph: a new directed
// graph whose vertices are the indices of the strongly connected components (as
// ordered by [Graph.TarjanSCC]) with an edge i->j whenever some edge of the
// original graph runs from component i to a different component j. The
// condensation is always acyclic. The second result maps each original vertex
// to its component index.
func (g *Graph) Condensation() (*Graph, map[int]int) {
	sccs := g.TarjanSCC()
	compOf := make(map[int]int)
	for i, comp := range sccs {
		for _, v := range comp {
			compOf[v] = i
		}
	}
	c := NewDirected()
	for i := range sccs {
		c.AddVertex(i)
	}
	for u, nbrs := range g.adj {
		for v := range nbrs {
			if compOf[u] != compOf[v] {
				c.AddEdge(compOf[u], compOf[v])
			}
		}
	}
	return c, compOf
}
