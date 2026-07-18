package graph

import "sort"

// BFS returns the vertices reachable from start in breadth-first order,
// including start itself. Neighbors are visited in ascending order. It returns
// [ErrVertexNotFound] if start is absent.
func (g *Graph) BFS(start int) ([]int, error) {
	if !g.HasVertex(start) {
		return nil, ErrVertexNotFound
	}
	visited := map[int]bool{start: true}
	order := []int{start}
	queue := []int{start}
	for len(queue) > 0 {
		u := queue[0]
		queue = queue[1:]
		for _, v := range graphSortedKeys(g.adj[u]) {
			if !visited[v] {
				visited[v] = true
				order = append(order, v)
				queue = append(queue, v)
			}
		}
	}
	return order, nil
}

// BFSDistances returns, for every vertex reachable from start, the minimum
// number of edges from start to that vertex (start maps to 0). Unreachable
// vertices are omitted. It returns [ErrVertexNotFound] if start is absent.
func (g *Graph) BFSDistances(start int) (map[int]int, error) {
	if !g.HasVertex(start) {
		return nil, ErrVertexNotFound
	}
	dist := map[int]int{start: 0}
	queue := []int{start}
	for len(queue) > 0 {
		u := queue[0]
		queue = queue[1:]
		for _, v := range graphSortedKeys(g.adj[u]) {
			if _, seen := dist[v]; !seen {
				dist[v] = dist[u] + 1
				queue = append(queue, v)
			}
		}
	}
	return dist, nil
}

// ShortestPathBFS returns a path with the fewest edges from start to goal in an
// unweighted sense, and whether such a path exists. Edge weights are ignored.
// The returned path begins with start and ends with goal.
func (g *Graph) ShortestPathBFS(start, goal int) ([]int, bool) {
	if !g.HasVertex(start) || !g.HasVertex(goal) {
		return nil, false
	}
	prev := map[int]int{start: start}
	queue := []int{start}
	for len(queue) > 0 {
		u := queue[0]
		queue = queue[1:]
		if u == goal {
			return graphReconstruct(prev, start, goal), true
		}
		for _, v := range graphSortedKeys(g.adj[u]) {
			if _, seen := prev[v]; !seen {
				prev[v] = u
				queue = append(queue, v)
			}
		}
	}
	return nil, false
}

// DFS returns the vertices reachable from start in depth-first preorder,
// including start itself. Neighbors are visited in ascending order. It returns
// [ErrVertexNotFound] if start is absent.
func (g *Graph) DFS(start int) ([]int, error) {
	return g.DFSPreorder(start)
}

// DFSPreorder returns the vertices reachable from start in depth-first
// preorder: each vertex is recorded when first visited. Neighbors are visited
// in ascending order. It returns [ErrVertexNotFound] if start is absent.
func (g *Graph) DFSPreorder(start int) ([]int, error) {
	if !g.HasVertex(start) {
		return nil, ErrVertexNotFound
	}
	visited := make(map[int]bool)
	var order []int
	var dfs func(u int)
	dfs = func(u int) {
		visited[u] = true
		order = append(order, u)
		for _, v := range graphSortedKeys(g.adj[u]) {
			if !visited[v] {
				dfs(v)
			}
		}
	}
	dfs(start)
	return order, nil
}

// DFSPostorder returns the vertices reachable from start in depth-first
// postorder: each vertex is recorded after all its descendants. Neighbors are
// visited in ascending order. It returns [ErrVertexNotFound] if start is absent.
func (g *Graph) DFSPostorder(start int) ([]int, error) {
	if !g.HasVertex(start) {
		return nil, ErrVertexNotFound
	}
	visited := make(map[int]bool)
	var order []int
	var dfs func(u int)
	dfs = func(u int) {
		visited[u] = true
		for _, v := range graphSortedKeys(g.adj[u]) {
			if !visited[v] {
				dfs(v)
			}
		}
		order = append(order, u)
	}
	dfs(start)
	return order, nil
}

// ReachableFrom returns the set of vertices reachable from start (including
// start) as a map. It returns [ErrVertexNotFound] if start is absent.
func (g *Graph) ReachableFrom(start int) (map[int]bool, error) {
	order, err := g.BFS(start)
	if err != nil {
		return nil, err
	}
	set := make(map[int]bool, len(order))
	for _, v := range order {
		set[v] = true
	}
	return set, nil
}

// graphUndirectedAdj returns an undirected view of the adjacency relation: for
// every edge u->v both u->v and v->u are present. It is used by algorithms that
// treat a directed graph as its underlying (weakly connected) undirected graph.
func (g *Graph) graphUndirectedAdj() map[int]map[int]bool {
	adj := make(map[int]map[int]bool, len(g.adj))
	for u := range g.adj {
		adj[u] = make(map[int]bool)
	}
	for u, nbrs := range g.adj {
		for v := range nbrs {
			adj[u][v] = true
			adj[v][u] = true
		}
	}
	return adj
}

// ConnectedComponents returns the connected components of the graph. For a
// directed graph the components are the weakly connected components (edge
// orientation is ignored). Each component is a slice of vertices in ascending
// order, and the components are ordered by their smallest vertex.
func (g *Graph) ConnectedComponents() [][]int {
	adj := g.graphUndirectedAdj()
	visited := make(map[int]bool)
	var comps [][]int
	for _, s := range g.Vertices() {
		if visited[s] {
			continue
		}
		var comp []int
		queue := []int{s}
		visited[s] = true
		for len(queue) > 0 {
			u := queue[0]
			queue = queue[1:]
			comp = append(comp, u)
			for _, v := range graphSortedKeys(adj[u]) {
				if !visited[v] {
					visited[v] = true
					queue = append(queue, v)
				}
			}
		}
		sort.Ints(comp)
		comps = append(comps, comp)
	}
	sort.Slice(comps, func(i, j int) bool { return comps[i][0] < comps[j][0] })
	return comps
}

// NumConnectedComponents returns the number of connected components (weakly
// connected for directed graphs). An empty graph has zero components.
func (g *Graph) NumConnectedComponents() int {
	return len(g.ConnectedComponents())
}

// IsConnected reports whether the graph is connected (weakly connected for
// directed graphs). The empty graph and single-vertex graph are considered
// connected.
func (g *Graph) IsConnected() bool {
	return g.NumConnectedComponents() <= 1
}

// TopologicalSort returns a topological ordering of a directed acyclic graph
// using Kahn's algorithm. Among vertices whose prerequisites are all satisfied
// it always emits the smallest identifier first, giving a deterministic order.
// It returns ErrNotDAG if the graph is undirected or contains a cycle.
func (g *Graph) TopologicalSort() ([]int, error) {
	if !g.directed {
		return nil, ErrNotDAG
	}
	indeg := make(map[int]int)
	for _, v := range g.Vertices() {
		indeg[v] = 0
	}
	for _, nbrs := range g.adj {
		for v := range nbrs {
			indeg[v]++
		}
	}
	// Seed with all zero in-degree vertices, kept sorted.
	var ready []int
	for _, v := range g.Vertices() {
		if indeg[v] == 0 {
			ready = append(ready, v)
		}
	}
	var order []int
	for len(ready) > 0 {
		// Pop the smallest ready vertex.
		u := ready[0]
		ready = ready[1:]
		order = append(order, u)
		for _, v := range graphSortedKeys(g.adj[u]) {
			indeg[v]--
			if indeg[v] == 0 {
				// Insert v keeping ready sorted.
				i := sort.SearchInts(ready, v)
				ready = append(ready, 0)
				copy(ready[i+1:], ready[i:])
				ready[i] = v
			}
		}
	}
	if len(order) != g.NumVertices() {
		return nil, ErrNotDAG
	}
	return order, nil
}

// HasCycle reports whether the graph contains a cycle. For directed graphs it
// detects directed cycles; for undirected graphs it detects any cycle,
// disregarding trivial back-edges to the immediate parent. A self-loop always
// counts as a cycle.
func (g *Graph) HasCycle() bool {
	if g.directed {
		const (
			white = 0
			gray  = 1
			black = 2
		)
		color := make(map[int]int)
		var dfs func(u int) bool
		dfs = func(u int) bool {
			color[u] = gray
			for _, v := range graphSortedKeys(g.adj[u]) {
				switch color[v] {
				case gray:
					return true
				case white:
					if dfs(v) {
						return true
					}
				}
			}
			color[u] = black
			return false
		}
		for _, s := range g.Vertices() {
			if color[s] == white && dfs(s) {
				return true
			}
		}
		return false
	}
	// Undirected: DFS tracking parent; self-loops and parallel handling.
	visited := make(map[int]bool)
	var dfs func(u, parent int) bool
	dfs = func(u, parent int) bool {
		visited[u] = true
		for _, v := range graphSortedKeys(g.adj[u]) {
			if v == u {
				return true // self-loop
			}
			if !visited[v] {
				if dfs(v, u) {
					return true
				}
			} else if v != parent {
				return true
			}
		}
		return false
	}
	for _, s := range g.Vertices() {
		if !visited[s] && dfs(s, -1) {
			return true
		}
	}
	return false
}

// IsDAG reports whether the graph is a directed acyclic graph.
func (g *Graph) IsDAG() bool {
	return g.directed && !g.HasCycle()
}

// IsTree reports whether an undirected graph is a tree: connected and acyclic
// with exactly NumVertices-1 edges. The empty graph is not a tree; a single
// vertex is. Directed graphs are never reported as trees by this method.
func (g *Graph) IsTree() bool {
	if g.directed || g.NumVertices() == 0 {
		return false
	}
	return g.IsConnected() && g.NumEdges() == g.NumVertices()-1
}

// IsBipartite reports whether the graph is bipartite (2-colorable), ignoring
// edge orientation. A graph with no edges is bipartite.
func (g *Graph) IsBipartite() bool {
	_, ok := g.TwoColoring()
	return ok
}

// TwoColoring attempts to 2-color the graph so that adjacent vertices receive
// different colors (0 or 1), ignoring edge orientation. It returns the coloring
// and true if the graph is bipartite, or nil and false otherwise.
func (g *Graph) TwoColoring() (map[int]int, bool) {
	adj := g.graphUndirectedAdj()
	color := make(map[int]int)
	for _, s := range g.Vertices() {
		if _, seen := color[s]; seen {
			continue
		}
		color[s] = 0
		queue := []int{s}
		for len(queue) > 0 {
			u := queue[0]
			queue = queue[1:]
			for _, v := range graphSortedKeys(adj[u]) {
				if v == u {
					return nil, false // self-loop is not 2-colorable
				}
				if c, seen := color[v]; !seen {
					color[v] = 1 - color[u]
					queue = append(queue, v)
				} else if c == color[u] {
					return nil, false
				}
			}
		}
	}
	return color, true
}

// graphReconstruct rebuilds a path from start to goal using a predecessor map
// where prev[start] == start. It returns nil if goal is unreachable.
func graphReconstruct(prev map[int]int, start, goal int) []int {
	if _, ok := prev[goal]; !ok {
		return nil
	}
	var rev []int
	for at := goal; ; at = prev[at] {
		rev = append(rev, at)
		if at == start {
			break
		}
	}
	for i, j := 0, len(rev)-1; i < j; i, j = i+1, j-1 {
		rev[i], rev[j] = rev[j], rev[i]
	}
	return rev
}
