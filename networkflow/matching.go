package networkflow

import "sort"

// BipartiteGraph is an undirected bipartite graph with a left vertex set
// {0..L-1} and a right vertex set {0..R-1}. Edges join a left vertex to a right
// vertex. It is the input type for maximum-cardinality matching.
type BipartiteGraph struct {
	nl  int
	nr  int
	adj [][]int
	set [][]bool
}

// NewBipartiteGraph returns an empty bipartite graph with nl left vertices and
// nr right vertices. It panics if either count is negative.
func NewBipartiteGraph(nl, nr int) *BipartiteGraph {
	if nl < 0 || nr < 0 {
		panic("networkflow: negative vertex count")
	}
	set := make([][]bool, nl)
	for i := range set {
		set[i] = make([]bool, nr)
	}
	return &BipartiteGraph{nl: nl, nr: nr, adj: make([][]int, nl), set: set}
}

// LeftSize returns the number of left vertices.
func (g *BipartiteGraph) LeftSize() int { return g.nl }

// RightSize returns the number of right vertices.
func (g *BipartiteGraph) RightSize() int { return g.nr }

// AddEdge adds an edge between left vertex u and right vertex v. Duplicate
// edges are ignored. It panics on an out-of-range vertex.
func (g *BipartiteGraph) AddEdge(u, v int) {
	if u < 0 || u >= g.nl || v < 0 || v >= g.nr {
		panic(ErrInvalidVertex)
	}
	if g.set[u][v] {
		return
	}
	g.set[u][v] = true
	g.adj[u] = append(g.adj[u], v)
}

// HasEdge reports whether an edge joins left vertex u and right vertex v.
func (g *BipartiteGraph) HasEdge(u, v int) bool {
	if u < 0 || u >= g.nl || v < 0 || v >= g.nr {
		return false
	}
	return g.set[u][v]
}

// LeftNeighbors returns the right vertices adjacent to left vertex u, sorted.
func (g *BipartiteGraph) LeftNeighbors(u int) []int {
	out := make([]int, len(g.adj[u]))
	copy(out, g.adj[u])
	sort.Ints(out)
	return out
}

// NumEdges returns the total number of edges.
func (g *BipartiteGraph) NumEdges() int {
	c := 0
	for u := 0; u < g.nl; u++ {
		c += len(g.adj[u])
	}
	return c
}

// Clone returns a deep copy of the bipartite graph.
func (g *BipartiteGraph) Clone() *BipartiteGraph {
	c := NewBipartiteGraph(g.nl, g.nr)
	for u := 0; u < g.nl; u++ {
		c.adj[u] = make([]int, len(g.adj[u]))
		copy(c.adj[u], g.adj[u])
		copy(c.set[u], g.set[u])
	}
	return c
}

// MatchingResult describes a matching in a bipartite graph.
type MatchingResult struct {
	// Size is the number of matched pairs.
	Size int
	// MatchL maps each left vertex to its matched right vertex, or -1.
	MatchL []int
	// MatchR maps each right vertex to its matched left vertex, or -1.
	MatchR []int
}

// Pairs returns the matched (left, right) pairs sorted by left vertex.
func (m *MatchingResult) Pairs() [][2]int {
	var out [][2]int
	for u, v := range m.MatchL {
		if v >= 0 {
			out = append(out, [2]int{u, v})
		}
	}
	return out
}

// IsLeftMatched reports whether left vertex u is matched.
func (m *MatchingResult) IsLeftMatched(u int) bool { return m.MatchL[u] >= 0 }

// IsRightMatched reports whether right vertex v is matched.
func (m *MatchingResult) IsRightMatched(v int) bool { return m.MatchR[v] >= 0 }

// IsPerfect reports whether every vertex on the smaller side is matched, i.e.
// the matching size equals min(LeftSize, RightSize).
func (m *MatchingResult) IsPerfect() bool {
	return m.Size == min(len(m.MatchL), len(m.MatchR))
}

// hopcroftKarp computes a maximum-cardinality matching and returns matchL and
// matchR arrays (values of -1 mean unmatched).
func hopcroftKarp(g *BipartiteGraph) ([]int, []int) {
	const inf = 1 << 30
	matchL := make([]int, g.nl)
	matchR := make([]int, g.nr)
	for i := range matchL {
		matchL[i] = -1
	}
	for i := range matchR {
		matchR[i] = -1
	}
	dist := make([]int, g.nl)

	bfs := func() bool {
		queue := make([]int, 0, g.nl)
		for u := 0; u < g.nl; u++ {
			if matchL[u] == -1 {
				dist[u] = 0
				queue = append(queue, u)
			} else {
				dist[u] = inf
			}
		}
		found := false
		for len(queue) > 0 {
			u := queue[0]
			queue = queue[1:]
			for _, v := range g.adj[u] {
				w := matchR[v]
				if w == -1 {
					found = true
				} else if dist[w] == inf {
					dist[w] = dist[u] + 1
					queue = append(queue, w)
				}
			}
		}
		return found
	}

	var dfs func(u int) bool
	dfs = func(u int) bool {
		for _, v := range g.adj[u] {
			w := matchR[v]
			if w == -1 || (dist[w] == dist[u]+1 && dfs(w)) {
				matchL[u] = v
				matchR[v] = u
				return true
			}
		}
		dist[u] = inf
		return false
	}

	for bfs() {
		for u := 0; u < g.nl; u++ {
			if matchL[u] == -1 {
				dfs(u)
			}
		}
	}
	return matchL, matchR
}

// kuhn computes a maximum-cardinality matching with the simple Kuhn
// augmenting-path method and returns matchL and matchR arrays.
func kuhn(g *BipartiteGraph) ([]int, []int) {
	matchL := make([]int, g.nl)
	matchR := make([]int, g.nr)
	for i := range matchL {
		matchL[i] = -1
	}
	for i := range matchR {
		matchR[i] = -1
	}
	var used []bool
	var try func(u int) bool
	try = func(u int) bool {
		for _, v := range g.adj[u] {
			if used[v] {
				continue
			}
			used[v] = true
			if matchR[v] == -1 || try(matchR[v]) {
				matchL[u] = v
				matchR[v] = u
				return true
			}
		}
		return false
	}
	for u := 0; u < g.nl; u++ {
		used = make([]bool, g.nr)
		try(u)
	}
	return matchL, matchR
}

// resultFrom packages matchL and matchR into a MatchingResult.
func resultFrom(matchL, matchR []int) *MatchingResult {
	size := 0
	for _, v := range matchL {
		if v >= 0 {
			size++
		}
	}
	return &MatchingResult{Size: size, MatchL: matchL, MatchR: matchR}
}

// HopcroftKarp returns a maximum-cardinality matching of a bipartite graph
// computed with the Hopcroft-Karp algorithm in O(E*sqrt(V)) time. The input is
// left unchanged.
func HopcroftKarp(g *BipartiteGraph) *MatchingResult {
	l, r := hopcroftKarp(g)
	return resultFrom(l, r)
}

// KuhnMatching returns a maximum-cardinality matching computed with the simpler
// Kuhn augmenting-path method in O(V*E) time. The input is left unchanged.
func KuhnMatching(g *BipartiteGraph) *MatchingResult {
	l, r := kuhn(g)
	return resultFrom(l, r)
}

// MaximumMatching returns a maximum-cardinality matching using the default
// engine (Hopcroft-Karp).
func MaximumMatching(g *BipartiteGraph) *MatchingResult { return HopcroftKarp(g) }

// MaximumMatchingSize returns just the size of a maximum-cardinality matching.
func MaximumMatchingSize(g *BipartiteGraph) int { return HopcroftKarp(g).Size }

// MinimumVertexCover returns a minimum vertex cover of a bipartite graph as two
// sets (left cover vertices, right cover vertices) via Konig's theorem applied
// to a maximum matching. The total size equals the maximum matching size. The
// input is left unchanged.
func MinimumVertexCover(g *BipartiteGraph) (left, right []int) {
	m := HopcroftKarp(g)
	// Z = vertices reachable from unmatched left vertices via alternating paths.
	visitedL := make([]bool, g.nl)
	visitedR := make([]bool, g.nr)
	var dfs func(u int)
	dfs = func(u int) {
		visitedL[u] = true
		for _, v := range g.adj[u] {
			// Follow only non-matching edges L->R.
			if m.MatchL[u] == v {
				continue
			}
			if !visitedR[v] {
				visitedR[v] = true
				w := m.MatchR[v]
				if w != -1 && !visitedL[w] {
					dfs(w)
				}
			}
		}
	}
	for u := 0; u < g.nl; u++ {
		if m.MatchL[u] == -1 {
			dfs(u)
		}
	}
	// Cover = (left not in Z) union (right in Z).
	for u := 0; u < g.nl; u++ {
		if !visitedL[u] {
			left = append(left, u)
		}
	}
	for v := 0; v < g.nr; v++ {
		if visitedR[v] {
			right = append(right, v)
		}
	}
	return left, right
}

// MaximumIndependentSet returns a maximum independent set of a bipartite graph
// as (left, right) vertex lists. It is the complement of a minimum vertex cover
// (Konig's theorem). The input is left unchanged.
func MaximumIndependentSet(g *BipartiteGraph) (left, right []int) {
	cl, cr := MinimumVertexCover(g)
	inCL := make([]bool, g.nl)
	inCR := make([]bool, g.nr)
	for _, u := range cl {
		inCL[u] = true
	}
	for _, v := range cr {
		inCR[v] = true
	}
	for u := 0; u < g.nl; u++ {
		if !inCL[u] {
			left = append(left, u)
		}
	}
	for v := 0; v < g.nr; v++ {
		if !inCR[v] {
			right = append(right, v)
		}
	}
	return left, right
}

// HasPerfectMatching reports whether the bipartite graph has a matching that
// saturates the smaller side.
func HasPerfectMatching(g *BipartiteGraph) bool { return HopcroftKarp(g).IsPerfect() }
