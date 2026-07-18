package graph

// GreedyColoring assigns each vertex a color (a non-negative integer) so that no
// edge joins two vertices of the same color, using the greedy first-fit
// heuristic. Vertices are processed in ascending order and each receives the
// smallest color not used by an already-colored neighbor, ignoring edge
// orientation. The number of colors used is at most one more than the maximum
// degree; the result is a valid coloring but not necessarily optimal. A
// self-loop makes proper coloring impossible and is reported by the coloring
// still being returned, so callers that may have self-loops should validate with
// [Graph.IsProperColoring].
func (g *Graph) GreedyColoring() map[int]int {
	adj := g.graphUndirectedAdj()
	color := make(map[int]int, g.NumVertices())
	for _, v := range g.Vertices() {
		used := make(map[int]bool)
		for _, u := range graphSortedKeys(adj[v]) {
			if c, ok := color[u]; ok {
				used[c] = true
			}
		}
		c := 0
		for used[c] {
			c++
		}
		color[v] = c
	}
	return color
}

// NumColors returns the number of distinct colors that appear in a coloring.
func NumColors(color map[int]int) int {
	seen := make(map[int]bool)
	for _, c := range color {
		seen[c] = true
	}
	return len(seen)
}

// IsProperColoring reports whether color assigns different colors to the two
// endpoints of every edge, ignoring edge orientation. A graph with a self-loop
// can never be properly colored.
func (g *Graph) IsProperColoring(color map[int]int) bool {
	for u, nbrs := range g.adj {
		cu, ok := color[u]
		if !ok {
			return false
		}
		for v := range nbrs {
			cv, ok := color[v]
			if !ok {
				return false
			}
			if cu == cv {
				return false
			}
		}
	}
	return true
}

// ChromaticNumberUpperBound returns an easily computed upper bound on the
// chromatic number: one more than the maximum vertex degree (Brooks' bound in
// its simplest form). The true chromatic number is at most this value.
func (g *Graph) ChromaticNumberUpperBound() int {
	maxDeg := 0
	for _, v := range g.Vertices() {
		if d, err := g.Degree(v); err == nil && d > maxDeg {
			maxDeg = d
		}
	}
	return maxDeg + 1
}
