package graph

import (
	"math"
	"sort"
)

// EdmondsKarp computes the maximum flow from source to sink in a directed graph,
// treating each edge weight as its capacity, using the Edmonds-Karp algorithm
// (BFS-augmenting-path Ford-Fulkerson). It returns the value of the maximum
// flow. It returns [ErrVertexNotFound] if source or sink is absent, or
// [ErrNegativeWeight] if any capacity is negative. When source equals sink the
// flow is zero.
//
// The graph is not modified. An undirected graph is treated as having, for each
// edge {u,v}, forward capacity in both directions equal to the edge weight.
func (g *Graph) EdmondsKarp(source, sink int) (float64, error) {
	if !g.HasVertex(source) || !g.HasVertex(sink) {
		return 0, ErrVertexNotFound
	}
	// Build residual capacity map, ensuring reverse arcs exist.
	residual := make(map[int]map[int]float64, g.NumVertices())
	ensure := func(u int) {
		if residual[u] == nil {
			residual[u] = make(map[int]float64)
		}
	}
	for _, v := range g.Vertices() {
		ensure(v)
	}
	for u, nbrs := range g.adj {
		for v, w := range nbrs {
			if w < 0 {
				return 0, ErrNegativeWeight
			}
			residual[u][v] += w
			if _, ok := residual[v][u]; !ok {
				residual[v][u] = 0 // ensure a reverse residual arc exists
			}
		}
	}
	if source == sink {
		return 0, nil
	}

	maxFlow := 0.0
	for {
		// BFS to find an augmenting path in the residual graph.
		parent := map[int]int{source: source}
		queue := []int{source}
		found := false
		for len(queue) > 0 && !found {
			u := queue[0]
			queue = queue[1:]
			for _, v := range graphSortedKeys(residual[u]) {
				if _, seen := parent[v]; !seen && residual[u][v] > 0 {
					parent[v] = u
					if v == sink {
						found = true
						break
					}
					queue = append(queue, v)
				}
			}
		}
		if !found {
			break
		}
		// Bottleneck along the path.
		bottleneck := math.Inf(1)
		for v := sink; v != source; v = parent[v] {
			u := parent[v]
			if residual[u][v] < bottleneck {
				bottleneck = residual[u][v]
			}
		}
		// Augment.
		for v := sink; v != source; v = parent[v] {
			u := parent[v]
			residual[u][v] -= bottleneck
			residual[v][u] += bottleneck
		}
		maxFlow += bottleneck
	}
	return maxFlow, nil
}

// MinCut computes a minimum s-t cut of a directed graph whose edge weights are
// capacities. It returns the cut value (equal to the maximum flow), the set of
// vertices on the source side of the cut, and any error from the underlying
// [Graph.EdmondsKarp] computation. Vertices reachable from source in the final
// residual graph form the source side.
func (g *Graph) MinCut(source, sink int) (value float64, sourceSide []int, err error) {
	// Re-run max flow but keep the residual to find reachable set.
	if !g.HasVertex(source) || !g.HasVertex(sink) {
		return 0, nil, ErrVertexNotFound
	}
	residual := make(map[int]map[int]float64, g.NumVertices())
	for _, v := range g.Vertices() {
		residual[v] = make(map[int]float64)
	}
	for u, nbrs := range g.adj {
		for v, w := range nbrs {
			if w < 0 {
				return 0, nil, ErrNegativeWeight
			}
			residual[u][v] += w
			if _, ok := residual[v][u]; !ok {
				residual[v][u] = 0
			}
		}
	}
	total := 0.0
	if source != sink {
		for {
			parent := map[int]int{source: source}
			queue := []int{source}
			found := false
			for len(queue) > 0 && !found {
				u := queue[0]
				queue = queue[1:]
				for _, v := range graphSortedKeys(residual[u]) {
					if _, seen := parent[v]; !seen && residual[u][v] > 0 {
						parent[v] = u
						if v == sink {
							found = true
							break
						}
						queue = append(queue, v)
					}
				}
			}
			if !found {
				break
			}
			bottleneck := math.Inf(1)
			for v := sink; v != source; v = parent[v] {
				if residual[parent[v]][v] < bottleneck {
					bottleneck = residual[parent[v]][v]
				}
			}
			for v := sink; v != source; v = parent[v] {
				u := parent[v]
				residual[u][v] -= bottleneck
				residual[v][u] += bottleneck
			}
			total += bottleneck
		}
	}
	// Reachable set from source in residual graph.
	reach := map[int]bool{source: true}
	queue := []int{source}
	for len(queue) > 0 {
		u := queue[0]
		queue = queue[1:]
		for _, v := range graphSortedKeys(residual[u]) {
			if !reach[v] && residual[u][v] > 0 {
				reach[v] = true
				queue = append(queue, v)
			}
		}
	}
	for _, v := range g.Vertices() {
		if reach[v] {
			sourceSide = append(sourceSide, v)
		}
	}
	return total, sourceSide, nil
}

// MaxBipartiteMatching computes a maximum-cardinality matching between the left
// and right vertex sets of a bipartite graph using Kuhn's augmenting-path
// algorithm. Only edges from a left vertex to a right vertex are used
// (orientation is ignored, but at least one such adjacency must exist). It
// returns a map from each matched left vertex to its matched right vertex and
// the size of the matching. Left vertices are processed in ascending order for
// deterministic results.
func (g *Graph) MaxBipartiteMatching(left, right []int) (matching map[int]int, size int) {
	rightSet := make(map[int]bool, len(right))
	for _, r := range right {
		rightSet[r] = true
	}
	adj := g.graphUndirectedAdj()
	matchR := make(map[int]int) // right vertex -> left vertex
	leftSorted := append([]int(nil), left...)
	sort.Ints(leftSorted)

	var tryKuhn func(u int, visited map[int]bool) bool
	tryKuhn = func(u int, visited map[int]bool) bool {
		for _, v := range graphSortedKeys(adj[u]) {
			if !rightSet[v] || visited[v] {
				continue
			}
			visited[v] = true
			if w, ok := matchR[v]; !ok || tryKuhn(w, visited) {
				matchR[v] = u
				return true
			}
		}
		return false
	}

	for _, u := range leftSorted {
		visited := make(map[int]bool)
		tryKuhn(u, visited)
	}

	matching = make(map[int]int)
	for r, l := range matchR {
		matching[l] = r
	}
	return matching, len(matchR)
}
