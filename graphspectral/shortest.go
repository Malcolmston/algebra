package graphspectral

import "math"

// BFSDistances returns the number of edges on a shortest (unweighted) path from
// source to every vertex, ignoring edge weights. Unreachable vertices are marked
// with +Inf and the source with 0.
func BFSDistances(g *Graph, source int) []float64 {
	dist := make([]float64, g.n)
	for i := range dist {
		dist[i] = math.Inf(1)
	}
	if source < 0 || source >= g.n {
		return dist
	}
	dist[source] = 0
	queue := []int{source}
	for len(queue) > 0 {
		v := queue[0]
		queue = queue[1:]
		for w := 0; w < g.n; w++ {
			if w != v && g.adj.At(v, w) != 0 && math.IsInf(dist[w], 1) {
				dist[w] = dist[v] + 1
				queue = append(queue, w)
			}
		}
	}
	return dist
}

// DijkstraDistances returns the weighted shortest-path distance from source to
// every vertex, treating each edge weight as a positive length. Unreachable
// vertices are marked with +Inf. Negative weights are not supported and yield
// undefined results.
func DijkstraDistances(g *Graph, source int) []float64 {
	dist := make([]float64, g.n)
	for i := range dist {
		dist[i] = math.Inf(1)
	}
	if source < 0 || source >= g.n {
		return dist
	}
	visited := make([]bool, g.n)
	dist[source] = 0
	for iter := 0; iter < g.n; iter++ {
		u := -1
		best := math.Inf(1)
		for i := 0; i < g.n; i++ {
			if !visited[i] && dist[i] < best {
				best = dist[i]
				u = i
			}
		}
		if u == -1 {
			break
		}
		visited[u] = true
		for w := 0; w < g.n; w++ {
			wt := g.adj.At(u, w)
			if w == u || wt == 0 {
				continue
			}
			if nd := dist[u] + wt; nd < dist[w] {
				dist[w] = nd
			}
		}
	}
	return dist
}

// AllPairsShortestPaths returns the matrix of weighted shortest-path distances
// between every ordered pair of vertices, computed by running Dijkstra from each
// source. Unreachable pairs hold +Inf.
func AllPairsShortestPaths(g *Graph) *Matrix {
	m := NewMatrix(g.n, g.n)
	for s := 0; s < g.n; s++ {
		d := DijkstraDistances(g, s)
		for t := 0; t < g.n; t++ {
			m.Set(s, t, d[t])
		}
	}
	return m
}

// Diameter returns the greatest weighted shortest-path distance between any two
// vertices of a connected graph. It returns +Inf if the graph is disconnected.
func Diameter(g *Graph) float64 {
	var d float64
	for s := 0; s < g.n; s++ {
		dist := DijkstraDistances(g, s)
		for t := 0; t < g.n; t++ {
			if math.IsInf(dist[t], 1) {
				return math.Inf(1)
			}
			if dist[t] > d {
				d = dist[t]
			}
		}
	}
	return d
}
