package graph

import (
	"container/heap"
	"math"
	"sort"
)

// graphPQItem is an entry in the priority queue used by Dijkstra and A*.
type graphPQItem struct {
	vertex   int
	priority float64
}

// graphPQ is a min-heap of graphPQItem ordered by priority, breaking ties by
// ascending vertex identifier for determinism.
type graphPQ []graphPQItem

// Len reports the number of items in the queue, implementing heap.Interface.
func (pq graphPQ) Len() int { return len(pq) }

// Less reports whether item i orders before item j, implementing heap.Interface.
// Items are ordered by ascending priority, breaking ties by ascending vertex
// identifier.
func (pq graphPQ) Less(i, j int) bool {
	if pq[i].priority != pq[j].priority {
		return pq[i].priority < pq[j].priority
	}
	return pq[i].vertex < pq[j].vertex
}

// Swap exchanges items i and j, implementing heap.Interface.
func (pq graphPQ) Swap(i, j int) { pq[i], pq[j] = pq[j], pq[i] }

// Push appends x to the queue, implementing heap.Interface.
func (pq *graphPQ) Push(x any) { *pq = append(*pq, x.(graphPQItem)) }

// Pop removes and returns the last item of the queue, implementing heap.Interface.
func (pq *graphPQ) Pop() any {
	old := *pq
	n := len(old)
	it := old[n-1]
	*pq = old[:n-1]
	return it
}

// Dijkstra computes single-source shortest-path distances from start over
// non-negative edge weights. It returns a map of reachable vertices to their
// distance from start and a predecessor map suitable for [Graph.PathTo].
// Unreachable vertices are omitted. It returns [ErrVertexNotFound] if start is
// absent, or [ErrNegativeWeight] if any reachable edge weight is negative.
func (g *Graph) Dijkstra(start int) (dist map[int]float64, prev map[int]int, err error) {
	if !g.HasVertex(start) {
		return nil, nil, ErrVertexNotFound
	}
	dist = map[int]float64{start: 0}
	prev = map[int]int{start: start}
	pq := &graphPQ{{vertex: start, priority: 0}}
	settled := make(map[int]bool)
	for pq.Len() > 0 {
		item := heap.Pop(pq).(graphPQItem)
		u := item.vertex
		if settled[u] {
			continue
		}
		settled[u] = true
		for _, v := range graphSortedKeys(g.adj[u]) {
			w := g.adj[u][v]
			if w < 0 {
				return nil, nil, ErrNegativeWeight
			}
			nd := dist[u] + w
			if d, ok := dist[v]; !ok || nd < d {
				dist[v] = nd
				prev[v] = u
				heap.Push(pq, graphPQItem{vertex: v, priority: nd})
			}
		}
	}
	return dist, prev, nil
}

// DijkstraPath returns the minimum-weight path from start to goal over
// non-negative edge weights, together with its total weight. The boolean result
// reports whether goal is reachable. It returns an error under the same
// conditions as [Graph.Dijkstra].
func (g *Graph) DijkstraPath(start, goal int) (path []int, weight float64, err error) {
	dist, prev, err := g.Dijkstra(start)
	if err != nil {
		return nil, 0, err
	}
	d, ok := dist[goal]
	if !ok {
		return nil, 0, nil
	}
	return graphReconstruct(prev, start, goal), d, nil
}

// PathTo reconstructs the path from start to goal from a predecessor map as
// returned by [Graph.Dijkstra], [Graph.BellmanFord] or [Graph.AStar]. The start
// vertex must map to itself in prev. It returns nil if goal is unreachable.
func (g *Graph) PathTo(prev map[int]int, start, goal int) []int {
	return graphReconstruct(prev, start, goal)
}

// BellmanFord computes single-source shortest-path distances from start and
// tolerates negative edge weights. It returns a map of reachable vertices to
// their distance and a predecessor map suitable for [Graph.PathTo]. It returns
// [ErrVertexNotFound] if start is absent, or [ErrNegativeCycle] if a
// negative-weight cycle is reachable from start.
func (g *Graph) BellmanFord(start int) (dist map[int]float64, prev map[int]int, err error) {
	if !g.HasVertex(start) {
		return nil, nil, ErrVertexNotFound
	}
	dist = map[int]float64{start: 0}
	prev = map[int]int{start: start}
	edges := g.graphDirectedEdges()
	n := g.NumVertices()
	for i := 0; i < n-1; i++ {
		changed := false
		for _, e := range edges {
			du, ok := dist[e.From]
			if !ok {
				continue
			}
			nd := du + e.Weight
			if d, ok := dist[e.To]; !ok || nd < d {
				dist[e.To] = nd
				prev[e.To] = e.From
				changed = true
			}
		}
		if !changed {
			break
		}
	}
	// One more relaxation detects a reachable negative cycle.
	for _, e := range edges {
		du, ok := dist[e.From]
		if !ok {
			continue
		}
		if d, ok := dist[e.To]; !ok || du+e.Weight < d {
			return nil, nil, ErrNegativeCycle
		}
	}
	return dist, prev, nil
}

// BellmanFordPath returns the minimum-weight path from start to goal, tolerating
// negative edge weights, together with its total weight. The boolean result
// reports whether goal is reachable. It returns an error under the same
// conditions as [Graph.BellmanFord].
func (g *Graph) BellmanFordPath(start, goal int) (path []int, weight float64, err error) {
	dist, prev, err := g.BellmanFord(start)
	if err != nil {
		return nil, 0, err
	}
	d, ok := dist[goal]
	if !ok {
		return nil, 0, nil
	}
	return graphReconstruct(prev, start, goal), d, nil
}

// FloydWarshall computes all-pairs shortest-path distances. It returns a nested
// map where dist[u][v] is the shortest distance from u to v, or
// math.Inf(1) if v is unreachable from u; dist[u][u] is 0 unless a negative
// self-reachable cycle exists. It returns [ErrNegativeCycle] if the graph
// contains a negative-weight cycle.
func (g *Graph) FloydWarshall() (map[int]map[int]float64, error) {
	verts := g.Vertices()
	dist := make(map[int]map[int]float64, len(verts))
	for _, u := range verts {
		row := make(map[int]float64, len(verts))
		for _, v := range verts {
			if u == v {
				row[v] = 0
			} else {
				row[v] = math.Inf(1)
			}
		}
		dist[u] = row
	}
	for u, nbrs := range g.adj {
		for v, w := range nbrs {
			if w < dist[u][v] {
				dist[u][v] = w
			}
		}
	}
	for _, k := range verts {
		for _, i := range verts {
			if math.IsInf(dist[i][k], 1) {
				continue
			}
			for _, j := range verts {
				if nd := dist[i][k] + dist[k][j]; nd < dist[i][j] {
					dist[i][j] = nd
				}
			}
		}
	}
	for _, v := range verts {
		if dist[v][v] < 0 {
			return nil, ErrNegativeCycle
		}
	}
	return dist, nil
}

// AStar finds a minimum-weight path from start to goal over non-negative edge
// weights, guided by the heuristic h, which must estimate the remaining cost
// from a vertex to goal and be admissible (never overestimate) for the result
// to be optimal. It returns the path, its total weight, and whether goal was
// reached. It returns [ErrVertexNotFound] if start or goal is absent, or
// [ErrNegativeWeight] if a relaxed edge weight is negative.
func (g *Graph) AStar(start, goal int, h func(int) float64) (path []int, weight float64, err error) {
	if !g.HasVertex(start) || !g.HasVertex(goal) {
		return nil, 0, ErrVertexNotFound
	}
	gScore := map[int]float64{start: 0}
	prev := map[int]int{start: start}
	pq := &graphPQ{{vertex: start, priority: h(start)}}
	settled := make(map[int]bool)
	for pq.Len() > 0 {
		item := heap.Pop(pq).(graphPQItem)
		u := item.vertex
		if settled[u] {
			continue
		}
		if u == goal {
			return graphReconstruct(prev, start, goal), gScore[goal], nil
		}
		settled[u] = true
		for _, v := range graphSortedKeys(g.adj[u]) {
			w := g.adj[u][v]
			if w < 0 {
				return nil, 0, ErrNegativeWeight
			}
			nd := gScore[u] + w
			if d, ok := gScore[v]; !ok || nd < d {
				gScore[v] = nd
				prev[v] = u
				heap.Push(pq, graphPQItem{vertex: v, priority: nd + h(v)})
			}
		}
	}
	return nil, 0, nil
}

// graphDirectedEdges returns every directed arc of the graph. For an undirected
// graph each edge yields two arcs (u->v and v->u). Arcs are sorted by
// (From, To) for deterministic relaxation order.
func (g *Graph) graphDirectedEdges() []Edge {
	var edges []Edge
	for u, nbrs := range g.adj {
		for v, w := range nbrs {
			edges = append(edges, Edge{From: u, To: v, Weight: w})
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
