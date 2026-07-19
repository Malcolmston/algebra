package graph

import (
	"container/heap"
	"errors"
	"sort"
)

// ErrNotUndirected reports that an operation requiring an undirected graph (such
// as [Graph.Kruskal] or [Graph.Prim]) was applied to a directed graph.
var ErrNotUndirected = errors.New("graph: operation requires an undirected graph")

// UnionFind is a disjoint-set (union-find) data structure over a fixed set of
// integer elements, with union by rank and path compression. It answers
// connectivity queries in near-constant amortized time and is used by
// [Graph.Kruskal].
type UnionFind struct {
	parent map[int]int
	rank   map[int]int
	count  int
}

// NewUnionFind returns a union-find structure whose given elements each start in
// their own singleton set. Duplicate elements are ignored.
func NewUnionFind(elements []int) *UnionFind {
	uf := &UnionFind{parent: make(map[int]int), rank: make(map[int]int)}
	for _, e := range elements {
		if _, ok := uf.parent[e]; !ok {
			uf.parent[e] = e
			uf.rank[e] = 0
			uf.count++
		}
	}
	return uf
}

// MakeSet adds x as a new singleton set if it is not already present.
func (uf *UnionFind) MakeSet(x int) {
	if _, ok := uf.parent[x]; !ok {
		uf.parent[x] = x
		uf.rank[x] = 0
		uf.count++
	}
}

// Find returns the canonical representative of the set containing x, adding x as
// a new singleton if it was not present. Path compression is applied.
func (uf *UnionFind) Find(x int) int {
	if _, ok := uf.parent[x]; !ok {
		uf.MakeSet(x)
		return x
	}
	root := x
	for uf.parent[root] != root {
		root = uf.parent[root]
	}
	for uf.parent[x] != root {
		uf.parent[x], x = root, uf.parent[x]
	}
	return root
}

// Union merges the sets containing x and y and reports whether a merge occurred
// (that is, whether they were previously in different sets). Absent elements are
// added automatically.
func (uf *UnionFind) Union(x, y int) bool {
	rx, ry := uf.Find(x), uf.Find(y)
	if rx == ry {
		return false
	}
	if uf.rank[rx] < uf.rank[ry] {
		rx, ry = ry, rx
	}
	uf.parent[ry] = rx
	if uf.rank[rx] == uf.rank[ry] {
		uf.rank[rx]++
	}
	uf.count--
	return true
}

// Connected reports whether x and y are in the same set.
func (uf *UnionFind) Connected(x, y int) bool {
	return uf.Find(x) == uf.Find(y)
}

// Count returns the number of disjoint sets currently maintained.
func (uf *UnionFind) Count() int { return uf.count }

// Kruskal computes a minimum spanning tree (or minimum spanning forest, if the
// graph is disconnected) of an undirected graph using Kruskal's algorithm. It
// returns the chosen edges and their total weight. Edges are considered in
// non-decreasing weight order, breaking ties by (From, To) for determinism. It
// returns [ErrNotUndirected] for a directed graph.
func (g *Graph) Kruskal() (mst []Edge, total float64, err error) {
	if g.directed {
		return nil, 0, ErrNotUndirected
	}
	edges := g.Edges()
	sort.SliceStable(edges, func(i, j int) bool {
		return edges[i].Weight < edges[j].Weight
	})
	uf := NewUnionFind(g.Vertices())
	for _, e := range edges {
		if e.From == e.To {
			continue // self-loops never belong to a spanning tree
		}
		if uf.Union(e.From, e.To) {
			mst = append(mst, e)
			total += e.Weight
		}
	}
	return mst, total, nil
}

// graphPrimItem is an entry in the priority queue used by Prim's algorithm.
type graphPrimItem struct {
	from, to int
	weight   float64
}

// graphPrimPQ is a min-heap of candidate tree edges ordered by weight, breaking
// ties by (to, from) for determinism.
type graphPrimPQ []graphPrimItem

// Len reports the number of items in the queue, implementing heap.Interface.
func (pq graphPrimPQ) Len() int { return len(pq) }

// Less reports whether item i orders before item j, implementing heap.Interface.
// Items are ordered by ascending weight, breaking ties by ascending to and then
// ascending from vertex.
func (pq graphPrimPQ) Less(i, j int) bool {
	if pq[i].weight != pq[j].weight {
		return pq[i].weight < pq[j].weight
	}
	if pq[i].to != pq[j].to {
		return pq[i].to < pq[j].to
	}
	return pq[i].from < pq[j].from
}

// Swap exchanges items i and j, implementing heap.Interface.
func (pq graphPrimPQ) Swap(i, j int) { pq[i], pq[j] = pq[j], pq[i] }

// Push appends x to the queue, implementing heap.Interface.
func (pq *graphPrimPQ) Push(x any) { *pq = append(*pq, x.(graphPrimItem)) }

// Pop removes and returns the last item of the queue, implementing heap.Interface.
func (pq *graphPrimPQ) Pop() any {
	old := *pq
	n := len(old)
	it := old[n-1]
	*pq = old[:n-1]
	return it
}

// Prim computes a minimum spanning tree of the connected component containing
// start in an undirected graph using Prim's algorithm. It returns the chosen
// edges and their total weight. It returns [ErrNotUndirected] for a directed
// graph or [ErrVertexNotFound] if start is absent. To span every component of a
// possibly disconnected graph, use [Graph.Kruskal].
func (g *Graph) Prim(start int) (mst []Edge, total float64, err error) {
	if g.directed {
		return nil, 0, ErrNotUndirected
	}
	if !g.HasVertex(start) {
		return nil, 0, ErrVertexNotFound
	}
	inTree := map[int]bool{start: true}
	pq := &graphPrimPQ{}
	push := func(u int) {
		for _, v := range graphSortedKeys(g.adj[u]) {
			if !inTree[v] {
				heap.Push(pq, graphPrimItem{from: u, to: v, weight: g.adj[u][v]})
			}
		}
	}
	push(start)
	for pq.Len() > 0 {
		it := heap.Pop(pq).(graphPrimItem)
		if inTree[it.to] {
			continue
		}
		inTree[it.to] = true
		mst = append(mst, Edge{From: it.from, To: it.to, Weight: it.weight})
		total += it.weight
		push(it.to)
	}
	sort.Slice(mst, func(i, j int) bool {
		if mst[i].From != mst[j].From {
			return mst[i].From < mst[j].From
		}
		return mst[i].To < mst[j].To
	})
	return mst, total, nil
}
