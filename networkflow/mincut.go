package networkflow

import (
	"math"
	"sort"
)

// MinCutResult describes a cut: a bipartition of the vertices and the total
// capacity of the edges crossing from one side to the other.
type MinCutResult struct {
	// Value is the total capacity (or weight) crossing the cut.
	Value float64
	// S is the set of vertices on the source (or first) side of the cut.
	S []int
	// T is the complementary set of vertices.
	T []int
}

// Contains reports whether vertex v lies on the S side of the cut.
func (c *MinCutResult) Contains(v int) bool {
	for _, u := range c.S {
		if u == v {
			return true
		}
	}
	return false
}

// SizeS returns the number of vertices on the S side.
func (c *MinCutResult) SizeS() int { return len(c.S) }

// SizeT returns the number of vertices on the T side.
func (c *MinCutResult) SizeT() int { return len(c.T) }

// MinCutST computes a minimum s-t cut of a directed [FlowNetwork] via the
// max-flow/min-cut theorem: the cut value equals the maximum flow, and the S
// side is the set of vertices reachable from s in the residual network. The
// input network is left unchanged.
func MinCutST(g *FlowNetwork, s, t int) *MinCutResult {
	res := DinicResult(g, s, t)
	seen := res.Residual.ReachableInResidual(s)
	var S, T []int
	for v := 0; v < g.n; v++ {
		if seen[v] {
			S = append(S, v)
		} else {
			T = append(T, v)
		}
	}
	return &MinCutResult{Value: float64(res.Value), S: S, T: T}
}

// MinCutSTEdges returns the directed caller edges crossing the minimum s-t cut,
// whose capacities sum to the cut value.
func MinCutSTEdges(g *FlowNetwork, s, t int) []Edge {
	return DinicResult(g, s, t).MinCutEdges()
}

// MinCutSTValue returns just the value of the minimum s-t cut, equal to the
// maximum flow.
func MinCutSTValue(g *FlowNetwork, s, t int) int64 { return Dinic(g, s, t) }

// WeightedGraph is an undirected graph on n vertices with real, non-negative
// edge weights stored as a symmetric adjacency matrix. It is the input type for
// the Stoer-Wagner global minimum cut.
type WeightedGraph struct {
	n int
	w [][]float64
}

// NewWeightedGraph returns an empty undirected weighted graph on n vertices. It
// panics if n is negative.
func NewWeightedGraph(n int) *WeightedGraph {
	if n < 0 {
		panic("networkflow: negative vertex count")
	}
	w := make([][]float64, n)
	for i := range w {
		w[i] = make([]float64, n)
	}
	return &WeightedGraph{n: n, w: w}
}

// NumVertices returns the number of vertices.
func (g *WeightedGraph) NumVertices() int { return g.n }

// AddEdge adds weight to the undirected edge {u,v} (accumulating on parallel
// additions). Self-loops are ignored. It panics on an out-of-range vertex.
func (g *WeightedGraph) AddEdge(u, v int, weight float64) {
	if u < 0 || u >= g.n || v < 0 || v >= g.n {
		panic(ErrInvalidVertex)
	}
	if u == v {
		return
	}
	g.w[u][v] += weight
	g.w[v][u] += weight
}

// SetWeight sets the weight of the undirected edge {u,v} to weight, replacing
// any previous value. It panics on an out-of-range vertex.
func (g *WeightedGraph) SetWeight(u, v int, weight float64) {
	if u < 0 || u >= g.n || v < 0 || v >= g.n {
		panic(ErrInvalidVertex)
	}
	if u == v {
		return
	}
	g.w[u][v] = weight
	g.w[v][u] = weight
}

// Weight returns the weight of the undirected edge {u,v} (0 if absent).
func (g *WeightedGraph) Weight(u, v int) float64 { return g.w[u][v] }

// WeightedDegree returns the sum of the weights of all edges incident to v.
func (g *WeightedGraph) WeightedDegree(v int) float64 {
	var d float64
	for u := 0; u < g.n; u++ {
		d += g.w[v][u]
	}
	return d
}

// TotalWeight returns the sum of all edge weights (each undirected edge counted
// once).
func (g *WeightedGraph) TotalWeight() float64 {
	var s float64
	for i := 0; i < g.n; i++ {
		for j := i + 1; j < g.n; j++ {
			s += g.w[i][j]
		}
	}
	return s
}

// Clone returns a deep copy of the weighted graph.
func (g *WeightedGraph) Clone() *WeightedGraph {
	c := NewWeightedGraph(g.n)
	for i := 0; i < g.n; i++ {
		copy(c.w[i], g.w[i])
	}
	return c
}

// StoerWagner computes a global minimum cut of an undirected weighted graph
// with the Stoer-Wagner algorithm in O(V^3) time. It returns the cut value and
// the vertex partition. For fewer than two vertices the cut value is zero. Edge
// weights must be non-negative. The input graph is left unchanged.
func StoerWagner(g *WeightedGraph) *MinCutResult {
	n := g.n
	if n < 2 {
		var all []int
		for v := 0; v < n; v++ {
			all = append(all, v)
		}
		return &MinCutResult{Value: 0, S: all, T: nil}
	}

	// Working weight matrix over "merged" super-vertices.
	w := make([][]float64, n)
	for i := range w {
		w[i] = make([]float64, n)
		copy(w[i], g.w[i])
	}
	// vertices[i] is the list of original vertices merged into super-vertex i.
	vertices := make([][]int, n)
	for i := range vertices {
		vertices[i] = []int{i}
	}
	active := make([]bool, n)
	for i := range active {
		active[i] = true
	}

	best := math.Inf(1)
	var bestSet []int

	for phase := 0; phase < n-1; phase++ {
		added := make([]bool, n)
		weights := make([]float64, n)
		prev, last := -1, -1
		// Number of still-active super-vertices remaining this phase.
		remaining := 0
		for i := 0; i < n; i++ {
			if active[i] {
				remaining++
			}
		}
		for k := 0; k < remaining; k++ {
			sel := -1
			for i := 0; i < n; i++ {
				if active[i] && !added[i] && (sel == -1 || weights[i] > weights[sel]) {
					sel = i
				}
			}
			if sel == -1 {
				break
			}
			added[sel] = true
			prev = last
			last = sel
			for i := 0; i < n; i++ {
				if active[i] && !added[i] {
					weights[i] += w[sel][i]
				}
			}
		}
		// Cut-of-the-phase: separates last from the rest.
		if weights[last] < best {
			best = weights[last]
			bestSet = append([]int(nil), vertices[last]...)
		}
		// Merge last into prev.
		if prev >= 0 {
			vertices[prev] = append(vertices[prev], vertices[last]...)
			for i := 0; i < n; i++ {
				if active[i] && i != prev && i != last {
					w[prev][i] += w[last][i]
					w[i][prev] += w[i][last]
				}
			}
			active[last] = false
		}
	}

	inS := make(map[int]bool, len(bestSet))
	for _, v := range bestSet {
		inS[v] = true
	}
	var S, T []int
	for v := 0; v < n; v++ {
		if inS[v] {
			S = append(S, v)
		} else {
			T = append(T, v)
		}
	}
	sort.Ints(S)
	sort.Ints(T)
	return &MinCutResult{Value: best, S: S, T: T}
}

// GlobalMinCut is an alias for [StoerWagner] returning only the cut value.
func GlobalMinCut(g *WeightedGraph) float64 { return StoerWagner(g).Value }
