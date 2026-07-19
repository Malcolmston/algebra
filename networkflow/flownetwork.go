package networkflow

import (
	"errors"
	"fmt"
	"sort"
)

// Sentinel errors returned by constructors and validators.
var (
	// ErrInvalidVertex indicates a vertex index outside [0, NumVertices).
	ErrInvalidVertex = errors.New("networkflow: invalid vertex index")
	// ErrNegativeCapacity indicates an edge capacity below zero.
	ErrNegativeCapacity = errors.New("networkflow: negative capacity")
	// ErrSourceEqualsSink indicates that a source and sink coincide where a
	// distinct pair is required.
	ErrSourceEqualsSink = errors.New("networkflow: source equals sink")
	// ErrDimensionMismatch indicates a non-rectangular or wrongly shaped matrix.
	ErrDimensionMismatch = errors.New("networkflow: dimension mismatch")
	// ErrNoPerfectMatching indicates that no perfect assignment exists.
	ErrNoPerfectMatching = errors.New("networkflow: no perfect matching exists")
)

// Edge is a read-only view of a directed edge in a [FlowNetwork].
type Edge struct {
	// From is the tail vertex.
	From int
	// To is the head vertex.
	To int
	// Cap is the edge capacity.
	Cap int64
	// Flow is the flow currently routed along the edge.
	Flow int64
}

// Residual returns the remaining capacity Cap-Flow of the edge.
func (e Edge) Residual() int64 { return e.Cap - e.Flow }

// IsSaturated reports whether the edge carries flow equal to its capacity.
func (e Edge) IsSaturated() bool { return e.Flow >= e.Cap }

// String renders the edge as "from->to cap/flow".
func (e Edge) String() string {
	return fmt.Sprintf("%d->%d %d/%d", e.From, e.To, e.Flow, e.Cap)
}

// fedge is the internal storage for a directed edge. Edges are always added in
// pairs (forward, reverse) so that the reverse edge of index i is i^1.
type fedge struct {
	from int
	to   int
	cap  int64
	flow int64
}

// FlowNetwork is a directed graph with non-negative integer edge capacities.
// Each edge is stored with an anti-parallel residual edge, so a network can be
// used directly as its own residual graph while an algorithm mutates the flow
// values. Vertices are the integers 0..NumVertices-1.
type FlowNetwork struct {
	n     int
	edges []fedge
	adj   [][]int
}

// NewFlowNetwork returns an empty flow network on n vertices (labelled
// 0..n-1). It panics if n is negative.
func NewFlowNetwork(n int) *FlowNetwork {
	if n < 0 {
		panic("networkflow: negative vertex count")
	}
	return &FlowNetwork{n: n, adj: make([][]int, n)}
}

// NumVertices returns the number of vertices in the network.
func (g *FlowNetwork) NumVertices() int { return g.n }

// NumEdges returns the number of directed edges the caller added; residual
// twin edges are not counted.
func (g *FlowNetwork) NumEdges() int { return len(g.edges) / 2 }

// validVertex reports whether v is a legal vertex index.
func (g *FlowNetwork) validVertex(v int) bool { return v >= 0 && v < g.n }

// AddEdge adds a directed edge from->to with the given capacity and returns its
// edge id. It panics on an out-of-range vertex or a negative capacity. Parallel
// edges are allowed and are kept distinct.
func (g *FlowNetwork) AddEdge(from, to int, capacity int64) int {
	if !g.validVertex(from) || !g.validVertex(to) {
		panic(ErrInvalidVertex)
	}
	if capacity < 0 {
		panic(ErrNegativeCapacity)
	}
	id := len(g.edges)
	g.edges = append(g.edges, fedge{from, to, capacity, 0})
	g.edges = append(g.edges, fedge{to, from, 0, 0})
	g.adj[from] = append(g.adj[from], id)
	g.adj[to] = append(g.adj[to], id+1)
	return id
}

// AddUndirectedEdge adds an undirected edge {u,v} with the given capacity in
// both directions and returns its edge id. It panics on an out-of-range vertex
// or a negative capacity.
func (g *FlowNetwork) AddUndirectedEdge(u, v int, capacity int64) int {
	if !g.validVertex(u) || !g.validVertex(v) {
		panic(ErrInvalidVertex)
	}
	if capacity < 0 {
		panic(ErrNegativeCapacity)
	}
	id := len(g.edges)
	g.edges = append(g.edges, fedge{u, v, capacity, 0})
	g.edges = append(g.edges, fedge{v, u, capacity, 0})
	g.adj[u] = append(g.adj[u], id)
	g.adj[v] = append(g.adj[v], id+1)
	return id
}

// Edge returns a read-only view of the directed edge with the given id (as
// returned by [FlowNetwork.AddEdge]).
func (g *FlowNetwork) Edge(id int) Edge {
	e := g.edges[id]
	return Edge{e.from, e.to, e.cap, e.flow}
}

// Edges returns a read-only view of every caller-added directed edge, in the
// order the edges were added. Residual twin edges are omitted.
func (g *FlowNetwork) Edges() []Edge {
	out := make([]Edge, 0, len(g.edges)/2)
	for i := 0; i < len(g.edges); i += 2 {
		e := g.edges[i]
		out = append(out, Edge{e.from, e.to, e.cap, e.flow})
	}
	return out
}

// OutEdges returns the ids of every residual arc leaving v (both caller edges
// and reverse twins). Use [FlowNetwork.Edge] to inspect each.
func (g *FlowNetwork) OutEdges(v int) []int {
	ids := make([]int, len(g.adj[v]))
	copy(ids, g.adj[v])
	return ids
}

// Neighbors returns the distinct head vertices reachable from v along
// positive-capacity caller edges, sorted ascending.
func (g *FlowNetwork) Neighbors(v int) []int {
	seen := map[int]bool{}
	for _, id := range g.adj[v] {
		e := g.edges[id]
		if id%2 == 0 && e.cap > 0 {
			seen[e.to] = true
		}
	}
	out := make([]int, 0, len(seen))
	for u := range seen {
		out = append(out, u)
	}
	sort.Ints(out)
	return out
}

// OutDegree returns the number of caller edges leaving v.
func (g *FlowNetwork) OutDegree(v int) int {
	c := 0
	for _, id := range g.adj[v] {
		if id%2 == 0 {
			c++
		}
	}
	return c
}

// InDegree returns the number of caller edges entering v.
func (g *FlowNetwork) InDegree(v int) int {
	c := 0
	for _, id := range g.adj[v] {
		if id%2 == 1 {
			c++
		}
	}
	return c
}

// Capacity returns the total capacity of all caller edges from u to v.
func (g *FlowNetwork) Capacity(u, v int) int64 {
	var c int64
	for i := 0; i < len(g.edges); i += 2 {
		if g.edges[i].from == u && g.edges[i].to == v {
			c += g.edges[i].cap
		}
	}
	return c
}

// FlowOn returns the net flow of all caller edges from u to v.
func (g *FlowNetwork) FlowOn(u, v int) int64 {
	var f int64
	for i := 0; i < len(g.edges); i += 2 {
		if g.edges[i].from == u && g.edges[i].to == v {
			f += g.edges[i].flow
		}
	}
	return f
}

// TotalCapacity returns the sum of all caller-edge capacities.
func (g *FlowNetwork) TotalCapacity() int64 {
	var c int64
	for i := 0; i < len(g.edges); i += 2 {
		c += g.edges[i].cap
	}
	return c
}

// ResetFlow zeroes the flow on every edge, leaving capacities untouched.
func (g *FlowNetwork) ResetFlow() {
	for i := range g.edges {
		g.edges[i].flow = 0
	}
}

// Clone returns a deep copy of the network, including current flow values.
func (g *FlowNetwork) Clone() *FlowNetwork {
	c := &FlowNetwork{n: g.n, adj: make([][]int, g.n)}
	c.edges = make([]fedge, len(g.edges))
	copy(c.edges, g.edges)
	for v := range g.adj {
		c.adj[v] = make([]int, len(g.adj[v]))
		copy(c.adj[v], g.adj[v])
	}
	return c
}

// Validate reports the first structural problem it finds, or nil if the network
// is well formed (all endpoints in range and all capacities non-negative).
func (g *FlowNetwork) Validate() error {
	for i := 0; i < len(g.edges); i += 2 {
		e := g.edges[i]
		if !g.validVertex(e.from) || !g.validVertex(e.to) {
			return ErrInvalidVertex
		}
		if e.cap < 0 {
			return ErrNegativeCapacity
		}
	}
	return nil
}

// ResidualCapacity returns the residual capacity along the edge with the given
// id, i.e. Cap-Flow (which for a reverse twin equals the flow pushed forward).
func (g *FlowNetwork) ResidualCapacity(id int) int64 {
	return g.edges[id].cap - g.edges[id].flow
}

// ResidualEdges returns, for every residual arc with positive residual
// capacity, an [Edge] whose Cap field holds the residual capacity and whose
// Flow field is zero. Reverse twins with pushed flow are therefore included.
func (g *FlowNetwork) ResidualEdges() []Edge {
	var out []Edge
	for i := range g.edges {
		e := g.edges[i]
		if r := e.cap - e.flow; r > 0 {
			out = append(out, Edge{e.from, e.to, r, 0})
		}
	}
	return out
}

// ReachableInResidual returns a boolean slice marking every vertex reachable
// from src along arcs with positive residual capacity (a BFS of the residual
// graph).
func (g *FlowNetwork) ReachableInResidual(src int) []bool {
	seen := make([]bool, g.n)
	if !g.validVertex(src) {
		return seen
	}
	seen[src] = true
	queue := []int{src}
	for len(queue) > 0 {
		v := queue[0]
		queue = queue[1:]
		for _, id := range g.adj[v] {
			e := g.edges[id]
			if !seen[e.to] && e.cap-e.flow > 0 {
				seen[e.to] = true
				queue = append(queue, e.to)
			}
		}
	}
	return seen
}

// FlowValue returns the value of the current flow measured as the net flow out
// of source s (out-flow minus in-flow). For a valid s-t flow this equals the
// net flow into the sink.
func (g *FlowNetwork) FlowValue(s int) int64 {
	var v int64
	for i := 0; i < len(g.edges); i += 2 {
		if g.edges[i].from == s {
			v += g.edges[i].flow
		}
		if g.edges[i].to == s {
			v -= g.edges[i].flow
		}
	}
	return v
}

// IsFeasibleFlow reports whether the current flow respects capacities on every
// edge and conserves flow at every vertex other than s and t.
func (g *FlowNetwork) IsFeasibleFlow(s, t int) bool {
	net := make([]int64, g.n)
	for i := 0; i < len(g.edges); i += 2 {
		e := g.edges[i]
		if e.flow < 0 || e.flow > e.cap {
			return false
		}
		net[e.from] += e.flow
		net[e.to] -= e.flow
	}
	for v := 0; v < g.n; v++ {
		if v == s || v == t {
			continue
		}
		if net[v] != 0 {
			return false
		}
	}
	return true
}
