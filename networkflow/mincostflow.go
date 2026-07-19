package networkflow

import "math"

// mcedge is the internal storage for a capacitated, costed directed edge.
type mcedge struct {
	from int
	to   int
	cap  int64
	flow int64
	cost int64
}

// MinCostNetwork is a directed graph with non-negative integer capacities and
// integer per-unit costs on each edge. It supports minimum-cost flow via
// successive shortest augmenting paths. Vertices are 0..NumVertices-1.
type MinCostNetwork struct {
	n     int
	edges []mcedge
	adj   [][]int
}

// NewMinCostNetwork returns an empty min-cost network on n vertices. It panics
// if n is negative.
func NewMinCostNetwork(n int) *MinCostNetwork {
	if n < 0 {
		panic("networkflow: negative vertex count")
	}
	return &MinCostNetwork{n: n, adj: make([][]int, n)}
}

// NumVertices returns the number of vertices.
func (g *MinCostNetwork) NumVertices() int { return g.n }

// NumEdges returns the number of caller edges (residual twins excluded).
func (g *MinCostNetwork) NumEdges() int { return len(g.edges) / 2 }

// AddEdge adds a directed edge from->to with the given capacity and per-unit
// cost, returning its edge id. The reverse residual twin carries the negated
// cost. It panics on an out-of-range vertex or negative capacity.
func (g *MinCostNetwork) AddEdge(from, to int, capacity, cost int64) int {
	if from < 0 || from >= g.n || to < 0 || to >= g.n {
		panic(ErrInvalidVertex)
	}
	if capacity < 0 {
		panic(ErrNegativeCapacity)
	}
	id := len(g.edges)
	g.edges = append(g.edges, mcedge{from, to, capacity, 0, cost})
	g.edges = append(g.edges, mcedge{to, from, 0, 0, -cost})
	g.adj[from] = append(g.adj[from], id)
	g.adj[to] = append(g.adj[to], id+1)
	return id
}

// Edge returns a read-only [Edge] view (capacity and flow) of the caller edge
// with the given id. Use [MinCostNetwork.EdgeCost] for its cost.
func (g *MinCostNetwork) Edge(id int) Edge {
	e := g.edges[id]
	return Edge{e.from, e.to, e.cap, e.flow}
}

// EdgeCost returns the per-unit cost of the edge with the given id.
func (g *MinCostNetwork) EdgeCost(id int) int64 { return g.edges[id].cost }

// Edges returns read-only [Edge] views of every caller edge in insertion order.
func (g *MinCostNetwork) Edges() []Edge {
	out := make([]Edge, 0, len(g.edges)/2)
	for i := 0; i < len(g.edges); i += 2 {
		e := g.edges[i]
		out = append(out, Edge{e.from, e.to, e.cap, e.flow})
	}
	return out
}

// ResetFlow zeroes the flow on every edge.
func (g *MinCostNetwork) ResetFlow() {
	for i := range g.edges {
		g.edges[i].flow = 0
	}
}

// Clone returns a deep copy of the network, including current flow values.
func (g *MinCostNetwork) Clone() *MinCostNetwork {
	c := &MinCostNetwork{n: g.n, adj: make([][]int, g.n)}
	c.edges = make([]mcedge, len(g.edges))
	copy(c.edges, g.edges)
	for v := range g.adj {
		c.adj[v] = make([]int, len(g.adj[v]))
		copy(c.adj[v], g.adj[v])
	}
	return c
}

// TotalCost returns the total cost of the current flow, summed as flow*cost
// over caller edges.
func (g *MinCostNetwork) TotalCost() int64 {
	var c int64
	for i := 0; i < len(g.edges); i += 2 {
		c += g.edges[i].flow * g.edges[i].cost
	}
	return c
}

// pathCost sums the per-unit cost along the augmenting path recorded in prevE
// from t back to s.
func (g *MinCostNetwork) pathCost(prevE []int, s, t int) int64 {
	var c int64
	for v := t; v != s; {
		id := prevE[v]
		c += g.edges[id].cost
		v = g.edges[id].from
	}
	return c
}

// pushAlong augments the path recorded in prevE by amount, updating both twins.
func (g *MinCostNetwork) pushAlong(prevE []int, s, t int, amount int64) {
	for v := t; v != s; {
		id := prevE[v]
		g.edges[id].flow += amount
		g.edges[id^1].flow -= amount
		v = g.edges[id].from
	}
}

// bottleneck returns the minimum residual capacity along the path in prevE,
// capped at limit.
func (g *MinCostNetwork) bottleneck(prevE []int, s, t int, limit int64) int64 {
	push := limit
	for v := t; v != s; {
		id := prevE[v]
		if r := g.edges[id].cap - g.edges[id].flow; r < push {
			push = r
		}
		v = g.edges[id].from
	}
	return push
}

// mcmfSPFA runs successive shortest augmenting paths, using SPFA
// (queue-based Bellman-Ford) shortest paths, until either maxFlow units are
// shipped or the sink is unreachable. It mutates g and returns (flow, cost).
func mcmfSPFA(g *MinCostNetwork, s, t int, maxFlow int64) (int64, int64) {
	var flow, cost int64
	n := g.n
	for flow < maxFlow {
		dist := make([]int64, n)
		inq := make([]bool, n)
		prevE := make([]int, n)
		for i := range dist {
			dist[i] = math.MaxInt64
			prevE[i] = -1
		}
		dist[s] = 0
		queue := []int{s}
		inq[s] = true
		for len(queue) > 0 {
			v := queue[0]
			queue = queue[1:]
			inq[v] = false
			for _, id := range g.adj[v] {
				e := g.edges[id]
				if e.cap-e.flow > 0 && dist[v] != math.MaxInt64 && dist[v]+e.cost < dist[e.to] {
					dist[e.to] = dist[v] + e.cost
					prevE[e.to] = id
					if !inq[e.to] {
						inq[e.to] = true
						queue = append(queue, e.to)
					}
				}
			}
		}
		if dist[t] == math.MaxInt64 {
			break
		}
		push := g.bottleneck(prevE, s, t, maxFlow-flow)
		g.pushAlong(prevE, s, t, push)
		flow += push
		cost += push * g.pathCost(prevE, s, t)
	}
	return flow, cost
}

// mcmfDijkstra runs successive shortest augmenting paths using Johnson
// potentials and Dijkstra with reduced costs. An initial Bellman-Ford pass sets
// the potentials so that negative edge costs are handled. It mutates g and
// returns (flow, cost).
func mcmfDijkstra(g *MinCostNetwork, s, t int, maxFlow int64) (int64, int64) {
	n := g.n
	h := make([]int64, n)
	for i := range h {
		h[i] = math.MaxInt64
	}
	h[s] = 0
	for iter := 0; iter < n-1; iter++ {
		changed := false
		for v := 0; v < n; v++ {
			if h[v] == math.MaxInt64 {
				continue
			}
			for _, id := range g.adj[v] {
				e := g.edges[id]
				if e.cap-e.flow > 0 && h[v]+e.cost < h[e.to] {
					h[e.to] = h[v] + e.cost
					changed = true
				}
			}
		}
		if !changed {
			break
		}
	}
	for i := range h {
		if h[i] == math.MaxInt64 {
			h[i] = 0
		}
	}

	var flow, cost int64
	for flow < maxFlow {
		dist := make([]int64, n)
		prevE := make([]int, n)
		visited := make([]bool, n)
		for i := range dist {
			dist[i] = math.MaxInt64
			prevE[i] = -1
		}
		dist[s] = 0
		for {
			u := -1
			for i := 0; i < n; i++ {
				if !visited[i] && dist[i] != math.MaxInt64 && (u == -1 || dist[i] < dist[u]) {
					u = i
				}
			}
			if u == -1 {
				break
			}
			visited[u] = true
			for _, id := range g.adj[u] {
				e := g.edges[id]
				if e.cap-e.flow > 0 && !visited[e.to] {
					nd := dist[u] + e.cost + h[u] - h[e.to]
					if nd < dist[e.to] {
						dist[e.to] = nd
						prevE[e.to] = id
					}
				}
			}
		}
		if dist[t] == math.MaxInt64 {
			break
		}
		for i := 0; i < n; i++ {
			if dist[i] != math.MaxInt64 {
				h[i] += dist[i]
			}
		}
		push := g.bottleneck(prevE, s, t, maxFlow-flow)
		g.pushAlong(prevE, s, t, push)
		flow += push
		cost += push * g.pathCost(prevE, s, t)
	}
	return flow, cost
}

// MCMFResult carries the outcome of a minimum-cost flow computation.
type MCMFResult struct {
	// Flow is the value of the flow shipped from source to sink.
	Flow int64
	// Cost is the total cost of that flow.
	Cost int64
	// Source is the flow source.
	Source int
	// Sink is the flow sink.
	Sink int
	// Residual is the network with the computed flow installed.
	Residual *MinCostNetwork
}

// FlowValue returns the value of the flow.
func (r *MCMFResult) FlowValue() int64 { return r.Flow }

// TotalCost returns the total cost of the flow.
func (r *MCMFResult) TotalCost() int64 { return r.Cost }

// FlowEdges returns every caller edge carrying strictly positive flow.
func (r *MCMFResult) FlowEdges() []Edge {
	var out []Edge
	for _, e := range r.Residual.Edges() {
		if e.Flow > 0 {
			out = append(out, e)
		}
	}
	return out
}

// checkMCST panics if s or t is out of range for g.
func checkMCST(g *MinCostNetwork, s, t int) {
	if s < 0 || s >= g.n || t < 0 || t >= g.n {
		panic(ErrInvalidVertex)
	}
}

// MinCostMaxFlow returns the maximum flow from s to t of least total cost,
// found by successive shortest paths with SPFA (Bellman-Ford). It tolerates
// negative edge costs but assumes no negative-cost cycle exists in the residual
// network. The input is left unchanged.
func MinCostMaxFlow(g *MinCostNetwork, s, t int) (flow, cost int64) {
	checkMCST(g, s, t)
	if s == t {
		return 0, 0
	}
	return mcmfSPFA(g.Clone(), s, t, math.MaxInt64)
}

// MinCostFlow ships up to want units from s to t at least total cost, returning
// the flow actually shipped and its cost. If less than want units can be sent
// it returns the maximum feasible flow. The input is left unchanged.
func MinCostFlow(g *MinCostNetwork, s, t int, want int64) (flow, cost int64) {
	checkMCST(g, s, t)
	if s == t || want <= 0 {
		return 0, 0
	}
	return mcmfSPFA(g.Clone(), s, t, want)
}

// MinCostMaxFlowDijkstra returns the minimum-cost maximum flow using Johnson
// potentials with Dijkstra. It is typically faster than [MinCostMaxFlow] on
// large graphs; the initial potentials are set with Bellman-Ford so negative
// edge costs are handled. The input is left unchanged.
func MinCostMaxFlowDijkstra(g *MinCostNetwork, s, t int) (flow, cost int64) {
	checkMCST(g, s, t)
	if s == t {
		return 0, 0
	}
	return mcmfDijkstra(g.Clone(), s, t, math.MaxInt64)
}

// MinCostMaxFlowResult runs SPFA-based min-cost max-flow and returns the full
// [MCMFResult], including the residual network.
func MinCostMaxFlowResult(g *MinCostNetwork, s, t int) *MCMFResult {
	checkMCST(g, s, t)
	c := g.Clone()
	var flow, cost int64
	if s != t {
		flow, cost = mcmfSPFA(c, s, t, math.MaxInt64)
	}
	return &MCMFResult{Flow: flow, Cost: cost, Source: s, Sink: t, Residual: c}
}
