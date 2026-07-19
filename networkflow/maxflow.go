package networkflow

import "math"

// MaxFlowResult carries the outcome of a maximum-flow computation together with
// the residual network, from which the minimum s-t cut and a flow decomposition
// can be recovered.
type MaxFlowResult struct {
	// Value is the value of the maximum flow.
	Value int64
	// Source is the flow source.
	Source int
	// Sink is the flow sink.
	Sink int
	// Residual is the network with the maximum flow installed on its edges.
	Residual *FlowNetwork
}

// FlowValue returns the value of the maximum flow.
func (r *MaxFlowResult) FlowValue() int64 { return r.Value }

// FlowEdges returns every caller edge that carries strictly positive flow.
func (r *MaxFlowResult) FlowEdges() []Edge {
	var out []Edge
	for _, e := range r.Residual.Edges() {
		if e.Flow > 0 {
			out = append(out, e)
		}
	}
	return out
}

// SourceSide returns the set S of vertices reachable from the source in the
// residual network; it is the source side of the induced minimum cut.
func (r *MaxFlowResult) SourceSide() []int {
	seen := r.Residual.ReachableInResidual(r.Source)
	var s []int
	for v, ok := range seen {
		if ok {
			s = append(s, v)
		}
	}
	return s
}

// SinkSide returns the complement of [MaxFlowResult.SourceSide]: the vertices
// not reachable from the source in the residual network.
func (r *MaxFlowResult) SinkSide() []int {
	seen := r.Residual.ReachableInResidual(r.Source)
	var t []int
	for v, ok := range seen {
		if !ok {
			t = append(t, v)
		}
	}
	return t
}

// MinCutEdges returns the caller edges that cross the induced minimum cut, i.e.
// edges whose tail is on the source side and whose head is on the sink side.
// Their capacities sum to [MaxFlowResult.Value].
func (r *MaxFlowResult) MinCutEdges() []Edge {
	seen := r.Residual.ReachableInResidual(r.Source)
	var out []Edge
	for _, e := range r.Residual.Edges() {
		if seen[e.From] && !seen[e.To] && e.Cap > 0 {
			out = append(out, e)
		}
	}
	return out
}

// MinCutValue returns the capacity of the induced minimum cut, which equals the
// maximum-flow value.
func (r *MaxFlowResult) MinCutValue() int64 {
	var c int64
	for _, e := range r.MinCutEdges() {
		c += e.Cap
	}
	return c
}

// edmondsKarp augments along shortest paths (by edge count) found with BFS and
// returns the value of the maximum flow, mutating g in place.
func edmondsKarp(g *FlowNetwork, s, t int) int64 {
	var total int64
	for {
		prev := make([]int, g.n)
		for i := range prev {
			prev[i] = -1
		}
		prev[s] = -2
		queue := []int{s}
		for len(queue) > 0 && prev[t] == -1 {
			v := queue[0]
			queue = queue[1:]
			for _, id := range g.adj[v] {
				e := g.edges[id]
				if prev[e.to] == -1 && e.cap-e.flow > 0 {
					prev[e.to] = id
					queue = append(queue, e.to)
				}
			}
		}
		if prev[t] == -1 {
			break
		}
		bottleneck := int64(math.MaxInt64)
		for v := t; v != s; {
			id := prev[v]
			if r := g.edges[id].cap - g.edges[id].flow; r < bottleneck {
				bottleneck = r
			}
			v = g.edges[id].from
		}
		for v := t; v != s; {
			id := prev[v]
			g.edges[id].flow += bottleneck
			g.edges[id^1].flow -= bottleneck
			v = g.edges[id].from
		}
		total += bottleneck
	}
	return total
}

// dinic computes the maximum flow with Dinic's algorithm (BFS level graph plus
// DFS blocking flows), mutating g in place, and returns the flow value.
func dinic(g *FlowNetwork, s, t int) int64 {
	level := make([]int, g.n)
	iter := make([]int, g.n)

	bfs := func() bool {
		for i := range level {
			level[i] = -1
		}
		level[s] = 0
		queue := []int{s}
		for len(queue) > 0 {
			v := queue[0]
			queue = queue[1:]
			for _, id := range g.adj[v] {
				e := g.edges[id]
				if level[e.to] < 0 && e.cap-e.flow > 0 {
					level[e.to] = level[v] + 1
					queue = append(queue, e.to)
				}
			}
		}
		return level[t] >= 0
	}

	var dfs func(v int, pushed int64) int64
	dfs = func(v int, pushed int64) int64 {
		if v == t {
			return pushed
		}
		for ; iter[v] < len(g.adj[v]); iter[v]++ {
			id := g.adj[v][iter[v]]
			e := g.edges[id]
			if level[v]+1 == level[e.to] && e.cap-e.flow > 0 {
				d := dfs(e.to, min64(pushed, e.cap-e.flow))
				if d > 0 {
					g.edges[id].flow += d
					g.edges[id^1].flow -= d
					return d
				}
			}
		}
		return 0
	}

	var total int64
	for bfs() {
		for i := range iter {
			iter[i] = 0
		}
		for {
			f := dfs(s, math.MaxInt64)
			if f == 0 {
				break
			}
			total += f
		}
	}
	return total
}

// pushRelabel computes the maximum flow with the preflow-push method, mutating g
// in place. When highest is true it selects the active vertex of greatest
// height; otherwise it uses a FIFO queue. It returns the flow value.
func pushRelabel(g *FlowNetwork, s, t int, highest bool) int64 {
	if s == t {
		return 0
	}
	n := g.n
	height := make([]int, n)
	excess := make([]int64, n)
	cur := make([]int, n)
	height[s] = n

	push := func(id int) {
		e := &g.edges[id]
		d := excess[e.from]
		if r := e.cap - e.flow; r < d {
			d = r
		}
		e.flow += d
		g.edges[id^1].flow -= d
		excess[e.from] -= d
		excess[e.to] += d
	}
	relabel := func(v int) {
		mh := math.MaxInt
		for _, id := range g.adj[v] {
			e := g.edges[id]
			if e.cap-e.flow > 0 && height[e.to]+1 < mh {
				mh = height[e.to] + 1
			}
		}
		if mh < math.MaxInt {
			height[v] = mh
		}
	}

	// Initial preflow: saturate all source arcs.
	for _, id := range g.adj[s] {
		e := &g.edges[id]
		if d := e.cap - e.flow; d > 0 {
			e.flow += d
			g.edges[id^1].flow -= d
			excess[e.to] += d
			excess[s] -= d
		}
	}

	active := func(v int) bool { return v != s && v != t && excess[v] > 0 }
	discharge := func(v int) {
		for excess[v] > 0 {
			if cur[v] == len(g.adj[v]) {
				relabel(v)
				cur[v] = 0
				continue
			}
			id := g.adj[v][cur[v]]
			e := g.edges[id]
			if e.cap-e.flow > 0 && height[v] == height[e.to]+1 {
				push(id)
			} else {
				cur[v]++
			}
		}
	}

	if highest {
		for {
			v := -1
			for u := 0; u < n; u++ {
				if active(u) && (v == -1 || height[u] > height[v]) {
					v = u
				}
			}
			if v == -1 {
				break
			}
			discharge(v)
		}
	} else {
		inQ := make([]bool, n)
		var queue []int
		enqueue := func(v int) {
			if !inQ[v] && active(v) {
				inQ[v] = true
				queue = append(queue, v)
			}
		}
		for u := 0; u < n; u++ {
			enqueue(u)
		}
		for len(queue) > 0 {
			v := queue[0]
			queue = queue[1:]
			inQ[v] = false
			for excess[v] > 0 {
				if cur[v] == len(g.adj[v]) {
					relabel(v)
					cur[v] = 0
					continue
				}
				id := g.adj[v][cur[v]]
				e := g.edges[id]
				if e.cap-e.flow > 0 && height[v] == height[e.to]+1 {
					push(id)
					if e.to != s && e.to != t {
						enqueue(e.to)
					}
				} else {
					cur[v]++
				}
			}
		}
	}
	return excess[t]
}

// checkST panics with a descriptive error if s or t is out of range.
func checkST(g *FlowNetwork, s, t int) {
	if !g.validVertex(s) || !g.validVertex(t) {
		panic(ErrInvalidVertex)
	}
}

// EdmondsKarp returns the value of the maximum flow from s to t using the
// Edmonds-Karp algorithm. The input network is left unchanged. It panics if s
// or t is out of range.
func EdmondsKarp(g *FlowNetwork, s, t int) int64 {
	checkST(g, s, t)
	if s == t {
		return 0
	}
	return edmondsKarp(g.Clone(), s, t)
}

// EdmondsKarpResult runs Edmonds-Karp and returns a [MaxFlowResult] holding the
// flow value and the saturated residual network. The input is left unchanged.
func EdmondsKarpResult(g *FlowNetwork, s, t int) *MaxFlowResult {
	checkST(g, s, t)
	c := g.Clone()
	var v int64
	if s != t {
		v = edmondsKarp(c, s, t)
	}
	return &MaxFlowResult{Value: v, Source: s, Sink: t, Residual: c}
}

// Dinic returns the value of the maximum flow from s to t using Dinic's
// algorithm. The input network is left unchanged. It panics if s or t is out of
// range.
func Dinic(g *FlowNetwork, s, t int) int64 {
	checkST(g, s, t)
	if s == t {
		return 0
	}
	return dinic(g.Clone(), s, t)
}

// DinicResult runs Dinic's algorithm and returns a [MaxFlowResult] holding the
// flow value and the saturated residual network. The input is left unchanged.
func DinicResult(g *FlowNetwork, s, t int) *MaxFlowResult {
	checkST(g, s, t)
	c := g.Clone()
	var v int64
	if s != t {
		v = dinic(c, s, t)
	}
	return &MaxFlowResult{Value: v, Source: s, Sink: t, Residual: c}
}

// PushRelabel returns the value of the maximum flow from s to t using the FIFO
// preflow-push method. The input network is left unchanged. It panics if s or t
// is out of range.
func PushRelabel(g *FlowNetwork, s, t int) int64 {
	checkST(g, s, t)
	if s == t {
		return 0
	}
	return pushRelabel(g.Clone(), s, t, false)
}

// PushRelabelHighestLabel returns the maximum flow from s to t using preflow
// push with the highest-label selection rule. The input is left unchanged.
func PushRelabelHighestLabel(g *FlowNetwork, s, t int) int64 {
	checkST(g, s, t)
	if s == t {
		return 0
	}
	return pushRelabel(g.Clone(), s, t, true)
}

// PushRelabelResult runs FIFO preflow push and returns a [MaxFlowResult]
// holding the flow value and the saturated residual network.
func PushRelabelResult(g *FlowNetwork, s, t int) *MaxFlowResult {
	checkST(g, s, t)
	c := g.Clone()
	var v int64
	if s != t {
		v = pushRelabel(c, s, t, false)
	}
	return &MaxFlowResult{Value: v, Source: s, Sink: t, Residual: c}
}

// MaxFlow returns the value of the maximum flow from s to t using the default
// engine (Dinic). The input network is left unchanged.
func MaxFlow(g *FlowNetwork, s, t int) int64 { return Dinic(g, s, t) }

// MaxFlowResultOf runs the default engine (Dinic) and returns the full
// [MaxFlowResult]. The input network is left unchanged.
func MaxFlowResultOf(g *FlowNetwork, s, t int) *MaxFlowResult { return DinicResult(g, s, t) }

// min64 returns the smaller of two int64 values.
func min64(a, b int64) int64 {
	if a < b {
		return a
	}
	return b
}

// max64 returns the larger of two int64 values.
func max64(a, b int64) int64 {
	if a > b {
		return a
	}
	return b
}

// abs64 returns the absolute value of an int64.
func abs64(a int64) int64 {
	if a < 0 {
		return -a
	}
	return a
}
