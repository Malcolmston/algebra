package networkflow

import "math"

// FlowPath is a single source-to-sink path (or cycle) carrying a fixed amount
// of flow, produced by [FlowDecomposition].
type FlowPath struct {
	// Vertices lists the vertices along the path in order. For a cycle the
	// first and last vertices coincide.
	Vertices []int
	// Flow is the amount of flow carried along the path.
	Flow int64
	// Cycle reports whether the path is a cycle rather than an s-t path.
	Cycle bool
}

// Length returns the number of edges in the path.
func (p FlowPath) Length() int {
	if len(p.Vertices) == 0 {
		return 0
	}
	return len(p.Vertices) - 1
}

// FlowDecomposition decomposes the flow currently installed in g into
// source-to-sink paths and, if any circulation remains, cycles. The sum of the
// path flows equals the value of the flow. The input network is not modified.
func FlowDecomposition(g *FlowNetwork, s, t int) []FlowPath {
	// Work on positive net flow of caller edges only.
	type arc struct {
		to   int
		flow int64
	}
	out := make([][]arc, g.n)
	for i := 0; i < len(g.edges); i += 2 {
		e := g.edges[i]
		if e.flow > 0 {
			out[e.from] = append(out[e.from], arc{e.to, e.flow})
		}
	}
	ptr := make([]int, g.n)

	// residualOut returns the next arc index from v with remaining flow, or -1.
	nextArc := func(v int) int {
		for ptr[v] < len(out[v]) {
			if out[v][ptr[v]].flow > 0 {
				return ptr[v]
			}
			ptr[v]++
		}
		return -1
	}

	var paths []FlowPath
	// Extract s-t paths.
	for {
		if nextArc(s) == -1 {
			break
		}
		path := []int{s}
		onPath := make([]int, g.n)
		for i := range onPath {
			onPath[i] = -1
		}
		onPath[s] = 0
		v := s
		bottleneck := int64(math.MaxInt64)
		cyclePaths := false
		for v != t {
			ai := nextArc(v)
			if ai == -1 {
				cyclePaths = true
				break
			}
			a := out[v][ai]
			if a.flow < bottleneck {
				bottleneck = a.flow
			}
			path = append(path, a.to)
			if onPath[a.to] != -1 {
				// Found a cycle within the walk; handle below.
				cyclePaths = true
				break
			}
			onPath[a.to] = len(path) - 1
			v = a.to
		}
		if cyclePaths {
			break
		}
		// Subtract bottleneck along the path.
		for i := 0; i+1 < len(path); i++ {
			u := path[i]
			for j := range out[u] {
				if out[u][j].to == path[i+1] && out[u][j].flow > 0 {
					out[u][j].flow -= bottleneck
					break
				}
			}
		}
		paths = append(paths, FlowPath{Vertices: path, Flow: bottleneck})
		for i := range ptr {
			ptr[i] = 0
		}
	}

	// Extract remaining circulation as cycles.
	for v := 0; v < g.n; v++ {
		for {
			ai := nextArc(v)
			if ai == -1 {
				break
			}
			// Walk until we revisit a vertex.
			pos := make(map[int]int)
			walk := []int{v}
			pos[v] = 0
			cur := v
			ok := true
			for {
				a := nextArc(cur)
				if a == -1 {
					ok = false
					break
				}
				nxt := out[cur][a].to
				if p, seen := pos[nxt]; seen {
					// cycle from p..end
					cyc := append([]int{}, walk[p:]...)
					cyc = append(cyc, nxt)
					var bn int64 = math.MaxInt64
					for i := 0; i+1 < len(cyc); i++ {
						uu := cyc[i]
						for j := range out[uu] {
							if out[uu][j].to == cyc[i+1] && out[uu][j].flow > 0 {
								if out[uu][j].flow < bn {
									bn = out[uu][j].flow
								}
								break
							}
						}
					}
					for i := 0; i+1 < len(cyc); i++ {
						uu := cyc[i]
						for j := range out[uu] {
							if out[uu][j].to == cyc[i+1] && out[uu][j].flow > 0 {
								out[uu][j].flow -= bn
								break
							}
						}
					}
					paths = append(paths, FlowPath{Vertices: cyc, Flow: bn, Cycle: true})
					break
				}
				walk = append(walk, nxt)
				pos[nxt] = len(walk) - 1
				cur = nxt
			}
			if !ok {
				break
			}
			for i := range ptr {
				ptr[i] = 0
			}
		}
	}
	return paths
}

// GomoryHuTree is a weighted tree on the same vertex set as a source
// [FlowNetwork] such that, for every pair of vertices u and v, the minimum
// weight of an edge on the tree path from u to v equals the minimum u-v cut in
// the original graph. It is produced by [GomoryHu].
type GomoryHuTree struct {
	n      int
	parent []int
	weight []int64
}

// NumVertices returns the number of vertices in the tree.
func (t *GomoryHuTree) NumVertices() int { return t.n }

// Parent returns the tree parent of v, or -1 for the root (vertex 0).
func (t *GomoryHuTree) Parent(v int) int { return t.parent[v] }

// ParentWeight returns the weight of the tree edge from v to its parent, i.e.
// the minimum cut separating v from its parent. It is 0 for the root.
func (t *GomoryHuTree) ParentWeight(v int) int64 { return t.weight[v] }

// TreeEdges returns the n-1 tree edges as {u, v, weight} triples with u the
// child and v its parent.
func (t *GomoryHuTree) TreeEdges() [][3]int64 {
	var out [][3]int64
	for v := 1; v < t.n; v++ {
		out = append(out, [3]int64{int64(v), int64(t.parent[v]), t.weight[v]})
	}
	return out
}

// MinCut returns the minimum u-v cut value, computed as the minimum edge weight
// on the tree path between u and v. For u == v it returns 0.
func (t *GomoryHuTree) MinCut(u, v int) int64 {
	if u == v {
		return 0
	}
	// Depth of each node.
	depth := func(x int) int {
		d := 0
		for x != 0 {
			x = t.parent[x]
			d++
		}
		return d
	}
	du, dv := depth(u), depth(v)
	best := int64(math.MaxInt64)
	for du > dv {
		if t.weight[u] < best {
			best = t.weight[u]
		}
		u = t.parent[u]
		du--
	}
	for dv > du {
		if t.weight[v] < best {
			best = t.weight[v]
		}
		v = t.parent[v]
		dv--
	}
	for u != v {
		if t.weight[u] < best {
			best = t.weight[u]
		}
		if t.weight[v] < best {
			best = t.weight[v]
		}
		u = t.parent[u]
		v = t.parent[v]
	}
	return best
}

// AllPairsMinCut returns the full symmetric matrix of pairwise minimum cut
// values derived from the tree.
func (t *GomoryHuTree) AllPairsMinCut() [][]int64 {
	m := make([][]int64, t.n)
	for i := range m {
		m[i] = make([]int64, t.n)
	}
	for i := 0; i < t.n; i++ {
		for j := i + 1; j < t.n; j++ {
			c := t.MinCut(i, j)
			m[i][j] = c
			m[j][i] = c
		}
	}
	return m
}

// GomoryHu builds a Gomory-Hu cut tree of an undirected flow network using
// Gusfield's algorithm, which performs n-1 maximum-flow computations. The
// network should have been built with [FlowNetwork.AddUndirectedEdge] so that
// capacities are symmetric. The input is left unchanged. A network with fewer
// than two vertices yields a trivial tree.
func GomoryHu(g *FlowNetwork) *GomoryHuTree {
	n := g.n
	parent := make([]int, n)
	weight := make([]int64, n)
	for i := range parent {
		parent[i] = 0
	}
	for i := 1; i < n; i++ {
		p := parent[i]
		res := DinicResult(g, i, p)
		f := res.Value
		seen := res.Residual.ReachableInResidual(i)
		for j := i + 1; j < n; j++ {
			if seen[j] && parent[j] == p {
				parent[j] = i
			}
		}
		weight[i] = f
		if seen[parent[p]] {
			// Gusfield's re-rooting fix-up: swap i and its parent p.
			parent[i] = parent[p]
			parent[p] = i
			weight[i] = weight[p]
			weight[p] = f
		}
	}
	return &GomoryHuTree{n: n, parent: parent, weight: weight}
}
