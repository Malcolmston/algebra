package networkflow

import (
	"fmt"
	"math"
	"testing"
)

const ftol = 1e-9

// clrsNetwork builds the classic 6-vertex maximum-flow instance from CLRS whose
// maximum flow from 0 to 5 is 23.
func clrsNetwork() *FlowNetwork {
	g := NewFlowNetwork(6)
	g.AddEdge(0, 1, 16)
	g.AddEdge(0, 2, 13)
	g.AddEdge(1, 2, 10)
	g.AddEdge(2, 1, 4)
	g.AddEdge(1, 3, 12)
	g.AddEdge(3, 2, 9)
	g.AddEdge(2, 4, 14)
	g.AddEdge(4, 3, 7)
	g.AddEdge(3, 5, 20)
	g.AddEdge(4, 5, 4)
	return g
}

func TestMaxFlowEngines(t *testing.T) {
	engines := []struct {
		name string
		fn   func(*FlowNetwork, int, int) int64
	}{
		{"EdmondsKarp", EdmondsKarp},
		{"Dinic", Dinic},
		{"PushRelabel", PushRelabel},
		{"PushRelabelHighestLabel", PushRelabelHighestLabel},
		{"MaxFlow", MaxFlow},
	}
	for _, e := range engines {
		g := clrsNetwork()
		if got := e.fn(g, 0, 5); got != 23 {
			t.Errorf("%s: max flow = %d, want 23", e.name, got)
		}
		// The input network must be untouched.
		if fv := g.FlowValue(0); fv != 0 {
			t.Errorf("%s: input network mutated, flow value = %d", e.name, fv)
		}
	}
}

func TestMaxFlowSmallCases(t *testing.T) {
	tests := []struct {
		name  string
		build func() *FlowNetwork
		s, t  int
		want  int64
	}{
		{
			name: "single edge",
			build: func() *FlowNetwork {
				g := NewFlowNetwork(2)
				g.AddEdge(0, 1, 5)
				return g
			},
			s: 0, t: 1, want: 5,
		},
		{
			name: "series bottleneck",
			build: func() *FlowNetwork {
				g := NewFlowNetwork(3)
				g.AddEdge(0, 1, 10)
				g.AddEdge(1, 2, 3)
				return g
			},
			s: 0, t: 2, want: 3,
		},
		{
			name: "parallel edges",
			build: func() *FlowNetwork {
				g := NewFlowNetwork(2)
				g.AddEdge(0, 1, 3)
				g.AddEdge(0, 1, 4)
				return g
			},
			s: 0, t: 1, want: 7,
		},
		{
			name: "diamond",
			build: func() *FlowNetwork {
				g := NewFlowNetwork(4)
				g.AddEdge(0, 1, 3)
				g.AddEdge(0, 2, 2)
				g.AddEdge(1, 3, 2)
				g.AddEdge(2, 3, 3)
				return g
			},
			s: 0, t: 3, want: 4,
		},
		{
			name: "disconnected",
			build: func() *FlowNetwork {
				g := NewFlowNetwork(4)
				g.AddEdge(0, 1, 5)
				g.AddEdge(2, 3, 5)
				return g
			},
			s: 0, t: 3, want: 0,
		},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			for _, fn := range []func(*FlowNetwork, int, int) int64{EdmondsKarp, Dinic, PushRelabel, PushRelabelHighestLabel} {
				if got := fn(tc.build(), tc.s, tc.t); got != tc.want {
					t.Errorf("max flow = %d, want %d", got, tc.want)
				}
			}
		})
	}
}

func TestMaxFlowResultMinCut(t *testing.T) {
	res := DinicResult(clrsNetwork(), 0, 5)
	if res.Value != 23 {
		t.Fatalf("value = %d, want 23", res.Value)
	}
	if res.MinCutValue() != 23 {
		t.Errorf("min cut value = %d, want 23", res.MinCutValue())
	}
	// Feasibility of the produced flow.
	if !res.Residual.IsFeasibleFlow(0, 5) {
		t.Errorf("produced flow is infeasible")
	}
	// Source side and sink side partition all vertices.
	if len(res.SourceSide())+len(res.SinkSide()) != 6 {
		t.Errorf("cut sides do not partition vertices")
	}
	// Min cut edges capacity sums to the flow value.
	var sum int64
	for _, e := range res.MinCutEdges() {
		sum += e.Cap
	}
	if sum != 23 {
		t.Errorf("min cut edge capacities sum to %d, want 23", sum)
	}
}

func TestMinCutST(t *testing.T) {
	c := MinCutST(clrsNetwork(), 0, 5)
	if math.Abs(c.Value-23) > ftol {
		t.Errorf("min s-t cut = %v, want 23", c.Value)
	}
	if c.SizeS()+c.SizeT() != 6 {
		t.Errorf("partition size mismatch")
	}
	if !c.Contains(0) {
		t.Errorf("source should be on S side")
	}
}

func TestStoerWagner(t *testing.T) {
	tests := []struct {
		name  string
		build func() *WeightedGraph
		want  float64
	}{
		{
			name: "triangle unit weights",
			build: func() *WeightedGraph {
				g := NewWeightedGraph(3)
				g.AddEdge(0, 1, 1)
				g.AddEdge(1, 2, 1)
				g.AddEdge(0, 2, 1)
				return g
			},
			want: 2,
		},
		{
			name: "path graph",
			build: func() *WeightedGraph {
				g := NewWeightedGraph(4)
				g.AddEdge(0, 1, 5)
				g.AddEdge(1, 2, 1)
				g.AddEdge(2, 3, 5)
				return g
			},
			want: 1,
		},
		{
			name: "stoer-wagner paper",
			build: func() *WeightedGraph {
				g := NewWeightedGraph(8)
				edges := [][3]int{
					{0, 1, 2}, {0, 4, 3}, {1, 2, 3}, {1, 4, 2}, {1, 5, 2},
					{2, 3, 4}, {2, 6, 2}, {3, 6, 2}, {3, 7, 2}, {4, 5, 3},
					{5, 6, 1}, {6, 7, 3},
				}
				for _, e := range edges {
					g.AddEdge(e[0], e[1], float64(e[2]))
				}
				return g
			},
			want: 4,
		},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			c := StoerWagner(tc.build())
			if math.Abs(c.Value-tc.want) > ftol {
				t.Errorf("global min cut = %v, want %v", c.Value, tc.want)
			}
			if c.SizeS() == 0 || c.SizeT() == 0 {
				t.Errorf("degenerate partition: |S|=%d |T|=%d", c.SizeS(), c.SizeT())
			}
		})
	}
}

func TestEdgeConnectivity(t *testing.T) {
	// A 4-cycle with unit weights has edge connectivity 2.
	g := NewWeightedGraph(4)
	g.AddEdge(0, 1, 1)
	g.AddEdge(1, 2, 1)
	g.AddEdge(2, 3, 1)
	g.AddEdge(3, 0, 1)
	if got := EdgeConnectivity(g); math.Abs(got-2) > ftol {
		t.Errorf("edge connectivity = %v, want 2", got)
	}
}

func TestBipartiteMatching(t *testing.T) {
	tests := []struct {
		name  string
		build func() *BipartiteGraph
		want  int
	}{
		{
			name: "perfect 3x3",
			build: func() *BipartiteGraph {
				g := NewBipartiteGraph(3, 3)
				g.AddEdge(0, 0)
				g.AddEdge(0, 1)
				g.AddEdge(1, 0)
				g.AddEdge(1, 2)
				g.AddEdge(2, 1)
				return g
			},
			want: 3,
		},
		{
			name: "star limits to 1",
			build: func() *BipartiteGraph {
				g := NewBipartiteGraph(3, 3)
				g.AddEdge(0, 0)
				g.AddEdge(1, 0)
				g.AddEdge(2, 0)
				return g
			},
			want: 1,
		},
		{
			name: "path",
			build: func() *BipartiteGraph {
				g := NewBipartiteGraph(2, 2)
				g.AddEdge(0, 0)
				g.AddEdge(1, 0)
				g.AddEdge(1, 1)
				return g
			},
			want: 2,
		},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			hk := HopcroftKarp(tc.build())
			ku := KuhnMatching(tc.build())
			if hk.Size != tc.want {
				t.Errorf("Hopcroft-Karp size = %d, want %d", hk.Size, tc.want)
			}
			if ku.Size != tc.want {
				t.Errorf("Kuhn size = %d, want %d", ku.Size, tc.want)
			}
			// Konig: min vertex cover size equals matching size.
			cl, cr := MinimumVertexCover(tc.build())
			if len(cl)+len(cr) != tc.want {
				t.Errorf("Konig cover size = %d, want %d", len(cl)+len(cr), tc.want)
			}
			// Every edge must be covered.
			for _, e := range tc.build().Edges() {
				covered := false
				for _, u := range cl {
					if u == e[0] {
						covered = true
					}
				}
				for _, v := range cr {
					if v == e[1] {
						covered = true
					}
				}
				if !covered {
					t.Errorf("edge %v not covered by vertex cover", e)
				}
			}
		})
	}
}

func TestMatchingValidity(t *testing.T) {
	g := NewBipartiteGraph(3, 3)
	g.AddEdge(0, 0)
	g.AddEdge(0, 1)
	g.AddEdge(1, 0)
	g.AddEdge(1, 2)
	g.AddEdge(2, 1)
	m := HopcroftKarp(g)
	seenR := map[int]bool{}
	for _, p := range m.Pairs() {
		if !g.HasEdge(p[0], p[1]) {
			t.Errorf("matched pair %v is not an edge", p)
		}
		if seenR[p[1]] {
			t.Errorf("right vertex %d matched twice", p[1])
		}
		seenR[p[1]] = true
	}
	if !m.IsPerfect() {
		t.Errorf("expected perfect matching")
	}
}

func TestHungarianAssignment(t *testing.T) {
	tests := []struct {
		name     string
		cost     [][]float64
		wantCost float64
	}{
		{
			name:     "3x3",
			cost:     [][]float64{{4, 1, 3}, {2, 0, 5}, {3, 2, 2}},
			wantCost: 5,
		},
		{
			name:     "identity-ish",
			cost:     [][]float64{{1, 4, 5}, {5, 1, 4}, {4, 5, 1}},
			wantCost: 3,
		},
		{
			name:     "2x2",
			cost:     [][]float64{{1, 2}, {2, 1}},
			wantCost: 2,
		},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			m, err := CostMatrixFrom(tc.cost)
			if err != nil {
				t.Fatal(err)
			}
			res := HungarianMinCost(m)
			if math.Abs(res.Cost-tc.wantCost) > ftol {
				t.Errorf("min assignment cost = %v, want %v", res.Cost, tc.wantCost)
			}
			if !IsValidAssignment(m, res.Assignment) {
				t.Errorf("assignment %v is not valid", res.Assignment)
			}
			if math.Abs(AssignmentCost(m, res.Assignment)-tc.wantCost) > ftol {
				t.Errorf("recomputed cost mismatch")
			}
		})
	}
}

func TestMaxWeightAssignment(t *testing.T) {
	// Maximizing over this matrix picks the anti-diagonal 3+3+3 = 9.
	w := [][]float64{{1, 1, 3}, {1, 3, 1}, {3, 1, 1}}
	m, _ := CostMatrixFrom(w)
	res := MaxWeightAssignment(m)
	if math.Abs(res.Cost-9) > ftol {
		t.Errorf("max weight = %v, want 9", res.Cost)
	}
}

func TestRectangularAssignment(t *testing.T) {
	// 2 rows, 3 columns: assign each row to a distinct column minimizing cost.
	cost := [][]float64{{9, 2, 7}, {6, 4, 3}}
	m, _ := CostMatrixFrom(cost)
	res := HungarianMinCost(m)
	// Best: row0->col1 (2), row1->col2 (3) = 5.
	if math.Abs(res.Cost-5) > ftol {
		t.Errorf("rectangular min cost = %v, want 5", res.Cost)
	}
	if res.Size() != 2 {
		t.Errorf("assigned %d rows, want 2", res.Size())
	}

	// 3 rows, 2 columns: only two rows can be assigned.
	cost2 := [][]float64{{9, 2}, {6, 4}, {7, 3}}
	m2, _ := CostMatrixFrom(cost2)
	res2 := HungarianMinCost(m2)
	if res2.Size() != 2 {
		t.Errorf("assigned %d rows, want 2", res2.Size())
	}
	if !IsValidAssignment(m2, res2.Assignment) {
		t.Errorf("invalid assignment %v", res2.Assignment)
	}
}

func TestMinCostMaxFlow(t *testing.T) {
	build := func() *MinCostNetwork {
		g := NewMinCostNetwork(4)
		g.AddEdge(0, 1, 1, 1)
		g.AddEdge(0, 2, 1, 2)
		g.AddEdge(1, 3, 1, 1)
		g.AddEdge(2, 3, 1, 1)
		g.AddEdge(1, 2, 1, 1)
		return g
	}
	f, c := MinCostMaxFlow(build(), 0, 3)
	if f != 2 || c != 5 {
		t.Errorf("SPFA min-cost max-flow = (%d,%d), want (2,5)", f, c)
	}
	f2, c2 := MinCostMaxFlowDijkstra(build(), 0, 3)
	if f2 != 2 || c2 != 5 {
		t.Errorf("Dijkstra min-cost max-flow = (%d,%d), want (2,5)", f2, c2)
	}
	// Bounded flow of 1 unit takes the single cheapest path (cost 2).
	f3, c3 := MinCostFlow(build(), 0, 3, 1)
	if f3 != 1 || c3 != 2 {
		t.Errorf("bounded min-cost flow = (%d,%d), want (1,2)", f3, c3)
	}
}

func TestMinCostNegativeEdges(t *testing.T) {
	// A negative-cost edge should be preferred when it is cheaper overall.
	g := NewMinCostNetwork(3)
	g.AddEdge(0, 1, 1, -5)
	g.AddEdge(1, 2, 1, 2)
	g.AddEdge(0, 2, 1, 10)
	f, c := MinCostMaxFlow(g, 0, 2)
	// Max flow is 2 (both edges into 2 have unit capacity): path 0-1-2 (cost -3)
	// plus path 0-2 (cost 10) = 7.
	if f != 2 || c != 7 {
		t.Errorf("min-cost max-flow with negative edge = (%d,%d), want (2,7)", f, c)
	}
}

func TestMCMFResult(t *testing.T) {
	g := NewMinCostNetwork(4)
	g.AddEdge(0, 1, 2, 1)
	g.AddEdge(1, 3, 2, 1)
	g.AddEdge(0, 2, 2, 2)
	g.AddEdge(2, 3, 2, 1)
	res := MinCostMaxFlowResult(g, 0, 3)
	if res.FlowValue() != 4 {
		t.Errorf("flow = %d, want 4", res.FlowValue())
	}
	if res.TotalCost() != res.Residual.TotalCost() {
		t.Errorf("cost mismatch between result and residual")
	}
	if len(res.FlowEdges()) == 0 {
		t.Errorf("expected some flow edges")
	}
}

func TestGomoryHu(t *testing.T) {
	// Undirected network; the tree's pairwise min cuts must match direct max
	// flows computed independently.
	g := NewFlowNetwork(6)
	g.AddUndirectedEdge(0, 1, 3)
	g.AddUndirectedEdge(0, 2, 2)
	g.AddUndirectedEdge(1, 2, 1)
	g.AddUndirectedEdge(1, 3, 4)
	g.AddUndirectedEdge(2, 4, 2)
	g.AddUndirectedEdge(3, 4, 2)
	g.AddUndirectedEdge(3, 5, 3)
	g.AddUndirectedEdge(4, 5, 3)

	tree := GomoryHu(g)
	for u := 0; u < 6; u++ {
		for v := u + 1; v < 6; v++ {
			want := Dinic(g, u, v)
			got := tree.MinCut(u, v)
			if got != want {
				t.Errorf("min cut (%d,%d): tree=%d, direct=%d", u, v, got, want)
			}
		}
	}
	if len(tree.TreeEdges()) != 5 {
		t.Errorf("Gomory-Hu tree should have 5 edges, got %d", len(tree.TreeEdges()))
	}
	// All-pairs matrix must be symmetric and consistent with MinCut.
	apc := tree.AllPairsMinCut()
	for u := 0; u < 6; u++ {
		for v := 0; v < 6; v++ {
			if apc[u][v] != apc[v][u] {
				t.Errorf("all-pairs matrix not symmetric at (%d,%d)", u, v)
			}
		}
	}
}

func TestFlowDecomposition(t *testing.T) {
	g := clrsNetwork()
	res := DinicResult(g, 0, 5)
	paths := FlowDecomposition(res.Residual, 0, 5)
	if PathFlowTotal(paths) != res.Value {
		t.Errorf("decomposed flow total = %d, want %d", PathFlowTotal(paths), res.Value)
	}
	for _, p := range paths {
		if p.Cycle {
			continue
		}
		if p.Vertices[0] != 0 || p.Vertices[len(p.Vertices)-1] != 5 {
			t.Errorf("path does not run from source to sink: %v", p.Vertices)
		}
		if p.Flow <= 0 {
			t.Errorf("non-positive path flow %d", p.Flow)
		}
	}
}

func TestFlowNetworkFromMatrix(t *testing.T) {
	cap := [][]int64{
		{0, 10, 5, 0},
		{0, 0, 15, 0},
		{0, 0, 0, 10},
		{0, 0, 0, 0},
	}
	g, err := FlowNetworkFromMatrix(cap)
	if err != nil {
		t.Fatal(err)
	}
	if g.NumEdges() != 4 {
		t.Errorf("edges = %d, want 4", g.NumEdges())
	}
	if got := Dinic(g, 0, 3); got != 10 {
		t.Errorf("max flow = %d, want 10", got)
	}
	// Round-trip the capacity matrix.
	m := g.CapacityMatrix()
	if m[0][1] != 10 || m[1][2] != 15 || m[2][3] != 10 {
		t.Errorf("capacity matrix round-trip failed: %v", m)
	}
	if _, err := FlowNetworkFromMatrix([][]int64{{0, 1}}); err != ErrDimensionMismatch {
		t.Errorf("expected ErrDimensionMismatch for non-square matrix")
	}
}

func TestResidualAndReach(t *testing.T) {
	g := clrsNetwork()
	Dinic(g, 0, 5) // does not mutate g (works on a clone)
	if g.FlowValue(0) != 0 {
		t.Errorf("Dinic mutated input")
	}
	// Reachability along positive-capacity edges from the source.
	r := g.Reachable(0)
	if len(r) != 6 {
		t.Errorf("all vertices should be reachable from source, got %d", len(r))
	}
}

func TestValidate(t *testing.T) {
	g := NewFlowNetwork(3)
	g.AddEdge(0, 1, 4)
	g.AddEdge(1, 2, 2)
	if err := g.Validate(); err != nil {
		t.Errorf("valid network reported error: %v", err)
	}
}

func TestPanicOnBadVertex(t *testing.T) {
	defer func() {
		if r := recover(); r == nil {
			t.Errorf("expected panic on out-of-range vertex")
		}
	}()
	g := NewFlowNetwork(2)
	g.AddEdge(0, 5, 1)
}

func ExampleMaxFlow() {
	// The classic CLRS network: maximum flow from vertex 0 to vertex 5.
	g := NewFlowNetwork(6)
	g.AddEdge(0, 1, 16)
	g.AddEdge(0, 2, 13)
	g.AddEdge(1, 2, 10)
	g.AddEdge(2, 1, 4)
	g.AddEdge(1, 3, 12)
	g.AddEdge(3, 2, 9)
	g.AddEdge(2, 4, 14)
	g.AddEdge(4, 3, 7)
	g.AddEdge(3, 5, 20)
	g.AddEdge(4, 5, 4)

	fmt.Println(MaxFlow(g, 0, 5))
	// Output: 23
}

func ExampleHungarianMinCost() {
	m, _ := CostMatrixFrom([][]float64{
		{4, 1, 3},
		{2, 0, 5},
		{3, 2, 2},
	})
	res := HungarianMinCost(m)
	fmt.Printf("cost=%.0f assignment=%v\n", res.Cost, res.Assignment)
	// Output: cost=5 assignment=[1 0 2]
}
