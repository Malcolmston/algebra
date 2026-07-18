package graph

import "testing"

// clrsFlowNetwork is the classic CLRS max-flow example whose maximum flow from
// s=0 to t=5 is 23.
func clrsFlowNetwork() *Graph {
	g := NewDirected()
	g.AddWeightedEdge(0, 1, 16)
	g.AddWeightedEdge(0, 2, 13)
	g.AddWeightedEdge(1, 2, 10)
	g.AddWeightedEdge(2, 1, 4)
	g.AddWeightedEdge(1, 3, 12)
	g.AddWeightedEdge(3, 2, 9)
	g.AddWeightedEdge(2, 4, 14)
	g.AddWeightedEdge(4, 3, 7)
	g.AddWeightedEdge(3, 5, 20)
	g.AddWeightedEdge(4, 5, 4)
	return g
}

func TestEdmondsKarp(t *testing.T) {
	g := clrsFlowNetwork()
	flow, err := g.EdmondsKarp(0, 5)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(flow, 23) {
		t.Fatalf("max flow = %v; want 23", flow)
	}
}

func TestMinCutEqualsMaxFlow(t *testing.T) {
	g := clrsFlowNetwork()
	flow, _ := g.EdmondsKarp(0, 5)
	cut, side, err := g.MinCut(0, 5)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(cut, flow) {
		t.Fatalf("min cut = %v; want %v (max-flow min-cut theorem)", cut, flow)
	}
	// Source must be on the source side, sink must not.
	inSide := map[int]bool{}
	for _, v := range side {
		inSide[v] = true
	}
	if !inSide[0] || inSide[5] {
		t.Fatalf("cut partition wrong: %v", side)
	}
}

func TestEdmondsKarpSameSourceSink(t *testing.T) {
	g := clrsFlowNetwork()
	if f, _ := g.EdmondsKarp(3, 3); f != 0 {
		t.Fatalf("flow with source==sink = %v; want 0", f)
	}
}

func TestMaxBipartiteMatching(t *testing.T) {
	// Left {0,1,2}, right {3,4,5}, admits a perfect matching of size 3.
	g := New()
	g.AddEdge(0, 3)
	g.AddEdge(0, 4)
	g.AddEdge(1, 4)
	g.AddEdge(1, 5)
	g.AddEdge(2, 3)
	m, size := g.MaxBipartiteMatching([]int{0, 1, 2}, []int{3, 4, 5})
	if size != 3 {
		t.Fatalf("matching size = %d; want 3", size)
	}
	// Every matched pair must be a real edge, and right endpoints unique.
	usedRight := map[int]bool{}
	for l, r := range m {
		if !g.HasEdge(l, r) {
			t.Fatalf("matched non-edge %d-%d", l, r)
		}
		if usedRight[r] {
			t.Fatalf("right vertex %d matched twice", r)
		}
		usedRight[r] = true
	}
}

func TestMaxBipartiteMatchingLimited(t *testing.T) {
	// left 0:{3,4}, 1:{3}, 2:{4} -> maximum matching is only 2.
	g := New()
	g.AddEdge(0, 3)
	g.AddEdge(0, 4)
	g.AddEdge(1, 3)
	g.AddEdge(2, 4)
	_, size := g.MaxBipartiteMatching([]int{0, 1, 2}, []int{3, 4})
	if size != 2 {
		t.Fatalf("matching size = %d; want 2", size)
	}
}
