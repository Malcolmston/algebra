package graph

import "testing"

func TestGreedyColoring(t *testing.T) {
	// 4-cycle is bipartite; greedy first-fit should use 2 colors.
	g := New()
	g.AddEdge(0, 1)
	g.AddEdge(1, 2)
	g.AddEdge(2, 3)
	g.AddEdge(3, 0)
	color := g.GreedyColoring()
	if !g.IsProperColoring(color) {
		t.Fatalf("coloring not proper: %v", color)
	}
	if n := NumColors(color); n != 2 {
		t.Fatalf("NumColors = %d; want 2", n)
	}
}

func TestGreedyColoringTriangle(t *testing.T) {
	// Triangle (K3) needs 3 colors.
	g := New()
	g.AddEdge(0, 1)
	g.AddEdge(1, 2)
	g.AddEdge(2, 0)
	color := g.GreedyColoring()
	if !g.IsProperColoring(color) {
		t.Fatalf("coloring not proper: %v", color)
	}
	if n := NumColors(color); n != 3 {
		t.Fatalf("NumColors = %d; want 3", n)
	}
	if ub := g.ChromaticNumberUpperBound(); ub < 3 {
		t.Fatalf("upper bound = %d; must be >= chromatic number 3", ub)
	}
}

func TestImproperColoringDetected(t *testing.T) {
	g := New()
	g.AddEdge(0, 1)
	bad := map[int]int{0: 0, 1: 0}
	if g.IsProperColoring(bad) {
		t.Fatal("equal-colored adjacent vertices must be improper")
	}
}
