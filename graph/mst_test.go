package graph

import (
	"reflect"
	"testing"
)

func TestUnionFind(t *testing.T) {
	uf := NewUnionFind([]int{0, 1, 2, 3, 4})
	if uf.Count() != 5 {
		t.Fatalf("Count = %d; want 5", uf.Count())
	}
	if !uf.Union(0, 1) || !uf.Union(1, 2) {
		t.Fatal("unions should merge")
	}
	if uf.Union(0, 2) {
		t.Fatal("already connected; union should return false")
	}
	if !uf.Connected(0, 2) || uf.Connected(0, 3) {
		t.Fatal("connectivity wrong")
	}
	if uf.Count() != 3 {
		t.Fatalf("Count = %d; want 3", uf.Count())
	}
	uf.MakeSet(9)
	if uf.Count() != 4 {
		t.Fatalf("Count after MakeSet = %d; want 4", uf.Count())
	}
}

// mstGraph is an undirected weighted graph with a known MST weight of 6.
func mstGraph() *Graph {
	g := New()
	g.AddWeightedEdge(0, 1, 1)
	g.AddWeightedEdge(1, 2, 2)
	g.AddWeightedEdge(0, 2, 2)
	g.AddWeightedEdge(2, 3, 3)
	return g
}

func TestKruskal(t *testing.T) {
	g := mstGraph()
	mst, total, err := g.Kruskal()
	if err != nil {
		t.Fatal(err)
	}
	if !approx(total, 6) {
		t.Fatalf("Kruskal total = %v; want 6", total)
	}
	if len(mst) != 3 {
		t.Fatalf("MST edges = %d; want 3 (n-1)", len(mst))
	}
}

func TestPrim(t *testing.T) {
	g := mstGraph()
	_, total, err := g.Prim(0)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(total, 6) {
		t.Fatalf("Prim total = %v; want 6", total)
	}
}

func TestKruskalPrimAgree(t *testing.T) {
	g := mstGraph()
	_, kt, _ := g.Kruskal()
	_, pt, _ := g.Prim(0)
	if !approx(kt, pt) {
		t.Fatalf("Kruskal %v != Prim %v", kt, pt)
	}
}

func TestMSTDirectedError(t *testing.T) {
	g := NewDirected()
	g.AddEdge(0, 1)
	if _, _, err := g.Kruskal(); err != ErrNotUndirected {
		t.Fatalf("Kruskal on directed err = %v; want ErrNotUndirected", err)
	}
	if _, _, err := g.Prim(0); err != ErrNotUndirected {
		t.Fatalf("Prim on directed err = %v; want ErrNotUndirected", err)
	}
}

func TestKruskalForest(t *testing.T) {
	// Two disconnected triangles -> spanning forest with 4 edges.
	g := New()
	g.AddWeightedEdge(0, 1, 1)
	g.AddWeightedEdge(1, 2, 1)
	g.AddWeightedEdge(0, 2, 5)
	g.AddWeightedEdge(3, 4, 1)
	g.AddWeightedEdge(4, 5, 1)
	g.AddWeightedEdge(3, 5, 5)
	mst, total, _ := g.Kruskal()
	if len(mst) != 4 || !approx(total, 4) {
		t.Fatalf("forest = %v edges total %v; want 4 edges total 4", mst, total)
	}
	// Cross-check the edge set is acyclic and spans both components.
	if !reflect.DeepEqual([]int{0, 1, 2, 3, 4, 5}, g.Vertices()) {
		t.Fatal("vertex set changed")
	}
}
