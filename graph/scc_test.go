package graph

import (
	"reflect"
	"testing"
)

func sccGraph() *Graph {
	g := NewDirected()
	g.AddEdge(0, 1)
	g.AddEdge(1, 2)
	g.AddEdge(2, 0) // SCC {0,1,2}
	g.AddEdge(2, 3)
	g.AddEdge(3, 4)
	g.AddEdge(4, 3) // SCC {3,4}
	return g
}

func TestTarjanSCC(t *testing.T) {
	sccs := sccGraph().TarjanSCC()
	want := [][]int{{0, 1, 2}, {3, 4}}
	if !reflect.DeepEqual(sccs, want) {
		t.Fatalf("TarjanSCC = %v; want %v", sccs, want)
	}
}

func TestKosarajuSCC(t *testing.T) {
	sccs := sccGraph().KosarajuSCC()
	want := [][]int{{0, 1, 2}, {3, 4}}
	if !reflect.DeepEqual(sccs, want) {
		t.Fatalf("KosarajuSCC = %v; want %v", sccs, want)
	}
}

func TestTarjanKosarajuAgree(t *testing.T) {
	g := sccGraph()
	if !reflect.DeepEqual(g.TarjanSCC(), g.KosarajuSCC()) {
		t.Fatal("Tarjan and Kosaraju disagree")
	}
}

func TestTransitiveClosure(t *testing.T) {
	g := NewDirected()
	g.AddEdge(0, 1)
	g.AddEdge(1, 2)
	cl := g.TransitiveClosure()
	if !cl[0][1] || !cl[0][2] || !cl[1][2] {
		t.Fatalf("closure missing reachability: %v", cl)
	}
	if cl[2][0] {
		t.Fatal("2 should not reach 0")
	}
	if cl[0][0] {
		t.Fatal("0 has no cycle, should not reach itself")
	}
}

func TestTransitiveClosureCycle(t *testing.T) {
	g := NewDirected()
	g.AddEdge(0, 1)
	g.AddEdge(1, 0)
	cl := g.TransitiveClosure()
	if !cl[0][0] || !cl[1][1] {
		t.Fatal("vertices on a cycle should reach themselves")
	}
}

func TestCondensation(t *testing.T) {
	g := sccGraph()
	c, compOf := g.Condensation()
	if compOf[0] != compOf[1] || compOf[0] != compOf[2] {
		t.Fatal("0,1,2 should share a component")
	}
	if compOf[3] != compOf[4] || compOf[0] == compOf[3] {
		t.Fatal("component assignment wrong")
	}
	if !c.IsDAG() {
		t.Fatal("condensation must be a DAG")
	}
	if c.NumVertices() != 2 {
		t.Fatalf("condensation vertices = %d; want 2", c.NumVertices())
	}
}
