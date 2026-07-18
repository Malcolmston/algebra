package graph

import (
	"errors"
	"reflect"
	"testing"
)

func equalInts(a, b []int) bool { return reflect.DeepEqual(a, b) }

func TestBasicConstruction(t *testing.T) {
	g := New()
	g.AddWeightedEdge(0, 1, 2.5)
	g.AddEdge(1, 2)
	if g.Directed() {
		t.Fatal("New should be undirected")
	}
	if !g.HasVertex(0) || !g.HasVertex(2) {
		t.Fatal("missing vertices")
	}
	if !g.HasEdge(1, 0) {
		t.Fatal("undirected edge should be symmetric")
	}
	if w, ok := g.Weight(0, 1); !ok || w != 2.5 {
		t.Fatalf("weight = %v, %v; want 2.5, true", w, ok)
	}
	if g.NumVertices() != 3 {
		t.Fatalf("NumVertices = %d; want 3", g.NumVertices())
	}
	if g.NumEdges() != 2 {
		t.Fatalf("NumEdges = %d; want 2", g.NumEdges())
	}
	if !equalInts(g.Vertices(), []int{0, 1, 2}) {
		t.Fatalf("Vertices = %v", g.Vertices())
	}
}

func TestDirectedDegrees(t *testing.T) {
	g := NewDirected()
	g.AddEdge(0, 1)
	g.AddEdge(0, 2)
	g.AddEdge(2, 0)
	if d, _ := g.OutDegree(0); d != 2 {
		t.Fatalf("OutDegree(0) = %d; want 2", d)
	}
	if d, _ := g.InDegree(0); d != 1 {
		t.Fatalf("InDegree(0) = %d; want 1", d)
	}
	if d, _ := g.Degree(0); d != 3 {
		t.Fatalf("Degree(0) = %d; want 3", d)
	}
	if g.NumEdges() != 3 {
		t.Fatalf("NumEdges = %d; want 3", g.NumEdges())
	}
}

func TestRemoveAndErrors(t *testing.T) {
	g := New()
	g.AddEdge(0, 1)
	g.AddEdge(1, 2)
	g.RemoveEdge(0, 1)
	if g.HasEdge(0, 1) {
		t.Fatal("edge should be removed")
	}
	g.RemoveVertex(2)
	if g.HasVertex(2) {
		t.Fatal("vertex should be removed")
	}
	if g.HasEdge(1, 2) {
		t.Fatal("incident edge should be removed")
	}
	if _, err := g.Neighbors(99); !errors.Is(err, ErrVertexNotFound) {
		t.Fatalf("Neighbors(absent) err = %v; want ErrVertexNotFound", err)
	}
}

func TestReverseAndAdjacencyMatrix(t *testing.T) {
	g := NewDirected()
	g.AddWeightedEdge(0, 1, 3)
	r := g.Reverse()
	if !r.HasEdge(1, 0) || r.HasEdge(0, 1) {
		t.Fatal("reverse failed")
	}
	m, order := g.AdjacencyMatrix()
	if !equalInts(order, []int{0, 1}) {
		t.Fatalf("order = %v", order)
	}
	if m[0][1] != 3 || m[1][0] != 0 {
		t.Fatalf("adjacency matrix = %v", m)
	}
}

func TestDegreeSequence(t *testing.T) {
	// Path 0-1-2-3: degrees 1,2,2,1 -> sorted desc 2,2,1,1.
	g := New()
	g.AddEdge(0, 1)
	g.AddEdge(1, 2)
	g.AddEdge(2, 3)
	if got := g.DegreeSequence(); !equalInts(got, []int{2, 2, 1, 1}) {
		t.Fatalf("DegreeSequence = %v; want [2 2 1 1]", got)
	}
}
