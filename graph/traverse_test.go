package graph

import (
	"reflect"
	"testing"
)

func TestBFSDFSOrder(t *testing.T) {
	g := New()
	g.AddEdge(0, 1)
	g.AddEdge(0, 2)
	g.AddEdge(1, 3)
	g.AddEdge(2, 3)
	bfs, _ := g.BFS(0)
	if !reflect.DeepEqual(bfs, []int{0, 1, 2, 3}) {
		t.Fatalf("BFS = %v; want [0 1 2 3]", bfs)
	}
	dfs, _ := g.DFSPreorder(0)
	if !reflect.DeepEqual(dfs, []int{0, 1, 3, 2}) {
		t.Fatalf("DFSPreorder = %v; want [0 1 3 2]", dfs)
	}
	post, _ := g.DFSPostorder(0)
	if !reflect.DeepEqual(post, []int{2, 3, 1, 0}) {
		t.Fatalf("DFSPostorder = %v; want [2 3 1 0]", post)
	}
}

func TestBFSDistances(t *testing.T) {
	g := New()
	g.AddEdge(0, 1)
	g.AddEdge(1, 2)
	g.AddEdge(2, 3)
	d, _ := g.BFSDistances(0)
	want := map[int]int{0: 0, 1: 1, 2: 2, 3: 3}
	if !reflect.DeepEqual(d, want) {
		t.Fatalf("BFSDistances = %v; want %v", d, want)
	}
}

func TestShortestPathBFS(t *testing.T) {
	g := New()
	g.AddEdge(0, 1)
	g.AddEdge(1, 2)
	g.AddEdge(0, 2)
	g.AddEdge(2, 3)
	path, ok := g.ShortestPathBFS(0, 3)
	if !ok || !reflect.DeepEqual(path, []int{0, 2, 3}) {
		t.Fatalf("ShortestPathBFS = %v, %v; want [0 2 3], true", path, ok)
	}
	if _, ok := g.ShortestPathBFS(0, 99); ok {
		t.Fatal("path to absent vertex should not exist")
	}
}

func TestConnectedComponents(t *testing.T) {
	g := New()
	g.AddEdge(0, 1)
	g.AddEdge(2, 3)
	g.AddVertex(4)
	comps := g.ConnectedComponents()
	want := [][]int{{0, 1}, {2, 3}, {4}}
	if !reflect.DeepEqual(comps, want) {
		t.Fatalf("ConnectedComponents = %v; want %v", comps, want)
	}
	if g.NumConnectedComponents() != 3 || g.IsConnected() {
		t.Fatal("connectivity mismatch")
	}
}

func TestTopologicalSort(t *testing.T) {
	g := NewDirected()
	g.AddEdge(0, 1)
	g.AddEdge(0, 2)
	g.AddEdge(1, 3)
	g.AddEdge(2, 3)
	order, err := g.TopologicalSort()
	if err != nil {
		t.Fatal(err)
	}
	if !reflect.DeepEqual(order, []int{0, 1, 2, 3}) {
		t.Fatalf("TopologicalSort = %v; want [0 1 2 3]", order)
	}
	// Verify it is a valid topological order in general.
	pos := map[int]int{}
	for i, v := range order {
		pos[v] = i
	}
	for _, e := range g.Edges() {
		if pos[e.From] >= pos[e.To] {
			t.Fatalf("edge %d->%d violates order", e.From, e.To)
		}
	}
}

func TestTopologicalSortCycle(t *testing.T) {
	g := NewDirected()
	g.AddEdge(0, 1)
	g.AddEdge(1, 2)
	g.AddEdge(2, 0)
	if _, err := g.TopologicalSort(); err == nil {
		t.Fatal("expected error on cyclic graph")
	}
}

func TestHasCycle(t *testing.T) {
	dag := NewDirected()
	dag.AddEdge(0, 1)
	dag.AddEdge(1, 2)
	if dag.HasCycle() || !dag.IsDAG() {
		t.Fatal("DAG reported cyclic")
	}
	cyc := NewDirected()
	cyc.AddEdge(0, 1)
	cyc.AddEdge(1, 0)
	if !cyc.HasCycle() || cyc.IsDAG() {
		t.Fatal("cycle not detected")
	}
	tree := New()
	tree.AddEdge(0, 1)
	tree.AddEdge(1, 2)
	if tree.HasCycle() {
		t.Fatal("tree should be acyclic")
	}
	loop := New()
	loop.AddEdge(0, 1)
	loop.AddEdge(1, 2)
	loop.AddEdge(2, 0)
	if !loop.HasCycle() {
		t.Fatal("undirected cycle not detected")
	}
}

func TestIsTree(t *testing.T) {
	tree := New()
	tree.AddEdge(0, 1)
	tree.AddEdge(0, 2)
	tree.AddEdge(2, 3)
	if !tree.IsTree() {
		t.Fatal("should be a tree")
	}
	tree.AddEdge(1, 3) // add a cycle
	if tree.IsTree() {
		t.Fatal("cyclic graph is not a tree")
	}
}

func TestBipartite(t *testing.T) {
	even := New() // 4-cycle is bipartite
	even.AddEdge(0, 1)
	even.AddEdge(1, 2)
	even.AddEdge(2, 3)
	even.AddEdge(3, 0)
	if !even.IsBipartite() {
		t.Fatal("even cycle should be bipartite")
	}
	color, ok := even.TwoColoring()
	if !ok || color[0] == color[1] {
		t.Fatal("two-coloring failed")
	}
	odd := New() // triangle is not bipartite
	odd.AddEdge(0, 1)
	odd.AddEdge(1, 2)
	odd.AddEdge(2, 0)
	if odd.IsBipartite() {
		t.Fatal("odd cycle should not be bipartite")
	}
}
