package graph

import (
	"errors"
	"math"
	"math/rand"
	"reflect"
	"testing"
)

const tol = 1e-9

func approx(a, b float64) bool { return math.Abs(a-b) <= tol }

// dijkstraGraph builds a small directed weighted graph with a known shortest
// path 0->2->1->3 of total weight 4.
func dijkstraGraph() *Graph {
	g := NewDirected()
	g.AddWeightedEdge(0, 1, 4)
	g.AddWeightedEdge(0, 2, 1)
	g.AddWeightedEdge(2, 1, 2)
	g.AddWeightedEdge(1, 3, 1)
	g.AddWeightedEdge(2, 3, 5)
	return g
}

func TestDijkstra(t *testing.T) {
	g := dijkstraGraph()
	dist, _, err := g.Dijkstra(0)
	if err != nil {
		t.Fatal(err)
	}
	want := map[int]float64{0: 0, 1: 3, 2: 1, 3: 4}
	for v, d := range want {
		if !approx(dist[v], d) {
			t.Fatalf("dist[%d] = %v; want %v", v, dist[v], d)
		}
	}
	path, w, err := g.DijkstraPath(0, 3)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(w, 4) || !reflect.DeepEqual(path, []int{0, 2, 1, 3}) {
		t.Fatalf("DijkstraPath = %v, %v; want [0 2 1 3], 4", path, w)
	}
}

func TestDijkstraNegativeWeight(t *testing.T) {
	g := NewDirected()
	g.AddWeightedEdge(0, 1, -1)
	if _, _, err := g.Dijkstra(0); !errors.Is(err, ErrNegativeWeight) {
		t.Fatalf("err = %v; want ErrNegativeWeight", err)
	}
}

func TestBellmanFord(t *testing.T) {
	g := NewDirected()
	g.AddWeightedEdge(0, 1, 1)
	g.AddWeightedEdge(1, 2, -2)
	g.AddWeightedEdge(0, 2, 4)
	dist, _, err := g.BellmanFord(0)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(dist[2], -1) {
		t.Fatalf("dist[2] = %v; want -1", dist[2])
	}
	path, w, err := g.BellmanFordPath(0, 2)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(w, -1) || !reflect.DeepEqual(path, []int{0, 1, 2}) {
		t.Fatalf("BellmanFordPath = %v, %v; want [0 1 2], -1", path, w)
	}
}

func TestBellmanFordNegativeCycle(t *testing.T) {
	g := NewDirected()
	g.AddWeightedEdge(0, 1, 1)
	g.AddWeightedEdge(1, 0, -2)
	if _, _, err := g.BellmanFord(0); !errors.Is(err, ErrNegativeCycle) {
		t.Fatalf("err = %v; want ErrNegativeCycle", err)
	}
}

func TestDijkstraVsBellmanFord(t *testing.T) {
	// On a non-negative graph both must agree.
	g := dijkstraGraph()
	dd, _, _ := g.Dijkstra(0)
	bd, _, _ := g.BellmanFord(0)
	for v := range dd {
		if !approx(dd[v], bd[v]) {
			t.Fatalf("mismatch at %d: dijkstra %v bellman %v", v, dd[v], bd[v])
		}
	}
}

func TestFloydWarshall(t *testing.T) {
	g := dijkstraGraph()
	dist, err := g.FloydWarshall()
	if err != nil {
		t.Fatal(err)
	}
	if !approx(dist[0][3], 4) || !approx(dist[0][1], 3) || !approx(dist[2][3], 3) {
		t.Fatalf("FloydWarshall distances wrong: %v", dist)
	}
	if !math.IsInf(dist[3][0], 1) {
		t.Fatalf("dist[3][0] = %v; want +Inf", dist[3][0])
	}
	// Cross-check against Dijkstra from each source.
	for _, s := range g.Vertices() {
		dd, _, _ := g.Dijkstra(s)
		for _, v := range g.Vertices() {
			if d, ok := dd[v]; ok {
				if !approx(d, dist[s][v]) {
					t.Fatalf("floyd/dijkstra mismatch %d->%d: %v vs %v", s, v, dist[s][v], d)
				}
			} else if !math.IsInf(dist[s][v], 1) {
				t.Fatalf("floyd says reachable %d->%d but dijkstra does not", s, v)
			}
		}
	}
}

func TestAStar(t *testing.T) {
	g := dijkstraGraph()
	// A zero heuristic makes A* equivalent to Dijkstra and stays admissible.
	path, w, err := g.AStar(0, 3, func(int) float64 { return 0 })
	if err != nil {
		t.Fatal(err)
	}
	if !approx(w, 4) || !reflect.DeepEqual(path, []int{0, 2, 1, 3}) {
		t.Fatalf("AStar = %v, %v; want [0 2 1 3], 4", path, w)
	}
}

// gridGraph builds an n-by-n directed grid where each cell connects to its right
// and bottom neighbor with unit weight. It is used for benchmarking.
func gridGraph(n int) *Graph {
	g := NewDirected()
	id := func(r, c int) int { return r*n + c }
	rng := rand.New(rand.NewSource(1))
	for r := 0; r < n; r++ {
		for c := 0; c < n; c++ {
			if c+1 < n {
				g.AddWeightedEdge(id(r, c), id(r, c+1), 1+rng.Float64())
			}
			if r+1 < n {
				g.AddWeightedEdge(id(r, c), id(r+1, c), 1+rng.Float64())
			}
		}
	}
	return g
}

func BenchmarkDijkstra(b *testing.B) {
	g := gridGraph(60)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if _, _, err := g.Dijkstra(0); err != nil {
			b.Fatal(err)
		}
	}
}
