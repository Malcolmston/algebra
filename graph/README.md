# graph — implements graph data structures and classical graph algorithms using only the Go

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/graph.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/graph)

Package graph implements graph data structures and classical graph algorithms
using only the Go standard library.

The central type is `Graph`, an adjacency-list graph whose vertices are
arbitrary non-negative or negative integer identifiers. A graph may be
directed or undirected and every edge carries a float64 weight (unweighted
edges use a weight of 1). Construct graphs with `New` or `NewDirected` and
populate them with `Graph.AddEdge` and `Graph.AddWeightedEdge`.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/graph"
```

## Usage

This is the package's own `ExampleGraph_AStar`, so it compiles and its output is
asserted on every `go test ./graph/`.

```go
g := graph.New()
	for i := 0; i < 3; i++ {
		g.AddWeightedEdge(i, i+1, 1)
	}
	h := func(v int) float64 { return math.Abs(float64(3 - v)) }

	path, weight, err := g.AStar(0, 3, h)
	fmt.Println(path, weight, err)

	path, weight, _ = g.AStar(0, 3, nil)
	fmt.Println(path, weight)
```

```
[0 1 2 3] 3 <nil>
[0 1 2 3] 3
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func NumColors(color map[int]int) int` | NumColors returns the number of distinct colors that appear in a coloring. |

### Types

| Type | What it is |
| --- | --- |
| `Edge` | Edge is a weighted connection from vertex From to vertex To. |
| `Graph` | Graph is an adjacency-list graph with integer vertices and float64 edge weights. |
| `UnionFind` | UnionFind is a disjoint-set (union-find) data structure over a fixed set of integer elements, with union by rank and path compression. |

<details>
<summary><code>Graph</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func New() *Graph` | New returns an empty undirected graph. |
| `func NewDirected() *Graph` | NewDirected returns an empty directed graph. |
| `func (g *Graph) AStar(start, goal int, h func(int) float64) (path []int, weight float64, err error)` | AStar finds a minimum-weight path from start to goal over non-negative edge weights, guided by the heuristic h, which estimates the remaining cost… |
| `func (g *Graph) AddEdge(u, v int)` | AddEdge adds an edge between u and v with weight 1, creating either vertex if necessary. |
| `func (g *Graph) AddVertex(v int)` | AddVertex adds v to the graph if it is not already present. |
| `func (g *Graph) AddWeightedEdge(u, v int, w float64)` | AddWeightedEdge adds an edge from u to v with the given weight, creating either vertex if necessary. |
| `func (g *Graph) AdjacencyMatrix() (matrix [][]float64, order []int)` | AdjacencyMatrix returns the weighted adjacency matrix of the graph together with the vertex ordering (ascending) that indexes its rows and columns. |
| `func (g *Graph) BFS(start int) ([]int, error)` | BFS returns the vertices reachable from start in breadth-first order, including start itself. |
| `func (g *Graph) BFSDistances(start int) (map[int]int, error)` | BFSDistances returns, for every vertex reachable from start, the minimum number of edges from start to that vertex (start maps to 0). |
| `func (g *Graph) BellmanFord(start int) (dist map[int]float64, prev map[int]int, err error)` | BellmanFord computes single-source shortest-path distances from start and tolerates negative edge weights. |
| `func (g *Graph) BellmanFordPath(start, goal int) (path []int, weight float64, err error)` | BellmanFordPath returns the minimum-weight path from start to goal, tolerating negative edge weights, together with its total weight. |
| `func (g *Graph) ChromaticNumberUpperBound() int` | ChromaticNumberUpperBound returns an easily computed upper bound on the chromatic number: one more than the maximum vertex degree (Brooks' bound in… |
| `func (g *Graph) Condensation() (*Graph, map[int]int)` | Condensation returns the condensation of a directed graph: a new directed graph whose vertices are the indices of the strongly connected components… |
| `func (g *Graph) ConnectedComponents() [][]int` | ConnectedComponents returns the connected components of the graph. |
| `func (g *Graph) Copy() *Graph` | Copy returns a deep copy of the graph. |
| `func (g *Graph) DFS(start int) ([]int, error)` | DFS returns the vertices reachable from start in depth-first preorder, including start itself. |
| `func (g *Graph) DFSPostorder(start int) ([]int, error)` | DFSPostorder returns the vertices reachable from start in depth-first postorder: each vertex is recorded after all its descendants. |
| `func (g *Graph) DFSPreorder(start int) ([]int, error)` | DFSPreorder returns the vertices reachable from start in depth-first preorder: each vertex is recorded when first visited. |
| `func (g *Graph) Degree(v int) (int, error)` | Degree returns the degree of v: for an undirected graph the number of incident edges (a self-loop counts twice), and for a directed graph the sum of… |
| `func (g *Graph) DegreeSequence() []int` | DegreeSequence returns the degrees of all vertices sorted in non-increasing order. |
| `func (g *Graph) Dijkstra(start int) (dist map[int]float64, prev map[int]int, err error)` | Dijkstra computes single-source shortest-path distances from start over non-negative edge weights. |
| `func (g *Graph) DijkstraPath(start, goal int) (path []int, weight float64, err error)` | DijkstraPath returns the minimum-weight path from start to goal over non-negative edge weights, together with its total weight. |
| `func (g *Graph) Directed() bool` | Directed reports whether the graph is directed. |
| `func (g *Graph) Edges() []Edge` | Edges returns all edges of the graph. |
| `func (g *Graph) EdmondsKarp(source, sink int) (float64, error)` | EdmondsKarp computes the maximum flow from source to sink in a directed graph, treating each edge weight as its capacity, using the Edmonds-Karp… |
| `func (g *Graph) FloydWarshall() (map[int]map[int]float64, error)` | FloydWarshall computes all-pairs shortest-path distances. |
| `func (g *Graph) GreedyColoring() map[int]int` | GreedyColoring assigns each vertex a color (a non-negative integer) so that no edge joins two vertices of the same color, using the greedy first-fit… |
| `func (g *Graph) HasCycle() bool` | HasCycle reports whether the graph contains a cycle. |
| `func (g *Graph) HasEdge(u, v int) bool` | HasEdge reports whether an edge from u to v exists. |
| `func (g *Graph) HasVertex(v int) bool` | HasVertex reports whether v is present in the graph. |
| `func (g *Graph) InDegree(v int) (int, error)` | InDegree returns the number of edges entering v. |
| `func (g *Graph) IsBipartite() bool` | IsBipartite reports whether the graph is bipartite (2-colorable), ignoring edge orientation. |
| `func (g *Graph) IsConnected() bool` | IsConnected reports whether the graph is connected (weakly connected for directed graphs). |
| `func (g *Graph) IsDAG() bool` | IsDAG reports whether the graph is a directed acyclic graph. |
| `func (g *Graph) IsProperColoring(color map[int]int) bool` | IsProperColoring reports whether color assigns different colors to the two endpoints of every edge, ignoring edge orientation. |
| `func (g *Graph) IsTree() bool` | IsTree reports whether an undirected graph is a tree: connected and acyclic with exactly NumVertices-1 edges. |
| `func (g *Graph) KosarajuSCC() [][]int` | KosarajuSCC returns the strongly connected components of a directed graph using Kosaraju's two-pass algorithm. |
| `func (g *Graph) Kruskal() (mst []Edge, total float64, err error)` | Kruskal computes a minimum spanning tree (or minimum spanning forest, if the graph is disconnected) of an undirected graph using Kruskal's algorithm. |
| `func (g *Graph) MaxBipartiteMatching(left, right []int) (matching map[int]int, size int)` | MaxBipartiteMatching computes a maximum-cardinality matching between the left and right vertex sets of a bipartite graph using Kuhn's augmenting-path… |
| `func (g *Graph) MinCut(source, sink int) (value float64, sourceSide []int, err error)` | MinCut computes a minimum s-t cut of a directed graph whose edge weights are capacities. |
| `func (g *Graph) Neighbors(v int) ([]int, error)` | Neighbors returns the out-neighbors of v in ascending order. |
| `func (g *Graph) NumConnectedComponents() int` | NumConnectedComponents returns the number of connected components (weakly connected for directed graphs). |
| `func (g *Graph) NumEdges() int` | NumEdges returns the number of edges (the size of the graph). |
| `func (g *Graph) NumVertices() int` | NumVertices returns the number of vertices (the order of the graph). |
| `func (g *Graph) OutDegree(v int) (int, error)` | OutDegree returns the number of edges leaving v. |
| `func (g *Graph) PathTo(prev map[int]int, start, goal int) []int` | PathTo reconstructs the path from start to goal from a predecessor map as returned by `Graph.Dijkstra`, `Graph.BellmanFord` or `Graph.AStar`. |
| `func (g *Graph) Prim(start int) (mst []Edge, total float64, err error)` | Prim computes a minimum spanning tree of the connected component containing start in an undirected graph using Prim's algorithm. |
| `func (g *Graph) ReachableFrom(start int) (map[int]bool, error)` | ReachableFrom returns the set of vertices reachable from start (including start) as a map. |
| `func (g *Graph) RemoveEdge(u, v int)` | RemoveEdge removes the edge from u to v if present. |
| `func (g *Graph) RemoveVertex(v int)` | RemoveVertex removes v and all edges incident to it. |
| `func (g *Graph) Reverse() *Graph` | Reverse returns a new directed graph with every edge reversed. |
| `func (g *Graph) ShortestPathBFS(start, goal int) ([]int, bool)` | ShortestPathBFS returns a path with the fewest edges from start to goal in an unweighted sense, and whether such a path exists. |
| `func (g *Graph) TarjanSCC() [][]int` | TarjanSCC returns the strongly connected components of a directed graph using Tarjan's algorithm. |
| `func (g *Graph) TopologicalSort() ([]int, error)` | TopologicalSort returns a topological ordering of a directed acyclic graph using Kahn's algorithm. |
| `func (g *Graph) TransitiveClosure() map[int]map[int]bool` | TransitiveClosure returns the reachability relation of the graph as a nested map: closure`u``v` is true when there is a directed path (of one or more… |
| `func (g *Graph) TwoColoring() (map[int]int, bool)` | TwoColoring attempts to 2-color the graph so that adjacent vertices receive different colors (0 or 1), ignoring edge orientation. |
| `func (g *Graph) Vertices() []int` | Vertices returns all vertices in ascending order. |
| `func (g *Graph) Weight(u, v int) (float64, bool)` | Weight returns the weight of the edge from u to v and whether it exists. |

</details>

<details>
<summary><code>UnionFind</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewUnionFind(elements []int) *UnionFind` | NewUnionFind returns a union-find structure whose given elements each start in their own singleton set. |
| `func (uf *UnionFind) Connected(x, y int) bool` | Connected reports whether x and y are in the same set. |
| `func (uf *UnionFind) Count() int` | Count returns the number of disjoint sets currently maintained. |
| `func (uf *UnionFind) Find(x int) int` | Find returns the canonical representative of the set containing x, adding x as a new singleton if it was not present. |
| `func (uf *UnionFind) MakeSet(x int)` | MakeSet adds x as a new singleton set if it is not already present. |
| `func (uf *UnionFind) Union(x, y int) bool` | Union merges the sets containing x and y and reports whether a merge occurred (that is, whether they were previously in different sets). |

</details>

### Variables

`ErrNegativeCycle`, `ErrNegativeWeight`, `ErrNotDAG`, `ErrNotUndirected`, `ErrVertexNotFound`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/graph).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
