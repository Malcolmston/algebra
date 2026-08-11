# networkflow — implements network-flow and matching algorithms on directed and undirected graphs using

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/networkflow.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/networkflow)

Package networkflow implements network-flow and matching algorithms on
directed and undirected graphs using only the Go standard library.

The package is organized around a small number of concrete graph types and the
classical algorithms that run on them.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/networkflow"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func AssignmentCost(m *CostMatrix, assignment []int) float64` | AssignmentCost returns the total cost of the given row-to-column assignment under matrix m. |
| `func Dinic(g *FlowNetwork, s, t int) int64` | Dinic returns the value of the maximum flow from s to t using Dinic's algorithm. |
| `func EdgeConnectivity(g *WeightedGraph) float64` | EdgeConnectivity returns the weighted edge connectivity of an undirected graph, i.e. |
| `func EdmondsKarp(g *FlowNetwork, s, t int) int64` | EdmondsKarp returns the value of the maximum flow from s to t using the Edmonds-Karp algorithm. |
| `func GlobalMinCut(g *WeightedGraph) float64` | GlobalMinCut is an alias for `StoerWagner` returning only the cut value. |
| `func HasPerfectMatching(g *BipartiteGraph) bool` | HasPerfectMatching reports whether the bipartite graph has a matching that saturates the smaller side. |
| `func IsValidAssignment(m *CostMatrix, assignment []int) bool` | IsValidAssignment reports whether assignment maps rows to distinct, in-range columns (allowing -1 for unassigned rows). |
| `func MaxFlow(g *FlowNetwork, s, t int) int64` | MaxFlow returns the value of the maximum flow from s to t using the default engine (Dinic). |
| `func MaximumIndependentSet(g *BipartiteGraph) (left, right []int)` | MaximumIndependentSet returns a maximum independent set of a bipartite graph as (left, right) vertex lists. |
| `func MaximumMatchingSize(g *BipartiteGraph) int` | MaximumMatchingSize returns just the size of a maximum-cardinality matching. |
| `func MinCostFlow(g *MinCostNetwork, s, t int, want int64) (flow, cost int64)` | MinCostFlow ships up to want units from s to t at least total cost, returning the flow actually shipped and its cost. |
| `func MinCostMaxFlow(g *MinCostNetwork, s, t int) (flow, cost int64)` | MinCostMaxFlow returns the maximum flow from s to t of least total cost, found by successive shortest paths with SPFA (Bellman-Ford). |
| `func MinCostMaxFlowDijkstra(g *MinCostNetwork, s, t int) (flow, cost int64)` | MinCostMaxFlowDijkstra returns the minimum-cost maximum flow using Johnson potentials with Dijkstra. |
| `func MinCutSTValue(g *FlowNetwork, s, t int) int64` | MinCutSTValue returns just the value of the minimum s-t cut, equal to the maximum flow. |
| `func MinimumVertexCover(g *BipartiteGraph) (left, right []int)` | MinimumVertexCover returns a minimum vertex cover of a bipartite graph as two sets (left cover vertices, right cover vertices) via Konig's theorem… |
| `func PathFlowTotal(paths []FlowPath) int64` | PathFlowTotal sums the flow over the given decomposition paths that are not cycles; for a valid decomposition this equals the flow value. |
| `func PushRelabel(g *FlowNetwork, s, t int) int64` | PushRelabel returns the value of the maximum flow from s to t using the FIFO preflow-push method. |
| `func PushRelabelHighestLabel(g *FlowNetwork, s, t int) int64` | PushRelabelHighestLabel returns the maximum flow from s to t using preflow push with the highest-label selection rule. |

### Types

| Type | What it is |
| --- | --- |
| `AssignmentResult` | AssignmentResult describes an optimal assignment. |
| `BipartiteGraph` | BipartiteGraph is an undirected bipartite graph with a left vertex set {0..L-1} and a right vertex set {0..R-1}. |
| `CostMatrix` | CostMatrix is a rectangular matrix of assignment costs with rows (workers) and columns (jobs). |
| `Edge` | Edge is a read-only view of a directed edge in a `FlowNetwork`. |
| `FlowNetwork` | FlowNetwork is a directed graph with non-negative integer edge capacities. |
| `FlowPath` | FlowPath is a single source-to-sink path (or cycle) carrying a fixed amount of flow, produced by `FlowDecomposition`. |
| `GomoryHuTree` | GomoryHuTree is a weighted tree on the same vertex set as a source `FlowNetwork` such that, for every pair of vertices u and v, the minimum weight of… |
| `MCMFResult` | MCMFResult carries the outcome of a minimum-cost flow computation. |
| `MatchingResult` | MatchingResult describes a matching in a bipartite graph. |
| `MaxFlowResult` | MaxFlowResult carries the outcome of a maximum-flow computation together with the residual network, from which the minimum s-t cut and a flow… |
| `MinCostNetwork` | MinCostNetwork is a directed graph with non-negative integer capacities and integer per-unit costs on each edge. |
| `MinCutResult` | MinCutResult describes a cut: a bipartition of the vertices and the total capacity of the edges crossing from one side to the other. |
| `WeightedGraph` | WeightedGraph is an undirected graph on n vertices with real, non-negative edge weights stored as a symmetric adjacency matrix. |

<details>
<summary><code>AssignmentResult</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func HungarianMinCost(m *CostMatrix) *AssignmentResult` | HungarianMinCost solves the minimum-cost assignment problem for the given cost matrix using the Hungarian (Kuhn-Munkres) algorithm in O(n^3) time. |
| `func KuhnMunkres(m *CostMatrix) *AssignmentResult` | KuhnMunkres is an alias for `HungarianMinCost`. |
| `func MaxWeightAssignment(m *CostMatrix) *AssignmentResult` | MaxWeightAssignment solves the maximum-weight assignment problem: it assigns rows to distinct columns to maximize total weight, by negating the… |
| `func (a *AssignmentResult) Pairs() [][2]int` | Pairs returns the assigned (row, column) pairs. |
| `func (a *AssignmentResult) Size() int` | Size returns the number of assigned rows. |

</details>

<details>
<summary><code>BipartiteGraph</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func BipartiteFromMatrix(adj [][]bool) (*BipartiteGraph, error)` | BipartiteFromMatrix builds a `BipartiteGraph` from a boolean adjacency matrix where adj`i``j` true means left i is joined to right j. |
| `func NewBipartiteGraph(nl, nr int) *BipartiteGraph` | NewBipartiteGraph returns an empty bipartite graph with nl left vertices and nr right vertices. |
| `func (g *BipartiteGraph) AddEdge(u, v int)` | AddEdge adds an edge between left vertex u and right vertex v. |
| `func (g *BipartiteGraph) Clone() *BipartiteGraph` | Clone returns a deep copy of the bipartite graph. |
| `func (g *BipartiteGraph) Edges() [][2]int` | Edges returns every edge of the bipartite graph as {left, right} pairs, sorted lexicographically. |
| `func (g *BipartiteGraph) HasEdge(u, v int) bool` | HasEdge reports whether an edge joins left vertex u and right vertex v. |
| `func (g *BipartiteGraph) LeftNeighbors(u int) []int` | LeftNeighbors returns the right vertices adjacent to left vertex u, sorted. |
| `func (g *BipartiteGraph) LeftSize() int` | LeftSize returns the number of left vertices. |
| `func (g *BipartiteGraph) NumEdges() int` | NumEdges returns the total number of edges. |
| `func (g *BipartiteGraph) RightNeighbors(v int) []int` | RightNeighbors returns the left vertices adjacent to right vertex v, sorted. |
| `func (g *BipartiteGraph) RightSize() int` | RightSize returns the number of right vertices. |

</details>

<details>
<summary><code>CostMatrix</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func CostMatrixFrom(m [][]float64) (*CostMatrix, error)` | CostMatrixFrom wraps a rectangular slice as a `CostMatrix` (copying it). |
| `func NewCostMatrix(rows, cols int) *CostMatrix` | NewCostMatrix returns a rows-by-cols cost matrix initialized to zero. |
| `func (m *CostMatrix) At(i, j int) float64` | At returns the cost of assigning row i to column j. |
| `func (m *CostMatrix) Clone() *CostMatrix` | Clone returns a deep copy of the cost matrix. |
| `func (m *CostMatrix) Cols() int` | Cols returns the number of columns (jobs). |
| `func (m *CostMatrix) Negate() *CostMatrix` | Negate returns a new cost matrix with every entry negated; solving a minimum-cost assignment on it maximizes the original weights. |
| `func (m *CostMatrix) Rows() int` | Rows returns the number of rows (workers). |
| `func (m *CostMatrix) Set(i, j int, cost float64)` | Set sets the cost of assigning row i to column j. |
| `func (m *CostMatrix) Transpose() *CostMatrix` | Transpose returns a new cost matrix with rows and columns exchanged. |

</details>

<details>
<summary><code>Edge</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func MinCutSTEdges(g *FlowNetwork, s, t int) []Edge` | MinCutSTEdges returns the directed caller edges crossing the minimum s-t cut, whose capacities sum to the cut value. |
| `func (e Edge) IsSaturated() bool` | IsSaturated reports whether the edge carries flow equal to its capacity. |
| `func (e Edge) Residual() int64` | Residual returns the remaining capacity Cap-Flow of the edge. |
| `func (e Edge) String() string` | String renders the edge as "from->to cap/flow". |

</details>

<details>
<summary><code>FlowNetwork</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func FlowNetworkFromMatrix(cap [][]int64) (*FlowNetwork, error)` | FlowNetworkFromMatrix builds a directed `FlowNetwork` from a square capacity matrix: an edge i->j is added for every strictly positive cap`i``j`. |
| `func NewFlowNetwork(n int) *FlowNetwork` | NewFlowNetwork returns an empty flow network on n vertices (labelled 0..n-1). |
| `func (g *FlowNetwork) AddEdge(from, to int, capacity int64) int` | AddEdge adds a directed edge from->to with the given capacity and returns its edge id. |
| `func (g *FlowNetwork) AddUndirectedEdge(u, v int, capacity int64) int` | AddUndirectedEdge adds an undirected edge {u,v} with the given capacity in both directions and returns its edge id. |
| `func (g *FlowNetwork) Capacity(u, v int) int64` | Capacity returns the total capacity of every arc from u to v. |
| `func (g *FlowNetwork) CapacityMatrix() [][]int64` | CapacityMatrix returns the n-by-n matrix whose (i,j) entry is the total capacity of caller edges from i to j. |
| `func (g *FlowNetwork) Clone() *FlowNetwork` | Clone returns a deep copy of the network, including current flow values. |
| `func (g *FlowNetwork) Edge(id int) Edge` | Edge returns a read-only view of the directed edge with the given id (as returned by `FlowNetwork.AddEdge`). |
| `func (g *FlowNetwork) Edges() []Edge` | Edges returns a read-only view of every caller-added directed edge, in the order the edges were added. |
| `func (g *FlowNetwork) FlowOn(u, v int) int64` | FlowOn returns the net flow of all caller edges from u to v. |
| `func (g *FlowNetwork) FlowValue(s int) int64` | FlowValue returns the value of the current flow measured as the net flow out of source s (out-flow minus in-flow). |
| `func (g *FlowNetwork) HasEdge(u, v int) bool` | HasEdge reports whether at least one caller edge runs from u to v. |
| `func (g *FlowNetwork) InDegree(v int) int` | InDegree returns the number of caller edges entering v. |
| `func (g *FlowNetwork) IsFeasibleFlow(s, t int) bool` | IsFeasibleFlow reports whether the current flow respects capacities on every edge and conserves flow at every vertex other than s and t. |
| `func (g *FlowNetwork) Neighbors(v int) []int` | Neighbors returns the distinct head vertices reachable from v along positive-capacity edges, sorted ascending. |
| `func (g *FlowNetwork) NumEdges() int` | NumEdges returns the number of directed edges the caller added; residual twin edges are not counted. |
| `func (g *FlowNetwork) NumVertices() int` | NumVertices returns the number of vertices in the network. |
| `func (g *FlowNetwork) OutDegree(v int) int` | OutDegree returns the number of caller edges leaving v. |
| `func (g *FlowNetwork) OutEdges(v int) []int` | OutEdges returns the ids of every residual arc leaving v (both caller edges and reverse twins). |
| `func (g *FlowNetwork) Reachable(src int) []int` | Reachable returns the vertices reachable from src along positive-capacity caller edges (ignoring flow), sorted ascending. |
| `func (g *FlowNetwork) ReachableInResidual(src int) []bool` | ReachableInResidual returns a boolean slice marking every vertex reachable from src along arcs with positive residual capacity (a BFS of the residual… |
| `func (g *FlowNetwork) ResetFlow()` | ResetFlow zeroes the flow on every edge, leaving capacities untouched. |
| `func (g *FlowNetwork) ResidualCapacity(id int) int64` | ResidualCapacity returns the residual capacity along the edge with the given id, i.e. |
| `func (g *FlowNetwork) ResidualEdges() []Edge` | ResidualEdges returns, for every residual arc with positive residual capacity, an `Edge` whose Cap field holds the residual capacity and whose Flow… |
| `func (g *FlowNetwork) TotalCapacity() int64` | TotalCapacity returns the sum of all caller-edge capacities. |
| `func (g *FlowNetwork) Validate() error` | Validate reports the first structural problem it finds, or nil if the network is well formed (all endpoints in range and all capacities non-negative). |

</details>

<details>
<summary><code>FlowPath</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func FlowDecomposition(g *FlowNetwork, s, t int) []FlowPath` | FlowDecomposition decomposes the flow currently installed in g into source-to-sink paths and, if any circulation remains, cycles. |
| `func (p FlowPath) Length() int` | Length returns the number of edges in the path. |

</details>

<details>
<summary><code>GomoryHuTree</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func GomoryHu(g *FlowNetwork) *GomoryHuTree` | GomoryHu builds a Gomory-Hu cut tree of an undirected flow network using Gusfield's algorithm, which performs n-1 maximum-flow computations. |
| `func (t *GomoryHuTree) AllPairsMinCut() [][]int64` | AllPairsMinCut returns the full symmetric matrix of pairwise minimum cut values derived from the tree. |
| `func (t *GomoryHuTree) MinCut(u, v int) int64` | MinCut returns the minimum u-v cut value, computed as the minimum edge weight on the tree path between u and v. |
| `func (t *GomoryHuTree) NumVertices() int` | NumVertices returns the number of vertices in the tree. |
| `func (t *GomoryHuTree) Parent(v int) int` | Parent returns the tree parent of v. |
| `func (t *GomoryHuTree) ParentWeight(v int) int64` | ParentWeight returns the weight of the tree edge from v to its parent, i.e. |
| `func (t *GomoryHuTree) Root() int` | Root returns the root vertex of the Gomory-Hu tree, which is always 0. |
| `func (t *GomoryHuTree) TreeEdges() [][3]int64` | TreeEdges returns the n-1 tree edges as {u, v, weight} triples with u the child and v its parent. |

</details>

<details>
<summary><code>MCMFResult</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func MinCostMaxFlowResult(g *MinCostNetwork, s, t int) *MCMFResult` | MinCostMaxFlowResult runs SPFA-based min-cost max-flow and returns the full `MCMFResult`, including the residual network. |
| `func (r *MCMFResult) FlowEdges() []Edge` | FlowEdges returns every caller edge carrying strictly positive flow. |
| `func (r *MCMFResult) FlowValue() int64` | FlowValue returns the value of the flow. |
| `func (r *MCMFResult) TotalCost() int64` | TotalCost returns the total cost of the flow. |

</details>

<details>
<summary><code>MatchingResult</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func HopcroftKarp(g *BipartiteGraph) *MatchingResult` | HopcroftKarp returns a maximum-cardinality matching of a bipartite graph computed with the Hopcroft-Karp algorithm in O(E*sqrt(V)) time. |
| `func KuhnMatching(g *BipartiteGraph) *MatchingResult` | KuhnMatching returns a maximum-cardinality matching computed with the simpler Kuhn augmenting-path method in O(V*E) time. |
| `func MaximumMatching(g *BipartiteGraph) *MatchingResult` | MaximumMatching returns a maximum-cardinality matching using the default engine (Hopcroft-Karp). |
| `func (m *MatchingResult) IsLeftMatched(u int) bool` | IsLeftMatched reports whether left vertex u is matched. |
| `func (m *MatchingResult) IsPerfect() bool` | IsPerfect reports whether every vertex on the smaller side is matched, i.e. |
| `func (m *MatchingResult) IsRightMatched(v int) bool` | IsRightMatched reports whether right vertex v is matched. |
| `func (m *MatchingResult) Pairs() [][2]int` | Pairs returns the matched (left, right) pairs sorted by left vertex. |
| `func (m *MatchingResult) UnmatchedLeft() []int` | UnmatchedLeft returns the left vertices left unmatched. |
| `func (m *MatchingResult) UnmatchedRight() []int` | UnmatchedRight returns the right vertices left unmatched. |

</details>

<details>
<summary><code>MaxFlowResult</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func DinicResult(g *FlowNetwork, s, t int) *MaxFlowResult` | DinicResult runs Dinic's algorithm and returns a `MaxFlowResult` holding the flow value and the saturated residual network. |
| `func EdmondsKarpResult(g *FlowNetwork, s, t int) *MaxFlowResult` | EdmondsKarpResult runs Edmonds-Karp and returns a `MaxFlowResult` holding the flow value and the saturated residual network. |
| `func MaxFlowResultOf(g *FlowNetwork, s, t int) *MaxFlowResult` | MaxFlowResultOf runs the default engine (Dinic) and returns the full `MaxFlowResult`. |
| `func PushRelabelResult(g *FlowNetwork, s, t int) *MaxFlowResult` | PushRelabelResult runs FIFO preflow push and returns a `MaxFlowResult` holding the flow value and the saturated residual network. |
| `func (r *MaxFlowResult) FlowEdges() []Edge` | FlowEdges returns every caller edge that carries strictly positive flow. |
| `func (r *MaxFlowResult) FlowValue() int64` | FlowValue returns the value of the maximum flow. |
| `func (r *MaxFlowResult) MinCutEdges() []Edge` | MinCutEdges returns the caller edges that cross the induced minimum cut, i.e. |
| `func (r *MaxFlowResult) MinCutValue() int64` | MinCutValue returns the capacity of the induced minimum cut, which equals the maximum-flow value. |
| `func (r *MaxFlowResult) SinkSide() []int` | SinkSide returns the complement of `MaxFlowResult.SourceSide`: the vertices not reachable from the source in the residual network. |
| `func (r *MaxFlowResult) SourceSide() []int` | SourceSide returns the set S of vertices reachable from the source in the residual network; it is the source side of the induced minimum cut. |

</details>

<details>
<summary><code>MinCostNetwork</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewMinCostNetwork(n int) *MinCostNetwork` | NewMinCostNetwork returns an empty min-cost network on n vertices. |
| `func (g *MinCostNetwork) AddEdge(from, to int, capacity, cost int64) int` | AddEdge adds a directed edge from->to with the given capacity and per-unit cost, returning its edge id. |
| `func (g *MinCostNetwork) Clone() *MinCostNetwork` | Clone returns a deep copy of the network, including current flow values. |
| `func (g *MinCostNetwork) Edge(id int) Edge` | Edge returns a read-only `Edge` view (capacity and flow) of the caller edge with the given id. |
| `func (g *MinCostNetwork) EdgeCost(id int) int64` | EdgeCost returns the per-unit cost of the edge with the given id. |
| `func (g *MinCostNetwork) Edges() []Edge` | Edges returns read-only `Edge` views of every caller edge in insertion order. |
| `func (g *MinCostNetwork) InDegree(v int) int` | InDegree returns the number of caller edges entering v. |
| `func (g *MinCostNetwork) NumEdges() int` | NumEdges returns the number of caller edges (residual twins excluded). |
| `func (g *MinCostNetwork) NumVertices() int` | NumVertices returns the number of vertices. |
| `func (g *MinCostNetwork) OutDegree(v int) int` | OutDegree returns the number of caller edges leaving v. |
| `func (g *MinCostNetwork) ResetFlow()` | ResetFlow zeroes the flow on every edge. |
| `func (g *MinCostNetwork) TotalCapacity() int64` | TotalCapacity returns the sum of all caller-edge capacities. |
| `func (g *MinCostNetwork) TotalCost() int64` | TotalCost returns the total cost of the current flow, summed as flow*cost over caller edges. |

</details>

<details>
<summary><code>MinCutResult</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func MinCut(g *FlowNetwork, s, t int) *MinCutResult` | MinCut is an alias for `MinCutST`. |
| `func MinCutST(g *FlowNetwork, s, t int) *MinCutResult` | MinCutST computes a minimum s-t cut of a directed `FlowNetwork` via the max-flow/min-cut theorem: the cut value equals the maximum flow, and the S… |
| `func StoerWagner(g *WeightedGraph) *MinCutResult` | StoerWagner computes a global minimum cut of an undirected weighted graph with the Stoer-Wagner algorithm in O(V^3) time. |
| `func (c *MinCutResult) Contains(v int) bool` | Contains reports whether vertex v lies on the S side of the cut. |
| `func (c *MinCutResult) SVertices() []int` | SVertices returns the S side of the cut. |
| `func (c *MinCutResult) SizeS() int` | SizeS returns the number of vertices on the S side. |
| `func (c *MinCutResult) SizeT() int` | SizeT returns the number of vertices on the T side. |
| `func (c *MinCutResult) TVertices() []int` | TVertices returns the T side of the cut. |

</details>

<details>
<summary><code>WeightedGraph</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewWeightedGraph(n int) *WeightedGraph` | NewWeightedGraph returns an empty undirected weighted graph on n vertices. |
| `func WeightedGraphFromMatrix(w [][]float64) (*WeightedGraph, error)` | WeightedGraphFromMatrix builds an undirected `WeightedGraph` from a symmetric weight matrix (only the upper triangle is read). |
| `func (g *WeightedGraph) AddEdge(u, v int, weight float64)` | AddEdge adds weight to the undirected edge {u,v} (accumulating on parallel additions). |
| `func (g *WeightedGraph) Clone() *WeightedGraph` | Clone returns a deep copy of the weighted graph. |
| `func (g *WeightedGraph) Edges() [][3]float64` | Edges returns every undirected edge with non-zero weight as {u, v, weight}, with u < v. |
| `func (g *WeightedGraph) NumVertices() int` | NumVertices returns the number of vertices. |
| `func (g *WeightedGraph) SetWeight(u, v int, weight float64)` | SetWeight sets the weight of the undirected edge {u,v} to weight, replacing any previous value. |
| `func (g *WeightedGraph) TotalWeight() float64` | TotalWeight returns the sum of all edge weights (each undirected edge counted once). |
| `func (g *WeightedGraph) Weight(u, v int) float64` | Weight returns the weight of the undirected edge {u,v} (0 if absent). |
| `func (g *WeightedGraph) WeightedDegree(v int) float64` | WeightedDegree returns the sum of the weights of all edges incident to v. |

</details>

### Variables

`ErrInvalidVertex`, `ErrNegativeCapacity`, `ErrSourceEqualsSink`, `ErrDimensionMismatch`, `ErrNoPerfectMatching`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/networkflow).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
