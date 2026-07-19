// Package networkflow implements network-flow and matching algorithms on
// directed and undirected graphs using only the Go standard library.
//
// The package is organized around a small number of concrete graph types and
// the classical algorithms that run on them.
//
// A [FlowNetwork] is a directed graph with non-negative integer edge
// capacities. Every edge is stored together with an anti-parallel residual
// edge, so the network doubles as its own residual graph while an algorithm
// runs. Maximum flow can be computed with three interchangeable engines:
// [EdmondsKarp] (shortest-augmenting-path Ford-Fulkerson, O(V*E^2)), [Dinic]
// (blocking flows on the level graph, O(V^2*E), and O(E*sqrt(V)) on unit
// networks), and [PushRelabel] (the FIFO preflow-push method). Each engine has
// a plain integer-valued form and a Result form that also exposes the residual
// network, the flow decomposition, and the induced minimum s-t cut.
//
// Minimum cuts come in two flavours. The minimum s-t cut is read directly off
// a maximum flow by [MinCutST] via the max-flow/min-cut theorem. The global
// minimum cut of an undirected weighted graph is computed on a [WeightedGraph]
// with the Stoer-Wagner algorithm ([StoerWagner]), which needs no flow
// computation at all. The [GomoryHuTree] built by [GomoryHu] compresses the
// minimum s-t cut between every pair of vertices into a single weighted tree
// using n-1 maximum-flow computations (Gusfield's algorithm).
//
// A [MinCostNetwork] adds a per-unit cost to every edge. Minimum-cost flows are
// found by successive shortest augmenting paths, either with a
// Bellman-Ford/SPFA shortest-path search that tolerates negative edge costs
// ([MinCostMaxFlow]) or with Johnson potentials and Dijkstra for graphs whose
// reduced costs stay non-negative ([MinCostMaxFlowDijkstra]).
//
// Matching lives on a [BipartiteGraph]. Maximum-cardinality matchings are found
// with the Hopcroft-Karp algorithm ([HopcroftKarp]) or the simpler Kuhn
// augmenting-path method ([KuhnMatching]); Konig's theorem then yields a
// minimum vertex cover and a maximum independent set for free. Maximum-weight
// (or minimum-cost) perfect assignment on a [CostMatrix] is solved in O(n^3) by
// the Hungarian / Kuhn-Munkres algorithm ([HungarianMinCost],
// [MaxWeightAssignment]).
//
// All randomness, where used, flows through a caller-supplied seed and the
// deterministic generator in math/rand; the package performs no wall-clock or
// cryptographic sampling and never touches global state.
package networkflow
