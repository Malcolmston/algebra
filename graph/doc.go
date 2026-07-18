// Package graph implements graph data structures and classical graph
// algorithms using only the Go standard library.
//
// The central type is [Graph], an adjacency-list graph whose vertices are
// arbitrary non-negative or negative integer identifiers. A graph may be
// directed or undirected and every edge carries a float64 weight (unweighted
// edges use a weight of 1). Construct graphs with [New] or [NewDirected] and
// populate them with [Graph.AddEdge] and [Graph.AddWeightedEdge].
//
// The package provides:
//
//   - Traversal and structure: [Graph.BFS], [Graph.DFS], preorder and
//     postorder walks, [Graph.ConnectedComponents], [Graph.TopologicalSort],
//     [Graph.HasCycle], [Graph.IsBipartite] and [Graph.IsTree].
//   - Shortest paths: [Graph.Dijkstra], [Graph.BellmanFord],
//     [Graph.FloydWarshall] and [Graph.AStar].
//   - Minimum spanning trees: [Graph.Kruskal] and [Graph.Prim], backed by the
//     [UnionFind] disjoint-set structure.
//   - Strong connectivity: [Graph.TarjanSCC], [Graph.KosarajuSCC] and
//     [Graph.TransitiveClosure].
//   - Coloring and flows: [Graph.GreedyColoring], [Graph.EdmondsKarp] max-flow
//     and [Graph.MaxBipartiteMatching].
//
// All algorithms are deterministic: neighbor iteration and every returned
// ordering break ties by ascending vertex identifier, so repeated runs on the
// same input yield identical results. The package has no third-party
// dependencies.
package graph
