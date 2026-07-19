// Package voronoi implements planar computational geometry centred on the
// Delaunay triangulation and its dual, the Voronoi diagram.
//
// The package is self-contained and depends only on the Go standard library.
// All coordinates are float64 and every routine is deterministic: identical
// inputs yield identical outputs.
//
// The core algorithms are:
//
//   - The incremental Bowyer-Watson algorithm for the Delaunay triangulation
//     of a planar point set ([Triangulate]).
//   - Construction of the dual Voronoi diagram from that triangulation
//     ([Voronoi]), including bounded and unbounded cells.
//   - Andrew's monotone-chain convex hull ([ConvexHull]).
//   - Alpha shapes and alpha complexes derived from the Delaunay triangulation
//     ([AlphaShape]).
//   - Nearest-neighbour, k-nearest, closest-pair, largest-empty-circle and
//     minimum-enclosing-circle queries.
//   - Point location within a triangulation ([Triangulation.Locate]).
//
// Supporting value types include [Point] (a location in the plane), [Circle],
// [Rect], [Edge] and [Triangle]. Numerically delicate predicates
// (orientation, in-circle) are computed as signed determinants so callers can
// interpret the sign directly, and approximate comparisons accept an explicit
// epsilon.
package voronoi

// Eps is the default absolute tolerance used by the approximate predicates and
// comparisons in this package. It is chosen to sit comfortably above
// double-precision rounding noise for coordinates of moderate magnitude.
const Eps = 1e-9
