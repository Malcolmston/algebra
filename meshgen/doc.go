// Package meshgen implements two- and three-dimensional mesh generation and
// iso-contouring using only the Go standard library.
//
// The package is self-contained and deterministic: identical inputs (including
// any caller-supplied random seed) always produce identical outputs. All
// coordinates are float64.
//
// The main capabilities are:
//
//   - Planar Delaunay triangulation via the incremental Bowyer-Watson
//     algorithm ([Triangulate], [Triangulation]).
//   - Constrained Delaunay triangulation by edge recovery
//     ([Triangulation.InsertConstraint], [TriangulateConstrained]).
//   - Ruppert-style Delaunay refinement to a minimum angle
//     ([Triangulation.Refine], [RefineToMinAngle]).
//   - Laplacian and centroidal (Lloyd) mesh smoothing ([LaplacianSmooth],
//     [LloydSmooth]).
//   - Triangle quality metrics: angles, radius ratios, aspect ratio,
//     shape regularity and edge-length statistics ([TriangleQuality]).
//   - Marching squares for planar iso-contours ([MarchingSquares]) and
//     marching (tetrahedra) cubes for 3-D iso-surfaces ([MarchingCubes]).
//   - Regular grid and quadtree background meshing ([GridMesh],
//     [Quadtree]).
//   - Mesh connectivity: edges, boundary detection, vertex and triangle
//     adjacency, and Euler characteristics ([Mesh]).
//   - Area and quality statistics over a whole mesh ([MeshStats]).
//
// Supporting value types include [Vec2] and [Vec3] (points/vectors in the
// plane and in space), [Edge], [Tri] and [Segment]. Numerically delicate
// geometric predicates (orientation, in-circle) are exposed as signed
// determinants so callers can interpret the sign directly.
package meshgen

// Eps is the default absolute tolerance used by the approximate predicates and
// comparisons in this package.
const Eps = 1e-9
