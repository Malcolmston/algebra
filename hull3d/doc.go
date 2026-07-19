// Package hull3d implements three-dimensional computational geometry over
// float64 coordinates using only the Go standard library.
//
// The package centres on convex polytopes in R^3 and the algorithms that build
// and query them:
//
//   - Convex hulls of finite point clouds via an incremental / QuickHull style
//     construction ([ConvexHull], [QuickHull]) with a gift-wrapping variant
//     ([GiftWrapHull]) provided as an independent cross-check.
//   - Delaunay tetrahedralization ([Delaunay]) using the incremental
//     Bowyer–Watson algorithm and its dual, the 3-D Voronoi diagram
//     ([VoronoiCells]).
//   - Half-space intersection ([HalfSpaceIntersection]) computed as the polar
//     dual of a convex hull.
//   - Minkowski sums and differences of convex polytopes ([MinkowskiSum],
//     [MinkowskiDifference]).
//   - Collision and proximity queries between convex bodies via the
//     Gilbert–Johnson–Keerthi distance algorithm ([GJKDistance],
//     [GJKIntersect]) and Expanding Polytope Algorithm penetration depth
//     ([EPAPenetration]).
//   - Metric and combinatorial queries on polytopes: volume, surface area and
//     centroid ([Polytope.Volume], [Polytope.SurfaceArea], [Polytope.Centroid]),
//     face / edge / vertex enumeration, point containment, and supporting and
//     separating planes.
//
// The fundamental value type is [Vec3], a three-component vector used
// interchangeably as a point or a displacement, equipped with a full algebra of
// linear operations. [Plane] represents an oriented plane (the boundary of a
// half-space); [Polytope] is a boundary representation storing vertices together
// with triangular faces given as oriented index triples.
//
// Volumes, areas and centroids of closed triangulated surfaces are obtained from
// the divergence theorem: the signed volume is one sixth of the sum over
// outward-oriented triangles of the scalar triple product of their vertices.
// Orientation predicates ([Orient3D], [InSphere]) use the sign of a determinant;
// they take or use an explicit epsilon so callers control the trade-off between
// robustness and strictness.
//
// All computation is deterministic. Where randomness is useful (for symbolic
// perturbation or shuffling of insertion order) the caller supplies a
// [math/rand.Rand] seeded however they wish; the package never reads the clock
// or a cryptographic source.
package hull3d
