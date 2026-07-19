// Package polyhedra provides tools for three-dimensional convex geometry with a
// focus on polyhedra: the five Platonic solids, the thirteen Archimedean
// solids, general triangulated and polygonal meshes, and 3-D convex hulls.
//
// The package is organised around a small number of value types. [Vec3] is the
// core three-component vector, used interchangeably as a point or a
// displacement. [Polyhedron] is a boundary representation (a set of vertices
// together with polygonal faces given as ordered index loops) on top of which
// the package computes surface area, enclosed volume, centroid, face normals
// and areas, edge enumeration, the Euler characteristic, and the (polar/canonical)
// dual polyhedron. [PlatonicSolid] and [ArchimedeanSolid] capture the
// combinatorial and metric data of the regular and semiregular convex solids
// and can materialise them as [Polyhedron] meshes with exact vertex
// coordinates.
//
// Volumes and areas of closed meshes are obtained from the divergence theorem:
// the signed volume is one sixth of the sum over outward-oriented triangles of
// the scalar triple product of their vertices, and face normals are computed
// with Newell's method so that they are well defined even for slightly
// non-planar polygons. The convex hull of a point cloud is computed with an
// incremental algorithm ([ConvexHull]); a gift-wrapping variant
// ([GiftWrapHull]) is provided as an independent cross-check.
//
// All computation is performed in float64 using only the Go standard library.
// Routines are deterministic. Predicates that must tolerate floating-point
// noise take or use an explicit epsilon so that callers control the trade-off
// between robustness and strictness.
package polyhedra
