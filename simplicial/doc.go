// Package simplicial implements abstract simplicial complexes and the core
// machinery of computational topology over the Go standard library alone.
//
// The fundamental combinatorial object is a [Simplex]: a finite, non-empty set
// of integer vertices, stored in ascending order. A k-simplex has k+1 vertices
// and dimension k. Simplices are glued together into an abstract
// [Complex], a downward-closed family of simplices — whenever a simplex
// belongs to the complex so does every one of its faces. From a complex the
// package builds the chain complex of boundary maps
//
//	… → C_{k+1} → C_k → C_{k-1} → … → C_0 → 0
//
// in three flavours of coefficients: the two-element field GF(2) ([GF2Matrix]),
// the rationals Q ([RatMatrix]) and the integers Z ([IntMatrix]). Over a field
// the ranks of consecutive boundary maps give the Betti numbers
//
//	b_k = dim ker ∂_k − rank ∂_{k+1},
//
// and their alternating sum reproduces the Euler characteristic. Over Z the
// Smith normal form ([IntMatrix.SmithNormalForm]) additionally exposes the
// torsion of the homology groups: the elementary divisors of ∂_{k+1} greater
// than one are exactly the torsion coefficients of H_k. All of this is wrapped
// by [Complex.BettiNumbers], [Complex.HomologyZ] and friends.
//
// The package also turns geometry into topology. A [PointCloud] in Euclidean
// space induces the [VietorisRips] and [Cech] complexes at a chosen scale,
// the latter using an exact minimal-enclosing-ball test ([MinimalEnclosingBall]).
// Letting the scale grow produces a [Filtration], whose evolving homology is
// summarised by persistent homology: the standard matrix-reduction algorithm
// ([PersistentHomology]) returns [PersistencePair]s — a birth/death time for
// every topological feature — from which barcodes and persistence diagrams are
// read off.
//
// A small zoo of named complexes ([SphereComplex], [TorusComplex],
// [ProjectivePlaneComplex], [Cone], [Suspension], …) is
// provided so that the homology routines can be exercised against spaces with
// known invariants.
//
// Everything is deterministic and depends only on the standard library:
// integer and rational arithmetic use math/big, so ranks, invariant factors
// and torsion coefficients are computed exactly with no floating-point error.
package simplicial
