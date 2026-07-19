// Package matroids implements matroid theory and matroid-optimization
// algorithms using only the Go standard library.
//
// A matroid is a combinatorial structure that abstracts the notion of linear
// independence. Every matroid in this package has a finite ground set whose
// elements are the integers 0, 1, ..., n-1, where n is reported by
// [Matroid.Size]. The single defining operation is the rank function
// [Matroid.Rank], which maps a subset of the ground set to a non-negative
// integer satisfying the matroid rank axioms (normalisation, monotonicity,
// unit increase and submodularity). All higher-level notions — independence,
// bases, circuits, closure, flats, loops, coloops and connectivity — are
// derived from the rank oracle by the generic free functions in this package,
// so they work uniformly across every matroid representation.
//
// The package provides several concrete constructions:
//
//   - [UniformMatroid] (and the special cases [FreeMatroid] and
//     [TrivialMatroid]) in which a set is independent exactly when it is small
//     enough.
//   - [PartitionMatroid], whose independence is a per-block capacity
//     constraint; this doubles as the "transversal / partition" building block
//     used by matroid union.
//   - [GraphicMatroid], the cycle matroid of a multigraph whose independent
//     sets are the forests, backed by a union-find structure.
//   - [LinearMatroid], the column matroid of an exact rational matrix (built on
//     math/big), and [BinaryMatroid], the column matroid of a 0/1 matrix over
//     the two-element field GF(2).
//   - [TransversalMatroid], whose independent sets are the partial transversals
//     of a family of sets, computed by bipartite matching.
//   - [ExplicitMatroid], a matroid given directly by its independent sets,
//     bases or circuits, with axiom validation.
//
// Generic operations available on any [Matroid] include rank, closure and
// [Flats], [Bases], [Circuits], [Loops], [Coloops] and [Girth]; matroid
// [Duality]; the minor operations [Deletion], [Contraction] and [Minor];
// [DirectSum]; [Components] and [Connectivity]; the [Greedy] algorithm for
// maximum-weight independent sets; [Union] of several matroids; and both
// unweighted [Intersection] and [WeightedIntersection] of two matroids via
// augmenting paths in the exchange graph.
//
// No function in this package performs I/O, uses cgo, or depends on wall-clock
// time. Any randomness is driven by a caller-supplied *math/rand.Rand.
package matroids
