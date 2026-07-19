// Package tropical implements algebra over the tropical (min-plus and
// max-plus) semirings and the numerical objects built on top of them.
//
// A semiring replaces ordinary addition and multiplication with two new
// operations. In the min-plus semiring the "sum" of two numbers is their
// minimum and the "product" is their ordinary sum; its neutral element for
// addition (the tropical zero) is +Inf and its neutral element for
// multiplication (the tropical one) is 0. The max-plus semiring is the
// order-dual: addition is the maximum, the tropical zero is -Inf and the
// tropical one is still 0. Both are idempotent, commutative, complete
// dioids and the two are exchanged by negating every entry.
//
// Tropical linear algebra is unexpectedly useful. Tropical matrix
// multiplication over min-plus is exactly the shortest-path relaxation, so a
// tropical matrix power gives shortest walks of a bounded length and the
// Kleene star gives all-pairs shortest paths. Over max-plus the same objects
// describe the longest paths and the timing behaviour of discrete-event
// systems: the unique eigenvalue of an irreducible matrix is the maximum
// cycle mean, computed here with Karp's algorithm. A tropical polynomial is a
// piecewise-linear convex (min-plus) or concave (max-plus) function whose
// corners are its roots; their multiplicities are read off the Newton
// polygon. The tropical determinant coincides with the tropical permanent and
// solves the optimal assignment problem, and residuation gives the greatest
// (max-plus) or least (min-plus) solution of a tropical linear system as well
// as the least solution x = A*b of the affine fixed-point equation x = Ax + b.
//
// This package depends only on the Go standard library. Scalars are ordinary
// float64 values together with the appropriate infinity for the tropical zero.
// Routines are deterministic: identical inputs always produce identical
// outputs. Where a computation can diverge (a negative cycle for min-plus, a
// positive cycle for max-plus) the corresponding function reports it instead
// of returning a meaningless number.
package tropical
