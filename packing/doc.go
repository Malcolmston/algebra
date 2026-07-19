// Package packing implements classical packing and covering algorithms.
//
// The package collects several loosely related families of problems that all
// ask how densely a set of objects can be arranged, or how few resources are
// needed to cover a demand:
//
//   - One dimensional bin packing. The online heuristics next-fit, first-fit,
//     best-fit and worst-fit together with their offline "decreasing" variants
//     (items sorted by non-increasing size before packing), plus the classical
//     Martello-Toth lower bounds L1 and L2 and the known worst-case
//     approximation ratios.
//
//   - Sphere packing in Euclidean space. Center densities, packing densities,
//     kissing numbers, covering radii and covering thickness for the standard
//     families of lattices: the integer lattice Z^n, the root lattices A_n,
//     D_n, E6, E7, E8, their duals where relevant, and the Leech lattice in
//     dimension 24.
//
//   - Circle packing. Best known arrangements of equal circles inside a unit
//     square and inside a circle, the hexagonal (densest planar) packing and
//     its covering counterpart, and the associated densities.
//
//   - Greedy approximation. The greedy 0/1 and fractional knapsack heuristics
//     and the greedy set-cover approximation together with its logarithmic
//     performance guarantee.
//
// Everything is implemented with the Go standard library only. Numeric
// quantities are float64 unless the exact value is integral. Functions that
// consult tables of best known values document the source of those values and
// expose a companion predicate reporting whether the requested instance is
// covered by the table.
//
// Geometric conventions for lattices follow Conway and Sloane, "Sphere
// Packings, Lattices and Groups" (SPLAG). Root lattices are scaled so that the
// minimal (shortest non-zero) vectors have squared length 2, hence minimal
// distance sqrt(2) and packing radius 1/sqrt(2); the Leech lattice is scaled to
// minimal squared length 4, minimal distance 2 and packing radius 1. The center
// density is delta = rho^n / V where rho is the packing radius and V is the
// covolume (the volume of a fundamental domain), and the packing density is
// Delta = delta * V_n where V_n is the volume of the unit n-ball. The covering
// thickness is Theta = V_n * R^n / V where R is the covering radius.
package packing
