// Package settheory implements finite set theory, binary relations, order
// theory and lattices over integer and string domains, built entirely on the
// Go standard library.
//
// The package is self-contained and deterministic: identical inputs always
// yield identical outputs, iteration-order sensitivity is eliminated by
// returning sorted slices from every accessor, and no third-party dependency
// is required.
//
// The following families are covered:
//
//   - finite sets: an IntSet and a StringSet backed by hash maps, providing
//     membership, cardinality, union, intersection, difference, symmetric
//     difference, subset/superset/disjointness tests, Cartesian product and
//     power set enumeration;
//   - binary relations on an integer domain: composition, inverse, the
//     reflexive/symmetric/transitive property tests, and the reflexive,
//     symmetric, transitive (via Warshall's algorithm) and equivalence
//     closures;
//   - partitions: equivalence classes induced by an equivalence relation, the
//     relation induced by a partition, refinement testing, and the Bell and
//     Stirling-of-the-second-kind counting numbers;
//   - partial orders (posets): comparability, covering (Hasse) edges, a
//     deterministic topological order (linear extension), minimal/maximal and
//     least/greatest elements, upper and lower bounds, chain height;
//   - lattices: pairwise meet (greatest lower bound) and join (least upper
//     bound), lattice detection, and the divisibility lattice of a natural
//     number.
//
// All construction functions validate their inputs and report violations as
// errors rather than panicking.
package settheory

import "fmt"

// settheoryErrorf builds a formatted error. It centralizes error creation so
// that every exported routine reports failures in a consistent style.
func settheoryErrorf(format string, args ...any) error {
	return fmt.Errorf(format, args...)
}
