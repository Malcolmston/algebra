// Package seq implements classical integer sequences and figurate numbers.
//
// The package is a self-contained catalogue of the recurrences, closed forms
// and iterative procedures that recur throughout elementary number theory and
// recreational mathematics. It is organised into three broad families:
//
//   - Linear recurrences and their relatives: the Fibonacci numbers computed by
//     the fast-doubling method, together with the Lucas, Pell, Pell-Lucas,
//     Jacobsthal, Jacobsthal-Lucas, Tribonacci, Tetranacci, Padovan and Perrin
//     sequences, a fully general constant-coefficient LinearRecurrence engine
//     and the Zeckendorf representation.
//
//   - Figurate numbers: the polygonal numbers (triangular, square, pentagonal,
//     hexagonal, heptagonal, octagonal and the general s-gonal family), the
//     centered polygonal numbers, the three-dimensional figurate numbers
//     (tetrahedral, square-pyramidal), the four-dimensional pentatope numbers
//     and assorted companions (pronic, star, cube and gnomonic numbers), each
//     paired where natural with a membership predicate.
//
//   - Digit- and iteration-driven sequences: the Collatz trajectory, happy
//     numbers, Kaprekar numbers and the Kaprekar 6174 routine, Recaman's
//     sequence and the look-and-say sequence.
//
// Machine-word routines return uint64 (or int64 where negative terms arise) and
// document the index beyond which they overflow; arbitrary-precision *big.Int
// variants are provided for the Fibonacci-family sequences where exact values
// far beyond the 64-bit range are wanted. Every routine is deterministic and
// depends only on the Go standard library.
package seq

import "math"

// IntSqrt returns the integer square root of v, that is the largest integer r
// with r*r <= v. It is exact for every uint64 input.
func IntSqrt(v uint64) uint64 {
	if v == 0 {
		return 0
	}
	r := uint64(math.Sqrt(float64(v)))
	// Correct any floating-point rounding error in either direction. The
	// guards keep the (r+1)*(r+1) probe from overflowing near 2^64.
	for r > 0 && r*r > v {
		r--
	}
	for r < 0xFFFFFFFF && (r+1)*(r+1) <= v {
		r++
	}
	return r
}

// IsPerfectSquare reports whether v is a perfect square.
func IsPerfectSquare(v uint64) bool {
	r := IntSqrt(v)
	return r*r == v
}
