// Package discretemath implements a broad collection of discrete-mathematics
// primitives built entirely on the Go standard library.
//
// The package is deliberately self-contained and deterministic: identical
// inputs always yield identical outputs, no global state is mutated after
// initialization, and no third-party dependency is required.
//
// The following families are covered:
//
//   - bit tricks: population count, parity, Gray-code encode/decode, Morton
//     (Z-order) interleaving in two and three dimensions, bit reversal,
//     rotation, power-of-two tests, integer base-2 logarithms, leading and
//     trailing zero counts (the latter implemented with a De Bruijn table),
//     and single-bit set/clear/toggle/test helpers;
//   - Hamming distance over machine words, byte slices and strings;
//   - De Bruijn sequences B(k, n) over integer alphabets and over arbitrary
//     rune alphabets;
//   - base conversion for any radix in the range 2..36, in both signed and
//     unsigned forms;
//   - Roman-numeral encoding and decoding (1..3999);
//   - English cardinal and ordinal number-to-words conversion;
//   - set operations over comparable element types: union, intersection,
//     difference, symmetric difference, subset and disjointness tests, and the
//     power set of a slice;
//   - Cartesian products of two heterogeneous slices and of an arbitrary number
//     of homogeneous slices;
//   - run-length encoding and decoding of slices and of ASCII strings;
//   - edit distances: Levenshtein and the optimal-string-alignment variant of
//     Damerau-Levenshtein, plus the longest common subsequence.
//
// Every exported identifier is documented and every routine is validated in the
// accompanying tests against closed-form or reference values.
package discretemath
