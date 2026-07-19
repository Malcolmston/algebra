// Package designs implements combinatorial designs and finite geometry.
//
// The package provides exact, integer-arithmetic constructions and analysers
// for the classical objects of design theory:
//
//   - Balanced incomplete block designs: the parameter quintuple
//     [BIBDParams] (v, b, r, k, lambda) with the necessary divisibility
//     conditions, Fisher's inequality, symmetric designs and their
//     complements, derived and residual designs; a general [Design] type that
//     stores an incidence structure and recovers its parameters, replication
//     numbers, concurrence (pairwise lambda) matrix, dual and complement.
//   - Latin squares and mutually orthogonal Latin squares ([LatinSquare],
//     [MOLS]), Graeco-Latin squares and orthogonality testing, built directly
//     and from finite fields.
//   - Steiner triple and quadruple systems ([SteinerTripleSystem],
//     [SteinerQuadrupleSystem]) via the Bose construction, Stinson's
//     hill-climbing algorithm and the Boolean construction on GF(2)^k.
//   - Resolvable designs, parallel classes, one-factorizations of complete
//     graphs and round-robin tournament schedules ([OneFactorization]).
//   - Hadamard matrices ([HadamardMatrix]) by the Sylvester and Paley
//     constructions, conference matrices, and the symmetric 2-designs they
//     produce.
//   - Difference sets ([DifferenceSet]) in cyclic groups, the Paley
//     (quadratic-residue) and Singer planar families, and development into
//     symmetric designs.
//   - Finite projective and affine planes ([ProjectivePlane], [AffinePlane])
//     of prime-power order with full incidence, joins, meets and parallel
//     classes.
//
// Finite fields GF(q) of arbitrary prime-power order are provided by
// [GaloisField], built on an irreducible polynomial found over the prime
// field.
//
// All computation is exact (int, and math/big only where large exponents
// arise). The package depends only on the Go standard library. Every routine
// that requires randomness takes a caller-supplied *math/rand.Rand so that
// results are fully reproducible from a seed; no global, time-based or
// cryptographic source of randomness is used.
package designs
