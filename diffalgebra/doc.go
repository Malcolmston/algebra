// Package diffalgebra provides symbolic tools for differential algebra and
// linear ordinary differential equations over the rational numbers.
//
// The package is built up from a small tower of exact algebraic types and a
// large collection of genuinely distinct operations on them:
//
//   - Exact rational numbers via math/big and dense univariate polynomials
//     [Poly] over Q with the full complement of ring operations: addition,
//     multiplication, Euclidean division, greatest common divisors, extended
//     GCD, resultants, discriminants, square-free factorisation, rational-root
//     isolation, composition and formal differentiation and integration.
//
//   - Rational functions [RatFunc] = Poly/Poly forming the differential field
//     Q(x) with the standard derivation d/dx ([Derivation]), logarithmic
//     derivatives, partial-fraction decomposition and exact arithmetic.
//
//   - Linear differential operators [Operator] with rational-function
//     coefficients, forming the non-commutative Ore ring Q(x)[D] with the
//     relation D*a = a*D + a'. Operators can be added, composed, applied to
//     functions, raised to powers and reduced to their formal adjoint and
//     symbol.
//
//   - Wronskians ([WronskianPoly], [WronskianRatFunc]) and the induced
//     linear-independence tests, built on exact determinants over Q(x).
//
//   - Symbolic solution of constant-coefficient linear ODEs and recurrences
//     ([SolveLinearConstantODE], [SolveLinearRecurrence]) via numerically
//     computed characteristic roots grouped into exact multiplicities, together
//     with initial-value fitting by exact/complex linear algebra.
//
//   - The variation-of-parameters construction ([VariationOfParameters]) that
//     turns a fundamental system into the integrands of a particular solution.
//
//   - Elementary integration of rational functions by Hermite reduction
//     ([HermiteReduce]) and the Rothstein-Trager resultant method
//     ([IntegrateRational]), plus the Risch structure-theorem heuristics for
//     simple exponential integrands ([RischExpIntegrate]).
//
//   - The Kovacic algorithm ([Kovacic]) for detecting Liouvillian solutions of
//     second-order linear ODEs y” = r y with r in Q(x), including reduction of
//     a general second-order equation to normal form ([ReduceToNormalForm]).
//
// Everything is implemented with the Go standard library only. Symbolic data is
// exact over Q; numerical routines (root finding, initial-value fitting) accept
// a caller-supplied seed for reproducibility and never read the wall clock.
package diffalgebra
