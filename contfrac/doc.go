// Package contfrac is a self-contained toolkit for continued fractions and
// Diophantine approximation, built entirely on the Go standard library
// (math, math/big, sort, errors, fmt, strings).
//
// It is a sibling subpackage of github.com/malcolmston/algebra but does not
// depend on it. Everything is expressed in terms of a small set of value
// types defined here:
//
//   - [Frac] is an exact rational with int64 numerator and denominator, with
//     the usual field of arithmetic and ordering operations plus conversions
//     to and from continued fractions.
//   - [CF] is a finite continued fraction [a0; a1, a2, ...] stored as a slice
//     of partial quotients. The Euclidean/floor convention is used, so a0 may
//     be zero or negative but every later term is a positive integer.
//   - [QuadraticSurd] represents (P + sqrt(D))/Q and [PeriodicCF] represents an
//     eventually periodic continued fraction; the two convert to one another,
//     realising the classical theorem of Lagrange.
//
// # Contents
//
// Core expansion and evaluation: [FromRational], [FromFloat], [FromBigFloat],
// [Convergents], [BestApproximation], [Continuant].
//
// Quadratic irrationals: [SqrtCF], [SqrtCFPeriod], [QuadraticSurd.CF],
// [PeriodicCF.Surd], and the Pell-equation solvers [PellFundamental],
// [PellSolutions], [PellNegative].
//
// Trees and sequences of fractions: [Mediant], [SternBrocotPath],
// [SternBrocotFromPath], [FareySequence], [FareySuccessor], [PathToCF].
//
// Decompositions and constants: [EgyptianFraction], [ECF], [PiCF],
// [GoldenRatioCF].
//
// All floating-point routines document their tolerance and every exact routine
// is closed over the integers or math/big, so results are reproducible.
package contfrac
