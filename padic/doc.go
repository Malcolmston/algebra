// Package padic is a self-contained toolkit for p-adic numbers and p-adic
// arithmetic, built entirely on the Go standard library (math, math/big,
// math/rand, sort, errors, fmt, strings).
//
// It is a sibling subpackage of github.com/malcolmston/algebra but does not
// depend on it. Every value is expressed in terms of the type [Padic], a
// p-adic number carried to a finite, tracked precision.
//
// # Representation
//
// A non-zero p-adic number x is written uniquely as
//
//	x = p^v * u
//
// where v = [Padic.Valuation] is an integer (possibly negative) and u is a
// p-adic unit, i.e. an integer coprime to p. This package stores u modulo
// p^r, where r = [Padic.RelativePrecision] is the number of known p-adic
// digits of the unit. The absolute precision, [Padic.AbsolutePrecision], is
// v + r: the element is known modulo p^(v+r). The p-adic zero is represented
// as "zero to some absolute precision" and carries no unit; [Padic.IsZero]
// detects it.
//
// Arithmetic tracks precision honestly. Adding two numbers yields a result
// whose absolute precision is the smaller of the inputs; a subtraction that
// cancels leading digits correctly reports the reduced relative precision
// (catastrophic cancellation is visible, not hidden).
//
// # Contents
//
// Construction: [New], [Zero], [One], [FromInt], [FromInt64], [FromBigInt],
// [FromRat], [FromRational], [Uniformizer], [TeichmullerUnit].
//
// Field arithmetic: [Padic.Add], [Padic.Sub], [Padic.Neg], [Padic.Mul],
// [Padic.Div], [Padic.Inv], [Padic.Pow], [Padic.Square].
//
// Valuation and absolute value: [ValuationInt], [ValuationRat],
// [Padic.Valuation], [Padic.AbsValue], [AbsValueRat], [Padic.NormFloat].
//
// Teichmuller theory: [Teichmuller], [TeichmullerRep], [Padic.Teichmuller],
// [IsTeichmuller].
//
// Hensel lifting of polynomial roots: [HenselLift], [HenselLiftPadic],
// [SimpleRootsModP], [PadicRoots].
//
// p-adic square roots, logarithm and exponential: [Padic.Sqrt],
// [Padic.IsSquare], [SqrtInt], [Padic.Log], [Padic.Exp].
//
// Newton polygons: [NewtonPolygonFromInts], [NewtonPolygonFromRats],
// [NewtonPolygonFromPadics] and the methods on [NewtonPolygon].
//
// p-adic expansion of rationals: [ExpandRational], [Padic.Digits],
// [DigitsToRat].
//
// Convergence helpers: [StrassmannBound], [SeriesValuationAt],
// [ConvergesInUnitBall].
//
// Number-theoretic utilities: [IsPrime], [NextPrime], [LegendreSymbol],
// [JacobiSymbol], [SqrtModP], [PowMod], [InvMod].
package padic
