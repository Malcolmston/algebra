// Package complexanalysis provides numerical tools for complex analysis using
// only the Go standard library.
//
// It is organised into a few areas:
//
//   - Elementary and inverse functions on complex128 with well-defined
//     principal branches (see elementary.go): [Sqrt], [Cbrt], [Exp], [Log],
//     [Pow], the trigonometric and hyperbolic families, and their inverses.
//   - Special functions (see special.go): the complex [Gamma] and [LogGamma]
//     via the Lanczos approximation, [Digamma], [Beta], the error functions
//     [Erf]/[Erfc], the Riemann [Zeta] function via Borwein's algorithm, and
//     the combinatorial [Factorial], [RisingFactorial] and [Binomial].
//   - Numerical contour integration and its consequences (see contour.go):
//     [IntegrateCircle], [IntegrateSegment], [IntegratePolygon] and the general
//     [ContourIntegral]; [Residue], [ResidueSimplePole] and [ResidueOrderM];
//     the [CauchyIntegralValue] and [CauchyDerivative] formulas; the
//     [WindingNumber]; and the [ArgumentPrinciple] with [CountZeros].
//   - Conformal maps (see conformal.go): the [Mobius] transformation type with
//     composition, inversion and fixed points, together with the [CrossRatio],
//     the [CayleyTransform] and the [JoukowskiMap].
//   - Laurent and Taylor coefficient extraction and numeric analytic
//     continuation (see laurent.go): [LaurentCoefficient], [TaylorCoefficient],
//     [PowerSeriesEval] and [AnalyticContinuation].
//
// Every function that samples an analytic function does so by evaluating a
// user-supplied [Function]. Contour methods use the composite trapezoidal rule
// on a circular contour, which converges geometrically for functions that are
// analytic in an annulus, so modest sample counts (a few hundred points) give
// near machine-precision results.
//
// The package is deterministic and depends only on math and math/cmplx.
package complexanalysis

// Function is a complex-valued function of a single complex variable. It is the
// common argument type for the numerical routines in this package.
type Function func(complex128) complex128
