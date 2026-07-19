// Package quadrature implements numerical integration ("quadrature") in pure
// Go using only the standard library.
//
// The package is organized around three complementary ideas that together
// cover the classical toolbox of one- and multi-dimensional numerical
// integration:
//
//   - Gaussian quadrature. Nodes and weights for the Gauss rules associated
//     with the classical weight functions are generated from the three-term
//     recurrence of the underlying orthogonal polynomials via the
//     Golub-Welsch algorithm (an eigenvalue problem for a symmetric
//     tridiagonal Jacobi matrix). Legendre, Chebyshev (four kinds), Hermite
//     (physicists' and probabilists'), Laguerre (including the generalized
//     form) and Jacobi families are provided, together with the Gauss-Lobatto
//     and Gauss-Radau variants that fix one or both endpoints and an embedded
//     Gauss-Kronrod pair for error estimation.
//
//   - Interpolatory and iterative rules. The Newton-Cotes family (midpoint,
//     trapezoid, Simpson, Simpson's 3/8 and Boole) is provided in single-panel
//     and composite forms, with exact rational weight generation for arbitrary
//     order. Romberg extrapolation, adaptive Simpson and adaptive
//     Gauss-Kronrod refinement, the tanh-sinh (double-exponential) rule for
//     endpoint singularities and semi-infinite/infinite ranges, and the
//     Clenshaw-Curtis and Fejer rules built on Chebyshev points round out the
//     one-dimensional methods.
//
//   - Multi-dimensional integration. Tensor-product Gauss and Newton-Cotes
//     rules over rectangular boxes, together with plain, stratified and
//     quasi-random (Halton) Monte-Carlo estimators for higher dimensions.
//
// A [Rule] value bundles a slice of nodes with the matching weights and can be
// evaluated against a function, rescaled from the canonical interval to an
// arbitrary [a, b], or composed into product rules. Deterministic routines are
// free of global state; the Monte-Carlo routines take an explicit seed so that
// results are reproducible.
//
// All abscissae for the finite Gauss and interpolatory rules are returned in
// ascending order. Unless otherwise noted, weight-function rules such as
// [GaussHermite] or [GaussChebyshev1] return the nodes and weights for the
// integral of f against that weight; the corresponding IntegrateGaussXxx
// helpers apply them.
package quadrature
