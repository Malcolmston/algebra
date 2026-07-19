// Package approxtheory implements the core tools of classical approximation
// theory using nothing but the Go standard library.
//
// The package is organised around several themes:
//
//   - Chebyshev series: computing interpolation coefficients from function
//     samples (a discrete cosine transform), evaluating a series with the
//     numerically stable Clenshaw recurrence, differentiating and integrating
//     in coefficient space, finding the real roots on an interval and doing
//     ordinary arithmetic on series.
//
//   - Remez exchange: computing the best (minimax) polynomial approximation of
//     a continuous function on an interval, together with the equioscillating
//     error.
//
//   - Pade approximants: turning a truncated Taylor series into a rational
//     function that usually reproduces the original function far better than
//     the polynomial partial sum.
//
//   - Polynomial fitting: ordinary and weighted least-squares fits together
//     with a discrete minimax (Chebyshev) fit obtained by an exchange
//     algorithm.
//
//   - Rational interpolation: Thiele's continued-fraction interpolant built
//     from reciprocal differences.
//
//   - Bernstein polynomials and Bezier evaluation, degree elevation and the
//     Bernstein approximation operator.
//
//   - Error analysis: Lebesgue functions and constants, node polynomials and
//     the associated interpolation error bounds, and assorted error norms.
//
// Unless otherwise noted, polynomials are represented in the monomial basis by
// a slice of coefficients in ascending order, so that coeffs[i] multiplies
// x**i. Chebyshev series carry their own domain [A, B]; all other routines act
// on the interval the caller supplies.
//
// The implementations favour clarity and numerical soundness over raw speed
// and depend only on the packages math, sort, errors, fmt and strings.
package approxtheory
