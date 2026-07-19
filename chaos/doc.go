// Package chaos implements tools for the study of chaos and nonlinear
// dynamics using only the Go standard library.
//
// The package collects, in a single self-contained place, the computational
// primitives that recur throughout the study of deterministic chaos:
//
//   - Canonical systems: the one-dimensional logistic, tent, sine, Gauss and
//     cubic maps ([Map1D]); the two-dimensional Henon and standard (Chirikov)
//     maps ([Map2D]); and the continuous Lorenz and Rossler flows
//     ([Flow]) integrated with fixed-step Runge-Kutta schemes.
//
//   - Fixed points and linear stability: root finding for map and flow fixed
//     points, numerical Jacobians, eigenvalues of small matrices, and the
//     classification of equilibria into nodes, saddles, foci and centres.
//
//   - Lyapunov exponents: the largest exponent by the Benettin
//     single-trajectory method and by direct trajectory separation, and the
//     full spectrum by the QR (Gram-Schmidt) reorthonormalisation method for
//     both maps and flows.
//
//   - Bifurcation analysis: bifurcation-diagram sampling for parameterised
//     one-dimensional families, period-doubling cascade detection and the
//     estimation of the Feigenbaum constants delta and alpha.
//
//   - Poincare sections and return maps: hyperplane crossings of a flow with
//     linear interpolation to the section, and first-return maps.
//
//   - Fractal dimensions: box-counting (capacity) dimension, the
//     Grassberger-Procaccia correlation dimension, information dimension and
//     the Kaplan-Yorke (Lyapunov) dimension.
//
// Vectors are represented by the [Vec] type (a slice of float64) and matrices
// by the [Mat] type (a slice of rows). Discrete maps are ordinary Go closures
// of type [Map1D] or [MapN]; continuous systems are described by a right-hand
// side of type [Field] wrapped in a [Flow].
//
// Randomness, where used (for example to seed ensembles of initial
// conditions), is always drawn from a caller-supplied *math/rand.Rand so that
// every routine is reproducible. The implementations depend only on the
// packages math, math/cmplx, sort, errors, fmt and strings and favour
// clarity and numerical soundness over raw speed.
package chaos
