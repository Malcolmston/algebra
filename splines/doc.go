// Package splines implements splines and parametric curves using nothing but
// the Go standard library.
//
// The package provides a broad, self-contained toolkit for computer-aided
// geometric design and numerical interpolation:
//
//   - Cubic splines with natural, clamped, not-a-knot and periodic end
//     conditions, together with value, derivative and integral evaluation.
//   - Piecewise cubic Hermite interpolation, Catmull-Rom / cardinal splines,
//     shape-preserving monotone interpolation (Fritsch-Carlson) and Akima
//     splines.
//   - Bezier curves evaluated with de Casteljau's algorithm, plus degree
//     elevation, subdivision, derivatives and the Bernstein basis.
//   - B-spline curves evaluated with de Boor's algorithm, knot insertion
//     (Boehm), the Cox-de Boor basis and derivatives.
//   - NURBS curves (rational B-splines) with control point weights.
//   - Arc-length computation and arc-length (unit-speed) reparameterisation.
//   - Tensor-product Bezier and B-spline surface evaluation.
//
// Points and vectors are represented by the [Vec] type, an alias-free slice of
// float64 coordinates, so every routine works uniformly in any dimension.
// Scalar (one dimensional) interpolation is available through dedicated
// float64 helpers that avoid the small allocation overhead of [Vec].
//
// All algorithms are implemented from first principles with no third-party
// dependencies and no cgo.
package splines
