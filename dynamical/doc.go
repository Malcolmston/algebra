// Package dynamical provides primitives for the study of discrete and
// continuous dynamical systems, implemented with the Go standard library only.
//
// The package is organized around a handful of small concepts:
//
//   - Iterated one-dimensional maps ([Map1D]) such as the logistic, tent,
//     sine, Gauss, cubic and circle maps, together with helpers to build
//     parameterized families as closures.
//   - Iterated two-dimensional maps ([Map2D]) over [Point2D] values, most
//     notably the Henon map.
//   - Orbits and trajectories: forward iteration with optional transient
//     discarding, and single n-th iterate evaluation.
//   - Fixed points and their linear stability, including closed-form fixed
//     points and multipliers for the logistic and tent maps.
//   - Lyapunov exponents estimated both from the analytic derivative and by
//     the trajectory-separation (renormalization) method, in one and two
//     dimensions.
//   - Bifurcation-diagram sampling for parameterized one-dimensional families.
//   - Period detection for periodic and eventually-periodic orbits.
//   - Cobweb (staircase) diagram segment data for visualizing 1-D iteration.
//   - Newton's method in the real and complex plane, and the basin-of-
//     attraction grid produced by Newton iteration of a complex map.
//   - Continuous flows in three dimensions ([Flow3D]) integrated with a
//     fixed-step classical Runge-Kutta (RK4) scheme, with the Lorenz and
//     Rossler systems provided as concrete examples.
//
// All routines are deterministic: given the same inputs they produce the same
// outputs, and no global random state is used.
package dynamical
