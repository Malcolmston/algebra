// Package odesolvers provides numerical integrators for systems of
// first-order ordinary differential equations (ODEs) of the form
//
//	y'(t) = f(t, y),   y(t0) = y0,
//
// implemented using only the Go standard library.
//
// The package collects the classical families of ODE methods behind a small,
// uniform interface built around the [Field] function type and the [Solution]
// container:
//
//   - Explicit Runge-Kutta methods described by a [ButcherTableau]: forward
//     Euler, the explicit midpoint and Heun methods, Ralston's method, the
//     classical fourth-order RK4 and the 3/8-rule, and SSPRK3.
//   - Embedded adaptive Runge-Kutta pairs with automatic step-size control:
//     Heun-Euler 2(1), Bogacki-Shampine 3(2), Runge-Kutta-Fehlberg RKF45,
//     Cash-Karp, Dormand-Prince DOPRI5 5(4) and the high-order Fehlberg
//     RKF78 7(8) pair.
//   - Implicit / stiff solvers built on a generic implicit Runge-Kutta stage
//     solver: backward Euler, the trapezoidal rule, the implicit midpoint and
//     Gauss-Legendre methods, and the Radau IIA methods of orders 3 and 5.
//   - Fixed-step backward-differentiation formulas (BDF1 through BDF6).
//   - Linear multistep methods: Adams-Bashforth, Adams-Moulton and their
//     Adams-Bashforth-Moulton predictor-corrector combination.
//   - Symplectic integrators for separable second-order systems
//     q” = a(t, q): the velocity- and position-Verlet schemes, the
//     leapfrog method and the fourth-order Yoshida composition.
//   - Event / root detection during integration ([ScanEvents]).
//   - Boundary-value problems via single and multiple shooting.
//
// All solvers are deterministic. Where a pseudo-random source is required the
// caller supplies a seed and a math/rand source is used; crypto/rand and the
// wall clock are never consulted.
//
// The [Field] type returns a freshly allocated derivative slice on each call
// and never mutates its input, which keeps the method implementations simple
// and free of hidden aliasing.
package odesolvers
