// Package autodiff implements automatic differentiation in pure Go using only
// the standard library.
//
// Automatic differentiation (AD) evaluates the exact derivative of a function
// specified by ordinary program code, up to floating-point round-off. Unlike
// finite differences it introduces no truncation error and needs no step size,
// and unlike symbolic differentiation it does not build an expression graph
// that grows with the algebra. The package provides three complementary
// engines:
//
//   - Forward mode via [Dual] numbers. A dual number carries a value together
//     with its derivative along a chosen direction. Evaluating a function on a
//     seeded dual number yields the value and the directional derivative in a
//     single pass. Forward mode is efficient when the number of inputs is small
//     relative to the number of outputs and underlies [Derivative], [Gradient],
//     [Jacobian], [DirectionalDerivative] and [PartialDerivative].
//
//   - Second-order forward mode via [HyperDual] numbers. A hyper-dual number
//     carries a value, two first-order derivative slots and one mixed
//     second-order slot, so a single evaluation returns both first and second
//     derivatives with no cancellation error. This underlies [SecondDerivative],
//     [Hessian], [GradientHessian] and [HessianVectorProduct].
//
//   - Reverse mode via a small [Tape]. The tape records every elementary
//     operation as it executes, then propagates adjoints backward from the
//     output to every input in one sweep, yielding the whole gradient of a
//     scalar function at the cost of roughly one function evaluation regardless
//     of the number of inputs.
//
// All routines are deterministic and free of global mutable state (each [Tape]
// is independent), so results are reproducible across runs and safe to use from
// distinct goroutines that operate on distinct tapes.
//
// # Seeding
//
// Forward-mode derivatives are obtained by seeding: the derivative slot of the
// independent variable is set to one (see [Variable]) while constants keep a
// zero slot (see [Constant]). The chain rule then transports the seed through
// every elementary operation, and the derivative appears in the output's slot.
package autodiff
