package odesolvers

import "errors"

// Sentinel errors returned by the solvers. Callers may test for them with
// errors.Is.
var (
	// ErrDimensionMismatch reports inconsistent vector or matrix dimensions.
	ErrDimensionMismatch = errors.New("odesolvers: dimension mismatch")
	// ErrMaxSteps reports that an integrator reached its step-count limit
	// before covering the requested interval.
	ErrMaxSteps = errors.New("odesolvers: maximum number of steps exceeded")
	// ErrStepTooSmall reports that adaptive step-size control drove the step
	// below the allowed minimum, usually a sign of a discontinuity or an
	// excessively tight tolerance.
	ErrStepTooSmall = errors.New("odesolvers: required step size underflowed the minimum")
	// ErrSingularMatrix reports a singular (or numerically singular) linear
	// system encountered during an implicit solve.
	ErrSingularMatrix = errors.New("odesolvers: singular matrix")
	// ErrNoConvergence reports that an iterative solver (Newton, shooting)
	// failed to converge within its iteration budget.
	ErrNoConvergence = errors.New("odesolvers: iteration failed to converge")
	// ErrInvalidOrder reports an out-of-range method order request.
	ErrInvalidOrder = errors.New("odesolvers: invalid method order")
	// ErrInvalidInput reports malformed or inconsistent solver inputs.
	ErrInvalidInput = errors.New("odesolvers: invalid input")
)
