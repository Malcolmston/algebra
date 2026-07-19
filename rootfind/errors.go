package rootfind

import "errors"

// ErrNoConvergence is returned by an iterative solver that failed to reach the
// requested tolerance within its iteration budget.
var ErrNoConvergence = errors.New("rootfind: iteration did not converge")

// ErrNoBracket is returned by a bracketing solver when the supplied endpoints
// do not straddle a root, i.e. the function has the same sign at both ends.
var ErrNoBracket = errors.New("rootfind: endpoints do not bracket a root")

// ErrZeroDerivative is returned by a derivative-based solver when the derivative
// vanishes at the current iterate, preventing a well-defined update step.
var ErrZeroDerivative = errors.New("rootfind: derivative vanished during iteration")

// ErrDegreeTooLow is returned when a routine requires a polynomial of at least
// a certain degree (for example a nonconstant polynomial) and it is not met.
var ErrDegreeTooLow = errors.New("rootfind: polynomial degree too low for this operation")

// ErrZeroPolynomial is returned when an operation that is undefined on the zero
// polynomial (such as division or making a polynomial monic) receives one.
var ErrZeroPolynomial = errors.New("rootfind: operation undefined for the zero polynomial")

// ErrEmptyInterval is returned by interval routines when the supplied interval
// is degenerate or has its endpoints reversed.
var ErrEmptyInterval = errors.New("rootfind: invalid or empty interval")

// ErrBadInput is returned for malformed arguments that are not covered by a more
// specific sentinel error, such as a negative tolerance or iteration count.
var ErrBadInput = errors.New("rootfind: invalid input")
