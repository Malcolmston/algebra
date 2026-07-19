package chaos

import "errors"

var (
	// ErrDimensionMismatch is returned when two vectors or matrices have
	// incompatible shapes for the requested operation.
	ErrDimensionMismatch = errors.New("chaos: dimension mismatch")

	// ErrEmpty is returned when an operation requires a non-empty input but
	// received one of length zero.
	ErrEmpty = errors.New("chaos: empty input")

	// ErrNonSquare is returned when a square matrix was required.
	ErrNonSquare = errors.New("chaos: matrix is not square")

	// ErrSingular is returned when a matrix is singular to working precision.
	ErrSingular = errors.New("chaos: matrix is singular")

	// ErrNoConvergence is returned when an iterative routine fails to reach
	// the requested tolerance within the allotted iterations.
	ErrNoConvergence = errors.New("chaos: iteration did not converge")

	// ErrBadParameter is returned when a numeric parameter is out of its
	// valid range (for example a negative step count).
	ErrBadParameter = errors.New("chaos: invalid parameter")

	// ErrNoCrossing is returned when a Poincare section produced no crossings.
	ErrNoCrossing = errors.New("chaos: no section crossings found")
)
