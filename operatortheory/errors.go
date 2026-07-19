package operatortheory

import "errors"

// Sentinel errors returned throughout the package. Callers may test for these
// with errors.Is.
var (
	// ErrDimensionMismatch reports that two operands have incompatible shapes.
	ErrDimensionMismatch = errors.New("operatortheory: dimension mismatch")
	// ErrNotSquare reports that a square matrix was required but not provided.
	ErrNotSquare = errors.New("operatortheory: matrix is not square")
	// ErrNotHermitian reports that a Hermitian matrix was required but not
	// provided.
	ErrNotHermitian = errors.New("operatortheory: matrix is not Hermitian")
	// ErrNotNormal reports that a normal matrix was required but not provided.
	ErrNotNormal = errors.New("operatortheory: matrix is not normal")
	// ErrSingular reports that a matrix is singular (or numerically so).
	ErrSingular = errors.New("operatortheory: matrix is singular")
	// ErrEmpty reports that an empty matrix or vector was supplied.
	ErrEmpty = errors.New("operatortheory: empty input")
	// ErrOutOfRange reports an index outside the valid range.
	ErrOutOfRange = errors.New("operatortheory: index out of range")
	// ErrInvalidArgument reports a value outside its permitted domain.
	ErrInvalidArgument = errors.New("operatortheory: invalid argument")
	// ErrNoConvergence reports that an iterative method failed to converge
	// within its iteration budget.
	ErrNoConvergence = errors.New("operatortheory: iteration did not converge")
)
