package graphspectral

import "errors"

// Sentinel errors returned throughout the package. Callers may test for these
// with errors.Is.
var (
	// ErrDimensionMismatch reports that two operands have incompatible shapes.
	ErrDimensionMismatch = errors.New("graphspectral: dimension mismatch")
	// ErrNotSquare reports that a square matrix was required but not provided.
	ErrNotSquare = errors.New("graphspectral: matrix is not square")
	// ErrNotSymmetric reports that a symmetric matrix was required but not provided.
	ErrNotSymmetric = errors.New("graphspectral: matrix is not symmetric")
	// ErrSingular reports that a matrix is singular (or numerically so).
	ErrSingular = errors.New("graphspectral: matrix is singular")
	// ErrEmpty reports that an empty matrix, vector or graph was supplied.
	ErrEmpty = errors.New("graphspectral: empty input")
	// ErrOutOfRange reports an index outside the valid range.
	ErrOutOfRange = errors.New("graphspectral: index out of range")
	// ErrInvalidArgument reports a value outside its permitted domain.
	ErrInvalidArgument = errors.New("graphspectral: invalid argument")
	// ErrNotConnected reports that a connected graph was required but the graph
	// has more than one connected component.
	ErrNotConnected = errors.New("graphspectral: graph is not connected")
	// ErrNoConvergence reports that an iterative method failed to converge
	// within its iteration budget.
	ErrNoConvergence = errors.New("graphspectral: iteration did not converge")
)
