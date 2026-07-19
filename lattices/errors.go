package lattices

import "errors"

// ErrDimMismatch is returned when two vectors or matrices have incompatible
// dimensions for the requested operation.
var ErrDimMismatch = errors.New("lattices: dimension mismatch")

// ErrEmpty is returned when an operation requires a non-empty basis, vector or
// matrix but received an empty one.
var ErrEmpty = errors.New("lattices: empty input")

// ErrNotSquare is returned when a square matrix is required but the argument is
// rectangular.
var ErrNotSquare = errors.New("lattices: matrix is not square")

// ErrSingular is returned when a matrix is singular and therefore cannot be
// inverted or used to solve a linear system.
var ErrSingular = errors.New("lattices: matrix is singular")

// ErrNotFullRank is returned when a basis is expected to be full rank (its
// vectors linearly independent) but is rank deficient.
var ErrNotFullRank = errors.New("lattices: basis is not full rank")

// ErrRankMismatch is returned when the number of basis vectors does not match
// the ambient dimension for an operation that requires a square basis.
var ErrRankMismatch = errors.New("lattices: basis rank does not match ambient dimension")

// ErrBadParameter is returned when a numeric parameter is outside its valid
// range (for example an LLL delta outside (1/4, 1]).
var ErrBadParameter = errors.New("lattices: parameter out of range")

// ErrNoSolution is returned when an enumeration or search finds no vector
// satisfying the requested constraints.
var ErrNoSolution = errors.New("lattices: no solution found")
