package liealgebra

import "errors"

// Sentinel errors returned throughout the package.
var (
	// ErrDim indicates a shape or dimension mismatch between operands.
	ErrDim = errors.New("liealgebra: dimension mismatch")
	// ErrNotSquare indicates that a square matrix was required.
	ErrNotSquare = errors.New("liealgebra: matrix is not square")
	// ErrSingular indicates a singular (non-invertible) linear system.
	ErrSingular = errors.New("liealgebra: singular matrix")
	// ErrRank indicates a rank-deficient or linearly dependent basis.
	ErrRank = errors.New("liealgebra: rank deficient basis")
	// ErrType indicates an unknown or unsupported Dynkin/Cartan type.
	ErrType = errors.New("liealgebra: unknown Dynkin type")
	// ErrRange indicates an out-of-range rank or index argument.
	ErrRange = errors.New("liealgebra: argument out of range")
)
