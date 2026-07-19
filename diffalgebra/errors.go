package diffalgebra

import "errors"

// Sentinel errors returned throughout the package.
var (
	// ErrDivByZero indicates division by a zero polynomial or rational function.
	ErrDivByZero = errors.New("diffalgebra: division by zero")
	// ErrDim indicates a shape or dimension mismatch between operands.
	ErrDim = errors.New("diffalgebra: dimension mismatch")
	// ErrNotSquare indicates that a square matrix was required.
	ErrNotSquare = errors.New("diffalgebra: matrix is not square")
	// ErrSingular indicates a singular (non-invertible) linear system.
	ErrSingular = errors.New("diffalgebra: singular linear system")
	// ErrEmpty indicates that an empty argument was supplied where at least one
	// element is required.
	ErrEmpty = errors.New("diffalgebra: empty argument")
	// ErrDegree indicates an unsupported or inconsistent polynomial degree.
	ErrDegree = errors.New("diffalgebra: unsupported degree")
	// ErrNoSolution indicates that no solution of the requested kind exists.
	ErrNoSolution = errors.New("diffalgebra: no solution")
	// ErrNonRational indicates that a residue or exponent was not rational and
	// could therefore not be expressed exactly over Q.
	ErrNonRational = errors.New("diffalgebra: non-rational value")
	// ErrConverge indicates that a numerical iteration failed to converge.
	ErrConverge = errors.New("diffalgebra: failed to converge")
)
