package exterioralgebra

import "errors"

// ErrDim is returned when two operands live in exterior algebras of different
// ambient dimension, or when a supplied dimension is invalid.
var ErrDim = errors.New("exterioralgebra: incompatible or invalid dimension")

// ErrIndex is returned when a basis index is negative, is not less than the
// ambient dimension, or when a blade specification repeats an index.
var ErrIndex = errors.New("exterioralgebra: basis index out of range or repeated")

// ErrGrade is returned when a grade argument is negative or exceeds the ambient
// dimension, or when a homogeneous operation receives a mixed-grade Form.
var ErrGrade = errors.New("exterioralgebra: invalid or incompatible grade")

// ErrMap is returned when a polynomial map has the wrong number of components
// for the requested pullback, or when its components disagree on arity.
var ErrMap = errors.New("exterioralgebra: malformed polynomial map")
