package simplicial

import "errors"

// ErrEmptySimplex is returned when an operation requires a non-empty simplex
// but receives one with no vertices.
var ErrEmptySimplex = errors.New("simplicial: empty simplex")

// ErrDimension is returned when a dimension argument is negative or otherwise
// outside the range meaningful for the operand.
var ErrDimension = errors.New("simplicial: invalid dimension")

// ErrShape is returned when two matrices have incompatible shapes for the
// requested linear-algebra operation.
var ErrShape = errors.New("simplicial: incompatible matrix shape")

// ErrNotSquare is returned when a square matrix is required, for example by a
// determinant, but a rectangular one is supplied.
var ErrNotSquare = errors.New("simplicial: matrix is not square")

// ErrSingular is returned when an operation requires an invertible matrix but
// the operand is singular.
var ErrSingular = errors.New("simplicial: matrix is singular")

// ErrEmptyCloud is returned when a point-cloud operation requires at least one
// point but the cloud is empty.
var ErrEmptyCloud = errors.New("simplicial: empty point cloud")

// ErrDimMismatch is returned when two points, or a point and a metric, have
// mismatched ambient dimensions.
var ErrDimMismatch = errors.New("simplicial: point dimension mismatch")
