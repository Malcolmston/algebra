package matrix

import "github.com/malcolmston/algebra"

// Solve returns the unique solution x of the linear system a·x = b, where a is
// an m×n coefficient matrix and b is a length-m vector.
//
// The augmented matrix [a|b] is reduced to reduced row-echelon form with exact
// symbolic pivots and the outcome is classified:
//
//   - a unique solution is returned when the coefficient matrix has full column
//     rank and the system is consistent;
//   - [ErrInconsistent] is returned when the system has no solution;
//   - [ErrUnderdetermined] is returned when the system has infinitely many
//     solutions (fewer pivots than unknowns);
//   - [ErrDimension] is returned when len(b) != a.Rows().
func Solve(a *Matrix, b *Vector) (*Vector, error) {
	if a.rows != len(b.data) {
		return nil, ErrDimension
	}
	n := a.cols
	// Build the augmented matrix [a | b].
	aug := New(a.rows, n+1)
	for i := 0; i < a.rows; i++ {
		for j := 0; j < n; j++ {
			aug.data[i][j] = a.data[i][j]
		}
		aug.data[i][n] = b.data[i]
	}
	rref, _ := aug.RREF()

	// Detect inconsistency: a row that is all zero across the coefficient
	// columns but has a nonzero augmented entry.
	for i := 0; i < rref.rows; i++ {
		allZero := true
		for j := 0; j < n; j++ {
			if !isZeroExpr(rref.data[i][j]) {
				allZero = false
				break
			}
		}
		if allZero && !isZeroExpr(rref.data[i][n]) {
			return nil, ErrInconsistent
		}
	}

	// Count pivots in the coefficient part.
	pivots := 0
	for i := 0; i < rref.rows; i++ {
		for j := 0; j < n; j++ {
			if !isZeroExpr(rref.data[i][j]) {
				pivots++
				break
			}
		}
	}
	if pivots < n {
		return nil, ErrUnderdetermined
	}

	// Full column rank and consistent: with n pivots in RREF, the pivot for
	// column i sits at row i, so x_i is the augmented entry there.
	out := make([]algebra.Expr, n)
	for i := 0; i < n; i++ {
		out[i] = rref.data[i][n]
	}
	return &Vector{data: out}, nil
}
