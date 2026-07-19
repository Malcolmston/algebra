package liealgebra

import "math"

// luDecomp holds an LU factorisation with partial pivoting: the combined lower
// and upper factors, the pivot permutation and the sign of the permutation.
type luDecomp struct {
	n    int
	lu   []float64 // combined L (unit diagonal) and U, row-major
	piv  []int     // row permutation
	sign float64
}

// luFactor computes an LU decomposition with partial pivoting of a square
// matrix. It returns [ErrNotSquare] for non-square input and never fails on a
// singular matrix (the factorisation is still returned; detect singularity via
// a zero pivot in the U diagonal).
func luFactor(m *Matrix) (*luDecomp, error) {
	if !m.IsSquare() {
		return nil, ErrNotSquare
	}
	n := m.Rows
	lu := make([]float64, n*n)
	copy(lu, m.Data)
	piv := make([]int, n)
	for i := range piv {
		piv[i] = i
	}
	sign := 1.0
	for k := 0; k < n; k++ {
		// Find pivot.
		p := k
		max := math.Abs(lu[k*n+k])
		for i := k + 1; i < n; i++ {
			if v := math.Abs(lu[i*n+k]); v > max {
				max = v
				p = i
			}
		}
		if p != k {
			for j := 0; j < n; j++ {
				lu[k*n+j], lu[p*n+j] = lu[p*n+j], lu[k*n+j]
			}
			piv[k], piv[p] = piv[p], piv[k]
			sign = -sign
		}
		pivot := lu[k*n+k]
		if pivot == 0 {
			continue
		}
		for i := k + 1; i < n; i++ {
			lu[i*n+k] /= pivot
			f := lu[i*n+k]
			for j := k + 1; j < n; j++ {
				lu[i*n+j] -= f * lu[k*n+j]
			}
		}
	}
	return &luDecomp{n: n, lu: lu, piv: piv, sign: sign}, nil
}

// det returns the determinant implied by the factorisation.
func (d *luDecomp) det() float64 {
	det := d.sign
	for i := 0; i < d.n; i++ {
		det *= d.lu[i*d.n+i]
	}
	return det
}

// solveVec solves A x = b using the factorisation.
func (d *luDecomp) solveVec(b []float64) []float64 {
	n := d.n
	x := make([]float64, n)
	for i := 0; i < n; i++ {
		x[i] = b[d.piv[i]]
	}
	// Forward substitution (unit lower).
	for i := 0; i < n; i++ {
		s := x[i]
		for j := 0; j < i; j++ {
			s -= d.lu[i*n+j] * x[j]
		}
		x[i] = s
	}
	// Back substitution (upper).
	for i := n - 1; i >= 0; i-- {
		s := x[i]
		for j := i + 1; j < n; j++ {
			s -= d.lu[i*n+j] * x[j]
		}
		x[i] = s / d.lu[i*n+i]
	}
	return x
}

// Det returns the determinant of a square matrix.
func Det(m *Matrix) (float64, error) {
	d, err := luFactor(m)
	if err != nil {
		return 0, err
	}
	return d.det(), nil
}

// Inverse returns the matrix inverse, or [ErrSingular] if the matrix is not
// invertible.
func Inverse(m *Matrix) (*Matrix, error) {
	d, err := luFactor(m)
	if err != nil {
		return nil, err
	}
	if d.det() == 0 {
		return nil, ErrSingular
	}
	n := d.n
	inv := NewMatrix(n, n)
	e := make([]float64, n)
	for c := 0; c < n; c++ {
		for i := range e {
			e[i] = 0
		}
		e[c] = 1
		x := d.solveVec(e)
		for i := 0; i < n; i++ {
			inv.Data[i*n+c] = x[i]
		}
	}
	return inv, nil
}

// Solve solves the linear system A x = b for a square A. It returns
// [ErrSingular] if A is not invertible and [ErrDim] on a length mismatch.
func Solve(a *Matrix, b []float64) ([]float64, error) {
	if !a.IsSquare() {
		return nil, ErrNotSquare
	}
	if len(b) != a.Rows {
		return nil, ErrDim
	}
	d, err := luFactor(a)
	if err != nil {
		return nil, err
	}
	if d.det() == 0 {
		return nil, ErrSingular
	}
	return d.solveVec(b), nil
}

// SolveLeastSquares solves the (possibly overdetermined) system A x = b in the
// least-squares sense via the normal equations AᵀA x = Aᵀb. It requires A to
// have full column rank and returns [ErrRank] otherwise.
func SolveLeastSquares(a *Matrix, b []float64) ([]float64, error) {
	if len(b) != a.Rows {
		return nil, ErrDim
	}
	at := a.Transpose()
	ata, err := at.Mul(a)
	if err != nil {
		return nil, err
	}
	atb, err := at.MatVec(b)
	if err != nil {
		return nil, err
	}
	x, err := Solve(ata, atb)
	if err == ErrSingular {
		return nil, ErrRank
	}
	return x, err
}

// Rank returns the numerical rank of a matrix using Gaussian elimination with
// the given tolerance for treating a pivot as zero.
func Rank(m *Matrix, tol float64) int {
	a := m.Clone()
	n, cols := a.Rows, a.Cols
	rank := 0
	row := 0
	for col := 0; col < cols && row < n; col++ {
		// Find pivot in this column at or below row.
		p := row
		max := math.Abs(a.Data[row*cols+col])
		for i := row + 1; i < n; i++ {
			if v := math.Abs(a.Data[i*cols+col]); v > max {
				max = v
				p = i
			}
		}
		if max <= tol {
			continue
		}
		if p != row {
			for j := 0; j < cols; j++ {
				a.Data[row*cols+j], a.Data[p*cols+j] = a.Data[p*cols+j], a.Data[row*cols+j]
			}
		}
		pivot := a.Data[row*cols+col]
		for i := 0; i < n; i++ {
			if i == row {
				continue
			}
			f := a.Data[i*cols+col] / pivot
			for j := col; j < cols; j++ {
				a.Data[i*cols+j] -= f * a.Data[row*cols+j]
			}
		}
		row++
		rank++
	}
	return rank
}
