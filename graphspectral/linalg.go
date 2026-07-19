package graphspectral

import "math"

// Determinant returns the determinant of a square matrix, computed by Gaussian
// elimination with partial pivoting. It returns ErrNotSquare for a non-square
// matrix.
func Determinant(m *Matrix) (float64, error) {
	if !m.IsSquare() {
		return 0, ErrNotSquare
	}
	n := m.rows
	a := m.ToRows()
	det := 1.0
	for col := 0; col < n; col++ {
		// partial pivot
		piv := col
		max := math.Abs(a[col][col])
		for r := col + 1; r < n; r++ {
			if v := math.Abs(a[r][col]); v > max {
				max = v
				piv = r
			}
		}
		if max == 0 {
			return 0, nil
		}
		if piv != col {
			a[col], a[piv] = a[piv], a[col]
			det = -det
		}
		det *= a[col][col]
		for r := col + 1; r < n; r++ {
			f := a[r][col] / a[col][col]
			if f == 0 {
				continue
			}
			for c := col; c < n; c++ {
				a[r][c] -= f * a[col][c]
			}
		}
	}
	return det, nil
}

// SolveLinear solves the system A·x = b for a square matrix A, using Gaussian
// elimination with partial pivoting. It returns ErrNotSquare, ErrDimensionMismatch
// or ErrSingular as appropriate.
func SolveLinear(a *Matrix, b []float64) ([]float64, error) {
	if !a.IsSquare() {
		return nil, ErrNotSquare
	}
	n := a.rows
	if len(b) != n {
		return nil, ErrDimensionMismatch
	}
	m := a.ToRows()
	x := make([]float64, n)
	copy(x, b)
	for col := 0; col < n; col++ {
		piv := col
		max := math.Abs(m[col][col])
		for r := col + 1; r < n; r++ {
			if v := math.Abs(m[r][col]); v > max {
				max = v
				piv = r
			}
		}
		if max == 0 {
			return nil, ErrSingular
		}
		if piv != col {
			m[col], m[piv] = m[piv], m[col]
			x[col], x[piv] = x[piv], x[col]
		}
		for r := col + 1; r < n; r++ {
			f := m[r][col] / m[col][col]
			if f == 0 {
				continue
			}
			for c := col; c < n; c++ {
				m[r][c] -= f * m[col][c]
			}
			x[r] -= f * x[col]
		}
	}
	// back substitution
	for r := n - 1; r >= 0; r-- {
		s := x[r]
		for c := r + 1; c < n; c++ {
			s -= m[r][c] * x[c]
		}
		x[r] = s / m[r][r]
	}
	return x, nil
}

// Inverse returns the inverse of a square matrix using Gauss-Jordan elimination
// with partial pivoting. It returns ErrNotSquare or ErrSingular as appropriate.
func Inverse(a *Matrix) (*Matrix, error) {
	if !a.IsSquare() {
		return nil, ErrNotSquare
	}
	n := a.rows
	m := a.ToRows()
	inv := make([][]float64, n)
	for i := range inv {
		inv[i] = make([]float64, n)
		inv[i][i] = 1
	}
	for col := 0; col < n; col++ {
		piv := col
		max := math.Abs(m[col][col])
		for r := col + 1; r < n; r++ {
			if v := math.Abs(m[r][col]); v > max {
				max = v
				piv = r
			}
		}
		if max == 0 {
			return nil, ErrSingular
		}
		if piv != col {
			m[col], m[piv] = m[piv], m[col]
			inv[col], inv[piv] = inv[piv], inv[col]
		}
		d := m[col][col]
		for c := 0; c < n; c++ {
			m[col][c] /= d
			inv[col][c] /= d
		}
		for r := 0; r < n; r++ {
			if r == col {
				continue
			}
			f := m[r][col]
			if f == 0 {
				continue
			}
			for c := 0; c < n; c++ {
				m[r][c] -= f * m[col][c]
				inv[r][c] -= f * inv[col][c]
			}
		}
	}
	out, err := NewMatrixFromRows(inv)
	if err != nil {
		return nil, err
	}
	return out, nil
}

// Rank returns the numerical rank of a matrix: the number of pivots found during
// Gaussian elimination whose magnitude exceeds tol.
func Rank(m *Matrix, tol float64) int {
	a := m.ToRows()
	rows, cols := m.rows, m.cols
	rank := 0
	pr := 0
	for c := 0; c < cols && pr < rows; c++ {
		piv := pr
		max := math.Abs(a[pr][c])
		for r := pr + 1; r < rows; r++ {
			if v := math.Abs(a[r][c]); v > max {
				max = v
				piv = r
			}
		}
		if max <= tol {
			continue
		}
		a[pr], a[piv] = a[piv], a[pr]
		for r := 0; r < rows; r++ {
			if r == pr {
				continue
			}
			f := a[r][c] / a[pr][c]
			for cc := c; cc < cols; cc++ {
				a[r][cc] -= f * a[pr][cc]
			}
		}
		pr++
		rank++
	}
	return rank
}
