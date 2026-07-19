package lattices

import (
	"fmt"
	"math"
	"strings"
)

// Matrix is a dense real matrix stored row-major with float64 entries.
type Matrix struct {
	rows, cols int
	data       [][]float64
}

// NewMatrix builds a rows-by-cols matrix from the given row slices. Each row
// must have length cols. It panics if the shape is inconsistent.
func NewMatrix(rows [][]float64) Matrix {
	r := len(rows)
	c := 0
	if r > 0 {
		c = len(rows[0])
	}
	d := make([][]float64, r)
	for i := range rows {
		if len(rows[i]) != c {
			panic("lattices: ragged matrix")
		}
		d[i] = make([]float64, c)
		copy(d[i], rows[i])
	}
	return Matrix{rows: r, cols: c, data: d}
}

// ZeroMatrix returns an r-by-c matrix of zeros.
func ZeroMatrix(r, c int) Matrix {
	d := make([][]float64, r)
	for i := range d {
		d[i] = make([]float64, c)
	}
	return Matrix{rows: r, cols: c, data: d}
}

// IdentityMatrix returns the n-by-n identity matrix.
func IdentityMatrix(n int) Matrix {
	m := ZeroMatrix(n, n)
	for i := 0; i < n; i++ {
		m.data[i][i] = 1
	}
	return m
}

// DiagMatrix returns the square matrix with the given diagonal entries.
func DiagMatrix(diag ...float64) Matrix {
	n := len(diag)
	m := ZeroMatrix(n, n)
	for i := 0; i < n; i++ {
		m.data[i][i] = diag[i]
	}
	return m
}

// Rows returns the number of rows of m.
func (m Matrix) Rows() int { return m.rows }

// Cols returns the number of columns of m.
func (m Matrix) Cols() int { return m.cols }

// At returns the entry in row i and column j.
func (m Matrix) At(i, j int) float64 { return m.data[i][j] }

// Set stores x in row i and column j.
func (m Matrix) Set(i, j int, x float64) { m.data[i][j] = x }

// Clone returns an independent copy of m.
func (m Matrix) Clone() Matrix {
	d := make([][]float64, m.rows)
	for i := range m.data {
		d[i] = make([]float64, m.cols)
		copy(d[i], m.data[i])
	}
	return Matrix{rows: m.rows, cols: m.cols, data: d}
}

// Row returns a copy of row i as a Vec.
func (m Matrix) Row(i int) Vec {
	v := make(Vec, m.cols)
	copy(v, m.data[i])
	return v
}

// Col returns a copy of column j as a Vec.
func (m Matrix) Col(j int) Vec {
	v := make(Vec, m.rows)
	for i := 0; i < m.rows; i++ {
		v[i] = m.data[i][j]
	}
	return v
}

// IsSquare reports whether m has equal numbers of rows and columns.
func (m Matrix) IsSquare() bool { return m.rows == m.cols }

// IsSymmetric reports whether m is square and equal to its transpose within the
// tolerance tol.
func (m Matrix) IsSymmetric(tol float64) bool {
	if !m.IsSquare() {
		return false
	}
	for i := 0; i < m.rows; i++ {
		for j := i + 1; j < m.cols; j++ {
			if math.Abs(m.data[i][j]-m.data[j][i]) > tol {
				return false
			}
		}
	}
	return true
}

// Transpose returns the transpose of m.
func (m Matrix) Transpose() Matrix {
	t := ZeroMatrix(m.cols, m.rows)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			t.data[j][i] = m.data[i][j]
		}
	}
	return t
}

// Add returns m+n. It panics if the shapes differ.
func (m Matrix) Add(n Matrix) Matrix {
	m.mustSameShape(n)
	r := ZeroMatrix(m.rows, m.cols)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			r.data[i][j] = m.data[i][j] + n.data[i][j]
		}
	}
	return r
}

// Sub returns m-n. It panics if the shapes differ.
func (m Matrix) Sub(n Matrix) Matrix {
	m.mustSameShape(n)
	r := ZeroMatrix(m.rows, m.cols)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			r.data[i][j] = m.data[i][j] - n.data[i][j]
		}
	}
	return r
}

// Scale returns the scalar multiple s*m.
func (m Matrix) Scale(s float64) Matrix {
	r := ZeroMatrix(m.rows, m.cols)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			r.data[i][j] = s * m.data[i][j]
		}
	}
	return r
}

// Mul returns the matrix product m*n. It panics if the inner dimensions do not
// agree.
func (m Matrix) Mul(n Matrix) Matrix {
	if m.cols != n.rows {
		panic(ErrDimMismatch)
	}
	r := ZeroMatrix(m.rows, n.cols)
	for i := 0; i < m.rows; i++ {
		for k := 0; k < m.cols; k++ {
			a := m.data[i][k]
			if a == 0 {
				continue
			}
			for j := 0; j < n.cols; j++ {
				r.data[i][j] += a * n.data[k][j]
			}
		}
	}
	return r
}

// MulVec returns the matrix-vector product m*v. It panics if the dimensions do
// not agree.
func (m Matrix) MulVec(v Vec) Vec {
	if m.cols != len(v) {
		panic(ErrDimMismatch)
	}
	r := make(Vec, m.rows)
	for i := 0; i < m.rows; i++ {
		var s float64
		for j := 0; j < m.cols; j++ {
			s += m.data[i][j] * v[j]
		}
		r[i] = s
	}
	return r
}

// Trace returns the sum of the diagonal entries of a square matrix. It returns
// ErrNotSquare for a rectangular matrix.
func (m Matrix) Trace() (float64, error) {
	if !m.IsSquare() {
		return 0, ErrNotSquare
	}
	var s float64
	for i := 0; i < m.rows; i++ {
		s += m.data[i][i]
	}
	return s, nil
}

// Det returns the determinant of a square matrix via LU decomposition with
// partial pivoting. It returns ErrNotSquare for a rectangular matrix.
func (m Matrix) Det() (float64, error) {
	if !m.IsSquare() {
		return 0, ErrNotSquare
	}
	n := m.rows
	a := m.Clone().data
	det := 1.0
	for k := 0; k < n; k++ {
		// partial pivot
		p := k
		max := math.Abs(a[k][k])
		for i := k + 1; i < n; i++ {
			if v := math.Abs(a[i][k]); v > max {
				max = v
				p = i
			}
		}
		if max == 0 {
			return 0, nil
		}
		if p != k {
			a[k], a[p] = a[p], a[k]
			det = -det
		}
		det *= a[k][k]
		for i := k + 1; i < n; i++ {
			f := a[i][k] / a[k][k]
			for j := k; j < n; j++ {
				a[i][j] -= f * a[k][j]
			}
		}
	}
	return det, nil
}

// Inverse returns the inverse of a square matrix by Gauss-Jordan elimination
// with partial pivoting. It returns ErrNotSquare for a rectangular matrix and
// ErrSingular when the matrix is (numerically) singular.
func (m Matrix) Inverse() (Matrix, error) {
	if !m.IsSquare() {
		return Matrix{}, ErrNotSquare
	}
	n := m.rows
	a := m.Clone().data
	inv := IdentityMatrix(n).data
	for k := 0; k < n; k++ {
		p := k
		max := math.Abs(a[k][k])
		for i := k + 1; i < n; i++ {
			if v := math.Abs(a[i][k]); v > max {
				max = v
				p = i
			}
		}
		if max < 1e-300 {
			return Matrix{}, ErrSingular
		}
		if p != k {
			a[k], a[p] = a[p], a[k]
			inv[k], inv[p] = inv[p], inv[k]
		}
		piv := a[k][k]
		for j := 0; j < n; j++ {
			a[k][j] /= piv
			inv[k][j] /= piv
		}
		for i := 0; i < n; i++ {
			if i == k {
				continue
			}
			f := a[i][k]
			if f == 0 {
				continue
			}
			for j := 0; j < n; j++ {
				a[i][j] -= f * a[k][j]
				inv[i][j] -= f * inv[k][j]
			}
		}
	}
	return Matrix{rows: n, cols: n, data: inv}, nil
}

// Solve returns the solution x of m*x = b for a square matrix m using Gaussian
// elimination with partial pivoting. It returns ErrNotSquare, ErrDimMismatch or
// ErrSingular as appropriate.
func (m Matrix) Solve(b Vec) (Vec, error) {
	if !m.IsSquare() {
		return nil, ErrNotSquare
	}
	n := m.rows
	if len(b) != n {
		return nil, ErrDimMismatch
	}
	a := m.Clone().data
	x := make(Vec, n)
	copy(x, b)
	for k := 0; k < n; k++ {
		p := k
		max := math.Abs(a[k][k])
		for i := k + 1; i < n; i++ {
			if v := math.Abs(a[i][k]); v > max {
				max = v
				p = i
			}
		}
		if max < 1e-300 {
			return nil, ErrSingular
		}
		if p != k {
			a[k], a[p] = a[p], a[k]
			x[k], x[p] = x[p], x[k]
		}
		for i := k + 1; i < n; i++ {
			f := a[i][k] / a[k][k]
			for j := k; j < n; j++ {
				a[i][j] -= f * a[k][j]
			}
			x[i] -= f * x[k]
		}
	}
	for i := n - 1; i >= 0; i-- {
		s := x[i]
		for j := i + 1; j < n; j++ {
			s -= a[i][j] * x[j]
		}
		x[i] = s / a[i][i]
	}
	return x, nil
}

// Rank returns the numerical rank of m computed by Gaussian elimination, using
// tol as the threshold below which a pivot is treated as zero.
func (m Matrix) Rank(tol float64) int {
	a := m.Clone().data
	rows, cols := m.rows, m.cols
	rank := 0
	for col := 0; col < cols && rank < rows; col++ {
		p := rank
		max := math.Abs(a[rank][col])
		for i := rank + 1; i < rows; i++ {
			if v := math.Abs(a[i][col]); v > max {
				max = v
				p = i
			}
		}
		if max <= tol {
			continue
		}
		a[rank], a[p] = a[p], a[rank]
		piv := a[rank][col]
		for i := 0; i < rows; i++ {
			if i == rank {
				continue
			}
			f := a[i][col] / piv
			for j := col; j < cols; j++ {
				a[i][j] -= f * a[rank][j]
			}
		}
		rank++
	}
	return rank
}

// Equal reports whether m and n have the same shape and entries.
func (m Matrix) Equal(n Matrix) bool {
	if m.rows != n.rows || m.cols != n.cols {
		return false
	}
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			if m.data[i][j] != n.data[i][j] {
				return false
			}
		}
	}
	return true
}

// ApproxEqual reports whether m and n have the same shape and every entry
// agrees within the absolute tolerance tol.
func (m Matrix) ApproxEqual(n Matrix, tol float64) bool {
	if m.rows != n.rows || m.cols != n.cols {
		return false
	}
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			if math.Abs(m.data[i][j]-n.data[i][j]) > tol {
				return false
			}
		}
	}
	return true
}

// String renders m with one row per line.
func (m Matrix) String() string {
	var sb strings.Builder
	for i := 0; i < m.rows; i++ {
		parts := make([]string, m.cols)
		for j := 0; j < m.cols; j++ {
			parts[j] = fmt.Sprintf("%g", m.data[i][j])
		}
		sb.WriteString("[" + strings.Join(parts, " ") + "]")
		if i < m.rows-1 {
			sb.WriteByte('\n')
		}
	}
	return sb.String()
}

func (m Matrix) mustSameShape(n Matrix) {
	if m.rows != n.rows || m.cols != n.cols {
		panic(ErrDimMismatch)
	}
}
