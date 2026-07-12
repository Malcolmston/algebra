package matrix

import (
	"errors"
	"strings"

	"github.com/malcolmston/algebra"
)

// Sentinel errors returned by the package.
var (
	// ErrDimension reports that the shapes of the operands are incompatible.
	ErrDimension = errors.New("matrix: incompatible dimensions")
	// ErrNotSquare reports that an operation requires a square matrix.
	ErrNotSquare = errors.New("matrix: matrix is not square")
	// ErrSingular reports that a matrix has a zero determinant and cannot be
	// inverted.
	ErrSingular = errors.New("matrix: matrix is singular")
	// ErrInconsistent reports that a linear system has no solution.
	ErrInconsistent = errors.New("matrix: system is inconsistent (no solution)")
	// ErrUnderdetermined reports that a linear system has infinitely many
	// solutions.
	ErrUnderdetermined = errors.New("matrix: system is underdetermined (infinitely many solutions)")
	// ErrUnsupported reports that an operation is not implemented for the given
	// input (for example Eigenvalues of a matrix larger than 3×3).
	ErrUnsupported = errors.New("matrix: unsupported operation for this input")
)

// zero is the canonical zero expression used to fill empty entries.
func zero() algebra.Expr { return algebra.Int(0) }

// one is the canonical unit expression.
func one() algebra.Expr { return algebra.Int(1) }

// simp is a short alias for algebra.Simplify.
func simp(e algebra.Expr) algebra.Expr { return algebra.Simplify(e) }

// Matrix is a two-dimensional, rows×cols array of algebra.Expr entries. The
// zero value is not usable; construct matrices with [New], [FromInts],
// [FromExpr], [FromFloats], [Identity], [Zeros], [Ones] or [Diag].
type Matrix struct {
	rows, cols int
	data       [][]algebra.Expr
}

// New returns a rows×cols matrix with every entry initialized to 0. It panics
// if rows or cols is negative.
func New(rows, cols int) *Matrix {
	if rows < 0 || cols < 0 {
		panic("matrix: negative dimension")
	}
	data := make([][]algebra.Expr, rows)
	for i := range data {
		data[i] = make([]algebra.Expr, cols)
		for j := range data[i] {
			data[i][j] = zero()
		}
	}
	return &Matrix{rows: rows, cols: cols, data: data}
}

// FromInts builds a matrix from a rectangular slice of int64 values. It panics
// if the rows have differing lengths.
func FromInts(vals [][]int64) *Matrix {
	r := len(vals)
	c := 0
	if r > 0 {
		c = len(vals[0])
	}
	m := New(r, c)
	for i := 0; i < r; i++ {
		if len(vals[i]) != c {
			panic("matrix: ragged input to FromInts")
		}
		for j := 0; j < c; j++ {
			m.data[i][j] = algebra.Int(vals[i][j])
		}
	}
	return m
}

// FromFloats builds a matrix from a rectangular slice of float64 values, storing
// each as an inexact algebra.Flt literal. It panics on ragged input.
func FromFloats(vals [][]float64) *Matrix {
	r := len(vals)
	c := 0
	if r > 0 {
		c = len(vals[0])
	}
	m := New(r, c)
	for i := 0; i < r; i++ {
		if len(vals[i]) != c {
			panic("matrix: ragged input to FromFloats")
		}
		for j := 0; j < c; j++ {
			m.data[i][j] = algebra.Flt(vals[i][j])
		}
	}
	return m
}

// FromExpr builds a matrix from a rectangular slice of expressions. Each entry
// is simplified. A nil entry is treated as 0. It panics on ragged input.
func FromExpr(vals [][]algebra.Expr) *Matrix {
	r := len(vals)
	c := 0
	if r > 0 {
		c = len(vals[0])
	}
	m := New(r, c)
	for i := 0; i < r; i++ {
		if len(vals[i]) != c {
			panic("matrix: ragged input to FromExpr")
		}
		for j := 0; j < c; j++ {
			if vals[i][j] == nil {
				m.data[i][j] = zero()
			} else {
				m.data[i][j] = simp(vals[i][j])
			}
		}
	}
	return m
}

// Identity returns the n×n identity matrix. It panics if n is negative.
func Identity(n int) *Matrix {
	m := New(n, n)
	for i := 0; i < n; i++ {
		m.data[i][i] = one()
	}
	return m
}

// Zeros returns a rows×cols matrix of zeros. It is a synonym for [New].
func Zeros(rows, cols int) *Matrix { return New(rows, cols) }

// Ones returns a rows×cols matrix with every entry equal to 1.
func Ones(rows, cols int) *Matrix {
	m := New(rows, cols)
	for i := 0; i < rows; i++ {
		for j := 0; j < cols; j++ {
			m.data[i][j] = one()
		}
	}
	return m
}

// Diag returns the square diagonal matrix whose main diagonal holds the given
// expressions and whose off-diagonal entries are 0.
func Diag(diag ...algebra.Expr) *Matrix {
	n := len(diag)
	m := New(n, n)
	for i := 0; i < n; i++ {
		if diag[i] == nil {
			m.data[i][i] = zero()
		} else {
			m.data[i][i] = simp(diag[i])
		}
	}
	return m
}

// Rows returns the number of rows.
func (m *Matrix) Rows() int { return m.rows }

// Cols returns the number of columns.
func (m *Matrix) Cols() int { return m.cols }

// IsSquare reports whether the matrix has equal row and column counts.
func (m *Matrix) IsSquare() bool { return m.rows == m.cols }

// At returns the entry at row i, column j (both 0-based). It panics if the
// indices are out of range.
func (m *Matrix) At(i, j int) algebra.Expr {
	m.checkIndex(i, j)
	return m.data[i][j]
}

// Set stores the simplified value of e at row i, column j (both 0-based). A nil
// value stores 0. It panics if the indices are out of range.
func (m *Matrix) Set(i, j int, e algebra.Expr) {
	m.checkIndex(i, j)
	if e == nil {
		m.data[i][j] = zero()
		return
	}
	m.data[i][j] = simp(e)
}

func (m *Matrix) checkIndex(i, j int) {
	if i < 0 || i >= m.rows || j < 0 || j >= m.cols {
		panic("matrix: index out of range")
	}
}

// Row returns a copy of row i as a [Vector]. It panics if i is out of range.
func (m *Matrix) Row(i int) *Vector {
	if i < 0 || i >= m.rows {
		panic("matrix: row index out of range")
	}
	out := make([]algebra.Expr, m.cols)
	copy(out, m.data[i])
	return &Vector{data: out}
}

// Col returns a copy of column j as a [Vector]. It panics if j is out of range.
func (m *Matrix) Col(j int) *Vector {
	if j < 0 || j >= m.cols {
		panic("matrix: column index out of range")
	}
	out := make([]algebra.Expr, m.rows)
	for i := 0; i < m.rows; i++ {
		out[i] = m.data[i][j]
	}
	return &Vector{data: out}
}

// Clone returns a deep copy of the matrix. Because expressions are immutable,
// only the backing slices are duplicated.
func (m *Matrix) Clone() *Matrix {
	out := New(m.rows, m.cols)
	for i := 0; i < m.rows; i++ {
		copy(out.data[i], m.data[i])
	}
	return out
}

// Equal reports whether m and n have the same shape and structurally equal
// entries. Entries are compared after simplification, so mathematically equal
// expressions built through the public API compare equal.
func (m *Matrix) Equal(n *Matrix) bool {
	if n == nil || m.rows != n.rows || m.cols != n.cols {
		return false
	}
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			if !simp(m.data[i][j]).Equal(simp(n.data[i][j])) {
				return false
			}
		}
	}
	return true
}

// String renders the matrix as an aligned, bracketed grid. Each column is
// padded to the width of its widest entry so the output lines up.
func (m *Matrix) String() string {
	if m.rows == 0 || m.cols == 0 {
		return "[ ]"
	}
	cells := make([][]string, m.rows)
	widths := make([]int, m.cols)
	for i := 0; i < m.rows; i++ {
		cells[i] = make([]string, m.cols)
		for j := 0; j < m.cols; j++ {
			s := m.data[i][j].String()
			cells[i][j] = s
			if len(s) > widths[j] {
				widths[j] = len(s)
			}
		}
	}
	var b strings.Builder
	for i := 0; i < m.rows; i++ {
		b.WriteString("[ ")
		for j := 0; j < m.cols; j++ {
			pad := widths[j] - len(cells[i][j])
			b.WriteString(strings.Repeat(" ", pad))
			b.WriteString(cells[i][j])
			if j != m.cols-1 {
				b.WriteString("  ")
			}
		}
		b.WriteString(" ]")
		if i != m.rows-1 {
			b.WriteByte('\n')
		}
	}
	return b.String()
}
