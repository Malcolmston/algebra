package lattices

import (
	"math/big"
	"strings"
)

// IntMatrix is a dense matrix of exact integers (*big.Int) stored row-major. It
// is used for exact integer lattice computations such as the Hermite normal
// form and integer determinants.
type IntMatrix struct {
	rows, cols int
	data       [][]*big.Int
}

// NewIntMatrix builds an integer matrix from int64 rows. Each row must have
// equal length. It panics on a ragged input.
func NewIntMatrix(rows [][]int64) IntMatrix {
	r := len(rows)
	c := 0
	if r > 0 {
		c = len(rows[0])
	}
	d := make([][]*big.Int, r)
	for i := range rows {
		if len(rows[i]) != c {
			panic("lattices: ragged matrix")
		}
		d[i] = make([]*big.Int, c)
		for j := range rows[i] {
			d[i][j] = big.NewInt(rows[i][j])
		}
	}
	return IntMatrix{rows: r, cols: c, data: d}
}

// ZeroIntMatrix returns an r-by-c integer matrix of zeros.
func ZeroIntMatrix(r, c int) IntMatrix {
	d := make([][]*big.Int, r)
	for i := range d {
		d[i] = make([]*big.Int, c)
		for j := range d[i] {
			d[i][j] = new(big.Int)
		}
	}
	return IntMatrix{rows: r, cols: c, data: d}
}

// IdentityIntMatrix returns the n-by-n integer identity matrix.
func IdentityIntMatrix(n int) IntMatrix {
	m := ZeroIntMatrix(n, n)
	for i := 0; i < n; i++ {
		m.data[i][i].SetInt64(1)
	}
	return m
}

// IntMatrixFromBasis converts a basis with integer entries to an IntMatrix. It
// returns ErrBadParameter if any entry is not (within rounding) an integer.
func IntMatrixFromBasis(b Basis) (IntMatrix, error) {
	m := ZeroIntMatrix(len(b), b.Dim())
	for i := range b {
		for j, x := range b[i] {
			r := new(big.Rat).SetFloat64(x)
			if r == nil || !r.IsInt() {
				return IntMatrix{}, ErrBadParameter
			}
			m.data[i][j].Set(r.Num())
		}
	}
	return m, nil
}

// Rows returns the number of rows of m.
func (m IntMatrix) Rows() int { return m.rows }

// Cols returns the number of columns of m.
func (m IntMatrix) Cols() int { return m.cols }

// At returns a copy of the entry in row i and column j.
func (m IntMatrix) At(i, j int) *big.Int { return new(big.Int).Set(m.data[i][j]) }

// Set stores a copy of x in row i and column j.
func (m IntMatrix) Set(i, j int, x *big.Int) { m.data[i][j].Set(x) }

// SetInt64 stores the int64 value x in row i and column j.
func (m IntMatrix) SetInt64(i, j int, x int64) { m.data[i][j].SetInt64(x) }

// Clone returns an independent deep copy of m.
func (m IntMatrix) Clone() IntMatrix {
	d := make([][]*big.Int, m.rows)
	for i := range m.data {
		d[i] = make([]*big.Int, m.cols)
		for j := range m.data[i] {
			d[i][j] = new(big.Int).Set(m.data[i][j])
		}
	}
	return IntMatrix{rows: m.rows, cols: m.cols, data: d}
}

// IsSquare reports whether m has equal numbers of rows and columns.
func (m IntMatrix) IsSquare() bool { return m.rows == m.cols }

// Transpose returns the transpose of m.
func (m IntMatrix) Transpose() IntMatrix {
	t := ZeroIntMatrix(m.cols, m.rows)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			t.data[j][i].Set(m.data[i][j])
		}
	}
	return t
}

// Mul returns the exact product m*n. It panics if the inner dimensions differ.
func (m IntMatrix) Mul(n IntMatrix) IntMatrix {
	if m.cols != n.rows {
		panic(ErrDimMismatch)
	}
	r := ZeroIntMatrix(m.rows, n.cols)
	tmp := new(big.Int)
	for i := 0; i < m.rows; i++ {
		for k := 0; k < m.cols; k++ {
			a := m.data[i][k]
			if a.Sign() == 0 {
				continue
			}
			for j := 0; j < n.cols; j++ {
				tmp.Mul(a, n.data[k][j])
				r.data[i][j].Add(r.data[i][j], tmp)
			}
		}
	}
	return r
}

// Det returns the exact integer determinant of a square matrix using the
// fraction-free Bareiss algorithm. It returns ErrNotSquare for a rectangular
// matrix.
func (m IntMatrix) Det() (*big.Int, error) {
	if !m.IsSquare() {
		return nil, ErrNotSquare
	}
	n := m.rows
	a := m.Clone().data
	sign := int64(1)
	prev := big.NewInt(1)
	num := new(big.Int)
	term := new(big.Int)
	for k := 0; k < n-1; k++ {
		if a[k][k].Sign() == 0 {
			swapped := false
			for p := k + 1; p < n; p++ {
				if a[p][k].Sign() != 0 {
					a[k], a[p] = a[p], a[k]
					sign = -sign
					swapped = true
					break
				}
			}
			if !swapped {
				return new(big.Int), nil
			}
		}
		for i := k + 1; i < n; i++ {
			for j := k + 1; j < n; j++ {
				num.Mul(a[i][j], a[k][k])
				term.Mul(a[i][k], a[k][j])
				num.Sub(num, term)
				a[i][j].Quo(num, prev)
			}
		}
		prev.Set(a[k][k])
	}
	det := new(big.Int).Set(a[n-1][n-1])
	if sign < 0 {
		det.Neg(det)
	}
	return det, nil
}

// Rat returns the exact rational version of m.
func (m IntMatrix) Rat() RatMatrix {
	r := ZeroRatMatrix(m.rows, m.cols)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			r.data[i][j].SetInt(m.data[i][j])
		}
	}
	return r
}

// Equal reports whether m and n have the same shape and equal entries.
func (m IntMatrix) Equal(n IntMatrix) bool {
	if m.rows != n.rows || m.cols != n.cols {
		return false
	}
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			if m.data[i][j].Cmp(n.data[i][j]) != 0 {
				return false
			}
		}
	}
	return true
}

// String renders m with one row per line.
func (m IntMatrix) String() string {
	var sb strings.Builder
	for i := 0; i < m.rows; i++ {
		parts := make([]string, m.cols)
		for j := 0; j < m.cols; j++ {
			parts[j] = m.data[i][j].String()
		}
		sb.WriteString("[" + strings.Join(parts, " ") + "]")
		if i < m.rows-1 {
			sb.WriteByte('\n')
		}
	}
	return sb.String()
}
