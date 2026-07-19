package lattices

import (
	"math/big"
	"strings"
)

// RatMatrix is a dense matrix of exact rationals stored row-major.
type RatMatrix struct {
	rows, cols int
	data       [][]*big.Rat
}

// NewRatMatrix builds a rows-by-cols rational matrix from integer rows. Each
// row must have equal length. It panics on a ragged input.
func NewRatMatrix(rows [][]int64) RatMatrix {
	r := len(rows)
	c := 0
	if r > 0 {
		c = len(rows[0])
	}
	d := make([][]*big.Rat, r)
	for i := range rows {
		if len(rows[i]) != c {
			panic("lattices: ragged matrix")
		}
		d[i] = make([]*big.Rat, c)
		for j := range rows[i] {
			d[i][j] = new(big.Rat).SetInt64(rows[i][j])
		}
	}
	return RatMatrix{rows: r, cols: c, data: d}
}

// ZeroRatMatrix returns an r-by-c rational matrix of zeros.
func ZeroRatMatrix(r, c int) RatMatrix {
	d := make([][]*big.Rat, r)
	for i := range d {
		d[i] = make([]*big.Rat, c)
		for j := range d[i] {
			d[i][j] = new(big.Rat)
		}
	}
	return RatMatrix{rows: r, cols: c, data: d}
}

// IdentityRatMatrix returns the n-by-n rational identity matrix.
func IdentityRatMatrix(n int) RatMatrix {
	m := ZeroRatMatrix(n, n)
	for i := 0; i < n; i++ {
		m.data[i][i].SetInt64(1)
	}
	return m
}

// RatMatrixFromMatrix converts a float64 Matrix to an exact RatMatrix.
func RatMatrixFromMatrix(m Matrix) RatMatrix {
	r := ZeroRatMatrix(m.rows, m.cols)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			r.data[i][j] = ratFromFloat(m.data[i][j])
		}
	}
	return r
}

// Rows returns the number of rows of m.
func (m RatMatrix) Rows() int { return m.rows }

// Cols returns the number of columns of m.
func (m RatMatrix) Cols() int { return m.cols }

// At returns a copy of the entry in row i and column j.
func (m RatMatrix) At(i, j int) *big.Rat { return new(big.Rat).Set(m.data[i][j]) }

// Set stores a copy of x in row i and column j.
func (m RatMatrix) Set(i, j int, x *big.Rat) { m.data[i][j].Set(x) }

// Clone returns an independent deep copy of m.
func (m RatMatrix) Clone() RatMatrix {
	d := make([][]*big.Rat, m.rows)
	for i := range m.data {
		d[i] = make([]*big.Rat, m.cols)
		for j := range m.data[i] {
			d[i][j] = new(big.Rat).Set(m.data[i][j])
		}
	}
	return RatMatrix{rows: m.rows, cols: m.cols, data: d}
}

// IsSquare reports whether m has equal numbers of rows and columns.
func (m RatMatrix) IsSquare() bool { return m.rows == m.cols }

// Transpose returns the transpose of m.
func (m RatMatrix) Transpose() RatMatrix {
	t := ZeroRatMatrix(m.cols, m.rows)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			t.data[j][i].Set(m.data[i][j])
		}
	}
	return t
}

// Mul returns the exact product m*n. It panics if the inner dimensions differ.
func (m RatMatrix) Mul(n RatMatrix) RatMatrix {
	if m.cols != n.rows {
		panic(ErrDimMismatch)
	}
	r := ZeroRatMatrix(m.rows, n.cols)
	tmp := new(big.Rat)
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

// Add returns the exact sum m+n. It panics if the shapes differ.
func (m RatMatrix) Add(n RatMatrix) RatMatrix {
	if m.rows != n.rows || m.cols != n.cols {
		panic(ErrDimMismatch)
	}
	r := ZeroRatMatrix(m.rows, m.cols)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			r.data[i][j].Add(m.data[i][j], n.data[i][j])
		}
	}
	return r
}

// Sub returns the exact difference m-n. It panics if the shapes differ.
func (m RatMatrix) Sub(n RatMatrix) RatMatrix {
	if m.rows != n.rows || m.cols != n.cols {
		panic(ErrDimMismatch)
	}
	r := ZeroRatMatrix(m.rows, m.cols)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			r.data[i][j].Sub(m.data[i][j], n.data[i][j])
		}
	}
	return r
}

// Det returns the exact determinant of a square rational matrix by fraction
// (exact) Gaussian elimination. It returns ErrNotSquare for a rectangular
// matrix.
func (m RatMatrix) Det() (*big.Rat, error) {
	if !m.IsSquare() {
		return nil, ErrNotSquare
	}
	n := m.rows
	a := m.Clone().data
	det := new(big.Rat).SetInt64(1)
	tmp := new(big.Rat)
	for k := 0; k < n; k++ {
		if a[k][k].Sign() == 0 {
			swapped := false
			for i := k + 1; i < n; i++ {
				if a[i][k].Sign() != 0 {
					a[k], a[i] = a[i], a[k]
					det.Neg(det)
					swapped = true
					break
				}
			}
			if !swapped {
				return new(big.Rat), nil
			}
		}
		det.Mul(det, a[k][k])
		for i := k + 1; i < n; i++ {
			if a[i][k].Sign() == 0 {
				continue
			}
			f := new(big.Rat).Quo(a[i][k], a[k][k])
			for j := k; j < n; j++ {
				tmp.Mul(f, a[k][j])
				a[i][j].Sub(a[i][j], tmp)
			}
		}
	}
	return det, nil
}

// Inverse returns the exact inverse of a square rational matrix by Gauss-Jordan
// elimination. It returns ErrNotSquare or ErrSingular as appropriate.
func (m RatMatrix) Inverse() (RatMatrix, error) {
	if !m.IsSquare() {
		return RatMatrix{}, ErrNotSquare
	}
	n := m.rows
	a := m.Clone().data
	inv := IdentityRatMatrix(n).data
	tmp := new(big.Rat)
	for k := 0; k < n; k++ {
		if a[k][k].Sign() == 0 {
			swapped := false
			for i := k + 1; i < n; i++ {
				if a[i][k].Sign() != 0 {
					a[k], a[i] = a[i], a[k]
					inv[k], inv[i] = inv[i], inv[k]
					swapped = true
					break
				}
			}
			if !swapped {
				return RatMatrix{}, ErrSingular
			}
		}
		piv := new(big.Rat).Set(a[k][k])
		for j := 0; j < n; j++ {
			a[k][j].Quo(a[k][j], piv)
			inv[k][j].Quo(inv[k][j], piv)
		}
		for i := 0; i < n; i++ {
			if i == k || a[i][k].Sign() == 0 {
				continue
			}
			f := new(big.Rat).Set(a[i][k])
			for j := 0; j < n; j++ {
				tmp.Mul(f, a[k][j])
				a[i][j].Sub(a[i][j], tmp)
				tmp.Mul(f, inv[k][j])
				inv[i][j].Sub(inv[i][j], tmp)
			}
		}
	}
	return RatMatrix{rows: n, cols: n, data: inv}, nil
}

// Float returns the float64 approximation of m as a Matrix.
func (m RatMatrix) Float() Matrix {
	r := ZeroMatrix(m.rows, m.cols)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			f, _ := m.data[i][j].Float64()
			r.data[i][j] = f
		}
	}
	return r
}

// Equal reports whether m and n have the same shape and exactly equal entries.
func (m RatMatrix) Equal(n RatMatrix) bool {
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

// String renders m with one row per line using rational entries.
func (m RatMatrix) String() string {
	var sb strings.Builder
	for i := 0; i < m.rows; i++ {
		parts := make([]string, m.cols)
		for j := 0; j < m.cols; j++ {
			parts[j] = m.data[i][j].RatString()
		}
		sb.WriteString("[" + strings.Join(parts, " ") + "]")
		if i < m.rows-1 {
			sb.WriteByte('\n')
		}
	}
	return sb.String()
}
