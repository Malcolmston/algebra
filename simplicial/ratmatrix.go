package simplicial

import (
	"math/big"
	"strings"
)

// RatMatrix is a dense matrix over the field of rationals Q, with entries stored
// exactly as *big.Rat. A nil entry is treated as zero on read but the
// constructors always populate every cell.
type RatMatrix struct {
	rows, cols int
	data       [][]*big.Rat
}

// NewRatMatrix returns a rows×cols zero matrix over Q.
func NewRatMatrix(rows, cols int) *RatMatrix {
	if rows < 0 || cols < 0 {
		panic("simplicial: negative matrix dimension")
	}
	d := make([][]*big.Rat, rows)
	for i := range d {
		d[i] = make([]*big.Rat, cols)
		for j := range d[i] {
			d[i][j] = new(big.Rat)
		}
	}
	return &RatMatrix{rows: rows, cols: cols, data: d}
}

// RatIdentity returns the n×n identity matrix over Q.
func RatIdentity(n int) *RatMatrix {
	m := NewRatMatrix(n, n)
	for i := 0; i < n; i++ {
		m.data[i][i].SetInt64(1)
	}
	return m
}

// Rows returns the number of rows.
func (m *RatMatrix) Rows() int { return m.rows }

// Cols returns the number of columns.
func (m *RatMatrix) Cols() int { return m.cols }

// At returns a copy of the entry in row i, column j.
func (m *RatMatrix) At(i, j int) *big.Rat { return new(big.Rat).Set(m.data[i][j]) }

// Set stores a copy of v in row i, column j.
func (m *RatMatrix) Set(i, j int, v *big.Rat) { m.data[i][j].Set(v) }

// SetInt stores the integer v in row i, column j.
func (m *RatMatrix) SetInt(i, j int, v int64) { m.data[i][j].SetInt64(v) }

// Clone returns an independent copy of the matrix.
func (m *RatMatrix) Clone() *RatMatrix {
	out := NewRatMatrix(m.rows, m.cols)
	for i := range m.data {
		for j := range m.data[i] {
			out.data[i][j].Set(m.data[i][j])
		}
	}
	return out
}

// Equal reports whether m and n have the same shape and equal entries.
func (m *RatMatrix) Equal(n *RatMatrix) bool {
	if m.rows != n.rows || m.cols != n.cols {
		return false
	}
	for i := range m.data {
		for j := range m.data[i] {
			if m.data[i][j].Cmp(n.data[i][j]) != 0 {
				return false
			}
		}
	}
	return true
}

// IsZero reports whether every entry is zero.
func (m *RatMatrix) IsZero() bool {
	for i := range m.data {
		for j := range m.data[i] {
			if m.data[i][j].Sign() != 0 {
				return false
			}
		}
	}
	return true
}

// Transpose returns the transpose of the matrix.
func (m *RatMatrix) Transpose() *RatMatrix {
	out := NewRatMatrix(m.cols, m.rows)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			out.data[j][i].Set(m.data[i][j])
		}
	}
	return out
}

// Add returns m+n. The shapes must match; otherwise nil and [ErrShape].
func (m *RatMatrix) Add(n *RatMatrix) (*RatMatrix, error) {
	if m.rows != n.rows || m.cols != n.cols {
		return nil, ErrShape
	}
	out := NewRatMatrix(m.rows, m.cols)
	for i := range m.data {
		for j := range m.data[i] {
			out.data[i][j].Add(m.data[i][j], n.data[i][j])
		}
	}
	return out, nil
}

// Sub returns m−n. The shapes must match; otherwise nil and [ErrShape].
func (m *RatMatrix) Sub(n *RatMatrix) (*RatMatrix, error) {
	if m.rows != n.rows || m.cols != n.cols {
		return nil, ErrShape
	}
	out := NewRatMatrix(m.rows, m.cols)
	for i := range m.data {
		for j := range m.data[i] {
			out.data[i][j].Sub(m.data[i][j], n.data[i][j])
		}
	}
	return out, nil
}

// Mul returns the product m·n. The inner dimensions must agree; otherwise nil
// and [ErrShape].
func (m *RatMatrix) Mul(n *RatMatrix) (*RatMatrix, error) {
	if m.cols != n.rows {
		return nil, ErrShape
	}
	out := NewRatMatrix(m.rows, n.cols)
	tmp := new(big.Rat)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < n.cols; j++ {
			acc := new(big.Rat)
			for k := 0; k < m.cols; k++ {
				tmp.Mul(m.data[i][k], n.data[k][j])
				acc.Add(acc, tmp)
			}
			out.data[i][j] = acc
		}
	}
	return out, nil
}

// Scale returns the matrix with every entry multiplied by s.
func (m *RatMatrix) Scale(s *big.Rat) *RatMatrix {
	out := NewRatMatrix(m.rows, m.cols)
	for i := range m.data {
		for j := range m.data[i] {
			out.data[i][j].Mul(m.data[i][j], s)
		}
	}
	return out
}

// rowEchelon returns a reduced copy together with the pivot columns.
func (m *RatMatrix) rowEchelon() (*RatMatrix, []int) {
	r := m.Clone()
	var pivots []int
	row := 0
	for col := 0; col < r.cols && row < r.rows; col++ {
		sel := -1
		for i := row; i < r.rows; i++ {
			if r.data[i][col].Sign() != 0 {
				sel = i
				break
			}
		}
		if sel < 0 {
			continue
		}
		r.data[row], r.data[sel] = r.data[sel], r.data[row]
		// normalize pivot row to leading 1
		inv := new(big.Rat).Inv(r.data[row][col])
		for j := col; j < r.cols; j++ {
			r.data[row][j].Mul(r.data[row][j], inv)
		}
		// eliminate the column in all other rows
		factor := new(big.Rat)
		term := new(big.Rat)
		for i := 0; i < r.rows; i++ {
			if i == row || r.data[i][col].Sign() == 0 {
				continue
			}
			factor.Set(r.data[i][col])
			for j := col; j < r.cols; j++ {
				term.Mul(factor, r.data[row][j])
				r.data[i][j].Sub(r.data[i][j], term)
			}
		}
		pivots = append(pivots, col)
		row++
	}
	return r, pivots
}

// RowEchelon returns the reduced row-echelon form of the matrix over Q.
func (m *RatMatrix) RowEchelon() *RatMatrix {
	r, _ := m.rowEchelon()
	return r
}

// Rank returns the rank of the matrix over Q.
func (m *RatMatrix) Rank() int {
	_, pivots := m.rowEchelon()
	return len(pivots)
}

// Nullity returns cols − rank, the dimension of the kernel over Q.
func (m *RatMatrix) Nullity() int { return m.cols - m.Rank() }

// Determinant returns the determinant of a square matrix over Q. It returns
// [ErrNotSquare] for a rectangular matrix.
func (m *RatMatrix) Determinant() (*big.Rat, error) {
	if m.rows != m.cols {
		return nil, ErrNotSquare
	}
	r := m.Clone()
	n := m.rows
	det := new(big.Rat).SetInt64(1)
	tmp := new(big.Rat)
	term := new(big.Rat)
	for col := 0; col < n; col++ {
		sel := -1
		for i := col; i < n; i++ {
			if r.data[i][col].Sign() != 0 {
				sel = i
				break
			}
		}
		if sel < 0 {
			return new(big.Rat), nil
		}
		if sel != col {
			r.data[col], r.data[sel] = r.data[sel], r.data[col]
			det.Neg(det)
		}
		det.Mul(det, r.data[col][col])
		inv := new(big.Rat).Inv(r.data[col][col])
		for i := col + 1; i < n; i++ {
			if r.data[i][col].Sign() == 0 {
				continue
			}
			tmp.Mul(r.data[i][col], inv)
			for j := col; j < n; j++ {
				term.Mul(tmp, r.data[col][j])
				r.data[i][j].Sub(r.data[i][j], term)
			}
		}
	}
	return det, nil
}

// KernelBasis returns a basis of the null space {x : m·x = 0} over Q. Each basis
// vector is a length-cols slice of *big.Rat.
func (m *RatMatrix) KernelBasis() [][]*big.Rat {
	r, pivots := m.rowEchelon()
	pivotRowOf := make(map[int]int, len(pivots))
	for row, col := range pivots {
		pivotRowOf[col] = row
	}
	var basis [][]*big.Rat
	for free := 0; free < m.cols; free++ {
		if _, ok := pivotRowOf[free]; ok {
			continue
		}
		vec := make([]*big.Rat, m.cols)
		for i := range vec {
			vec[i] = new(big.Rat)
		}
		vec[free].SetInt64(1)
		for col, row := range pivotRowOf {
			vec[col].Neg(r.data[row][free])
		}
		basis = append(basis, vec)
	}
	return basis
}

// String returns a multi-line rendering of the matrix.
func (m *RatMatrix) String() string {
	var b strings.Builder
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			if j > 0 {
				b.WriteByte(' ')
			}
			b.WriteString(m.data[i][j].RatString())
		}
		b.WriteByte('\n')
	}
	return b.String()
}
