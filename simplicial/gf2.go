package simplicial

import "strings"

// GF2Matrix is a dense matrix over the two-element field GF(2) = {0,1}, with
// arithmetic performed modulo 2. Entries are stored as bytes that are always 0
// or 1.
type GF2Matrix struct {
	rows, cols int
	data       [][]uint8
}

// NewGF2Matrix returns a rows×cols zero matrix over GF(2). It panics if either
// dimension is negative.
func NewGF2Matrix(rows, cols int) *GF2Matrix {
	if rows < 0 || cols < 0 {
		panic("simplicial: negative matrix dimension")
	}
	d := make([][]uint8, rows)
	for i := range d {
		d[i] = make([]uint8, cols)
	}
	return &GF2Matrix{rows: rows, cols: cols, data: d}
}

// GF2Identity returns the n×n identity matrix over GF(2).
func GF2Identity(n int) *GF2Matrix {
	m := NewGF2Matrix(n, n)
	for i := 0; i < n; i++ {
		m.data[i][i] = 1
	}
	return m
}

// Rows returns the number of rows.
func (m *GF2Matrix) Rows() int { return m.rows }

// Cols returns the number of columns.
func (m *GF2Matrix) Cols() int { return m.cols }

// At returns the entry in row i, column j as 0 or 1.
func (m *GF2Matrix) At(i, j int) int { return int(m.data[i][j]) }

// Set stores v (reduced modulo 2) in row i, column j.
func (m *GF2Matrix) Set(i, j, v int) {
	if v&1 == 0 {
		m.data[i][j] = 0
	} else {
		m.data[i][j] = 1
	}
}

// Clone returns an independent copy of the matrix.
func (m *GF2Matrix) Clone() *GF2Matrix {
	out := NewGF2Matrix(m.rows, m.cols)
	for i := range m.data {
		copy(out.data[i], m.data[i])
	}
	return out
}

// Equal reports whether m and n have the same shape and entries.
func (m *GF2Matrix) Equal(n *GF2Matrix) bool {
	if m.rows != n.rows || m.cols != n.cols {
		return false
	}
	for i := range m.data {
		for j := range m.data[i] {
			if m.data[i][j] != n.data[i][j] {
				return false
			}
		}
	}
	return true
}

// IsZero reports whether every entry is zero.
func (m *GF2Matrix) IsZero() bool {
	for i := range m.data {
		for j := range m.data[i] {
			if m.data[i][j] != 0 {
				return false
			}
		}
	}
	return true
}

// Row returns a copy of row i as a slice of 0/1 ints.
func (m *GF2Matrix) Row(i int) []int {
	out := make([]int, m.cols)
	for j := range out {
		out[j] = int(m.data[i][j])
	}
	return out
}

// Col returns a copy of column j as a slice of 0/1 ints.
func (m *GF2Matrix) Col(j int) []int {
	out := make([]int, m.rows)
	for i := range out {
		out[i] = int(m.data[i][j])
	}
	return out
}

// SwapRows exchanges rows i and j in place.
func (m *GF2Matrix) SwapRows(i, j int) { m.data[i], m.data[j] = m.data[j], m.data[i] }

// AddRow adds (XORs) row src into row dst in place.
func (m *GF2Matrix) AddRow(dst, src int) {
	for j := 0; j < m.cols; j++ {
		m.data[dst][j] ^= m.data[src][j]
	}
}

// Transpose returns the transpose of the matrix.
func (m *GF2Matrix) Transpose() *GF2Matrix {
	out := NewGF2Matrix(m.cols, m.rows)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			out.data[j][i] = m.data[i][j]
		}
	}
	return out
}

// Add returns the entrywise sum (XOR) m+n over GF(2). The matrices must have the
// same shape; otherwise nil and [ErrShape] are returned.
func (m *GF2Matrix) Add(n *GF2Matrix) (*GF2Matrix, error) {
	if m.rows != n.rows || m.cols != n.cols {
		return nil, ErrShape
	}
	out := NewGF2Matrix(m.rows, m.cols)
	for i := range m.data {
		for j := range m.data[i] {
			out.data[i][j] = m.data[i][j] ^ n.data[i][j]
		}
	}
	return out, nil
}

// Mul returns the matrix product m·n over GF(2). The inner dimensions must
// agree; otherwise nil and [ErrShape] are returned.
func (m *GF2Matrix) Mul(n *GF2Matrix) (*GF2Matrix, error) {
	if m.cols != n.rows {
		return nil, ErrShape
	}
	out := NewGF2Matrix(m.rows, n.cols)
	for i := 0; i < m.rows; i++ {
		for k := 0; k < m.cols; k++ {
			if m.data[i][k] == 0 {
				continue
			}
			for j := 0; j < n.cols; j++ {
				out.data[i][j] ^= n.data[k][j]
			}
		}
	}
	return out, nil
}

// rowEchelon reduces a copy to row echelon form and returns the reduced matrix
// together with the list of pivot columns, one per pivot row.
func (m *GF2Matrix) rowEchelon() (*GF2Matrix, []int) {
	r := m.Clone()
	var pivots []int
	row := 0
	for col := 0; col < r.cols && row < r.rows; col++ {
		// find a pivot in this column at or below `row`
		sel := -1
		for i := row; i < r.rows; i++ {
			if r.data[i][col] != 0 {
				sel = i
				break
			}
		}
		if sel < 0 {
			continue
		}
		r.SwapRows(row, sel)
		for i := 0; i < r.rows; i++ {
			if i != row && r.data[i][col] != 0 {
				r.AddRow(i, row)
			}
		}
		pivots = append(pivots, col)
		row++
	}
	return r, pivots
}

// RowEchelon returns the reduced row-echelon form of the matrix over GF(2).
func (m *GF2Matrix) RowEchelon() *GF2Matrix {
	r, _ := m.rowEchelon()
	return r
}

// Rank returns the rank of the matrix over GF(2).
func (m *GF2Matrix) Rank() int {
	_, pivots := m.rowEchelon()
	return len(pivots)
}

// Nullity returns the dimension of the kernel of the matrix acting on column
// vectors, namely cols − rank.
func (m *GF2Matrix) Nullity() int { return m.cols - m.Rank() }

// KernelBasis returns a basis of the null space {x : m·x = 0} over GF(2). Each
// basis vector is returned as a length-cols slice of 0/1 ints; the number of
// vectors equals the nullity.
func (m *GF2Matrix) KernelBasis() [][]int {
	r, pivots := m.rowEchelon()
	pivotSet := make(map[int]int, len(pivots)) // col -> pivot row
	for row, col := range pivots {
		pivotSet[col] = row
	}
	var basis [][]int
	for free := 0; free < m.cols; free++ {
		if _, isPivot := pivotSet[free]; isPivot {
			continue
		}
		vec := make([]int, m.cols)
		vec[free] = 1
		for col, row := range pivotSet {
			if r.data[row][free] != 0 {
				vec[col] = 1
			}
		}
		basis = append(basis, vec)
	}
	return basis
}

// String returns a multi-line 0/1 rendering of the matrix.
func (m *GF2Matrix) String() string {
	var b strings.Builder
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			if j > 0 {
				b.WriteByte(' ')
			}
			if m.data[i][j] != 0 {
				b.WriteByte('1')
			} else {
				b.WriteByte('0')
			}
		}
		b.WriteByte('\n')
	}
	return b.String()
}
