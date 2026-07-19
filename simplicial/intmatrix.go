package simplicial

import (
	"math/big"
	"strings"
)

// IntMatrix is a dense matrix over the integers Z, with entries stored exactly
// as *big.Int. It is the setting for the Smith normal form, from which the
// torsion of homology groups is read.
type IntMatrix struct {
	rows, cols int
	data       [][]*big.Int
}

// NewIntMatrix returns a rows×cols zero matrix over Z.
func NewIntMatrix(rows, cols int) *IntMatrix {
	if rows < 0 || cols < 0 {
		panic("simplicial: negative matrix dimension")
	}
	d := make([][]*big.Int, rows)
	for i := range d {
		d[i] = make([]*big.Int, cols)
		for j := range d[i] {
			d[i][j] = new(big.Int)
		}
	}
	return &IntMatrix{rows: rows, cols: cols, data: d}
}

// IntIdentity returns the n×n identity matrix over Z.
func IntIdentity(n int) *IntMatrix {
	m := NewIntMatrix(n, n)
	for i := 0; i < n; i++ {
		m.data[i][i].SetInt64(1)
	}
	return m
}

// IntMatrixFromInt64 builds an r×c integer matrix from a row-major slice of
// int64 values. It panics if len(vals) != r*c.
func IntMatrixFromInt64(r, c int, vals []int64) *IntMatrix {
	if len(vals) != r*c {
		panic("simplicial: value count does not match shape")
	}
	m := NewIntMatrix(r, c)
	for i := 0; i < r; i++ {
		for j := 0; j < c; j++ {
			m.data[i][j].SetInt64(vals[i*c+j])
		}
	}
	return m
}

// Rows returns the number of rows.
func (m *IntMatrix) Rows() int { return m.rows }

// Cols returns the number of columns.
func (m *IntMatrix) Cols() int { return m.cols }

// At returns a copy of the entry in row i, column j.
func (m *IntMatrix) At(i, j int) *big.Int { return new(big.Int).Set(m.data[i][j]) }

// Set stores a copy of v in row i, column j.
func (m *IntMatrix) Set(i, j int, v *big.Int) { m.data[i][j].Set(v) }

// SetInt stores the integer v in row i, column j.
func (m *IntMatrix) SetInt(i, j int, v int64) { m.data[i][j].SetInt64(v) }

// Clone returns an independent copy of the matrix.
func (m *IntMatrix) Clone() *IntMatrix {
	out := NewIntMatrix(m.rows, m.cols)
	for i := range m.data {
		for j := range m.data[i] {
			out.data[i][j].Set(m.data[i][j])
		}
	}
	return out
}

// Equal reports whether m and n have the same shape and entries.
func (m *IntMatrix) Equal(n *IntMatrix) bool {
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
func (m *IntMatrix) IsZero() bool {
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
func (m *IntMatrix) Transpose() *IntMatrix {
	out := NewIntMatrix(m.cols, m.rows)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			out.data[j][i].Set(m.data[i][j])
		}
	}
	return out
}

// Add returns m+n. The shapes must match; otherwise nil and [ErrShape].
func (m *IntMatrix) Add(n *IntMatrix) (*IntMatrix, error) {
	if m.rows != n.rows || m.cols != n.cols {
		return nil, ErrShape
	}
	out := NewIntMatrix(m.rows, m.cols)
	for i := range m.data {
		for j := range m.data[i] {
			out.data[i][j].Add(m.data[i][j], n.data[i][j])
		}
	}
	return out, nil
}

// Sub returns m−n. The shapes must match; otherwise nil and [ErrShape].
func (m *IntMatrix) Sub(n *IntMatrix) (*IntMatrix, error) {
	if m.rows != n.rows || m.cols != n.cols {
		return nil, ErrShape
	}
	out := NewIntMatrix(m.rows, m.cols)
	for i := range m.data {
		for j := range m.data[i] {
			out.data[i][j].Sub(m.data[i][j], n.data[i][j])
		}
	}
	return out, nil
}

// Mul returns the product m·n. The inner dimensions must agree; otherwise nil
// and [ErrShape].
func (m *IntMatrix) Mul(n *IntMatrix) (*IntMatrix, error) {
	if m.cols != n.rows {
		return nil, ErrShape
	}
	out := NewIntMatrix(m.rows, n.cols)
	tmp := new(big.Int)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < n.cols; j++ {
			acc := out.data[i][j]
			for k := 0; k < m.cols; k++ {
				tmp.Mul(m.data[i][k], n.data[k][j])
				acc.Add(acc, tmp)
			}
		}
	}
	return out, nil
}

// Neg returns −m.
func (m *IntMatrix) Neg() *IntMatrix {
	out := NewIntMatrix(m.rows, m.cols)
	for i := range m.data {
		for j := range m.data[i] {
			out.data[i][j].Neg(m.data[i][j])
		}
	}
	return out
}

// ToRat returns the matrix reinterpreted over the rationals.
func (m *IntMatrix) ToRat() *RatMatrix {
	out := NewRatMatrix(m.rows, m.cols)
	for i := range m.data {
		for j := range m.data[i] {
			out.data[i][j].SetInt(m.data[i][j])
		}
	}
	return out
}

// Rank returns the rank of the integer matrix, computed over Q (equivalently the
// number of non-zero invariant factors of its Smith normal form).
func (m *IntMatrix) Rank() int { return m.ToRat().Rank() }

// Determinant returns the determinant of a square integer matrix. It returns
// [ErrNotSquare] for a rectangular matrix.
func (m *IntMatrix) Determinant() (*big.Int, error) {
	if m.rows != m.cols {
		return nil, ErrNotSquare
	}
	d, err := m.ToRat().Determinant()
	if err != nil {
		return nil, err
	}
	// The determinant of an integer matrix is an integer.
	return new(big.Int).Set(d.Num()), nil
}

// --- Smith normal form -----------------------------------------------------

// SmithResult holds the Smith normal form decomposition U·A·V = D of an integer
// matrix A, where U and V are unimodular (invertible over Z) and D is diagonal
// with non-negative entries d_1 | d_2 | … | d_r each dividing the next.
type SmithResult struct {
	// D is the diagonal Smith normal form of the original matrix.
	D *IntMatrix
	// U is the unimodular left transform.
	U *IntMatrix
	// V is the unimodular right transform.
	V *IntMatrix
}

// Diagonal returns the diagonal entries of D, i.e. the invariant factors,
// padded with nothing beyond the smaller dimension. Trailing zeros indicate the
// rank deficiency.
func (s *SmithResult) Diagonal() []*big.Int {
	n := s.D.rows
	if s.D.cols < n {
		n = s.D.cols
	}
	out := make([]*big.Int, n)
	for i := 0; i < n; i++ {
		out[i] = new(big.Int).Set(s.D.data[i][i])
	}
	return out
}

// Rank returns the number of non-zero diagonal entries of D.
func (s *SmithResult) Rank() int {
	r := 0
	for _, d := range s.Diagonal() {
		if d.Sign() != 0 {
			r++
		}
	}
	return r
}

// InvariantFactors returns the non-zero diagonal entries of D, the invariant
// factors of the matrix.
func (s *SmithResult) InvariantFactors() []*big.Int {
	var out []*big.Int
	for _, d := range s.Diagonal() {
		if d.Sign() != 0 {
			out = append(out, d)
		}
	}
	return out
}

// TorsionFactors returns the invariant factors of D that are strictly greater
// than one; these are the elementary divisors giving rise to torsion.
func (s *SmithResult) TorsionFactors() []*big.Int {
	var out []*big.Int
	one := big.NewInt(1)
	for _, d := range s.InvariantFactors() {
		if d.Cmp(one) > 0 {
			out = append(out, new(big.Int).Set(d))
		}
	}
	return out
}

func swapRows(a *IntMatrix, i, j int) { a.data[i], a.data[j] = a.data[j], a.data[i] }

func swapCols(a *IntMatrix, i, j int) {
	for r := 0; r < a.rows; r++ {
		a.data[r][i], a.data[r][j] = a.data[r][j], a.data[r][i]
	}
}

func negateRow(a *IntMatrix, i int) {
	for j := 0; j < a.cols; j++ {
		a.data[i][j].Neg(a.data[i][j])
	}
}

func negateCol(a *IntMatrix, j int) {
	for i := 0; i < a.rows; i++ {
		a.data[i][j].Neg(a.data[i][j])
	}
}

// addRowMul performs row[dst] += q * row[src].
func addRowMul(a *IntMatrix, dst, src int, q *big.Int) {
	t := new(big.Int)
	for j := 0; j < a.cols; j++ {
		t.Mul(q, a.data[src][j])
		a.data[dst][j].Add(a.data[dst][j], t)
	}
}

// addColMul performs col[dst] += q * col[src].
func addColMul(a *IntMatrix, dst, src int, q *big.Int) {
	t := new(big.Int)
	for i := 0; i < a.rows; i++ {
		t.Mul(q, a.data[i][src])
		a.data[i][dst].Add(a.data[i][dst], t)
	}
}

// SmithNormalForm computes the Smith normal form of the matrix and returns U, D
// and V with U·A·V = D. The original matrix is not modified.
func (m *IntMatrix) SmithNormalForm() *SmithResult {
	a := m.Clone()
	u := IntIdentity(m.rows)
	v := IntIdentity(m.cols)

	minDim := m.rows
	if m.cols < minDim {
		minDim = m.cols
	}

	q := new(big.Int)
	for t := 0; t < minDim; t++ {
		for {
			// 1. find a non-zero pivot in submatrix [t:,t:] and bring to (t,t).
			pr, pc := -1, -1
			var best *big.Int
			for i := t; i < a.rows; i++ {
				for j := t; j < a.cols; j++ {
					if a.data[i][j].Sign() == 0 {
						continue
					}
					abs := new(big.Int).Abs(a.data[i][j])
					if best == nil || abs.Cmp(best) < 0 {
						best = abs
						pr, pc = i, j
					}
				}
			}
			if pr < 0 {
				break // submatrix is zero; done with this and later pivots
			}
			if pr != t {
				swapRows(a, t, pr)
				swapRows(u, t, pr)
			}
			if pc != t {
				swapCols(a, t, pc)
				swapCols(v, t, pc)
			}

			// 2. reduce column t below/above pivot.
			changed := false
			for i := 0; i < a.rows; i++ {
				if i == t || a.data[i][t].Sign() == 0 {
					continue
				}
				q.Quo(a.data[i][t], a.data[t][t])
				q.Neg(q)
				addRowMul(a, i, t, q)
				addRowMul(u, i, t, q)
				if a.data[i][t].Sign() != 0 {
					// remainder non-zero: promote and repeat
					swapRows(a, i, t)
					swapRows(u, i, t)
					changed = true
				}
			}
			if changed {
				continue
			}

			// 3. reduce row t left/right of pivot.
			for j := 0; j < a.cols; j++ {
				if j == t || a.data[t][j].Sign() == 0 {
					continue
				}
				q.Quo(a.data[t][j], a.data[t][t])
				q.Neg(q)
				addColMul(a, j, t, q)
				addColMul(v, j, t, q)
				if a.data[t][j].Sign() != 0 {
					swapCols(a, j, t)
					swapCols(v, j, t)
					changed = true
				}
			}
			if changed {
				continue
			}

			// 4. ensure a[t][t] divides every remaining entry; if not, fold the
			// offending row into row t and repeat.
			divisibilityOK := true
			for i := t + 1; i < a.rows && divisibilityOK; i++ {
				for j := t + 1; j < a.cols; j++ {
					r := new(big.Int)
					r.Rem(a.data[i][j], a.data[t][t])
					if r.Sign() != 0 {
						addRowMul(a, t, i, big.NewInt(1))
						addRowMul(u, t, i, big.NewInt(1))
						divisibilityOK = false
						break
					}
				}
			}
			if divisibilityOK {
				break
			}
		}
		// make the pivot non-negative
		if a.data[t][t].Sign() < 0 {
			negateRow(a, t)
			negateRow(u, t)
		}
	}

	return &SmithResult{D: a, U: u, V: v}
}

// SmithNormalFormDiagonal returns just the diagonal (invariant factors) of the
// Smith normal form.
func (m *IntMatrix) SmithNormalFormDiagonal() []*big.Int {
	return m.SmithNormalForm().Diagonal()
}

// InvariantFactors returns the non-zero invariant factors of the matrix.
func (m *IntMatrix) InvariantFactors() []*big.Int {
	return m.SmithNormalForm().InvariantFactors()
}

// String returns a multi-line rendering of the matrix.
func (m *IntMatrix) String() string {
	var b strings.Builder
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			if j > 0 {
				b.WriteByte(' ')
			}
			b.WriteString(m.data[i][j].String())
		}
		b.WriteByte('\n')
	}
	return b.String()
}
