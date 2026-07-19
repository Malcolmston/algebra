package galois

import "math/big"

// MatModP is a dense matrix with entries in the prime field GF(p). Entries are
// stored in Data[row][col], each reduced into [0, p).
type MatModP struct {
	P    *big.Int
	Rows int
	Cols int
	Data [][]*big.Int
}

// NewMatModP returns a rows×cols zero matrix over GF(p).
func NewMatModP(p *big.Int, rows, cols int) *MatModP {
	d := make([][]*big.Int, rows)
	for i := range d {
		d[i] = make([]*big.Int, cols)
		for j := range d[i] {
			d[i][j] = big.NewInt(0)
		}
	}
	return &MatModP{P: clone(p), Rows: rows, Cols: cols, Data: d}
}

// IdentityModP returns the n×n identity matrix over GF(p).
func IdentityModP(p *big.Int, n int) *MatModP {
	m := NewMatModP(p, n, n)
	for i := 0; i < n; i++ {
		m.Data[i][i] = big.NewInt(1)
	}
	return m
}

// Get returns the entry at (i, j).
func (m *MatModP) Get(i, j int) *big.Int { return clone(m.Data[i][j]) }

// Set stores v (reduced mod p) at (i, j).
func (m *MatModP) Set(i, j int, v *big.Int) { m.Data[i][j] = new(big.Int).Mod(v, m.P) }

// Clone returns an independent copy of the matrix.
func (m *MatModP) Clone() *MatModP {
	c := NewMatModP(m.P, m.Rows, m.Cols)
	for i := 0; i < m.Rows; i++ {
		for j := 0; j < m.Cols; j++ {
			c.Data[i][j] = clone(m.Data[i][j])
		}
	}
	return c
}

// Transpose returns the transpose of the matrix.
func (m *MatModP) Transpose() *MatModP {
	t := NewMatModP(m.P, m.Cols, m.Rows)
	for i := 0; i < m.Rows; i++ {
		for j := 0; j < m.Cols; j++ {
			t.Data[j][i] = clone(m.Data[i][j])
		}
	}
	return t
}

// Mul returns the matrix product m*other over GF(p). The inner dimensions must
// agree; a mismatch returns nil.
func (m *MatModP) Mul(other *MatModP) *MatModP {
	if m.Cols != other.Rows {
		return nil
	}
	r := NewMatModP(m.P, m.Rows, other.Cols)
	for i := 0; i < m.Rows; i++ {
		for j := 0; j < other.Cols; j++ {
			acc := big.NewInt(0)
			for k := 0; k < m.Cols; k++ {
				acc.Add(acc, new(big.Int).Mul(m.Data[i][k], other.Data[k][j]))
			}
			r.Data[i][j] = acc.Mod(acc, m.P)
		}
	}
	return r
}

// Sub returns the entrywise difference m-other over GF(p). A shape mismatch
// returns nil.
func (m *MatModP) Sub(other *MatModP) *MatModP {
	if m.Rows != other.Rows || m.Cols != other.Cols {
		return nil
	}
	r := NewMatModP(m.P, m.Rows, m.Cols)
	for i := 0; i < m.Rows; i++ {
		for j := 0; j < m.Cols; j++ {
			r.Data[i][j] = SubMod(m.Data[i][j], other.Data[i][j], m.P)
		}
	}
	return r
}

// RowReduce returns the reduced row-echelon form of the matrix together with
// its rank. The receiver is left unchanged.
func (m *MatModP) RowReduce() (*MatModP, int) {
	a := m.Clone()
	rank := 0
	for col := 0; col < a.Cols && rank < a.Rows; col++ {
		// find a pivot in this column at or below rank.
		pivot := -1
		for r := rank; r < a.Rows; r++ {
			if a.Data[r][col].Sign() != 0 {
				pivot = r
				break
			}
		}
		if pivot == -1 {
			continue
		}
		a.Data[rank], a.Data[pivot] = a.Data[pivot], a.Data[rank]
		inv, _ := InvMod(a.Data[rank][col], a.P)
		for j := 0; j < a.Cols; j++ {
			a.Data[rank][j] = MulMod(a.Data[rank][j], inv, a.P)
		}
		for r := 0; r < a.Rows; r++ {
			if r == rank || a.Data[r][col].Sign() == 0 {
				continue
			}
			factor := clone(a.Data[r][col])
			for j := 0; j < a.Cols; j++ {
				t := MulMod(factor, a.Data[rank][j], a.P)
				a.Data[r][j] = SubMod(a.Data[r][j], t, a.P)
			}
		}
		rank++
	}
	return a, rank
}

// Rank returns the rank of the matrix over GF(p).
func (m *MatModP) Rank() int {
	_, r := m.RowReduce()
	return r
}

// NullSpace returns a basis for the right null space {x : m·x = 0} over GF(p),
// each basis vector represented as a length-Cols slice.
func (m *MatModP) NullSpace() [][]*big.Int {
	rref, _ := m.RowReduce()
	// identify pivot columns.
	pivotCol := make([]int, 0, rref.Rows)
	isPivot := make([]bool, m.Cols)
	for r := 0; r < rref.Rows; r++ {
		lead := -1
		for c := 0; c < rref.Cols; c++ {
			if rref.Data[r][c].Sign() != 0 {
				lead = c
				break
			}
		}
		if lead >= 0 {
			pivotCol = append(pivotCol, lead)
			isPivot[lead] = true
		}
	}
	var basis [][]*big.Int
	for free := 0; free < m.Cols; free++ {
		if isPivot[free] {
			continue
		}
		vec := make([]*big.Int, m.Cols)
		for i := range vec {
			vec[i] = big.NewInt(0)
		}
		vec[free] = big.NewInt(1)
		for r, pc := range pivotCol {
			// pivot row r has the equation: x[pc] + sum_{free} coeff*x[free] = 0.
			vec[pc] = NegMod(rref.Data[r][free], m.P)
		}
		basis = append(basis, vec)
	}
	return basis
}
