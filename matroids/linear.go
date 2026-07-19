package matroids

import (
	"math/big"
)

// LinearMatroid is the column matroid of an exact rational matrix. The
// ground-set elements are the columns; a set of columns is independent exactly
// when the columns are linearly independent over the rationals. Ranks are
// computed by exact Gaussian elimination using math/big rationals, so results
// are free of floating-point error.
type LinearMatroid struct {
	rows int
	cols int
	data [][]*big.Rat // data[i][j] = entry in row i, column j
}

// NewLinearMatroid builds a linear matroid from a rows×cols matrix of int64
// entries given in row-major order. The ground set is the set of columns. It
// panics if any row length differs from the first row's length.
func NewLinearMatroid(matrix [][]int64) *LinearMatroid {
	rows := len(matrix)
	cols := 0
	if rows > 0 {
		cols = len(matrix[0])
	}
	data := make([][]*big.Rat, rows)
	for i := range matrix {
		if len(matrix[i]) != cols {
			panic("matroids: ragged matrix")
		}
		data[i] = make([]*big.Rat, cols)
		for j := range matrix[i] {
			data[i][j] = new(big.Rat).SetInt64(matrix[i][j])
		}
	}
	return &LinearMatroid{rows: rows, cols: cols, data: data}
}

// NewLinearMatroidRat builds a linear matroid from a rows×cols matrix of
// *big.Rat entries in row-major order. The matrix is copied. It panics on a
// ragged matrix.
func NewLinearMatroidRat(matrix [][]*big.Rat) *LinearMatroid {
	rows := len(matrix)
	cols := 0
	if rows > 0 {
		cols = len(matrix[0])
	}
	data := make([][]*big.Rat, rows)
	for i := range matrix {
		if len(matrix[i]) != cols {
			panic("matroids: ragged matrix")
		}
		data[i] = make([]*big.Rat, cols)
		for j := range matrix[i] {
			data[i][j] = new(big.Rat).Set(matrix[i][j])
		}
	}
	return &LinearMatroid{rows: rows, cols: cols, data: data}
}

// Size returns the number of columns (ground-set elements).
func (m *LinearMatroid) Size() int { return m.cols }

// NumRows returns the number of rows of the matrix.
func (m *LinearMatroid) NumRows() int { return m.rows }

// Entry returns a copy of the matrix entry in row i, column j.
func (m *LinearMatroid) Entry(i, j int) *big.Rat { return new(big.Rat).Set(m.data[i][j]) }

// Rank returns the rank of the submatrix formed by the columns in set, computed
// by exact Gaussian elimination.
func (m *LinearMatroid) Rank(set []int) int {
	cols := distinctInRange(set, m.cols)
	if len(cols) == 0 || m.rows == 0 {
		return 0
	}
	// Build a working copy: rows × len(cols).
	work := make([][]*big.Rat, m.rows)
	for i := 0; i < m.rows; i++ {
		work[i] = make([]*big.Rat, len(cols))
		for k, j := range cols {
			work[i][k] = new(big.Rat).Set(m.data[i][j])
		}
	}
	return ratMatrixRank(work, m.rows, len(cols))
}

// distinctInRange returns the distinct elements of set that fall in [0, n),
// sorted.
func distinctInRange(set []int, n int) []int {
	seen := make(map[int]bool, len(set))
	var out []int
	for _, e := range set {
		if e >= 0 && e < n && !seen[e] {
			seen[e] = true
			out = append(out, e)
		}
	}
	// keep sorted for determinism
	for i := 1; i < len(out); i++ {
		for j := i; j > 0 && out[j-1] > out[j]; j-- {
			out[j-1], out[j] = out[j], out[j-1]
		}
	}
	return out
}

// ratMatrixRank computes the rank of an r×c matrix of rationals by Gaussian
// elimination. The matrix is modified in place.
func ratMatrixRank(a [][]*big.Rat, r, c int) int {
	rank := 0
	for col := 0; col < c && rank < r; col++ {
		// find pivot
		pivot := -1
		for i := rank; i < r; i++ {
			if a[i][col].Sign() != 0 {
				pivot = i
				break
			}
		}
		if pivot == -1 {
			continue
		}
		a[rank], a[pivot] = a[pivot], a[rank]
		pivInv := new(big.Rat).Inv(a[rank][col])
		for j := col; j < c; j++ {
			a[rank][j].Mul(a[rank][j], pivInv)
		}
		for i := 0; i < r; i++ {
			if i == rank || a[i][col].Sign() == 0 {
				continue
			}
			factor := new(big.Rat).Set(a[i][col])
			for j := col; j < c; j++ {
				t := new(big.Rat).Mul(factor, a[rank][j])
				a[i][j].Sub(a[i][j], t)
			}
		}
		rank++
	}
	return rank
}

// BinaryMatroid is the column matroid of a 0/1 matrix over the two-element
// field GF(2). The ground-set elements are the columns; a set of columns is
// independent when they are linearly independent over GF(2). Columns are stored
// as big.Int bit vectors (bit i set means row i has a 1), so any number of rows
// is supported.
type BinaryMatroid struct {
	rows int
	cols int
	col  []*big.Int // col[j] holds the bits of column j
}

// NewBinaryMatroid builds a GF(2) linear matroid from a rows×cols 0/1 matrix in
// row-major order. Non-zero entries are treated as 1 modulo 2. It panics on a
// ragged matrix.
func NewBinaryMatroid(matrix [][]int) *BinaryMatroid {
	rows := len(matrix)
	cols := 0
	if rows > 0 {
		cols = len(matrix[0])
	}
	col := make([]*big.Int, cols)
	for j := 0; j < cols; j++ {
		col[j] = new(big.Int)
	}
	for i := 0; i < rows; i++ {
		if len(matrix[i]) != cols {
			panic("matroids: ragged matrix")
		}
		for j := 0; j < cols; j++ {
			if matrix[i][j]%2 != 0 {
				col[j].SetBit(col[j], i, 1)
			}
		}
	}
	return &BinaryMatroid{rows: rows, cols: cols, col: col}
}

// Size returns the number of columns (ground-set elements).
func (m *BinaryMatroid) Size() int { return m.cols }

// NumRows returns the number of rows of the matrix.
func (m *BinaryMatroid) NumRows() int { return m.rows }

// Rank returns the GF(2) rank of the columns in set, computed by Gaussian
// elimination over GF(2). Each accumulated basis vector is kept in reduced
// form with a distinct leading (highest set) bit.
func (m *BinaryMatroid) Rank(set []int) int {
	cols := distinctInRange(set, m.cols)
	// basis[k] is a reduced independent vector; leadBit[k] is its highest set
	// bit, used as the pivot.
	var basis []*big.Int
	var leadBit []int
	for _, j := range cols {
		v := new(big.Int).Set(m.col[j])
		for {
			if v.Sign() == 0 {
				break
			}
			lead := v.BitLen() - 1
			// reduce by any basis vector with the same leading bit
			reduced := false
			for k := range basis {
				if leadBit[k] == lead {
					v.Xor(v, basis[k])
					reduced = true
					break
				}
			}
			if !reduced {
				basis = append(basis, v)
				leadBit = append(leadBit, lead)
				break
			}
		}
	}
	return len(basis)
}
