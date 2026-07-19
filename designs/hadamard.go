package designs

import "errors"

// HadamardMatrix is a square matrix with entries in {-1,+1} whose rows are
// pairwise orthogonal, i.e. H*H^T = n*I where n is the order.
type HadamardMatrix [][]int

// Order returns the side length of the matrix.
func (h HadamardMatrix) Order() int { return len(h) }

// IsHadamard reports whether the matrix has entries in {-1,+1} and satisfies
// H*H^T = n*I, the defining property of a Hadamard matrix.
func (h HadamardMatrix) IsHadamard() bool {
	n := len(h)
	if n == 0 {
		return false
	}
	for i := 0; i < n; i++ {
		if len(h[i]) != n {
			return false
		}
		for j := 0; j < n; j++ {
			if h[i][j] != 1 && h[i][j] != -1 {
				return false
			}
		}
	}
	for i := 0; i < n; i++ {
		for j := i; j < n; j++ {
			dot := 0
			for c := 0; c < n; c++ {
				dot += h[i][c] * h[j][c]
			}
			if i == j {
				if dot != n {
					return false
				}
			} else if dot != 0 {
				return false
			}
		}
	}
	return true
}

// Transpose returns the transpose of the matrix, which is again Hadamard.
func (h HadamardMatrix) Transpose() HadamardMatrix {
	n := len(h)
	t := make(HadamardMatrix, n)
	for i := 0; i < n; i++ {
		t[i] = make([]int, n)
		for j := 0; j < n; j++ {
			t[i][j] = h[j][i]
		}
	}
	return t
}

// Row returns a copy of row i.
func (h HadamardMatrix) Row(i int) []int { return append([]int(nil), h[i]...) }

// IsNormalized reports whether the first row and first column consist entirely
// of +1 entries.
func (h HadamardMatrix) IsNormalized() bool {
	n := len(h)
	for i := 0; i < n; i++ {
		if h[0][i] != 1 || h[i][0] != 1 {
			return false
		}
	}
	return true
}

// Normalize returns an equivalent Hadamard matrix whose first row and first
// column are all +1, obtained by negating the rows and columns that begin with
// -1.
func (h HadamardMatrix) Normalize() HadamardMatrix {
	n := len(h)
	out := make(HadamardMatrix, n)
	for i := 0; i < n; i++ {
		out[i] = append([]int(nil), h[i]...)
	}
	// Fix columns using the first row.
	for j := 0; j < n; j++ {
		if out[0][j] == -1 {
			for i := 0; i < n; i++ {
				out[i][j] = -out[i][j]
			}
		}
	}
	// Fix rows using the first column.
	for i := 0; i < n; i++ {
		if out[i][0] == -1 {
			for j := 0; j < n; j++ {
				out[i][j] = -out[i][j]
			}
		}
	}
	return out
}

// KroneckerProduct returns the Kronecker (tensor) product of two Hadamard
// matrices, itself a Hadamard matrix of order equal to the product of the two
// orders.
func KroneckerProduct(a, b HadamardMatrix) HadamardMatrix {
	na, nb := len(a), len(b)
	n := na * nb
	out := make(HadamardMatrix, n)
	for i := 0; i < n; i++ {
		out[i] = make([]int, n)
	}
	for ia := 0; ia < na; ia++ {
		for ja := 0; ja < na; ja++ {
			for ib := 0; ib < nb; ib++ {
				for jb := 0; jb < nb; jb++ {
					out[ia*nb+ib][ja*nb+jb] = a[ia][ja] * b[ib][jb]
				}
			}
		}
	}
	return out
}

// SylvesterHadamard returns the Sylvester Hadamard matrix of order 2**k, built
// recursively from H(1)=[1] by H(2m)=[[H,H],[H,-H]]. It reports an error when
// k<0.
func SylvesterHadamard(k int) (HadamardMatrix, error) {
	if k < 0 {
		return nil, errors.New("designs: exponent must be non-negative")
	}
	h := HadamardMatrix{{1}}
	base := HadamardMatrix{{1, 1}, {1, -1}}
	for i := 0; i < k; i++ {
		h = KroneckerProduct(h, base)
	}
	return h, nil
}

// QuadraticCharacter returns the quadratic character chi(x) of the element x in
// GF(q): 0 when x is zero, +1 when x is a non-zero square, and -1 when x is a
// non-square. It reports an error when q is not an odd prime power.
func QuadraticCharacter(q, x int) (int, error) {
	f, err := NewGaloisField(q)
	if err != nil {
		return 0, err
	}
	if f.Characteristic() == 2 {
		return 0, errors.New("designs: quadratic character requires odd characteristic")
	}
	x = ((x % q) + q) % q
	if x == 0 {
		return 0, nil
	}
	if f.Pow(x, (q-1)/2) == 1 {
		return 1, nil
	}
	return -1, nil
}

// JacobsthalMatrix returns the q-by-q Jacobsthal matrix Q of the field GF(q),
// with Q[i][j] = chi(i-j) where chi is the quadratic character. It reports an
// error when q is not an odd prime power.
func JacobsthalMatrix(q int) ([][]int, error) {
	f, err := NewGaloisField(q)
	if err != nil {
		return nil, err
	}
	if f.Characteristic() == 2 {
		return nil, errors.New("designs: Jacobsthal matrix requires odd characteristic")
	}
	chi := func(x int) int {
		if x == 0 {
			return 0
		}
		if f.Pow(x, (q-1)/2) == 1 {
			return 1
		}
		return -1
	}
	m := make([][]int, q)
	for i := 0; i < q; i++ {
		m[i] = make([]int, q)
		for j := 0; j < q; j++ {
			m[i][j] = chi(f.Sub(i, j))
		}
	}
	return m, nil
}

// PaleyConstructionI returns the Paley type-I Hadamard matrix of order q+1,
// defined for a prime power q congruent to 3 modulo 4. It reports an error when
// q is not such a prime power.
func PaleyConstructionI(q int) (HadamardMatrix, error) {
	if !IsPrimePower(q) || q%4 != 3 {
		return nil, errors.New("designs: Paley I requires a prime power q = 3 (mod 4)")
	}
	Q, err := JacobsthalMatrix(q)
	if err != nil {
		return nil, err
	}
	n := q + 1
	h := make(HadamardMatrix, n)
	for i := 0; i < n; i++ {
		h[i] = make([]int, n)
	}
	for j := 0; j < n; j++ {
		h[0][j] = 1
		h[j][0] = 1
	}
	for i := 0; i < q; i++ {
		for j := 0; j < q; j++ {
			if i == j {
				h[i+1][j+1] = -1
			} else {
				h[i+1][j+1] = Q[i][j]
			}
		}
	}
	return h, nil
}

// SymmetricConferenceMatrix returns the symmetric conference matrix of order
// q+1 built from GF(q) for a prime power q congruent to 1 modulo 4. It has a
// zero diagonal, entries +/-1 off the diagonal, and satisfies C*C^T = q*I. It
// reports an error when q is not such a prime power.
func SymmetricConferenceMatrix(q int) ([][]int, error) {
	if !IsPrimePower(q) || q%4 != 1 {
		return nil, errors.New("designs: symmetric conference matrix requires q = 1 (mod 4)")
	}
	Q, err := JacobsthalMatrix(q)
	if err != nil {
		return nil, err
	}
	n := q + 1
	c := make([][]int, n)
	for i := 0; i < n; i++ {
		c[i] = make([]int, n)
	}
	for j := 1; j < n; j++ {
		c[0][j] = 1
		c[j][0] = 1
	}
	for i := 0; i < q; i++ {
		for j := 0; j < q; j++ {
			c[i+1][j+1] = Q[i][j]
		}
	}
	return c, nil
}

// PaleyConstructionII returns the Paley type-II Hadamard matrix of order
// 2(q+1), defined for a prime power q congruent to 1 modulo 4. It reports an
// error when q is not such a prime power.
func PaleyConstructionII(q int) (HadamardMatrix, error) {
	c, err := SymmetricConferenceMatrix(q)
	if err != nil {
		return nil, err
	}
	m := len(c)
	n := 2 * m
	h := make(HadamardMatrix, n)
	for i := 0; i < n; i++ {
		h[i] = make([]int, n)
	}
	place := func(bi, bj int, blk [2][2]int) {
		h[2*bi][2*bj] = blk[0][0]
		h[2*bi][2*bj+1] = blk[0][1]
		h[2*bi+1][2*bj] = blk[1][0]
		h[2*bi+1][2*bj+1] = blk[1][1]
	}
	zero := [2][2]int{{1, -1}, {-1, -1}}
	plus := [2][2]int{{1, 1}, {1, -1}}
	minus := [2][2]int{{-1, -1}, {-1, 1}}
	for i := 0; i < m; i++ {
		for j := 0; j < m; j++ {
			switch {
			case c[i][j] == 0:
				place(i, j, zero)
			case c[i][j] == 1:
				place(i, j, plus)
			default:
				place(i, j, minus)
			}
		}
	}
	return h, nil
}

// HadamardOrderValid reports whether n is a possible Hadamard order: 1, 2, or a
// positive multiple of 4.
func HadamardOrderValid(n int) bool {
	return n == 1 || n == 2 || (n > 0 && n%4 == 0)
}

// HadamardToDesign converts a normalized Hadamard matrix of order 4t into the
// symmetric 2-(4t-1, 2t-1, t-1) design obtained by deleting the first row and
// column and reading each remaining row as a block indicator (+1 -> present).
// It reports an error when the matrix is not a Hadamard matrix of order
// divisible by 4.
func HadamardToDesign(h HadamardMatrix) (*Design, error) {
	n := len(h)
	if !h.IsHadamard() {
		return nil, errors.New("designs: not a Hadamard matrix")
	}
	if n%4 != 0 {
		return nil, errors.New("designs: order must be a multiple of 4")
	}
	hn := h.Normalize()
	v := n - 1
	blocks := make([][]int, 0, v)
	for i := 1; i < n; i++ {
		var b []int
		for j := 1; j < n; j++ {
			if hn[i][j] == 1 {
				b = append(b, j-1)
			}
		}
		blocks = append(blocks, b)
	}
	return NewDesign(v, blocks)
}
