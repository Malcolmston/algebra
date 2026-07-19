package codingtheory

import "errors"

// ErrRank is returned when a generator matrix is not full row rank.
var ErrRank = errors.New("codingtheory: generator matrix is not full rank")

// ErrShape is returned when a matrix argument is empty or ragged.
var ErrShape = errors.New("codingtheory: empty or ragged matrix")

// LinearCode is a binary linear [n,k] block code over GF(2) described by a
// k-by-n generator matrix G and an (n-k)-by-n parity-check matrix H with
// G*Hᵀ = 0. Messages of length k are encoded as m*G; a received word r is a
// code word exactly when H*rᵀ = 0.
type LinearCode struct {
	n, k int
	g    [][]int // k x n generator
	h    [][]int // (n-k) x n parity check
}

// rrefGF2 returns the reduced row-echelon form of a copy of m over GF(2) along
// with the list of pivot column indices.
func rrefGF2(m [][]int) ([][]int, []int) {
	rows := len(m)
	cols := len(m[0])
	a := make([][]int, rows)
	for i := range m {
		a[i] = append([]int(nil), m[i]...)
	}
	var pivots []int
	r := 0
	for c := 0; c < cols && r < rows; c++ {
		pivot := -1
		for i := r; i < rows; i++ {
			if a[i][c] != 0 {
				pivot = i
				break
			}
		}
		if pivot == -1 {
			continue
		}
		a[r], a[pivot] = a[pivot], a[r]
		for i := 0; i < rows; i++ {
			if i != r && a[i][c] != 0 {
				for j := 0; j < cols; j++ {
					a[i][j] ^= a[r][j]
				}
			}
		}
		pivots = append(pivots, c)
		r++
	}
	return a, pivots
}

// nullSpaceGF2 returns a basis of the right null space {x : m*x = 0} of a GF(2)
// matrix m with n columns.
func nullSpaceGF2(m [][]int, n int) [][]int {
	if len(m) == 0 {
		basis := make([][]int, n)
		for i := 0; i < n; i++ {
			basis[i] = make([]int, n)
			basis[i][i] = 1
		}
		return basis
	}
	a, pivots := rrefGF2(m)
	isPivot := make([]bool, n)
	for _, p := range pivots {
		isPivot[p] = true
	}
	// Map each pivot column to its row.
	pivotRow := make([]int, n)
	for i, p := range pivots {
		pivotRow[p] = i
	}
	var basis [][]int
	for free := 0; free < n; free++ {
		if isPivot[free] {
			continue
		}
		vec := make([]int, n)
		vec[free] = 1
		for _, p := range pivots {
			// pivot variable equals sum of free variables in its row.
			if a[pivotRow[p]][free] != 0 {
				vec[p] = 1
			}
		}
		basis = append(basis, vec)
	}
	return basis
}

// NewLinearCode builds a linear code from a k-by-n generator matrix G. G must
// be non-empty, rectangular and full row rank; the parity-check matrix is
// derived as a basis of the null space of G. It returns ErrShape or ErrRank on
// invalid input.
func NewLinearCode(g [][]int) (*LinearCode, error) {
	if len(g) == 0 || len(g[0]) == 0 {
		return nil, ErrShape
	}
	n := len(g[0])
	for _, row := range g {
		if len(row) != n {
			return nil, ErrShape
		}
		if !validBits(row) {
			return nil, ErrBit
		}
	}
	k := len(g)
	if RankGF2(g) != k {
		return nil, ErrRank
	}
	gc := make([][]int, k)
	for i := range g {
		gc[i] = append([]int(nil), g[i]...)
	}
	h := nullSpaceGF2(g, n)
	return &LinearCode{n: n, k: k, g: gc, h: h}, nil
}

// NewLinearCodeGH builds a linear code from explicit generator and parity-check
// matrices. It verifies that shapes are consistent and G*Hᵀ = 0.
func NewLinearCodeGH(g, h [][]int) (*LinearCode, error) {
	if len(g) == 0 || len(g[0]) == 0 || len(h) == 0 {
		return nil, ErrShape
	}
	n := len(g[0])
	for _, row := range g {
		if len(row) != n || !validBits(row) {
			return nil, ErrShape
		}
	}
	for _, row := range h {
		if len(row) != n || !validBits(row) {
			return nil, ErrShape
		}
	}
	// Check orthogonality G*Hᵀ = 0.
	for _, gr := range g {
		for _, hr := range h {
			if DotGF2(gr, hr) != 0 {
				return nil, errors.New("codingtheory: G and H are not orthogonal")
			}
		}
	}
	gc := make([][]int, len(g))
	for i := range g {
		gc[i] = append([]int(nil), g[i]...)
	}
	hc := make([][]int, len(h))
	for i := range h {
		hc[i] = append([]int(nil), h[i]...)
	}
	return &LinearCode{n: n, k: len(g), g: gc, h: hc}, nil
}

// N returns the code length n.
func (c *LinearCode) N() int { return c.n }

// K returns the message dimension k.
func (c *LinearCode) K() int { return c.k }

// Redundancy returns the number of parity checks n-k.
func (c *LinearCode) Redundancy() int { return c.n - c.k }

// Rate returns the code rate k/n.
func (c *LinearCode) Rate() float64 { return float64(c.k) / float64(c.n) }

// Generator returns a copy of the generator matrix.
func (c *LinearCode) Generator() [][]int { return copyMat(c.g) }

// ParityCheck returns a copy of the parity-check matrix.
func (c *LinearCode) ParityCheck() [][]int { return copyMat(c.h) }

// Encode encodes a length-k message into a length-n code word as m*G. It
// returns ErrLength if the message length is not k.
func (c *LinearCode) Encode(msg []int) ([]int, error) {
	if len(msg) != c.k {
		return nil, ErrLength
	}
	if !validBits(msg) {
		return nil, ErrBit
	}
	return VecMatGF2(msg, c.g), nil
}

// Syndrome returns the syndrome H*rᵀ of a received length-n word. It returns
// ErrLength if the word length is not n.
func (c *LinearCode) Syndrome(word []int) ([]int, error) {
	if len(word) != c.n {
		return nil, ErrLength
	}
	return MatVecGF2(c.h, word), nil
}

// IsCodeword reports whether a length-n word is a valid code word (zero
// syndrome).
func (c *LinearCode) IsCodeword(word []int) bool {
	if len(word) != c.n {
		return false
	}
	s := MatVecGF2(c.h, word)
	return HammingWeight(s) == 0
}

// Codewords returns all 2^k code words in message-counting order. It is
// intended for small codes; it allocates 2^k slices.
func (c *LinearCode) Codewords() [][]int {
	total := 1 << uint(c.k)
	out := make([][]int, total)
	for i := 0; i < total; i++ {
		msg := make([]int, c.k)
		for b := 0; b < c.k; b++ {
			msg[c.k-1-b] = (i >> uint(b)) & 1
		}
		out[i] = VecMatGF2(msg, c.g)
	}
	return out
}

// MinimumDistance returns the minimum Hamming distance of the code, computed as
// the minimum weight over all non-zero code words. It is exponential in k and
// meant for small codes.
func (c *LinearCode) MinimumDistance() int {
	best := c.n + 1
	total := 1 << uint(c.k)
	for i := 1; i < total; i++ {
		msg := make([]int, c.k)
		for b := 0; b < c.k; b++ {
			msg[c.k-1-b] = (i >> uint(b)) & 1
		}
		w := HammingWeight(VecMatGF2(msg, c.g))
		if w < best {
			best = w
		}
	}
	return best
}

// SyndromeTable builds a standard-array coset-leader table mapping each
// syndrome (packed as a big-endian integer) to a minimum-weight error pattern.
// It is exponential in n and intended for small codes.
func (c *LinearCode) SyndromeTable() map[uint64][]int {
	table := make(map[uint64][]int)
	total := 1 << uint(c.n)
	for e := 0; e < total; e++ {
		err := make([]int, c.n)
		for b := 0; b < c.n; b++ {
			err[c.n-1-b] = (e >> uint(b)) & 1
		}
		s := BitsToUint(MatVecGF2(c.h, err))
		if cur, ok := table[s]; !ok || HammingWeight(err) < HammingWeight(cur) {
			table[s] = err
		}
	}
	return table
}

// DecodeSyndrome performs syndrome (coset-leader) decoding of a received word
// using a table built by SyndromeTable, returning the most likely code word.
// It returns ErrLength on a bad word length.
func (c *LinearCode) DecodeSyndrome(word []int, table map[uint64][]int) ([]int, error) {
	if len(word) != c.n {
		return nil, ErrLength
	}
	s := BitsToUint(MatVecGF2(c.h, word))
	err := table[s]
	if err == nil {
		return append([]int(nil), word...), nil
	}
	return XORVectors(word, err), nil
}

// copyMat returns a deep copy of a GF(2) matrix.
func copyMat(m [][]int) [][]int {
	out := make([][]int, len(m))
	for i := range m {
		out[i] = append([]int(nil), m[i]...)
	}
	return out
}
