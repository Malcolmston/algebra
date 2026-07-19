package codingtheory

import "math/rand"

// This file provides helpers for low-density parity-check (LDPC) codes: a code
// wrapper around a sparse parity-check matrix, degree/density statistics, a
// Gallager regular-ensemble construction, and Gallager's hard-decision
// bit-flipping decoder.

// LDPCCode wraps a binary parity-check matrix H (m rows by n columns) over
// GF(2). A length-n word c is a code word exactly when H*cᵀ = 0.
type LDPCCode struct {
	n int
	m int
	h [][]int
	// adjacency: for each check, the variable indices it touches, and vice versa
	checkNbrs [][]int
	varNbrs   [][]int
}

// NewLDPCCode builds an LDPC code from a rectangular binary parity-check matrix.
// It returns ErrShape on an empty or ragged matrix and ErrBit on non-bit entries.
func NewLDPCCode(h [][]int) (*LDPCCode, error) {
	if len(h) == 0 || len(h[0]) == 0 {
		return nil, ErrShape
	}
	n := len(h[0])
	for _, row := range h {
		if len(row) != n {
			return nil, ErrShape
		}
		if !validBits(row) {
			return nil, ErrBit
		}
	}
	m := len(h)
	code := &LDPCCode{n: n, m: m, h: copyMat(h)}
	code.checkNbrs = make([][]int, m)
	code.varNbrs = make([][]int, n)
	for i := 0; i < m; i++ {
		for j := 0; j < n; j++ {
			if h[i][j] != 0 {
				code.checkNbrs[i] = append(code.checkNbrs[i], j)
				code.varNbrs[j] = append(code.varNbrs[j], i)
			}
		}
	}
	return code, nil
}

// N returns the code length (number of variable nodes).
func (c *LDPCCode) N() int { return c.n }

// M returns the number of parity checks (check nodes).
func (c *LDPCCode) M() int { return c.m }

// ParityCheck returns a copy of the parity-check matrix.
func (c *LDPCCode) ParityCheck() [][]int { return copyMat(c.h) }

// Syndrome returns the length-m syndrome H*wᵀ of a received length-n word.
func (c *LDPCCode) Syndrome(word []int) ([]int, error) {
	if len(word) != c.n {
		return nil, ErrLength
	}
	return MatVecGF2(c.h, word), nil
}

// IsCodeword reports whether a length-n word satisfies all parity checks.
func (c *LDPCCode) IsCodeword(word []int) bool {
	if len(word) != c.n {
		return false
	}
	for i := 0; i < c.m; i++ {
		s := 0
		for _, j := range c.checkNbrs[i] {
			s ^= word[j]
		}
		if s != 0 {
			return false
		}
	}
	return true
}

// Density returns the fraction of ones in the parity-check matrix.
func (c *LDPCCode) Density() float64 {
	ones := 0
	for j := range c.varNbrs {
		ones += len(c.varNbrs[j])
	}
	return float64(ones) / float64(c.n*c.m)
}

// VariableDegrees returns the column weights (variable-node degrees) of H.
func (c *LDPCCode) VariableDegrees() []int {
	out := make([]int, c.n)
	for j := range c.varNbrs {
		out[j] = len(c.varNbrs[j])
	}
	return out
}

// CheckDegrees returns the row weights (check-node degrees) of H.
func (c *LDPCCode) CheckDegrees() []int {
	out := make([]int, c.m)
	for i := range c.checkNbrs {
		out[i] = len(c.checkNbrs[i])
	}
	return out
}

// IsRegular reports whether every variable node has the same degree and every
// check node has the same degree.
func (c *LDPCCode) IsRegular() bool {
	vd := c.VariableDegrees()
	cd := c.CheckDegrees()
	for _, d := range vd {
		if d != vd[0] {
			return false
		}
	}
	for _, d := range cd {
		if d != cd[0] {
			return false
		}
	}
	return true
}

// BitFlipDecode runs Gallager's hard-decision bit-flipping algorithm for up to
// maxIter iterations, flipping in each round every bit that fails a majority of
// its parity checks. It returns the decoded word, the number of iterations
// used, and whether all checks were satisfied. The received word must have
// length n.
func (c *LDPCCode) BitFlipDecode(received []int, maxIter int) (decoded []int, iters int, ok bool, err error) {
	if len(received) != c.n {
		return nil, 0, false, ErrLength
	}
	if !validBits(received) {
		return nil, 0, false, ErrBit
	}
	word := append([]int(nil), received...)
	for it := 0; it < maxIter; it++ {
		// compute all check values
		checkVal := make([]int, c.m)
		unsatisfied := 0
		for i := 0; i < c.m; i++ {
			s := 0
			for _, j := range c.checkNbrs[i] {
				s ^= word[j]
			}
			checkVal[i] = s
			if s != 0 {
				unsatisfied++
			}
		}
		if unsatisfied == 0 {
			return word, it, true, nil
		}
		// for each variable, count how many of its checks are unsatisfied
		flippedAny := false
		bestVar := -1
		bestCount := 0
		for j := 0; j < c.n; j++ {
			cnt := 0
			for _, i := range c.varNbrs[j] {
				if checkVal[i] != 0 {
					cnt++
				}
			}
			if 2*cnt > len(c.varNbrs[j]) {
				word[j] ^= 1
				flippedAny = true
			}
			if cnt > bestCount {
				bestCount = cnt
				bestVar = j
			}
		}
		if !flippedAny {
			// no strict majority anywhere: flip the single most-suspect bit to
			// escape the stall, mirroring Gallager's serial variant.
			if bestVar >= 0 {
				word[bestVar] ^= 1
			} else {
				break
			}
		}
	}
	ok = c.IsCodeword(word)
	return word, maxIter, ok, nil
}

// GallagerLDPC constructs a regular (wc, wr) Gallager parity-check matrix of
// length n. It requires n to be a multiple of wr and returns an m-by-n matrix
// with m = wc*n/wr, every column of weight wc and every row of weight wr. The
// seed makes the random column permutations of the sub-bands reproducible.
func GallagerLDPC(n, wc, wr int, seed int64) ([][]int, error) {
	if wr < 1 || wc < 1 || n%wr != 0 || wc > n {
		return nil, ErrFieldParam
	}
	bandRows := n / wr
	m := wc * bandRows
	h := make([][]int, m)
	for i := range h {
		h[i] = make([]int, n)
	}
	// First band: row r has ones in columns [r*wr, (r+1)*wr).
	for r := 0; r < bandRows; r++ {
		for c := r * wr; c < (r+1)*wr; c++ {
			h[r][c] = 1
		}
	}
	rng := rand.New(rand.NewSource(seed))
	// Subsequent bands are column permutations of the first band.
	for band := 1; band < wc; band++ {
		perm := rng.Perm(n)
		for r := 0; r < bandRows; r++ {
			dst := band*bandRows + r
			for c := 0; c < n; c++ {
				h[dst][c] = h[r][perm[c]]
			}
		}
	}
	return h, nil
}
