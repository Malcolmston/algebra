package codingtheory

// This file implements the classic Hamming(7,4) code and general Hamming codes
// Hamming(2^r-1, 2^r-1-r), together with their single-error-correcting decoders
// and the extended (SECDED) variants.

// Hamming74Generator returns the 4-by-7 generator matrix of the systematic
// Hamming(7,4) code with layout [data | parity]. Row i encodes the i-th data
// bit together with the parity bits it participates in.
func Hamming74Generator() [][]int {
	return [][]int{
		{1, 0, 0, 0, 1, 1, 0},
		{0, 1, 0, 0, 1, 0, 1},
		{0, 0, 1, 0, 0, 1, 1},
		{0, 0, 0, 1, 1, 1, 1},
	}
}

// Hamming74ParityCheck returns the 3-by-7 parity-check matrix matching
// Hamming74Generator.
func Hamming74ParityCheck() [][]int {
	return [][]int{
		{1, 1, 0, 1, 1, 0, 0},
		{1, 0, 1, 1, 0, 1, 0},
		{0, 1, 1, 1, 0, 0, 1},
	}
}

// Hamming74Encode encodes four data bits into a seven-bit systematic
// Hamming(7,4) code word. It returns ErrLength if len(data) != 4.
func Hamming74Encode(data []int) ([]int, error) {
	if len(data) != 4 {
		return nil, ErrLength
	}
	if !validBits(data) {
		return nil, ErrBit
	}
	d0, d1, d2, d3 := data[0], data[1], data[2], data[3]
	p0 := d0 ^ d1 ^ d3
	p1 := d0 ^ d2 ^ d3
	p2 := d1 ^ d2 ^ d3
	return []int{d0, d1, d2, d3, p0, p1, p2}, nil
}

// Hamming74Syndrome returns the three-bit syndrome of a received Hamming(7,4)
// word. A zero syndrome indicates no detected single-bit error.
func Hamming74Syndrome(word []int) ([]int, error) {
	if len(word) != 7 {
		return nil, ErrLength
	}
	return MatVecGF2(Hamming74ParityCheck(), word), nil
}

// Hamming74Decode corrects up to one bit error in a received seven-bit word and
// returns the four recovered data bits along with the corrected code word and
// the (0-based) error position, which is -1 when the syndrome is zero.
func Hamming74Decode(word []int) (data, corrected []int, errPos int, err error) {
	if len(word) != 7 {
		return nil, nil, -1, ErrLength
	}
	if !validBits(word) {
		return nil, nil, -1, ErrBit
	}
	corrected = append([]int(nil), word...)
	s := MatVecGF2(Hamming74ParityCheck(), corrected)
	errPos = -1
	if HammingWeight(s) != 0 {
		h := Hamming74ParityCheck()
		for col := 0; col < 7; col++ {
			match := true
			for r := 0; r < 3; r++ {
				if h[r][col] != s[r] {
					match = false
					break
				}
			}
			if match {
				corrected[col] ^= 1
				errPos = col
				break
			}
		}
	}
	data = []int{corrected[0], corrected[1], corrected[2], corrected[3]}
	return data, corrected, errPos, nil
}

// HammingCode is a general binary Hamming code of length n = 2^r-1 with
// dimension k = n-r and minimum distance 3, correcting any single bit error.
// The parity-check matrix has as its columns all non-zero r-bit vectors in
// ascending numeric order, so the syndrome of a single error equals the binary
// index (1-based) of the erroneous position.
type HammingCode struct {
	r int // number of parity bits
	n int // 2^r - 1
	k int // n - r
	h [][]int
}

// NewHammingCode returns the Hamming code with r parity bits (r >= 2), of
// length 2^r-1 and dimension 2^r-1-r. It returns ErrFieldParam for r < 2.
func NewHammingCode(r int) (*HammingCode, error) {
	if r < 2 || r > 20 {
		return nil, ErrFieldParam
	}
	n := (1 << uint(r)) - 1
	h := make([][]int, r)
	for i := range h {
		h[i] = make([]int, n)
	}
	for col := 1; col <= n; col++ {
		for bit := 0; bit < r; bit++ {
			h[r-1-bit][col-1] = (col >> uint(bit)) & 1
		}
	}
	return &HammingCode{r: r, n: n, k: n - r, h: h}, nil
}

// N returns the code length 2^r-1.
func (c *HammingCode) N() int { return c.n }

// K returns the message dimension 2^r-1-r.
func (c *HammingCode) K() int { return c.k }

// R returns the number of parity bits r.
func (c *HammingCode) R() int { return c.r }

// ParityCheck returns a copy of the r-by-n parity-check matrix.
func (c *HammingCode) ParityCheck() [][]int { return copyMat(c.h) }

// Syndrome returns the r-bit syndrome H*wᵀ of a received length-n word.
func (c *HammingCode) Syndrome(word []int) ([]int, error) {
	if len(word) != c.n {
		return nil, ErrLength
	}
	return MatVecGF2(c.h, word), nil
}

// SyndromePosition returns the 0-based position of a single bit error implied
// by a received word, or -1 if the syndrome is zero. The syndrome, read as a
// big-endian binary number, is the 1-based column index of the error.
func (c *HammingCode) SyndromePosition(word []int) (int, error) {
	s, err := c.Syndrome(word)
	if err != nil {
		return -1, err
	}
	val := int(BitsToUint(s))
	if val == 0 {
		return -1, nil
	}
	return val - 1, nil
}

// Decode corrects up to one bit error in a received length-n word and returns
// the corrected code word together with the error position (-1 if none).
func (c *HammingCode) Decode(word []int) (corrected []int, errPos int, err error) {
	if len(word) != c.n {
		return nil, -1, ErrLength
	}
	if !validBits(word) {
		return nil, -1, ErrBit
	}
	corrected = append([]int(nil), word...)
	pos, _ := c.SyndromePosition(corrected)
	if pos >= 0 {
		corrected[pos] ^= 1
	}
	return corrected, pos, nil
}

// LinearCode converts the Hamming code into a general LinearCode built from its
// parity-check matrix.
func (c *HammingCode) LinearCode() (*LinearCode, error) {
	g := nullSpaceGF2(c.h, c.n)
	return NewLinearCodeGH(g, c.h)
}

// ExtendedHamming74Encode encodes four data bits into an eight-bit extended
// Hamming(8,4) SECDED code word by appending an overall even-parity bit to the
// Hamming(7,4) code word.
func ExtendedHamming74Encode(data []int) ([]int, error) {
	cw, err := Hamming74Encode(data)
	if err != nil {
		return nil, err
	}
	overall := 0
	for _, b := range cw {
		overall ^= b
	}
	return append(cw, overall), nil
}

// ExtendedHamming74Decode decodes an eight-bit extended Hamming code word. It
// corrects any single error and flags (via doubleError) an uncorrectable
// double error detected by the overall parity bit.
func ExtendedHamming74Decode(word []int) (data []int, doubleError bool, err error) {
	if len(word) != 8 {
		return nil, false, ErrLength
	}
	if !validBits(word) {
		return nil, false, ErrBit
	}
	inner := word[:7]
	overall := 0
	for _, b := range word {
		overall ^= b
	}
	s := MatVecGF2(Hamming74ParityCheck(), inner)
	synZero := HammingWeight(s) == 0
	switch {
	case synZero && overall == 0:
		return []int{word[0], word[1], word[2], word[3]}, false, nil
	case overall == 1:
		// odd overall parity: a single error occurred (possibly on the parity
		// bit); correct it via the inner syndrome.
		d, _, _, e := Hamming74Decode(append([]int(nil), inner...))
		return d, false, e
	default:
		// overall parity even but non-zero syndrome: double error detected.
		return []int{word[0], word[1], word[2], word[3]}, true, nil
	}
}
