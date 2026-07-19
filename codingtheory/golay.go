package codingtheory

// This file implements the perfect binary Golay code [23,12,7] and the
// extended Golay code [24,12,8]. The [23,12] code is realised as a cyclic code
// with generator polynomial g(x) = x^11+x^9+x^7+x^6+x^5+x+1 (the integer 2787).
// Because the code is perfect, every one of the 2^11 syndromes corresponds to a
// unique error pattern of weight at most three, which we tabulate once for
// exact three-error decoding.

// Golay23Generator is the generator polynomial of the cyclic Golay(23,12) code,
// x^11+x^9+x^7+x^6+x^5+x+1, packed as the integer 2787.
const Golay23Generator = 2787

// GolayCode is the perfect binary Golay code together with a complete
// syndrome-to-error table for exact correction of up to three bit errors.
type GolayCode struct {
	cc    *CyclicCode
	table map[int]int // syndrome polynomial -> error pattern (packed int)
}

// NewGolayCode returns the perfect binary Golay(23,12) code with its decoding
// table precomputed.
func NewGolayCode() *GolayCode {
	cc, err := NewCyclicCode(23, Golay23Generator)
	if err != nil {
		panic("codingtheory: internal Golay generator is invalid: " + err.Error())
	}
	table := make(map[int]int, 2048)
	n := 23
	// weight 0
	table[0] = 0
	// weight 1..3 error patterns
	for i := 0; i < n; i++ {
		e1 := 1 << uint(i)
		table[cc.SyndromePoly(e1)] = e1
		for j := i + 1; j < n; j++ {
			e2 := e1 | 1<<uint(j)
			table[cc.SyndromePoly(e2)] = e2
			for k := j + 1; k < n; k++ {
				e3 := e2 | 1<<uint(k)
				table[cc.SyndromePoly(e3)] = e3
			}
		}
	}
	return &GolayCode{cc: cc, table: table}
}

// N returns the code length 23.
func (g *GolayCode) N() int { return 23 }

// K returns the message dimension 12.
func (g *GolayCode) K() int { return 12 }

// MinimumDistance returns the minimum distance 7 of the Golay(23,12) code.
func (g *GolayCode) MinimumDistance() int { return 7 }

// Cyclic returns the underlying cyclic code.
func (g *GolayCode) Cyclic() *CyclicCode { return g.cc }

// Encode systematically encodes a 12-bit message into a 23-bit Golay code word.
func (g *GolayCode) Encode(msg []int) ([]int, error) { return g.cc.Encode(msg) }

// Syndrome returns the 11-bit syndrome of a received 23-bit word.
func (g *GolayCode) Syndrome(word []int) ([]int, error) { return g.cc.Syndrome(word) }

// Decode corrects up to three bit errors in a received 23-bit word and returns
// the corrected code word and the number of errors corrected.
func (g *GolayCode) Decode(word []int) (corrected []int, numErrors int, err error) {
	if len(word) != 23 {
		return nil, 0, ErrLength
	}
	if !validBits(word) {
		return nil, 0, ErrBit
	}
	r := BitsToPoly(word)
	syn := g.cc.SyndromePoly(r)
	e := g.table[syn]
	corrected = PolyToBits(r^e, 23)
	return corrected, HammingWeightUint(uint64(e)), nil
}

// DecodeMessage decodes a received 23-bit word and returns the 12 message bits.
func (g *GolayCode) DecodeMessage(word []int) (msg []int, numErrors int, err error) {
	corrected, n, err := g.Decode(word)
	if err != nil {
		return nil, 0, err
	}
	msg, _ = g.cc.DecodeMessage(corrected)
	return msg, n, nil
}

// ExtendedGolayEncode encodes a 12-bit message into a 24-bit extended Golay
// code word by appending an overall even-parity bit to the Golay(23,12) code
// word.
func ExtendedGolayEncode(g *GolayCode, msg []int) ([]int, error) {
	cw, err := g.Encode(msg)
	if err != nil {
		return nil, err
	}
	parity := 0
	for _, b := range cw {
		parity ^= b
	}
	return append(cw, parity), nil
}

// ExtendedGolayDecode decodes a 24-bit extended Golay code word. It corrects up
// to three errors using the inner Golay code and reports the number of errors
// corrected in the 23-bit prefix.
func ExtendedGolayDecode(g *GolayCode, word []int) (msg []int, numErrors int, err error) {
	if len(word) != 24 {
		return nil, 0, ErrLength
	}
	if !validBits(word) {
		return nil, 0, ErrBit
	}
	return g.DecodeMessage(word[:23])
}
