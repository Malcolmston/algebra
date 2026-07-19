package codingtheory

// This file implements binary BCH codes over GF(2^m). A narrow-sense
// t-error-correcting BCH code of length n = 2^m-1 has generator polynomial
// g(x) = lcm of the minimal polynomials of alpha^1 … alpha^(2t); decoding uses
// syndromes over GF(2^m), the Berlekamp-Massey algorithm and a Chien search.

// BCHCode is a binary BCH code of length n = 2^m-1 designed to correct t
// errors. It is a cyclic code whose generator polynomial is the least common
// multiple of the minimal polynomials of the roots.
type BCHCode struct {
	field *Field
	n     int
	k     int
	t     int
	b     int // first root exponent (1 for narrow-sense)
	cc    *CyclicCode
}

// NewBCHCode returns the narrow-sense binary BCH code over the given field that
// corrects up to t errors. It returns ErrFieldParam if t is too large for the
// field.
func NewBCHCode(field *Field, t int) (*BCHCode, error) {
	return NewBCHCodeB(field, t, 1)
}

// NewBCHCodeB returns a binary BCH code over the given field correcting up to t
// errors with first root exponent b (b=1 gives the narrow-sense code). The
// generator has 2t consecutive roots alpha^b … alpha^(b+2t-1).
func NewBCHCodeB(field *Field, t, b int) (*BCHCode, error) {
	n := field.Order()
	if t < 1 || 2*t >= n || b < 1 {
		return nil, ErrFieldParam
	}
	// g = lcm of minimal polynomials of alpha^(b) .. alpha^(b+2t-1)
	seen := map[int]bool{}
	gen := 1
	for i := 0; i < 2*t; i++ {
		e := b + i
		mp := field.MinimalPoly(e)
		if mp == 0 {
			return nil, ErrFieldParam
		}
		// deduplicate identical minimal polynomials via cyclotomic coset leader
		leader := cosetLeader(e, n)
		if seen[leader] {
			continue
		}
		seen[leader] = true
		gen = gf2Lcm(gen, mp)
	}
	cc, err := NewCyclicCode(n, gen)
	if err != nil {
		return nil, err
	}
	return &BCHCode{field: field, n: n, k: cc.K(), t: t, b: b, cc: cc}, nil
}

// cosetLeader returns the smallest element of the cyclotomic coset of e modulo
// n under multiplication by two.
func cosetLeader(e, n int) int {
	e = ((e % n) + n) % n
	min := e
	s := e
	for {
		s = (s * 2) % n
		if s == e {
			break
		}
		if s < min {
			min = s
		}
	}
	return min
}

// gf2Lcm returns the least common multiple of two GF(2) polynomials.
func gf2Lcm(a, b int) int {
	if a == 0 || b == 0 {
		return 0
	}
	g := GF2PolyGCD(a, b)
	return GF2PolyMul(GF2PolyDiv(a, g), b)
}

// Field returns the underlying field.
func (c *BCHCode) Field() *Field { return c.field }

// N returns the code length 2^m-1.
func (c *BCHCode) N() int { return c.n }

// K returns the message dimension.
func (c *BCHCode) K() int { return c.k }

// T returns the designed error-correction capability.
func (c *BCHCode) T() int { return c.t }

// DesignedDistance returns the designed minimum distance 2t+1.
func (c *BCHCode) DesignedDistance() int { return 2*c.t + 1 }

// Generator returns the generator polynomial packed as a GF(2) polynomial int.
func (c *BCHCode) Generator() int { return c.cc.Generator() }

// Cyclic returns the underlying cyclic code.
func (c *BCHCode) Cyclic() *CyclicCode { return c.cc }

// Encode systematically encodes a length-k message into a length-n code word.
func (c *BCHCode) Encode(msg []int) ([]int, error) { return c.cc.Encode(msg) }

// Syndromes returns the 2t syndrome values S_j = R(alpha^(b+j)) for
// j = 0 … 2t-1 of a received length-n binary word.
func (c *BCHCode) Syndromes(received []int) []int {
	f := c.field
	s := make([]int, 2*c.t)
	for j := 0; j < 2*c.t; j++ {
		exp := c.b + j
		var acc int
		for i, bit := range received {
			if bit&1 != 0 {
				acc ^= f.Exp(exp * i)
			}
		}
		s[j] = acc
	}
	return s
}

// Decode corrects up to t bit errors in a received length-n word, returning the
// corrected code word and the number of corrected bits. It returns
// ErrTooManyErrors on an uncorrectable pattern.
func (c *BCHCode) Decode(received []int) (corrected []int, numErrors int, err error) {
	if len(received) != c.n {
		return nil, 0, ErrLength
	}
	if !validBits(received) {
		return nil, 0, ErrBit
	}
	f := c.field
	s := c.Syndromes(received)
	if HammingWeight(s) == 0 {
		return append([]int(nil), received...), 0, nil
	}
	lambda := f.BerlekampMassey(s)
	positions := f.ChienSearch(lambda, c.n)
	nerr := len(lambda) - 1
	if len(positions) != nerr || nerr == 0 || nerr > c.t {
		return nil, 0, ErrTooManyErrors
	}
	corrected = append([]int(nil), received...)
	for _, pos := range positions {
		corrected[pos] ^= 1
	}
	// verify all syndromes now vanish
	if HammingWeight(c.Syndromes(corrected)) != 0 {
		return nil, 0, ErrTooManyErrors
	}
	return corrected, len(positions), nil
}

// DecodeMessage decodes a received word and returns the k message bits.
func (c *BCHCode) DecodeMessage(received []int) (msg []int, numErrors int, err error) {
	corrected, n, err := c.Decode(received)
	if err != nil {
		return nil, 0, err
	}
	m, _ := c.cc.DecodeMessage(corrected)
	return m, n, nil
}
