package codingtheory

import "errors"

// ErrTooManyErrors is returned when a decoder detects more symbol errors than
// it can correct.
var ErrTooManyErrors = errors.New("codingtheory: too many errors to correct")

// ReedSolomon is a Reed-Solomon code over GF(2^m). It has length n = 2^m-1,
// dimension k = n - nsym and can correct up to nsym/2 symbol errors. The nsym
// generator roots are consecutive powers alpha^fcr … alpha^(fcr+nsym-1) of the
// field's primitive element.
//
// Messages and code words are big-endian []int symbol slices: index 0 holds
// the first (highest-order) symbol. Systematic code words carry the k message
// symbols first followed by nsym parity symbols.
type ReedSolomon struct {
	field *Field
	n     int
	k     int
	nsym  int
	fcr   int
	genLE []int // generator polynomial, little-endian (index = power of x)
}

// NewReedSolomon returns a Reed-Solomon code over the given field with nsym
// parity symbols and first consecutive root exponent fcr=0. It returns
// ErrFieldParam if nsym is not in [1, 2^m-1].
func NewReedSolomon(field *Field, nsym int) (*ReedSolomon, error) {
	return NewReedSolomonFCR(field, nsym, 0)
}

// NewReedSolomonFCR returns a Reed-Solomon code with nsym parity symbols and a
// configurable first-consecutive-root exponent fcr.
func NewReedSolomonFCR(field *Field, nsym, fcr int) (*ReedSolomon, error) {
	n := field.Order()
	if nsym < 1 || nsym >= n {
		return nil, ErrFieldParam
	}
	rs := &ReedSolomon{field: field, n: n, k: n - nsym, nsym: nsym, fcr: fcr}
	rs.genLE = rs.generatorPoly()
	return rs, nil
}

// generatorPoly builds g(x) = prod_{i=0}^{nsym-1} (x - alpha^(fcr+i)) as a
// little-endian monic polynomial.
func (rs *ReedSolomon) generatorPoly() []int {
	f := rs.field
	g := []int{1}
	for i := 0; i < rs.nsym; i++ {
		root := f.Exp(rs.fcr + i)
		// multiply g by (x - root) = (x + root) = little-endian [root, 1]
		g = f.PolyMul(g, []int{root, 1})
	}
	return g
}

// Field returns the underlying field.
func (rs *ReedSolomon) Field() *Field { return rs.field }

// N returns the code length 2^m-1.
func (rs *ReedSolomon) N() int { return rs.n }

// K returns the message dimension n-nsym.
func (rs *ReedSolomon) K() int { return rs.k }

// NSym returns the number of parity symbols.
func (rs *ReedSolomon) NSym() int { return rs.nsym }

// FCR returns the first-consecutive-root exponent.
func (rs *ReedSolomon) FCR() int { return rs.fcr }

// Correction returns the maximum number of symbol errors the code can correct,
// nsym/2.
func (rs *ReedSolomon) Correction() int { return rs.nsym / 2 }

// GeneratorPoly returns the generator polynomial as a big-endian symbol slice
// (index 0 is the highest-degree coefficient).
func (rs *ReedSolomon) GeneratorPoly() []int { return reverseInts(rs.genLE) }

// Encode systematically encodes a length-k message into a length-n code word.
// It returns ErrLength if the message length is not k.
func (rs *ReedSolomon) Encode(msg []int) ([]int, error) {
	if len(msg) != rs.k {
		return nil, ErrLength
	}
	f := rs.field
	// little-endian message poly: msg[0] is highest-degree symbol.
	mLE := reverseInts(msg)
	// shift up by nsym: multiply by x^nsym
	shifted := make([]int, len(mLE)+rs.nsym)
	copy(shifted[rs.nsym:], mLE)
	parity := f.PolyMod(shifted, rs.genLE)
	cwLE := f.PolyAdd(shifted, parity)
	// pad to length n
	full := make([]int, rs.n)
	copy(full, cwLE)
	return reverseInts(full), nil
}

// Syndromes returns the nsym syndrome values S_j = R(alpha^(fcr+j)) of a
// received length-n word, in order j = 0 … nsym-1.
func (rs *ReedSolomon) Syndromes(received []int) []int {
	f := rs.field
	rLE := reverseInts(received)
	s := make([]int, rs.nsym)
	for j := 0; j < rs.nsym; j++ {
		s[j] = f.PolyEval(rLE, f.Exp(rs.fcr+j))
	}
	return s
}

// Decode corrects up to nsym/2 symbol errors in a received length-n word,
// returning the corrected code word and the number of corrected symbols. It
// returns ErrTooManyErrors when the error pattern is uncorrectable.
func (rs *ReedSolomon) Decode(received []int) (corrected []int, numErrors int, err error) {
	if len(received) != rs.n {
		return nil, 0, ErrLength
	}
	f := rs.field
	s := rs.Syndromes(received)
	if HammingWeight(s) == 0 {
		return append([]int(nil), received...), 0, nil
	}
	lambda := f.BerlekampMassey(s)
	positions := f.ChienSearch(lambda, rs.n)
	nerr := len(lambda) - 1
	if len(positions) != nerr || nerr == 0 || nerr > rs.Correction() {
		return nil, 0, ErrTooManyErrors
	}
	// Error evaluator Omega(x) = S(x)*Lambda(x) mod x^nsym.
	sPoly := append([]int(nil), s...) // little-endian S_0 + S_1 x + ...
	omega := f.PolyMul(sPoly, lambda)
	if len(omega) > rs.nsym {
		omega = omega[:rs.nsym]
	}
	lambdaDeriv := f.PolyDeriv(lambda)
	rLE := reverseInts(received)
	out := make([]int, len(rLE))
	copy(out, rLE)
	for _, pos := range positions {
		x := f.Exp(pos)  // X_l = alpha^pos
		xinv := f.Inv(x) // X_l^{-1}
		num := f.PolyEval(omega, xinv)
		den := f.PolyEval(lambdaDeriv, xinv)
		if den == 0 {
			return nil, 0, ErrTooManyErrors
		}
		// Forney: magnitude = X^(1-fcr) * Omega(X^-1) / Lambda'(X^-1)
		mag := f.Mul(f.Div(num, den), f.Pow(x, 1-rs.fcr))
		out[pos] ^= mag
	}
	return reverseInts(out), len(positions), nil
}

// DecodeMessage decodes a received word and returns the k recovered message
// symbols.
func (rs *ReedSolomon) DecodeMessage(received []int) (msg []int, numErrors int, err error) {
	corrected, n, err := rs.Decode(received)
	if err != nil {
		return nil, 0, err
	}
	return corrected[:rs.k], n, nil
}

// IsCodeword reports whether a received length-n word has all-zero syndromes.
func (rs *ReedSolomon) IsCodeword(word []int) bool {
	if len(word) != rs.n {
		return false
	}
	return HammingWeight(rs.Syndromes(word)) == 0
}

// reverseInts returns a reversed copy of s.
func reverseInts(s []int) []int {
	out := make([]int, len(s))
	for i := range s {
		out[len(s)-1-i] = s[i]
	}
	return out
}

// polyShiftLE multiplies a little-endian polynomial by x^m (prepends m zeros).
func polyShiftLE(p []int, m int) []int {
	out := make([]int, len(p)+m)
	copy(out[m:], p)
	return out
}
