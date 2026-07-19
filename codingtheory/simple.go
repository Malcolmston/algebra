package codingtheory

// This file implements a handful of elementary block codes: repetition codes,
// the single-parity-check code, and majority-logic helpers.

// RepetitionEncode repeats each message bit n times, producing a code word of
// length n*len(msg). It returns ErrBit on a non-bit input.
func RepetitionEncode(msg []int, n int) ([]int, error) {
	if n < 1 {
		return nil, ErrFieldParam
	}
	if !validBits(msg) {
		return nil, ErrBit
	}
	out := make([]int, 0, len(msg)*n)
	for _, b := range msg {
		for i := 0; i < n; i++ {
			out = append(out, b)
		}
	}
	return out, nil
}

// RepetitionDecode decodes a repetition-coded word by majority vote over each
// block of n symbols, recovering the original message. The word length must be
// a multiple of n.
func RepetitionDecode(word []int, n int) ([]int, error) {
	if n < 1 || len(word)%n != 0 {
		return nil, ErrLength
	}
	if !validBits(word) {
		return nil, ErrBit
	}
	blocks := len(word) / n
	out := make([]int, blocks)
	for b := 0; b < blocks; b++ {
		ones := 0
		for i := 0; i < n; i++ {
			ones += word[b*n+i]
		}
		if 2*ones > n {
			out[b] = 1
		}
	}
	return out, nil
}

// MajorityVote returns the majority bit of a bit slice; ties resolve to 0.
func MajorityVote(bits []int) int {
	ones := 0
	for _, b := range bits {
		ones += b & 1
	}
	if 2*ones > len(bits) {
		return 1
	}
	return 0
}

// ParityBit returns the even-parity bit of a bit slice (the exclusive-or of all
// entries), so that appending it makes the total parity even.
func ParityBit(bits []int) int {
	p := 0
	for _, b := range bits {
		p ^= b & 1
	}
	return p
}

// SingleParityEncode appends an even-parity bit to the message, producing a
// length-(k+1) code word of even weight.
func SingleParityEncode(msg []int) ([]int, error) {
	if !validBits(msg) {
		return nil, ErrBit
	}
	return append(append([]int(nil), msg...), ParityBit(msg)), nil
}

// SingleParityCheck reports whether a word has even parity (no detected single
// error).
func SingleParityCheck(word []int) bool { return ParityBit(word) == 0 }

// SingleParityDecode strips the trailing parity bit and reports whether the
// received word passed the even-parity check (errorDetected is true when an
// odd number of bit errors is present).
func SingleParityDecode(word []int) (msg []int, errorDetected bool, err error) {
	if len(word) < 1 {
		return nil, false, ErrLength
	}
	if !validBits(word) {
		return nil, false, ErrBit
	}
	return append([]int(nil), word[:len(word)-1]...), ParityBit(word) != 0, nil
}
