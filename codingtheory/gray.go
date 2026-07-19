package codingtheory

// This file implements binary reflected Gray codes and n-ary Gray codes.

// BinaryToGray converts a non-negative integer to its reflected binary Gray
// code, in which consecutive values differ in exactly one bit.
func BinaryToGray(n uint64) uint64 { return n ^ (n >> 1) }

// GrayToBinary converts a reflected binary Gray code back to the ordinary
// binary integer it encodes.
func GrayToBinary(g uint64) uint64 {
	var b uint64
	for g != 0 {
		b ^= g
		g >>= 1
	}
	return b
}

// GrayCodeBits returns the Gray code of value with a fixed width of nbits,
// most-significant bit first.
func GrayCodeBits(value uint64, nbits int) []int {
	return UintToBits(BinaryToGray(value), nbits)
}

// GrayBitsToBinary converts a fixed-width Gray-code bit slice (most-significant
// bit first) back to its integer value.
func GrayBitsToBinary(bits []int) uint64 {
	return GrayToBinary(BitsToUint(bits))
}

// GrayCodeSequence returns the length-2^n reflected Gray code sequence as a
// slice of integers, where consecutive entries differ in exactly one bit.
func GrayCodeSequence(n int) []uint64 {
	total := 1 << uint(n)
	out := make([]uint64, total)
	for i := 0; i < total; i++ {
		out[i] = BinaryToGray(uint64(i))
	}
	return out
}

// GrayCodeTable returns the length-2^n reflected Gray code sequence as bit
// slices of width n (most-significant bit first).
func GrayCodeTable(n int) [][]int {
	seq := GrayCodeSequence(n)
	out := make([][]int, len(seq))
	for i, g := range seq {
		out[i] = UintToBits(g, n)
	}
	return out
}

// NaryToGray converts a value expressed as base-radix digits (most-significant
// digit first) to its modular n-ary Gray code digits: consecutive values differ
// in exactly one digit, and that digit changes by one modulo radix. The input
// digits are not modified.
func NaryToGray(digits []int, radix int) []int {
	out := make([]int, len(digits))
	shift := 0
	// process from most significant to least significant
	for i := 0; i < len(digits); i++ {
		out[i] = (digits[i] + shift) % radix
		shift = (shift + radix - out[i]) % radix
	}
	return out
}

// GrayToNary inverts NaryToGray, converting modular n-ary Gray digits
// (most-significant first) back to ordinary base-radix digits.
func GrayToNary(gray []int, radix int) []int {
	out := make([]int, len(gray))
	shift := 0
	for i := 0; i < len(gray); i++ {
		out[i] = ((gray[i]-shift)%radix + radix) % radix
		shift = (shift + radix - gray[i]) % radix
	}
	return out
}
