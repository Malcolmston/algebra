package quasirandom

// PopCount returns the number of one bits in n.
func PopCount(n uint64) int {
	c := 0
	for n != 0 {
		n &= n - 1
		c++
	}
	return c
}

// BitLength returns the number of bits required to represent n, i.e. the index
// of the highest set bit plus one; BitLength(0)==0.
func BitLength(n uint64) int {
	c := 0
	for n != 0 {
		c++
		n >>= 1
	}
	return c
}

// IsPowerOfTwo reports whether n is a positive power of two.
func IsPowerOfTwo(n uint64) bool { return n != 0 && n&(n-1) == 0 }

// NextPowerOfTwo returns the smallest power of two greater than or equal to n,
// with NextPowerOfTwo(0)==1.
func NextPowerOfTwo(n uint64) uint64 {
	if n <= 1 {
		return 1
	}
	return uint64(1) << uint(BitLength(n-1))
}

// ReverseBits returns n with its low width bits reversed, the base-two analogue
// of ReverseDigits used implicitly by the base-two radical inverse.
func ReverseBits(n uint64, width int) uint64 {
	var r uint64
	for i := 0; i < width; i++ {
		r = (r << 1) | (n & 1)
		n >>= 1
	}
	return r
}
