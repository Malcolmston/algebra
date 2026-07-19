package quasirandom

// Digits returns the base-b representation of the non-negative integer n as a
// slice of digits ordered from least significant to most significant. The
// value zero yields an empty slice. It returns an error when base < 2 or n is
// negative.
func Digits(n uint64, base int) ([]int, error) {
	if base < 2 {
		return nil, ErrBadBase
	}
	b := uint64(base)
	var out []int
	for n > 0 {
		out = append(out, int(n%b))
		n /= b
	}
	return out, nil
}

// DigitsBig behaves like Digits but always returns a slice of length at least
// width, left-padding the most-significant end with zeros. A width of zero is
// equivalent to Digits.
func DigitsPadded(n uint64, base, width int) ([]int, error) {
	d, err := Digits(n, base)
	if err != nil {
		return nil, err
	}
	for len(d) < width {
		d = append(d, 0)
	}
	return d, nil
}

// FromDigits reconstructs the integer whose base-b digits (least significant
// first) are given. It returns an error when base < 2 or a digit is out of the
// range [0,base).
func FromDigits(digits []int, base int) (uint64, error) {
	if base < 2 {
		return 0, ErrBadBase
	}
	b := uint64(base)
	var n, p uint64 = 0, 1
	for _, d := range digits {
		if d < 0 || d >= base {
			return 0, ErrBadBase
		}
		n += uint64(d) * p
		p *= b
	}
	return n, nil
}

// DigitCount returns the number of base-b digits required to represent n, with
// zero requiring zero digits. It returns an error when base < 2.
func DigitCount(n uint64, base int) (int, error) {
	d, err := Digits(n, base)
	if err != nil {
		return 0, err
	}
	return len(d), nil
}

// DigitSum returns the sum of the base-b digits of n. It returns an error when
// base < 2.
func DigitSum(n uint64, base int) (int, error) {
	d, err := Digits(n, base)
	if err != nil {
		return 0, err
	}
	s := 0
	for _, x := range d {
		s += x
	}
	return s, nil
}

// ReverseDigits returns the integer obtained by reversing the base-b digit
// string of n over exactly width digits. It is the integer that the radical
// inverse scales by base^-width. It returns an error when base < 2.
func ReverseDigits(n uint64, base, width int) (uint64, error) {
	d, err := DigitsPadded(n, base, width)
	if err != nil {
		return 0, err
	}
	if len(d) > width {
		width = len(d)
	}
	rev := make([]int, width)
	for i := 0; i < width; i++ {
		rev[width-1-i] = d[i]
	}
	return FromDigits(rev, base)
}

// DigitAt returns the base-b digit of n at position pos (position zero is the
// least significant digit). Positions beyond the length of n yield zero. It
// returns an error when base < 2 or pos is negative.
func DigitAt(n uint64, base, pos int) (int, error) {
	if base < 2 {
		return 0, ErrBadBase
	}
	if pos < 0 {
		return 0, ErrNonPositive
	}
	b := uint64(base)
	for i := 0; i < pos && n > 0; i++ {
		n /= b
	}
	return int(n % b), nil
}

// GrayCode returns the reflected binary Gray code of n, namely n XOR (n>>1).
func GrayCode(n uint64) uint64 { return n ^ (n >> 1) }

// InverseGrayCode returns the integer whose Gray code is g, inverting GrayCode.
func InverseGrayCode(g uint64) uint64 {
	var n uint64
	for ; g != 0; g >>= 1 {
		n ^= g
	}
	return n
}

// TrailingZeroBits returns the number of trailing zero bits of n, i.e. the
// index (from zero) of the lowest set bit. For n==0 it returns 64.
func TrailingZeroBits(n uint64) int {
	if n == 0 {
		return 64
	}
	c := 0
	for n&1 == 0 {
		n >>= 1
		c++
	}
	return c
}

// LowestZeroBit returns the position (from zero) of the lowest zero bit of n,
// which drives the Antonov–Saleev Gray-code Sobol recurrence.
func LowestZeroBit(n uint64) int {
	c := 0
	for n&1 == 1 {
		n >>= 1
		c++
	}
	return c
}
