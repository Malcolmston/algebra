package discretemath

import "strings"

// discretemathDeBruijn64 is the 64-bit De Bruijn constant B(2, 6) used to build
// a perfect-hash table for isolating the index of the lowest set bit.
const discretemathDeBruijn64 = 0x03f79d71b4ca8b09

// discretemathDeBruijnTable64 maps the top six bits of (isolatedBit * constant)
// to that bit's position, giving a branch-free trailing-zero count.
var discretemathDeBruijnTable64 [64]int

func init() {
	for i := 0; i < 64; i++ {
		discretemathDeBruijnTable64[(uint64(1)<<uint(i)*discretemathDeBruijn64)>>58] = i
	}
}

// PopCount returns the number of one bits (the population count) in x.
//
// It uses the classic SWAR (SIMD-within-a-register) algorithm and runs in a
// fixed number of operations independent of the value of x.
func PopCount(x uint64) int {
	x = x - ((x >> 1) & 0x5555555555555555)
	x = (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333)
	x = (x + (x >> 4)) & 0x0f0f0f0f0f0f0f0f
	return int((x * 0x0101010101010101) >> 56)
}

// Parity returns the parity of x: 1 if the number of set bits is odd and 0 if
// it is even.
func Parity(x uint64) int {
	x ^= x >> 32
	x ^= x >> 16
	x ^= x >> 8
	x ^= x >> 4
	x ^= x >> 2
	x ^= x >> 1
	return int(x & 1)
}

// HammingDistance returns the number of bit positions at which a and b differ,
// which is the population count of their exclusive-or.
func HammingDistance(a, b uint64) int {
	return PopCount(a ^ b)
}

// HammingDistanceBytes returns the number of differing bits between two equal
// length byte slices. It returns an error if the slices differ in length.
func HammingDistanceBytes(a, b []byte) (int, error) {
	if len(a) != len(b) {
		return 0, discretemathErrorf("HammingDistanceBytes: length mismatch %d != %d", len(a), len(b))
	}
	d := 0
	for i := range a {
		d += PopCount(uint64(a[i] ^ b[i]))
	}
	return d, nil
}

// GrayEncode converts the binary value x into its reflected binary Gray code,
// in which successive integers differ in exactly one bit.
func GrayEncode(x uint64) uint64 {
	return x ^ (x >> 1)
}

// GrayDecode converts a reflected binary Gray code g back into the ordinary
// binary value it represents. It is the inverse of GrayEncode.
func GrayDecode(g uint64) uint64 {
	g ^= g >> 32
	g ^= g >> 16
	g ^= g >> 8
	g ^= g >> 4
	g ^= g >> 2
	g ^= g >> 1
	return g
}

// discretemathPart1By1 spreads the low 32 bits of x so that bit i moves to
// position 2*i, leaving the odd positions clear (used for 2-D Morton codes).
func discretemathPart1By1(x uint32) uint64 {
	v := uint64(x)
	v = (v | (v << 16)) & 0x0000ffff0000ffff
	v = (v | (v << 8)) & 0x00ff00ff00ff00ff
	v = (v | (v << 4)) & 0x0f0f0f0f0f0f0f0f
	v = (v | (v << 2)) & 0x3333333333333333
	v = (v | (v << 1)) & 0x5555555555555555
	return v
}

// discretemathCompact1By1 is the inverse of discretemathPart1By1: it gathers the
// even-position bits of v back into the low 32 bits.
func discretemathCompact1By1(v uint64) uint32 {
	v &= 0x5555555555555555
	v = (v | (v >> 1)) & 0x3333333333333333
	v = (v | (v >> 2)) & 0x0f0f0f0f0f0f0f0f
	v = (v | (v >> 4)) & 0x00ff00ff00ff00ff
	v = (v | (v >> 8)) & 0x0000ffff0000ffff
	v = (v | (v >> 16)) & 0x00000000ffffffff
	return uint32(v)
}

// MortonEncode2D interleaves the bits of x and y into a single Morton (Z-order)
// code, with the bits of x occupying the even positions and the bits of y the
// odd positions.
func MortonEncode2D(x, y uint32) uint64 {
	return discretemathPart1By1(x) | (discretemathPart1By1(y) << 1)
}

// MortonDecode2D recovers the two coordinates x and y from a 2-D Morton code m.
// It is the inverse of MortonEncode2D.
func MortonDecode2D(m uint64) (x, y uint32) {
	return discretemathCompact1By1(m), discretemathCompact1By1(m >> 1)
}

// discretemathPart1By2 spreads the low 21 bits of x so that bit i moves to
// position 3*i (used for 3-D Morton codes).
func discretemathPart1By2(x uint32) uint64 {
	v := uint64(x) & 0x1fffff
	v = (v | (v << 32)) & 0x1f00000000ffff
	v = (v | (v << 16)) & 0x1f0000ff0000ff
	v = (v | (v << 8)) & 0x100f00f00f00f00f
	v = (v | (v << 4)) & 0x10c30c30c30c30c3
	v = (v | (v << 2)) & 0x1249249249249249
	return v
}

// discretemathCompact1By2 is the inverse of discretemathPart1By2: it gathers
// every third bit of v back into the low 21 bits.
func discretemathCompact1By2(v uint64) uint32 {
	v &= 0x1249249249249249
	v = (v | (v >> 2)) & 0x10c30c30c30c30c3
	v = (v | (v >> 4)) & 0x100f00f00f00f00f
	v = (v | (v >> 8)) & 0x1f0000ff0000ff
	v = (v | (v >> 16)) & 0x1f00000000ffff
	v = (v | (v >> 32)) & 0x1fffff
	return uint32(v)
}

// MortonEncode3D interleaves the low 21 bits of x, y and z into a single 3-D
// Morton (Z-order) code. Only the low 21 bits of each coordinate are used, so
// the result fits in 63 bits.
func MortonEncode3D(x, y, z uint32) uint64 {
	return discretemathPart1By2(x) | (discretemathPart1By2(y) << 1) | (discretemathPart1By2(z) << 2)
}

// MortonDecode3D recovers the three coordinates x, y and z from a 3-D Morton
// code m. It is the inverse of MortonEncode3D.
func MortonDecode3D(m uint64) (x, y, z uint32) {
	return discretemathCompact1By2(m), discretemathCompact1By2(m >> 1), discretemathCompact1By2(m >> 2)
}

// ReverseBits returns x with the order of its 64 bits reversed, so that bit i of
// the result equals bit 63-i of x.
func ReverseBits(x uint64) uint64 {
	x = ((x & 0x5555555555555555) << 1) | ((x >> 1) & 0x5555555555555555)
	x = ((x & 0x3333333333333333) << 2) | ((x >> 2) & 0x3333333333333333)
	x = ((x & 0x0f0f0f0f0f0f0f0f) << 4) | ((x >> 4) & 0x0f0f0f0f0f0f0f0f)
	x = ((x & 0x00ff00ff00ff00ff) << 8) | ((x >> 8) & 0x00ff00ff00ff00ff)
	x = ((x & 0x0000ffff0000ffff) << 16) | ((x >> 16) & 0x0000ffff0000ffff)
	x = (x << 32) | (x >> 32)
	return x
}

// ReverseBits32 returns x with the order of its 32 bits reversed.
func ReverseBits32(x uint32) uint32 {
	x = ((x & 0x55555555) << 1) | ((x >> 1) & 0x55555555)
	x = ((x & 0x33333333) << 2) | ((x >> 2) & 0x33333333)
	x = ((x & 0x0f0f0f0f) << 4) | ((x >> 4) & 0x0f0f0f0f)
	x = ((x & 0x00ff00ff) << 8) | ((x >> 8) & 0x00ff00ff)
	x = (x << 16) | (x >> 16)
	return x
}

// RotateLeft returns x rotated left by k bits within a 64-bit word. Negative
// values of k rotate right; k is reduced modulo 64.
func RotateLeft(x uint64, k int) uint64 {
	n := uint(((k % 64) + 64) % 64)
	return x<<n | x>>(64-n)
}

// RotateRight returns x rotated right by k bits within a 64-bit word. Negative
// values of k rotate left; k is reduced modulo 64.
func RotateRight(x uint64, k int) uint64 {
	return RotateLeft(x, -k)
}

// TrailingZeros returns the number of consecutive zero bits at the low end of x,
// counted using a De Bruijn perfect-hash table. It returns 64 when x is zero.
func TrailingZeros(x uint64) int {
	if x == 0 {
		return 64
	}
	isolated := x & -x
	return discretemathDeBruijnTable64[(isolated*discretemathDeBruijn64)>>58]
}

// LeadingZeros returns the number of consecutive zero bits at the high end of x.
// It returns 64 when x is zero.
func LeadingZeros(x uint64) int {
	if x == 0 {
		return 64
	}
	x |= x >> 1
	x |= x >> 2
	x |= x >> 4
	x |= x >> 8
	x |= x >> 16
	x |= x >> 32
	return 64 - PopCount(x)
}

// Log2Floor returns the floor of the base-2 logarithm of x, equivalently the
// index of its highest set bit. It returns -1 when x is zero.
func Log2Floor(x uint64) int {
	if x == 0 {
		return -1
	}
	return 63 - LeadingZeros(x)
}

// Log2Ceil returns the ceiling of the base-2 logarithm of x, equivalently the
// smallest e such that 2**e >= x. It returns -1 when x is zero and 0 when x is
// one.
func Log2Ceil(x uint64) int {
	if x == 0 {
		return -1
	}
	if x == 1 {
		return 0
	}
	return Log2Floor(x-1) + 1
}

// IsPowerOfTwo reports whether x is a positive power of two (1, 2, 4, 8, ...).
func IsPowerOfTwo(x uint64) bool {
	return x != 0 && x&(x-1) == 0
}

// NextPowerOfTwo returns the smallest power of two that is greater than or equal
// to x. It returns 1 for x of 0 or 1. The result is undefined (overflows to 0)
// when x exceeds the largest representable power of two.
func NextPowerOfTwo(x uint64) uint64 {
	if x <= 1 {
		return 1
	}
	x--
	x |= x >> 1
	x |= x >> 2
	x |= x >> 4
	x |= x >> 8
	x |= x >> 16
	x |= x >> 32
	return x + 1
}

// PrevPowerOfTwo returns the largest power of two that is less than or equal to
// x. It returns 0 when x is zero.
func PrevPowerOfTwo(x uint64) uint64 {
	if x == 0 {
		return 0
	}
	return uint64(1) << uint(Log2Floor(x))
}

// SetBit returns x with bit i set to one. Bits are numbered from zero at the
// least-significant end.
func SetBit(x uint64, i uint) uint64 {
	return x | (uint64(1) << i)
}

// ClearBit returns x with bit i set to zero. Bits are numbered from zero at the
// least-significant end.
func ClearBit(x uint64, i uint) uint64 {
	return x &^ (uint64(1) << i)
}

// ToggleBit returns x with bit i flipped. Bits are numbered from zero at the
// least-significant end.
func ToggleBit(x uint64, i uint) uint64 {
	return x ^ (uint64(1) << i)
}

// TestBit reports whether bit i of x is set. Bits are numbered from zero at the
// least-significant end.
func TestBit(x uint64, i uint) bool {
	return x&(uint64(1)<<i) != 0
}

// LowestSetBit returns x with all bits cleared except its lowest set bit. It
// returns 0 when x is zero.
func LowestSetBit(x uint64) uint64 {
	return x & -x
}

// ClearLowestSetBit returns x with its lowest set bit cleared. It returns 0 when
// x is zero.
func ClearLowestSetBit(x uint64) uint64 {
	return x & (x - 1)
}

// BitString returns the base-2 representation of x, left-padded with zeros to at
// least width characters. A width of zero or less yields the minimal-length
// representation.
func BitString(x uint64, width int) string {
	var buf [64]byte
	n := 0
	for i := 63; i >= 0; i-- {
		if x&(uint64(1)<<uint(i)) != 0 {
			buf[n] = '1'
		} else {
			buf[n] = '0'
		}
		n++
	}
	// Trim leading zeros down to the minimal representation.
	start := 0
	for start < 63 && buf[start] == '0' {
		start++
	}
	s := string(buf[start:])
	if pad := width - len(s); pad > 0 {
		s = strings.Repeat("0", pad) + s
	}
	return s
}
