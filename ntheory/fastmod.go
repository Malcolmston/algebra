package ntheory

import "math/bits"

// This file provides stdlib-only (math/bits) fast modular primitives for
// unsigned 64-bit integers. They give exact uint64 fast paths where the rest of
// the package drops down to math/big, and they are the foundation the other
// performance-sensitive number-theory routines (primality, factorization,
// primitive roots, modular square roots) build on.
//
// All results are normalized to the half-open range [0, m). Every routine is
// deterministic and panics on a non-positive modulus rather than returning a
// silently wrong value.

// MulModU64 returns a*b mod m computed from the full 128-bit product with
// math/bits (bits.Mul64 followed by bits.Div64), so no intermediate value
// overflows and no math/big allocation is required.
//
// m must be positive; MulModU64 panics if m == 0. The high word of the product
// is reduced modulo m first (so it is strictly less than m) to keep bits.Div64
// from overflowing. The result lies in [0, m).
func MulModU64(a, b, m uint64) uint64 {
	if m == 0 {
		panic("ntheory: MulModU64 requires a positive modulus")
	}
	hi, lo := bits.Mul64(a, b)
	// Ensure the high word is < m so bits.Div64 cannot overflow. This is exactly
	// the "if the high word >= m it reduces first" step.
	hi %= m
	_, rem := bits.Div64(hi, lo, m)
	return rem
}

// AddModU64 returns (a + b) mod m without overflowing, even when a + b does not
// fit in a uint64. m must be positive; AddModU64 panics if m == 0. The result
// lies in [0, m).
func AddModU64(a, b, m uint64) uint64 {
	if m == 0 {
		panic("ntheory: AddModU64 requires a positive modulus")
	}
	a %= m
	b %= m
	// If b >= m-a then a+b >= m and the true result is a+b-m == b-(m-a), computed
	// without ever forming the overflowing sum a+b.
	if b >= m-a {
		return b - (m - a)
	}
	return a + b
}

// SubModU64 returns (a - b) mod m without underflowing. m must be positive;
// SubModU64 panics if m == 0. The result lies in [0, m).
func SubModU64(a, b, m uint64) uint64 {
	if m == 0 {
		panic("ntheory: SubModU64 requires a positive modulus")
	}
	a %= m
	b %= m
	if a >= b {
		return a - b
	}
	return m - (b - a)
}

// ModPowU64 returns base**exp mod m by square-and-multiply. base is reduced
// modulo m first. m must be positive; ModPowU64 panics if m == 0, and returns 0
// when m == 1.
//
// When m is odd it multiplies through a Montgomery context (division-free REDC
// in the inner loop); for even m it falls back to square-and-multiply over
// [MulModU64]. The result lies in [0, m).
func ModPowU64(base, exp, m uint64) uint64 {
	if m == 0 {
		panic("ntheory: ModPowU64 requires a positive modulus")
	}
	if m == 1 {
		return 0
	}
	base %= m
	if m&1 == 1 {
		// m is odd and, since m != 1, m >= 3: Montgomery applies.
		return NewMontgomery(m).PowMont(base, exp)
	}
	result := uint64(1) % m
	b := base
	e := exp
	for e > 0 {
		if e&1 == 1 {
			result = MulModU64(result, b, m)
		}
		b = MulModU64(b, b, m)
		e >>= 1
	}
	return result
}

// ntheoryMontInverse returns n**-1 mod 2**64 for odd n, using Newton iteration.
// Each step x <- x*(2 - n*x) doubles the number of correct low bits; starting
// from x == 1 (correct modulo 2) six iterations are enough to reach 64 bits. All
// arithmetic wraps modulo 2**64, which is exactly the ring we invert in.
func ntheoryMontInverse(n uint64) uint64 {
	x := uint64(1)
	for i := 0; i < 6; i++ {
		x *= 2 - n*x
	}
	return x
}

// Montgomery holds the precomputed constants for Montgomery reduction with a
// fixed odd modulus. It lets repeated modular multiplications avoid division by
// working in the Montgomery domain, where reduction (REDC) uses only
// multiplications, additions and a single conditional subtraction.
type Montgomery struct {
	n    uint64 // the odd modulus
	nInv uint64 // -n**-1 mod 2**64, via Newton iteration
	r2   uint64 // 2**128 mod n, used to map values into the Montgomery domain
	one  uint64 // R mod n == 2**64 mod n, the Montgomery form of 1
}

// NewMontgomery returns a Montgomery context for the odd modulus m. It panics if
// m is even or m < 3, since Montgomery reduction requires an odd modulus and
// R == 2**64 to be coprime to it.
func NewMontgomery(m uint64) *Montgomery {
	if m < 3 || m&1 == 0 {
		panic("ntheory: NewMontgomery requires an odd modulus >= 3")
	}
	inv := ntheoryMontInverse(m) // n**-1 mod 2**64
	one := (0 - m) % m           // 2**64 mod m == R mod n
	mo := &Montgomery{
		n:    m,
		nInv: -inv, // -n**-1 mod 2**64
		one:  one,
	}
	mo.r2 = MulModU64(one, one, m) // (2**64)**2 mod m == 2**128 mod m
	return mo
}

// Modulus returns the odd modulus this context reduces by.
func (mo *Montgomery) Modulus() uint64 {
	return mo.n
}

// ntheoryRedc performs the Montgomery reduction REDC of the 128-bit value
// hi:lo, returning (hi:lo) * R**-1 mod n normalized to [0, n). It uses only
// multiply/add/subtract, no division.
func (mo *Montgomery) ntheoryRedc(hi, lo uint64) uint64 {
	m := lo * mo.nInv
	mHi, mLo := bits.Mul64(m, mo.n)
	// Add the 128-bit values hi:lo and mHi:mLo. By construction the low 64 bits
	// sum to zero modulo 2**64, so only the carry out and the high words matter.
	_, carry := bits.Add64(lo, mLo, 0)
	res, carry2 := bits.Add64(hi, mHi, carry)
	if carry2 != 0 || res >= mo.n {
		res -= mo.n
	}
	return res
}

// ToMont maps a into the Montgomery domain, returning a*R mod n. a is reduced
// modulo n first.
func (mo *Montgomery) ToMont(a uint64) uint64 {
	return mo.MulMont(a%mo.n, mo.r2)
}

// FromMont maps a Montgomery-domain value back to normal form, returning the
// Montgomery reduction a*R**-1 mod n.
func (mo *Montgomery) FromMont(a uint64) uint64 {
	return mo.ntheoryRedc(0, a)
}

// MulMont returns REDC(a*b), the Montgomery product of a and b. Both operands
// must already be in Montgomery form; the result is in Montgomery form and lies
// in [0, n).
func (mo *Montgomery) MulMont(a, b uint64) uint64 {
	hi, lo := bits.Mul64(a, b)
	return mo.ntheoryRedc(hi, lo)
}

// PowMont returns base**exp mod n. base is given in NORMAL (non-Montgomery)
// form and the result is returned in normal form; the exponentiation itself runs
// entirely in the Montgomery domain via square-and-multiply over [Montgomery.MulMont].
func (mo *Montgomery) PowMont(base, exp uint64) uint64 {
	bm := mo.ToMont(base % mo.n)
	result := mo.one
	e := exp
	for e > 0 {
		if e&1 == 1 {
			result = mo.MulMont(result, bm)
		}
		bm = mo.MulMont(bm, bm)
		e >>= 1
	}
	return mo.FromMont(result)
}

// ntheoryReciprocal2by1 returns the Moeller-Granlund reciprocal
// floor((2**128 - 1) / d) - 2**64 of a normalized divisor d (that is, d with its
// most significant bit set, 2**63 <= d < 2**64). This precomputed reciprocal
// drives the division-free 2-by-1 division ntheoryDiv2by1.
func ntheoryReciprocal2by1(d uint64) uint64 {
	// floor((2**128 - 1)/d) - 2**64 == floor((2**128 - 1 - 2**64*d)/d), whose
	// 128-bit numerator is (^d):(^0) with high word ^d < d (no overflow).
	q, _ := bits.Div64(^d, ^uint64(0), d)
	return q
}

// ntheoryDiv2by1 divides the 128-bit value u1:u0 by the normalized divisor d
// (2**63 <= d < 2**64), returning the 64-bit quotient and remainder. It requires
// u1 < d and uses the precomputed reciprocal v from ntheoryReciprocal2by1, so it
// performs no hardware division. This is Algorithm 4 of Moeller and Granlund,
// "Improved Division by Invariant Integers".
func ntheoryDiv2by1(u1, u0, d, v uint64) (q, r uint64) {
	q1, q0 := bits.Mul64(v, u1)
	var c uint64
	q0, c = bits.Add64(q0, u0, 0)
	q1, _ = bits.Add64(q1, u1, c)
	q1++
	r = u0 - q1*d
	if r > q0 {
		q1--
		r += d
	}
	if r >= d {
		q1++
		r -= d
	}
	return q1, r
}

// Barrett holds a precomputed reciprocal for repeated reduction by a fixed
// modulus. Because it does not require the modulus to be odd, it is the natural
// choice for even moduli, where Montgomery reduction does not apply. The
// modulus is normalized (shifted so its top bit is set) and reduced with a
// Moeller-Granlund reciprocal, so the hot path uses only multiply/add/subtract.
type Barrett struct {
	m     uint64 // the modulus, as supplied
	d     uint64 // normalized modulus m << shift (top bit set)
	mu    uint64 // reciprocal of d, from ntheoryReciprocal2by1 (floor(2**128/d)-like)
	shift uint   // normalization shift == bits.LeadingZeros64(m)
}

// NewBarrett returns a Barrett context for the modulus m. m must be positive;
// NewBarrett panics if m == 0.
func NewBarrett(m uint64) *Barrett {
	if m == 0 {
		panic("ntheory: NewBarrett requires a positive modulus")
	}
	shift := uint(bits.LeadingZeros64(m))
	d := m << shift
	return &Barrett{
		m:     m,
		d:     d,
		mu:    ntheoryReciprocal2by1(d),
		shift: shift,
	}
}

// Modulus returns the modulus this context reduces by.
func (b *Barrett) Modulus() uint64 {
	return b.m
}

// Reduce returns x mod m for the context's modulus m. The result lies in
// [0, m).
func (b *Barrett) Reduce(x uint64) uint64 {
	// Normalize the dividend by the same shift as the divisor: (x << shift) mod d
	// equals (x mod m) << shift. The high word x >> (64-shift) is < 2**shift <= d,
	// satisfying the u1 < d precondition of ntheoryDiv2by1.
	var u1, u0 uint64
	if b.shift == 0 {
		u1, u0 = 0, x
	} else {
		u1 = x >> (64 - b.shift)
		u0 = x << b.shift
	}
	_, r := ntheoryDiv2by1(u1, u0, b.d, b.mu)
	return r >> b.shift
}

// MulMod returns (a*c) mod m for the context's modulus m. It forms the exact
// 128-bit product with bits.Mul64 and reduces it with the precomputed
// reciprocal. The result lies in [0, m).
func (b *Barrett) MulMod(a, c uint64) uint64 {
	hi, lo := bits.Mul64(a, c)
	return b.ntheoryReduce128(hi, lo)
}

// ntheoryReduce128 returns the 128-bit value hi:lo modulo the context's
// modulus m, in [0, m). It first reduces the high word so that hi < m, then
// performs one normalized 2-by-1 division.
func (b *Barrett) ntheoryReduce128(hi, lo uint64) uint64 {
	hi = b.Reduce(hi) // now hi < m, so (hi:lo) < m*2**64
	var u1, u0 uint64
	if b.shift == 0 {
		u1, u0 = hi, lo
	} else {
		u1 = (hi << b.shift) | (lo >> (64 - b.shift))
		u0 = lo << b.shift
	}
	_, r := ntheoryDiv2by1(u1, u0, b.d, b.mu)
	return r >> b.shift
}
