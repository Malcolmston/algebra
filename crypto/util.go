package crypto

import (
	"math/big"
	"math/rand"
)

// Shared constants used throughout the package. They are never mutated; callers
// that need a modifiable copy must allocate their own via new(big.Int).Set.
var (
	cryptoZero  = big.NewInt(0)
	cryptoOne   = big.NewInt(1)
	cryptoTwo   = big.NewInt(2)
	cryptoThree = big.NewInt(3)
)

// cryptoRandBelow returns a uniformly distributed integer in the half-open
// interval [0, n) using rng as the source of randomness. It panics if n is not
// positive or if rng is nil. The result is produced by rejection sampling over
// whole bytes, so the distribution is exactly uniform.
func cryptoRandBelow(rng *rand.Rand, n *big.Int) *big.Int {
	if rng == nil {
		panic("crypto: nil random source")
	}
	if n.Sign() <= 0 {
		panic("crypto: cryptoRandBelow requires n > 0")
	}
	bits := n.BitLen()
	byteLen := (bits + 7) / 8
	excess := uint(byteLen*8 - bits)
	buf := make([]byte, byteLen)
	x := new(big.Int)
	for {
		for i := range buf {
			buf[i] = byte(rng.Intn(256))
		}
		// Clear the surplus high bits so the candidate rarely exceeds n; this
		// keeps the expected number of rejections below two.
		if excess > 0 {
			buf[0] &= byte(0xff >> excess)
		}
		x.SetBytes(buf)
		if x.Cmp(n) < 0 {
			return x
		}
	}
}

// cryptoRandRange returns a uniformly distributed integer in the closed
// interval [lo, hi]. It panics if lo > hi.
func cryptoRandRange(rng *rand.Rand, lo, hi *big.Int) *big.Int {
	if lo.Cmp(hi) > 0 {
		panic("crypto: cryptoRandRange requires lo <= hi")
	}
	span := new(big.Int).Sub(hi, lo)
	span.Add(span, cryptoOne)
	r := cryptoRandBelow(rng, span)
	return r.Add(r, lo)
}

// cryptoRandBits returns a uniformly distributed non-negative integer whose bit
// length is at most bits (the top bit may be zero). It panics for bits <= 0.
func cryptoRandBits(rng *rand.Rand, bits int) *big.Int {
	if bits <= 0 {
		panic("crypto: cryptoRandBits requires bits > 0")
	}
	n := new(big.Int).Lsh(cryptoOne, uint(bits))
	return cryptoRandBelow(rng, n)
}
