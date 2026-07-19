package groebner

import (
	"math/big"
)

// NewRat returns the rational number a/b as a *big.Rat. It panics only if b is
// zero, matching the behaviour of big.NewRat.
func NewRat(a, b int64) *big.Rat {
	return big.NewRat(a, b)
}

// RatFromInt returns the integer a as a *big.Rat.
func RatFromInt(a int64) *big.Rat {
	return big.NewRat(a, 1)
}

// RatFromString parses a rational number from its decimal or fractional string
// representation (for example "3/4", "-2", or "1.5"). The boolean result is
// false when the string cannot be parsed.
func RatFromString(s string) (*big.Rat, bool) {
	r := new(big.Rat)
	_, ok := r.SetString(s)
	return r, ok
}

// RatIsZero reports whether the rational r is zero.
func RatIsZero(r *big.Rat) bool { return r.Sign() == 0 }

// RatIsOne reports whether the rational r equals 1.
func RatIsOne(r *big.Rat) bool { return r.Cmp(bigOne) == 0 }

// CloneRat returns an independent copy of r.
func CloneRat(r *big.Rat) *big.Rat { return new(big.Rat).Set(r) }

var (
	bigZero = big.NewRat(0, 1)
	bigOne  = big.NewRat(1, 1)
)

func ratAdd(a, b *big.Rat) *big.Rat { return new(big.Rat).Add(a, b) }
func ratSub(a, b *big.Rat) *big.Rat { return new(big.Rat).Sub(a, b) }
func ratMul(a, b *big.Rat) *big.Rat { return new(big.Rat).Mul(a, b) }
func ratDiv(a, b *big.Rat) *big.Rat { return new(big.Rat).Quo(a, b) }
func ratNeg(a *big.Rat) *big.Rat    { return new(big.Rat).Neg(a) }

// RatToFloat returns the nearest float64 to the rational r.
func RatToFloat(r *big.Rat) float64 {
	f, _ := r.Float64()
	return f
}
