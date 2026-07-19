package diffalgebra

import (
	"math/big"
	"strings"
)

// ---------------------------------------------------------------------------
// Small helpers around *big.Rat. These are the exact scalar arithmetic that
// the polynomial and rational-function layers are built on.
// ---------------------------------------------------------------------------

// cloneRat returns a fresh copy of r so callers never alias a shared *big.Rat.
func cloneRat(r *big.Rat) *big.Rat { return new(big.Rat).Set(r) }

// ratInt returns the rational number n/1.
func ratInt(n int64) *big.Rat { return big.NewRat(n, 1) }

// ratZero reports whether r is exactly zero.
func ratZero(r *big.Rat) bool { return r.Sign() == 0 }

// ratAdd returns a+b as a new value.
func ratAdd(a, b *big.Rat) *big.Rat { return new(big.Rat).Add(a, b) }

// ratSub returns a-b as a new value.
func ratSub(a, b *big.Rat) *big.Rat { return new(big.Rat).Sub(a, b) }

// ratMul returns a*b as a new value.
func ratMul(a, b *big.Rat) *big.Rat { return new(big.Rat).Mul(a, b) }

// ratNeg returns -a as a new value.
func ratNeg(a *big.Rat) *big.Rat { return new(big.Rat).Neg(a) }

// ratInv returns 1/a as a new value; a must be non-zero.
func ratInv(a *big.Rat) *big.Rat { return new(big.Rat).Inv(a) }

// ratDiv returns a/b as a new value; b must be non-zero.
func ratDiv(a, b *big.Rat) *big.Rat { return new(big.Rat).Quo(a, b) }

// ratPow raises base to the integer power n (n may be negative).
func ratPow(base *big.Rat, n int) *big.Rat {
	if n < 0 {
		return ratInv(ratPow(base, -n))
	}
	result := ratInt(1)
	b := cloneRat(base)
	for n > 0 {
		if n&1 == 1 {
			result = ratMul(result, b)
		}
		b = ratMul(b, b)
		n >>= 1
	}
	return result
}

// RatFromInt returns the rational number n/1.
func RatFromInt(n int64) *big.Rat { return ratInt(n) }

// RatFromFrac returns the rational number p/q. It panics only if q is zero,
// matching the behaviour of math/big.NewRat.
func RatFromFrac(p, q int64) *big.Rat { return big.NewRat(p, q) }

// RatFromFloat returns the exact rational value of the IEEE-754 float f. It
// returns nil if f is not finite.
func RatFromFloat(f float64) *big.Rat {
	r := new(big.Rat).SetFloat64(f)
	if r == nil {
		return nil
	}
	return r
}

// RatFromString parses a rational number written as "p/q", a decimal, or an
// integer. It reports whether parsing succeeded.
func RatFromString(s string) (*big.Rat, bool) {
	r, ok := new(big.Rat).SetString(strings.TrimSpace(s))
	return r, ok
}

// RatToFloat returns the nearest float64 to r.
func RatToFloat(r *big.Rat) float64 {
	f, _ := r.Float64()
	return f
}

// RatString formats r as "p/q" (or "p" when the denominator is one).
func RatString(r *big.Rat) string { return r.RatString() }
