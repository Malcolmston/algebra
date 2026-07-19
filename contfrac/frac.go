package contfrac

import (
	"fmt"
	"math"
	"math/big"
)

// Frac is an exact rational number with an int64 numerator and denominator.
// The zero value Frac{} is 0/0 and is not a valid number; use [NewFrac] or the
// package constructors to build values. Most methods assume a non-zero
// denominator and, for canonical results, a reduced representation with a
// positive denominator (which [NewFrac] and [Frac.Reduce] guarantee).
type Frac struct {
	Num int64 // numerator
	Den int64 // denominator
}

// NewFrac returns the fraction n/d reduced to lowest terms with a positive
// denominator. It panics if d == 0.
func NewFrac(n, d int64) Frac {
	if d == 0 {
		panic("contfrac: NewFrac requires a non-zero denominator")
	}
	p, q := ReduceFraction(n, d)
	return Frac{p, q}
}

// Reduce returns an equivalent fraction in lowest terms with a positive
// denominator. It panics if the denominator is zero.
func (f Frac) Reduce() Frac {
	return NewFrac(f.Num, f.Den)
}

// Float returns the value of f as a float64.
func (f Frac) Float() float64 {
	return float64(f.Num) / float64(f.Den)
}

// Rat returns f as a *big.Rat.
func (f Frac) Rat() *big.Rat {
	return big.NewRat(f.Num, f.Den)
}

// Add returns the reduced sum f + g.
func (f Frac) Add(g Frac) Frac {
	return NewFrac(f.Num*g.Den+g.Num*f.Den, f.Den*g.Den)
}

// Sub returns the reduced difference f - g.
func (f Frac) Sub(g Frac) Frac {
	return NewFrac(f.Num*g.Den-g.Num*f.Den, f.Den*g.Den)
}

// Mul returns the reduced product f * g.
func (f Frac) Mul(g Frac) Frac {
	return NewFrac(f.Num*g.Num, f.Den*g.Den)
}

// Div returns the reduced quotient f / g. It panics if g is zero.
func (f Frac) Div(g Frac) Frac {
	if g.Num == 0 {
		panic("contfrac: division by zero fraction")
	}
	return NewFrac(f.Num*g.Den, f.Den*g.Num)
}

// Neg returns -f.
func (f Frac) Neg() Frac {
	return Frac{-f.Num, f.Den}
}

// Inv returns the reciprocal 1/f. It panics if f is zero.
func (f Frac) Inv() Frac {
	if f.Num == 0 {
		panic("contfrac: reciprocal of zero")
	}
	return NewFrac(f.Den, f.Num)
}

// Abs returns |f|.
func (f Frac) Abs() Frac {
	if f.Num < 0 {
		return Frac{-f.Num, f.Den}
	}
	return f
}

// Sign returns -1, 0 or +1 according to the sign of f (assuming a positive
// denominator, as produced by the constructors).
func (f Frac) Sign() int {
	switch {
	case f.Num < 0:
		return -1
	case f.Num > 0:
		return 1
	default:
		return 0
	}
}

// Cmp compares f and g, returning -1 if f < g, 0 if they are equal and +1 if
// f > g. The comparison is exact for inputs whose cross products fit in int64.
func (f Frac) Cmp(g Frac) int {
	lhs := f.Num * g.Den
	rhs := g.Num * f.Den
	// Normalise for denominator signs (constructors keep them positive, but be
	// robust to hand-built values).
	if (f.Den < 0) != (g.Den < 0) {
		lhs, rhs = rhs, lhs
	}
	switch {
	case lhs < rhs:
		return -1
	case lhs > rhs:
		return 1
	default:
		return 0
	}
}

// Equal reports whether f and g represent the same rational number.
func (f Frac) Equal(g Frac) bool {
	return f.Cmp(g) == 0
}

// Less reports whether f < g.
func (f Frac) Less(g Frac) bool {
	return f.Cmp(g) < 0
}

// IsZero reports whether f == 0.
func (f Frac) IsZero() bool {
	return f.Num == 0
}

// IsInteger reports whether f is an integer.
func (f Frac) IsInteger() bool {
	r := f.Reduce()
	return r.Den == 1
}

// Floor returns the greatest integer not exceeding f.
func (f Frac) Floor() int64 {
	r := f.Reduce()
	return floorDiv(r.Num, r.Den)
}

// Ceil returns the least integer not less than f.
func (f Frac) Ceil() int64 {
	r := f.Reduce()
	return -floorDiv(-r.Num, r.Den)
}

// Round returns f rounded to the nearest integer, rounding halves away from
// zero.
func (f Frac) Round() int64 {
	r := f.Reduce()
	if r.Num >= 0 {
		return floorDiv(2*r.Num+r.Den, 2*r.Den)
	}
	return -((-2*r.Num + r.Den) / (2 * r.Den))
}

// Pow returns f raised to the integer power n (which may be negative). It
// panics for a negative power of zero.
func (f Frac) Pow(n int) Frac {
	if n < 0 {
		return f.Inv().Pow(-n)
	}
	result := Frac{1, 1}
	base := f
	for n > 0 {
		if n&1 == 1 {
			result = result.Mul(base)
		}
		base = base.Mul(base)
		n >>= 1
	}
	return result.Reduce()
}

// Mediant returns the mediant (f.Num+g.Num)/(f.Den+g.Den) of f and g. The
// mediant of two reduced fractions lies strictly between them.
func (f Frac) Mediant(g Frac) Frac {
	return Frac{f.Num + g.Num, f.Den + g.Den}
}

// String renders f as "p/q", or just "p" when the denominator is 1.
func (f Frac) String() string {
	r := f.Reduce()
	if r.Den == 1 {
		return fmt.Sprintf("%d", r.Num)
	}
	return fmt.Sprintf("%d/%d", r.Num, r.Den)
}

// FracFromFloat returns the exact Frac whose value equals x, provided x is a
// dyadic rational representable in an int64 numerator and denominator, together
// with ok. When x cannot be represented exactly (it needs more than 62 bits or
// is non-finite), ok is false. For a best approximation with a bounded
// denominator use [BestApproximationFrac] instead.
func FracFromFloat(x float64) (Frac, bool) {
	if math.IsNaN(x) || math.IsInf(x, 0) {
		return Frac{}, false
	}
	if x == math.Trunc(x) && math.Abs(x) < 9.2e18 {
		return Frac{int64(x), 1}, true
	}
	den := int64(1)
	y := x
	for i := 0; i < 62; i++ {
		if y == math.Trunc(y) {
			if math.Abs(y) >= 9.2e18 {
				return Frac{}, false
			}
			return NewFrac(int64(y), den), true
		}
		y *= 2
		den <<= 1
		if den <= 0 {
			return Frac{}, false
		}
	}
	return Frac{}, false
}
