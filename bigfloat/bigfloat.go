package bigfloat

import (
	"errors"
	"fmt"
	"math"
	"math/big"
	"strings"
)

// Domain errors returned by the partial functions in this package.
var (
	// ErrDomain reports that an argument lies outside a function's domain.
	ErrDomain = errors.New("bigfloat: argument out of domain")
	// ErrNegative reports that a function received a negative argument where a
	// non-negative one was required.
	ErrNegative = errors.New("bigfloat: negative argument")
	// ErrPole reports that a function was evaluated at one of its poles.
	ErrPole = errors.New("bigfloat: function pole")
	// ErrParse reports that a string could not be parsed as a big.Float.
	ErrParse = errors.New("bigfloat: cannot parse number")
)

const defaultGuard = 64

// working returns the internal working precision used to evaluate a result at
// the requested precision prec. A guard band is always added; prec == 0 is
// treated as double precision.
func working(prec uint) uint {
	if prec == 0 {
		prec = 53
	}
	return prec + defaultGuard
}

// newF allocates a new zero *big.Float with the given precision.
func newF(prec uint) *big.Float { return new(big.Float).SetPrec(prec) }

// intF returns n as a *big.Float with the given precision.
func intF(prec uint, n int64) *big.Float { return new(big.Float).SetPrec(prec).SetInt64(n) }

// oneF returns 1 as a *big.Float with the given precision.
func oneF(prec uint) *big.Float { return intF(prec, 1) }

// clone copies x into a fresh *big.Float of the given precision.
func clone(prec uint, x *big.Float) *big.Float { return new(big.Float).SetPrec(prec).Set(x) }

// roundTo returns x rounded to prec bits in a new value.
func roundTo(prec uint, x *big.Float) *big.Float { return new(big.Float).SetPrec(prec).Set(x) }

// pow2 returns 2**n as a *big.Float with the given precision.
func pow2(prec uint, n int) *big.Float {
	return new(big.Float).SetPrec(prec).SetMantExp(oneF(prec), n)
}

// mulPow2 returns x * 2**n. SetMantExp(x, n) computes x*2**n directly because it
// treats its first argument as an unnormalised mantissa.
func mulPow2(x *big.Float, n int) *big.Float {
	return new(big.Float).SetPrec(x.Prec()).SetMantExp(x, n)
}

// tiny reports whether term is negligible relative to sum at prec bits, i.e.
// adding term no longer affects the leading prec bits of sum.
func tiny(term, sum *big.Float, prec uint) bool {
	if term.Sign() == 0 {
		return true
	}
	if sum.Sign() == 0 {
		return false
	}
	return term.MantExp(nil) < sum.MantExp(nil)-int(prec)-4
}

// seriesIters returns a reduction/step count on the order of sqrt(prec), used
// to balance argument-halving passes against series length.
func seriesIters(prec uint) int { return int(math.Sqrt(float64(prec))) + 2 }

// nearestInt rounds x to the nearest integer (ties away from zero).
func nearestInt(x *big.Float) *big.Int {
	f := new(big.Float).SetPrec(x.Prec() + 4).Set(x)
	half := big.NewFloat(0.5)
	if f.Signbit() {
		f.Sub(f, half)
	} else {
		f.Add(f, half)
	}
	i, _ := f.Int(nil)
	return i
}

// -----------------------------------------------------------------------------
// Construction and conversion.
// -----------------------------------------------------------------------------

// New returns a new zero value with the given precision in bits.
func New(prec uint) *big.Float { return newF(prec) }

// Zero returns 0 with the given precision.
func Zero(prec uint) *big.Float { return newF(prec) }

// One returns 1 with the given precision.
func One(prec uint) *big.Float { return intF(prec, 1) }

// Two returns 2 with the given precision.
func Two(prec uint) *big.Float { return intF(prec, 2) }

// Half returns 1/2 with the given precision.
func Half(prec uint) *big.Float { return new(big.Float).SetPrec(prec).SetFloat64(0.5) }

// FromInt returns the integer n as a *big.Float with the given precision.
func FromInt(n int64, prec uint) *big.Float { return intF(prec, n) }

// FromBig returns the big.Int n as a *big.Float with the given precision.
func FromBig(n *big.Int, prec uint) *big.Float { return new(big.Float).SetPrec(prec).SetInt(n) }

// FromRat returns the rational q as a *big.Float with the given precision.
func FromRat(q *big.Rat, prec uint) *big.Float { return new(big.Float).SetPrec(prec).SetRat(q) }

// FromFloat64 returns the float64 x as a *big.Float with the given precision.
func FromFloat64(x float64, prec uint) *big.Float {
	return new(big.Float).SetPrec(prec).SetFloat64(x)
}

// Parse parses s (in any base-10 float syntax accepted by big.Float.SetString)
// into a *big.Float with the given precision. It returns ErrParse on failure.
func Parse(s string, prec uint) (*big.Float, error) {
	z := new(big.Float).SetPrec(prec)
	if _, ok := z.SetString(strings.TrimSpace(s)); !ok {
		return nil, fmt.Errorf("%w: %q", ErrParse, s)
	}
	return z, nil
}

// MustParse is like Parse but panics on a parse error. It is intended for
// package-level initialisers and tests with known-good literals.
func MustParse(s string, prec uint) *big.Float {
	z, err := Parse(s, prec)
	if err != nil {
		panic(err)
	}
	return z
}

// Float64 returns the value of x as a float64, together with the accuracy of
// the conversion.
func Float64(x *big.Float) (float64, big.Accuracy) { return x.Float64() }

// String formats x using big.Float's default text representation.
func String(x *big.Float) string { return x.Text('g', -1) }

// Text formats x with the given format byte and number of digits, matching
// big.Float.Text.
func Text(x *big.Float, format byte, digits int) string { return x.Text(format, digits) }

// -----------------------------------------------------------------------------
// Basic arithmetic (each rounds to the requested precision).
// -----------------------------------------------------------------------------

// Add returns x + y rounded to prec bits.
func Add(x, y *big.Float, prec uint) *big.Float { return new(big.Float).SetPrec(prec).Add(x, y) }

// Sub returns x - y rounded to prec bits.
func Sub(x, y *big.Float, prec uint) *big.Float { return new(big.Float).SetPrec(prec).Sub(x, y) }

// Mul returns x * y rounded to prec bits.
func Mul(x, y *big.Float, prec uint) *big.Float { return new(big.Float).SetPrec(prec).Mul(x, y) }

// Div returns x / y rounded to prec bits.
func Div(x, y *big.Float, prec uint) *big.Float { return new(big.Float).SetPrec(prec).Quo(x, y) }

// Recip returns 1 / x rounded to prec bits.
func Recip(x *big.Float, prec uint) *big.Float {
	return new(big.Float).SetPrec(prec).Quo(oneF(prec), x)
}

// Neg returns -x rounded to prec bits.
func Neg(x *big.Float, prec uint) *big.Float { return new(big.Float).SetPrec(prec).Neg(x) }

// Abs returns |x| rounded to prec bits.
func Abs(x *big.Float, prec uint) *big.Float { return new(big.Float).SetPrec(prec).Abs(x) }

// Square returns x*x rounded to prec bits.
func Square(x *big.Float, prec uint) *big.Float { return new(big.Float).SetPrec(prec).Mul(x, x) }

// Cube returns x*x*x rounded to prec bits.
func Cube(x *big.Float, prec uint) *big.Float {
	z := new(big.Float).SetPrec(prec).Mul(x, x)
	return z.Mul(z, x)
}

// FMA returns x*y + z rounded to prec bits (computed at higher internal
// precision to reduce intermediate rounding).
func FMA(x, y, z *big.Float, prec uint) *big.Float {
	wp := working(prec)
	t := new(big.Float).SetPrec(wp).Mul(x, y)
	t.Add(t, z)
	return roundTo(prec, t)
}

// Copysign returns a value with the magnitude of x and the sign of y.
func Copysign(x, y *big.Float, prec uint) *big.Float {
	z := new(big.Float).SetPrec(prec).Abs(x)
	if y.Signbit() {
		z.Neg(z)
	}
	return z
}

// Dim returns max(x-y, 0) rounded to prec bits.
func Dim(x, y *big.Float, prec uint) *big.Float {
	z := new(big.Float).SetPrec(prec).Sub(x, y)
	if z.Sign() < 0 {
		z.SetInt64(0)
	}
	return z
}

// Hypot returns sqrt(x*x + y*y), computed without undue overflow, to prec bits.
func Hypot(x, y *big.Float, prec uint) *big.Float {
	wp := working(prec)
	ax := new(big.Float).SetPrec(wp).Abs(x)
	ay := new(big.Float).SetPrec(wp).Abs(y)
	if ax.Sign() == 0 {
		return roundTo(prec, ay)
	}
	if ay.Sign() == 0 {
		return roundTo(prec, ax)
	}
	// Scale by the larger to avoid overflow/underflow.
	big1, small := ax, ay
	if ax.Cmp(ay) < 0 {
		big1, small = ay, ax
	}
	r := new(big.Float).SetPrec(wp).Quo(small, big1)
	r.Mul(r, r)
	r.Add(r, oneF(wp))
	r.Sqrt(r)
	r.Mul(r, big1)
	return roundTo(prec, r)
}

// -----------------------------------------------------------------------------
// Comparison and predicates.
// -----------------------------------------------------------------------------

// Cmp compares x and y and returns -1, 0, or +1.
func Cmp(x, y *big.Float) int { return x.Cmp(y) }

// Sign returns -1, 0, or +1 according to the sign of x.
func Sign(x *big.Float) int { return x.Sign() }

// Signbit reports whether x is negative or negative zero.
func Signbit(x *big.Float) bool { return x.Signbit() }

// IsZero reports whether x is zero.
func IsZero(x *big.Float) bool { return x.Sign() == 0 }

// IsInf reports whether x is an infinity.
func IsInf(x *big.Float) bool { return x.IsInf() }

// IsInteger reports whether x has no fractional part.
func IsInteger(x *big.Float) bool { return x.IsInt() }

// Equal reports whether x and y are exactly equal.
func Equal(x, y *big.Float) bool { return x.Cmp(y) == 0 }

// AlmostEqual reports whether |x-y| <= tol.
func AlmostEqual(x, y, tol *big.Float) bool {
	d := new(big.Float).SetPrec(maxPrec(x, y)).Sub(x, y)
	d.Abs(d)
	return d.Cmp(tol) <= 0
}

func maxPrec(x, y *big.Float) uint {
	if x.Prec() > y.Prec() {
		return x.Prec()
	}
	return y.Prec()
}

// Min returns the smaller of x and y (rounded to prec bits).
func Min(x, y *big.Float, prec uint) *big.Float {
	if x.Cmp(y) <= 0 {
		return roundTo(prec, x)
	}
	return roundTo(prec, y)
}

// Max returns the larger of x and y (rounded to prec bits).
func Max(x, y *big.Float, prec uint) *big.Float {
	if x.Cmp(y) >= 0 {
		return roundTo(prec, x)
	}
	return roundTo(prec, y)
}

// Clamp returns x confined to the interval [lo, hi], rounded to prec bits.
func Clamp(x, lo, hi *big.Float, prec uint) *big.Float {
	if x.Cmp(lo) < 0 {
		return roundTo(prec, lo)
	}
	if x.Cmp(hi) > 0 {
		return roundTo(prec, hi)
	}
	return roundTo(prec, x)
}

// -----------------------------------------------------------------------------
// Rounding to integral values.
// -----------------------------------------------------------------------------

// Trunc returns the integer part of x (rounding toward zero), to prec bits.
func Trunc(x *big.Float, prec uint) *big.Float {
	i, _ := x.Int(nil)
	return new(big.Float).SetPrec(prec).SetInt(i)
}

// Floor returns the greatest integer <= x, to prec bits.
func Floor(x *big.Float, prec uint) *big.Float {
	i, acc := x.Int(nil)
	if x.Sign() < 0 && acc != big.Exact {
		i.Sub(i, big.NewInt(1))
	}
	return new(big.Float).SetPrec(prec).SetInt(i)
}

// Ceil returns the least integer >= x, to prec bits.
func Ceil(x *big.Float, prec uint) *big.Float {
	i, acc := x.Int(nil)
	if x.Sign() > 0 && acc != big.Exact {
		i.Add(i, big.NewInt(1))
	}
	return new(big.Float).SetPrec(prec).SetInt(i)
}

// Round returns x rounded to the nearest integer, ties away from zero.
func Round(x *big.Float, prec uint) *big.Float {
	return new(big.Float).SetPrec(prec).SetInt(nearestInt(x))
}

// Frac returns the fractional part x - Trunc(x), to prec bits.
func Frac(x *big.Float, prec uint) *big.Float {
	wp := working(prec)
	i, _ := x.Int(nil)
	t := new(big.Float).SetPrec(wp).SetInt(i)
	return roundTo(prec, t.Sub(x, t))
}

// Ldexp returns x * 2**n, to prec bits.
func Ldexp(x *big.Float, n int, prec uint) *big.Float {
	return new(big.Float).SetPrec(prec).SetMantExp(x, n)
}

// Frexp breaks x into a normalised fraction f in [0.5,1) and exponent e such
// that x == f * 2**e. The fraction is returned to prec bits.
func Frexp(x *big.Float, prec uint) (frac *big.Float, exp int) {
	m := new(big.Float).SetPrec(prec)
	exp = x.MantExp(m)
	return m, exp
}

// Mod returns the floating-point remainder x - y*Trunc(x/y), to prec bits.
func Mod(x, y *big.Float, prec uint) *big.Float {
	wp := working(prec) + uint(absInt(x.MantExp(nil)))
	q := new(big.Float).SetPrec(wp).Quo(x, y)
	i, _ := q.Int(nil)
	t := new(big.Float).SetPrec(wp).SetInt(i)
	t.Mul(t, y)
	t.Sub(x, t)
	return roundTo(prec, t)
}

// Remainder returns the IEEE 754 remainder x - y*Round(x/y), to prec bits.
func Remainder(x, y *big.Float, prec uint) *big.Float {
	wp := working(prec) + uint(absInt(x.MantExp(nil)))
	q := new(big.Float).SetPrec(wp).Quo(x, y)
	n := nearestInt(q)
	t := new(big.Float).SetPrec(wp).SetInt(n)
	t.Mul(t, y)
	t.Sub(x, t)
	return roundTo(prec, t)
}

func absInt(n int) int {
	if n < 0 {
		return -n
	}
	return n
}
