package bigfloat

import (
	"fmt"
	"math/big"
)

// Exp returns e**x to prec bits.
func Exp(x *big.Float, prec uint) *big.Float {
	wp := expGuard(x, prec)
	return roundTo(prec, bfExp(clone(wp, x), wp))
}

// Expm1 returns e**x - 1 to prec bits, accurate even when x is small.
func Expm1(x *big.Float, prec uint) *big.Float {
	wp := expGuard(x, prec) + 8
	e := bfExp(clone(wp, x), wp)
	e.Sub(e, oneF(wp))
	return roundTo(prec, e)
}

// Exp2 returns 2**x to prec bits.
func Exp2(x *big.Float, prec uint) *big.Float {
	wp := expGuard(x, prec)
	arg := new(big.Float).SetPrec(wp).Mul(clone(wp, x), bfLn2(wp))
	return roundTo(prec, bfExp(arg, wp))
}

// Exp10 returns 10**x to prec bits.
func Exp10(x *big.Float, prec uint) *big.Float {
	wp := expGuard(x, prec)
	arg := new(big.Float).SetPrec(wp).Mul(clone(wp, x), bfLn(intF(wp, 10), wp))
	return roundTo(prec, bfExp(arg, wp))
}

// Log returns the natural logarithm of x to prec bits. It returns ErrDomain if
// x <= 0.
func Log(x *big.Float, prec uint) (*big.Float, error) {
	if x.Sign() <= 0 {
		return nil, fmt.Errorf("%w: Log of non-positive %s", ErrDomain, String(x))
	}
	wp := working(prec)
	return roundTo(prec, bfLn(clone(wp, x), wp)), nil
}

// Log2 returns the base-2 logarithm of x to prec bits. It returns ErrDomain if
// x <= 0.
func Log2(x *big.Float, prec uint) (*big.Float, error) {
	if x.Sign() <= 0 {
		return nil, fmt.Errorf("%w: Log2 of non-positive %s", ErrDomain, String(x))
	}
	wp := working(prec)
	r := new(big.Float).SetPrec(wp).Quo(bfLn(clone(wp, x), wp), bfLn2(wp))
	return roundTo(prec, r), nil
}

// Log10 returns the base-10 logarithm of x to prec bits. It returns ErrDomain
// if x <= 0.
func Log10(x *big.Float, prec uint) (*big.Float, error) {
	if x.Sign() <= 0 {
		return nil, fmt.Errorf("%w: Log10 of non-positive %s", ErrDomain, String(x))
	}
	wp := working(prec)
	r := new(big.Float).SetPrec(wp).Quo(bfLn(clone(wp, x), wp), bfLn(intF(wp, 10), wp))
	return roundTo(prec, r), nil
}

// Log1p returns log(1+x) to prec bits, accurate for small x. It returns
// ErrDomain if x <= -1.
func Log1p(x *big.Float, prec uint) (*big.Float, error) {
	wp := working(prec) + 8
	arg := new(big.Float).SetPrec(wp).Add(oneF(wp), clone(wp, x))
	if arg.Sign() <= 0 {
		return nil, fmt.Errorf("%w: Log1p of %s", ErrDomain, String(x))
	}
	return roundTo(prec, bfLn(arg, wp)), nil
}

// LogBase returns the logarithm of x in the given base to prec bits. It returns
// ErrDomain if x <= 0 or base <= 0 or base == 1.
func LogBase(x, base *big.Float, prec uint) (*big.Float, error) {
	if x.Sign() <= 0 || base.Sign() <= 0 || base.Cmp(oneF(base.Prec())) == 0 {
		return nil, fmt.Errorf("%w: LogBase(%s, %s)", ErrDomain, String(x), String(base))
	}
	wp := working(prec)
	r := new(big.Float).SetPrec(wp).Quo(bfLn(clone(wp, x), wp), bfLn(clone(wp, base), wp))
	return roundTo(prec, r), nil
}

// Sqrt returns the square root of x to prec bits. It returns ErrNegative if
// x < 0.
func Sqrt(x *big.Float, prec uint) (*big.Float, error) {
	if x.Sign() < 0 {
		return nil, fmt.Errorf("%w: Sqrt of %s", ErrNegative, String(x))
	}
	return new(big.Float).SetPrec(prec).Sqrt(x), nil
}

// Rsqrt returns the reciprocal square root 1/sqrt(x) to prec bits. It returns
// ErrDomain if x <= 0.
func Rsqrt(x *big.Float, prec uint) (*big.Float, error) {
	if x.Sign() <= 0 {
		return nil, fmt.Errorf("%w: Rsqrt of %s", ErrDomain, String(x))
	}
	wp := working(prec)
	s := new(big.Float).SetPrec(wp).Sqrt(x)
	return roundTo(prec, new(big.Float).SetPrec(wp).Quo(oneF(wp), s)), nil
}

// Cbrt returns the real cube root of x to prec bits, for any real x (negative
// arguments yield negative roots).
func Cbrt(x *big.Float, prec uint) *big.Float {
	if x.Sign() == 0 {
		return newF(prec)
	}
	wp := working(prec)
	ax := new(big.Float).SetPrec(wp).Abs(x)
	l := bfLn(ax, wp)
	l.Quo(l, intF(wp, 3))
	r := bfExp(l, wp)
	if x.Sign() < 0 {
		r.Neg(r)
	}
	return roundTo(prec, r)
}

// Root returns the real n-th root of x to prec bits. It returns ErrDomain if
// n == 0, or if x < 0 with n even.
func Root(x *big.Float, n int, prec uint) (*big.Float, error) {
	if n == 0 {
		return nil, fmt.Errorf("%w: 0th root", ErrDomain)
	}
	if x.Sign() == 0 {
		return newF(prec), nil
	}
	neg := x.Sign() < 0
	if neg && n%2 == 0 {
		return nil, fmt.Errorf("%w: even root of negative %s", ErrDomain, String(x))
	}
	wp := working(prec)
	ax := new(big.Float).SetPrec(wp).Abs(x)
	l := bfLn(ax, wp)
	l.Quo(l, intF(wp, int64(n)))
	r := bfExp(l, wp)
	if neg {
		r.Neg(r)
	}
	return roundTo(prec, r), nil
}

// Pow returns x**y to prec bits, defined as exp(y*log(x)) for x > 0. Special
// cases: x**0 == 1 for any x; 0**y == 0 for y > 0. Negative x is supported only
// when y is an integer, via PowInt. It returns ErrDomain otherwise.
func Pow(x, y *big.Float, prec uint) (*big.Float, error) {
	if y.Sign() == 0 {
		return oneF(prec), nil
	}
	if x.Sign() == 0 {
		if y.Sign() > 0 {
			return newF(prec), nil
		}
		return nil, fmt.Errorf("%w: 0 raised to non-positive power", ErrDomain)
	}
	if x.Sign() < 0 {
		if !y.IsInt() {
			return nil, fmt.Errorf("%w: negative base with non-integer exponent", ErrDomain)
		}
		yi, _ := y.Int(nil)
		return PowInt(x, yi.Int64(), prec), nil
	}
	wp := working(prec)
	l := bfLn(clone(wp, x), wp)
	l.Mul(l, clone(wp, y))
	// Extra guard against cancellation when the product is large.
	if e := l.MantExp(nil); e > 0 {
		wp2 := wp + uint(e)
		l = new(big.Float).SetPrec(wp2).Mul(bfLn(clone(wp2, x), wp2), clone(wp2, y))
		return roundTo(prec, bfExp(l, wp2)), nil
	}
	return roundTo(prec, bfExp(l, wp)), nil
}

// PowInt returns x**n for integer n to prec bits, via binary exponentiation.
// It is exact up to the final rounding and valid for any real x (with x != 0
// when n < 0).
func PowInt(x *big.Float, n int64, prec uint) *big.Float {
	wp := working(prec) + 8
	if n == 0 {
		return oneF(prec)
	}
	neg := n < 0
	if neg {
		n = -n
	}
	result := oneF(wp)
	base := clone(wp, x)
	for n > 0 {
		if n&1 == 1 {
			result.Mul(result, base)
		}
		n >>= 1
		if n > 0 {
			base.Mul(base, base)
		}
	}
	if neg {
		result.Quo(oneF(wp), result)
	}
	return roundTo(prec, result)
}

// Sigmoid returns the logistic function 1/(1+e**-x) to prec bits.
func Sigmoid(x *big.Float, prec uint) *big.Float {
	wp := expGuard(x, prec) + 4
	negx := new(big.Float).SetPrec(wp).Neg(clone(wp, x))
	d := bfExp(negx, wp)
	d.Add(d, oneF(wp))
	return roundTo(prec, new(big.Float).SetPrec(wp).Quo(oneF(wp), d))
}

// Logit returns log(x/(1-x)) to prec bits, the inverse of Sigmoid. It returns
// ErrDomain unless 0 < x < 1.
func Logit(x *big.Float, prec uint) (*big.Float, error) {
	wp := working(prec) + 8
	if x.Sign() <= 0 || x.Cmp(oneF(x.Prec())) >= 0 {
		return nil, fmt.Errorf("%w: Logit of %s", ErrDomain, String(x))
	}
	num := clone(wp, x)
	den := new(big.Float).SetPrec(wp).Sub(oneF(wp), clone(wp, x))
	r := num.Quo(num, den)
	return roundTo(prec, bfLn(r, wp)), nil
}

// Softplus returns log(1+e**x) to prec bits, computed in a numerically stable
// way for large |x|.
func Softplus(x *big.Float, prec uint) *big.Float {
	wp := expGuard(x, prec) + 8
	// softplus(x) = max(x,0) + log(1+e**-|x|).
	ax := new(big.Float).SetPrec(wp).Abs(clone(wp, x))
	ax.Neg(ax)
	e := bfExp(ax, wp)
	e.Add(e, oneF(wp))
	r := bfLn(e, wp)
	if x.Sign() > 0 {
		r.Add(r, clone(wp, x))
	}
	return roundTo(prec, r)
}

// LogAddExp returns log(e**x + e**y) to prec bits, computed stably.
func LogAddExp(x, y *big.Float, prec uint) *big.Float {
	wp := working(prec) + 8
	hi, lo := clone(wp, x), clone(wp, y)
	if hi.Cmp(lo) < 0 {
		hi, lo = lo, hi
	}
	d := new(big.Float).SetPrec(wp).Sub(lo, hi) // <= 0
	e := bfExp(d, wp)
	e.Add(e, oneF(wp))
	r := bfLn(e, wp)
	r.Add(r, hi)
	return roundTo(prec, r)
}

// Sinc returns the unnormalised cardinal sine sin(x)/x to prec bits, with
// Sinc(0) == 1.
func Sinc(x *big.Float, prec uint) *big.Float {
	if x.Sign() == 0 {
		return oneF(prec)
	}
	wp := expGuard(x, prec) + 4
	halfPi := mulPow2(bfPi(wp), -1)
	sin, _ := bfSinCos(clone(wp, x), halfPi, wp)
	r := new(big.Float).SetPrec(wp).Quo(sin, clone(wp, x))
	return roundTo(prec, r)
}
