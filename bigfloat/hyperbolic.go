package bigfloat

import (
	"fmt"
	"math/big"
)

// sinhCosh returns sinh(x) and cosh(x) at prec bits from a single pair of
// exponentials.
func sinhCosh(x *big.Float, prec uint) (sinh, cosh *big.Float) {
	ex := bfExp(clone(prec, x), prec)
	enx := new(big.Float).SetPrec(prec).Quo(oneF(prec), ex)
	sinh = new(big.Float).SetPrec(prec).Sub(ex, enx)
	sinh = mulPow2(sinh, -1)
	cosh = new(big.Float).SetPrec(prec).Add(ex, enx)
	cosh = mulPow2(cosh, -1)
	return sinh, cosh
}

// Sinh returns the hyperbolic sine of x to prec bits.
func Sinh(x *big.Float, prec uint) *big.Float {
	wp := expGuard(x, prec) + 8
	s, _ := sinhCosh(x, wp)
	return roundTo(prec, s)
}

// Cosh returns the hyperbolic cosine of x to prec bits.
func Cosh(x *big.Float, prec uint) *big.Float {
	wp := expGuard(x, prec) + 8
	_, c := sinhCosh(x, wp)
	return roundTo(prec, c)
}

// Tanh returns the hyperbolic tangent of x to prec bits.
func Tanh(x *big.Float, prec uint) *big.Float {
	wp := expGuard(x, prec) + 8
	s, c := sinhCosh(x, wp)
	return roundTo(prec, new(big.Float).SetPrec(wp).Quo(s, c))
}

// Coth returns the hyperbolic cotangent of x to prec bits. It returns ErrPole
// if x == 0.
func Coth(x *big.Float, prec uint) (*big.Float, error) {
	if x.Sign() == 0 {
		return nil, fmt.Errorf("%w: Coth at 0", ErrPole)
	}
	wp := expGuard(x, prec) + 8
	s, c := sinhCosh(x, wp)
	return roundTo(prec, new(big.Float).SetPrec(wp).Quo(c, s)), nil
}

// Sech returns the hyperbolic secant 1/cosh(x) to prec bits.
func Sech(x *big.Float, prec uint) *big.Float {
	wp := expGuard(x, prec) + 8
	_, c := sinhCosh(x, wp)
	return roundTo(prec, new(big.Float).SetPrec(wp).Quo(oneF(wp), c))
}

// Csch returns the hyperbolic cosecant 1/sinh(x) to prec bits. It returns
// ErrPole if x == 0.
func Csch(x *big.Float, prec uint) (*big.Float, error) {
	if x.Sign() == 0 {
		return nil, fmt.Errorf("%w: Csch at 0", ErrPole)
	}
	wp := expGuard(x, prec) + 8
	s, _ := sinhCosh(x, wp)
	return roundTo(prec, new(big.Float).SetPrec(wp).Quo(oneF(wp), s)), nil
}

// Asinh returns the inverse hyperbolic sine log(x + sqrt(x^2+1)) to prec bits,
// for any real x.
func Asinh(x *big.Float, prec uint) *big.Float {
	wp := working(prec) + uint(absInt(x.MantExp(nil))) + 8
	ax := new(big.Float).SetPrec(wp).Abs(clone(wp, x))
	d := new(big.Float).SetPrec(wp).Mul(ax, ax)
	d.Add(d, oneF(wp))
	d.Sqrt(d)
	d.Add(d, ax)
	r := bfLn(d, wp)
	if x.Sign() < 0 {
		r.Neg(r)
	}
	return roundTo(prec, r)
}

// Acosh returns the inverse hyperbolic cosine log(x + sqrt(x^2-1)) to prec bits.
// It returns ErrDomain if x < 1.
func Acosh(x *big.Float, prec uint) (*big.Float, error) {
	if x.Cmp(oneF(x.Prec())) < 0 {
		return nil, fmt.Errorf("%w: Acosh of %s", ErrDomain, String(x))
	}
	wp := working(prec) + uint(absInt(x.MantExp(nil))) + 8
	d := new(big.Float).SetPrec(wp).Mul(clone(wp, x), clone(wp, x))
	d.Sub(d, oneF(wp))
	d.Sqrt(d)
	d.Add(d, clone(wp, x))
	return roundTo(prec, bfLn(d, wp)), nil
}

// Atanh returns the inverse hyperbolic tangent (1/2)*log((1+x)/(1-x)) to prec
// bits. It returns ErrDomain unless -1 < x < 1.
func Atanh(x *big.Float, prec uint) (*big.Float, error) {
	if new(big.Float).Abs(x).Cmp(oneF(x.Prec())) >= 0 {
		return nil, fmt.Errorf("%w: Atanh of %s", ErrDomain, String(x))
	}
	// atanh(x) = (1/2) log((1+x)/(1-x)); the log form stays accurate as
	// |x| -> 1, where the atanh power series would converge too slowly.
	wp := working(prec) + 8
	num := new(big.Float).SetPrec(wp).Add(oneF(wp), clone(wp, x))
	den := new(big.Float).SetPrec(wp).Sub(oneF(wp), clone(wp, x))
	r := bfLn(num.Quo(num, den), wp)
	return roundTo(prec, mulPow2(r, -1)), nil
}

// Acoth returns the inverse hyperbolic cotangent (1/2)*log((x+1)/(x-1)) to prec
// bits. It returns ErrDomain unless |x| > 1.
func Acoth(x *big.Float, prec uint) (*big.Float, error) {
	if new(big.Float).Abs(x).Cmp(oneF(x.Prec())) <= 0 {
		return nil, fmt.Errorf("%w: Acoth of %s", ErrDomain, String(x))
	}
	// acoth(x) = (1/2) log((x+1)/(x-1)).
	wp := working(prec) + 8
	num := new(big.Float).SetPrec(wp).Add(clone(wp, x), oneF(wp))
	den := new(big.Float).SetPrec(wp).Sub(clone(wp, x), oneF(wp))
	r := bfLn(num.Quo(num, den), wp)
	return roundTo(prec, mulPow2(r, -1)), nil
}

// Asech returns the inverse hyperbolic secant acosh(1/x) to prec bits. It
// returns ErrDomain unless 0 < x <= 1.
func Asech(x *big.Float, prec uint) (*big.Float, error) {
	if x.Sign() <= 0 || x.Cmp(oneF(x.Prec())) > 0 {
		return nil, fmt.Errorf("%w: Asech of %s", ErrDomain, String(x))
	}
	wp := working(prec) + 8
	inv := new(big.Float).SetPrec(wp).Quo(oneF(wp), clone(wp, x))
	return Acosh(inv, prec)
}

// Acsch returns the inverse hyperbolic cosecant asinh(1/x) to prec bits. It
// returns ErrDomain if x == 0.
func Acsch(x *big.Float, prec uint) (*big.Float, error) {
	if x.Sign() == 0 {
		return nil, fmt.Errorf("%w: Acsch at 0", ErrDomain)
	}
	wp := working(prec) + 8
	inv := new(big.Float).SetPrec(wp).Quo(oneF(wp), clone(wp, x))
	return Asinh(inv, prec), nil
}
