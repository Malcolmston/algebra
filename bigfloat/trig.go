package bigfloat

import (
	"fmt"
	"math/big"
)

// halfPiAt returns pi/2 at prec bits (a small internal convenience).
func halfPiAt(prec uint) *big.Float { return mulPow2(bfPi(prec), -1) }

// Sin returns the sine of x (radians) to prec bits.
func Sin(x *big.Float, prec uint) *big.Float {
	wp := expGuard(x, prec)
	sin, _ := bfSinCos(clone(wp, x), halfPiAt(wp), wp)
	return roundTo(prec, sin)
}

// Cos returns the cosine of x (radians) to prec bits.
func Cos(x *big.Float, prec uint) *big.Float {
	wp := expGuard(x, prec)
	_, cos := bfSinCos(clone(wp, x), halfPiAt(wp), wp)
	return roundTo(prec, cos)
}

// SinCos returns both the sine and cosine of x (radians) to prec bits.
func SinCos(x *big.Float, prec uint) (sin, cos *big.Float) {
	wp := expGuard(x, prec)
	s, c := bfSinCos(clone(wp, x), halfPiAt(wp), wp)
	return roundTo(prec, s), roundTo(prec, c)
}

// Tan returns the tangent of x (radians) to prec bits.
func Tan(x *big.Float, prec uint) *big.Float {
	wp := expGuard(x, prec) + 8
	s, c := bfSinCos(clone(wp, x), halfPiAt(wp), wp)
	return roundTo(prec, new(big.Float).SetPrec(wp).Quo(s, c))
}

// Cot returns the cotangent of x (radians) to prec bits.
func Cot(x *big.Float, prec uint) *big.Float {
	wp := expGuard(x, prec) + 8
	s, c := bfSinCos(clone(wp, x), halfPiAt(wp), wp)
	return roundTo(prec, new(big.Float).SetPrec(wp).Quo(c, s))
}

// Sec returns the secant 1/cos(x) to prec bits.
func Sec(x *big.Float, prec uint) *big.Float {
	wp := expGuard(x, prec) + 8
	_, c := bfSinCos(clone(wp, x), halfPiAt(wp), wp)
	return roundTo(prec, new(big.Float).SetPrec(wp).Quo(oneF(wp), c))
}

// Csc returns the cosecant 1/sin(x) to prec bits.
func Csc(x *big.Float, prec uint) *big.Float {
	wp := expGuard(x, prec) + 8
	s, _ := bfSinCos(clone(wp, x), halfPiAt(wp), wp)
	return roundTo(prec, new(big.Float).SetPrec(wp).Quo(oneF(wp), s))
}

// Sinpi returns sin(pi*x) to prec bits, reducing x modulo 2 exactly so it stays
// accurate for large x.
func Sinpi(x *big.Float, prec uint) *big.Float {
	wp := working(prec) + 8
	n := nearestInt(clone(wp, x))
	t := new(big.Float).SetPrec(wp).Sub(clone(wp, x), new(big.Float).SetPrec(wp).SetInt(n))
	arg := new(big.Float).SetPrec(wp).Mul(bfPi(wp), t)
	sin, _ := bfSinCos(arg, halfPiAt(wp), wp)
	if n.Bit(0) == 1 {
		sin.Neg(sin)
	}
	return roundTo(prec, sin)
}

// Cospi returns cos(pi*x) to prec bits, reducing x modulo 2 exactly.
func Cospi(x *big.Float, prec uint) *big.Float {
	wp := working(prec) + 8
	n := nearestInt(clone(wp, x))
	t := new(big.Float).SetPrec(wp).Sub(clone(wp, x), new(big.Float).SetPrec(wp).SetInt(n))
	arg := new(big.Float).SetPrec(wp).Mul(bfPi(wp), t)
	_, cos := bfSinCos(arg, halfPiAt(wp), wp)
	if n.Bit(0) == 1 {
		cos.Neg(cos)
	}
	return roundTo(prec, cos)
}

// Tanpi returns tan(pi*x) to prec bits.
func Tanpi(x *big.Float, prec uint) *big.Float {
	wp := working(prec) + 8
	s := Sinpi(clone(wp, x), wp)
	c := Cospi(clone(wp, x), wp)
	return roundTo(prec, new(big.Float).SetPrec(wp).Quo(s, c))
}

// Atan returns the arctangent of x to prec bits, in (-pi/2, pi/2).
func Atan(x *big.Float, prec uint) *big.Float {
	wp := working(prec)
	return roundTo(prec, bfAtan(clone(wp, x), bfPi(wp), wp))
}

// Atan2 returns the angle of the point (x, y) in (-pi, pi], to prec bits.
func Atan2(y, x *big.Float, prec uint) *big.Float {
	wp := working(prec)
	pi := bfPi(wp)
	switch {
	case x.Sign() == 0 && y.Sign() == 0:
		return newF(prec)
	case x.Sign() == 0:
		hp := mulPow2(clone(wp, pi), -1)
		if y.Sign() < 0 {
			hp.Neg(hp)
		}
		return roundTo(prec, hp)
	}
	r := new(big.Float).SetPrec(wp).Quo(clone(wp, y), clone(wp, x))
	a := bfAtan(r, pi, wp)
	if x.Sign() < 0 {
		if y.Sign() >= 0 {
			a.Add(a, pi)
		} else {
			a.Sub(a, pi)
		}
	}
	return roundTo(prec, a)
}

// Acot returns the arccotangent atan(1/x) to prec bits (principal value in
// (0, pi)); Acot(0) == pi/2.
func Acot(x *big.Float, prec uint) *big.Float {
	wp := working(prec)
	pi := bfPi(wp)
	if x.Sign() == 0 {
		return roundTo(prec, mulPow2(pi, -1))
	}
	inv := new(big.Float).SetPrec(wp).Quo(oneF(wp), clone(wp, x))
	a := bfAtan(inv, pi, wp)
	if x.Sign() < 0 {
		a.Add(a, pi)
	}
	return roundTo(prec, a)
}

// Asin returns the arcsine of x to prec bits, in [-pi/2, pi/2]. It returns
// ErrDomain unless -1 <= x <= 1.
func Asin(x *big.Float, prec uint) (*big.Float, error) {
	wp := working(prec) + 8
	one := oneF(wp)
	ax := new(big.Float).SetPrec(wp).Abs(clone(wp, x))
	if ax.Cmp(one) > 0 {
		return nil, fmt.Errorf("%w: Asin of %s", ErrDomain, String(x))
	}
	if ax.Cmp(one) == 0 {
		hp := mulPow2(bfPi(wp), -1)
		if x.Sign() < 0 {
			hp.Neg(hp)
		}
		return roundTo(prec, hp), nil
	}
	// asin(x) = atan(x / sqrt(1-x^2)).
	d := new(big.Float).SetPrec(wp).Mul(clone(wp, x), clone(wp, x))
	d.Sub(one, d)
	d.Sqrt(d)
	r := new(big.Float).SetPrec(wp).Quo(clone(wp, x), d)
	return roundTo(prec, bfAtan(r, bfPi(wp), wp)), nil
}

// Acos returns the arccosine of x to prec bits, in [0, pi]. It returns
// ErrDomain unless -1 <= x <= 1.
func Acos(x *big.Float, prec uint) (*big.Float, error) {
	wp := working(prec) + 8
	as, err := Asin(x, wp)
	if err != nil {
		return nil, fmt.Errorf("%w: Acos of %s", ErrDomain, String(x))
	}
	hp := mulPow2(bfPi(wp), -1)
	hp.Sub(hp, as)
	return roundTo(prec, hp), nil
}

// Asec returns the arcsecant acos(1/x) to prec bits. It returns ErrDomain if
// |x| < 1.
func Asec(x *big.Float, prec uint) (*big.Float, error) {
	wp := working(prec) + 8
	if new(big.Float).Abs(x).Cmp(oneF(x.Prec())) < 0 {
		return nil, fmt.Errorf("%w: Asec of %s", ErrDomain, String(x))
	}
	inv := new(big.Float).SetPrec(wp).Quo(oneF(wp), clone(wp, x))
	return Acos(inv, prec)
}

// Acsc returns the arccosecant asin(1/x) to prec bits. It returns ErrDomain if
// |x| < 1.
func Acsc(x *big.Float, prec uint) (*big.Float, error) {
	wp := working(prec) + 8
	if new(big.Float).Abs(x).Cmp(oneF(x.Prec())) < 0 {
		return nil, fmt.Errorf("%w: Acsc of %s", ErrDomain, String(x))
	}
	inv := new(big.Float).SetPrec(wp).Quo(oneF(wp), clone(wp, x))
	return Asin(inv, prec)
}

// Deg2Rad converts x degrees to radians to prec bits.
func Deg2Rad(x *big.Float, prec uint) *big.Float {
	wp := working(prec)
	r := new(big.Float).SetPrec(wp).Mul(clone(wp, x), bfPi(wp))
	r.Quo(r, intF(wp, 180))
	return roundTo(prec, r)
}

// Rad2Deg converts x radians to degrees to prec bits.
func Rad2Deg(x *big.Float, prec uint) *big.Float {
	wp := working(prec)
	r := new(big.Float).SetPrec(wp).Mul(clone(wp, x), intF(wp, 180))
	r.Quo(r, bfPi(wp))
	return roundTo(prec, r)
}

// Versin returns the versed sine 1 - cos(x) to prec bits.
func Versin(x *big.Float, prec uint) *big.Float {
	wp := expGuard(x, prec) + 4
	_, c := bfSinCos(clone(wp, x), halfPiAt(wp), wp)
	r := new(big.Float).SetPrec(wp).Sub(oneF(wp), c)
	return roundTo(prec, r)
}

// Coversin returns the coversed sine 1 - sin(x) to prec bits.
func Coversin(x *big.Float, prec uint) *big.Float {
	wp := expGuard(x, prec) + 4
	s, _ := bfSinCos(clone(wp, x), halfPiAt(wp), wp)
	r := new(big.Float).SetPrec(wp).Sub(oneF(wp), s)
	return roundTo(prec, r)
}

// Haversin returns the haversine (1 - cos(x))/2 to prec bits.
func Haversin(x *big.Float, prec uint) *big.Float {
	wp := expGuard(x, prec) + 4
	_, c := bfSinCos(clone(wp, x), halfPiAt(wp), wp)
	r := new(big.Float).SetPrec(wp).Sub(oneF(wp), c)
	return roundTo(prec, mulPow2(r, -1))
}

// Exsec returns the exsecant sec(x) - 1 to prec bits.
func Exsec(x *big.Float, prec uint) *big.Float {
	wp := expGuard(x, prec) + 8
	_, c := bfSinCos(clone(wp, x), halfPiAt(wp), wp)
	r := new(big.Float).SetPrec(wp).Quo(oneF(wp), c)
	r.Sub(r, oneF(wp))
	return roundTo(prec, r)
}
