package bigfloat

import (
	"fmt"
	"math"
	"math/big"
)

// AGM returns the arithmetic-geometric mean of a and b to prec bits. It
// requires a >= 0 and b >= 0 and returns ErrNegative otherwise.
func AGM(a, b *big.Float, prec uint) (*big.Float, error) {
	if a.Sign() < 0 || b.Sign() < 0 {
		return nil, fmt.Errorf("%w: AGM requires non-negative arguments", ErrNegative)
	}
	wp := working(prec) + 16
	an := clone(wp, a)
	bn := clone(wp, b)
	for i := 0; i < int(wp); i++ {
		nextA := new(big.Float).SetPrec(wp).Add(an, bn)
		nextA = mulPow2(nextA, -1)
		nextB := new(big.Float).SetPrec(wp).Mul(an, bn)
		nextB.Sqrt(nextB)
		diff := new(big.Float).SetPrec(wp).Sub(nextA, nextB)
		an, bn = nextA, nextB
		if tiny(diff, an, prec) {
			break
		}
	}
	return roundTo(prec, an), nil
}

// EllipticK returns the complete elliptic integral of the first kind K(k) to
// prec bits, using K(k) = pi / (2 * AGM(1, sqrt(1-k^2))). It returns ErrDomain
// unless -1 < k < 1.
func EllipticK(k *big.Float, prec uint) (*big.Float, error) {
	if new(big.Float).Abs(k).Cmp(oneF(k.Prec())) >= 0 {
		return nil, fmt.Errorf("%w: EllipticK requires |k| < 1", ErrDomain)
	}
	wp := working(prec) + 16
	kp := new(big.Float).SetPrec(wp).Mul(clone(wp, k), clone(wp, k))
	kp.Sub(oneF(wp), kp)
	kp.Sqrt(kp) // sqrt(1-k^2)
	m, err := AGM(oneF(wp), kp, wp)
	if err != nil {
		return nil, err
	}
	res := new(big.Float).SetPrec(wp).Quo(bfPi(wp), mulPow2(m, 1))
	return roundTo(prec, res), nil
}

// EllipticE returns the complete elliptic integral of the second kind E(k) to
// prec bits, via the descending AGM sequence
//
//	E(k) = K(k) * (1 - sum_{n>=0} 2^{n-1} c_n^2),
//
// with a_0 = 1, b_0 = sqrt(1-k^2), c_0 = k. It returns ErrDomain unless
// -1 < k < 1.
func EllipticE(k *big.Float, prec uint) (*big.Float, error) {
	if new(big.Float).Abs(k).Cmp(oneF(k.Prec())) >= 0 {
		return nil, fmt.Errorf("%w: EllipticE requires |k| < 1", ErrDomain)
	}
	wp := working(prec) + 16
	a := oneF(wp)
	b := new(big.Float).SetPrec(wp).Mul(clone(wp, k), clone(wp, k))
	b.Sub(oneF(wp), b)
	b.Sqrt(b)
	c := new(big.Float).SetPrec(wp).Abs(clone(wp, k)) // c_0 = |k|
	// sum starts with 2^{-1} c_0^2.
	sum := new(big.Float).SetPrec(wp).Mul(c, c)
	sum = mulPow2(sum, -1)
	for n := 1; n < int(wp); n++ {
		nextA := mulPow2(new(big.Float).SetPrec(wp).Add(a, b), -1)
		nextB := new(big.Float).SetPrec(wp).Mul(a, b)
		nextB.Sqrt(nextB)
		cn := mulPow2(new(big.Float).SetPrec(wp).Sub(a, b), -1)
		a, b, c = nextA, nextB, cn
		csq := new(big.Float).SetPrec(wp).Mul(c, c)
		term := mulPow2(csq, n-1) // 2^{n-1} c_n^2
		sum.Add(sum, term)
		if tiny(term, sum, prec) {
			break
		}
	}
	kk := new(big.Float).SetPrec(wp).Quo(bfPi(wp), mulPow2(a, 1)) // K(k)
	factor := new(big.Float).SetPrec(wp).Sub(oneF(wp), sum)
	return roundTo(prec, new(big.Float).SetPrec(wp).Mul(kk, factor)), nil
}

// LambertW returns the principal branch W0(x) of the Lambert W function to prec
// bits, i.e. the solution w of w*e^w = x with w >= -1. It is computed by
// Halley's (Newton-type) iteration. It returns ErrDomain if x < -1/e.
func LambertW(x *big.Float, prec uint) (*big.Float, error) {
	wp := working(prec) + 16
	// Domain check against -1/e.
	minv := new(big.Float).SetPrec(wp).Quo(oneF(wp), bfExp(oneF(wp), wp))
	minv.Neg(minv)
	if x.Cmp(minv) < 0 {
		return nil, fmt.Errorf("%w: LambertW requires x >= -1/e", ErrDomain)
	}
	if x.Sign() == 0 {
		return newF(prec), nil
	}
	// Initial guess from float64.
	xf, _ := x.Float64()
	var w0 float64
	switch {
	case xf < -0.3:
		// near the branch point; expand in p = sqrt(2(e x + 1)).
		p := math.Sqrt(2 * (math.E*xf + 1))
		w0 = -1 + p - p*p/3
	case xf < 10:
		w0 = math.Log1p(xf) // decent for moderate positive x
		if w0 == 0 {
			w0 = xf
		}
	default:
		l := math.Log(xf)
		w0 = l - math.Log(l)
	}
	w := new(big.Float).SetPrec(wp).SetFloat64(w0)
	xw := clone(wp, x)
	for i := 0; i < int(wp); i++ {
		ew := bfExp(w, wp)
		we := new(big.Float).SetPrec(wp).Mul(w, ew) // w e^w
		f := new(big.Float).SetPrec(wp).Sub(we, xw) // w e^w - x
		// Halley: w -= f / (e^w (w+1) - (w+2) f / (2 w + 2)).
		wp1 := new(big.Float).SetPrec(wp).Add(w, oneF(wp))
		denom := new(big.Float).SetPrec(wp).Mul(ew, wp1)
		wp2 := new(big.Float).SetPrec(wp).Add(w, intF(wp, 2))
		num2 := new(big.Float).SetPrec(wp).Mul(wp2, f)
		den2 := mulPow2(wp1, 1) // 2(w+1)
		corr := new(big.Float).SetPrec(wp).Quo(num2, den2)
		denom.Sub(denom, corr)
		dw := new(big.Float).SetPrec(wp).Quo(f, denom)
		w.Sub(w, dw)
		if tiny(dw, w, prec) {
			break
		}
	}
	return roundTo(prec, w), nil
}
