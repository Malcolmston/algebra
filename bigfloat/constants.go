package bigfloat

import "math/big"

// bfPi returns pi at prec bits via the Machin formula
// pi = 16*atan(1/5) - 4*atan(1/239).
func bfPi(prec uint) *big.Float {
	fifth := new(big.Float).SetPrec(prec).Quo(oneF(prec), intF(prec, 5))
	c239 := new(big.Float).SetPrec(prec).Quo(oneF(prec), intF(prec, 239))
	a := bfAtanCore(fifth, prec)
	b := bfAtanCore(c239, prec)
	a.Mul(a, intF(prec, 16))
	b.Mul(b, intF(prec, 4))
	return a.Sub(a, b)
}

// Pi returns the ratio of a circle's circumference to its diameter, pi, to prec
// bits.
func Pi(prec uint) *big.Float {
	wp := working(prec)
	return roundTo(prec, bfPi(wp))
}

// TwoPi returns 2*pi to prec bits.
func TwoPi(prec uint) *big.Float {
	wp := working(prec)
	return roundTo(prec, mulPow2(bfPi(wp), 1))
}

// HalfPi returns pi/2 to prec bits.
func HalfPi(prec uint) *big.Float {
	wp := working(prec)
	return roundTo(prec, mulPow2(bfPi(wp), -1))
}

// QuarterPi returns pi/4 to prec bits.
func QuarterPi(prec uint) *big.Float {
	wp := working(prec)
	return roundTo(prec, mulPow2(bfPi(wp), -2))
}

// InvPi returns 1/pi to prec bits.
func InvPi(prec uint) *big.Float {
	wp := working(prec)
	return roundTo(prec, new(big.Float).SetPrec(wp).Quo(oneF(wp), bfPi(wp)))
}

// TwoInvPi returns 2/pi to prec bits.
func TwoInvPi(prec uint) *big.Float {
	wp := working(prec)
	r := new(big.Float).SetPrec(wp).Quo(intF(wp, 2), bfPi(wp))
	return roundTo(prec, r)
}

// Degree returns the size of one degree in radians, pi/180, to prec bits.
func Degree(prec uint) *big.Float {
	wp := working(prec)
	r := new(big.Float).SetPrec(wp).Quo(bfPi(wp), intF(wp, 180))
	return roundTo(prec, r)
}

// E returns Euler's number e = exp(1) to prec bits.
func E(prec uint) *big.Float {
	wp := working(prec)
	return roundTo(prec, bfExp(oneF(wp), wp))
}

// InvE returns 1/e to prec bits.
func InvE(prec uint) *big.Float {
	wp := working(prec)
	one := oneF(wp)
	neg := new(big.Float).SetPrec(wp).Neg(one)
	return roundTo(prec, bfExp(neg, wp))
}

// Ln2 returns the natural logarithm of 2 to prec bits.
func Ln2(prec uint) *big.Float {
	wp := working(prec)
	return roundTo(prec, bfLn2(wp))
}

// Ln10 returns the natural logarithm of 10 to prec bits.
func Ln10(prec uint) *big.Float {
	wp := working(prec)
	return roundTo(prec, bfLn(intF(wp, 10), wp))
}

// LnPi returns the natural logarithm of pi to prec bits.
func LnPi(prec uint) *big.Float {
	wp := working(prec)
	return roundTo(prec, bfLn(bfPi(wp), wp))
}

// Ln2Pi returns the natural logarithm of 2*pi to prec bits.
func Ln2Pi(prec uint) *big.Float {
	wp := working(prec)
	return roundTo(prec, bfLn(mulPow2(bfPi(wp), 1), wp))
}

// Log2E returns log base 2 of e, i.e. 1/ln(2), to prec bits.
func Log2E(prec uint) *big.Float {
	wp := working(prec)
	r := new(big.Float).SetPrec(wp).Quo(oneF(wp), bfLn2(wp))
	return roundTo(prec, r)
}

// Log10E returns log base 10 of e, i.e. 1/ln(10), to prec bits.
func Log10E(prec uint) *big.Float {
	wp := working(prec)
	r := new(big.Float).SetPrec(wp).Quo(oneF(wp), bfLn(intF(wp, 10), wp))
	return roundTo(prec, r)
}

// bfSqrtInt returns sqrt(n) at prec bits.
func bfSqrtInt(n int64, prec uint) *big.Float {
	return new(big.Float).SetPrec(prec).Sqrt(intF(prec, n))
}

// Sqrt2 returns the square root of 2 to prec bits.
func Sqrt2(prec uint) *big.Float { wp := working(prec); return roundTo(prec, bfSqrtInt(2, wp)) }

// Sqrt3 returns the square root of 3 to prec bits.
func Sqrt3(prec uint) *big.Float { wp := working(prec); return roundTo(prec, bfSqrtInt(3, wp)) }

// Sqrt5 returns the square root of 5 to prec bits.
func Sqrt5(prec uint) *big.Float { wp := working(prec); return roundTo(prec, bfSqrtInt(5, wp)) }

// InvSqrt2 returns 1/sqrt(2) to prec bits.
func InvSqrt2(prec uint) *big.Float {
	wp := working(prec)
	r := new(big.Float).SetPrec(wp).Quo(oneF(wp), bfSqrtInt(2, wp))
	return roundTo(prec, r)
}

// SqrtPi returns sqrt(pi) to prec bits.
func SqrtPi(prec uint) *big.Float {
	wp := working(prec)
	r := new(big.Float).SetPrec(wp).Sqrt(bfPi(wp))
	return roundTo(prec, r)
}

// SqrtTwoPi returns sqrt(2*pi) to prec bits.
func SqrtTwoPi(prec uint) *big.Float {
	wp := working(prec)
	r := new(big.Float).SetPrec(wp).Sqrt(mulPow2(bfPi(wp), 1))
	return roundTo(prec, r)
}

// InvSqrtPi returns 1/sqrt(pi) to prec bits.
func InvSqrtPi(prec uint) *big.Float {
	wp := working(prec)
	s := new(big.Float).SetPrec(wp).Sqrt(bfPi(wp))
	r := new(big.Float).SetPrec(wp).Quo(oneF(wp), s)
	return roundTo(prec, r)
}

// GoldenRatio returns the golden ratio (1+sqrt(5))/2 to prec bits.
func GoldenRatio(prec uint) *big.Float {
	wp := working(prec)
	r := new(big.Float).SetPrec(wp).Add(oneF(wp), bfSqrtInt(5, wp))
	return roundTo(prec, mulPow2(r, -1))
}

// InvGoldenRatio returns 1/phi = (sqrt(5)-1)/2 to prec bits.
func InvGoldenRatio(prec uint) *big.Float {
	wp := working(prec)
	r := new(big.Float).SetPrec(wp).Sub(bfSqrtInt(5, wp), oneF(wp))
	return roundTo(prec, mulPow2(r, -1))
}

// SilverRatio returns the silver ratio 1+sqrt(2) to prec bits.
func SilverRatio(prec uint) *big.Float {
	wp := working(prec)
	r := new(big.Float).SetPrec(wp).Add(oneF(wp), bfSqrtInt(2, wp))
	return roundTo(prec, r)
}

// EulerGamma returns the Euler-Mascheroni constant gamma to prec bits, computed
// with the Brent-McMillan series
//
//	gamma = S(n)/I(n) - ln(n),
//
// where I(n) = sum (n^k/k!)^2 and S(n) = sum (n^k/k!)^2 * H_k.
func EulerGamma(prec uint) *big.Float {
	wp := working(prec) + 16
	// Choose n so that the truncation error ~ pi*exp(-4n) is below 2^-wp.
	n := int(float64(wp)/5.0) + 6
	nf := intF(wp, int64(n))
	n2 := new(big.Float).SetPrec(wp).Mul(nf, nf)

	term := oneF(wp) // (n^k/k!)^2, k=0
	iSum := oneF(wp) // I(n)
	sSum := newF(wp) // S(n), H_0 = 0
	harm := newF(wp) // H_k
	for k := 1; k < 8*n+50; k++ {
		kf := intF(wp, int64(k))
		// term *= n^2 / k^2.
		term.Mul(term, n2)
		k2 := new(big.Float).SetPrec(wp).Mul(kf, kf)
		term.Quo(term, k2)
		// H_k = H_{k-1} + 1/k.
		harm.Add(harm, new(big.Float).SetPrec(wp).Quo(oneF(wp), kf))
		iSum.Add(iSum, term)
		sSum.Add(sSum, new(big.Float).SetPrec(wp).Mul(term, harm))
		if k > 2*n && tiny(term, iSum, wp) {
			break
		}
	}
	res := new(big.Float).SetPrec(wp).Quo(sSum, iSum)
	res.Sub(res, bfLn(nf, wp))
	return roundTo(prec, res)
}

// Catalan returns Catalan's constant G to prec bits using the rapidly
// convergent identity
//
//	G = (pi/8)*ln(2+sqrt(3)) + (3/8) * sum_{k>=0} 1/((2k+1)^2 * C(2k,k)).
func Catalan(prec uint) *big.Float {
	wp := working(prec) + 16
	// First piece: (pi/8) * ln(2+sqrt(3)).
	pi := bfPi(wp)
	root := new(big.Float).SetPrec(wp).Add(intF(wp, 2), bfSqrtInt(3, wp))
	first := new(big.Float).SetPrec(wp).Mul(pi, bfLn(root, wp))
	first = mulPow2(first, -3) // divide by 8
	// Second piece: (3/8) * sum t_k, t_0 = 1, t_k/t_{k-1} = k(2k-1)/(2(2k+1)^2).
	t := oneF(wp)
	sum := oneF(wp)
	for k := 1; k < 4*int(wp)+20; k++ {
		num := intF(wp, int64(k*(2*k-1)))
		den := intF(wp, int64(2*(2*k+1)*(2*k+1)))
		t = new(big.Float).SetPrec(wp).Mul(t, num)
		t.Quo(t, den)
		sum.Add(sum, t)
		if tiny(t, sum, wp) {
			break
		}
	}
	second := new(big.Float).SetPrec(wp).Mul(sum, intF(wp, 3))
	second = mulPow2(second, -3)
	res := new(big.Float).SetPrec(wp).Add(first, second)
	return roundTo(prec, res)
}

// Apery returns Apery's constant zeta(3) to prec bits using
//
//	zeta(3) = (5/2) * sum_{n>=1} (-1)^{n-1} / (n^3 * C(2n,n)).
func Apery(prec uint) *big.Float {
	wp := working(prec) + 16
	sum := newF(wp)
	// binom(2n,n) via incremental ratio: c_n = c_{n-1} * (2(2n-1)/n).
	c := big.NewInt(1) // C(0,0)
	for n := 1; n < 4*int(wp)+20; n++ {
		c.Mul(c, big.NewInt(int64(2*(2*n-1))))
		c.Div(c, big.NewInt(int64(n)))
		// term = (-1)^{n-1} / (n^3 * c).
		den := new(big.Float).SetPrec(wp).SetInt(c)
		den.Mul(den, intF(wp, int64(n)*int64(n)*int64(n)))
		term := new(big.Float).SetPrec(wp).Quo(oneF(wp), den)
		if n%2 == 0 {
			term.Neg(term)
		}
		sum.Add(sum, term)
		if tiny(term, sum, wp) {
			break
		}
	}
	res := new(big.Float).SetPrec(wp).Mul(sum, new(big.Float).SetPrec(wp).SetFloat64(2.5))
	return roundTo(prec, res)
}

// Zeta3 is an alias for Apery; it returns zeta(3) to prec bits.
func Zeta3(prec uint) *big.Float { return Apery(prec) }

// PlasticNumber returns the plastic number, the real root of x^3 = x + 1, to
// prec bits. It is found by Newton's method on f(x) = x^3 - x - 1.
func PlasticNumber(prec uint) *big.Float {
	wp := working(prec) + 8
	x := new(big.Float).SetPrec(wp).SetFloat64(1.3247179572447458)
	for i := 0; i < int(wp); i++ {
		x2 := new(big.Float).SetPrec(wp).Mul(x, x)
		f := new(big.Float).SetPrec(wp).Mul(x2, x)
		f.Sub(f, x)
		f.Sub(f, oneF(wp))
		fp := new(big.Float).SetPrec(wp).Mul(x2, intF(wp, 3))
		fp.Sub(fp, oneF(wp))
		dx := new(big.Float).SetPrec(wp).Quo(f, fp)
		x.Sub(x, dx)
		if tiny(dx, x, prec) {
			break
		}
	}
	return roundTo(prec, x)
}
