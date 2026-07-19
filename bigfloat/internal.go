package bigfloat

import "math/big"

// This file holds the unexported series kernels shared by the exported
// elementary functions. Every kernel works at the precision it is handed (the
// caller is responsible for adding a guard band) and assumes its argument lies
// in the reduced range documented on the function.

// twoAtanh returns 2*atanh(y) = log((1+y)/(1-y)) for |y| < 1, evaluated at prec
// bits. It powers both natural logarithms and the inverse hyperbolic tangent.
func twoAtanh(y *big.Float, prec uint) *big.Float {
	sum := clone(prec, y)
	term := clone(prec, y)
	y2 := new(big.Float).SetPrec(prec).Mul(y, y)
	tmp := newF(prec)
	for k := 1; k < int(prec)+10; k++ {
		term.Mul(term, y2)
		tmp.Quo(term, intF(prec, int64(2*k+1)))
		sum.Add(sum, tmp)
		if tiny(tmp, sum, prec) {
			break
		}
	}
	return sum.Mul(sum, intF(prec, 2))
}

// bfLn2 returns log(2) at prec bits via 2*atanh(1/3).
func bfLn2(prec uint) *big.Float {
	third := new(big.Float).SetPrec(prec).Quo(oneF(prec), intF(prec, 3))
	return twoAtanh(third, prec)
}

// bfLn returns log(x) for x > 0 at prec bits.
func bfLn(x *big.Float, prec uint) *big.Float {
	// x = mant * 2**e with mant in [0.5,1).
	mant := newF(prec)
	e := x.MantExp(mant)
	// Bring mant close to 1 with a few square roots so the atanh argument is
	// small: log(mant) = 2**s * log(mant**(1/2**s)).
	s := 0
	nine10 := new(big.Float).SetPrec(prec).SetFloat64(0.9)
	for mant.Cmp(nine10) < 0 {
		mant.Sqrt(mant)
		s++
	}
	// y = (mant-1)/(mant+1); log(mant**(1/2**s)) = twoAtanh(y).
	num := new(big.Float).SetPrec(prec).Sub(mant, oneF(prec))
	den := new(big.Float).SetPrec(prec).Add(mant, oneF(prec))
	y := num.Quo(num, den)
	lm := twoAtanh(y, prec)
	lm = mulPow2(lm, s) // undo the square-root reductions
	// log(x) = log(mant) + e*log(2).
	res := new(big.Float).SetPrec(prec).SetInt64(int64(e))
	res.Mul(res, bfLn2(prec))
	res.Add(res, lm)
	return res
}

// bfExpSmall returns exp(s) for |s| small (well under 1) at prec bits, by
// direct Taylor summation.
func bfExpSmall(s *big.Float, prec uint) *big.Float {
	sum := oneF(prec)
	term := oneF(prec)
	for k := 1; k < int(prec)+20; k++ {
		term = new(big.Float).SetPrec(prec).Mul(term, s)
		term.Quo(term, intF(prec, int64(k)))
		sum.Add(sum, term)
		if tiny(term, sum, prec) {
			break
		}
	}
	return sum
}

// bfExp returns exp(x) at prec bits for any finite x.
func bfExp(x *big.Float, prec uint) *big.Float {
	if x.Sign() == 0 {
		return oneF(prec)
	}
	ln2 := bfLn2(prec)
	// k = round(x/ln2); r = x - k*ln2, |r| <= ln2/2.
	kf := new(big.Float).SetPrec(prec).Quo(x, ln2)
	ki := nearestInt(kf)
	kln2 := new(big.Float).SetPrec(prec).SetInt(ki)
	kln2.Mul(kln2, ln2)
	r := new(big.Float).SetPrec(prec).Sub(x, kln2)
	// Scale r down by 2**t, sum the series, then square back t times.
	t := seriesIters(prec)
	s := new(big.Float).SetPrec(prec).SetMantExp(r, -t)
	e := bfExpSmall(s, prec)
	for i := 0; i < t; i++ {
		e.Mul(e, e)
	}
	// exp(x) = exp(r) * 2**k.
	k := int(ki.Int64())
	return new(big.Float).SetPrec(prec).SetMantExp(e, k)
}

// bfAtanCore returns atan(a) for 0 <= a <= 1 at prec bits.
func bfAtanCore(a *big.Float, prec uint) *big.Float {
	if a.Sign() == 0 {
		return newF(prec)
	}
	x := clone(prec, a)
	// Argument halving: atan(x) = 2*atan(x/(1+sqrt(1+x^2))).
	m := 0
	thr := pow2(prec, -6)
	one := oneF(prec)
	for x.Cmp(thr) > 0 {
		x2 := new(big.Float).SetPrec(prec).Mul(x, x)
		x2.Add(x2, one)
		x2.Sqrt(x2)
		x2.Add(x2, one)
		x.Quo(x, x2)
		m++
	}
	// Taylor series for the reduced argument.
	negx2 := new(big.Float).SetPrec(prec).Mul(x, x)
	negx2.Neg(negx2)
	term := clone(prec, x)
	sum := clone(prec, x)
	tmp := newF(prec)
	for k := 1; k < int(prec)+10; k++ {
		term.Mul(term, negx2)
		tmp.Quo(term, intF(prec, int64(2*k+1)))
		sum.Add(sum, tmp)
		if tiny(tmp, sum, prec) {
			break
		}
	}
	return mulPow2(sum, m)
}

// bfAtan returns atan(x) at prec bits for any finite x. It needs pi for the
// large-argument reflection and therefore takes it as a parameter to avoid
// recomputation.
func bfAtan(x *big.Float, pi *big.Float, prec uint) *big.Float {
	if x.Sign() == 0 {
		return newF(prec)
	}
	neg := x.Signbit()
	ax := new(big.Float).SetPrec(prec).Abs(x)
	var res *big.Float
	if ax.Cmp(oneF(prec)) > 0 {
		inv := new(big.Float).SetPrec(prec).Quo(oneF(prec), ax)
		res = new(big.Float).SetPrec(prec).Mul(pi, new(big.Float).SetFloat64(0.5))
		res.Sub(res, bfAtanCore(inv, prec))
	} else {
		res = bfAtanCore(ax, prec)
	}
	if neg {
		res.Neg(res)
	}
	return res
}

// bfSinCosSmall returns sin(s) and cos(s) for |s| small at prec bits.
func bfSinCosSmall(s *big.Float, prec uint) (sin, cos *big.Float) {
	s2 := new(big.Float).SetPrec(prec).Mul(s, s)
	negs2 := new(big.Float).SetPrec(prec).Neg(s2)
	// sin.
	sin = clone(prec, s)
	tsin := clone(prec, s)
	for n := 1; n < int(prec)+20; n++ {
		d := intF(prec, int64((2*n)*(2*n+1)))
		tsin.Mul(tsin, negs2)
		tsin.Quo(tsin, d)
		sin.Add(sin, tsin)
		if tiny(tsin, sin, prec) {
			break
		}
	}
	// cos.
	cos = oneF(prec)
	tcos := oneF(prec)
	for n := 1; n < int(prec)+20; n++ {
		d := intF(prec, int64((2*n-1)*(2*n)))
		tcos.Mul(tcos, negs2)
		tcos.Quo(tcos, d)
		cos.Add(cos, tcos)
		if tiny(tcos, cos, prec) {
			break
		}
	}
	return sin, cos
}

// bfSinCos returns sin(x) and cos(x) at prec bits for any finite x. halfPi must
// be pi/2 at (at least) prec bits.
func bfSinCos(x, halfPi *big.Float, prec uint) (sin, cos *big.Float) {
	// k = round(x/(pi/2)); r = x - k*(pi/2) in [-pi/4, pi/4].
	kf := new(big.Float).SetPrec(prec).Quo(x, halfPi)
	ki := nearestInt(kf)
	khp := new(big.Float).SetPrec(prec).SetInt(ki)
	khp.Mul(khp, halfPi)
	r := new(big.Float).SetPrec(prec).Sub(x, khp)
	// Scale down, series, double-angle back up.
	t := seriesIters(prec)
	s := new(big.Float).SetPrec(prec).SetMantExp(r, -t)
	sinr, cosr := bfSinCosSmall(s, prec)
	for i := 0; i < t; i++ {
		// sin(2u)=2 sin u cos u ; cos(2u)=cos^2 u - sin^2 u.
		twosc := new(big.Float).SetPrec(prec).Mul(sinr, cosr)
		twosc.Mul(twosc, intF(prec, 2))
		c2 := new(big.Float).SetPrec(prec).Mul(cosr, cosr)
		s2 := new(big.Float).SetPrec(prec).Mul(sinr, sinr)
		c2.Sub(c2, s2)
		sinr, cosr = twosc, c2
	}
	// Quadrant from k mod 4.
	q := int(new(big.Int).Mod(ki, big.NewInt(4)).Int64())
	switch q {
	case 0:
		sin, cos = sinr, cosr
	case 1:
		sin, cos = cosr, new(big.Float).SetPrec(prec).Neg(sinr)
	case 2:
		sin, cos = new(big.Float).SetPrec(prec).Neg(sinr), new(big.Float).SetPrec(prec).Neg(cosr)
	default: // 3
		sin, cos = new(big.Float).SetPrec(prec).Neg(cosr), sinr
	}
	return sin, cos
}

// expGuard returns the working precision for exp/trig of x, adding bits to
// cover cancellation in the range reduction of large arguments.
func expGuard(x *big.Float, prec uint) uint {
	wp := working(prec)
	if e := x.MantExp(nil); e > 0 {
		wp += uint(e)
	}
	return wp
}
