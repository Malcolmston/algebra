package bigfloat

import (
	"fmt"
	"math"
	"math/big"
)

// -----------------------------------------------------------------------------
// Integer combinatorics (exact).
// -----------------------------------------------------------------------------

// factorialBigInt computes n! as a *big.Int.
func factorialBigInt(n uint) *big.Int {
	r := big.NewInt(1)
	for i := int64(2); i <= int64(n); i++ {
		r.Mul(r, big.NewInt(i))
	}
	return r
}

// binomBigInt computes C(n,k) as a *big.Int for 0 <= k <= n.
func binomBigInt(n, k int) *big.Int {
	if k < 0 || k > n {
		return big.NewInt(0)
	}
	if k > n-k {
		k = n - k
	}
	r := big.NewInt(1)
	for i := 0; i < k; i++ {
		r.Mul(r, big.NewInt(int64(n-i)))
		r.Div(r, big.NewInt(int64(i+1)))
	}
	return r
}

// FactorialBig returns n! as an exact *big.Int.
func FactorialBig(n uint) *big.Int { return factorialBigInt(n) }

// Factorial returns n! rounded to prec bits.
func Factorial(n uint, prec uint) *big.Float {
	return new(big.Float).SetPrec(prec).SetInt(factorialBigInt(n))
}

// DoubleFactorialBig returns the double factorial n!! as an exact *big.Int.
func DoubleFactorialBig(n uint) *big.Int {
	r := big.NewInt(1)
	for i := int64(n); i > 1; i -= 2 {
		r.Mul(r, big.NewInt(i))
	}
	return r
}

// DoubleFactorial returns the double factorial n!! rounded to prec bits.
func DoubleFactorial(n uint, prec uint) *big.Float {
	return new(big.Float).SetPrec(prec).SetInt(DoubleFactorialBig(n))
}

// BinomialBig returns the binomial coefficient C(n,k) as an exact *big.Int.
func BinomialBig(n, k uint) *big.Int { return binomBigInt(int(n), int(k)) }

// Binomial returns the binomial coefficient C(n,k) rounded to prec bits.
func Binomial(n, k uint, prec uint) *big.Float {
	return new(big.Float).SetPrec(prec).SetInt(binomBigInt(int(n), int(k)))
}

// RisingFactorial returns the Pochhammer symbol x*(x+1)*...*(x+n-1) to prec
// bits.
func RisingFactorial(x *big.Float, n uint, prec uint) *big.Float {
	wp := working(prec) + 8
	r := oneF(wp)
	for i := uint(0); i < n; i++ {
		t := new(big.Float).SetPrec(wp).Add(clone(wp, x), intF(wp, int64(i)))
		r.Mul(r, t)
	}
	return roundTo(prec, r)
}

// Pochhammer is an alias for RisingFactorial.
func Pochhammer(x *big.Float, n uint, prec uint) *big.Float { return RisingFactorial(x, n, prec) }

// FallingFactorial returns x*(x-1)*...*(x-n+1) to prec bits.
func FallingFactorial(x *big.Float, n uint, prec uint) *big.Float {
	wp := working(prec) + 8
	r := oneF(wp)
	for i := uint(0); i < n; i++ {
		t := new(big.Float).SetPrec(wp).Sub(clone(wp, x), intF(wp, int64(i)))
		r.Mul(r, t)
	}
	return roundTo(prec, r)
}

// -----------------------------------------------------------------------------
// Bernoulli numbers.
// -----------------------------------------------------------------------------

// computeBernoulli returns B_0..B_upTo as exact rationals using the recurrence
// sum_{j=0}^{m} C(m+1,j) B_j = 0 for m >= 1.
func computeBernoulli(upTo int) []*big.Rat {
	b := make([]*big.Rat, upTo+1)
	b[0] = big.NewRat(1, 1)
	for m := 1; m <= upTo; m++ {
		s := new(big.Rat)
		for j := 0; j < m; j++ {
			c := new(big.Rat).SetInt(binomBigInt(m+1, j))
			s.Add(s, new(big.Rat).Mul(c, b[j]))
		}
		s.Neg(s)
		s.Quo(s, big.NewRat(int64(m+1), 1))
		b[m] = s
	}
	return b
}

// Bernoulli returns the n-th Bernoulli number B_n as an exact *big.Rat, using
// the convention B_1 = -1/2.
func Bernoulli(n uint) *big.Rat {
	b := computeBernoulli(int(n))
	return b[n]
}

// BernoulliFloat returns the n-th Bernoulli number rounded to prec bits.
func BernoulliFloat(n uint, prec uint) *big.Float {
	return new(big.Float).SetPrec(prec).SetRat(Bernoulli(n))
}

// -----------------------------------------------------------------------------
// Gamma family.
// -----------------------------------------------------------------------------

// lgammaPos returns log(Gamma(x)) for x > 0 at prec bits, via the Stirling
// asymptotic series after shifting the argument up to a large value.
func lgammaPos(x *big.Float, prec uint) *big.Float {
	xf, _ := x.Float64()
	zmin := float64(prec)
	N := int(math.Ceil(zmin - xf))
	if N < 0 {
		N = 0
	}
	// P = x(x+1)...(x+N-1); z = x+N.
	lnP := newF(prec)
	if N > 0 {
		P := oneF(prec)
		for j := 0; j < N; j++ {
			P.Mul(P, new(big.Float).SetPrec(prec).Add(clone(prec, x), intF(prec, int64(j))))
		}
		lnP = bfLn(P, prec)
	}
	z := new(big.Float).SetPrec(prec).Add(clone(prec, x), intF(prec, int64(N)))
	lnz := bfLn(z, prec)
	// (z-1/2)*ln z - z + 1/2*ln(2pi).
	lg := new(big.Float).SetPrec(prec).Sub(z, new(big.Float).SetPrec(prec).SetFloat64(0.5))
	lg.Mul(lg, lnz)
	lg.Sub(lg, z)
	ln2pi := bfLn(mulPow2(bfPi(prec), 1), prec)
	lg.Add(lg, mulPow2(ln2pi, -1))
	// Correction series: sum B_{2k}/((2k)(2k-1) z^{2k-1}).
	zinv := new(big.Float).SetPrec(prec).Quo(oneF(prec), z)
	zinv2 := new(big.Float).SetPrec(prec).Mul(zinv, zinv)
	zp := clone(prec, zinv) // z^{-(2k-1)} for k=1
	bern := computeBernoulli(2 * (int(prec) + 4))
	prevMag := math.MaxInt
	for k := 1; k <= int(prec); k++ {
		if 2*k >= len(bern) {
			break
		}
		coef := new(big.Float).SetPrec(prec).SetRat(bern[2*k])
		coef.Quo(coef, intF(prec, int64((2*k)*(2*k-1))))
		term := new(big.Float).SetPrec(prec).Mul(coef, zp)
		mag := term.MantExp(nil)
		if term.Sign() != 0 && mag > prevMag {
			break // asymptotic series diverging
		}
		prevMag = mag
		lg.Add(lg, term)
		if tiny(term, lg, prec) {
			break
		}
		zp.Mul(zp, zinv2)
	}
	lg.Sub(lg, lnP)
	return lg
}

// Gamma returns the gamma function Gamma(x) to prec bits, for any real x that
// is not a non-positive integer (a pole), for which it returns ErrPole.
func Gamma(x *big.Float, prec uint) (*big.Float, error) {
	wp := working(prec) + 16
	if x.IsInt() {
		if x.Sign() <= 0 {
			return nil, fmt.Errorf("%w: Gamma at non-positive integer %s", ErrPole, String(x))
		}
		// Gamma(n) = (n-1)!.
		xi, _ := x.Int(nil)
		n := xi.Int64()
		return new(big.Float).SetPrec(prec).SetInt(factorialBigInt(uint(n - 1))), nil
	}
	if x.Sign() > 0 {
		return roundTo(prec, bfExp(lgammaPos(clone(wp, x), wp), wp)), nil
	}
	// Reflection: Gamma(x) = pi / (sin(pi x) * Gamma(1-x)).
	omx := new(big.Float).SetPrec(wp).Sub(oneF(wp), clone(wp, x)) // 1-x > 1
	g1 := bfExp(lgammaPos(omx, wp), wp)
	sp := Sinpi(clone(wp, x), wp)
	den := new(big.Float).SetPrec(wp).Mul(sp, g1)
	res := new(big.Float).SetPrec(wp).Quo(bfPi(wp), den)
	return roundTo(prec, res), nil
}

// Lgamma returns the natural logarithm of the absolute value of Gamma(x) and
// the sign of Gamma(x) (+1 or -1), to prec bits. It returns ErrPole at the
// non-positive integers.
func Lgamma(x *big.Float, prec uint) (val *big.Float, sign int, err error) {
	wp := working(prec) + 16
	if x.IsInt() && x.Sign() <= 0 {
		return nil, 0, fmt.Errorf("%w: Lgamma at non-positive integer %s", ErrPole, String(x))
	}
	if x.Sign() > 0 {
		return roundTo(prec, lgammaPos(clone(wp, x), wp)), 1, nil
	}
	// ln|Gamma(x)| = ln(pi) - ln|sin(pi x)| - lgammaPos(1-x).
	omx := new(big.Float).SetPrec(wp).Sub(oneF(wp), clone(wp, x))
	lg1 := lgammaPos(omx, wp)
	sp := Sinpi(clone(wp, x), wp)
	asp := new(big.Float).SetPrec(wp).Abs(sp)
	v := new(big.Float).SetPrec(wp).Sub(bfLn(bfPi(wp), wp), bfLn(asp, wp))
	v.Sub(v, lg1)
	sign = 1
	if sp.Sign() < 0 {
		sign = -1
	}
	return roundTo(prec, v), sign, nil
}

// Digamma returns the digamma function psi(x) = Gamma'(x)/Gamma(x) to prec bits,
// for any real x that is not a non-positive integer, for which it returns
// ErrPole.
func Digamma(x *big.Float, prec uint) (*big.Float, error) {
	wp := working(prec) + 16
	if x.IsInt() && x.Sign() <= 0 {
		return nil, fmt.Errorf("%w: Digamma at non-positive integer %s", ErrPole, String(x))
	}
	if x.Sign() > 0 {
		return roundTo(prec, digammaPos(clone(wp, x), wp)), nil
	}
	// Reflection: psi(1-x) - psi(x) = pi*cot(pi x).
	omx := new(big.Float).SetPrec(wp).Sub(oneF(wp), clone(wp, x))
	p1 := digammaPos(omx, wp)
	s := Sinpi(clone(wp, x), wp)
	c := Cospi(clone(wp, x), wp)
	cot := new(big.Float).SetPrec(wp).Quo(c, s)
	cot.Mul(cot, bfPi(wp))
	res := new(big.Float).SetPrec(wp).Sub(p1, cot)
	return roundTo(prec, res), nil
}

// digammaPos returns psi(x) for x > 0 at prec bits, by upward recurrence into
// the Stirling regime.
func digammaPos(x *big.Float, prec uint) *big.Float {
	xf, _ := x.Float64()
	zmin := float64(prec)
	N := int(math.Ceil(zmin - xf))
	if N < 0 {
		N = 0
	}
	// psi(x) = psi(x+N) - sum_{j=0}^{N-1} 1/(x+j).
	corr := newF(prec)
	for j := 0; j < N; j++ {
		d := new(big.Float).SetPrec(prec).Add(clone(prec, x), intF(prec, int64(j)))
		corr.Add(corr, new(big.Float).SetPrec(prec).Quo(oneF(prec), d))
	}
	z := new(big.Float).SetPrec(prec).Add(clone(prec, x), intF(prec, int64(N)))
	// psi(z) = ln z - 1/(2z) - sum B_{2k}/(2k z^{2k}).
	res := bfLn(z, prec)
	half := new(big.Float).SetPrec(prec).Quo(oneF(prec), z)
	res.Sub(res, mulPow2(half, -1))
	zinv := new(big.Float).SetPrec(prec).Quo(oneF(prec), z)
	zinv2 := new(big.Float).SetPrec(prec).Mul(zinv, zinv)
	zp := clone(prec, zinv2) // z^{-2k} for k=1
	bern := computeBernoulli(2 * (int(prec) + 4))
	prevMag := math.MaxInt
	for k := 1; k <= int(prec); k++ {
		if 2*k >= len(bern) {
			break
		}
		coef := new(big.Float).SetPrec(prec).SetRat(bern[2*k])
		coef.Quo(coef, intF(prec, int64(2*k)))
		term := new(big.Float).SetPrec(prec).Mul(coef, zp)
		mag := term.MantExp(nil)
		if term.Sign() != 0 && mag > prevMag {
			break
		}
		prevMag = mag
		res.Sub(res, term)
		if tiny(term, res, prec) {
			break
		}
		zp.Mul(zp, zinv2)
	}
	res.Sub(res, corr)
	return res
}

// Beta returns the beta function B(a,b) = Gamma(a)Gamma(b)/Gamma(a+b) to prec
// bits, for a > 0 and b > 0. It returns ErrDomain otherwise.
func Beta(a, b *big.Float, prec uint) (*big.Float, error) {
	if a.Sign() <= 0 || b.Sign() <= 0 {
		return nil, fmt.Errorf("%w: Beta requires positive arguments", ErrDomain)
	}
	wp := working(prec) + 16
	la := lgammaPos(clone(wp, a), wp)
	lb := lgammaPos(clone(wp, b), wp)
	ab := new(big.Float).SetPrec(wp).Add(clone(wp, a), clone(wp, b))
	lab := lgammaPos(ab, wp)
	s := new(big.Float).SetPrec(wp).Add(la, lb)
	s.Sub(s, lab)
	return roundTo(prec, bfExp(s, wp)), nil
}

// LogBeta returns log(B(a,b)) to prec bits, for a > 0 and b > 0. It returns
// ErrDomain otherwise.
func LogBeta(a, b *big.Float, prec uint) (*big.Float, error) {
	if a.Sign() <= 0 || b.Sign() <= 0 {
		return nil, fmt.Errorf("%w: LogBeta requires positive arguments", ErrDomain)
	}
	wp := working(prec) + 16
	la := lgammaPos(clone(wp, a), wp)
	lb := lgammaPos(clone(wp, b), wp)
	ab := new(big.Float).SetPrec(wp).Add(clone(wp, a), clone(wp, b))
	lab := lgammaPos(ab, wp)
	s := new(big.Float).SetPrec(wp).Add(la, lb)
	s.Sub(s, lab)
	return roundTo(prec, s), nil
}

// BinomialReal returns the generalised binomial coefficient
// Gamma(x+1)/(Gamma(y+1)*Gamma(x-y+1)) to prec bits, for real x and y where all
// three gamma values are finite. It returns an error at a pole.
func BinomialReal(x, y *big.Float, prec uint) (*big.Float, error) {
	wp := working(prec) + 16
	gx, err := Gamma(new(big.Float).SetPrec(wp).Add(clone(wp, x), oneF(wp)), wp)
	if err != nil {
		return nil, err
	}
	gy, err := Gamma(new(big.Float).SetPrec(wp).Add(clone(wp, y), oneF(wp)), wp)
	if err != nil {
		return nil, err
	}
	xy := new(big.Float).SetPrec(wp).Sub(clone(wp, x), clone(wp, y))
	xy.Add(xy, oneF(wp))
	gxy, err := Gamma(xy, wp)
	if err != nil {
		return nil, err
	}
	r := new(big.Float).SetPrec(wp).Quo(gx, new(big.Float).SetPrec(wp).Mul(gy, gxy))
	return roundTo(prec, r), nil
}
