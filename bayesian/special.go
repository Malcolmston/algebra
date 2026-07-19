package bayesian

import (
	"math"
)

// LogGamma returns the natural logarithm of the absolute value of the Gamma
// function, ln|Γ(x)|. It is a thin, sign-discarding wrapper around the standard
// library's math.Lgamma and underlies most of the log-domain computations in
// this package.
func LogGamma(x float64) float64 {
	v, _ := math.Lgamma(x)
	return v
}

// GammaFunc returns the Gamma function Γ(x). It is a convenience wrapper around
// math.Gamma provided so callers of this package need not import math directly.
func GammaFunc(x float64) float64 {
	return math.Gamma(x)
}

// LogBeta returns the natural logarithm of the Beta function,
// ln B(a,b) = ln Γ(a) + ln Γ(b) − ln Γ(a+b), for a,b > 0.
func LogBeta(a, b float64) float64 {
	return LogGamma(a) + LogGamma(b) - LogGamma(a+b)
}

// BetaFunc returns the Beta function B(a,b) = Γ(a)Γ(b)/Γ(a+b) for a,b > 0.
func BetaFunc(a, b float64) float64 {
	return math.Exp(LogBeta(a, b))
}

// LogFactorial returns ln(n!) computed in the log-gamma domain as ln Γ(n+1).
// It returns NaN for negative n.
func LogFactorial(n int) float64 {
	if n < 0 {
		return math.NaN()
	}
	return LogGamma(float64(n) + 1)
}

// LogChoose returns the natural logarithm of the binomial coefficient
// ln C(n,k). It returns math.Inf(-1) when k is outside [0,n].
func LogChoose(n, k int) float64 {
	if k < 0 || k > n {
		return math.Inf(-1)
	}
	return LogFactorial(n) - LogFactorial(k) - LogFactorial(n-k)
}

// LogMultinomialCoefficient returns the natural logarithm of the multinomial
// coefficient (Σkᵢ)! / Πkᵢ! for the non-negative counts ks.
func LogMultinomialCoefficient(ks []int) float64 {
	total := 0
	den := 0.0
	for _, k := range ks {
		if k < 0 {
			return math.NaN()
		}
		total += k
		den += LogFactorial(k)
	}
	return LogFactorial(total) - den
}

// LogSumExp returns ln(Σ exp(xᵢ)) computed in a numerically stable way by
// factoring out the maximum element. It returns math.Inf(-1) for an empty slice.
func LogSumExp(xs []float64) float64 {
	if len(xs) == 0 {
		return math.Inf(-1)
	}
	max := math.Inf(-1)
	for _, x := range xs {
		if x > max {
			max = x
		}
	}
	if math.IsInf(max, -1) {
		return math.Inf(-1)
	}
	var sum float64
	for _, x := range xs {
		sum += math.Exp(x - max)
	}
	return max + math.Log(sum)
}

// Digamma returns the digamma function ψ(x) = d/dx ln Γ(x). It uses the
// recurrence ψ(x) = ψ(x+1) − 1/x to shift the argument into the asymptotic
// regime and then an asymptotic expansion. It is valid for x > 0.
func Digamma(x float64) float64 {
	var result float64
	for x < 6 {
		result -= 1 / x
		x++
	}
	f := 1 / (x * x)
	result += math.Log(x) - 1/(2*x) -
		f*(1.0/12-f*(1.0/120-f*(1.0/252-f*(1.0/240-f*(1.0/132)))))
	return result
}

// Trigamma returns the trigamma function ψ₁(x) = d²/dx² ln Γ(x), the derivative
// of the digamma function, for x > 0.
func Trigamma(x float64) float64 {
	var result float64
	for x < 6 {
		result += 1 / (x * x)
		x++
	}
	f := 1 / (x * x)
	result += 1/x + f/2 + (f/x)*(1.0/6-f*(1.0/30-f*(1.0/42-f*(1.0/30))))
	return result
}

// gammaSeries evaluates the regularized lower incomplete gamma function P(a,x)
// via its series representation, appropriate for x < a+1.
func gammaSeries(a, x float64) float64 {
	if x <= 0 {
		return 0
	}
	ap := a
	sum := 1 / a
	del := sum
	for n := 0; n < 500; n++ {
		ap++
		del *= x / ap
		sum += del
		if math.Abs(del) < math.Abs(sum)*1e-15 {
			break
		}
	}
	return sum * math.Exp(-x+a*math.Log(x)-LogGamma(a))
}

// gammaContinuedFraction evaluates the regularized upper incomplete gamma
// function Q(a,x) via a continued fraction, appropriate for x >= a+1.
func gammaContinuedFraction(a, x float64) float64 {
	const tiny = 1e-300
	b := x + 1 - a
	c := 1 / tiny
	d := 1 / b
	h := d
	for i := 1; i < 500; i++ {
		an := -float64(i) * (float64(i) - a)
		b += 2
		d = an*d + b
		if math.Abs(d) < tiny {
			d = tiny
		}
		c = b + an/c
		if math.Abs(c) < tiny {
			c = tiny
		}
		d = 1 / d
		del := d * c
		h *= del
		if math.Abs(del-1) < 1e-15 {
			break
		}
	}
	return math.Exp(-x+a*math.Log(x)-LogGamma(a)) * h
}

// RegularizedGammaP returns the regularized lower incomplete gamma function
// P(a,x) = γ(a,x)/Γ(a) for a > 0 and x >= 0. It equals the CDF of a Gamma(a,1)
// random variable at x.
func RegularizedGammaP(a, x float64) float64 {
	if x < 0 || a <= 0 {
		return math.NaN()
	}
	if x == 0 {
		return 0
	}
	if x < a+1 {
		return gammaSeries(a, x)
	}
	return 1 - gammaContinuedFraction(a, x)
}

// RegularizedGammaQ returns the regularized upper incomplete gamma function
// Q(a,x) = 1 − P(a,x) = Γ(a,x)/Γ(a) for a > 0 and x >= 0.
func RegularizedGammaQ(a, x float64) float64 {
	return 1 - RegularizedGammaP(a, x)
}

// InverseRegularizedGammaP returns the value x such that RegularizedGammaP(a,x)
// = p, for a > 0 and p in [0,1]. It uses a bracketing bisection combined with
// Newton refinement and is accurate to roughly 1e-12.
func InverseRegularizedGammaP(p, a float64) float64 {
	if a <= 0 || p < 0 || p > 1 {
		return math.NaN()
	}
	if p == 0 {
		return 0
	}
	if p == 1 {
		return math.Inf(1)
	}
	// Initial guess (Wilson–Hilferty).
	x := a
	{
		g := LogGamma(a)
		if a > 1 {
			pp := p
			if pp > 0.5 {
				pp = 1 - pp
			}
			t := math.Sqrt(-2 * math.Log(pp))
			xx := (2.30753 + t*0.27061) / (1 + t*(0.99229+t*0.04481))
			if p < 0.5 {
				xx = -xx
			}
			x = a * math.Pow(1-1/(9*a)-xx/(3*math.Sqrt(a)), 3)
			if x <= 0 {
				x = a
			}
		} else {
			t := 1 - a*(0.253+a*0.12)
			if p < t {
				x = math.Pow(p/t, 1/a)
			} else {
				x = 1 - math.Log(1-(p-t)/(1-t))
			}
		}
		_ = g
	}
	// Newton with a bisection fallback.
	lo, hi := 0.0, math.Inf(1)
	for i := 0; i < 200; i++ {
		f := RegularizedGammaP(a, x) - p
		if f > 0 {
			hi = x
		} else {
			lo = x
		}
		// derivative of P wrt x = pdf of Gamma(a,1)
		d := math.Exp(-x + (a-1)*math.Log(x) - LogGamma(a))
		var step float64
		if d > 0 {
			step = f / d
		}
		next := x - step
		if next <= lo || next >= hi || d == 0 {
			if math.IsInf(hi, 1) {
				next = 2 * x
			} else {
				next = 0.5 * (lo + hi)
			}
		}
		if math.Abs(next-x) < 1e-14*math.Abs(next)+1e-300 {
			return next
		}
		x = next
	}
	return x
}

// betaContinuedFraction evaluates the continued fraction used by the
// regularized incomplete beta function.
func betaContinuedFraction(x, a, b float64) float64 {
	const tiny = 1e-300
	qab := a + b
	qap := a + 1
	qam := a - 1
	c := 1.0
	d := 1 - qab*x/qap
	if math.Abs(d) < tiny {
		d = tiny
	}
	d = 1 / d
	h := d
	for m := 1; m < 500; m++ {
		fm := float64(m)
		m2 := 2 * fm
		aa := fm * (b - fm) * x / ((qam + m2) * (a + m2))
		d = 1 + aa*d
		if math.Abs(d) < tiny {
			d = tiny
		}
		c = 1 + aa/c
		if math.Abs(c) < tiny {
			c = tiny
		}
		d = 1 / d
		h *= d * c
		aa = -(a + fm) * (qab + fm) * x / ((a + m2) * (qap + m2))
		d = 1 + aa*d
		if math.Abs(d) < tiny {
			d = tiny
		}
		c = 1 + aa/c
		if math.Abs(c) < tiny {
			c = tiny
		}
		d = 1 / d
		del := d * c
		h *= del
		if math.Abs(del-1) < 1e-15 {
			break
		}
	}
	return h
}

// RegularizedIncompleteBeta returns the regularized incomplete beta function
// Iₓ(a,b) for x in [0,1] and a,b > 0. It equals the CDF of a Beta(a,b) random
// variable evaluated at x.
func RegularizedIncompleteBeta(x, a, b float64) float64 {
	if a <= 0 || b <= 0 {
		return math.NaN()
	}
	if x <= 0 {
		return 0
	}
	if x >= 1 {
		return 1
	}
	front := math.Exp(a*math.Log(x) + b*math.Log(1-x) - LogBeta(a, b))
	if x < (a+1)/(a+b+2) {
		return front * betaContinuedFraction(x, a, b) / a
	}
	return 1 - front*betaContinuedFraction(1-x, b, a)/b
}

// InverseRegularizedIncompleteBeta returns the value x in [0,1] such that
// RegularizedIncompleteBeta(x,a,b) = p, i.e. the quantile of the Beta(a,b)
// distribution. It uses bisection with Newton refinement.
func InverseRegularizedIncompleteBeta(p, a, b float64) float64 {
	if a <= 0 || b <= 0 || p < 0 || p > 1 {
		return math.NaN()
	}
	if p == 0 {
		return 0
	}
	if p == 1 {
		return 1
	}
	lo, hi := 0.0, 1.0
	x := 0.5
	lb := LogBeta(a, b)
	for i := 0; i < 200; i++ {
		f := RegularizedIncompleteBeta(x, a, b) - p
		if f > 0 {
			hi = x
		} else {
			lo = x
		}
		// pdf of Beta(a,b) at x
		var d float64
		if x > 0 && x < 1 {
			d = math.Exp((a-1)*math.Log(x) + (b-1)*math.Log(1-x) - lb)
		}
		next := x
		if d > 0 {
			next = x - f/d
		}
		if next <= lo || next >= hi {
			next = 0.5 * (lo + hi)
		}
		if math.Abs(next-x) < 1e-15 {
			return next
		}
		x = next
	}
	return x
}

// ErfInv returns the inverse error function, the value y such that erf(y) = x
// for x in (−1,1). It uses a rational approximation refined by a Newton step.
func ErfInv(x float64) float64 {
	if x <= -1 {
		return math.Inf(-1)
	}
	if x >= 1 {
		return math.Inf(1)
	}
	if x == 0 {
		return 0
	}
	// Winitzki approximation as the seed.
	a := 0.147
	ln := math.Log(1 - x*x)
	t := 2/(math.Pi*a) + ln/2
	y := math.Copysign(math.Sqrt(math.Sqrt(t*t-ln/a)-t), x)
	// One Newton step: f(y) = erf(y) − x, f'(y) = 2/√π e^{−y²}.
	for i := 0; i < 3; i++ {
		e := math.Erf(y) - x
		y -= e / (2 / math.Sqrt(math.Pi) * math.Exp(-y*y))
	}
	return y
}

// StdNormalQuantile returns the p-quantile of the standard normal distribution,
// the value z such that Φ(z) = p for p in (0,1).
func StdNormalQuantile(p float64) float64 {
	if p <= 0 {
		return math.Inf(-1)
	}
	if p >= 1 {
		return math.Inf(1)
	}
	return math.Sqrt2 * ErfInv(2*p-1)
}

// StdNormalCDF returns the standard normal cumulative distribution function
// Φ(z) = ½(1 + erf(z/√2)).
func StdNormalCDF(z float64) float64 {
	return 0.5 * math.Erfc(-z/math.Sqrt2)
}

// StdNormalPDF returns the standard normal probability density φ(z).
func StdNormalPDF(z float64) float64 {
	return math.Exp(-0.5*z*z) / math.Sqrt(2*math.Pi)
}
