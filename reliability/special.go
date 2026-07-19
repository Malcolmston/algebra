package reliability

import "math"

// normalCDF returns the cumulative distribution function of the standard
// normal distribution evaluated at x.
func normalCDF(x float64) float64 {
	return 0.5 * math.Erfc(-x/math.Sqrt2)
}

// normalPDF returns the probability density of the standard normal
// distribution evaluated at x.
func normalPDF(x float64) float64 {
	return math.Exp(-0.5*x*x) / math.Sqrt(2*math.Pi)
}

// normalQuantile returns the inverse of the standard-normal CDF (the probit
// function) for probability p in (0,1), using Acklam's rational approximation
// refined by one Halley step. It returns -Inf at 0 and +Inf at 1 and NaN
// outside [0,1].
func normalQuantile(p float64) float64 {
	if math.IsNaN(p) || p < 0 || p > 1 {
		return math.NaN()
	}
	if p == 0 {
		return math.Inf(-1)
	}
	if p == 1 {
		return math.Inf(1)
	}
	a := [6]float64{-3.969683028665376e+01, 2.209460984245205e+02,
		-2.759285104469687e+02, 1.383577518672690e+02,
		-3.066479806614716e+01, 2.506628277459239e+00}
	b := [5]float64{-5.447609879822406e+01, 1.615858368580409e+02,
		-1.556989798598866e+02, 6.680131188771972e+01,
		-1.328068155288572e+01}
	c := [6]float64{-7.784894002430293e-03, -3.223964580411365e-01,
		-2.400758277161838e+00, -2.549732539343734e+00,
		4.374664141464968e+00, 2.938163982698783e+00}
	d := [4]float64{7.784695709041462e-03, 3.224671290700398e-01,
		2.445134137142996e+00, 3.754408661907416e+00}
	const plow = 0.02425
	const phigh = 1 - plow
	var x float64
	switch {
	case p < plow:
		q := math.Sqrt(-2 * math.Log(p))
		x = (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q + c[5]) /
			((((d[0]*q+d[1])*q+d[2])*q+d[3])*q + 1)
	case p <= phigh:
		q := p - 0.5
		r := q * q
		x = (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r + a[5]) * q /
			(((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r + 1)
	default:
		q := math.Sqrt(-2 * math.Log(1-p))
		x = -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q + c[5]) /
			((((d[0]*q+d[1])*q+d[2])*q+d[3])*q + 1)
	}
	// One Halley refinement step for full double precision.
	e := normalCDF(x) - p
	u := e * math.Sqrt(2*math.Pi) * math.Exp(0.5*x*x)
	x = x - u/(1+0.5*x*u)
	return x
}

// regGammaP returns the regularized lower incomplete gamma function
// P(a,x) = γ(a,x)/Γ(a) for a>0 and x>=0.
func regGammaP(a, x float64) float64 {
	if x < 0 || a <= 0 || math.IsNaN(x) || math.IsNaN(a) {
		return math.NaN()
	}
	if x == 0 {
		return 0
	}
	if x < a+1 {
		return gammaSeries(a, x)
	}
	return 1 - gammaContFrac(a, x)
}

// regGammaQ returns the regularized upper incomplete gamma function
// Q(a,x) = Γ(a,x)/Γ(a) = 1 - P(a,x).
func regGammaQ(a, x float64) float64 {
	if x < 0 || a <= 0 || math.IsNaN(x) || math.IsNaN(a) {
		return math.NaN()
	}
	if x == 0 {
		return 1
	}
	if x < a+1 {
		return 1 - gammaSeries(a, x)
	}
	return gammaContFrac(a, x)
}

// gammaSeries evaluates P(a,x) via its power series, valid for x < a+1.
func gammaSeries(a, x float64) float64 {
	const maxIter = 500
	const eps = 1e-15
	ap := a
	sum := 1.0 / a
	del := sum
	for i := 0; i < maxIter; i++ {
		ap++
		del *= x / ap
		sum += del
		if math.Abs(del) < math.Abs(sum)*eps {
			break
		}
	}
	return sum * math.Exp(-x+a*math.Log(x)-lgamma(a))
}

// gammaContFrac evaluates Q(a,x) via a continued fraction, valid for x >= a+1.
func gammaContFrac(a, x float64) float64 {
	const maxIter = 500
	const eps = 1e-15
	const fpmin = 1e-300
	b := x + 1 - a
	c := 1.0 / fpmin
	d := 1.0 / b
	h := d
	for i := 1; i <= maxIter; i++ {
		an := -float64(i) * (float64(i) - a)
		b += 2
		d = an*d + b
		if math.Abs(d) < fpmin {
			d = fpmin
		}
		c = b + an/c
		if math.Abs(c) < fpmin {
			c = fpmin
		}
		d = 1.0 / d
		del := d * c
		h *= del
		if math.Abs(del-1) < eps {
			break
		}
	}
	return math.Exp(-x+a*math.Log(x)-lgamma(a)) * h
}

// lgamma returns the natural logarithm of the absolute value of the gamma
// function.
func lgamma(x float64) float64 {
	v, _ := math.Lgamma(x)
	return v
}

// invRegGammaP returns x such that regGammaP(a,x)=p, using bisection with a
// Wilson–Hilferty initial bracket. It is used for gamma-distribution
// quantiles.
func invRegGammaP(a, p float64) float64 {
	if math.IsNaN(a) || math.IsNaN(p) || a <= 0 || p < 0 || p > 1 {
		return math.NaN()
	}
	if p == 0 {
		return 0
	}
	if p == 1 {
		return math.Inf(1)
	}
	// Wilson–Hilferty approximation for an initial guess.
	z := normalQuantile(p)
	g := a * math.Pow(1-1.0/(9*a)+z/(3*math.Sqrt(a)), 3)
	if g <= 0 || math.IsNaN(g) {
		g = a
	}
	lo, hi := 0.0, g
	// Expand hi until it brackets p.
	for regGammaP(a, hi) < p {
		lo = hi
		hi *= 2
		if hi > 1e300 {
			return math.Inf(1)
		}
	}
	for i := 0; i < 200; i++ {
		mid := 0.5 * (lo + hi)
		if regGammaP(a, mid) < p {
			lo = mid
		} else {
			hi = mid
		}
		if hi-lo < 1e-14*math.Max(1, hi) {
			break
		}
	}
	return 0.5 * (lo + hi)
}

// simpson integrates f over [a,b] using the composite Simpson rule with n
// (rounded up to an even number) subintervals.
func simpson(f func(float64) float64, a, b float64, n int) float64 {
	if n < 2 {
		n = 2
	}
	if n%2 != 0 {
		n++
	}
	h := (b - a) / float64(n)
	sum := f(a) + f(b)
	for i := 1; i < n; i++ {
		x := a + float64(i)*h
		if i%2 == 0 {
			sum += 2 * f(x)
		} else {
			sum += 4 * f(x)
		}
	}
	return sum * h / 3
}

// integrateToInfinity approximates the integral of a non-negative,
// eventually-decaying reliability-like function f from a to +infinity. It
// integrates over successive panels of width step until the running increment
// becomes negligible relative to the accumulated area.
func integrateToInfinity(f func(float64) float64, a, step float64) float64 {
	if step <= 0 {
		step = 1
	}
	total := 0.0
	x := a
	for panel := 0; panel < 100000; panel++ {
		inc := simpson(f, x, x+step, 16)
		total += inc
		x += step
		if inc <= 1e-14*math.Max(1, total) && panel > 4 {
			break
		}
	}
	return total
}
