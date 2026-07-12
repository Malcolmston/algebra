package stats

import "math"

// gammaLn returns the natural logarithm of the absolute value of the gamma
// function, ln|Γ(x)|. It is a thin wrapper around math.Lgamma used throughout
// the package to keep factorials, binomial coefficients and the incomplete
// gamma/beta functions finite and accurate for large arguments.
func gammaLn(x float64) float64 {
	v, _ := math.Lgamma(x)
	return v
}

// regularizedGammaP returns the regularized lower incomplete gamma function
// P(a, x) = γ(a, x) / Γ(a), for a > 0 and x >= 0. It is the CDF building
// block for the gamma-family distributions (Gamma, ChiSquared, Poisson).
func regularizedGammaP(a, x float64) float64 {
	if a <= 0 || x < 0 {
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

// regularizedGammaQ returns the regularized upper incomplete gamma function
// Q(a, x) = 1 - P(a, x).
func regularizedGammaQ(a, x float64) float64 {
	return 1 - regularizedGammaP(a, x)
}

// gammaSeries evaluates P(a, x) by its power series, accurate for x < a+1.
func gammaSeries(a, x float64) float64 {
	const eps = 3e-14
	ap := a
	sum := 1 / a
	del := sum
	for n := 0; n < 500; n++ {
		ap++
		del *= x / ap
		sum += del
		if math.Abs(del) < math.Abs(sum)*eps {
			break
		}
	}
	return sum * math.Exp(-x+a*math.Log(x)-gammaLn(a))
}

// gammaContinuedFraction evaluates Q(a, x) by its continued fraction,
// accurate for x >= a+1 (Lentz's algorithm).
func gammaContinuedFraction(a, x float64) float64 {
	const eps = 3e-14
	const fpmin = 1e-300
	b := x + 1 - a
	c := 1 / fpmin
	d := 1 / b
	h := d
	for i := 1; i < 500; i++ {
		fi := float64(i)
		an := -fi * (fi - a)
		b += 2
		d = an*d + b
		if math.Abs(d) < fpmin {
			d = fpmin
		}
		c = b + an/c
		if math.Abs(c) < fpmin {
			c = fpmin
		}
		d = 1 / d
		del := d * c
		h *= del
		if math.Abs(del-1) < eps {
			break
		}
	}
	return math.Exp(-x+a*math.Log(x)-gammaLn(a)) * h
}

// regularizedIncompleteBeta returns the regularized incomplete beta function
// I_x(a, b) for a, b > 0 and x in [0, 1]. It is the CDF building block for the
// Student's t distribution.
func regularizedIncompleteBeta(a, b, x float64) float64 {
	if x <= 0 {
		return 0
	}
	if x >= 1 {
		return 1
	}
	lnFront := gammaLn(a+b) - gammaLn(a) - gammaLn(b) +
		a*math.Log(x) + b*math.Log(1-x)
	front := math.Exp(lnFront)
	if x < (a+1)/(a+b+2) {
		return front * betaContinuedFraction(a, b, x) / a
	}
	return 1 - front*betaContinuedFraction(b, a, 1-x)/b
}

// betaContinuedFraction evaluates the continued fraction used by
// regularizedIncompleteBeta (Lentz's algorithm).
func betaContinuedFraction(a, b, x float64) float64 {
	const eps = 3e-14
	const fpmin = 1e-300
	qab := a + b
	qap := a + 1
	qam := a - 1
	c := 1.0
	d := 1 - qab*x/qap
	if math.Abs(d) < fpmin {
		d = fpmin
	}
	d = 1 / d
	h := d
	for m := 1; m < 500; m++ {
		fm := float64(m)
		aa := fm * (b - fm) * x / ((qam + 2*fm) * (a + 2*fm))
		d = 1 + aa*d
		if math.Abs(d) < fpmin {
			d = fpmin
		}
		c = 1 + aa/c
		if math.Abs(c) < fpmin {
			c = fpmin
		}
		d = 1 / d
		h *= d * c
		aa = -(a + fm) * (qab + fm) * x / ((a + 2*fm) * (qap + 2*fm))
		d = 1 + aa*d
		if math.Abs(d) < fpmin {
			d = fpmin
		}
		c = 1 + aa/c
		if math.Abs(c) < fpmin {
			c = fpmin
		}
		d = 1 / d
		del := d * c
		h *= del
		if math.Abs(del-1) < eps {
			break
		}
	}
	return h
}

// normQuantile returns the quantile (inverse CDF) of the standard normal
// distribution for p in (0, 1), using Acklam's rational approximation
// refined by one Halley step. It underlies Normal.Quantile.
func normQuantile(p float64) float64 {
	if math.IsNaN(p) || p < 0 || p > 1 {
		return math.NaN()
	}
	if p == 0 {
		return math.Inf(-1)
	}
	if p == 1 {
		return math.Inf(1)
	}
	a := [...]float64{-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00}
	b := [...]float64{-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02, 6.680131188771972e+01, -1.328068155288572e+01}
	c := [...]float64{-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00}
	d := [...]float64{7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00}
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
	e := 0.5*math.Erfc(-x/math.Sqrt2) - p
	u := e * math.Sqrt(2*math.Pi) * math.Exp(x*x/2)
	x = x - u/(1+x*u/2)
	return x
}
