package analyticnt

import "math"

// RiemannR returns the Riemann prime-counting approximation R(x), evaluated via
// the Gram series R(x) = 1 + Σ_{k≥1} (ln x)^k / (k · k! · ζ(k+1)). For x >= 2 it
// is a markedly better estimate of π(x) than li(x).
func RiemannR(x float64) float64 {
	if x <= 1 {
		panic("analyticnt: RiemannR requires x > 1")
	}
	lnx := math.Log(x)
	sum := 1.0
	term := 1.0 // (ln x)^k / k! accumulated
	for k := 1; k <= 200; k++ {
		term *= lnx / float64(k)
		delta := term / (float64(k) * Zeta(float64(k+1)))
		sum += delta
		if math.Abs(delta) < 1e-16*math.Abs(sum) && k > 5 {
			break
		}
	}
	return sum
}

// RiemannRPrime returns the derivative R'(x) via the series
// R'(x) = Σ_{k≥1} (ln x)^{k-1} / (x · k! · ζ(k+1)), useful for inverting R.
func RiemannRPrime(x float64) float64 {
	if x <= 1 {
		panic("analyticnt: RiemannRPrime requires x > 1")
	}
	lnx := math.Log(x)
	sum := 0.0
	term := 1.0 // (ln x)^{k-1}/(k-1)! -> handled incrementally
	pow := 1.0  // (ln x)^{k-1}
	fact := 1.0 // k!
	for k := 1; k <= 200; k++ {
		fact *= float64(k)
		delta := pow / (fact * Zeta(float64(k+1)))
		sum += delta
		pow *= lnx
		_ = term
		if math.Abs(delta) < 1e-18*(math.Abs(sum)+1) && k > 5 {
			break
		}
	}
	return sum / x
}

// RiemannRScaled returns R(x) computed with the Möbius-inversion form
// R(x) = Σ_{n≥1} μ(n)/n · li(x^{1/n}). It agrees with the Gram-series RiemannR
// but exposes the alternative definition.
func RiemannRScaled(x float64) float64 {
	if x <= 1 {
		panic("analyticnt: RiemannRScaled requires x > 1")
	}
	sum := 0.0
	for n := 1; n <= 60; n++ {
		mu := MobiusMu(int64(n))
		if mu == 0 {
			continue
		}
		root := math.Pow(x, 1.0/float64(n))
		if root <= 1.0000001 {
			break
		}
		sum += float64(mu) / float64(n) * Li(root)
	}
	return sum
}

// PrimePiApprox returns the recommended real-valued estimate of π(x), namely
// the Riemann R function for x >= 2 and 0 below that.
func PrimePiApprox(x float64) float64 {
	if x < 2 {
		return 0
	}
	return RiemannR(x)
}

// PrimePiLi returns the offset logarithmic-integral estimate Li(x) of π(x) for
// x >= 2.
func PrimePiLi(x float64) float64 {
	if x < 2 {
		return 0
	}
	return LiOffset(x)
}

// PrimePiLegendreApprox returns Legendre's classical estimate x/(ln x − B) with
// B = 1.08366 (Legendre's constant). It is historically important but less
// accurate than li or R.
func PrimePiLegendreApprox(x float64) float64 {
	if x < 2 {
		return 0
	}
	return x / (math.Log(x) - 1.08366)
}

// PrimePiGauss returns Gauss's original estimate x/ln x of π(x).
func PrimePiGauss(x float64) float64 {
	if x < 2 {
		return 0
	}
	return x / math.Log(x)
}

// NthPrimeApprox returns an asymptotic estimate of the n-th prime using the
// expansion p_n ≈ n(ln n + ln ln n − 1 + (ln ln n − 2)/ln n). n must be >= 6
// for the expansion to be meaningful; smaller n fall back to a coarse bound.
func NthPrimeApprox(n int) float64 {
	if n < 1 {
		panic("analyticnt: NthPrimeApprox requires n >= 1")
	}
	fn := float64(n)
	if n < 6 {
		return fn * (math.Log(fn) + math.Log(math.Log(fn+2)))
	}
	l := math.Log(fn)
	ll := math.Log(l)
	return fn * (l + ll - 1 + (ll-2)/l - (ll*ll-6*ll+11)/(2*l*l))
}

// NthPrimeInverseR returns the n-th prime estimated by inverting the Riemann R
// function, i.e. solving R(x) = n by Newton iteration. This is one of the most
// accurate closed-form estimators available.
func NthPrimeInverseR(n int) float64 {
	if n < 1 {
		panic("analyticnt: NthPrimeInverseR requires n >= 1")
	}
	x := NthPrimeApprox(n)
	if x < 2 {
		x = 2.5
	}
	for i := 0; i < 100; i++ {
		f := RiemannR(x) - float64(n)
		fp := RiemannRPrime(x)
		if fp == 0 {
			break
		}
		nx := x - f/fp
		if math.Abs(nx-x) < 1e-6*x {
			x = nx
			break
		}
		x = nx
	}
	return x
}
