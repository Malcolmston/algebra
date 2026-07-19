package analyticnt

import "math"

// PiError returns the signed error li(x) − π(x) of the offset-logarithmic
// integral estimate of the prime-counting function at integer x. For all
// tabulated x this is positive (the first sign change occurs near the enormous
// Skewes number).
func PiError(x int64) float64 {
	return LiOffset(float64(x)) - float64(PrimePi(x))
}

// PiErrorRiemann returns the signed error R(x) − π(x) of the Riemann estimate.
func PiErrorRiemann(x int64) float64 {
	return RiemannR(float64(x)) - float64(PrimePi(x))
}

// PsiError returns the normalized error (ψ(x) − x)/√x of the second Chebyshev
// function, the quantity whose boundedness is equivalent to the Riemann
// hypothesis.
func PsiError(x float64) float64 {
	return (ChebyshevPsi(x) - x) / math.Sqrt(x)
}

// ThetaError returns the normalized error (θ(x) − x)/√x of the first Chebyshev
// function.
func ThetaError(x float64) float64 {
	return (ChebyshevTheta(x) - x) / math.Sqrt(x)
}

// PrimeCountingBias returns li(x) − π(x) as a float and reports whether the
// value is positive (li over-estimates π), illustrating the Chebyshev bias.
func PrimeCountingBias(x int64) (bias float64, liOverestimates bool) {
	bias = PiError(x)
	return bias, bias > 0
}

// RiemannHypothesisBound returns the RH-conditional error bound for |π(x) −
// li(x)|, namely (1/(8π))·√x·ln x, which holds for x >= 2657 under the Riemann
// hypothesis (Schoenfeld).
func RiemannHypothesisBound(x float64) float64 {
	return math.Sqrt(x) * math.Log(x) / (8 * math.Pi)
}

// PsiSchoenfeldBound returns Schoenfeld's RH-conditional bound |ψ(x) − x| <
// (1/(8π))·√x·(ln x)², valid for x >= 73.2.
func PsiSchoenfeldBound(x float64) float64 {
	l := math.Log(x)
	return math.Sqrt(x) * l * l / (8 * math.Pi)
}

// PrimeDensity returns the prime density 1/ln x, the probability that a number
// near x is prime according to the prime number theorem.
func PrimeDensity(x float64) float64 {
	if x <= 1 {
		panic("analyticnt: PrimeDensity requires x > 1")
	}
	return 1 / math.Log(x)
}

// ExplicitPsi returns Riemann's explicit-formula approximation to ψ(x) using the
// first numZeros nontrivial zeros ρ = 1/2 + iγ:
// ψ(x) ≈ x − Σ_ρ x^ρ/ρ − ln(2π) − ½ln(1 − x^{-2}). The oscillatory zero sum is
// taken over conjugate pairs so the result is real.
func ExplicitPsi(x float64, numZeros int) float64 {
	if x <= 1 {
		panic("analyticnt: ExplicitPsi requires x > 1")
	}
	sum := x - math.Log(2*math.Pi) - 0.5*math.Log(1-1/(x*x))
	gammas := ZetaZeros(numZeros)
	lnx := math.Log(x)
	for _, g := range gammas {
		// x^ρ/ρ + x^{ρ̄}/ρ̄ = 2·Re(x^ρ/ρ); ρ = 1/2 + iγ.
		// x^ρ = √x·(cos(γ ln x) + i sin(γ ln x)); 1/ρ = (1/2 − iγ)/(1/4+γ²).
		sx := math.Sqrt(x)
		cosT := math.Cos(g * lnx)
		sinT := math.Sin(g * lnx)
		denom := 0.25 + g*g
		// Re( (cosT + i sinT) · (1/2 − iγ) ) / denom
		re := (cosT*0.5 + sinT*g) / denom
		sum -= 2 * sx * re
	}
	return sum
}

// SkewesNumberLog returns the natural logarithm of the classical Skewes number,
// the (very large) bound below which li(x) − π(x) was proven to change sign; it
// is included as a documented constant, ln(Sk) ≈ 7.7·10^{2}? Actually returns
// the widely cited exponent e^{e^{e^{79}}} is intractable, so this returns the
// modern Bays–Hudson crossover estimate ≈ 1.397·10^{316}.
func SkewesNumberLog() float64 {
	// ln(1.397e316) = ln 1.397 + 316 ln 10.
	return math.Log(1.397) + 316*math.Log(10)
}

// CramerGap returns Cramér's heuristic estimate (ln p)² for the gap following a
// prime p, the conjectured order of the maximal prime gap near p.
func CramerGap(p float64) float64 {
	l := math.Log(p)
	return l * l
}
