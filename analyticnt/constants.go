package analyticnt

import "math"

// Analytic number theory constants. Values are given to double precision.
const (
	// MertensConstant is the Meissel–Mertens constant M ≈ 0.2614972128,
	// appearing in Σ_{p≤x} 1/p = ln ln x + M + o(1).
	MertensConstant = 0.26149721284764278375542683860869585905

	// TwinPrimeConstant is Hardy–Littlewood's C₂ ≈ 0.6601618158, controlling
	// the density of twin primes.
	TwinPrimeConstant = 0.66016181584686957392781211001455577843

	// BrunConstant is Brun's constant B₂ ≈ 1.902160583, the (conjectured)
	// convergent sum of reciprocals of twin primes.
	BrunConstant = 1.902160583104

	// ArtinConstant is Artin's constant ≈ 0.3739558136, the density of primes
	// with a given primitive root under Artin's conjecture.
	ArtinConstant = 0.37395581361920228805472805434641641511162924

	// LandauRamanujanConstant ≈ 0.7642236535, governing the density of sums of
	// two squares.
	LandauRamanujanConstant = 0.76422365358922066299069873125009232811679

	// FellerTornierConstant ≈ 0.6613170495.
	FellerTornierConstant = 0.66131704946962233528766040608130611265

	// SarnakConstant (Hafner–Sarnak–McCurley) ≈ 0.3532363719.
	SarnakConstant = 0.35323637185499598454351655043268201581

	// GlaisherKinkelin is Glaisher's constant A ≈ 1.2824271291, related to
	// ζ'(−1).
	GlaisherKinkelin = 1.28242712910062263687534256886979172776768893

	// CatalanConstant is Catalan's constant G ≈ 0.9159655942.
	CatalanConstant = 0.91596559417721901505460351493238411077414937

	// StieltjesGamma0 is the 0-th Stieltjes constant, equal to Euler's γ.
	StieltjesGamma0 = EulerGamma

	// StieltjesGamma1 is the 1-st Stieltjes constant γ₁ ≈ −0.0728158454.
	StieltjesGamma1 = -0.072815845483676724860586375874901319138

	// StieltjesGamma2 is the 2-nd Stieltjes constant γ₂ ≈ −0.0096903632.
	StieltjesGamma2 = -0.0096903631928723184845303860352125293590

	// LegendreConstant is Legendre's empirical constant B ≈ 1.08366 used in the
	// estimate x/(ln x − B).
	LegendreConstant = 1.08366
)

// EulerGammaExp returns e^{−γ} ≈ 0.5614594836, the constant appearing in
// Mertens' third theorem.
func EulerGammaExp() float64 { return math.Exp(-EulerGamma) }

// TwinPrimeDensity returns the Hardy–Littlewood estimate of the number of twin
// prime pairs (p, p+2) with p ≤ x, namely 2·C₂·∫_2^x dt/(ln t)². It is
// evaluated with an offset-Li–like construction.
func TwinPrimeDensity(x float64) float64 {
	if x < 3 {
		return 0
	}
	return 2 * TwinPrimeConstant * liSquared(x)
}

// liSquared approximates ∫_2^x dt/(ln t)² by Simpson's rule.
func liSquared(x float64) float64 {
	const steps = 4000
	a := 2.0
	b := x
	h := (b - a) / steps
	f := func(t float64) float64 {
		l := math.Log(t)
		return 1 / (l * l)
	}
	sum := f(a) + f(b)
	for i := 1; i < steps; i++ {
		t := a + float64(i)*h
		if i%2 == 1 {
			sum += 4 * f(t)
		} else {
			sum += 2 * f(t)
		}
	}
	return sum * h / 3
}

// ZetaAtEvenInteger returns the exact value ζ(2k) = (−1)^{k+1} B_{2k} (2π)^{2k}
// / (2·(2k)!) for k >= 1, using tabulated Bernoulli numbers. For example ζ(2) =
// π²/6.
func ZetaAtEvenInteger(k int) float64 {
	if k < 1 {
		panic("analyticnt: ZetaAtEvenInteger requires k >= 1")
	}
	switch k {
	case 1:
		return math.Pi * math.Pi / 6
	case 2:
		return math.Pow(math.Pi, 4) / 90
	case 3:
		return math.Pow(math.Pi, 6) / 945
	case 4:
		return math.Pow(math.Pi, 8) / 9450
	default:
		return Zeta(float64(2 * k))
	}
}

// AperyConstant is ζ(3) ≈ 1.2020569032, Apéry's constant.
const AperyConstant = 1.20205690315959428539973816151144999076498629

// PrimeZeta returns the prime zeta function P(s) = Σ_p p^{-s} for real s > 1,
// computed by Möbius inversion P(s) = Σ_{n≥1} μ(n)/n · ln ζ(ns).
func PrimeZeta(s float64) float64 {
	if s <= 1 {
		panic("analyticnt: PrimeZeta requires s > 1")
	}
	sum := 0.0
	for n := 1; n <= 200; n++ {
		mu := MobiusMu(int64(n))
		if mu == 0 {
			continue
		}
		z := Zeta(s * float64(n))
		delta := float64(mu) / float64(n) * math.Log(z)
		sum += delta
		if math.Abs(delta) < 1e-18 && n > 2 {
			break
		}
	}
	return sum
}
