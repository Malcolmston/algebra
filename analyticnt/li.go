package analyticnt

import "math"

// EulerGamma is the Euler–Mascheroni constant γ.
const EulerGamma = 0.57721566490153286060651209008240243104215933593992

// E1 returns the exponential integral E1(x) = ∫_x^∞ e^{-t}/t dt for x > 0.
// It uses a power series for small x and a continued fraction for large x.
func E1(x float64) float64 {
	if x <= 0 {
		panic("analyticnt: E1 requires x > 0")
	}
	if x < 1 {
		// Series: E1(x) = -γ - ln x - Σ_{k≥1} (-1)^k x^k /(k·k!).
		sum := 0.0
		term := 1.0
		for k := 1; k <= 100; k++ {
			term *= -x / float64(k)
			delta := term / float64(k)
			sum += delta
			if math.Abs(delta) < 1e-18*math.Abs(sum) {
				break
			}
		}
		return -EulerGamma - math.Log(x) - sum
	}
	// Lentz continued fraction: E1(x) = e^{-x} · 1/(x+1-1²/(x+3-2²/(x+5-…))).
	const tiny = 1e-300
	b := x + 1
	c := 1 / tiny
	d := 1 / b
	f := d
	for i := 1; i <= 200; i++ {
		fi := float64(i)
		a := -fi * fi
		b += 2
		d = a*d + b
		if math.Abs(d) < tiny {
			d = tiny
		}
		c = b + a/c
		if math.Abs(c) < tiny {
			c = tiny
		}
		d = 1 / d
		delta := c * d
		f *= delta
		if math.Abs(delta-1) < 1e-17 {
			break
		}
	}
	return f * math.Exp(-x)
}

// Ei returns the exponential integral Ei(x) = -PV ∫_{-x}^∞ e^{-t}/t dt, the
// principal value for x > 0 and the ordinary integral for x < 0. Ei has a
// logarithmic singularity at x = 0.
func Ei(x float64) float64 {
	if x == 0 {
		return math.Inf(-1)
	}
	if x < 0 {
		return -E1(-x)
	}
	if x < 40 {
		// Power series: Ei(x) = γ + ln x + Σ_{k≥1} x^k/(k·k!).
		sum := 0.0
		term := 1.0
		for k := 1; k <= 500; k++ {
			term *= x / float64(k)
			delta := term / float64(k)
			sum += delta
			if math.Abs(delta) < 1e-18*math.Abs(sum) {
				break
			}
		}
		return EulerGamma + math.Log(x) + sum
	}
	// Asymptotic series: Ei(x) ~ e^x/x · Σ_{k≥0} k!/x^k, truncated at the
	// smallest term.
	sum := 1.0
	term := 1.0
	prev := math.Inf(1)
	for k := 1; k <= 60; k++ {
		term *= float64(k) / x
		if term > prev {
			break
		}
		prev = term
		sum += term
	}
	return math.Exp(x) / x * sum
}

// Li returns the logarithmic integral li(x) = PV ∫_0^x dt/ln t for x > 0, x ≠ 1.
// It is computed as Ei(ln x). li(1) is a singularity and li near 1 is large in
// magnitude.
func Li(x float64) float64 {
	if x <= 0 {
		panic("analyticnt: Li requires x > 0")
	}
	if x == 1 {
		return math.Inf(-1)
	}
	return Ei(math.Log(x))
}

// SoldnerConstant is the Ramanujan–Soldner constant μ ≈ 1.4513692349, the
// unique positive root of the logarithmic integral li(μ) = 0.
const SoldnerConstant = 1.45136923488338105028396848589202744949303228

// LiOffset returns the offset logarithmic integral Li(x) = li(x) − li(2) =
// ∫_2^x dt/ln t. This is the standard smooth estimate for the prime-counting
// function π(x). x must be >= 2.
func LiOffset(x float64) float64 {
	if x < 2 {
		panic("analyticnt: LiOffset requires x >= 2")
	}
	if x == 2 {
		return 0
	}
	return Li(x) - liOf2
}

// liOf2 is li(2), the additive constant relating li and the offset Li.
var liOf2 = Ei(math.Log(2))

// LiValueAt2 returns li(2) ≈ 1.045163780117493.
func LiValueAt2() float64 { return liOf2 }
