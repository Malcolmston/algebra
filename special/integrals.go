package special

import (
	"math"
	"math/cmplx"
)

// The functions in this file implement a broad collection of special
// mathematical functions using only the Go standard library. All routines
// aim for correctness against known closed-form and reference values.

// specialSqrtPi is √π. The Euler–Mascheroni constant specialEulerGamma is
// declared once in bessel.go and shared across the package.
const specialSqrtPi = 1.7724538509055160272981674833411451827975 // sqrt(π)

// -----------------------------------------------------------------------------
// Error function family
// -----------------------------------------------------------------------------

// Erf returns the error function of x,
//
//	erf(x) = (2/√π) ∫₀ˣ e^{-t²} dt.
//
// It delegates to the standard library implementation.
func Erf(x float64) float64 {
	return math.Erf(x)
}

// Erfc returns the complementary error function of x, erfc(x) = 1 - erf(x).
// It is accurate even when erf(x) is very close to 1.
func Erfc(x float64) float64 {
	return math.Erfc(x)
}

// Erfcx returns the scaled complementary error function,
//
//	erfcx(x) = e^{x²} · erfc(x).
//
// The scaling removes the exponential decay of erfc for large positive x,
// avoiding underflow.
func Erfcx(x float64) float64 {
	if x < 0 {
		return 2*math.Exp(x*x) - Erfcx(-x)
	}
	if x < 25 {
		return math.Exp(x*x) * math.Erfc(x)
	}
	// Asymptotic continued fraction for large x: erfcx(x) ≈ 1/(x√π) · (1 - 1/(2x²) + ...)
	inv := 1.0 / (x * specialSqrtPi)
	x2 := x * x
	return inv * (1 - 0.5/x2 + 0.75/(x2*x2) - 1.875/(x2*x2*x2))
}

// Erfi returns the imaginary error function,
//
//	erfi(x) = -i·erf(i·x) = (2/√π) ∫₀ˣ e^{t²} dt.
//
// It grows super-exponentially and is related to Dawson's function by
// erfi(x) = (2/√π)·e^{x²}·D(x).
func Erfi(x float64) float64 {
	return 2 / specialSqrtPi * math.Exp(x*x) * Dawson(x)
}

// InverseErf returns the inverse error function, the value y such that
// erf(y) = x for x in the open interval (-1, 1). It returns ±Inf at ±1 and
// NaN outside [-1, 1].
func InverseErf(x float64) float64 {
	if x <= -1 {
		if x == -1 {
			return math.Inf(-1)
		}
		return math.NaN()
	}
	if x >= 1 {
		if x == 1 {
			return math.Inf(1)
		}
		return math.NaN()
	}
	if x == 0 {
		return 0
	}
	// Winitzki initial approximation followed by Halley refinement.
	w := -math.Log((1 - x) * (1 + x))
	var y float64
	if w < 5 {
		w -= 2.5
		y = 2.81022636e-08
		y = 3.43273939e-07 + y*w
		y = -3.5233877e-06 + y*w
		y = -4.39150654e-06 + y*w
		y = 0.00021858087 + y*w
		y = -0.00125372503 + y*w
		y = -0.00417768164 + y*w
		y = 0.246640727 + y*w
		y = 1.50140941 + y*w
	} else {
		w = math.Sqrt(w) - 3
		y = -0.000200214257
		y = 0.000100950558 + y*w
		y = 0.00134934322 + y*w
		y = -0.00367342844 + y*w
		y = 0.00573950773 + y*w
		y = -0.0076224613 + y*w
		y = 0.00943887047 + y*w
		y = 1.00167406 + y*w
		y = 2.83297682 + y*w
	}
	y *= x
	// Two Halley iterations for full double precision.
	for i := 0; i < 2; i++ {
		err := math.Erf(y) - x
		y -= err / (2/specialSqrtPi*math.Exp(-y*y) - y*err)
	}
	return y
}

// InverseErfc returns the inverse complementary error function, the value y
// such that erfc(y) = x for x in the open interval (0, 2). It returns
// ±Inf at the endpoints and NaN outside [0, 2].
func InverseErfc(x float64) float64 {
	if x <= 0 {
		if x == 0 {
			return math.Inf(1)
		}
		return math.NaN()
	}
	if x >= 2 {
		if x == 2 {
			return math.Inf(-1)
		}
		return math.NaN()
	}
	return InverseErf(1 - x)
}

// Dawson returns Dawson's integral,
//
//	D(x) = e^{-x²} ∫₀ˣ e^{t²} dt.
//
// The implementation uses a rational/series scheme accurate to near machine
// precision across the whole real line.
func Dawson(x float64) float64 {
	if x < 0 {
		return -Dawson(-x)
	}
	if x == 0 {
		return 0
	}
	if x < 4 {
		// Power series: D(x) = Σ_{n≥0} (-1)^n 2^n x^{2n+1}/(2n+1)!! .
		x2 := x * x
		term := x
		sum := x
		n := 1
		for {
			term *= -2 * x2 / float64(2*n+1)
			sum += term
			if math.Abs(term) < 1e-18*math.Abs(sum) {
				break
			}
			n++
			if n > 400 {
				break
			}
		}
		return sum
	}
	// Asymptotic expansion for large x:
	//   D(x) ~ 1/(2x)·Σ_{k≥0} (2k-1)!!/(2x²)^k .
	x2 := x * x
	inv := 1.0 / (2 * x)
	sum := 1.0
	term := 1.0
	for k := 1; k <= 40; k++ {
		prev := term
		term *= float64(2*k-1) / (2 * x2)
		if term > prev { // asymptotic series past optimal truncation
			break
		}
		sum += term
		if term < 1e-18 {
			break
		}
	}
	return inv * sum
}

// -----------------------------------------------------------------------------
// Fresnel integrals
// -----------------------------------------------------------------------------

// FresnelS returns the Fresnel sine integral,
//
//	S(x) = ∫₀ˣ sin(π t²/2) dt.
func FresnelS(x float64) float64 {
	s, _ := Fresnel(x)
	return s
}

// FresnelC returns the Fresnel cosine integral,
//
//	C(x) = ∫₀ˣ cos(π t²/2) dt.
func FresnelC(x float64) float64 {
	_, c := Fresnel(x)
	return c
}

// Fresnel returns both Fresnel integrals S(x) and C(x) simultaneously,
//
//	S(x) = ∫₀ˣ sin(π t²/2) dt,   C(x) = ∫₀ˣ cos(π t²/2) dt.
//
// A power series is used for small |x| and an asymptotic auxiliary-function
// expansion for large |x|.
func Fresnel(x float64) (s, c float64) {
	if x < 0 {
		s, c = Fresnel(-x)
		return -s, -c
	}
	if x <= 1.5 {
		// Power series in x for small arguments.
		c = specialFresnelCSeries(x)
		s = specialFresnelSSeries(x)
		return s, c
	}
	// Complex continued fraction (modified Lentz) for the auxiliary integral,
	// accurate for all x > 1.5 (after Numerical Recipes).
	const tiny = 1e-300
	pix2 := math.Pi * x * x
	b := complex(1, -pix2)
	cc := complex(1/tiny, 0)
	d := 1 / b
	h := d
	n := -1.0
	for k := 2; k <= 400; k++ {
		n += 2
		a := complex(-n*(n+1), 0)
		b += complex(4, 0)
		d = a*d + b
		if cmplx.Abs(d) < tiny {
			d = complex(tiny, 0)
		}
		cc = b + a/cc
		if cmplx.Abs(cc) < tiny {
			cc = complex(tiny, 0)
		}
		d = 1 / d
		del := d * cc
		h *= del
		if math.Abs(real(del)-1)+math.Abs(imag(del)) < 1e-16 {
			break
		}
	}
	h *= complex(x, -x)
	sinT, cosT := math.Sincos(0.5 * pix2)
	cs := complex(0.5, 0.5) * (1 - complex(cosT, sinT)*h)
	return imag(cs), real(cs)
}

// specialFresnelCSeries evaluates C(x) by its Maclaurin series.
func specialFresnelCSeries(x float64) float64 {
	piBy2 := math.Pi / 2
	x2 := x * x
	sum := 0.0
	term := x
	n := 0
	for {
		den := float64(4*n + 1)
		add := term / den
		sum += add
		if math.Abs(add) < 1e-18*math.Abs(sum)+1e-300 {
			break
		}
		n++
		term *= -(piBy2 * piBy2) * x2 * x2 / float64((2*n)*(2*n-1))
		if n > 300 {
			break
		}
	}
	return sum
}

// specialFresnelSSeries evaluates S(x) by its Maclaurin series.
func specialFresnelSSeries(x float64) float64 {
	piBy2 := math.Pi / 2
	x2 := x * x
	sum := 0.0
	term := x * (piBy2 * x2)
	n := 0
	for {
		den := float64(4*n + 3)
		add := term / den
		sum += add
		if math.Abs(add) < 1e-18*math.Abs(sum)+1e-300 {
			break
		}
		n++
		term *= -(piBy2 * piBy2) * x2 * x2 / float64((2*n+1)*(2*n))
		if n > 300 {
			break
		}
	}
	return sum
}

// -----------------------------------------------------------------------------
// Exponential, logarithmic, sine and cosine integrals
// -----------------------------------------------------------------------------

// Ei returns the exponential integral,
//
//	Ei(x) = -∫_{-x}^∞ e^{-t}/t dt   (principal value for x > 0).
//
// It is defined for all real x ≠ 0; Ei(0) is -Inf.
func Ei(x float64) float64 {
	if x == 0 {
		return math.Inf(-1)
	}
	if x < 0 {
		return -E1(-x)
	}
	if x < 40 {
		// Series (all terms positive, no cancellation):
		//   Ei(x) = γ + ln x + Σ_{k≥1} x^k/(k·k!).
		sum := 0.0
		term := 1.0
		for k := 1; k <= 500; k++ {
			term *= x / float64(k)
			add := term / float64(k)
			sum += add
			if math.Abs(add) < 1e-18*math.Abs(sum) {
				break
			}
		}
		return specialEulerGamma + math.Log(x) + sum
	}
	// Asymptotic series for large x: Ei(x) ~ e^x/x · Σ k!/x^k.
	sum := 1.0
	term := 1.0
	for k := 1; k <= 60; k++ {
		prev := term
		term *= float64(k) / x
		if term > prev { // asymptotic series past optimal truncation
			break
		}
		sum += term
		if term < 1e-18*sum {
			break
		}
	}
	return math.Exp(x) / x * sum
}

// E1 returns the exponential integral E₁,
//
//	E₁(x) = ∫₁^∞ e^{-x t}/t dt = ∫_x^∞ e^{-t}/t dt,   x > 0.
//
// For x ≤ 0 the value is defined via analytic continuation as -Ei(-x).
func E1(x float64) float64 {
	if x == 0 {
		return math.Inf(1)
	}
	if x < 0 {
		return -Ei(-x)
	}
	if x <= 1 {
		// Series: E1(x) = -γ - ln x + Σ_{k≥1} (-1)^{k+1} x^k/(k·k!)
		sum := 0.0
		term := 1.0
		for k := 1; k <= 100; k++ {
			term *= -x / float64(k)
			add := -term / float64(k)
			sum += add
			if math.Abs(add) < 1e-18*math.Abs(sum) {
				break
			}
		}
		return -specialEulerGamma - math.Log(x) + sum
	}
	// Continued fraction (Lentz) for x > 1.
	const tiny = 1e-300
	b := x + 1
	c := 1.0 / tiny
	d := 1.0 / b
	h := d
	for i := 1; i <= 200; i++ {
		a := -float64(i * i)
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
		del := c * d
		h *= del
		if math.Abs(del-1) < 1e-16 {
			break
		}
	}
	return math.Exp(-x) * h
}

// En returns the generalized exponential integral of order n,
//
//	Eₙ(x) = ∫₁^∞ e^{-x t}/tⁿ dt,   x ≥ 0, n ≥ 0.
//
// E₀(x) = e^{-x}/x and E₁ coincides with the two-argument E1.
func En(n int, x float64) float64 {
	if n < 0 || x < 0 || (x == 0 && n <= 1) {
		return math.NaN()
	}
	if n == 0 {
		return math.Exp(-x) / x
	}
	if x == 0 {
		return 1.0 / float64(n-1)
	}
	if n == 1 {
		return E1(x)
	}
	const tiny = 1e-300
	if x > 1 {
		// Continued fraction.
		b := x + float64(n)
		c := 1.0 / tiny
		d := 1.0 / b
		h := d
		for i := 1; i <= 200; i++ {
			a := -float64(i * (n + i - 1))
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
			del := c * d
			h *= del
			if math.Abs(del-1) < 1e-16 {
				break
			}
		}
		return h * math.Exp(-x)
	}
	// Series representation for x ≤ 1.
	var ans float64
	if n == 1 {
		ans = -math.Log(x) - specialEulerGamma
	} else {
		ans = 1.0 / float64(n-1)
	}
	fact := 1.0
	for i := 1; i <= 200; i++ {
		fact *= -x / float64(i)
		var del float64
		if i != n-1 {
			del = -fact / float64(i-(n-1))
		} else {
			psi := -specialEulerGamma
			for ii := 1; ii <= n-1; ii++ {
				psi += 1.0 / float64(ii)
			}
			del = fact * (-math.Log(x) + psi)
		}
		ans += del
		if math.Abs(del) < math.Abs(ans)*1e-16 {
			break
		}
	}
	return ans
}

// Li returns the logarithmic integral,
//
//	li(x) = ∫₀ˣ dt/ln t   (principal value for x > 1),
//
// which equals Ei(ln x). It is defined for x > 0, x ≠ 1.
func Li(x float64) float64 {
	if x <= 0 {
		return math.NaN()
	}
	if x == 1 {
		return math.Inf(-1)
	}
	return Ei(math.Log(x))
}

// Si returns the sine integral,
//
//	Si(x) = ∫₀ˣ sin(t)/t dt.
func Si(x float64) float64 {
	si, _ := SinCosIntegral(x)
	return si
}

// Ci returns the cosine integral,
//
//	Ci(x) = γ + ln x + ∫₀ˣ (cos(t)-1)/t dt,   x > 0.
func Ci(x float64) float64 {
	_, ci := SinCosIntegral(x)
	return ci
}

// SinCosIntegral returns the sine integral Si(x) and cosine integral Ci(x)
// simultaneously. For x < 0, Si is odd (Si(-x) = -Si(x)) and Ci(x) is
// evaluated as Ci(|x|).
func SinCosIntegral(x float64) (si, ci float64) {
	neg := false
	if x < 0 {
		x = -x
		neg = true
	}
	if x == 0 {
		return 0, math.Inf(-1)
	}
	if x <= 2 {
		// Power series.
		// Si(x) = Σ_{k≥0} (-1)^k x^{2k+1}/((2k+1)(2k+1)!)
		// Ci(x) = γ + ln x + Σ_{k≥1} (-1)^k x^{2k}/((2k)(2k)!)
		x2 := x * x
		term := x
		si = x
		for k := 1; k <= 200; k++ {
			term *= -x2 / float64((2*k)*(2*k+1))
			add := term / float64(2*k+1)
			si += add
			if math.Abs(add) < 1e-18*math.Abs(si) {
				break
			}
		}
		termc := 1.0
		sumc := 0.0
		for k := 1; k <= 200; k++ {
			termc *= -x2 / float64((2*k-1)*(2*k))
			add := termc / float64(2*k)
			sumc += add
			if math.Abs(add) < 1e-18*(math.Abs(sumc)+1) {
				break
			}
		}
		ci = specialEulerGamma + math.Log(x) + sumc
	} else {
		// Evaluate the complex continued fraction for the auxiliary integral
		//   ∫_x^∞ e^{it}/t dt = -Ci(x) + i(Si(x) - π/2)
		// via a modified Lentz algorithm (after Numerical Recipes).
		const tiny = 1e-300
		b := complex(1, x)
		c := complex(1/tiny, 0)
		d := 1 / b
		h := d
		for i := 2; i <= 400; i++ {
			a := -complex(float64((i-1)*(i-1)), 0)
			b += complex(2, 0)
			d = a*d + b
			if cmplx.Abs(d) < tiny {
				d = complex(tiny, 0)
			}
			c = b + a/c
			if cmplx.Abs(c) < tiny {
				c = complex(tiny, 0)
			}
			d = 1 / d
			del := d * c
			h *= del
			if math.Abs(real(del)-1)+math.Abs(imag(del)) < 1e-16 {
				break
			}
		}
		sinx, cosx := math.Sincos(x)
		h = complex(cosx, -sinx) * h
		ci = -real(h)
		si = math.Pi/2 + imag(h)
	}
	if neg {
		si = -si
	}
	return si, ci
}

// Shi returns the hyperbolic sine integral,
//
//	Shi(x) = ∫₀ˣ sinh(t)/t dt.
func Shi(x float64) float64 {
	shi, _ := HypSinCosIntegral(x)
	return shi
}

// Chi returns the hyperbolic cosine integral,
//
//	Chi(x) = γ + ln x + ∫₀ˣ (cosh(t)-1)/t dt,   x > 0.
func Chi(x float64) float64 {
	_, chi := HypSinCosIntegral(x)
	return chi
}

// HypSinCosIntegral returns the hyperbolic sine integral Shi(x) and hyperbolic
// cosine integral Chi(x) simultaneously. Shi is odd; Chi is evaluated at |x|.
func HypSinCosIntegral(x float64) (shi, chi float64) {
	neg := false
	if x < 0 {
		x = -x
		neg = true
	}
	if x == 0 {
		return 0, math.Inf(-1)
	}
	x2 := x * x
	// Shi(x) = Σ_{k≥0} x^{2k+1}/((2k+1)(2k+1)!)
	term := x
	shi = x
	for k := 1; k <= 200; k++ {
		term *= x2 / float64((2*k)*(2*k+1))
		add := term / float64(2*k+1)
		shi += add
		if math.Abs(add) < 1e-18*math.Abs(shi) {
			break
		}
	}
	// Chi(x) = γ + ln x + Σ_{k≥1} x^{2k}/((2k)(2k)!)
	termc := 1.0
	sumc := 0.0
	for k := 1; k <= 200; k++ {
		termc *= x2 / float64((2*k-1)*(2*k))
		add := termc / float64(2*k)
		sumc += add
		if math.Abs(add) < 1e-18*(math.Abs(sumc)+1) {
			break
		}
	}
	chi = specialEulerGamma + math.Log(x) + sumc
	if neg {
		shi = -shi
	}
	return shi, chi
}

// -----------------------------------------------------------------------------
// Complete and incomplete elliptic integrals
// -----------------------------------------------------------------------------

// EllipticK returns the complete elliptic integral of the first kind,
//
//	K(m) = ∫₀^{π/2} dθ / √(1 - m sin²θ),
//
// with parameter m = k². It is defined for m < 1.
func EllipticK(m float64) float64 {
	if m >= 1 {
		if m == 1 {
			return math.Inf(1)
		}
		return math.NaN()
	}
	return specialRF(0, 1-m, 1)
}

// EllipticE returns the complete elliptic integral of the second kind,
//
//	E(m) = ∫₀^{π/2} √(1 - m sin²θ) dθ,
//
// with parameter m = k². It is defined for m ≤ 1.
func EllipticE(m float64) float64 {
	if m > 1 {
		return math.NaN()
	}
	if m == 1 {
		return 1
	}
	y := 1 - m
	return specialRF(0, y, 1) - m/3*specialRD(0, y, 1)
}

// EllipticF returns the incomplete elliptic integral of the first kind,
//
//	F(φ, m) = ∫₀^φ dθ / √(1 - m sin²θ),
//
// with amplitude φ (radians) and parameter m = k².
func EllipticF(phi, m float64) float64 {
	s := math.Sin(phi)
	c := math.Cos(phi)
	s2 := s * s
	return s * specialRF(c*c, 1-m*s2, 1)
}

// EllipticEInc returns the incomplete elliptic integral of the second kind,
//
//	E(φ, m) = ∫₀^φ √(1 - m sin²θ) dθ,
//
// with amplitude φ (radians) and parameter m = k².
func EllipticEInc(phi, m float64) float64 {
	s := math.Sin(phi)
	c := math.Cos(phi)
	s2 := s * s
	cc := c * c
	q := 1 - m*s2
	return s*specialRF(cc, q, 1) - m*s2*s/3*specialRD(cc, q, 1)
}

// EllipticPi returns the complete elliptic integral of the third kind,
//
//	Π(n, m) = ∫₀^{π/2} dθ / ((1 - n sin²θ)√(1 - m sin²θ)),
//
// with characteristic n and parameter m = k².
func EllipticPi(n, m float64) float64 {
	y := 1 - m
	return specialRF(0, y, 1) + n/3*specialRJ(0, y, 1, 1-n)
}

// EllipticPiInc returns the incomplete elliptic integral of the third kind,
//
//	Π(n; φ, m) = ∫₀^φ dθ / ((1 - n sin²θ)√(1 - m sin²θ)).
func EllipticPiInc(n, phi, m float64) float64 {
	s := math.Sin(phi)
	c := math.Cos(phi)
	s2 := s * s
	cc := c * c
	q := 1 - m*s2
	return s*specialRF(cc, q, 1) + n*s2*s/3*specialRJ(cc, q, 1, 1-n*s2)
}

// specialRF computes Carlson's symmetric elliptic integral of the first kind,
//
//	R_F(x, y, z) = ½ ∫₀^∞ dt / √((t+x)(t+y)(t+z)).
func specialRF(x, y, z float64) float64 {
	const errtol = 1e-13
	for i := 0; i < 100; i++ {
		lambda := math.Sqrt(x)*math.Sqrt(y) + math.Sqrt(y)*math.Sqrt(z) + math.Sqrt(z)*math.Sqrt(x)
		x = 0.25 * (x + lambda)
		y = 0.25 * (y + lambda)
		z = 0.25 * (z + lambda)
		mu := (x + y + z) / 3
		if mu == 0 {
			break
		}
		dx := math.Abs(1 - x/mu)
		dy := math.Abs(1 - y/mu)
		dz := math.Abs(1 - z/mu)
		if dx < errtol && dy < errtol && dz < errtol {
			e2 := (1-x/mu)*(1-y/mu) + (1-y/mu)*(1-z/mu) + (1-z/mu)*(1-x/mu)
			e3 := (1 - x/mu) * (1 - y/mu) * (1 - z/mu)
			return (1 - e2/10 + e3/14 + e2*e2/24 - 3*e2*e3/44) / math.Sqrt(mu)
		}
	}
	mu := (x + y + z) / 3
	return 1 / math.Sqrt(mu)
}

// specialRD computes Carlson's symmetric elliptic integral of the second kind,
//
//	R_D(x, y, z) = (3/2) ∫₀^∞ dt / ((t+z)√((t+x)(t+y)(t+z))).
func specialRD(x, y, z float64) float64 {
	const errtol = 1e-13
	sum := 0.0
	fac := 1.0
	for i := 0; i < 200; i++ {
		sx := math.Sqrt(x)
		sy := math.Sqrt(y)
		sz := math.Sqrt(z)
		lambda := sx*sy + sy*sz + sz*sx
		sum += fac / (sz * (z + lambda))
		fac *= 0.25
		x = 0.25 * (x + lambda)
		y = 0.25 * (y + lambda)
		z = 0.25 * (z + lambda)
		mu := (x + y + 3*z) / 5
		if mu == 0 {
			break
		}
		dx := math.Abs(1 - x/mu)
		dy := math.Abs(1 - y/mu)
		dz := math.Abs(1 - z/mu)
		if dx < errtol && dy < errtol && dz < errtol {
			ex := (1 - x/mu)
			ey := (1 - y/mu)
			ez := (1 - z/mu)
			ea := ex * ey
			eb := ez * ez
			ec := ea - eb
			ed := ea - 6*eb
			ee := ed + ec + ec
			s := ed*(-3.0/14+9.0/88*ed-4.5/26*ez*ee) +
				ez*(1.0/6*ee+ez*(-9.0/22*ec+3.0/26*ez*ea))
			return 3*sum + fac*(1+s)/(mu*math.Sqrt(mu))
		}
	}
	mu := (x + y + 3*z) / 5
	return 3*sum + fac/(mu*math.Sqrt(mu))
}

// specialRJ computes Carlson's symmetric elliptic integral of the third kind,
//
//	R_J(x, y, z, p) = (3/2) ∫₀^∞ dt / ((t+p)√((t+x)(t+y)(t+z))).
func specialRJ(x, y, z, p float64) float64 {
	const errtol = 1e-12
	sum := 0.0
	fac := 1.0
	for i := 0; i < 200; i++ {
		sx := math.Sqrt(x)
		sy := math.Sqrt(y)
		sz := math.Sqrt(z)
		lambda := sx*sy + sy*sz + sz*sx
		alpha := p*(sx+sy+sz) + sx*sy*sz
		alpha *= alpha
		beta := p * (p + lambda) * (p + lambda)
		sum += fac * specialRC(alpha, beta)
		fac *= 0.25
		x = 0.25 * (x + lambda)
		y = 0.25 * (y + lambda)
		z = 0.25 * (z + lambda)
		p = 0.25 * (p + lambda)
		mu := (x + y + z + 2*p) / 5
		if mu == 0 {
			break
		}
		dx := math.Abs(1 - x/mu)
		dy := math.Abs(1 - y/mu)
		dz := math.Abs(1 - z/mu)
		dp := math.Abs(1 - p/mu)
		if dx < errtol && dy < errtol && dz < errtol && dp < errtol {
			const (
				c1 = 3.0 / 14.0
				c2 = 1.0 / 3.0
				c3 = 3.0 / 22.0
				c4 = 3.0 / 26.0
			)
			c5 := 0.75 * c3
			c6 := 1.5 * c4
			c7 := 0.5 * c2
			c8 := c3 + c3
			ex := (1 - x/mu)
			ey := (1 - y/mu)
			ez := (1 - z/mu)
			ep := (1 - p/mu)
			ea := ex*(ey+ez) + ey*ez
			eb := ex * ey * ez
			ec := ep * ep
			ed := ea - 3*ec
			ee := eb + 2*ep*(ea-ec)
			s := 1 + ed*(-c1+c5*ed-c6*ep*ee) +
				eb*(c7+ep*(-c8+ep*c4)) +
				ep*ea*(c2-ep*c3) - c2*ep*ec
			return 3*sum + fac*s/(mu*math.Sqrt(mu))
		}
	}
	mu := (x + y + z + 2*p) / 5
	return 3*sum + fac/(mu*math.Sqrt(mu))
}

// specialRC computes Carlson's degenerate elliptic integral,
//
//	R_C(x, y) = ½ ∫₀^∞ dt / ((t+y)√(t+x)).
func specialRC(x, y float64) float64 {
	const errtol = 1e-13
	for i := 0; i < 100; i++ {
		lambda := 2*math.Sqrt(x)*math.Sqrt(y) + y
		x = 0.25 * (x + lambda)
		y = 0.25 * (y + lambda)
		mu := (x + 2*y) / 3
		if mu == 0 {
			break
		}
		s := (y - mu) / mu
		if math.Abs(s) < errtol {
			return (1 + s*s*(3.0/10+s*(1.0/7+s*(3.0/8+s*9.0/22)))) / math.Sqrt(mu)
		}
	}
	mu := (x + 2*y) / 3
	return 1 / math.Sqrt(mu)
}

// -----------------------------------------------------------------------------
// Zeta family
// -----------------------------------------------------------------------------

// Zeta returns the Riemann zeta function ζ(s) for real s ≠ 1.
//
//	ζ(s) = Σ_{n≥1} n^{-s}   (for s > 1, extended by analytic continuation).
//
// The reflection formula is used for s < 0.5.
func Zeta(s float64) float64 {
	if s == 1 {
		return math.Inf(1)
	}
	if s == 0 {
		return -0.5
	}
	if s < 0.5 {
		// Reflection: ζ(s) = 2^s π^{s-1} sin(πs/2) Γ(1-s) ζ(1-s).
		if s < 0 && s == math.Trunc(s) && int64(s)%2 == 0 {
			return 0 // trivial zeros at negative even integers
		}
		return math.Pow(2, s) * math.Pow(math.Pi, s-1) *
			math.Sin(math.Pi*s/2) * math.Gamma(1-s) * Zeta(1-s)
	}
	// ζ(s) = ζ(s, 1) via Euler–Maclaurin.
	return HurwitzZeta(s, 1)
}

// specialBernoulli2 lists the even-index Bernoulli numbers B₂, B₄, … used in
// the Euler–Maclaurin tail.
var specialBernoulli2 = []float64{
	1.0 / 6, -1.0 / 30, 1.0 / 42, -1.0 / 30, 5.0 / 66,
	-691.0 / 2730, 7.0 / 6, -3617.0 / 510,
}

// specialFactorial returns n! as a float64 for small n.
func specialFactorial(n int) float64 {
	r := 1.0
	for i := 2; i <= n; i++ {
		r *= float64(i)
	}
	return r
}

// HurwitzZeta returns the Hurwitz zeta function,
//
//	ζ(s, a) = Σ_{n≥0} (n + a)^{-s},   s > 1, a > 0,
//
// extended to s < 1 (s ≠ 1) by Euler–Maclaurin summation.
func HurwitzZeta(s, a float64) float64 {
	if s == 1 {
		return math.Inf(1)
	}
	const n = 24
	sum := 0.0
	for k := 0; k < n; k++ {
		sum += math.Pow(float64(k)+a, -s)
	}
	na := float64(n) + a
	// Boundary and integral terms.
	sum += math.Pow(na, -s) * 0.5
	sum += math.Pow(na, 1-s) / (s - 1)
	// Bernoulli tail: Σ_j B_{2j}/(2j)! · (s)_{2j-1} · na^{-s-2j+1}.
	poch := s                   // (s)_1
	napow := math.Pow(na, -s-1) // na^{-s-1}
	inv2 := 1.0 / (na * na)
	for j := 1; j <= len(specialBernoulli2); j++ {
		add := specialBernoulli2[j-1] / specialFactorial(2*j) * poch * napow
		sum += add
		if math.Abs(add) < 1e-19*math.Abs(sum) {
			break
		}
		poch *= (s + float64(2*j-1)) * (s + float64(2*j))
		napow *= inv2
	}
	return sum
}

// DirichletEta returns the Dirichlet eta function (alternating zeta),
//
//	η(s) = Σ_{n≥1} (-1)^{n-1} n^{-s} = (1 - 2^{1-s}) ζ(s).
func DirichletEta(s float64) float64 {
	if s == 1 {
		return math.Ln2
	}
	return (1 - math.Pow(2, 1-s)) * Zeta(s)
}

// DirichletBeta returns the Dirichlet beta function,
//
//	β(s) = Σ_{n≥0} (-1)^n (2n+1)^{-s}.
//
// It satisfies β(1) = π/4 and β(2) = Catalan's constant.
func DirichletBeta(s float64) float64 {
	if s == 1 {
		return math.Pi / 4
	}
	// β(s) = 4^{-s} (ζ(s,1/4) - ζ(s,3/4)).
	return math.Pow(4, -s) * (HurwitzZeta(s, 0.25) - HurwitzZeta(s, 0.75))
}

// -----------------------------------------------------------------------------
// Polylogarithm
// -----------------------------------------------------------------------------

// Li2 returns the dilogarithm,
//
//	Li₂(x) = Σ_{k≥1} x^k / k² = -∫₀ˣ ln(1-t)/t dt,
//
// for real x ≤ 1. Analytic continuation handles x < -1 and 0 < x ≤ 1.
func Li2(x float64) float64 {
	if x == 1 {
		return math.Pi * math.Pi / 6
	}
	if x == -1 {
		return -math.Pi * math.Pi / 12
	}
	if x > 1 {
		return math.NaN()
	}
	// Use identities to map x into [-1, 0.5] for fast convergence.
	if x < -1 {
		l := math.Log(-x)
		return -math.Pi*math.Pi/6 - 0.5*l*l - Li2(1/x)
	}
	if x < 0 {
		return 0.5*Li2(x*x) - Li2(-x)
	}
	if x <= 0.5 {
		return specialLi2Series(x)
	}
	// 0.5 < x < 1
	l := math.Log(x)
	return math.Pi*math.Pi/6 - l*math.Log(1-x) - specialLi2Series(1-x)
}

// specialLi2Series sums the defining power series of the dilogarithm; it is
// only used for arguments in [0, 0.5] where convergence is rapid.
func specialLi2Series(x float64) float64 {
	sum := 0.0
	term := 1.0
	for k := 1; k <= 1000; k++ {
		term *= x
		add := term / float64(k*k)
		sum += add
		if math.Abs(add) < 1e-18*math.Abs(sum)+1e-300 {
			break
		}
	}
	return sum
}

// Li3 returns the trilogarithm,
//
//	Li₃(x) = Σ_{k≥1} x^k / k³,
//
// for real x with |x| ≤ 1 (and x < 1). Values outside are returned via the
// series where it converges.
func Li3(x float64) float64 {
	if x == 1 {
		return specialApery
	}
	if math.Abs(x) > 1 {
		return math.NaN()
	}
	sum := 0.0
	term := 1.0
	for k := 1; k <= 100000; k++ {
		term *= x
		add := term / float64(k) / float64(k) / float64(k)
		sum += add
		if math.Abs(add) < 1e-18*math.Abs(sum)+1e-300 {
			break
		}
	}
	return sum
}

// specialApery is Apéry's constant ζ(3) = Li₃(1).
const specialApery = 1.2020569031595942853997381615114499907650

// Polylog returns the polylogarithm Liₛ(x) = Σ_{k≥1} x^k / k^s for integer
// order s and real |x| < 1 (and x = 1 when s > 1, giving ζ(s)). It sums the
// defining series directly.
func Polylog(s int, x float64) float64 {
	if x == 1 {
		if s > 1 {
			return Zeta(float64(s))
		}
		return math.Inf(1)
	}
	if math.Abs(x) > 1 {
		return math.NaN()
	}
	sum := 0.0
	term := 1.0
	for k := 1; k <= 200000; k++ {
		term *= x
		add := term / math.Pow(float64(k), float64(s))
		sum += add
		if math.Abs(add) < 1e-18*math.Abs(sum)+1e-300 {
			break
		}
	}
	return sum
}

// -----------------------------------------------------------------------------
// Lambert W
// -----------------------------------------------------------------------------

// LambertW returns the principal branch W₀(x) of the Lambert W function, the
// solution w of w·e^w = x for x ≥ -1/e. It returns NaN for x < -1/e.
func LambertW(x float64) float64 {
	const invE = 0.36787944117144232159552377016146086744581 // 1/e
	if x < -invE {
		return math.NaN()
	}
	if x == 0 {
		return 0
	}
	if x == -invE {
		return -1
	}
	// Initial guess.
	var w float64
	if x < 0 {
		// Series near the branch point at x = -1/e.
		p := math.Sqrt(2 * (math.E*x + 1))
		w = -1 + p - p*p/3 + 11.0/72*p*p*p
	} else {
		// log1p(x) is finite for all x ≥ 0 and a good starting point.
		w = math.Log1p(x)
	}
	// Halley iteration.
	for i := 0; i < 60; i++ {
		ew := math.Exp(w)
		f := w*ew - x
		wp1 := w + 1
		w -= f / (ew*wp1 - (w+2)*f/(2*wp1))
		if math.Abs(f) < 1e-16*(1+math.Abs(x)) {
			break
		}
	}
	return w
}

// LambertWm1 returns the secondary real branch W₋₁(x) of the Lambert W
// function for -1/e ≤ x < 0, the solution w ≤ -1 of w·e^w = x.
func LambertWm1(x float64) float64 {
	const invE = 0.36787944117144232159552377016146086744581
	if x < -invE || x >= 0 {
		if x == -invE {
			return -1
		}
		return math.NaN()
	}
	// Initial guess based on asymptotics of the -1 branch.
	l1 := math.Log(-x)
	l2 := math.Log(-l1)
	w := l1 - l2 + l2/l1
	for i := 0; i < 80; i++ {
		ew := math.Exp(w)
		f := w*ew - x
		wp1 := w + 1
		w -= f / (ew*wp1 - (w+2)*f/(2*wp1))
		if math.Abs(f) < 1e-16*(1+math.Abs(x)) {
			break
		}
	}
	return w
}

// -----------------------------------------------------------------------------
// Gamma family: digamma, polygamma, incomplete gamma and beta
// -----------------------------------------------------------------------------

// Digamma returns the digamma function ψ(x) = Γ'(x)/Γ(x), the logarithmic
// derivative of the gamma function, for real x that is not a non-positive
// integer.
func Digamma(x float64) float64 {
	if x <= 0 && x == math.Trunc(x) {
		return math.NaN()
	}
	result := 0.0
	// Reflection for negative x.
	if x < 0 {
		return Digamma(1-x) - math.Pi/math.Tan(math.Pi*x)
	}
	// Recurrence to push argument above 6.
	for x < 6 {
		result -= 1 / x
		x++
	}
	// Asymptotic expansion.
	inv := 1 / x
	inv2 := inv * inv
	result += math.Log(x) - 0.5*inv -
		inv2*(1.0/12-inv2*(1.0/120-inv2*(1.0/252-inv2*(1.0/240-inv2*(1.0/132)))))
	return result
}

// Trigamma returns the trigamma function ψ⁽¹⁾(x), the first derivative of the
// digamma function, for real x that is not a non-positive integer.
func Trigamma(x float64) float64 {
	if x <= 0 && x == math.Trunc(x) {
		return math.NaN()
	}
	if x < 0 {
		// Reflection: ψ'(1-x) + ψ'(x) = π²/sin²(πx).
		s := math.Sin(math.Pi * x)
		return math.Pi*math.Pi/(s*s) - Trigamma(1-x)
	}
	result := 0.0
	for x < 12 {
		result += 1 / (x * x)
		x++
	}
	inv := 1 / x
	inv2 := inv * inv
	// ψ'(x) ~ 1/x + 1/(2x²) + Σ_{k≥1} B_{2k} x^{-2k-1}
	//       = 1/x + 1/(2x²) + 1/6 x^{-3} - 1/30 x^{-5} + 1/42 x^{-7} - …
	result += inv + 0.5*inv2 +
		inv2*inv*(1.0/6-inv2*(1.0/30-inv2*(1.0/42-inv2*(1.0/30-inv2*(5.0/66)))))
	return result
}

// Polygamma returns the polygamma function ψ⁽ⁿ⁾(x), the n-th derivative of the
// digamma function, for integer order n ≥ 0 and real x > 0.
//
//	ψ⁽ⁿ⁾(x) = (-1)^{n+1} n! Σ_{k≥0} (x+k)^{-(n+1)}   for n ≥ 1.
func Polygamma(n int, x float64) float64 {
	if n < 0 {
		return math.NaN()
	}
	if n == 0 {
		return Digamma(x)
	}
	if x <= 0 && x == math.Trunc(x) {
		return math.NaN()
	}
	// ψ⁽ⁿ⁾(x) = (-1)^{n+1} n! ζ(n+1, x).
	sign := 1.0
	if n%2 == 0 {
		sign = -1.0
	}
	return sign * specialFactorial(n) * HurwitzZeta(float64(n+1), x)
}

// Beta returns the beta function B(a, b) = Γ(a)Γ(b)/Γ(a+b) for a, b > 0.
func Beta(a, b float64) float64 {
	lg, sign := math.Lgamma(a)
	lb, sb := math.Lgamma(b)
	lab, sab := math.Lgamma(a + b)
	return float64(sign*sb*sab) * math.Exp(lg+lb-lab)
}

// LogGamma returns the natural logarithm of the absolute value of the gamma
// function, ln|Γ(x)|.
func LogGamma(x float64) float64 {
	v, _ := math.Lgamma(x)
	return v
}

// GammaP returns the regularized lower incomplete gamma function,
//
//	P(a, x) = (1/Γ(a)) ∫₀ˣ t^{a-1} e^{-t} dt,   a > 0, x ≥ 0.
func GammaP(a, x float64) float64 {
	if x < 0 || a <= 0 {
		return math.NaN()
	}
	if x == 0 {
		return 0
	}
	if x < a+1 {
		return specialGammaSeries(a, x)
	}
	return 1 - specialGammaCF(a, x)
}

// GammaQ returns the regularized upper incomplete gamma function,
//
//	Q(a, x) = 1 - P(a, x) = (1/Γ(a)) ∫ₓ^∞ t^{a-1} e^{-t} dt.
func GammaQ(a, x float64) float64 {
	if x < 0 || a <= 0 {
		return math.NaN()
	}
	if x == 0 {
		return 1
	}
	if x < a+1 {
		return 1 - specialGammaSeries(a, x)
	}
	return specialGammaCF(a, x)
}

// specialGammaSeries evaluates P(a, x) by its power series, valid for x < a+1.
func specialGammaSeries(a, x float64) float64 {
	ap := a
	sum := 1 / a
	del := sum
	for i := 0; i < 500; i++ {
		ap++
		del *= x / ap
		sum += del
		if math.Abs(del) < math.Abs(sum)*1e-16 {
			break
		}
	}
	lg, _ := math.Lgamma(a)
	return sum * math.Exp(-x+a*math.Log(x)-lg)
}

// specialGammaCF evaluates Q(a, x) by a continued fraction, valid for x ≥ a+1.
func specialGammaCF(a, x float64) float64 {
	const tiny = 1e-300
	b := x + 1 - a
	c := 1.0 / tiny
	d := 1.0 / b
	h := d
	for i := 1; i <= 500; i++ {
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
		if math.Abs(del-1) < 1e-16 {
			break
		}
	}
	lg, _ := math.Lgamma(a)
	return math.Exp(-x+a*math.Log(x)-lg) * h
}

// LowerIncompleteGamma returns the (non-regularized) lower incomplete gamma
// function γ(a, x) = ∫₀ˣ t^{a-1} e^{-t} dt = P(a, x)·Γ(a).
func LowerIncompleteGamma(a, x float64) float64 {
	return GammaP(a, x) * math.Gamma(a)
}

// UpperIncompleteGamma returns the (non-regularized) upper incomplete gamma
// function Γ(a, x) = ∫ₓ^∞ t^{a-1} e^{-t} dt = Q(a, x)·Γ(a).
func UpperIncompleteGamma(a, x float64) float64 {
	return GammaQ(a, x) * math.Gamma(a)
}

// BetaInc returns the regularized incomplete beta function,
//
//	I_x(a, b) = (1/B(a,b)) ∫₀ˣ t^{a-1}(1-t)^{b-1} dt,
//
// for 0 ≤ x ≤ 1 and a, b > 0.
func BetaInc(a, b, x float64) float64 {
	if x < 0 || x > 1 {
		return math.NaN()
	}
	if x == 0 {
		return 0
	}
	if x == 1 {
		return 1
	}
	lbeta := specialLogBeta(a, b)
	front := math.Exp(math.Log(x)*a + math.Log(1-x)*b - lbeta)
	if x < (a+1)/(a+b+2) {
		return front * specialBetaCF(a, b, x) / a
	}
	return 1 - front*specialBetaCF(b, a, 1-x)/b
}

// specialLogBeta returns ln B(a, b).
func specialLogBeta(a, b float64) float64 {
	la, _ := math.Lgamma(a)
	lb, _ := math.Lgamma(b)
	lab, _ := math.Lgamma(a + b)
	return la + lb - lab
}

// specialBetaCF evaluates the continued fraction used by the regularized
// incomplete beta function (Lentz's algorithm).
func specialBetaCF(a, b, x float64) float64 {
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
	for m := 1; m <= 500; m++ {
		mf := float64(m)
		m2 := 2 * mf
		aa := mf * (b - mf) * x / ((qam + m2) * (a + m2))
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
		aa = -(a + mf) * (qab + mf) * x / ((a + m2) * (qap + m2))
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
		if math.Abs(del-1) < 1e-16 {
			break
		}
	}
	return h
}

// IncompleteBeta returns the (non-regularized) incomplete beta function
// B(x; a, b) = ∫₀ˣ t^{a-1}(1-t)^{b-1} dt = I_x(a, b)·B(a, b).
func IncompleteBeta(a, b, x float64) float64 {
	return BetaInc(a, b, x) * Beta(a, b)
}

// -----------------------------------------------------------------------------
// Bessel functions
// -----------------------------------------------------------------------------

// BesselJ returns the Bessel function of the first kind of integer order n,
// Jₙ(x).
func BesselJ(n int, x float64) float64 { return math.Jn(n, x) }

// BesselY returns the Bessel function of the second kind of integer order n,
// Yₙ(x), for x > 0.
func BesselY(n int, x float64) float64 { return math.Yn(n, x) }

// BesselI returns the modified Bessel function of the first kind of integer
// order n, Iₙ(x).
func BesselI(n int, x float64) float64 {
	if n < 0 {
		n = -n
	}
	if n == 0 {
		return BesselI0(x)
	}
	if n == 1 {
		return BesselI1(x)
	}
	if x == 0 {
		return 0
	}
	const acc = 40.0
	const bigno = 1e10
	const bigni = 1e-10
	tox := 2 / math.Abs(x)
	bip := 0.0
	ans := 0.0
	bi := 1.0
	m := 2 * (n + int(math.Sqrt(acc*float64(n))))
	for j := m; j > 0; j-- {
		bim := bip + float64(j)*tox*bi
		bip = bi
		bi = bim
		if math.Abs(bi) > bigno {
			ans *= bigni
			bi *= bigni
			bip *= bigni
		}
		if j == n {
			ans = bip
		}
	}
	ans *= BesselI0(x) / bi
	if x < 0 && n%2 == 1 {
		return -ans
	}
	return ans
}

// BesselK returns the modified Bessel function of the second kind of integer
// order n, Kₙ(x), for x > 0.
func BesselK(n int, x float64) float64 {
	if n < 0 {
		n = -n
	}
	if n == 0 {
		return BesselK0(x)
	}
	if n == 1 {
		return BesselK1(x)
	}
	tox := 2 / x
	bkm := BesselK0(x)
	bk := BesselK1(x)
	for j := 1; j < n; j++ {
		bkp := bkm + float64(j)*tox*bk
		bkm = bk
		bk = bkp
	}
	return bk
}

// -----------------------------------------------------------------------------
// Airy functions
// -----------------------------------------------------------------------------

// Airy returns the Airy functions Ai(x), Ai'(x), Bi(x) and Bi'(x)
// simultaneously. The implementation expresses the Airy functions in terms of
// modified Bessel functions (for x > 0) and ordinary Bessel functions (for
// x < 0), with a Maclaurin evaluation near the origin.
func Airy(x float64) (ai, aip, bi, bip float64) {
	// Values and derivatives at 0.
	const ai0 = 0.3550280538878172392600631860041831763979
	const aip0 = -0.2588194037928067984051835601892039634793
	const bi0 = 0.6149266274460007351509223690936135535960
	const bip0 = 0.4482883573538263579148237103988283908054

	absx := math.Abs(x)
	if absx < 1e-6 {
		// First-order Taylor about 0 (y'' = x·y ⇒ y''(0)=0, y'''(0)=y(0)).
		x2 := x * x
		ai = ai0 + aip0*x + ai0*x2*x/6
		aip = aip0 + ai0*x2/2
		bi = bi0 + bip0*x + bi0*x2*x/6
		bip = bip0 + bi0*x2/2
		return ai, aip, bi, bip
	}
	sqrt3 := math.Sqrt(3.0)
	z := 2.0 / 3.0 * absx * math.Sqrt(absx)
	rootx := math.Sqrt(absx)
	if x > 0 {
		k13 := specialBesselKfrac(1.0/3.0, z)
		k23 := specialBesselKfrac(2.0/3.0, z)
		i13 := specialBesselIfrac(1.0/3.0, z)
		im13 := specialBesselIfrac(-1.0/3.0, z)
		i23 := specialBesselIfrac(2.0/3.0, z)
		im23 := specialBesselIfrac(-2.0/3.0, z)
		ai = (1 / math.Pi) * math.Sqrt(absx/3) * k13
		aip = -(absx / (math.Pi * sqrt3)) * k23
		bi = math.Sqrt(absx/3) * (im13 + i13)
		bip = (absx / sqrt3) * (im23 + i23)
		return ai, aip, bi, bip
	}
	// x < 0: express via ordinary Bessel J of fractional order.
	j13 := specialBesselJfrac(1.0/3.0, z)
	jm13 := specialBesselJfrac(-1.0/3.0, z)
	j23 := specialBesselJfrac(2.0/3.0, z)
	jm23 := specialBesselJfrac(-2.0/3.0, z)
	ai = (rootx / 3) * (jm13 + j13)
	bi = (rootx / sqrt3) * (jm13 - j13)
	aip = (absx / 3) * (jm23 - j23)
	bip = (absx / sqrt3) * (jm23 + j23)
	return ai, aip, bi, bip
}

// specialBesselJfrac computes J_ν(x) for fractional order ν by ascending power
// series, adequate for the moderate arguments arising in the Airy evaluation.
func specialBesselJfrac(nu, x float64) float64 {
	// J_ν(x) = Σ_{m≥0} (-1)^m /(m! Γ(m+ν+1)) (x/2)^{2m+ν}
	halfx := x / 2
	sum := 0.0
	logHalf := math.Log(halfx)
	for m := 0; m < 200; m++ {
		lg1, _ := math.Lgamma(float64(m) + 1)
		lg2, _ := math.Lgamma(float64(m) + nu + 1)
		lterm := (2*float64(m)+nu)*logHalf - lg1 - lg2
		term := math.Exp(lterm)
		if m%2 == 1 {
			term = -term
		}
		sum += term
		if math.Abs(term) < 1e-18*math.Abs(sum)+1e-300 {
			break
		}
	}
	return sum
}

// specialBesselIfrac computes I_ν(x) for fractional order ν by ascending power
// series.
func specialBesselIfrac(nu, x float64) float64 {
	halfx := x / 2
	sum := 0.0
	logHalf := math.Log(halfx)
	for m := 0; m < 300; m++ {
		lg1, _ := math.Lgamma(float64(m) + 1)
		lg2, _ := math.Lgamma(float64(m) + nu + 1)
		lterm := (2*float64(m)+nu)*logHalf - lg1 - lg2
		sum += math.Exp(lterm)
		if math.Exp(lterm) < 1e-18*math.Abs(sum)+1e-300 {
			break
		}
	}
	return sum
}

// specialBesselKfrac computes K_ν(x) for fractional order ν using
// K_ν(x) = π/2 · (I_{-ν}(x) - I_ν(x)) / sin(νπ).
func specialBesselKfrac(nu, x float64) float64 {
	return math.Pi / 2 * (specialBesselIfrac(-nu, x) - specialBesselIfrac(nu, x)) / math.Sin(nu*math.Pi)
}

// -----------------------------------------------------------------------------
// Miscellaneous
// -----------------------------------------------------------------------------

// Sinc returns the unnormalized cardinal sine, sinc(x) = sin(x)/x, with the
// removable singularity Sinc(0) = 1 handled exactly.
func Sinc(x float64) float64 {
	if x == 0 {
		return 1
	}
	return math.Sin(x) / x
}

// SincNorm returns the normalized cardinal sine, sinc(x) = sin(πx)/(πx), with
// SincNorm(0) = 1.
func SincNorm(x float64) float64 {
	if x == 0 {
		return 1
	}
	px := math.Pi * x
	return math.Sin(px) / px
}

// Struve returns the Struve function Hₙ(x) of integer order n by its power
// series,
//
//	Hₙ(x) = Σ_{m≥0} (-1)^m /(Γ(m+3/2) Γ(m+n+3/2)) (x/2)^{2m+n+1}.
func Struve(n int, x float64) float64 {
	halfx := x / 2
	if halfx == 0 {
		if n == -1 {
			return 2 / math.Pi
		}
		return 0
	}
	logHalf := math.Log(math.Abs(halfx))
	sum := 0.0
	nf := float64(n)
	for m := 0; m < 300; m++ {
		lg1, _ := math.Lgamma(float64(m) + 1.5)
		lg2, _ := math.Lgamma(float64(m) + nf + 1.5)
		lterm := (2*float64(m)+nf+1)*logHalf - lg1 - lg2
		term := math.Exp(lterm)
		if m%2 == 1 {
			term = -term
		}
		sum += term
		if math.Abs(term) < 1e-18*math.Abs(sum)+1e-300 {
			break
		}
	}
	if halfx < 0 && (n+1)%2 != 0 {
		// account for sign of (x/2)^{n+1} when x<0
		sum = -sum
	}
	return sum
}

// Gamma returns the gamma function Γ(x); it delegates to the standard library
// and is provided for API completeness alongside the incomplete forms.
func Gamma(x float64) float64 { return math.Gamma(x) }

// ReciprocalGamma returns 1/Γ(x), which is entire and vanishes at the
// non-positive integers.
func ReciprocalGamma(x float64) float64 {
	if x <= 0 && x == math.Trunc(x) {
		return 0
	}
	lg, sign := math.Lgamma(x)
	return float64(sign) * math.Exp(-lg)
}
