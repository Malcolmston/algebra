package complexanalysis

import (
	"math"
	"math/cmplx"
)

// complexanalysisLanczosG is the Lanczos parameter g used together with the
// coefficient table below (g = 7, n = 9).
const complexanalysisLanczosG = 7.0

// complexanalysisLanczosCoef are the Lanczos g=7 coefficients.
var complexanalysisLanczosCoef = [9]float64{
	0.99999999999980993,
	676.5203681218851,
	-1259.1392167224028,
	771.32342877765313,
	-176.61502916214059,
	12.507343278686905,
	-0.13857109526572012,
	9.9843695780195716e-6,
	1.5056327351493116e-7,
}

// Gamma returns the value of the Euler gamma function at z, computed with the
// Lanczos approximation and the reflection formula for Re(z) < 0.5. It is
// accurate to roughly 1e-13 relative error away from the poles at the
// non-positive integers.
func Gamma(z complex128) complex128 {
	if real(z) < 0.5 {
		// Reflection: Gamma(z) = pi / (sin(pi z) Gamma(1-z)).
		return math.Pi / (cmplx.Sin(math.Pi*z) * Gamma(1-z))
	}
	z -= 1
	x := complex(complexanalysisLanczosCoef[0], 0)
	for i := 1; i < len(complexanalysisLanczosCoef); i++ {
		x += complex(complexanalysisLanczosCoef[i], 0) / (z + complex(float64(i), 0))
	}
	t := z + complex(complexanalysisLanczosG+0.5, 0)
	return complex(math.Sqrt(2*math.Pi), 0) * cmplx.Pow(t, z+0.5) * cmplx.Exp(-t) * x
}

// LogGamma returns the principal branch of the logarithm of the gamma function
// at z. For Re(z) >= 0.5 it uses the Lanczos series directly, avoiding the
// overflow that Gamma would suffer for large arguments; elsewhere it falls back
// to Log(Gamma(z)).
func LogGamma(z complex128) complex128 {
	if real(z) < 0.5 {
		return cmplx.Log(Gamma(z))
	}
	z -= 1
	x := complex(complexanalysisLanczosCoef[0], 0)
	for i := 1; i < len(complexanalysisLanczosCoef); i++ {
		x += complex(complexanalysisLanczosCoef[i], 0) / (z + complex(float64(i), 0))
	}
	t := z + complex(complexanalysisLanczosG+0.5, 0)
	return complex(0.5*math.Log(2*math.Pi), 0) + (z+0.5)*cmplx.Log(t) - t + cmplx.Log(x)
}

// Digamma returns the logarithmic derivative of the gamma function,
// psi(z) = Gamma'(z)/Gamma(z), using recurrence to raise the real part and then
// the standard asymptotic expansion.
func Digamma(z complex128) complex128 {
	var result complex128
	for real(z) < 6 {
		result -= 1 / z
		z += 1
	}
	r := 1 / z
	result += cmplx.Log(z) - 0.5*r
	r2 := r * r
	// Asymptotic Bernoulli terms: 1/12, -1/120, 1/252, -1/240, 1/132.
	result -= r2 * (1.0/12 - r2*(1.0/120-r2*(1.0/252-r2*(1.0/240-r2*(1.0/132)))))
	return result
}

// Beta returns the Euler beta function B(a, b) = Gamma(a)Gamma(b)/Gamma(a+b),
// evaluated through LogGamma to remain stable for moderately large arguments.
func Beta(a, b complex128) complex128 {
	return cmplx.Exp(LogGamma(a) + LogGamma(b) - LogGamma(a+b))
}

// Factorial returns z! defined as Gamma(z+1); for a non-negative integer z it
// equals the ordinary factorial.
func Factorial(z complex128) complex128 { return Gamma(z + 1) }

// RisingFactorial returns the Pochhammer symbol (z)_n = z(z+1)...(z+n-1),
// computed as Gamma(z+n)/Gamma(z). It returns 1 for n == 0.
func RisingFactorial(z complex128, n int) complex128 {
	if n == 0 {
		return 1
	}
	if n > 0 && n <= 32 {
		// Exact product for small integer n avoids branch issues near poles.
		p := complex(1, 0)
		for k := 0; k < n; k++ {
			p *= z + complex(float64(k), 0)
		}
		return p
	}
	return cmplx.Exp(LogGamma(z+complex(float64(n), 0)) - LogGamma(z))
}

// Binomial returns the generalized binomial coefficient C(z, k) =
// Gamma(z+1)/(Gamma(k+1)Gamma(z-k+1)) for a non-negative integer k, evaluated
// as a falling-factorial product so it is exact for polynomial z.
func Binomial(z complex128, k int) complex128 {
	if k < 0 {
		return 0
	}
	num := complex(1, 0)
	for i := 0; i < k; i++ {
		num *= z - complex(float64(i), 0)
	}
	return num / complex(complexanalysisFactorialFloat(k), 0)
}

// Erf returns the error function of z, 2/sqrt(pi) times the integral of
// exp(-t^2) from 0 to z, summed from its Taylor series. The series converges
// for every z but loses precision for large |z| (roughly |z| > 6); within that
// range it is accurate to about 1e-12.
func Erf(z complex128) complex128 {
	// erf(z) = 2/sqrt(pi) * sum_{n>=0} (-1)^n z^(2n+1) / (n! (2n+1)).
	const maxIter = 200
	z2 := z * z
	term := z // n = 0 term before the 1/(2n+1) factor and factorial
	sum := z
	for n := 1; n < maxIter; n++ {
		term *= -z2 / complex(float64(n), 0)
		add := term / complex(float64(2*n+1), 0)
		sum += add
		if cmplx.Abs(add) <= 1e-18*cmplx.Abs(sum) {
			break
		}
	}
	return complex(2/math.Sqrt(math.Pi), 0) * sum
}

// Erfc returns the complementary error function 1 - Erf(z).
func Erfc(z complex128) complex128 { return 1 - Erf(z) }

// complexanalysisEtaBorwein evaluates the Dirichlet eta function eta(s) with
// Borwein's algorithm using n terms.
func complexanalysisEtaBorwein(s complex128, n int) complex128 {
	// d_k = n * sum_{i=0}^{k} t_i, with t_0 = 1/n and the recurrence below.
	d := make([]float64, n+1)
	t := 1.0 / float64(n)
	sum := t
	d[0] = float64(n) * sum
	for i := 1; i <= n; i++ {
		fi := float64(i)
		t = t * 4.0 * (float64(n) + fi - 1.0) * (float64(n) - fi + 1.0) / (2.0 * fi * (2.0*fi - 1.0))
		sum += t
		d[i] = float64(n) * sum
	}
	var acc complex128
	for k := 0; k <= n-1; k++ {
		sign := 1.0
		if k%2 == 1 {
			sign = -1.0
		}
		base := complex(float64(k+1), 0)
		acc += complex(sign*(d[k]-d[n]), 0) / cmplx.Pow(base, s)
	}
	return -acc / complex(d[n], 0)
}

// Zeta returns the Riemann zeta function at s via the Dirichlet eta function
// and analytic continuation, eta(s) = (1 - 2^(1-s)) zeta(s). It is valid across
// the complex plane except at the pole s = 1, and is accurate to about 1e-12
// for moderate |Im(s)|.
func Zeta(s complex128) complex128 {
	if s == 1 {
		return cmplx.Inf()
	}
	eta := complexanalysisEtaBorwein(s, 40)
	return eta / (1 - cmplx.Pow(2, 1-s))
}
