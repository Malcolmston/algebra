package analyticnt

import (
	"math"
	"math/cmplx"
)

// bernoulliEven holds B_{2k} for k = 1..8 used by the Euler–Maclaurin zeta
// evaluation: B2, B4, …, B16.
var bernoulliEven = [...]float64{
	1.0 / 6.0,
	-1.0 / 30.0,
	1.0 / 42.0,
	-1.0 / 30.0,
	5.0 / 66.0,
	-691.0 / 2730.0,
	7.0 / 6.0,
	-3617.0 / 510.0,
}

// Zeta returns the Riemann zeta function ζ(s) for a real argument s ≠ 1, using
// the Euler–Maclaurin summation formula, which analytically continues ζ to the
// whole real line except the pole at s = 1.
func Zeta(s float64) float64 {
	if s == 1 {
		return math.Inf(1)
	}
	// Even negative integers are trivial zeros; return exact 0 to avoid noise.
	if s < 0 && s == math.Trunc(s) && math.Mod(s, 2) == 0 {
		return 0
	}
	const N = 20
	var sum float64
	for n := 1; n < N; n++ {
		sum += math.Pow(float64(n), -s)
	}
	nf := float64(N)
	sum += math.Pow(nf, -s) / 2
	sum += math.Pow(nf, 1-s) / (s - 1)
	// Bernoulli correction terms.
	term := math.Pow(nf, -s-1) * s
	fact := 1.0
	for k := 1; k <= len(bernoulliEven); k++ {
		fact = bernoulliEven[k-1]
		sum += fact / factorial(2*k) * term
		// Advance the falling-factorial product (s)(s+1)…(s+2k-1) and n^{-s-2k}.
		term *= (s + float64(2*k-1)) * (s + float64(2*k)) / (nf * nf)
	}
	return sum
}

// factorial returns m! as a float64 for small m.
func factorial(m int) float64 {
	f := 1.0
	for i := 2; i <= m; i++ {
		f *= float64(i)
	}
	return f
}

// ZetaComplex returns ζ(s) for a complex argument s ≠ 1 via Euler–Maclaurin
// summation. It is accurate throughout the critical strip for moderate |Im s|.
func ZetaComplex(s complex128) complex128 {
	if s == complex(1, 0) {
		return cmplx.Inf()
	}
	const N = 24
	var sum complex128
	for n := 1; n < N; n++ {
		sum += cmplx.Pow(complex(float64(n), 0), -s)
	}
	nf := complex(float64(N), 0)
	sum += cmplx.Pow(nf, -s) / 2
	sum += cmplx.Pow(nf, 1-s) / (s - 1)
	term := cmplx.Pow(nf, -s-1) * s
	n2 := nf * nf
	for k := 1; k <= len(bernoulliEven); k++ {
		coef := complex(bernoulliEven[k-1]/factorial(2*k), 0)
		sum += coef * term
		term *= (s + complex(float64(2*k-1), 0)) * (s + complex(float64(2*k), 0)) / n2
	}
	return sum
}

// DirichletEta returns the Dirichlet eta function η(s) = Σ (-1)^{n-1} n^{-s} =
// (1 − 2^{1−s}) ζ(s) for real s. At s = 1 the removable singularity gives
// η(1) = ln 2.
func DirichletEta(s float64) float64 {
	if s == 1 {
		return math.Ln2
	}
	return (1 - math.Pow(2, 1-s)) * Zeta(s)
}

// DirichletEtaComplex returns η(s) = (1 − 2^{1−s}) ζ(s) for complex s, with the
// removable value η(1) = ln 2.
func DirichletEtaComplex(s complex128) complex128 {
	if s == complex(1, 0) {
		return complex(math.Ln2, 0)
	}
	return (1 - cmplx.Pow(2, 1-s)) * ZetaComplex(s)
}

// ZetaPrime returns the derivative ζ'(s) for real s ≠ 1 via a central finite
// difference of the Euler–Maclaurin evaluation.
func ZetaPrime(s float64) float64 {
	const h = 1e-5
	return (Zeta(s+h) - Zeta(s-h)) / (2 * h)
}

// LogGammaComplex returns log Γ(z) (principal branch) for complex z with
// positive real part sufficiently large, using the Lanczos approximation. It is
// used by the Riemann–Siegel theta function.
func LogGammaComplex(z complex128) complex128 {
	// Lanczos g=7, n=9 coefficients.
	g := 7.0
	c := []float64{
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
	if real(z) < 0.5 {
		// Reflection: Γ(z)Γ(1−z) = π/sin(πz).
		return cmplx.Log(complex(math.Pi, 0)/cmplx.Sin(math.Pi*z)) - LogGammaComplex(1-z)
	}
	z -= 1
	x := complex(c[0], 0)
	for i := 1; i < len(c); i++ {
		x += complex(c[i], 0) / (z + complex(float64(i), 0))
	}
	t := z + complex(g+0.5, 0)
	return complex(0.5*math.Log(2*math.Pi), 0) + (z+0.5)*cmplx.Log(t) - t + cmplx.Log(x)
}

// RiemannSiegelTheta returns the Riemann–Siegel theta function
// θ(t) = arg Γ(1/4 + it/2) − (t/2) ln π, real for real t. It is central to the
// study of ζ(1/2 + it).
func RiemannSiegelTheta(t float64) float64 {
	z := complex(0.25, t/2)
	lg := LogGammaComplex(z)
	return imag(lg) - t/2*math.Log(math.Pi)
}

// RiemannSiegelThetaAsymptotic returns the standard asymptotic expansion of
// θ(t) for large t; it is faster than RiemannSiegelTheta but only accurate for
// t larger than roughly 10.
func RiemannSiegelThetaAsymptotic(t float64) float64 {
	return t/2*math.Log(t/(2*math.Pi)) - t/2 - math.Pi/8 +
		1/(48*t) + 7/(5760*math.Pow(t, 3))
}

// HardyZ returns Hardy's Z-function Z(t) = e^{iθ(t)} ζ(1/2 + it), which is real
// for real t and shares its real zeros with ζ on the critical line. It is
// evaluated here directly from ζ(1/2 + it).
func HardyZ(t float64) float64 {
	theta := RiemannSiegelTheta(t)
	zt := ZetaComplex(complex(0.5, t))
	rot := cmplx.Exp(complex(0, theta)) * zt
	return real(rot)
}

// RiemannSiegelZ is a synonym for HardyZ, the real-valued function whose sign
// changes locate the nontrivial zeros of ζ on the critical line.
func RiemannSiegelZ(t float64) float64 { return HardyZ(t) }

// GramPoint returns the n-th Gram point g_n, the solution of θ(g_n) = nπ, found
// by Newton iteration. Gram points interlace the zeros of ζ for small n.
func GramPoint(n int) float64 {
	target := float64(n) * math.Pi
	// Initial guess from the asymptotic inverse of θ.
	t := 2 * math.Pi * math.Exp(1+lambertW(float64(n)/math.E))
	if t < 10 {
		t = 10
	}
	for i := 0; i < 60; i++ {
		f := RiemannSiegelTheta(t) - target
		const h = 1e-4
		fp := (RiemannSiegelTheta(t+h) - RiemannSiegelTheta(t-h)) / (2 * h)
		if fp == 0 {
			break
		}
		nt := t - f/fp
		if math.Abs(nt-t) < 1e-9 {
			t = nt
			break
		}
		t = nt
	}
	return t
}

// lambertW returns the principal branch W0(x) for x >= -1/e via Newton
// iteration. It is a small internal helper.
func lambertW(x float64) float64 {
	if x < -1/math.E {
		return math.NaN()
	}
	w := 0.0
	if x > math.E {
		w = math.Log(x) - math.Log(math.Log(x))
	} else if x > 0 {
		w = x / math.E
	} else {
		w = 0
	}
	for i := 0; i < 60; i++ {
		ew := math.Exp(w)
		f := w*ew - x
		fp := ew*(w+1) - (w+2)*f/(2*w+2)
		if fp == 0 {
			break
		}
		nw := w - f/fp
		if math.Abs(nw-w) < 1e-14 {
			w = nw
			break
		}
		w = nw
	}
	return w
}

// ZetaZero returns the imaginary part of the n-th nontrivial zero of ζ on the
// critical line (n >= 1), located by bracketing sign changes of HardyZ and
// bisecting. The first zero is ≈ 14.134725.
func ZetaZero(n int) float64 {
	if n < 1 {
		panic("analyticnt: ZetaZero requires n >= 1")
	}
	zeros := ZetaZeros(n)
	return zeros[n-1]
}

// ZetaZeros returns the imaginary parts of the first n nontrivial zeros of ζ on
// the critical line, in increasing order. n must be >= 1.
func ZetaZeros(n int) []float64 {
	if n < 1 {
		panic("analyticnt: ZetaZeros requires n >= 1")
	}
	out := make([]float64, 0, n)
	const step = 0.05
	t := 1.0
	prev := HardyZ(t)
	for len(out) < n {
		t2 := t + step
		cur := HardyZ(t2)
		if prev == 0 {
			out = append(out, t)
		} else if (prev < 0) != (cur < 0) {
			root := bisectHardy(t, t2)
			out = append(out, root)
		}
		t = t2
		prev = cur
		if t > 1e6 {
			break
		}
	}
	return out
}

// bisectHardy refines a bracketed root of HardyZ in [a,b] by bisection.
func bisectHardy(a, b float64) float64 {
	fa := HardyZ(a)
	for i := 0; i < 100; i++ {
		m := (a + b) / 2
		fm := HardyZ(m)
		if fm == 0 || (b-a) < 1e-11 {
			return m
		}
		if (fa < 0) != (fm < 0) {
			b = m
		} else {
			a = m
			fa = fm
		}
	}
	return (a + b) / 2
}

// RiemannVonMangoldtN returns the counting function N(T), the approximate number
// of nontrivial zeros of ζ with imaginary part in (0, T], given by the
// Riemann–von Mangoldt formula θ(T)/π + 1.
func RiemannVonMangoldtN(T float64) float64 {
	if T <= 0 {
		return 0
	}
	return RiemannSiegelTheta(T)/math.Pi + 1
}
