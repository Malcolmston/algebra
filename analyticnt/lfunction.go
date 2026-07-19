package analyticnt

import (
	"math"
	"math/cmplx"
)

// HurwitzZeta returns the Hurwitz zeta function ζ(s, a) = Σ_{n≥0} (n+a)^{-s}
// for real s > 1 and a > 0, via Euler–Maclaurin summation (which also
// continues it to s < 1, s ≠ 1).
func HurwitzZeta(s, a float64) float64 {
	if a <= 0 {
		panic("analyticnt: HurwitzZeta requires a > 0")
	}
	const N = 20
	var sum float64
	for n := 0; n < N; n++ {
		sum += math.Pow(float64(n)+a, -s)
	}
	na := float64(N) + a
	sum += math.Pow(na, 1-s) / (s - 1)
	sum += math.Pow(na, -s) / 2
	term := s * math.Pow(na, -s-1)
	na2 := na * na
	for k := 1; k <= len(bernoulliEven); k++ {
		sum += bernoulliEven[k-1] / factorial(2*k) * term
		term *= (s + float64(2*k-1)) * (s + float64(2*k)) / na2
	}
	return sum
}

// HurwitzZetaComplex returns ζ(s, a) for complex s (s ≠ 1) and real a > 0.
func HurwitzZetaComplex(s complex128, a float64) complex128 {
	if a <= 0 {
		panic("analyticnt: HurwitzZetaComplex requires a > 0")
	}
	const N = 24
	var sum complex128
	for n := 0; n < N; n++ {
		sum += cmplx.Pow(complex(float64(n)+a, 0), -s)
	}
	na := complex(float64(N)+a, 0)
	sum += cmplx.Pow(na, 1-s) / (s - 1)
	sum += cmplx.Pow(na, -s) / 2
	term := s * cmplx.Pow(na, -s-1)
	na2 := na * na
	for k := 1; k <= len(bernoulliEven); k++ {
		sum += complex(bernoulliEven[k-1]/factorial(2*k), 0) * term
		term *= (s + complex(float64(2*k-1), 0)) * (s + complex(float64(2*k), 0)) / na2
	}
	return sum
}

// DirichletL returns the Dirichlet L-function L(s, χ) for a complex argument s
// (with s ≠ 1 when χ is principal), using the Hurwitz-zeta decomposition
// L(s, χ) = q^{-s} Σ_{a=1}^{q} χ(a) ζ(s, a/q).
func DirichletL(s complex128, chi DirichletCharacter) complex128 {
	q := chi.Q
	qf := complex(float64(q), 0)
	var sum complex128
	for a := int64(1); a <= q; a++ {
		ca := chi.Eval(a)
		if ca == 0 {
			continue
		}
		sum += ca * HurwitzZetaComplex(s, float64(a)/float64(q))
	}
	return cmplx.Pow(qf, -s) * sum
}

// LFunctionReal returns L(s, χ) for a real argument s, returning the value as a
// complex128 (which is real for real characters).
func LFunctionReal(s float64, chi DirichletCharacter) complex128 {
	return DirichletL(complex(s, 0), chi)
}

// DirichletLSeries returns the truncated Dirichlet series Σ_{n=1}^{terms}
// χ(n) n^{-s} for real s > 1. It converges slowly and is provided mainly for
// cross-checking DirichletL on the half-plane of absolute convergence.
func DirichletLSeries(s float64, chi DirichletCharacter, terms int) complex128 {
	var sum complex128
	for n := 1; n <= terms; n++ {
		c := chi.Eval(int64(n))
		if c == 0 {
			continue
		}
		sum += c * complex(math.Pow(float64(n), -s), 0)
	}
	return sum
}

// DirichletBeta returns the Dirichlet beta function β(s) = Σ_{n≥0} (−1)^n
// (2n+1)^{-s}, the L-function of the non-principal character modulo 4. β(2) is
// Catalan's constant.
func DirichletBeta(s float64) float64 {
	if s == 1 {
		// β is entire; the Hurwitz decomposition has a removable pole here.
		return math.Pi / 4
	}
	// β(s) = 4^{-s} (ζ(s, 1/4) − ζ(s, 3/4)).
	return math.Pow(4, -s) * (HurwitzZeta(s, 0.25) - HurwitzZeta(s, 0.75))
}

// DirichletLambda returns the Dirichlet lambda function λ(s) = Σ_{n≥0}
// (2n+1)^{-s} = (1 − 2^{-s}) ζ(s) for real s > 1.
func DirichletLambda(s float64) float64 {
	return (1 - math.Pow(2, -s)) * Zeta(s)
}

// GaussSum returns the Gauss sum τ(χ) = Σ_{a=1}^{q} χ(a) e^{2πi a/q} associated
// with a Dirichlet character χ modulo q.
func GaussSum(chi DirichletCharacter) complex128 {
	q := chi.Q
	var sum complex128
	for a := int64(1); a <= q; a++ {
		c := chi.Eval(a)
		if c == 0 {
			continue
		}
		sum += c * cmplx.Rect(1, 2*math.Pi*float64(a)/float64(q))
	}
	return sum
}

// ClassNumberDirichlet returns the class number h(−p) of the imaginary quadratic
// field Q(√−p) for a prime p ≡ 3 (mod 4), computed from the quadratic character
// via the analytic class-number formula h = (1/(2−(2|p))) Σ_{a=1}^{(p−1)/2}
// (a|p). It is an integer for these p.
func ClassNumberDirichlet(p int64) int64 {
	if p%4 != 3 || !IsPrime(p) {
		panic("analyticnt: ClassNumberDirichlet requires a prime p ≡ 3 (mod 4)")
	}
	var sum int64
	for a := int64(1); a <= (p-1)/2; a++ {
		sum += int64(LegendreSymbol(a, p))
	}
	denom := int64(2 - LegendreSymbol(2, p))
	return sum / denom
}
