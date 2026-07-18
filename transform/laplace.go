package transform

import (
	"math"
	"math/cmplx"
)

// Laplace numerically approximates the Laplace transform
//
//	F(s) = integral_0^upper f(t) e^{-s t} dt
//
// of the real function f at the complex point s, using composite Simpson
// quadrature with n subintervals over the truncated interval [0, upper]. For
// convergent transforms upper should be large enough that f(t)e^{-s t} has
// decayed to negligible size; n must be positive and is rounded up to the next
// even number.
func Laplace(f func(float64) float64, s complex128, upper float64, n int) complex128 {
	if n < 2 {
		n = 2
	}
	if n&1 == 1 {
		n++
	}
	h := upper / float64(n)
	g := func(t float64) complex128 {
		return complex(f(t), 0) * cmplx.Exp(-s*complex(t, 0))
	}
	sum := g(0) + g(upper)
	for i := 1; i < n; i++ {
		t := float64(i) * h
		if i&1 == 1 {
			sum += 4 * g(t)
		} else {
			sum += 2 * g(t)
		}
	}
	return sum * complex(h/3, 0)
}

// InverseLaplaceTalbot approximates the inverse Laplace transform f(t) of the
// transform F, supplied as a callable, using the fixed Talbot method of Abate
// and Valko. The contour is deformed into the left half-plane so that the
// integrand decays rapidly; m controls the number of terms (and thus the
// accuracy), with values around 20-40 giving several correct digits for
// smooth transforms. The argument t must be positive.
func InverseLaplaceTalbot(F func(complex128) complex128, t float64, m int) float64 {
	if m < 1 {
		m = 1
	}
	r := 2 * float64(m) / (5 * t)
	// k = 0 term.
	s0 := complex(r, 0)
	sum := 0.5 * real(cmplx.Exp(s0*complex(t, 0))*F(s0))
	for k := 1; k < m; k++ {
		theta := float64(k) * math.Pi / float64(m)
		cot := 1 / math.Tan(theta)
		s := complex(r*theta*cot, r*theta)
		sigma := theta + (theta*cot-1)*cot
		term := cmplx.Exp(s*complex(t, 0)) * F(s) * complex(1, sigma)
		sum += real(term)
	}
	return r / float64(m) * sum
}

// transformFactorial returns k! as a float64 for small non-negative k.
func transformFactorial(k int) float64 {
	if k < 0 {
		return math.NaN()
	}
	r := 1.0
	for i := 2; i <= k; i++ {
		r *= float64(i)
	}
	return r
}

// StehfestCoefficients returns the n Gaver-Stehfest weights V_1..V_n (indexed
// from zero) used by [InverseLaplaceStehfest]. n must be a positive even
// number; the function panics otherwise. The weights alternate in sign and
// grow rapidly, which is why the method is best used in the range n = 8..16.
func StehfestCoefficients(n int) []float64 {
	if n <= 0 || n&1 == 1 {
		panic("transform: StehfestCoefficients requires a positive even n")
	}
	nn := n / 2
	v := make([]float64, n)
	for k := 1; k <= n; k++ {
		lo := (k + 1) / 2
		hi := k
		if nn < hi {
			hi = nn
		}
		var sum float64
		for j := lo; j <= hi; j++ {
			num := math.Pow(float64(j), float64(nn)) * transformFactorial(2*j)
			den := transformFactorial(nn-j) * transformFactorial(j) *
				transformFactorial(j-1) * transformFactorial(k-j) *
				transformFactorial(2*j-k)
			sum += num / den
		}
		v[k-1] = transformSignPow(nn+k) * sum
	}
	return v
}

// InverseLaplaceStehfest approximates the inverse Laplace transform f(t) using
// the Gaver-Stehfest algorithm, which evaluates the real-argument transform F
// at n points:
//
//	f(t) ~= (ln 2 / t) * sum_{k=1}^{n} V_k * F(k ln 2 / t).
//
// n must be even (odd values are rounded up). The method needs only
// real-valued transform samples and works well for smooth, non-oscillatory
// functions; typical choices are n = 10..14.
func InverseLaplaceStehfest(F func(float64) float64, t float64, n int) float64 {
	if n < 2 {
		n = 2
	}
	if n&1 == 1 {
		n++
	}
	v := StehfestCoefficients(n)
	ln2 := math.Ln2
	var sum float64
	for k := 1; k <= n; k++ {
		sum += v[k-1] * F(float64(k)*ln2/t)
	}
	return ln2 / t * sum
}

// transformBinom returns the binomial coefficient C(n, k) as a float64.
func transformBinom(n, k int) float64 {
	if k < 0 || k > n {
		return 0
	}
	if k > n-k {
		k = n - k
	}
	r := 1.0
	for i := 0; i < k; i++ {
		r = r * float64(n-i) / float64(i+1)
	}
	return r
}

// InverseLaplaceEuler approximates the inverse Laplace transform f(t) using
// the Fourier-series method with Euler summation (Abate and Whitt). The
// Bromwich integral is written as an alternating series which is accelerated
// by averaging partial sums with binomial weights; m sets the number of
// Euler-averaged terms (around 15 is a good default). The transform F may be
// complex-valued and t must be positive.
func InverseLaplaceEuler(F func(complex128) complex128, t float64, m int) float64 {
	if m < 1 {
		m = 15
	}
	const a = 18.4
	const ntr = 15
	term := func(k int) float64 {
		s := complex(a/(2*t), math.Pi*float64(k)/t)
		re := real(F(s))
		switch {
		case k == 0:
			return re / 2
		case k&1 == 0:
			return re
		default:
			return -re
		}
	}
	sums := make([]float64, m+1)
	acc := 0.0
	for k := 0; k <= ntr; k++ {
		acc += term(k)
	}
	sums[0] = acc
	for j := 1; j <= m; j++ {
		acc += term(ntr + j)
		sums[j] = acc
	}
	var avg float64
	for j := 0; j <= m; j++ {
		avg += transformBinom(m, j) * sums[j]
	}
	avg *= math.Pow(2, -float64(m))
	return math.Exp(a/2) / t * avg
}
