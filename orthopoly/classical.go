// Package orthopoly implements the classical families of orthogonal
// polynomials together with the tools needed to work with them: pointwise
// evaluation via numerically stable three-term recurrences, derivatives,
// weight functions, squared norms, the roots (zeros) of each family, and the
// Gaussian quadrature nodes and weights derived from those roots.
//
// The families covered are Legendre, Chebyshev (first, second, third and
// fourth kinds), Hermite (both the physicists' Hn and the probabilists' Hen
// conventions), Laguerre and generalized (associated) Laguerre, along with the
// Gegenbauer (ultraspherical) and Jacobi polynomials that unify them.
//
// Evaluation routines are built on the forward recurrences for each family and
// are accurate across the natural domain of the polynomial. Roots and Gaussian
// quadrature nodes are computed either in closed form (Chebyshev) or by
// Newton's method with well-conditioned starting estimates (Legendre, Hermite,
// Laguerre). All routines use only the Go standard library and are
// deterministic.
package orthopoly

import "math"

// orthopolyHermitePiM4 is pi^(-1/4), the leading coefficient of the
// L2-normalized physicists' Hermite recurrence used by the Gauss-Hermite node
// solver.
const orthopolyHermitePiM4 = 0.7511255444649425

// orthopolyFactorialFloat returns n! as a float64. It is used for norm
// constants and overflows to +Inf for large n, which is acceptable for the
// moderate degrees used in practice.
func orthopolyFactorialFloat(n int) float64 {
	r := 1.0
	for i := 2; i <= n; i++ {
		r *= float64(i)
	}
	return r
}

// orthopolyRequireNonNeg panics with a package-qualified message when n is
// negative, guarding the degree argument shared by every family.
func orthopolyRequireNonNeg(fn string, n int) {
	if n < 0 {
		panic("orthopoly: " + fn + " requires n >= 0")
	}
}

// ---------------------------------------------------------------------------
// Legendre
// ---------------------------------------------------------------------------

// orthopolyLegendrePair returns Pn(x) and P(n-1)(x) evaluated with the standard
// three-term recurrence (n+1)P(n+1) = (2n+1)x Pn - n P(n-1).
func orthopolyLegendrePair(n int, x float64) (pn, pnm1 float64) {
	if n == 0 {
		return 1, 0
	}
	pnm1 = 1
	pn = x
	for k := 1; k < n; k++ {
		pnm1, pn = pn, ((2*float64(k)+1)*x*pn-float64(k)*pnm1)/float64(k+1)
	}
	return pn, pnm1
}

// LegendreP returns the Legendre polynomial P_n(x). The P_n are orthogonal on
// [-1, 1] with unit weight. n must be non-negative.
func LegendreP(n int, x float64) float64 {
	orthopolyRequireNonNeg("LegendreP", n)
	pn, _ := orthopolyLegendrePair(n, x)
	return pn
}

// LegendrePValues returns the slice [P_0(x), P_1(x), ..., P_n(x)] of Legendre
// polynomials up to degree n. n must be non-negative.
func LegendrePValues(n int, x float64) []float64 {
	orthopolyRequireNonNeg("LegendrePValues", n)
	out := make([]float64, n+1)
	out[0] = 1
	if n == 0 {
		return out
	}
	out[1] = x
	for k := 1; k < n; k++ {
		out[k+1] = ((2*float64(k)+1)*x*out[k] - float64(k)*out[k-1]) / float64(k+1)
	}
	return out
}

// LegendrePDerivative returns the derivative P_n'(x) of the Legendre polynomial.
// The endpoints x = ±1, where the closed form P_n'(±1) = (±1)^(n+1) n(n+1)/2
// applies, are handled without cancellation. n must be non-negative.
func LegendrePDerivative(n int, x float64) float64 {
	orthopolyRequireNonNeg("LegendrePDerivative", n)
	if n == 0 {
		return 0
	}
	if x == 1 {
		return float64(n) * float64(n+1) / 2
	}
	if x == -1 {
		s := 1.0
		if (n+1)%2 != 0 {
			s = -1
		}
		return s * float64(n) * float64(n+1) / 2
	}
	pn, pnm1 := orthopolyLegendrePair(n, x)
	return float64(n) * (x*pn - pnm1) / (x*x - 1)
}

// NormalizedLegendreP returns the L2-orthonormal Legendre polynomial
// sqrt((2n+1)/2) P_n(x), whose square integrates to 1 over [-1, 1].
func NormalizedLegendreP(n int, x float64) float64 {
	return math.Sqrt((2*float64(n)+1)/2) * LegendreP(n, x)
}

// ShiftedLegendreP returns the shifted Legendre polynomial P_n(2x-1), which is
// orthogonal on [0, 1] with unit weight. n must be non-negative.
func ShiftedLegendreP(n int, x float64) float64 {
	return LegendreP(n, 2*x-1)
}

// LegendrePNorm returns the squared L2 norm of P_n on [-1, 1], namely
// integral of P_n(x)^2 dx = 2/(2n+1).
func LegendrePNorm(n int) float64 {
	return 2 / (2*float64(n) + 1)
}

// orthopolyGaussLegendre computes the n-point Gauss-Legendre nodes and weights
// on [-1, 1] by Newton iteration on the roots of P_n.
func orthopolyGaussLegendre(n int) (nodes, weights []float64) {
	if n <= 0 {
		panic("orthopoly: Gauss-Legendre requires n >= 1")
	}
	nodes = make([]float64, n)
	weights = make([]float64, n)
	m := (n + 1) / 2
	for i := 0; i < m; i++ {
		z := math.Cos(math.Pi * (float64(i) + 0.75) / (float64(n) + 0.5))
		var pp float64
		for iter := 0; iter < 100; iter++ {
			p1, p2 := 1.0, 0.0
			for j := 0; j < n; j++ {
				p1, p2 = ((2*float64(j)+1)*z*p1-float64(j)*p2)/float64(j+1), p1
			}
			pp = float64(n) * (z*p1 - p2) / (z*z - 1)
			z1 := z
			z = z1 - p1/pp
			if math.Abs(z-z1) <= 1e-15 {
				break
			}
		}
		nodes[i] = -z
		nodes[n-1-i] = z
		w := 2 / ((1 - z*z) * pp * pp)
		weights[i] = w
		weights[n-1-i] = w
	}
	return nodes, weights
}

// LegendreZeros returns the n roots of the Legendre polynomial P_n in ascending
// order (the Gauss-Legendre nodes). n must be positive.
func LegendreZeros(n int) []float64 {
	nodes, _ := orthopolyGaussLegendre(n)
	return nodes
}

// LegendreWeights returns the n Gauss-Legendre quadrature weights corresponding
// to LegendreZeros(n). n must be positive.
func LegendreWeights(n int) []float64 {
	_, weights := orthopolyGaussLegendre(n)
	return weights
}

// ---------------------------------------------------------------------------
// Chebyshev, first kind
// ---------------------------------------------------------------------------

// ChebyshevFirst returns the Chebyshev polynomial of the first kind T_n(x),
// which satisfies T_n(cos θ) = cos(nθ) and is orthogonal on [-1, 1] with weight
// 1/sqrt(1-x^2). n must be non-negative.
func ChebyshevFirst(n int, x float64) float64 {
	orthopolyRequireNonNeg("ChebyshevFirst", n)
	if n == 0 {
		return 1
	}
	tkm1, tk := 1.0, x
	for k := 1; k < n; k++ {
		tkm1, tk = tk, 2*x*tk-tkm1
	}
	return tk
}

// ChebyshevFirstValues returns the slice [T_0(x), ..., T_n(x)] of first-kind
// Chebyshev polynomials up to degree n. n must be non-negative.
func ChebyshevFirstValues(n int, x float64) []float64 {
	orthopolyRequireNonNeg("ChebyshevFirstValues", n)
	out := make([]float64, n+1)
	out[0] = 1
	if n == 0 {
		return out
	}
	out[1] = x
	for k := 1; k < n; k++ {
		out[k+1] = 2*x*out[k] - out[k-1]
	}
	return out
}

// ChebyshevFirstDerivative returns the derivative T_n'(x) = n U_(n-1)(x), where
// U is the Chebyshev polynomial of the second kind. n must be non-negative.
func ChebyshevFirstDerivative(n int, x float64) float64 {
	orthopolyRequireNonNeg("ChebyshevFirstDerivative", n)
	if n == 0 {
		return 0
	}
	return float64(n) * ChebyshevSecond(n-1, x)
}

// ChebyshevFirstRoots returns the n roots of T_n, x_k = cos((2k-1)π/(2n)), in
// ascending order. These are the Chebyshev-Gauss nodes. n must be positive.
func ChebyshevFirstRoots(n int) []float64 {
	if n <= 0 {
		panic("orthopoly: ChebyshevFirstRoots requires n >= 1")
	}
	r := make([]float64, n)
	for k := 1; k <= n; k++ {
		r[n-k] = math.Cos(math.Pi * (2*float64(k) - 1) / (2 * float64(n)))
	}
	return r
}

// ChebyshevFirstExtrema returns the n+1 extrema of T_n on [-1, 1],
// x_k = cos(kπ/n) for k = 0..n, in ascending order. These are the
// Chebyshev-Gauss-Lobatto points. n must be positive.
func ChebyshevFirstExtrema(n int) []float64 {
	if n <= 0 {
		panic("orthopoly: ChebyshevFirstExtrema requires n >= 1")
	}
	r := make([]float64, n+1)
	for k := 0; k <= n; k++ {
		r[n-k] = math.Cos(math.Pi * float64(k) / float64(n))
	}
	return r
}

// ChebyshevFirstWeight returns the orthogonality weight of the first-kind
// Chebyshev polynomials, 1/sqrt(1-x^2), for x in (-1, 1). It returns +Inf at
// x = ±1.
func ChebyshevFirstWeight(x float64) float64 {
	return 1 / math.Sqrt(1-x*x)
}

// ChebyshevFirstNorm returns the squared weighted L2 norm of T_n on [-1, 1],
// namely π for n = 0 and π/2 for n > 0.
func ChebyshevFirstNorm(n int) float64 {
	if n == 0 {
		return math.Pi
	}
	return math.Pi / 2
}

// ShiftedChebyshevFirst returns the shifted first-kind Chebyshev polynomial
// T_n(2x-1), orthogonal on [0, 1] with weight 1/sqrt(x-x^2). n must be
// non-negative.
func ShiftedChebyshevFirst(n int, x float64) float64 {
	return ChebyshevFirst(n, 2*x-1)
}

// ClenshawChebyshevFirst evaluates the Chebyshev series sum_{k} c[k] T_k(x)
// using Clenshaw's recurrence, which is more stable than summing evaluated
// polynomials. An empty coefficient slice yields 0.
func ClenshawChebyshevFirst(c []float64, x float64) float64 {
	var b1, b2 float64
	for k := len(c) - 1; k >= 1; k-- {
		b1, b2 = 2*x*b1-b2+c[k], b1
	}
	if len(c) == 0 {
		return 0
	}
	return x*b1 - b2 + c[0]
}

// ---------------------------------------------------------------------------
// Chebyshev, second kind
// ---------------------------------------------------------------------------

// ChebyshevSecond returns the Chebyshev polynomial of the second kind U_n(x),
// which satisfies U_n(cos θ) sin θ = sin((n+1)θ) and is orthogonal on [-1, 1]
// with weight sqrt(1-x^2). n must be non-negative.
func ChebyshevSecond(n int, x float64) float64 {
	orthopolyRequireNonNeg("ChebyshevSecond", n)
	if n == 0 {
		return 1
	}
	ukm1, uk := 1.0, 2*x
	for k := 1; k < n; k++ {
		ukm1, uk = uk, 2*x*uk-ukm1
	}
	return uk
}

// ChebyshevSecondValues returns the slice [U_0(x), ..., U_n(x)] of second-kind
// Chebyshev polynomials up to degree n. n must be non-negative.
func ChebyshevSecondValues(n int, x float64) []float64 {
	orthopolyRequireNonNeg("ChebyshevSecondValues", n)
	out := make([]float64, n+1)
	out[0] = 1
	if n == 0 {
		return out
	}
	out[1] = 2 * x
	for k := 1; k < n; k++ {
		out[k+1] = 2*x*out[k] - out[k-1]
	}
	return out
}

// ChebyshevSecondDerivative returns the derivative U_n'(x). The endpoints
// x = ±1, where U_n'(±1) = (±1)^(n+1) n(n+1)(n+2)/3, are handled in closed
// form. n must be non-negative.
func ChebyshevSecondDerivative(n int, x float64) float64 {
	orthopolyRequireNonNeg("ChebyshevSecondDerivative", n)
	if n == 0 {
		return 0
	}
	if x == 1 {
		return float64(n) * float64(n+1) * float64(n+2) / 3
	}
	if x == -1 {
		s := 1.0
		if (n+1)%2 != 0 {
			s = -1
		}
		return s * float64(n) * float64(n+1) * float64(n+2) / 3
	}
	return (float64(n+1)*ChebyshevFirst(n+1, x) - x*ChebyshevSecond(n, x)) / (x*x - 1)
}

// ChebyshevSecondRoots returns the n roots of U_n, x_k = cos(kπ/(n+1)) for
// k = 1..n, in ascending order. n must be positive.
func ChebyshevSecondRoots(n int) []float64 {
	if n <= 0 {
		panic("orthopoly: ChebyshevSecondRoots requires n >= 1")
	}
	r := make([]float64, n)
	for k := 1; k <= n; k++ {
		r[n-k] = math.Cos(math.Pi * float64(k) / float64(n+1))
	}
	return r
}

// ChebyshevSecondWeight returns the orthogonality weight of the second-kind
// Chebyshev polynomials, sqrt(1-x^2), for x in [-1, 1].
func ChebyshevSecondWeight(x float64) float64 {
	return math.Sqrt(1 - x*x)
}

// ChebyshevSecondNorm returns the squared weighted L2 norm of U_n on [-1, 1],
// which is π/2 for every n >= 0.
func ChebyshevSecondNorm(n int) float64 {
	return math.Pi / 2
}

// ---------------------------------------------------------------------------
// Chebyshev, third and fourth kinds
// ---------------------------------------------------------------------------

// ChebyshevThird returns the Chebyshev polynomial of the third kind V_n(x), for
// which V_n(cos θ) = cos((n+1/2)θ)/cos(θ/2). It is orthogonal on [-1, 1] with
// weight sqrt((1+x)/(1-x)). n must be non-negative.
func ChebyshevThird(n int, x float64) float64 {
	orthopolyRequireNonNeg("ChebyshevThird", n)
	if n == 0 {
		return 1
	}
	vkm1, vk := 1.0, 2*x-1
	for k := 1; k < n; k++ {
		vkm1, vk = vk, 2*x*vk-vkm1
	}
	return vk
}

// ChebyshevFourth returns the Chebyshev polynomial of the fourth kind W_n(x),
// for which W_n(cos θ) = sin((n+1/2)θ)/sin(θ/2). It is orthogonal on [-1, 1]
// with weight sqrt((1-x)/(1+x)). n must be non-negative.
func ChebyshevFourth(n int, x float64) float64 {
	orthopolyRequireNonNeg("ChebyshevFourth", n)
	if n == 0 {
		return 1
	}
	wkm1, wk := 1.0, 2*x+1
	for k := 1; k < n; k++ {
		wkm1, wk = wk, 2*x*wk-wkm1
	}
	return wk
}

// ---------------------------------------------------------------------------
// Hermite, physicists' convention (H)
// ---------------------------------------------------------------------------

// HermiteH returns the physicists' Hermite polynomial H_n(x), orthogonal on
// the real line with weight exp(-x^2) and satisfying the recurrence
// H_(n+1) = 2x H_n - 2n H_(n-1). n must be non-negative.
func HermiteH(n int, x float64) float64 {
	orthopolyRequireNonNeg("HermiteH", n)
	if n == 0 {
		return 1
	}
	hkm1, hk := 1.0, 2*x
	for k := 1; k < n; k++ {
		hkm1, hk = hk, 2*x*hk-2*float64(k)*hkm1
	}
	return hk
}

// HermiteHValues returns the slice [H_0(x), ..., H_n(x)] of physicists' Hermite
// polynomials up to degree n. n must be non-negative.
func HermiteHValues(n int, x float64) []float64 {
	orthopolyRequireNonNeg("HermiteHValues", n)
	out := make([]float64, n+1)
	out[0] = 1
	if n == 0 {
		return out
	}
	out[1] = 2 * x
	for k := 1; k < n; k++ {
		out[k+1] = 2*x*out[k] - 2*float64(k)*out[k-1]
	}
	return out
}

// HermiteHDerivative returns the derivative H_n'(x) = 2n H_(n-1)(x).
// n must be non-negative.
func HermiteHDerivative(n int, x float64) float64 {
	orthopolyRequireNonNeg("HermiteHDerivative", n)
	if n == 0 {
		return 0
	}
	return 2 * float64(n) * HermiteH(n-1, x)
}

// HermiteHWeight returns the orthogonality weight exp(-x^2) of the physicists'
// Hermite polynomials.
func HermiteHWeight(x float64) float64 {
	return math.Exp(-x * x)
}

// HermiteHNorm returns the squared weighted L2 norm of H_n over the real line,
// namely sqrt(π) 2^n n!.
func HermiteHNorm(n int) float64 {
	return math.Sqrt(math.Pi) * math.Ldexp(orthopolyFactorialFloat(n), n)
}

// NormalizedHermiteH returns the L2-orthonormal physicists' Hermite coefficient
// H_n(x)/sqrt(sqrt(π) 2^n n!).
func NormalizedHermiteH(n int, x float64) float64 {
	return HermiteH(n, x) / math.Sqrt(HermiteHNorm(n))
}

// HermiteFunction returns the Hermite function psi_n(x) =
// exp(-x^2/2) H_n(x) / sqrt(2^n n! sqrt(π)), the orthonormal eigenfunctions of
// the quantum harmonic oscillator. It is evaluated stably via the normalized
// recurrence, avoiding overflow of H_n for large x. n must be non-negative.
func HermiteFunction(n int, x float64) float64 {
	orthopolyRequireNonNeg("HermiteFunction", n)
	base := math.Pow(math.Pi, -0.25) * math.Exp(-x*x/2)
	if n == 0 {
		return base
	}
	pkm1, pk := 0.0, base
	for k := 0; k < n; k++ {
		pkm1, pk = pk, x*math.Sqrt(2/float64(k+1))*pk-math.Sqrt(float64(k)/float64(k+1))*pkm1
	}
	return pk
}

// orthopolyGaussHermite computes the n-point Gauss-Hermite nodes and weights
// for weight exp(-x^2) by Newton iteration on the L2-normalized recurrence.
func orthopolyGaussHermite(n int) (nodes, weights []float64) {
	if n <= 0 {
		panic("orthopoly: Gauss-Hermite requires n >= 1")
	}
	nodes = make([]float64, n)
	weights = make([]float64, n)
	m := (n + 1) / 2
	var z float64
	for i := 0; i < m; i++ {
		switch i {
		case 0:
			z = math.Sqrt(float64(2*n+1)) - 1.85575*math.Pow(float64(2*n+1), -1.0/6.0)
		case 1:
			z -= 1.14 * math.Pow(float64(n), 0.426) / z
		case 2:
			z = 1.86*z - 0.86*nodes[n-1]
		case 3:
			z = 1.91*z - 0.91*nodes[n-2]
		default:
			z = 2*z - nodes[n-i+1]
		}
		var pp float64
		for iter := 0; iter < 100; iter++ {
			p1, p2 := orthopolyHermitePiM4, 0.0
			for j := 1; j <= n; j++ {
				p1, p2 = z*math.Sqrt(2/float64(j))*p1-math.Sqrt(float64(j-1)/float64(j))*p2, p1
			}
			pp = math.Sqrt(2*float64(n)) * p2
			z1 := z
			z = z1 - p1/pp
			if math.Abs(z-z1) <= 1e-14 {
				break
			}
		}
		w := 2 / (pp * pp)
		nodes[n-1-i] = z
		nodes[i] = -z
		weights[n-1-i] = w
		weights[i] = w
	}
	return nodes, weights
}

// HermiteHRoots returns the n roots of the physicists' Hermite polynomial H_n
// in ascending order (the Gauss-Hermite nodes). n must be positive.
func HermiteHRoots(n int) []float64 {
	nodes, _ := orthopolyGaussHermite(n)
	return nodes
}

// HermiteHWeights returns the n Gauss-Hermite quadrature weights corresponding
// to HermiteHRoots(n). n must be positive.
func HermiteHWeights(n int) []float64 {
	_, weights := orthopolyGaussHermite(n)
	return weights
}

// ---------------------------------------------------------------------------
// Hermite, probabilists' convention (He)
// ---------------------------------------------------------------------------

// HermiteHe returns the probabilists' Hermite polynomial He_n(x), orthogonal
// on the real line with weight exp(-x^2/2) and satisfying the recurrence
// He_(n+1) = x He_n - n He_(n-1). n must be non-negative.
func HermiteHe(n int, x float64) float64 {
	orthopolyRequireNonNeg("HermiteHe", n)
	if n == 0 {
		return 1
	}
	hkm1, hk := 1.0, x
	for k := 1; k < n; k++ {
		hkm1, hk = hk, x*hk-float64(k)*hkm1
	}
	return hk
}

// HermiteHeValues returns the slice [He_0(x), ..., He_n(x)] of probabilists'
// Hermite polynomials up to degree n. n must be non-negative.
func HermiteHeValues(n int, x float64) []float64 {
	orthopolyRequireNonNeg("HermiteHeValues", n)
	out := make([]float64, n+1)
	out[0] = 1
	if n == 0 {
		return out
	}
	out[1] = x
	for k := 1; k < n; k++ {
		out[k+1] = x*out[k] - float64(k)*out[k-1]
	}
	return out
}

// HermiteHeDerivative returns the derivative He_n'(x) = n He_(n-1)(x).
// n must be non-negative.
func HermiteHeDerivative(n int, x float64) float64 {
	orthopolyRequireNonNeg("HermiteHeDerivative", n)
	if n == 0 {
		return 0
	}
	return float64(n) * HermiteHe(n-1, x)
}

// HermiteHeWeight returns the orthogonality weight exp(-x^2/2) of the
// probabilists' Hermite polynomials.
func HermiteHeWeight(x float64) float64 {
	return math.Exp(-x * x / 2)
}

// HermiteHeNorm returns the squared weighted L2 norm of He_n over the real
// line, namely sqrt(2π) n!.
func HermiteHeNorm(n int) float64 {
	return math.Sqrt(2*math.Pi) * orthopolyFactorialFloat(n)
}

// NormalizedHermiteHe returns the L2-orthonormal probabilists' Hermite
// coefficient He_n(x)/sqrt(sqrt(2π) n!).
func NormalizedHermiteHe(n int, x float64) float64 {
	return HermiteHe(n, x) / math.Sqrt(HermiteHeNorm(n))
}

// HermiteHeRoots returns the n roots of the probabilists' Hermite polynomial
// He_n in ascending order. Because He_n(x) = 2^(-n/2) H_n(x/sqrt(2)), its roots
// are the Gauss-Hermite nodes scaled by sqrt(2). n must be positive.
func HermiteHeRoots(n int) []float64 {
	nodes, _ := orthopolyGaussHermite(n)
	r := make([]float64, n)
	for i, z := range nodes {
		r[i] = z * math.Sqrt2
	}
	return r
}

// ---------------------------------------------------------------------------
// Laguerre and generalized Laguerre
// ---------------------------------------------------------------------------

// GeneralizedLaguerreL returns the generalized (associated) Laguerre polynomial
// L_n^(alpha)(x), orthogonal on [0, ∞) with weight x^alpha exp(-x). It satisfies
// (n+1) L_(n+1) = (2n+1+alpha-x) L_n - (n+alpha) L_(n-1). alpha must exceed -1
// and n must be non-negative.
func GeneralizedLaguerreL(n int, alpha, x float64) float64 {
	orthopolyRequireNonNeg("GeneralizedLaguerreL", n)
	if n == 0 {
		return 1
	}
	lkm1, lk := 1.0, 1+alpha-x
	for k := 1; k < n; k++ {
		fk := float64(k)
		lkm1, lk = lk, ((2*fk+1+alpha-x)*lk-(fk+alpha)*lkm1)/(fk+1)
	}
	return lk
}

// GeneralizedLaguerreLValues returns the slice [L_0^(alpha)(x), ...,
// L_n^(alpha)(x)] up to degree n. alpha must exceed -1 and n must be
// non-negative.
func GeneralizedLaguerreLValues(n int, alpha, x float64) []float64 {
	orthopolyRequireNonNeg("GeneralizedLaguerreLValues", n)
	out := make([]float64, n+1)
	out[0] = 1
	if n == 0 {
		return out
	}
	out[1] = 1 + alpha - x
	for k := 1; k < n; k++ {
		fk := float64(k)
		out[k+1] = ((2*fk+1+alpha-x)*out[k] - (fk+alpha)*out[k-1]) / (fk + 1)
	}
	return out
}

// GeneralizedLaguerreLDerivative returns the derivative
// d/dx L_n^(alpha)(x) = -L_(n-1)^(alpha+1)(x). n must be non-negative.
func GeneralizedLaguerreLDerivative(n int, alpha, x float64) float64 {
	orthopolyRequireNonNeg("GeneralizedLaguerreLDerivative", n)
	if n == 0 {
		return 0
	}
	return -GeneralizedLaguerreL(n-1, alpha+1, x)
}

// GeneralizedLaguerreWeight returns the orthogonality weight x^alpha exp(-x) of
// the generalized Laguerre polynomials, for x >= 0.
func GeneralizedLaguerreWeight(alpha, x float64) float64 {
	return math.Pow(x, alpha) * math.Exp(-x)
}

// GeneralizedLaguerreNorm returns the squared weighted L2 norm of L_n^(alpha)
// on [0, ∞), namely Γ(n+alpha+1)/n!.
func GeneralizedLaguerreNorm(n int, alpha float64) float64 {
	lg, _ := math.Lgamma(float64(n) + alpha + 1)
	return math.Exp(lg) / orthopolyFactorialFloat(n)
}

// LaguerreL returns the (ordinary) Laguerre polynomial L_n(x) = L_n^(0)(x),
// orthogonal on [0, ∞) with weight exp(-x). n must be non-negative.
func LaguerreL(n int, x float64) float64 {
	return GeneralizedLaguerreL(n, 0, x)
}

// LaguerreLValues returns the slice [L_0(x), ..., L_n(x)] of ordinary Laguerre
// polynomials up to degree n. n must be non-negative.
func LaguerreLValues(n int, x float64) []float64 {
	return GeneralizedLaguerreLValues(n, 0, x)
}

// LaguerreLDerivative returns the derivative L_n'(x) = -L_(n-1)^(1)(x) of the
// ordinary Laguerre polynomial. n must be non-negative.
func LaguerreLDerivative(n int, x float64) float64 {
	return GeneralizedLaguerreLDerivative(n, 0, x)
}

// LaguerreWeightFunction returns the orthogonality weight exp(-x) of the
// ordinary Laguerre polynomials, for x >= 0.
func LaguerreWeightFunction(x float64) float64 {
	return math.Exp(-x)
}

// LaguerreNorm returns the squared weighted L2 norm of L_n on [0, ∞), which is
// 1 for every n >= 0.
func LaguerreNorm(n int) float64 {
	return 1
}

// orthopolyGaussLaguerre computes the n-point generalized Gauss-Laguerre nodes
// and weights for weight x^alpha exp(-x) by Newton iteration on the roots of
// L_n^(alpha).
func orthopolyGaussLaguerre(n int, alpha float64) (nodes, weights []float64) {
	if n <= 0 {
		panic("orthopoly: Gauss-Laguerre requires n >= 1")
	}
	nodes = make([]float64, n)
	weights = make([]float64, n)
	lgAlfN, _ := math.Lgamma(alpha + float64(n))
	lgN, _ := math.Lgamma(float64(n))
	scale := math.Exp(lgAlfN - lgN)
	var z float64
	for i := 0; i < n; i++ {
		switch i {
		case 0:
			z = (1 + alpha) * (3 + 0.92*alpha) / (1 + 2.4*float64(n) + 1.8*alpha)
		case 1:
			z += (15 + 6.25*alpha) / (1 + 0.9*alpha + 2.5*float64(n))
		default:
			ai := float64(i - 1)
			z += ((1+2.55*ai)/(1.9*ai) + 1.26*ai*alpha/(1+3.5*ai)) * (z - nodes[i-2]) / (1 + 0.3*alpha)
		}
		var pp, p2 float64
		for iter := 0; iter < 100; iter++ {
			p1 := 1.0
			p2 = 0.0
			for j := 1; j <= n; j++ {
				fj := float64(j)
				p1, p2 = ((2*fj-1+alpha-z)*p1-(fj-1+alpha)*p2)/fj, p1
			}
			pp = (float64(n)*p1 - (float64(n)+alpha)*p2) / z
			z1 := z
			z = z1 - p1/pp
			if math.Abs(z-z1) <= 1e-14 {
				break
			}
		}
		nodes[i] = z
		weights[i] = -scale / (pp * float64(n) * p2)
	}
	return nodes, weights
}

// LaguerreRoots returns the n roots of the ordinary Laguerre polynomial L_n in
// ascending order (the Gauss-Laguerre nodes). n must be positive.
func LaguerreRoots(n int) []float64 {
	nodes, _ := orthopolyGaussLaguerre(n, 0)
	return nodes
}

// LaguerreWeights returns the n Gauss-Laguerre quadrature weights corresponding
// to LaguerreRoots(n). n must be positive.
func LaguerreWeights(n int) []float64 {
	_, weights := orthopolyGaussLaguerre(n, 0)
	return weights
}

// GeneralizedLaguerreRoots returns the n roots of L_n^(alpha) in ascending
// order. alpha must exceed -1 and n must be positive.
func GeneralizedLaguerreRoots(n int, alpha float64) []float64 {
	nodes, _ := orthopolyGaussLaguerre(n, alpha)
	return nodes
}

// GeneralizedLaguerreWeights returns the n generalized Gauss-Laguerre weights
// corresponding to GeneralizedLaguerreRoots(n, alpha). alpha must exceed -1 and
// n must be positive.
func GeneralizedLaguerreWeights(n int, alpha float64) []float64 {
	_, weights := orthopolyGaussLaguerre(n, alpha)
	return weights
}

// ---------------------------------------------------------------------------
// Gegenbauer (ultraspherical) and Jacobi
// ---------------------------------------------------------------------------

// GegenbauerC returns the Gegenbauer (ultraspherical) polynomial C_n^(alpha)(x),
// orthogonal on [-1, 1] with weight (1-x^2)^(alpha-1/2). It generalizes the
// Legendre (alpha = 1/2) and second-kind Chebyshev (alpha = 1) polynomials and
// satisfies n C_n = 2(n+alpha-1) x C_(n-1) - (n+2alpha-2) C_(n-2). alpha should
// exceed -1/2 and n must be non-negative.
func GegenbauerC(n int, alpha, x float64) float64 {
	orthopolyRequireNonNeg("GegenbauerC", n)
	if n == 0 {
		return 1
	}
	ckm1, ck := 1.0, 2*alpha*x
	for k := 2; k <= n; k++ {
		fk := float64(k)
		ckm1, ck = ck, (2*(fk+alpha-1)*x*ck-(fk+2*alpha-2)*ckm1)/fk
	}
	return ck
}

// GegenbauerCValues returns the slice [C_0^(alpha)(x), ..., C_n^(alpha)(x)] up
// to degree n. alpha should exceed -1/2 and n must be non-negative.
func GegenbauerCValues(n int, alpha, x float64) []float64 {
	orthopolyRequireNonNeg("GegenbauerCValues", n)
	out := make([]float64, n+1)
	out[0] = 1
	if n == 0 {
		return out
	}
	out[1] = 2 * alpha * x
	for k := 2; k <= n; k++ {
		fk := float64(k)
		out[k] = (2*(fk+alpha-1)*x*out[k-1] - (fk+2*alpha-2)*out[k-2]) / fk
	}
	return out
}

// GegenbauerCDerivative returns the derivative
// d/dx C_n^(alpha)(x) = 2 alpha C_(n-1)^(alpha+1)(x). n must be non-negative.
func GegenbauerCDerivative(n int, alpha, x float64) float64 {
	orthopolyRequireNonNeg("GegenbauerCDerivative", n)
	if n == 0 {
		return 0
	}
	return 2 * alpha * GegenbauerC(n-1, alpha+1, x)
}

// JacobiP returns the Jacobi polynomial P_n^(a,b)(x), orthogonal on [-1, 1]
// with weight (1-x)^a (1+x)^b. It generalizes the Legendre (a = b = 0),
// Chebyshev, Gegenbauer and Laguerre limits. Both a and b must exceed -1 and n
// must be non-negative.
func JacobiP(n int, a, b, x float64) float64 {
	orthopolyRequireNonNeg("JacobiP", n)
	if n == 0 {
		return 1
	}
	pkm1 := 1.0
	pk := (a - b + (a+b+2)*x) / 2
	for k := 2; k <= n; k++ {
		fk := float64(k)
		c := 2 * fk * (fk + a + b) * (2*fk + a + b - 2)
		c1 := (2*fk + a + b - 1) * (a*a - b*b)
		c2 := (2*fk + a + b - 1) * (2*fk + a + b) * (2*fk + a + b - 2)
		c3 := 2 * (fk + a - 1) * (fk + b - 1) * (2*fk + a + b)
		pkm1, pk = pk, ((c1+c2*x)*pk-c3*pkm1)/c
	}
	return pk
}

// JacobiPValues returns the slice [P_0^(a,b)(x), ..., P_n^(a,b)(x)] up to degree
// n. Both a and b must exceed -1 and n must be non-negative.
func JacobiPValues(n int, a, b, x float64) []float64 {
	orthopolyRequireNonNeg("JacobiPValues", n)
	out := make([]float64, n+1)
	out[0] = 1
	if n == 0 {
		return out
	}
	out[1] = (a - b + (a+b+2)*x) / 2
	for k := 2; k <= n; k++ {
		fk := float64(k)
		c := 2 * fk * (fk + a + b) * (2*fk + a + b - 2)
		c1 := (2*fk + a + b - 1) * (a*a - b*b)
		c2 := (2*fk + a + b - 1) * (2*fk + a + b) * (2*fk + a + b - 2)
		c3 := 2 * (fk + a - 1) * (fk + b - 1) * (2*fk + a + b)
		out[k] = ((c1+c2*x)*out[k-1] - c3*out[k-2]) / c
	}
	return out
}

// JacobiPDerivative returns the derivative
// d/dx P_n^(a,b)(x) = (n+a+b+1)/2 * P_(n-1)^(a+1,b+1)(x). n must be non-negative.
func JacobiPDerivative(n int, a, b, x float64) float64 {
	orthopolyRequireNonNeg("JacobiPDerivative", n)
	if n == 0 {
		return 0
	}
	return (float64(n) + a + b + 1) / 2 * JacobiP(n-1, a+1, b+1, x)
}

// ---------------------------------------------------------------------------
// Quadrature rule helper
// ---------------------------------------------------------------------------

// QuadratureRule bundles a set of Gaussian quadrature nodes with their weights
// so that a weighted integral can be approximated as
// sum_i Weights[i] f(Nodes[i]). The weight function is implied by the family
// that produced the rule.
type QuadratureRule struct {
	// Nodes are the evaluation abscissae (the polynomial roots).
	Nodes []float64
	// Weights are the quadrature weights aligned with Nodes.
	Weights []float64
}

// Order returns the number of nodes in the rule.
func (q QuadratureRule) Order() int {
	return len(q.Nodes)
}

// Integrate approximates the weighted integral of f by evaluating f at each node
// and forming the weighted sum sum_i Weights[i] f(Nodes[i]).
func (q QuadratureRule) Integrate(f func(float64) float64) float64 {
	var s float64
	for i, x := range q.Nodes {
		s += q.Weights[i] * f(x)
	}
	return s
}

// NewGaussLegendreRule returns the n-point Gauss-Legendre rule for unit weight
// on [-1, 1]. The rule integrates polynomials up to degree 2n-1 exactly.
// n must be positive.
func NewGaussLegendreRule(n int) QuadratureRule {
	nodes, weights := orthopolyGaussLegendre(n)
	return QuadratureRule{Nodes: nodes, Weights: weights}
}

// NewGaussHermiteRule returns the n-point Gauss-Hermite rule for weight
// exp(-x^2) on the real line. n must be positive.
func NewGaussHermiteRule(n int) QuadratureRule {
	nodes, weights := orthopolyGaussHermite(n)
	return QuadratureRule{Nodes: nodes, Weights: weights}
}

// NewGaussLaguerreRule returns the n-point Gauss-Laguerre rule for weight
// exp(-x) on [0, ∞). n must be positive.
func NewGaussLaguerreRule(n int) QuadratureRule {
	nodes, weights := orthopolyGaussLaguerre(n, 0)
	return QuadratureRule{Nodes: nodes, Weights: weights}
}

// NewGeneralizedLaguerreRule returns the n-point generalized Gauss-Laguerre
// rule for weight x^alpha exp(-x) on [0, ∞). alpha must exceed -1 and n must be
// positive.
func NewGeneralizedLaguerreRule(n int, alpha float64) QuadratureRule {
	nodes, weights := orthopolyGaussLaguerre(n, alpha)
	return QuadratureRule{Nodes: nodes, Weights: weights}
}
