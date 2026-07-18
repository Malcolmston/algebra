package orthopoly

import (
	"math"
	"sort"
)

// -----------------------------------------------------------------------------
// Unexported helpers. Every helper in this file is prefixed with "orthopoly"
// so the two source files that make up this package can never collide on
// private helper names.
// -----------------------------------------------------------------------------

// orthopolyLgamma returns the natural logarithm of the absolute value of the
// gamma function at x.
func orthopolyLgamma(x float64) float64 {
	v, _ := math.Lgamma(x)
	return v
}

// orthopolyFactorial returns n! as a float64. It returns NaN for negative n.
func orthopolyFactorial(n int) float64 {
	if n < 0 {
		return math.NaN()
	}
	r := 1.0
	for i := 2; i <= n; i++ {
		r *= float64(i)
	}
	return r
}

// orthopolyBinomial returns the binomial coefficient C(n,k) as a float64.
func orthopolyBinomial(n, k int) float64 {
	if k < 0 || k > n || n < 0 {
		return 0
	}
	if k > n-k {
		k = n - k
	}
	r := 1.0
	for i := 0; i < k; i++ {
		r = r * float64(n-i) / float64(i+1)
	}
	return math.Round(r)
}

// orthopolySign returns |a| with the sign of b (the Fortran/EISPACK SIGN).
func orthopolySign(a, b float64) float64 {
	return math.Copysign(a, b)
}

// orthopolyPlgndr evaluates the associated Legendre function P_l^m(x) for the
// canonical range 0 <= m <= l including the Condon-Shortley phase (-1)^m. The
// argument x must lie in [-1, 1].
func orthopolyPlgndr(l, m int, x float64) float64 {
	if m < 0 || m > l {
		return math.NaN()
	}
	pmm := 1.0
	if m > 0 {
		somx2 := math.Sqrt((1 - x) * (1 + x))
		fact := 1.0
		for i := 1; i <= m; i++ {
			pmm *= -fact * somx2
			fact += 2.0
		}
	}
	if l == m {
		return pmm
	}
	pmmp1 := x * float64(2*m+1) * pmm
	if l == m+1 {
		return pmmp1
	}
	var pll float64
	for ll := m + 2; ll <= l; ll++ {
		pll = (x*float64(2*ll-1)*pmmp1 - float64(ll+m-1)*pmm) / float64(ll-m)
		pmm = pmmp1
		pmmp1 = pll
	}
	return pll
}

// orthopolyTQLI diagonalises a real symmetric tridiagonal matrix using the QL
// algorithm with implicit shifts (the EISPACK/Numerical Recipes tqli routine).
// On entry d holds the diagonal and e[1..n-1] the subdiagonal (e[0] arbitrary);
// z holds the first row of the accumulating eigenvector matrix (initialised to
// the first row of the identity). On exit d holds the eigenvalues and z the
// corresponding first eigenvector components.
func orthopolyTQLI(d, e, z []float64) {
	n := len(d)
	if n == 0 {
		return
	}
	for i := 1; i < n; i++ {
		e[i-1] = e[i]
	}
	e[n-1] = 0
	for l := 0; l < n; l++ {
		iter := 0
		for {
			var m int
			for m = l; m < n-1; m++ {
				dd := math.Abs(d[m]) + math.Abs(d[m+1])
				if math.Abs(e[m])+dd == dd {
					break
				}
			}
			if m == l {
				break
			}
			iter++
			if iter > 60 {
				break
			}
			g := (d[l+1] - d[l]) / (2 * e[l])
			r := math.Hypot(g, 1)
			g = d[m] - d[l] + e[l]/(g+orthopolySign(r, g))
			s, c := 1.0, 1.0
			p := 0.0
			var i int
			for i = m - 1; i >= l; i-- {
				f := s * e[i]
				b := c * e[i]
				r = math.Hypot(f, g)
				e[i+1] = r
				if r == 0 {
					d[i+1] -= p
					e[m] = 0
					break
				}
				s = f / r
				c = g / r
				g = d[i+1] - p
				r = (d[i]-g)*s + 2*c*b
				p = s * r
				d[i+1] = g + p
				g = c*r - b
				f = z[i+1]
				z[i+1] = s*z[i] + c*f
				z[i] = c*z[i] - s*f
			}
			if r == 0 && i >= l {
				continue
			}
			d[l] -= p
			e[l] = g
			e[m] = 0
		}
	}
}

// orthopolyGauss returns the Gauss quadrature nodes and weights defined by the
// monic three-term recurrence with diagonal coefficients a (length n) and
// subdiagonal squares beta (length n, where beta[i] with i >= 1 is the monic
// recurrence coefficient beta_i). mu0 is the integral of the weight function.
func orthopolyGauss(a, beta []float64, mu0 float64) ([]float64, []float64) {
	n := len(a)
	d := make([]float64, n)
	copy(d, a)
	e := make([]float64, n)
	for i := 1; i < n; i++ {
		e[i] = math.Sqrt(beta[i])
	}
	z := make([]float64, n)
	if n > 0 {
		z[0] = 1
	}
	orthopolyTQLI(d, e, z)
	type nw struct{ x, w float64 }
	pairs := make([]nw, n)
	for i := 0; i < n; i++ {
		pairs[i] = nw{d[i], mu0 * z[i] * z[i]}
	}
	sort.Slice(pairs, func(i, j int) bool { return pairs[i].x < pairs[j].x })
	nodes := make([]float64, n)
	weights := make([]float64, n)
	for i := range pairs {
		nodes[i] = pairs[i].x
		weights[i] = pairs[i].w
	}
	return nodes, weights
}

// -----------------------------------------------------------------------------
// Gegenbauer (ultraspherical) polynomials
// -----------------------------------------------------------------------------

// Gegenbauer evaluates the ultraspherical (Gegenbauer) polynomial C_n^alpha(x)
// using the stable three-term recurrence
//
//	n*C_n = 2*(n+alpha-1)*x*C_{n-1} - (n+2*alpha-2)*C_{n-2},
//
// with C_0 = 1 and C_1 = 2*alpha*x.
func Gegenbauer(n int, alpha, x float64) float64 {
	if n < 0 {
		return math.NaN()
	}
	if n == 0 {
		return 1
	}
	cm := 1.0
	c := 2 * alpha * x
	if n == 1 {
		return c
	}
	for k := 2; k <= n; k++ {
		kf := float64(k)
		next := (2*(kf+alpha-1)*x*c - (kf+2*alpha-2)*cm) / kf
		cm, c = c, next
	}
	return c
}

// GegenbauerDerivative returns the derivative d/dx C_n^alpha(x), using the
// identity d/dx C_n^alpha(x) = 2*alpha*C_{n-1}^{alpha+1}(x).
func GegenbauerDerivative(n int, alpha, x float64) float64 {
	if n <= 0 {
		return 0
	}
	return 2 * alpha * Gegenbauer(n-1, alpha+1, x)
}

// GegenbauerWeight returns the weight function (1-x^2)^(alpha-1/2) for which the
// Gegenbauer polynomials are orthogonal on [-1, 1].
func GegenbauerWeight(alpha, x float64) float64 {
	return math.Pow(1-x*x, alpha-0.5)
}

// GegenbauerNorm returns the squared L2 norm of C_n^alpha under the Gegenbauer
// weight, that is the value of the integral of (1-x^2)^(alpha-1/2)*C_n^alpha(x)^2
// over [-1, 1]. It is valid for alpha > -1/2 and alpha != 0.
func GegenbauerNorm(n int, alpha float64) float64 {
	ln := math.Ln2*(1-2*alpha) + math.Log(math.Pi) +
		orthopolyLgamma(float64(n)+2*alpha) -
		orthopolyLgamma(float64(n)+1) -
		2*orthopolyLgamma(alpha)
	return math.Exp(ln) / (float64(n) + alpha)
}

// -----------------------------------------------------------------------------
// Jacobi polynomials
// -----------------------------------------------------------------------------

// Jacobi evaluates the Jacobi polynomial P_n^{a,b}(x) on [-1, 1] using the
// classical three-term recurrence. The parameters a and b must be greater than
// -1 for the polynomials to be orthogonal, but the recurrence is evaluated for
// any real a, b.
func Jacobi(n int, a, b, x float64) float64 {
	if n < 0 {
		return math.NaN()
	}
	if n == 0 {
		return 1
	}
	p0 := 1.0
	p1 := 0.5 * (a - b + (a+b+2)*x)
	if n == 1 {
		return p1
	}
	for k := 2; k <= n; k++ {
		kf := float64(k)
		c := 2*kf + a + b
		c1 := 2 * kf * (kf + a + b) * (c - 2)
		c2 := (c - 1) * (c*(c-2)*x + a*a - b*b)
		c3 := 2 * (kf + a - 1) * (kf + b - 1) * c
		p2 := (c2*p1 - c3*p0) / c1
		p0, p1 = p1, p2
	}
	return p1
}

// JacobiDerivative returns the derivative d/dx P_n^{a,b}(x), using the identity
// d/dx P_n^{a,b}(x) = (n+a+b+1)/2 * P_{n-1}^{a+1,b+1}(x).
func JacobiDerivative(n int, a, b, x float64) float64 {
	if n <= 0 {
		return 0
	}
	return 0.5 * (float64(n) + a + b + 1) * Jacobi(n-1, a+1, b+1, x)
}

// JacobiWeight returns the weight function (1-x)^a*(1+x)^b for which the Jacobi
// polynomials are orthogonal on [-1, 1].
func JacobiWeight(a, b, x float64) float64 {
	return math.Pow(1-x, a) * math.Pow(1+x, b)
}

// JacobiNorm returns the squared L2 norm of P_n^{a,b} under the Jacobi weight,
// i.e. the integral of (1-x)^a*(1+x)^b*P_n^{a,b}(x)^2 over [-1, 1]. Valid for
// a > -1 and b > -1.
func JacobiNorm(n int, a, b float64) float64 {
	nf := float64(n)
	ln := (a+b+1)*math.Ln2 +
		orthopolyLgamma(nf+a+1) + orthopolyLgamma(nf+b+1) -
		orthopolyLgamma(nf+a+b+1) - orthopolyLgamma(nf+1)
	return math.Exp(ln) / (2*nf + a + b + 1)
}

// -----------------------------------------------------------------------------
// Legendre polynomials
// -----------------------------------------------------------------------------

// Legendre evaluates the Legendre polynomial P_n(x) using Bonnet's recurrence
// (n+1)*P_{n+1} = (2n+1)*x*P_n - n*P_{n-1}.
func Legendre(n int, x float64) float64 {
	if n < 0 {
		return math.NaN()
	}
	if n == 0 {
		return 1
	}
	p0 := 1.0
	p1 := x
	for k := 1; k < n; k++ {
		kf := float64(k)
		p2 := ((2*kf+1)*x*p1 - kf*p0) / (kf + 1)
		p0, p1 = p1, p2
	}
	return p1
}

// LegendreDerivative returns the derivative P_n'(x). Interior points use the
// identity (x^2-1)*P_n'(x) = n*(x*P_n(x) - P_{n-1}(x)); the endpoints x = +/-1
// use the closed forms P_n'(1) = n(n+1)/2 and P_n'(-1) = (-1)^(n-1)*n(n+1)/2.
func LegendreDerivative(n int, x float64) float64 {
	if n <= 0 {
		return 0
	}
	nf := float64(n)
	if math.Abs(x-1) < 1e-14 {
		return nf * (nf + 1) / 2
	}
	if math.Abs(x+1) < 1e-14 {
		s := 1.0
		if n%2 == 0 {
			s = -1
		}
		return s * nf * (nf + 1) / 2
	}
	return nf * (x*Legendre(n, x) - Legendre(n-1, x)) / (x*x - 1)
}

// LegendreWeight returns the (constant) weight function 1 for which the
// Legendre polynomials are orthogonal on [-1, 1]. It is provided for API
// symmetry with the other weighted families.
func LegendreWeight(x float64) float64 {
	return 1
}

// LegendreNorm returns the squared L2 norm of P_n on [-1, 1], namely 2/(2n+1).
func LegendreNorm(n int) float64 {
	return 2 / (2*float64(n) + 1)
}

// LegendreBasis returns the slice [P_0(x), P_1(x), ..., P_n(x)] evaluated in a
// single stable pass of Bonnet's recurrence.
func LegendreBasis(n int, x float64) []float64 {
	if n < 0 {
		return nil
	}
	out := make([]float64, n+1)
	out[0] = 1
	if n == 0 {
		return out
	}
	out[1] = x
	for k := 1; k < n; k++ {
		kf := float64(k)
		out[k+1] = ((2*kf+1)*x*out[k] - kf*out[k-1]) / (kf + 1)
	}
	return out
}

// LegendreRoots returns the n roots of P_n(x), which coincide with the
// Gauss-Legendre quadrature nodes on [-1, 1], sorted in increasing order.
func LegendreRoots(n int) []float64 {
	nodes, _ := GaussLegendre(n)
	return nodes
}

// -----------------------------------------------------------------------------
// Associated Legendre functions and spherical harmonics
// -----------------------------------------------------------------------------

// AssociatedLegendre evaluates the associated Legendre function P_l^m(x) for
// -l <= m <= l on the interval [-1, 1], including the Condon-Shortley phase
// (-1)^m. Negative orders use the identity
// P_l^{-m} = (-1)^m*(l-m)!/(l+m)!*P_l^m.
func AssociatedLegendre(l, m int, x float64) float64 {
	if l < 0 || m < -l || m > l {
		return math.NaN()
	}
	if m >= 0 {
		return orthopolyPlgndr(l, m, x)
	}
	mm := -m
	factor := orthopolyFactorial(l-mm) / orthopolyFactorial(l+mm)
	if mm%2 == 1 {
		factor = -factor
	}
	return factor * orthopolyPlgndr(l, mm, x)
}

// NormalizedAssociatedLegendre returns the fully normalized associated Legendre
// function, sqrt((2l+1)/2 * (l-|m|)!/(l+|m|)!) * P_l^{|m|}(x), whose square
// integrates to 1 over [-1, 1] against the Legendre weight.
func NormalizedAssociatedLegendre(l, m int, x float64) float64 {
	if l < 0 {
		return math.NaN()
	}
	if m < 0 {
		m = -m
	}
	if m > l {
		return math.NaN()
	}
	lnNorm := 0.5 * (math.Log(2*float64(l)+1) - math.Ln2 +
		orthopolyLgamma(float64(l-m+1)) - orthopolyLgamma(float64(l+m+1)))
	return math.Exp(lnNorm) * orthopolyPlgndr(l, m, x)
}

// SphericalHarmonicNormalization returns the normalization constant
// sqrt((2l+1)/(4*pi) * (l-|m|)!/(l+|m|)!) used in the real and complex
// spherical harmonics of degree l and order m.
func SphericalHarmonicNormalization(l, m int) float64 {
	if m < 0 {
		m = -m
	}
	lnNorm := 0.5 * (math.Log(2*float64(l)+1) - math.Log(4*math.Pi) +
		orthopolyLgamma(float64(l-m+1)) - orthopolyLgamma(float64(l+m+1)))
	return math.Exp(lnNorm)
}

// RealSphericalHarmonic evaluates the real spherical harmonic Y_l^m(theta, phi),
// where theta is the polar (colatitude) angle in [0, pi] and phi the azimuthal
// angle. The real basis is defined without the Condon-Shortley phase as
//
//	m > 0:  sqrt(2)*N_l^m*cos(m*phi)*P_l^m(cos theta),
//	m = 0:  N_l^0*P_l^0(cos theta),
//	m < 0:  sqrt(2)*N_l^{|m|}*sin(|m|*phi)*P_l^{|m|}(cos theta),
//
// with N_l^m the SphericalHarmonicNormalization and P_l^m the associated
// Legendre function taken without the (-1)^m phase.
func RealSphericalHarmonic(l, m int, theta, phi float64) float64 {
	if l < 0 || m < -l || m > l {
		return math.NaN()
	}
	ct := math.Cos(theta)
	am := m
	if am < 0 {
		am = -am
	}
	// orthopolyPlgndr carries the Condon-Shortley phase (-1)^m; remove it.
	p := orthopolyPlgndr(l, am, ct)
	if am%2 == 1 {
		p = -p
	}
	n := SphericalHarmonicNormalization(l, am)
	switch {
	case m > 0:
		return math.Sqrt2 * n * math.Cos(float64(m)*phi) * p
	case m < 0:
		return math.Sqrt2 * n * math.Sin(float64(am)*phi) * p
	default:
		return n * p
	}
}

// -----------------------------------------------------------------------------
// Chebyshev polynomials
// -----------------------------------------------------------------------------

// ChebyshevT evaluates the Chebyshev polynomial of the first kind T_n(x) using
// the recurrence T_{n+1} = 2*x*T_n - T_{n-1}.
func ChebyshevT(n int, x float64) float64 {
	if n < 0 {
		n = -n
	}
	if n == 0 {
		return 1
	}
	t0 := 1.0
	t1 := x
	for k := 1; k < n; k++ {
		t0, t1 = t1, 2*x*t1-t0
	}
	return t1
}

// ChebyshevU evaluates the Chebyshev polynomial of the second kind U_n(x) using
// the recurrence U_{n+1} = 2*x*U_n - U_{n-1}, with U_0 = 1 and U_1 = 2*x.
func ChebyshevU(n int, x float64) float64 {
	if n < 0 {
		return math.NaN()
	}
	if n == 0 {
		return 1
	}
	u0 := 1.0
	u1 := 2 * x
	for k := 1; k < n; k++ {
		u0, u1 = u1, 2*x*u1-u0
	}
	return u1
}

// ChebyshevTDerivative returns T_n'(x) via the identity T_n'(x) = n*U_{n-1}(x).
func ChebyshevTDerivative(n int, x float64) float64 {
	if n <= 0 {
		return 0
	}
	return float64(n) * ChebyshevU(n-1, x)
}

// ChebyshevUDerivative returns U_n'(x). Interior points use
// U_n'(x) = ((n+1)*T_{n+1}(x) - x*U_n(x))/(x^2-1); the endpoints x = +/-1 use
// the closed forms U_n'(1) = n(n+1)(n+2)/3 and
// U_n'(-1) = (-1)^(n+1)*n(n+1)(n+2)/3.
func ChebyshevUDerivative(n int, x float64) float64 {
	if n <= 0 {
		return 0
	}
	nf := float64(n)
	if math.Abs(x-1) < 1e-14 {
		return nf * (nf + 1) * (nf + 2) / 3
	}
	if math.Abs(x+1) < 1e-14 {
		s := 1.0
		if n%2 == 0 {
			s = -1
		}
		return s * nf * (nf + 1) * (nf + 2) / 3
	}
	return (float64(n+1)*ChebyshevT(n+1, x) - x*ChebyshevU(n, x)) / (x*x - 1)
}

// ChebyshevTWeight returns the weight 1/sqrt(1-x^2) for which the Chebyshev
// polynomials of the first kind are orthogonal on (-1, 1).
func ChebyshevTWeight(x float64) float64 {
	return 1 / math.Sqrt(1-x*x)
}

// ChebyshevUWeight returns the weight sqrt(1-x^2) for which the Chebyshev
// polynomials of the second kind are orthogonal on [-1, 1].
func ChebyshevUWeight(x float64) float64 {
	return math.Sqrt(1 - x*x)
}

// ChebyshevTRoots returns the n roots of T_n, x_k = cos((2k-1)*pi/(2n)) for
// k = 1..n, sorted in increasing order.
func ChebyshevTRoots(n int) []float64 {
	if n <= 0 {
		return nil
	}
	out := make([]float64, n)
	for k := 1; k <= n; k++ {
		out[n-k] = math.Cos(float64(2*k-1) * math.Pi / float64(2*n))
	}
	return out
}

// ChebyshevTExtrema returns the n+1 Chebyshev-Gauss-Lobatto points
// x_k = cos(k*pi/n) for k = 0..n, the extrema of T_n, sorted in increasing
// order.
func ChebyshevTExtrema(n int) []float64 {
	if n < 0 {
		return nil
	}
	out := make([]float64, n+1)
	for k := 0; k <= n; k++ {
		out[n-k] = math.Cos(float64(k) * math.Pi / float64(n))
	}
	return out
}

// ChebyshevTBasis returns the slice [T_0(x), T_1(x), ..., T_n(x)] evaluated in a
// single pass of the recurrence.
func ChebyshevTBasis(n int, x float64) []float64 {
	if n < 0 {
		return nil
	}
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

// -----------------------------------------------------------------------------
// Hermite polynomials
// -----------------------------------------------------------------------------

// HermitePhysicists evaluates the physicists' Hermite polynomial H_n(x) with the
// recurrence H_{n+1} = 2*x*H_n - 2*n*H_{n-1}, orthogonal with weight exp(-x^2).
func HermitePhysicists(n int, x float64) float64 {
	if n < 0 {
		return math.NaN()
	}
	if n == 0 {
		return 1
	}
	h0 := 1.0
	h1 := 2 * x
	for k := 1; k < n; k++ {
		kf := float64(k)
		h0, h1 = h1, 2*x*h1-2*kf*h0
	}
	return h1
}

// HermiteProbabilists evaluates the probabilists' Hermite polynomial He_n(x)
// with the recurrence He_{n+1} = x*He_n - n*He_{n-1}, orthogonal with weight
// exp(-x^2/2).
func HermiteProbabilists(n int, x float64) float64 {
	if n < 0 {
		return math.NaN()
	}
	if n == 0 {
		return 1
	}
	h0 := 1.0
	h1 := x
	for k := 1; k < n; k++ {
		kf := float64(k)
		h0, h1 = h1, x*h1-kf*h0
	}
	return h1
}

// HermitePhysicistsDerivative returns H_n'(x) = 2*n*H_{n-1}(x).
func HermitePhysicistsDerivative(n int, x float64) float64 {
	if n <= 0 {
		return 0
	}
	return 2 * float64(n) * HermitePhysicists(n-1, x)
}

// HermiteProbabilistsDerivative returns He_n'(x) = n*He_{n-1}(x).
func HermiteProbabilistsDerivative(n int, x float64) float64 {
	if n <= 0 {
		return 0
	}
	return float64(n) * HermiteProbabilists(n-1, x)
}

// HermiteWeight returns the weight exp(-x^2) for which the physicists' Hermite
// polynomials are orthogonal on the whole real line.
func HermiteWeight(x float64) float64 {
	return math.Exp(-x * x)
}

// HermitePhysicistsNorm returns the squared L2 norm of H_n against exp(-x^2),
// namely sqrt(pi)*2^n*n!.
func HermitePhysicistsNorm(n int) float64 {
	return math.Sqrt(math.Pi) * math.Pow(2, float64(n)) * orthopolyFactorial(n)
}

// HermitePhysicistsBasis returns [H_0(x), H_1(x), ..., H_n(x)] evaluated in a
// single pass of the recurrence.
func HermitePhysicistsBasis(n int, x float64) []float64 {
	if n < 0 {
		return nil
	}
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

// -----------------------------------------------------------------------------
// Laguerre polynomials
// -----------------------------------------------------------------------------

// Laguerre evaluates the Laguerre polynomial L_n(x) using the recurrence
// (n+1)*L_{n+1} = (2n+1-x)*L_n - n*L_{n-1}, orthogonal with weight exp(-x) on
// [0, inf).
func Laguerre(n int, x float64) float64 {
	return AssociatedLaguerre(n, 0, x)
}

// AssociatedLaguerre evaluates the generalized Laguerre polynomial L_n^alpha(x)
// using the recurrence
// (n+1)*L_{n+1}^alpha = (2n+1+alpha-x)*L_n^alpha - (n+alpha)*L_{n-1}^alpha,
// orthogonal with weight x^alpha*exp(-x) on [0, inf).
func AssociatedLaguerre(n int, alpha, x float64) float64 {
	if n < 0 {
		return math.NaN()
	}
	if n == 0 {
		return 1
	}
	l0 := 1.0
	l1 := 1 + alpha - x
	for k := 1; k < n; k++ {
		kf := float64(k)
		l2 := ((2*kf+1+alpha-x)*l1 - (kf+alpha)*l0) / (kf + 1)
		l0, l1 = l1, l2
	}
	return l1
}

// LaguerreDerivative returns L_n'(x) = -L_{n-1}^1(x).
func LaguerreDerivative(n int, x float64) float64 {
	if n <= 0 {
		return 0
	}
	return -AssociatedLaguerre(n-1, 1, x)
}

// AssociatedLaguerreDerivative returns d/dx L_n^alpha(x) = -L_{n-1}^{alpha+1}(x).
func AssociatedLaguerreDerivative(n int, alpha, x float64) float64 {
	if n <= 0 {
		return 0
	}
	return -AssociatedLaguerre(n-1, alpha+1, x)
}

// LaguerreWeight returns the weight x^alpha*exp(-x) for which the generalized
// Laguerre polynomials are orthogonal on [0, inf).
func LaguerreWeight(alpha, x float64) float64 {
	return math.Pow(x, alpha) * math.Exp(-x)
}

// LaguerreBasis returns [L_0(x), L_1(x), ..., L_n(x)] evaluated in a single pass
// of the recurrence.
func LaguerreBasis(n int, x float64) []float64 {
	if n < 0 {
		return nil
	}
	out := make([]float64, n+1)
	out[0] = 1
	if n == 0 {
		return out
	}
	out[1] = 1 - x
	for k := 1; k < n; k++ {
		kf := float64(k)
		out[k+1] = ((2*kf+1-x)*out[k] - kf*out[k-1]) / (kf + 1)
	}
	return out
}

// -----------------------------------------------------------------------------
// Bernoulli and Euler polynomials
// -----------------------------------------------------------------------------

// BernoulliNumber returns the n-th Bernoulli number B_n in the convention with
// B_1 = -1/2, computed from the standard recurrence.
func BernoulliNumber(n int) float64 {
	if n < 0 {
		return math.NaN()
	}
	b := make([]float64, n+1)
	b[0] = 1
	for m := 1; m <= n; m++ {
		var s float64
		for k := 0; k < m; k++ {
			s += orthopolyBinomial(m+1, k) * b[k]
		}
		b[m] = -s / float64(m+1)
	}
	return b[n]
}

// BernoulliPolynomial evaluates the Bernoulli polynomial B_n(x) via
// B_n(x) = sum_{k=0}^n C(n,k)*B_{n-k}*x^k, where B_j are the Bernoulli numbers.
func BernoulliPolynomial(n int, x float64) float64 {
	if n < 0 {
		return math.NaN()
	}
	b := make([]float64, n+1)
	b[0] = 1
	for m := 1; m <= n; m++ {
		var s float64
		for k := 0; k < m; k++ {
			s += orthopolyBinomial(m+1, k) * b[k]
		}
		b[m] = -s / float64(m+1)
	}
	var sum, xp float64
	xp = 1
	for k := 0; k <= n; k++ {
		sum += orthopolyBinomial(n, k) * b[n-k] * xp
		xp *= x
	}
	return sum
}

// EulerPolynomial evaluates the Euler polynomial E_n(x) using its relation to
// the Bernoulli polynomials,
// E_n(x) = 2/(n+1) * (B_{n+1}(x) - 2^{n+1}*B_{n+1}(x/2)).
func EulerPolynomial(n int, x float64) float64 {
	if n < 0 {
		return math.NaN()
	}
	return 2 / float64(n+1) * (BernoulliPolynomial(n+1, x) -
		math.Pow(2, float64(n+1))*BernoulliPolynomial(n+1, x/2))
}

// EulerNumber returns the n-th Euler number E_n = 2^n*E_n(1/2). The odd-indexed
// Euler numbers are zero.
func EulerNumber(n int) float64 {
	if n < 0 {
		return math.NaN()
	}
	if n%2 == 1 {
		return 0
	}
	return math.Pow(2, float64(n)) * EulerPolynomial(n, 0.5)
}

// -----------------------------------------------------------------------------
// Bernstein basis
// -----------------------------------------------------------------------------

// Bernstein evaluates the Bernstein basis polynomial
// b_{i,n}(x) = C(n,i)*x^i*(1-x)^(n-i) for 0 <= i <= n.
func Bernstein(i, n int, x float64) float64 {
	if i < 0 || i > n || n < 0 {
		return 0
	}
	return orthopolyBinomial(n, i) * math.Pow(x, float64(i)) * math.Pow(1-x, float64(n-i))
}

// BernsteinDerivative returns the derivative of the Bernstein basis polynomial
// b_{i,n}'(x) = n*(b_{i-1,n-1}(x) - b_{i,n-1}(x)).
func BernsteinDerivative(i, n int, x float64) float64 {
	if n <= 0 {
		return 0
	}
	return float64(n) * (Bernstein(i-1, n-1, x) - Bernstein(i, n-1, x))
}

// BernsteinBasis returns the full degree-n Bernstein basis
// [b_{0,n}(x), ..., b_{n,n}(x)] evaluated at x.
func BernsteinBasis(n int, x float64) []float64 {
	if n < 0 {
		return nil
	}
	out := make([]float64, n+1)
	for i := 0; i <= n; i++ {
		out[i] = Bernstein(i, n, x)
	}
	return out
}

// BernsteinBezier evaluates the Bezier curve with the given control ordinates at
// parameter x in [0, 1] using the numerically stable de Casteljau algorithm.
// The polynomial degree equals len(coeffs)-1.
func BernsteinBezier(coeffs []float64, x float64) float64 {
	n := len(coeffs)
	if n == 0 {
		return math.NaN()
	}
	tmp := make([]float64, n)
	copy(tmp, coeffs)
	for r := 1; r < n; r++ {
		for i := 0; i < n-r; i++ {
			tmp[i] = (1-x)*tmp[i] + x*tmp[i+1]
		}
	}
	return tmp[0]
}

// -----------------------------------------------------------------------------
// Zernike polynomials
// -----------------------------------------------------------------------------

// ZernikeRadial evaluates the Zernike radial polynomial R_n^m(rho) for
// 0 <= m <= n on the unit disk. It returns 0 when n-m is odd, per the standard
// definition.
func ZernikeRadial(n, m int, rho float64) float64 {
	if n < 0 || m < 0 || m > n {
		return math.NaN()
	}
	if (n-m)%2 != 0 {
		return 0
	}
	half := (n - m) / 2
	var sum float64
	for k := 0; k <= half; k++ {
		coeff := orthopolyFactorial(n-k) /
			(orthopolyFactorial(k) *
				orthopolyFactorial((n+m)/2-k) *
				orthopolyFactorial((n-m)/2-k))
		if k%2 == 1 {
			coeff = -coeff
		}
		sum += coeff * math.Pow(rho, float64(n-2*k))
	}
	return sum
}

// ZernikeRadialDerivative returns the derivative d/drho R_n^m(rho), obtained by
// differentiating the explicit power-series definition term by term.
func ZernikeRadialDerivative(n, m int, rho float64) float64 {
	if n < 0 || m < 0 || m > n {
		return math.NaN()
	}
	if (n-m)%2 != 0 {
		return 0
	}
	half := (n - m) / 2
	var sum float64
	for k := 0; k <= half; k++ {
		p := n - 2*k
		if p == 0 {
			continue
		}
		coeff := orthopolyFactorial(n-k) /
			(orthopolyFactorial(k) *
				orthopolyFactorial((n+m)/2-k) *
				orthopolyFactorial((n-m)/2-k))
		if k%2 == 1 {
			coeff = -coeff
		}
		sum += coeff * float64(p) * math.Pow(rho, float64(p-1))
	}
	return sum
}

// ZernikeNormalization returns the normalization factor that makes the Zernike
// polynomial orthonormal over the unit disk: sqrt(2*(n+1)) for m != 0 and
// sqrt(n+1) for m == 0.
func ZernikeNormalization(n, m int) float64 {
	if m == 0 {
		return math.Sqrt(float64(n + 1))
	}
	return math.Sqrt(float64(2 * (n + 1)))
}

// Zernike evaluates the (unnormalized) Zernike polynomial Z_n^m(rho, theta) on
// the unit disk. Positive m selects the cosine azimuthal factor, negative m the
// sine factor, and m == 0 the purely radial polynomial.
func Zernike(n, m int, rho, theta float64) float64 {
	if n < 0 || m < -n || m > n {
		return math.NaN()
	}
	am := m
	if am < 0 {
		am = -am
	}
	r := ZernikeRadial(n, am, rho)
	switch {
	case m > 0:
		return r * math.Cos(float64(m)*theta)
	case m < 0:
		return r * math.Sin(float64(am)*theta)
	default:
		return r
	}
}

// NollToNM converts a 1-based Noll single index j to the Zernike degree n and
// signed order m. Noll indexing assigns j = 1 to piston, j = 2, 3 to tilt, and
// so on.
func NollToNM(j int) (n, m int) {
	if j < 1 {
		return -1, 0
	}
	j1 := j - 1
	n = 0
	for j1 > n {
		n++
		j1 -= n
	}
	m = (n % 2) + 2*((j1+((n+1)%2))/2)
	if j%2 == 1 {
		m = -m
	}
	return n, m
}

// -----------------------------------------------------------------------------
// Gauss quadrature rules
// -----------------------------------------------------------------------------

// GolubWelsch computes Gauss quadrature nodes and weights from the monic
// three-term recurrence of an orthogonal polynomial family via the
// Golub-Welsch algorithm. alpha holds the n diagonal recurrence coefficients
// (a_0..a_{n-1}); beta holds the n subdiagonal coefficients where beta[i], for
// i >= 1, is the monic recurrence coefficient beta_i (beta[0] is ignored). mu0
// is the zeroth moment (the integral of the weight function). The returned
// nodes are sorted in increasing order.
func GolubWelsch(alpha, beta []float64, mu0 float64) (nodes, weights []float64) {
	if len(alpha) == 0 || len(alpha) != len(beta) {
		return nil, nil
	}
	return orthopolyGauss(alpha, beta, mu0)
}

// GaussLegendre returns the n-point Gauss-Legendre quadrature nodes and weights
// on [-1, 1] (weight function 1). The rule integrates polynomials up to degree
// 2n-1 exactly.
func GaussLegendre(n int) (nodes, weights []float64) {
	if n <= 0 {
		return nil, nil
	}
	a := make([]float64, n)
	b := make([]float64, n)
	for k := 1; k < n; k++ {
		kf := float64(k)
		b[k] = kf * kf / (4*kf*kf - 1)
	}
	return orthopolyGauss(a, b, 2)
}

// GaussHermite returns the n-point Gauss-Hermite quadrature nodes and weights on
// (-inf, inf) with weight exp(-x^2). The rule integrates polynomials up to
// degree 2n-1 exactly.
func GaussHermite(n int) (nodes, weights []float64) {
	if n <= 0 {
		return nil, nil
	}
	a := make([]float64, n)
	b := make([]float64, n)
	for k := 1; k < n; k++ {
		b[k] = float64(k) / 2
	}
	return orthopolyGauss(a, b, math.Sqrt(math.Pi))
}

// GaussLaguerre returns the n-point Gauss-Laguerre quadrature nodes and weights
// on [0, inf) with weight exp(-x). The rule integrates polynomials up to degree
// 2n-1 exactly.
func GaussLaguerre(n int) (nodes, weights []float64) {
	return GaussLaguerreGeneralized(n, 0)
}

// GaussLaguerreGeneralized returns the n-point generalized Gauss-Laguerre
// quadrature nodes and weights on [0, inf) with weight x^alpha*exp(-x), for
// alpha > -1.
func GaussLaguerreGeneralized(n int, alpha float64) (nodes, weights []float64) {
	if n <= 0 {
		return nil, nil
	}
	a := make([]float64, n)
	b := make([]float64, n)
	for k := 0; k < n; k++ {
		a[k] = 2*float64(k) + 1 + alpha
	}
	for k := 1; k < n; k++ {
		kf := float64(k)
		b[k] = kf * (kf + alpha)
	}
	mu0 := math.Gamma(1 + alpha)
	return orthopolyGauss(a, b, mu0)
}

// GaussJacobi returns the n-point Gauss-Jacobi quadrature nodes and weights on
// [-1, 1] with weight (1-x)^a*(1+x)^b, for a > -1 and b > -1.
func GaussJacobi(n int, a, b float64) (nodes, weights []float64) {
	if n <= 0 {
		return nil, nil
	}
	diag := make([]float64, n)
	off := make([]float64, n)
	ab := a + b
	diag[0] = (b - a) / (ab + 2)
	for k := 1; k < n; k++ {
		kf := float64(k)
		den := (2*kf + ab) * (2*kf + ab + 2)
		diag[k] = (b*b - a*a) / den
	}
	if n > 1 {
		off[1] = 4 * (a + 1) * (b + 1) / ((ab + 2) * (ab + 2) * (ab + 3))
	}
	for k := 2; k < n; k++ {
		kf := float64(k)
		num := 4 * kf * (kf + a) * (kf + b) * (kf + ab)
		den := (2*kf + ab) * (2*kf + ab) * (2*kf + ab + 1) * (2*kf + ab - 1)
		off[k] = num / den
	}
	mu0 := math.Exp((ab+1)*math.Ln2 +
		orthopolyLgamma(a+1) + orthopolyLgamma(b+1) - orthopolyLgamma(ab+2))
	return orthopolyGauss(diag, off, mu0)
}

// GaussGegenbauer returns the n-point Gauss-Gegenbauer quadrature nodes and
// weights on [-1, 1] with weight (1-x^2)^(alpha-1/2), for alpha > -1/2. It is
// the symmetric Gauss-Jacobi rule with a = b = alpha - 1/2.
func GaussGegenbauer(n int, alpha float64) (nodes, weights []float64) {
	return GaussJacobi(n, alpha-0.5, alpha-0.5)
}

// GaussChebyshevT returns the n-point Gauss-Chebyshev quadrature of the first
// kind on [-1, 1] with weight 1/sqrt(1-x^2). Nodes are x_k = cos((2k-1)*pi/(2n))
// and every weight equals pi/n. The nodes are returned in increasing order.
func GaussChebyshevT(n int) (nodes, weights []float64) {
	if n <= 0 {
		return nil, nil
	}
	nodes = make([]float64, n)
	weights = make([]float64, n)
	w := math.Pi / float64(n)
	for k := 1; k <= n; k++ {
		nodes[n-k] = math.Cos(float64(2*k-1) * math.Pi / float64(2*n))
		weights[n-k] = w
	}
	return nodes, weights
}

// GaussChebyshevU returns the n-point Gauss-Chebyshev quadrature of the second
// kind on [-1, 1] with weight sqrt(1-x^2). Nodes are x_k = cos(k*pi/(n+1)) with
// weights pi/(n+1)*sin^2(k*pi/(n+1)). The nodes are returned in increasing
// order.
func GaussChebyshevU(n int) (nodes, weights []float64) {
	if n <= 0 {
		return nil, nil
	}
	nodes = make([]float64, n)
	weights = make([]float64, n)
	for k := 1; k <= n; k++ {
		ang := float64(k) * math.Pi / float64(n+1)
		nodes[n-k] = math.Cos(ang)
		s := math.Sin(ang)
		weights[n-k] = math.Pi / float64(n+1) * s * s
	}
	return nodes, weights
}

// GaussLegendreIntegrate approximates the integral of f over [lo, hi] using the
// n-point Gauss-Legendre rule, mapping the standard interval [-1, 1] to
// [lo, hi] with the affine change of variables.
func GaussLegendreIntegrate(f func(float64) float64, lo, hi float64, n int) float64 {
	nodes, weights := GaussLegendre(n)
	half := (hi - lo) / 2
	mid := (hi + lo) / 2
	var sum float64
	for i := range nodes {
		sum += weights[i] * f(mid+half*nodes[i])
	}
	return half * sum
}

// GaussHermiteIntegrate approximates the integral of f(x)*exp(-x^2) over the
// whole real line using the n-point Gauss-Hermite rule. f should be the factor
// multiplying the Gaussian weight.
func GaussHermiteIntegrate(f func(float64) float64, n int) float64 {
	nodes, weights := GaussHermite(n)
	var sum float64
	for i := range nodes {
		sum += weights[i] * f(nodes[i])
	}
	return sum
}

// GaussLaguerreIntegrate approximates the integral of f(x)*exp(-x) over
// [0, inf) using the n-point Gauss-Laguerre rule. f should be the factor
// multiplying the exponential weight.
func GaussLaguerreIntegrate(f func(float64) float64, n int) float64 {
	nodes, weights := GaussLaguerre(n)
	var sum float64
	for i := range nodes {
		sum += weights[i] * f(nodes[i])
	}
	return sum
}
