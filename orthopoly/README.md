# orthopoly — implements the classical families of orthogonal polynomials together with the tools

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/orthopoly.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/orthopoly)

Package orthopoly implements the classical families of orthogonal polynomials
together with the tools needed to work with them: pointwise evaluation via
numerically stable three-term recurrences, derivatives, weight functions,
squared norms, the roots (zeros) of each family, and the Gaussian quadrature
nodes and weights derived from those roots.

The families covered are Legendre, Chebyshev (first, second, third and fourth
kinds), Hermite (both the physicists' Hn and the probabilists' Hen
conventions), Laguerre and generalized (associated) Laguerre, along with the
Gegenbauer (ultraspherical) and Jacobi polynomials that unify them.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/orthopoly"
```

## Usage

This is the package's own `ExampleBernsteinBasis`, so it compiles and its output is
asserted on every `go test ./orthopoly/`.

```go
b := orthopoly.BernsteinBasis(3, 0.5)
	sum := 0.0
	for i, v := range b {
		if i > 0 {
			fmt.Print(" ")
		}
		fmt.Printf("%.4f", v)
		sum += v
	}
	fmt.Printf("\nsum=%.4f\n", sum)
```

```
0.1250 0.3750 0.3750 0.1250
sum=1.0000
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func AssociatedLaguerre(n int, alpha, x float64) float64` | AssociatedLaguerre evaluates the generalized Laguerre polynomial L_n^alpha(x) using the recurrence (n+1)*L_{n+1}^alpha = (2n+1+alpha-x)*L_n^alpha -… |
| `func AssociatedLaguerreDerivative(n int, alpha, x float64) float64` | AssociatedLaguerreDerivative returns d/dx L_n^alpha(x) = -L_{n-1}^{alpha+1}(x). |
| `func AssociatedLegendre(l, m int, x float64) float64` | AssociatedLegendre evaluates the associated Legendre function P_l^m(x) for -l <= m <= l on the interval [-1, 1], including the Condon-Shortley phase… |
| `func BernoulliNumber(n int) float64` | BernoulliNumber returns the n-th Bernoulli number B_n in the convention with B_1 = -1/2, computed from the standard recurrence. |
| `func BernoulliPolynomial(n int, x float64) float64` | BernoulliPolynomial evaluates the Bernoulli polynomial B_n(x) via B_n(x) = sum_{k=0}^n C(n,k)*B_{n-k}*x^k, where B_j are the Bernoulli numbers. |
| `func Bernstein(i, n int, x float64) float64` | Bernstein evaluates the Bernstein basis polynomial b_{i,n}(x) = C(n,i)*x^i*(1-x)^(n-i) for 0 <= i <= n. |
| `func BernsteinBasis(n int, x float64) []float64` | BernsteinBasis returns the full degree-n Bernstein basis [b_{0,n}(x), ..., b_{n,n}(x)] evaluated at x. |
| `func BernsteinBezier(coeffs []float64, x float64) float64` | BernsteinBezier evaluates the Bezier curve with the given control ordinates at parameter x in [0, 1] using the numerically stable de Casteljau… |
| `func BernsteinDerivative(i, n int, x float64) float64` | BernsteinDerivative returns the derivative of the Bernstein basis polynomial b_{i,n}'(x) = n*(b_{i-1,n-1}(x) - b_{i,n-1}(x)). |
| `func ChebyshevFirst(n int, x float64) float64` | ChebyshevFirst returns the Chebyshev polynomial of the first kind T_n(x), which satisfies T_n(cos θ) = cos(nθ) and is orthogonal on [-1, 1] with… |
| `func ChebyshevFirstDerivative(n int, x float64) float64` | ChebyshevFirstDerivative returns the derivative T_n'(x) = n U_(n-1)(x), where U is the Chebyshev polynomial of the second kind. |
| `func ChebyshevFirstExtrema(n int) []float64` | ChebyshevFirstExtrema returns the n+1 extrema of T_n on [-1, 1], x_k = cos(kπ/n) for k = 0..n, in ascending order. |
| `func ChebyshevFirstNorm(n int) float64` | ChebyshevFirstNorm returns the squared weighted L2 norm of T_n on [-1, 1], namely π for n = 0 and π/2 for n > 0. |
| `func ChebyshevFirstRoots(n int) []float64` | ChebyshevFirstRoots returns the n roots of T_n, x_k = cos((2k-1)π/(2n)), in ascending order. |
| `func ChebyshevFirstValues(n int, x float64) []float64` | ChebyshevFirstValues returns the slice [T_0(x), ..., T_n(x)] of first-kind Chebyshev polynomials up to degree n. |
| `func ChebyshevFirstWeight(x float64) float64` | ChebyshevFirstWeight returns the orthogonality weight of the first-kind Chebyshev polynomials, 1/sqrt(1-x^2), for x in (-1, 1). |
| `func ChebyshevFourth(n int, x float64) float64` | ChebyshevFourth returns the Chebyshev polynomial of the fourth kind W_n(x), for which W_n(cos θ) = sin((n+1/2)θ)/sin(θ/2). |
| `func ChebyshevSecond(n int, x float64) float64` | ChebyshevSecond returns the Chebyshev polynomial of the second kind U_n(x), which satisfies U_n(cos θ) sin θ = sin((n+1)θ) and is orthogonal on… |
| `func ChebyshevSecondDerivative(n int, x float64) float64` | ChebyshevSecondDerivative returns the derivative U_n'(x). |
| `func ChebyshevSecondNorm(n int) float64` | ChebyshevSecondNorm returns the squared weighted L2 norm of U_n on [-1, 1], which is π/2 for every n >= 0. |
| `func ChebyshevSecondRoots(n int) []float64` | ChebyshevSecondRoots returns the n roots of U_n, x_k = cos(kπ/(n+1)) for k = 1..n, in ascending order. |
| `func ChebyshevSecondValues(n int, x float64) []float64` | ChebyshevSecondValues returns the slice [U_0(x), ..., U_n(x)] of second-kind Chebyshev polynomials up to degree n. |
| `func ChebyshevSecondWeight(x float64) float64` | ChebyshevSecondWeight returns the orthogonality weight of the second-kind Chebyshev polynomials, sqrt(1-x^2), for x in [-1, 1]. |
| `func ChebyshevT(n int, x float64) float64` | ChebyshevT evaluates the Chebyshev polynomial of the first kind T_n(x) using the recurrence T_{n+1} = 2*x*T_n - T_{n-1}. |
| `func ChebyshevTBasis(n int, x float64) []float64` | ChebyshevTBasis returns the slice [T_0(x), T_1(x), ..., T_n(x)] evaluated in a single pass of the recurrence. |
| `func ChebyshevTDerivative(n int, x float64) float64` | ChebyshevTDerivative returns T_n'(x) via the identity T_n'(x) = n*U_{n-1}(x). |
| `func ChebyshevTExtrema(n int) []float64` | ChebyshevTExtrema returns the n+1 Chebyshev-Gauss-Lobatto points x_k = cos(k*pi/n) for k = 0..n, the extrema of T_n, sorted in increasing order. |
| `func ChebyshevTRoots(n int) []float64` | ChebyshevTRoots returns the n roots of T_n, x_k = cos((2k-1)*pi/(2n)) for k = 1..n, sorted in increasing order. |
| `func ChebyshevTWeight(x float64) float64` | ChebyshevTWeight returns the weight 1/sqrt(1-x^2) for which the Chebyshev polynomials of the first kind are orthogonal on (-1, 1). |
| `func ChebyshevThird(n int, x float64) float64` | ChebyshevThird returns the Chebyshev polynomial of the third kind V_n(x), for which V_n(cos θ) = cos((n+1/2)θ)/cos(θ/2). |
| `func ChebyshevU(n int, x float64) float64` | ChebyshevU evaluates the Chebyshev polynomial of the second kind U_n(x) using the recurrence U_{n+1} = 2*x*U_n - U_{n-1}, with U_0 = 1 and U_1 = 2*x. |
| `func ChebyshevUDerivative(n int, x float64) float64` | ChebyshevUDerivative returns U_n'(x). |
| `func ChebyshevUWeight(x float64) float64` | ChebyshevUWeight returns the weight sqrt(1-x^2) for which the Chebyshev polynomials of the second kind are orthogonal on [-1, 1]. |
| `func ClenshawChebyshevFirst(c []float64, x float64) float64` | ClenshawChebyshevFirst evaluates the Chebyshev series sum_{k} c`k` T_k(x) using Clenshaw's recurrence, which is more stable than summing evaluated… |
| `func EulerNumber(n int) float64` | EulerNumber returns the n-th Euler number E_n = 2^n*E_n(1/2). |
| `func EulerPolynomial(n int, x float64) float64` | EulerPolynomial evaluates the Euler polynomial E_n(x) using its relation to the Bernoulli polynomials, E_n(x) = 2/(n+1) * (B_{n+1}(x) -… |
| `func GaussChebyshevT(n int) (nodes, weights []float64)` | GaussChebyshevT returns the n-point Gauss-Chebyshev quadrature of the first kind on [-1, 1] with weight 1/sqrt(1-x^2). |
| `func GaussChebyshevU(n int) (nodes, weights []float64)` | GaussChebyshevU returns the n-point Gauss-Chebyshev quadrature of the second kind on [-1, 1] with weight sqrt(1-x^2). |
| `func GaussGegenbauer(n int, alpha float64) (nodes, weights []float64)` | GaussGegenbauer returns the n-point Gauss-Gegenbauer quadrature nodes and weights on [-1, 1] with weight (1-x^2)^(alpha-1/2), for alpha > -1/2. |
| `func GaussHermite(n int) (nodes, weights []float64)` | GaussHermite returns the n-point Gauss-Hermite quadrature nodes and weights on (-inf, inf) with weight exp(-x^2). |
| `func GaussHermiteIntegrate(f func(float64) float64, n int) float64` | GaussHermiteIntegrate approximates the integral of f(x)*exp(-x^2) over the whole real line using the n-point Gauss-Hermite rule. |
| `func GaussJacobi(n int, a, b float64) (nodes, weights []float64)` | GaussJacobi returns the n-point Gauss-Jacobi quadrature nodes and weights on [-1, 1] with weight (1-x)^a*(1+x)^b, for a > -1 and b > -1. |
| `func GaussLaguerre(n int) (nodes, weights []float64)` | GaussLaguerre returns the n-point Gauss-Laguerre quadrature nodes and weights on [0, inf) with weight exp(-x). |
| `func GaussLaguerreGeneralized(n int, alpha float64) (nodes, weights []float64)` | GaussLaguerreGeneralized returns the n-point generalized Gauss-Laguerre quadrature nodes and weights on [0, inf) with weight x^alpha*exp(-x), for… |
| `func GaussLaguerreIntegrate(f func(float64) float64, n int) float64` | GaussLaguerreIntegrate approximates the integral of f(x)*exp(-x) over [0, inf) using the n-point Gauss-Laguerre rule. |
| `func GaussLegendre(n int) (nodes, weights []float64)` | GaussLegendre returns the n-point Gauss-Legendre quadrature nodes and weights on [-1, 1] (weight function 1). |
| `func GaussLegendreIntegrate(f func(float64) float64, lo, hi float64, n int) float64` | GaussLegendreIntegrate approximates the integral of f over [lo, hi] using the n-point Gauss-Legendre rule, mapping the standard interval [-1, 1] to… |
| `func Gegenbauer(n int, alpha, x float64) float64` | Gegenbauer evaluates the ultraspherical (Gegenbauer) polynomial C_n^alpha(x) using the stable three-term recurrence |
| `func GegenbauerC(n int, alpha, x float64) float64` | GegenbauerC returns the Gegenbauer (ultraspherical) polynomial C_n^(alpha)(x), orthogonal on [-1, 1] with weight (1-x^2)^(alpha-1/2). |
| `func GegenbauerCDerivative(n int, alpha, x float64) float64` | GegenbauerCDerivative returns the derivative d/dx C_n^(alpha)(x) = 2 alpha C_(n-1)^(alpha+1)(x). |
| `func GegenbauerCValues(n int, alpha, x float64) []float64` | GegenbauerCValues returns the slice [C_0^(alpha)(x), ..., C_n^(alpha)(x)] up to degree n. |
| `func GegenbauerDerivative(n int, alpha, x float64) float64` | GegenbauerDerivative returns the derivative d/dx C_n^alpha(x), using the identity d/dx C_n^alpha(x) = 2*alpha*C_{n-1}^{alpha+1}(x). |
| `func GegenbauerNorm(n int, alpha float64) float64` | GegenbauerNorm returns the squared L2 norm of C_n^alpha under the Gegenbauer weight, that is the value of the integral of… |
| `func GegenbauerWeight(alpha, x float64) float64` | GegenbauerWeight returns the weight function (1-x^2)^(alpha-1/2) for which the Gegenbauer polynomials are orthogonal on [-1, 1]. |
| `func GeneralizedLaguerreL(n int, alpha, x float64) float64` | GeneralizedLaguerreL returns the generalized (associated) Laguerre polynomial L_n^(alpha)(x), orthogonal on [0, ∞) with weight x^alpha exp(-x). |
| `func GeneralizedLaguerreLDerivative(n int, alpha, x float64) float64` | GeneralizedLaguerreLDerivative returns the derivative d/dx L_n^(alpha)(x) = -L_(n-1)^(alpha+1)(x). |
| `func GeneralizedLaguerreLValues(n int, alpha, x float64) []float64` | GeneralizedLaguerreLValues returns the slice [L_0^(alpha)(x), ..., L_n^(alpha)(x)] up to degree n. |
| `func GeneralizedLaguerreNorm(n int, alpha float64) float64` | GeneralizedLaguerreNorm returns the squared weighted L2 norm of L_n^(alpha) on [0, ∞), namely Γ(n+alpha+1)/n!. |
| `func GeneralizedLaguerreRoots(n int, alpha float64) []float64` | GeneralizedLaguerreRoots returns the n roots of L_n^(alpha) in ascending order. |
| `func GeneralizedLaguerreWeight(alpha, x float64) float64` | GeneralizedLaguerreWeight returns the orthogonality weight x^alpha exp(-x) of the generalized Laguerre polynomials, for x >= 0. |
| `func GeneralizedLaguerreWeights(n int, alpha float64) []float64` | GeneralizedLaguerreWeights returns the n generalized Gauss-Laguerre weights corresponding to GeneralizedLaguerreRoots(n, alpha). |
| `func GolubWelsch(alpha, beta []float64, mu0 float64) (nodes, weights []float64)` | GolubWelsch computes Gauss quadrature nodes and weights from the monic three-term recurrence of an orthogonal polynomial family via the Golub-Welsch… |
| `func HermiteFunction(n int, x float64) float64` | HermiteFunction returns the Hermite function psi_n(x) = exp(-x^2/2) H_n(x) / sqrt(2^n n! sqrt(π)), the orthonormal eigenfunctions of the quantum… |
| `func HermiteH(n int, x float64) float64` | HermiteH returns the physicists' Hermite polynomial H_n(x), orthogonal on the real line with weight exp(-x^2) and satisfying the recurrence H_(n+1) =… |
| `func HermiteHDerivative(n int, x float64) float64` | HermiteHDerivative returns the derivative H_n'(x) = 2n H_(n-1)(x). |
| `func HermiteHNorm(n int) float64` | HermiteHNorm returns the squared weighted L2 norm of H_n over the real line, namely sqrt(π) 2^n n!. |
| `func HermiteHRoots(n int) []float64` | HermiteHRoots returns the n roots of the physicists' Hermite polynomial H_n in ascending order (the Gauss-Hermite nodes). |
| `func HermiteHValues(n int, x float64) []float64` | HermiteHValues returns the slice [H_0(x), ..., H_n(x)] of physicists' Hermite polynomials up to degree n. |
| `func HermiteHWeight(x float64) float64` | HermiteHWeight returns the orthogonality weight exp(-x^2) of the physicists' Hermite polynomials. |
| `func HermiteHWeights(n int) []float64` | HermiteHWeights returns the n Gauss-Hermite quadrature weights corresponding to HermiteHRoots(n). |
| `func HermiteHe(n int, x float64) float64` | HermiteHe returns the probabilists' Hermite polynomial He_n(x), orthogonal on the real line with weight exp(-x^2/2) and satisfying the recurrence… |
| `func HermiteHeDerivative(n int, x float64) float64` | HermiteHeDerivative returns the derivative He_n'(x) = n He_(n-1)(x). |
| `func HermiteHeNorm(n int) float64` | HermiteHeNorm returns the squared weighted L2 norm of He_n over the real line, namely sqrt(2π) n!. |
| `func HermiteHeRoots(n int) []float64` | HermiteHeRoots returns the n roots of the probabilists' Hermite polynomial He_n in ascending order. |
| `func HermiteHeValues(n int, x float64) []float64` | HermiteHeValues returns the slice [He_0(x), ..., He_n(x)] of probabilists' Hermite polynomials up to degree n. |
| `func HermiteHeWeight(x float64) float64` | HermiteHeWeight returns the orthogonality weight exp(-x^2/2) of the probabilists' Hermite polynomials. |
| `func HermitePhysicists(n int, x float64) float64` | HermitePhysicists evaluates the physicists' Hermite polynomial H_n(x) with the recurrence H_{n+1} = 2*x*H_n - 2*n*H_{n-1}, orthogonal with weight… |
| `func HermitePhysicistsBasis(n int, x float64) []float64` | HermitePhysicistsBasis returns [H_0(x), H_1(x), ..., H_n(x)] evaluated in a single pass of the recurrence. |
| `func HermitePhysicistsDerivative(n int, x float64) float64` | HermitePhysicistsDerivative returns H_n'(x) = 2*n*H_{n-1}(x). |
| `func HermitePhysicistsNorm(n int) float64` | HermitePhysicistsNorm returns the squared L2 norm of H_n against exp(-x^2), namely sqrt(pi)*2^n*n!. |
| `func HermiteProbabilists(n int, x float64) float64` | HermiteProbabilists evaluates the probabilists' Hermite polynomial He_n(x) with the recurrence He_{n+1} = x*He_n - n*He_{n-1}, orthogonal with weight… |
| `func HermiteProbabilistsDerivative(n int, x float64) float64` | HermiteProbabilistsDerivative returns He_n'(x) = n*He_{n-1}(x). |
| `func HermiteWeight(x float64) float64` | HermiteWeight returns the weight exp(-x^2) for which the physicists' Hermite polynomials are orthogonal on the whole real line. |
| `func Jacobi(n int, a, b, x float64) float64` | Jacobi evaluates the Jacobi polynomial P_n^{a,b}(x) on [-1, 1] using the classical three-term recurrence. |
| `func JacobiDerivative(n int, a, b, x float64) float64` | JacobiDerivative returns the derivative d/dx P_n^{a,b}(x), using the identity d/dx P_n^{a,b}(x) = (n+a+b+1)/2 * P_{n-1}^{a+1,b+1}(x). |
| `func JacobiNorm(n int, a, b float64) float64` | JacobiNorm returns the squared L2 norm of P_n^{a,b} under the Jacobi weight, i.e. |
| `func JacobiP(n int, a, b, x float64) float64` | JacobiP returns the Jacobi polynomial P_n^(a,b)(x), orthogonal on [-1, 1] with weight (1-x)^a (1+x)^b. |
| `func JacobiPDerivative(n int, a, b, x float64) float64` | JacobiPDerivative returns the derivative d/dx P_n^(a,b)(x) = (n+a+b+1)/2 * P_(n-1)^(a+1,b+1)(x). |
| `func JacobiPValues(n int, a, b, x float64) []float64` | JacobiPValues returns the slice [P_0^(a,b)(x), ..., P_n^(a,b)(x)] up to degree n. |
| `func JacobiWeight(a, b, x float64) float64` | JacobiWeight returns the weight function (1-x)^a*(1+x)^b for which the Jacobi polynomials are orthogonal on [-1, 1]. |
| `func Laguerre(n int, x float64) float64` | Laguerre evaluates the Laguerre polynomial L_n(x) using the recurrence (n+1)*L_{n+1} = (2n+1-x)*L_n - n*L_{n-1}, orthogonal with weight exp(-x) on… |
| `func LaguerreBasis(n int, x float64) []float64` | LaguerreBasis returns [L_0(x), L_1(x), ..., L_n(x)] evaluated in a single pass of the recurrence. |
| `func LaguerreDerivative(n int, x float64) float64` | LaguerreDerivative returns L_n'(x) = -L_{n-1}^1(x). |
| `func LaguerreL(n int, x float64) float64` | LaguerreL returns the (ordinary) Laguerre polynomial L_n(x) = L_n^(0)(x), orthogonal on [0, ∞) with weight exp(-x). |
| `func LaguerreLDerivative(n int, x float64) float64` | LaguerreLDerivative returns the derivative L_n'(x) = -L_(n-1)^(1)(x) of the ordinary Laguerre polynomial. |
| `func LaguerreLValues(n int, x float64) []float64` | LaguerreLValues returns the slice [L_0(x), ..., L_n(x)] of ordinary Laguerre polynomials up to degree n. |
| `func LaguerreNorm(n int) float64` | LaguerreNorm returns the squared weighted L2 norm of L_n on [0, ∞), which is 1 for every n >= 0. |
| `func LaguerreRoots(n int) []float64` | LaguerreRoots returns the n roots of the ordinary Laguerre polynomial L_n in ascending order (the Gauss-Laguerre nodes). |
| `func LaguerreWeight(alpha, x float64) float64` | LaguerreWeight returns the weight x^alpha*exp(-x) for which the generalized Laguerre polynomials are orthogonal on [0, inf). |
| `func LaguerreWeightFunction(x float64) float64` | LaguerreWeightFunction returns the orthogonality weight exp(-x) of the ordinary Laguerre polynomials, for x >= 0. |
| `func LaguerreWeights(n int) []float64` | LaguerreWeights returns the n Gauss-Laguerre quadrature weights corresponding to LaguerreRoots(n). |
| `func Legendre(n int, x float64) float64` | Legendre evaluates the Legendre polynomial P_n(x) using Bonnet's recurrence (n+1)*P_{n+1} = (2n+1)*x*P_n - n*P_{n-1}. |
| `func LegendreBasis(n int, x float64) []float64` | LegendreBasis returns the slice [P_0(x), P_1(x), ..., P_n(x)] evaluated in a single stable pass of Bonnet's recurrence. |
| `func LegendreDerivative(n int, x float64) float64` | LegendreDerivative returns the derivative P_n'(x). |
| `func LegendreNorm(n int) float64` | LegendreNorm returns the squared L2 norm of P_n on [-1, 1], namely 2/(2n+1). |
| `func LegendreP(n int, x float64) float64` | LegendreP returns the Legendre polynomial P_n(x). |
| `func LegendrePDerivative(n int, x float64) float64` | LegendrePDerivative returns the derivative P_n'(x) of the Legendre polynomial. |
| `func LegendrePNorm(n int) float64` | LegendrePNorm returns the squared L2 norm of P_n on [-1, 1], namely integral of P_n(x)^2 dx = 2/(2n+1). |
| `func LegendrePValues(n int, x float64) []float64` | LegendrePValues returns the slice [P_0(x), P_1(x), ..., P_n(x)] of Legendre polynomials up to degree n. |
| `func LegendreRoots(n int) []float64` | LegendreRoots returns the n roots of P_n(x), which coincide with the Gauss-Legendre quadrature nodes on [-1, 1], sorted in increasing order. |
| `func LegendreWeight(x float64) float64` | LegendreWeight returns the (constant) weight function 1 for which the Legendre polynomials are orthogonal on [-1, 1]. |
| `func LegendreWeights(n int) []float64` | LegendreWeights returns the n Gauss-Legendre quadrature weights corresponding to LegendreZeros(n). |
| `func LegendreZeros(n int) []float64` | LegendreZeros returns the n roots of the Legendre polynomial P_n in ascending order (the Gauss-Legendre nodes). |
| `func NollToNM(j int) (n, m int)` | NollToNM converts a 1-based Noll single index j to the Zernike degree n and signed order m. |
| `func NormalizedAssociatedLegendre(l, m int, x float64) float64` | NormalizedAssociatedLegendre returns the fully normalized associated Legendre function, sqrt((2l+1)/2 * (l-\|m\|)!/(l+\|m\|)!) * P_l^{\|m\|}(x),… |
| `func NormalizedHermiteH(n int, x float64) float64` | NormalizedHermiteH returns the L2-orthonormal physicists' Hermite coefficient H_n(x)/sqrt(sqrt(π) 2^n n!). |
| `func NormalizedHermiteHe(n int, x float64) float64` | NormalizedHermiteHe returns the L2-orthonormal probabilists' Hermite coefficient He_n(x)/sqrt(sqrt(2π) n!). |
| `func NormalizedLegendreP(n int, x float64) float64` | NormalizedLegendreP returns the L2-orthonormal Legendre polynomial sqrt((2n+1)/2) P_n(x), whose square integrates to 1 over [-1, 1]. |
| `func RealSphericalHarmonic(l, m int, theta, phi float64) float64` | RealSphericalHarmonic evaluates the real spherical harmonic Y_l^m(theta, phi), where theta is the polar (colatitude) angle in [0, pi] and phi the… |
| `func ShiftedChebyshevFirst(n int, x float64) float64` | ShiftedChebyshevFirst returns the shifted first-kind Chebyshev polynomial T_n(2x-1), orthogonal on [0, 1] with weight 1/sqrt(x-x^2). |
| `func ShiftedLegendreP(n int, x float64) float64` | ShiftedLegendreP returns the shifted Legendre polynomial P_n(2x-1), which is orthogonal on [0, 1] with unit weight. |
| `func SphericalHarmonicNormalization(l, m int) float64` | SphericalHarmonicNormalization returns the normalization constant sqrt((2l+1)/(4*pi) * (l-\|m\|)!/(l+\|m\|)!) used in the real and complex spherical… |
| `func Zernike(n, m int, rho, theta float64) float64` | Zernike evaluates the (unnormalized) Zernike polynomial Z_n^m(rho, theta) on the unit disk. |
| `func ZernikeNormalization(n, m int) float64` | ZernikeNormalization returns the normalization factor that makes the Zernike polynomial orthonormal over the unit disk: sqrt(2*(n+1)) for m != 0 and… |
| `func ZernikeRadial(n, m int, rho float64) float64` | ZernikeRadial evaluates the Zernike radial polynomial R_n^m(rho) for 0 <= m <= n on the unit disk. |
| `func ZernikeRadialDerivative(n, m int, rho float64) float64` | ZernikeRadialDerivative returns the derivative d/drho R_n^m(rho), obtained by differentiating the explicit power-series definition term by term. |

### Types

| Type | What it is |
| --- | --- |
| `QuadratureRule` | QuadratureRule bundles a set of Gaussian quadrature nodes with their weights so that a weighted integral can be approximated as sum_i Weights`i`… |

<details>
<summary><code>QuadratureRule</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewGaussHermiteRule(n int) QuadratureRule` | NewGaussHermiteRule returns the n-point Gauss-Hermite rule for weight exp(-x^2) on the real line. |
| `func NewGaussLaguerreRule(n int) QuadratureRule` | NewGaussLaguerreRule returns the n-point Gauss-Laguerre rule for weight exp(-x) on [0, ∞). |
| `func NewGaussLegendreRule(n int) QuadratureRule` | NewGaussLegendreRule returns the n-point Gauss-Legendre rule for unit weight on [-1, 1]. |
| `func NewGeneralizedLaguerreRule(n int, alpha float64) QuadratureRule` | NewGeneralizedLaguerreRule returns the n-point generalized Gauss-Laguerre rule for weight x^alpha exp(-x) on [0, ∞). |
| `func (q QuadratureRule) Integrate(f func(float64) float64) float64` | Integrate approximates the weighted integral of f by evaluating f at each node and forming the weighted sum sum_i Weights`i` f(Nodes`i`). |
| `func (q QuadratureRule) Order() int` | Order returns the number of nodes in the rule. |

</details>

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/orthopoly).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
