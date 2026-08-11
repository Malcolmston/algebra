# special — implements special mathematical functions in pure Go

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/special.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/special)

Package special implements special mathematical functions in pure Go.

- Bessel functions of the first and second kind (J0, J1, Jn, Y0, Y1, Yn), the
modified Bessel functions (I0, I1, In, K0, K1, Kn), the spherical Bessel
functions (SphericalJn, SphericalYn) and the Struve functions (StruveH0,
StruveH1). - Airy functions Ai, Bi and their derivatives. - Elliptic integrals
and related quantities. - The error function family, including the Fresnel
integrals. - The exponential, sine and cosine integrals. - The Riemann zeta
family (zeta, eta, beta). - The Lambert W function. - The gamma family
(log-gamma, digamma, incomplete forms).

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/special"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func Airy(x float64) (ai, aip, bi, bip float64)` | Airy returns the Airy functions Ai(x), Ai'(x), Bi(x) and Bi'(x) simultaneously. |
| `func AiryAi(x float64) float64` | AiryAi returns the Airy function of the first kind, Ai(x). |
| `func AiryAiPrime(x float64) float64` | AiryAiPrime returns the derivative of the Airy function of the first kind, Ai'(x). |
| `func AiryBi(x float64) float64` | AiryBi returns the Airy function of the second kind, Bi(x). |
| `func AiryBiPrime(x float64) float64` | AiryBiPrime returns the derivative of the Airy function of the second kind, Bi'(x). |
| `func AiryZeroAi(s int) float64` | AiryZeroAi returns the s-th real zero (s >= 1) of the Airy function Ai, accurate to full double precision. |
| `func AiryZeroBi(s int) float64` | AiryZeroBi returns the s-th real zero (s >= 1) of the Airy function Bi, accurate to full double precision. |
| `func BesselI(n int, x float64) float64` | BesselI returns the modified Bessel function of the first kind of integer order n, Iₙ(x). |
| `func BesselI0(x float64) float64` | BesselI0 returns the modified Bessel function of the first kind of order zero, I0(x). |
| `func BesselI0Prime(x float64) float64` | BesselI0Prime returns the derivative of I0, which equals I1(x). |
| `func BesselI0e(x float64) float64` | BesselI0e returns the exponentially scaled I0, namely e^{-\|x\|} I0(x). |
| `func BesselI1(x float64) float64` | BesselI1 returns the modified Bessel function of the first kind of order one, I1(x). |
| `func BesselI1Prime(x float64) float64` | BesselI1Prime returns the derivative of I1(x). |
| `func BesselI1e(x float64) float64` | BesselI1e returns the exponentially scaled I1, namely e^{-\|x\|} I1(x). |
| `func BesselIn(n int, x float64) float64` | BesselIn returns the modified Bessel function of the first kind of integer order n, In(x). |
| `func BesselInPrime(n int, x float64) float64` | BesselInPrime returns the derivative of the modified Bessel function In(x). |
| `func BesselIne(n int, x float64) float64` | BesselIne returns the exponentially scaled In, namely e^{-\|x\|} In(x). |
| `func BesselInu(nu, x float64) float64` | BesselInu returns the modified Bessel function of the first kind of arbitrary real order nu, I_nu(x), for x >= 0. |
| `func BesselJ(n int, x float64) float64` | BesselJ returns the Bessel function of the first kind of integer order n, Jₙ(x). |
| `func BesselJ0(x float64) float64` | BesselJ0 returns the Bessel function of the first kind of order zero, J0(x). |
| `func BesselJ0Prime(x float64) float64` | BesselJ0Prime returns the derivative of J0, which equals -J1(x). |
| `func BesselJ1(x float64) float64` | BesselJ1 returns the Bessel function of the first kind of order one, J1(x). |
| `func BesselJ1Prime(x float64) float64` | BesselJ1Prime returns the derivative of J1(x). |
| `func BesselJZero(n, s int) float64` | BesselJZero returns the s-th positive zero (s >= 1) of the Bessel function J_n, accurate to full double precision. |
| `func BesselJn(n int, x float64) float64` | BesselJn returns the Bessel function of the first kind of integer order n, Jn(x). |
| `func BesselJnPrime(n int, x float64) float64` | BesselJnPrime returns the derivative of the Bessel function Jn(x). |
| `func BesselJnu(nu, x float64) float64` | BesselJnu returns the Bessel function of the first kind of arbitrary real order nu, J_nu(x), for x >= 0. |
| `func BesselK(n int, x float64) float64` | BesselK returns the modified Bessel function of the second kind of integer order n, Kₙ(x), for x > 0. |
| `func BesselK0(x float64) float64` | BesselK0 returns the modified Bessel function of the second kind of order zero, K0(x). |
| `func BesselK0Prime(x float64) float64` | BesselK0Prime returns the derivative of K0, which equals -K1(x). |
| `func BesselK0e(x float64) float64` | BesselK0e returns the exponentially scaled K0, namely e^{x} K0(x). |
| `func BesselK1(x float64) float64` | BesselK1 returns the modified Bessel function of the second kind of order one, K1(x). |
| `func BesselK1Prime(x float64) float64` | BesselK1Prime returns the derivative of K1(x). |
| `func BesselK1e(x float64) float64` | BesselK1e returns the exponentially scaled K1, namely e^{x} K1(x). |
| `func BesselKn(n int, x float64) float64` | BesselKn returns the modified Bessel function of the second kind of integer order n, Kn(x). |
| `func BesselKnPrime(n int, x float64) float64` | BesselKnPrime returns the derivative of the modified Bessel function Kn(x). |
| `func BesselKne(n int, x float64) float64` | BesselKne returns the exponentially scaled Kn, namely e^{x} Kn(x). |
| `func BesselWronskianIK(x float64) float64` | BesselWronskianIK returns the Wronskian I_n(x) K_n'(x) - I_n'(x) K_n(x), which equals -1/x for every order n. |
| `func BesselWronskianJY(x float64) float64` | BesselWronskianJY returns the Wronskian J_n(x) Y_n'(x) - J_n'(x) Y_n(x), which equals 2/(pi x) for every order n. |
| `func BesselY(n int, x float64) float64` | BesselY returns the Bessel function of the second kind of integer order n, Yₙ(x), for x > 0. |
| `func BesselY0(x float64) float64` | BesselY0 returns the Bessel function of the second kind of order zero, Y0(x). |
| `func BesselY0Prime(x float64) float64` | BesselY0Prime returns the derivative of Y0, which equals -Y1(x). |
| `func BesselY1(x float64) float64` | BesselY1 returns the Bessel function of the second kind of order one, Y1(x). |
| `func BesselY1Prime(x float64) float64` | BesselY1Prime returns the derivative of Y1(x). |
| `func BesselYZero(n, s int) float64` | BesselYZero returns the s-th positive zero (s >= 1) of the Bessel function Y_n, accurate to full double precision. |
| `func BesselYn(n int, x float64) float64` | BesselYn returns the Bessel function of the second kind of integer order n, Yn(x). |
| `func BesselYnPrime(n int, x float64) float64` | BesselYnPrime returns the derivative of the Bessel function Yn(x). |
| `func Beta(a, b float64) float64` | Beta returns the beta function B(a, b) = Γ(a)Γ(b)/Γ(a+b) for a, b > 0. |
| `func BetaInc(a, b, x float64) float64` | BetaInc returns the regularized incomplete beta function, |
| `func Chi(x float64) float64` | Chi returns the hyperbolic cosine integral, |
| `func Ci(x float64) float64` | Ci returns the cosine integral, |
| `func Dawson(x float64) float64` | Dawson returns Dawson's integral, |
| `func Digamma(x float64) float64` | Digamma returns the digamma function ψ(x) = Γ'(x)/Γ(x), the logarithmic derivative of the gamma function, for real x that is not a non-positive… |
| `func DirichletBeta(s float64) float64` | DirichletBeta returns the Dirichlet beta function, |
| `func DirichletEta(s float64) float64` | DirichletEta returns the Dirichlet eta function (alternating zeta), |
| `func E1(x float64) float64` | E1 returns the exponential integral E₁, |
| `func Ei(x float64) float64` | Ei returns the exponential integral, |
| `func EllipticE(m float64) float64` | EllipticE returns the complete elliptic integral of the second kind, |
| `func EllipticEInc(phi, m float64) float64` | EllipticEInc returns the incomplete elliptic integral of the second kind, |
| `func EllipticF(phi, m float64) float64` | EllipticF returns the incomplete elliptic integral of the first kind, |
| `func EllipticK(m float64) float64` | EllipticK returns the complete elliptic integral of the first kind, |
| `func EllipticPi(n, m float64) float64` | EllipticPi returns the complete elliptic integral of the third kind, |
| `func EllipticPiInc(n, phi, m float64) float64` | EllipticPiInc returns the incomplete elliptic integral of the third kind, |
| `func En(n int, x float64) float64` | En returns the generalized exponential integral of order n, |
| `func Erf(x float64) float64` | Erf returns the error function of x, |
| `func Erfc(x float64) float64` | Erfc returns the complementary error function of x, erfc(x) = 1 - erf(x). |
| `func Erfcx(x float64) float64` | Erfcx returns the scaled complementary error function, |
| `func Erfi(x float64) float64` | Erfi returns the imaginary error function, |
| `func Fresnel(x float64) (s, c float64)` | Fresnel returns both Fresnel integrals S(x) and C(x) simultaneously, |
| `func FresnelC(x float64) float64` | FresnelC returns the Fresnel cosine integral, |
| `func FresnelS(x float64) float64` | FresnelS returns the Fresnel sine integral, |
| `func Gamma(x float64) float64` | Gamma returns the gamma function Γ(x); it delegates to the standard library and is provided for API completeness alongside the incomplete forms. |
| `func GammaP(a, x float64) float64` | GammaP returns the regularized lower incomplete gamma function, |
| `func GammaQ(a, x float64) float64` | GammaQ returns the regularized upper incomplete gamma function, |
| `func HurwitzZeta(s, a float64) float64` | HurwitzZeta returns the Hurwitz zeta function, |
| `func HypSinCosIntegral(x float64) (shi, chi float64)` | HypSinCosIntegral returns the hyperbolic sine integral Shi(x) and hyperbolic cosine integral Chi(x) simultaneously. |
| `func IncompleteBeta(a, b, x float64) float64` | IncompleteBeta returns the (non-regularized) incomplete beta function B(x; a, b) = ∫₀ˣ t^{a-1}(1-t)^{b-1} dt = I_x(a, b)·B(a, b). |
| `func InverseErf(x float64) float64` | InverseErf returns the inverse error function, the value y such that erf(y) = x for x in the open interval (-1, 1). |
| `func InverseErfc(x float64) float64` | InverseErfc returns the inverse complementary error function, the value y such that erfc(y) = x for x in the open interval (0, 2). |
| `func KelvinBei(x float64) float64` | KelvinBei returns the Kelvin function bei(x) = Im[J0(x e^{3iπ/4})]. |
| `func KelvinBer(x float64) float64` | KelvinBer returns the Kelvin function ber(x) = Re[J0(x e^{3iπ/4})]. |
| `func KelvinKei(x float64) float64` | KelvinKei returns the Kelvin function kei(x) = Im[K0(x e^{iπ/4})]. |
| `func KelvinKer(x float64) float64` | KelvinKer returns the Kelvin function ker(x) = Re[K0(x e^{iπ/4})]. |
| `func LambertW(x float64) float64` | LambertW returns the principal branch W₀(x) of the Lambert W function, the solution w of w·e^w = x for x ≥ -1/e. |
| `func LambertWm1(x float64) float64` | LambertWm1 returns the secondary real branch W₋₁(x) of the Lambert W function for -1/e ≤ x < 0, the solution w ≤ -1 of w·e^w = x. |
| `func Li(x float64) float64` | Li returns the logarithmic integral, |
| `func Li2(x float64) float64` | Li2 returns the dilogarithm, |
| `func Li3(x float64) float64` | Li3 returns the trilogarithm, |
| `func LogGamma(x float64) float64` | LogGamma returns the natural logarithm of the absolute value of the gamma function, ln\|Γ(x)\|. |
| `func LowerIncompleteGamma(a, x float64) float64` | LowerIncompleteGamma returns the (non-regularized) lower incomplete gamma function γ(a, x) = ∫₀ˣ t^{a-1} e^{-t} dt = P(a, x)·Γ(a). |
| `func ModifiedSphericalI0(x float64) float64` | ModifiedSphericalI0 returns the modified spherical Bessel function of the first kind of order zero, i0(x) = sinh(x)/x. |
| `func ModifiedSphericalI1(x float64) float64` | ModifiedSphericalI1 returns the modified spherical Bessel function of the first kind of order one, i1(x) = (x cosh(x) - sinh(x)) / x^2. |
| `func ModifiedSphericalK0(x float64) float64` | ModifiedSphericalK0 returns the modified spherical Bessel function of the second kind of order zero, k0(x) = (pi/2) e^{-x}/x. |
| `func ModifiedSphericalK1(x float64) float64` | ModifiedSphericalK1 returns the modified spherical Bessel function of the second kind of order one, k1(x) = (pi/2) e^{-x} (1 + 1/x) / x. |
| `func Polygamma(n int, x float64) float64` | Polygamma returns the polygamma function ψ⁽ⁿ⁾(x), the n-th derivative of the digamma function, for integer order n ≥ 0 and real x > 0. |
| `func Polylog(s int, x float64) float64` | Polylog returns the polylogarithm Liₛ(x) = Σ_{k≥1} x^k / k^s for integer order s and real \|x\| < 1 (and x = 1 when s > 1, giving ζ(s)). |
| `func ReciprocalGamma(x float64) float64` | ReciprocalGamma returns 1/Γ(x), which is entire and vanishes at the non-positive integers. |
| `func RiccatiBesselChi(n int, x float64) float64` | RiccatiBesselChi returns the Riccati-Bessel function chi_n(x) = -x yn(x). |
| `func RiccatiBesselPsi(n int, x float64) float64` | RiccatiBesselPsi returns the Riccati-Bessel function psi_n(x) = x jn(x). |
| `func Shi(x float64) float64` | Shi returns the hyperbolic sine integral, |
| `func Si(x float64) float64` | Si returns the sine integral, |
| `func SinCosIntegral(x float64) (si, ci float64)` | SinCosIntegral returns the sine integral Si(x) and cosine integral Ci(x) simultaneously. |
| `func Sinc(x float64) float64` | Sinc returns the unnormalized cardinal sine, sinc(x) = sin(x)/x, with the removable singularity Sinc(0) = 1 handled exactly. |
| `func SincNorm(x float64) float64` | SincNorm returns the normalized cardinal sine, sinc(x) = sin(πx)/(πx), with SincNorm(0) = 1. |
| `func SphericalJ0(x float64) float64` | SphericalJ0 returns the spherical Bessel function of the first kind of order zero, j0(x) = sin(x)/x, with the removable singularity handled at x = 0. |
| `func SphericalJ1(x float64) float64` | SphericalJ1 returns the spherical Bessel function of the first kind of order one, j1(x) = sin(x)/x^2 - cos(x)/x. |
| `func SphericalJn(n int, x float64) float64` | SphericalJn returns the spherical Bessel function of the first kind of order n, jn(x) = sqrt(pi/(2x)) J_{n+1/2}(x). |
| `func SphericalJnPrime(n int, x float64) float64` | SphericalJnPrime returns the derivative of the spherical Bessel function jn(x), using jn'(x) = j_{n-1}(x) - (n+1)/x * jn(x). |
| `func SphericalY0(x float64) float64` | SphericalY0 returns the spherical Bessel function of the second kind of order zero, y0(x) = -cos(x)/x. |
| `func SphericalY1(x float64) float64` | SphericalY1 returns the spherical Bessel function of the second kind of order one, y1(x) = -cos(x)/x^2 - sin(x)/x. |
| `func SphericalYn(n int, x float64) float64` | SphericalYn returns the spherical Bessel function of the second kind of order n, yn(x). |
| `func SphericalYnPrime(n int, x float64) float64` | SphericalYnPrime returns the derivative of the spherical Bessel function yn(x), using yn'(x) = y_{n-1}(x) - (n+1)/x * yn(x). |
| `func Struve(n int, x float64) float64` | Struve returns the Struve function Hₙ(x) of integer order n. |
| `func StruveH0(x float64) float64` | StruveH0 returns the Struve function of order zero, H0(x). |
| `func StruveH1(x float64) float64` | StruveH1 returns the Struve function of order one, H1(x). |
| `func StruveL0(x float64) float64` | StruveL0 returns the modified Struve function of order zero, L0(x). |
| `func StruveL1(x float64) float64` | StruveL1 returns the modified Struve function of order one, L1(x). |
| `func Trigamma(x float64) float64` | Trigamma returns the trigamma function ψ⁽¹⁾(x), the first derivative of the digamma function, for real x that is not a non-positive integer. |
| `func UpperIncompleteGamma(a, x float64) float64` | UpperIncompleteGamma returns the (non-regularized) upper incomplete gamma function Γ(a, x) = ∫ₓ^∞ t^{a-1} e^{-t} dt = Q(a, x)·Γ(a). |
| `func Zeta(s float64) float64` | Zeta returns the Riemann zeta function ζ(s) for real s ≠ 1. |

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/special).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
