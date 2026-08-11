# quadrature — implements numerical integration ("quadrature") in pure Go using only the standard

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/quadrature.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/quadrature)

Package quadrature implements numerical integration ("quadrature") in pure Go
using only the standard library.

- Gaussian quadrature. Nodes and weights for the Gauss rules associated with
the classical weight functions are generated from the three-term recurrence of
the underlying orthogonal polynomials via the Golub-Welsch algorithm (an
eigenvalue problem for a symmetric tridiagonal Jacobi matrix). Legendre,
Chebyshev (four kinds), Hermite (physicists' and probabilists'), Laguerre
(including the generalized form) and Jacobi families are provided, together
with the Gauss-Lobatto and Gauss-Radau variants that fix one or both endpoints
and an embedded Gauss-Kronrod pair for error estimation.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/quadrature"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func AdaptiveLobatto(f Func, a, b, tol float64) float64` | AdaptiveLobatto approximates the integral of f over [a, b] to the requested absolute tolerance by recursively subdividing where a 3-point… |
| `func AdaptiveSimpson(f Func, a, b, tol float64) float64` | AdaptiveSimpson approximates the integral of f over [a, b] to the requested absolute tolerance using recursive adaptive Simpson quadrature with error… |
| `func AdaptiveTrapezoid(f Func, a, b, tol float64) float64` | AdaptiveTrapezoid approximates the integral of f over [a, b] to the requested absolute tolerance using recursive adaptive trapezoidal refinement. |
| `func Average(f Func, a, b, tol float64) float64` | Average returns the mean value of f over [a, b], namely the integral divided by the interval length, computed with adaptive Simpson quadrature. |
| `func Boole(f Func, a, b float64) float64` | Boole approximates the integral of f over [a, b] with Boole's rule using five equally spaced points, exact for polynomials up to degree five. |
| `func BooleSamples(ys []float64, h float64) float64` | BooleSamples integrates uniformly spaced samples ys with spacing h using the composite Boole rule. |
| `func ChebyshevT(n int, x float64) float64` | ChebyshevT evaluates the Chebyshev polynomial of the first kind of degree n at x using the recurrence T_{k+1} = 2x T_k - T_{k-1}. |
| `func ChebyshevU(n int, x float64) float64` | ChebyshevU evaluates the Chebyshev polynomial of the second kind of degree n at x using the recurrence U_{k+1} = 2x U_k - U_{k-1}. |
| `func ClenshawCurtis(f Func, a, b float64, n int) float64` | ClenshawCurtis approximates the integral of f over [a, b] with the (n+1)-point Clenshaw-Curtis rule. |
| `func ClenshawCurtisWeights(n int) (nodes, weights []float64)` | ClenshawCurtisWeights returns the n+1 Clenshaw-Curtis weights on [-1, 1] for the nodes cos(j*pi/n), j = 0..n (Trefethen's algorithm). |
| `func CompositeBoole(f Func, a, b float64, n int) float64` | CompositeBoole approximates the integral of f over [a, b] using the composite Boole rule. |
| `func CompositeGaussLegendre(f Func, a, b float64, n, panels int) float64` | CompositeGaussLegendre approximates the integral of f over [a, b] by splitting the interval into panels equal subintervals and applying the n-point… |
| `func CompositeMidpoint(f Func, a, b float64, n int) float64` | CompositeMidpoint approximates the integral of f over [a, b] using the midpoint rule on n equal subintervals. |
| `func CompositeSimpson(f Func, a, b float64, n int) float64` | CompositeSimpson approximates the integral of f over [a, b] using the composite Simpson's 1/3 rule. |
| `func CompositeSimpson38(f Func, a, b float64, n int) float64` | CompositeSimpson38 approximates the integral of f over [a, b] using the composite Simpson's 3/8 rule. |
| `func CompositeTrapezoid(f Func, a, b float64, n int) float64` | CompositeTrapezoid approximates the integral of f over [a, b] using the trapezoidal rule on n equal subintervals. |
| `func CumulativeTrapezoid(xs, ys []float64) []float64` | CumulativeTrapezoid returns the running integral of samples ys at abscissae xs computed with the trapezoidal rule. |
| `func DoubleAverage(f Func2, ax, bx, ay, by float64, n int) float64` | DoubleAverage returns the mean value of f over the rectangle [ax, bx] x [ay, by], computed with the n-point tensor Gauss-Legendre rule. |
| `func DoubleBoole(f Func2, ax, bx, ay, by float64, nx, ny int) float64` | DoubleBoole approximates the double integral of f over the rectangle using the composite Boole rule with nx and ny subintervals. |
| `func DoubleExponential(f Func, a, b, tol float64) float64` | DoubleExponential is an alias for TanhSinh, the name by which the method is also known. |
| `func DoubleGaussLegendre(f Func2, ax, bx, ay, by float64, n int) float64` | DoubleGaussLegendre approximates the double integral of f over the rectangle [ax, bx] x [ay, by] using the n-point Gauss-Legendre rule in each… |
| `func DoubleMidpoint(f Func2, ax, bx, ay, by float64, nx, ny int) float64` | DoubleMidpoint approximates the double integral of f over the rectangle using the composite midpoint rule with nx and ny subintervals. |
| `func DoubleSimpson(f Func2, ax, bx, ay, by float64, nx, ny int) float64` | DoubleSimpson approximates the double integral of f over the rectangle using the composite Simpson rule with nx and ny subintervals. |
| `func DoubleTrapezoid(f Func2, ax, bx, ay, by float64, nx, ny int) float64` | DoubleTrapezoid approximates the double integral of f over the rectangle [ax, bx] x [ay, by] using the composite trapezoidal rule with nx and ny… |
| `func ExpSinh(f Func, a, tol float64) float64` | ExpSinh approximates the integral of f over the semi-infinite interval [a, inf) using the exp-sinh double-exponential transformation x = a +… |
| `func ExpSinhLeft(f Func, b, tol float64) float64` | ExpSinhLeft approximates the integral of f over the semi-infinite interval (-inf, b] using the exp-sinh transformation x = b - exp((pi/2) sinh t). |
| `func ExpectationExponential(f Func, n int) float64` | ExpectationExponential approximates E[f(X)] for an exponential random variable X with rate 1 (density e^{-x} on [0, inf)) using the n-point… |
| `func ExpectationGaussHermite(f Func, n int) float64` | ExpectationGaussHermite approximates E[f(X)] for a standard normal random variable X ~ N(0, 1) using the n-point probabilists' Gauss-Hermite rule. |
| `func ExpectationGaussHermiteNormal(f Func, mean, sd float64, n int) float64` | ExpectationGaussHermiteNormal approximates E[f(X)] for a normal random variable X ~ N(mean, sd^2) using the n-point probabilists' Gauss-Hermite rule. |
| `func Fejer1Weights(n int) (nodes, weights []float64)` | Fejer1Weights returns the n nodes and weights of Fejer's first rule on [-1, 1], an open Chebyshev rule with nodes cos((2i-1)pi/(2n)). |
| `func Fejer2Weights(n int) (nodes, weights []float64)` | Fejer2Weights returns the n nodes and weights of Fejer's second rule on [-1, 1], an open Chebyshev rule with nodes cos(i*pi/(n+1)), i = 1..n. |
| `func GaussChebyshev1(n int) (nodes, weights []float64)` | GaussChebyshev1 returns the n nodes and weights of the Gauss-Chebyshev rule of the first kind, for the integral of f(x)/sqrt(1-x^2) over [-1, 1]. |
| `func GaussChebyshev2(n int) (nodes, weights []float64)` | GaussChebyshev2 returns the n nodes and weights of the Gauss-Chebyshev rule of the second kind, for the integral of f(x)*sqrt(1-x^2) over [-1, 1]. |
| `func GaussChebyshev3(n int) (nodes, weights []float64)` | GaussChebyshev3 returns the n nodes and weights of the Gauss-Chebyshev rule of the third kind, for the integral of f(x)*sqrt((1+x)/(1-x)) over [-1,… |
| `func GaussChebyshev4(n int) (nodes, weights []float64)` | GaussChebyshev4 returns the n nodes and weights of the Gauss-Chebyshev rule of the fourth kind, for the integral of f(x)*sqrt((1-x)/(1+x)) over [-1,… |
| `func GaussGegenbauer(n int, lambda float64) (nodes, weights []float64)` | GaussGegenbauer returns the n nodes and weights of the Gauss-Gegenbauer (ultraspherical) rule for the integral of f(x)*(1-x^2)^(lambda-1/2) over [-1,… |
| `func GaussHermite(n int) (nodes, weights []float64)` | GaussHermite returns the n nodes and weights of the Gauss-Hermite rule for the integral of f(x)*e^{-x^2} over the whole real line (physicists'… |
| `func GaussHermiteProb(n int) (nodes, weights []float64)` | GaussHermiteProb returns the n nodes and weights of the Gauss-Hermite rule for the integral of f(x)*e^{-x^2/2} over the whole real line… |
| `func GaussJacobi(n int, a, b float64) (nodes, weights []float64)` | GaussJacobi returns the n nodes and weights of the Gauss-Jacobi rule for the integral of f(x)*(1-x)^a*(1+x)^b over [-1, 1]. |
| `func GaussKronrod15() (nodes, kronrod, gauss []float64)` | GaussKronrod15 returns the abscissae and weights of the 15-point Gauss-Kronrod rule on [-1, 1] together with the weights of the embedded 7-point… |
| `func GaussKronrod15Eval(f Func, a, b float64) (value, absErr float64)` | GaussKronrod15Eval applies the 15-point Gauss-Kronrod rule to f over [a, b], returning the Kronrod estimate together with a QUADPACK-style estimate… |
| `func GaussLaguerre(n int) (nodes, weights []float64)` | GaussLaguerre returns the n nodes and weights of the Gauss-Laguerre rule for the integral of f(x)*e^{-x} over [0, inf). |
| `func GaussLaguerreGen(n int, a float64) (nodes, weights []float64)` | GaussLaguerreGen returns the n nodes and weights of the generalized Gauss-Laguerre rule for the integral of f(x)*x^a*e^{-x} over [0, inf). |
| `func GaussLegendre(n int) (nodes, weights []float64)` | GaussLegendre returns the n nodes and weights of the Gauss-Legendre rule on [-1, 1], which integrates polynomials up to degree 2n-1 exactly. |
| `func GaussLobatto(n int) (nodes, weights []float64)` | GaussLobatto returns the n nodes and weights of the Gauss-Lobatto-Legendre rule on [-1, 1]. |
| `func GaussRadau(n int) (nodes, weights []float64)` | GaussRadau returns the n nodes and weights of the Gauss-Radau-Legendre rule on [-1, 1] that fixes the left endpoint -1 as a node. |
| `func GaussRadauRight(n int) (nodes, weights []float64)` | GaussRadauRight returns the n nodes and weights of the Gauss-Radau-Legendre rule on [-1, 1] that fixes the right endpoint +1 as a node. |
| `func GolubWelsch(alpha, beta []float64) (nodes, weights []float64)` | GolubWelsch computes Gauss quadrature nodes and weights from the monic three-term recurrence coefficients of the associated orthogonal polynomials. |
| `func Halton(i, dim int) []float64` | Halton returns the i-th point (i >= 0) of the Halton low-discrepancy sequence in dim dimensions, using the first dim prime numbers as bases. |
| `func HaltonSequence(n, dim int) [][]float64` | HaltonSequence returns the first n points of the dim-dimensional Halton sequence, skipping the degenerate index 0. |
| `func HermiteH(n int, x float64) float64` | HermiteH evaluates the physicists' Hermite polynomial H_n at x using the recurrence H_{k+1} = 2x H_k - 2k H_{k-1}. |
| `func HermiteHe(n int, x float64) float64` | HermiteHe evaluates the probabilists' Hermite polynomial He_n at x using the recurrence He_{k+1} = x He_k - k He_{k-1}. |
| `func IntegrateBoxGaussLegendre(f FuncN, lower, upper []float64, nPerDim int) float64` | IntegrateBoxGaussLegendre approximates the integral of f over the box given by lower and upper corners using an nPerDim-point tensor Gauss-Legendre… |
| `func IntegrateFejer1(f Func, a, b float64, n int) float64` | IntegrateFejer1 approximates the integral of f over [a, b] with Fejer's first rule of n points. |
| `func IntegrateFejer2(f Func, a, b float64, n int) float64` | IntegrateFejer2 approximates the integral of f over [a, b] with Fejer's second rule of n points. |
| `func IntegrateGaussChebyshev1(f Func, n int) float64` | IntegrateGaussChebyshev1 approximates the integral of f(x)/sqrt(1-x^2) over [-1, 1] with the n-point first-kind Gauss-Chebyshev rule. |
| `func IntegrateGaussChebyshev2(f Func, n int) float64` | IntegrateGaussChebyshev2 approximates the integral of f(x)*sqrt(1-x^2) over [-1, 1] with the n-point second-kind Gauss-Chebyshev rule. |
| `func IntegrateGaussChebyshev3(f Func, n int) float64` | IntegrateGaussChebyshev3 approximates the integral of f(x)*sqrt((1+x)/(1-x)) over [-1, 1] with the n-point third-kind rule. |
| `func IntegrateGaussChebyshev4(f Func, n int) float64` | IntegrateGaussChebyshev4 approximates the integral of f(x)*sqrt((1-x)/(1+x)) over [-1, 1] with the n-point fourth-kind rule. |
| `func IntegrateGaussHermite(f Func, n int) float64` | IntegrateGaussHermite approximates the integral of f(x)*e^{-x^2} over the whole real line with the n-point Gauss-Hermite rule (physicists'… |
| `func IntegrateGaussHermiteFunc(f Func, n int) float64` | IntegrateGaussHermiteFunc approximates the integral of f (without any weight factor) over the whole real line by folding the Gaussian weight back… |
| `func IntegrateGaussHermiteProb(f Func, n int) float64` | IntegrateGaussHermiteProb approximates the integral of f(x)*e^{-x^2/2} over the whole real line with the n-point probabilists' Gauss-Hermite rule. |
| `func IntegrateGaussJacobi(f Func, n int, a, b float64) float64` | IntegrateGaussJacobi approximates the integral of f(x)*(1-x)^a*(1+x)^b over [-1, 1] with the n-point Gauss-Jacobi rule. |
| `func IntegrateGaussKronrod(f Func, a, b float64) float64` | IntegrateGaussKronrod is an alias for GaussKronrod15Eval returning only the integral estimate over [a, b]. |
| `func IntegrateGaussLaguerre(f Func, n int) float64` | IntegrateGaussLaguerre approximates the integral of f(x)*e^{-x} over [0, inf) with the n-point Gauss-Laguerre rule. |
| `func IntegrateGaussLaguerreFunc(f Func, n int) float64` | IntegrateGaussLaguerreFunc approximates the integral of f (without any weight factor) over [0, inf) by folding the e^{-x} weight back in: the summand… |
| `func IntegrateGaussLaguerreGen(f Func, n int, a float64) float64` | IntegrateGaussLaguerreGen approximates the integral of f(x)*x^a*e^{-x} over [0, inf) with the n-point generalized Gauss-Laguerre rule. |
| `func IntegrateGaussLegendre(f Func, a, b float64, n int) float64` | IntegrateGaussLegendre approximates the integral of f over [a, b] with the n-point Gauss-Legendre rule, mapping the canonical nodes onto [a, b]. |
| `func IntegrateGaussLobatto(f Func, a, b float64, n int) float64` | IntegrateGaussLobatto approximates the integral of f over [a, b] with the n-point Gauss-Lobatto rule. |
| `func IntegrateGaussRadau(f Func, a, b float64, n int) float64` | IntegrateGaussRadau approximates the integral of f over [a, b] with the n-point left Gauss-Radau rule. |
| `func IntegrateInfinite(f Func, tol float64) float64` | IntegrateInfinite approximates the integral of f over the entire real line using the substitution x = t/(1-t^2), t in (-1, 1), evaluated with the… |
| `func IntegrateNewtonCotesClosed(f Func, a, b float64, n int) float64` | IntegrateNewtonCotesClosed approximates the integral of f over [a, b] with a single closed (n+1)-point Newton-Cotes rule of order n. |
| `func IntegrateNewtonCotesOpen(f Func, a, b float64, n int) float64` | IntegrateNewtonCotesOpen approximates the integral of f over [a, b] with a single open n-point Newton-Cotes rule. |
| `func IntegrateProduct2(f Func2, rx, ry Rule) float64` | IntegrateProduct2 applies the tensor product of two one-dimensional rules to a bivariate function, evaluating f at every (node_i, node_j) pair… |
| `func IntegrateSamples(ys []float64, h float64) float64` | IntegrateSamples integrates uniformly spaced samples ys with spacing h, selecting the highest-order composite rule the sample count permits: Boole… |
| `func JacobiP(n int, a, b, x float64) float64` | JacobiP evaluates the Jacobi polynomial P_n^{(a,b)} at x using the standard three-term recurrence. |
| `func LaguerreL(n int, x float64) float64` | LaguerreL evaluates the Laguerre polynomial L_n at x using the recurrence (k+1) L_{k+1} = (2k+1-x) L_k - k L_{k-1}. |
| `func LaguerreLGen(n int, a, x float64) float64` | LaguerreLGen evaluates the generalized Laguerre polynomial L_n^{(a)} at x using the recurrence (k+1) L_{k+1} = (2k+1+a-x) L_k - (k+a) L_{k-1}. |
| `func LegendreP(n int, x float64) float64` | LegendreP evaluates the Legendre polynomial of degree n at x using the stable three-term recurrence. |
| `func LegendrePDeriv(n int, x float64) float64` | LegendrePDeriv evaluates the derivative of the Legendre polynomial of degree n at x. |
| `func Midpoint(f Func, a, b float64) float64` | Midpoint approximates the integral of f over [a, b] with the single-panel midpoint rule (b-a)*f((a+b)/2), exact for linear functions. |
| `func MidpointSamples(ys []float64, h float64) float64` | MidpointSamples integrates samples ys, interpreted as the values of f at the midpoints of n equal panels of width h, using the midpoint rule. |
| `func NewtonCotesClosed(n int) []float64` | NewtonCotesClosed returns the closed (n+1)-point Newton-Cotes weights as float64 values (see NewtonCotesClosedExact). |
| `func NewtonCotesClosedExact(n int) []*big.Rat` | NewtonCotesClosedExact returns the exact rational weights of the closed (n+1)-point Newton-Cotes rule on the reference interval [0, n], scaled so… |
| `func NewtonCotesOpen(n int) []float64` | NewtonCotesOpen returns the open n-point Newton-Cotes weights as float64 values (see NewtonCotesOpenExact). |
| `func NewtonCotesOpenExact(n int) []*big.Rat` | NewtonCotesOpenExact returns the exact rational weights of the open n-point Newton-Cotes rule with interior nodes at 1,2,...,n on the reference… |
| `func RecurrenceChebyshev1(n int) (alpha, beta []float64)` | RecurrenceChebyshev1 returns the monic recurrence coefficients for the Chebyshev polynomials of the first kind (weight 1/sqrt(1-x^2) on [-1, 1]). |
| `func RecurrenceChebyshev2(n int) (alpha, beta []float64)` | RecurrenceChebyshev2 returns the monic recurrence coefficients for the Chebyshev polynomials of the second kind (weight sqrt(1-x^2) on [-1, 1]). |
| `func RecurrenceHermite(n int) (alpha, beta []float64)` | RecurrenceHermite returns the monic recurrence coefficients for the physicists' Hermite polynomials (weight e^{-x^2} on the whole real line). |
| `func RecurrenceHermiteProb(n int) (alpha, beta []float64)` | RecurrenceHermiteProb returns the monic recurrence coefficients for the probabilists' Hermite polynomials (weight e^{-x^2/2} on the whole real line). |
| `func RecurrenceJacobi(n int, a, b float64) (alpha, beta []float64)` | RecurrenceJacobi returns the monic recurrence coefficients for the Jacobi polynomials with parameters a and b (weight (1-x)^a * (1+x)^b on [-1, 1]). |
| `func RecurrenceLaguerre(n int, a float64) (alpha, beta []float64)` | RecurrenceLaguerre returns the monic recurrence coefficients for the generalized Laguerre polynomials with parameter a (weight x^a * e^{-x} on [0,… |
| `func RecurrenceLegendre(n int) (alpha, beta []float64)` | RecurrenceLegendre returns the monic three-term recurrence coefficients for the Legendre polynomials (weight 1 on [-1, 1]) for an n-point rule. |
| `func Rescale(nodes, weights []float64, a, b float64) (rn, rw []float64)` | Rescale maps parallel node/weight slices from the canonical interval [-1, 1] to [a, b], returning fresh slices. |
| `func RichardsonExtrapolate(coarse, fine, ratio, order float64) float64` | RichardsonExtrapolate performs one step of Richardson extrapolation on two estimates A(h) and A(h/ratio) of the same quantity whose leading error… |
| `func RombergTable(f Func, a, b float64, k int) [][]float64` | RombergTable computes the full Romberg extrapolation table for the integral of f over [a, b] using k levels of trapezoidal refinement (rows 0..k-1). |
| `func RombergValue(f Func, a, b float64, k int) float64` | RombergValue approximates the integral of f over [a, b] with the k-level Romberg table, returning only the best estimate T`k-1``k-1`. |
| `func Simpson(f Func, a, b float64) float64` | Simpson approximates the integral of f over [a, b] with Simpson's 1/3 rule using the two endpoints and the midpoint, exact for cubics. |
| `func Simpson38(f Func, a, b float64) float64` | Simpson38 approximates the integral of f over [a, b] with Simpson's 3/8 rule using four equally spaced points, exact for cubics. |
| `func SimpsonSamples(ys []float64, h float64) float64` | SimpsonSamples integrates uniformly spaced samples ys with spacing h using the composite Simpson rule. |
| `func SinhSinh(f Func, tol float64) float64` | SinhSinh approximates the integral of f over the entire real line (-inf, inf) using the sinh-sinh double-exponential transformation x = sinh((pi/2)… |
| `func TanhSinh(f Func, a, b, tol float64) float64` | TanhSinh approximates the integral of f over the finite interval [a, b] using the tanh-sinh (double-exponential) rule. |
| `func Trapezoid(f Func, a, b float64) float64` | Trapezoid approximates the integral of f over [a, b] with the single-panel trapezoidal rule (b-a)/2*(f(a)+f(b)), exact for linear functions. |
| `func TrapezoidNonUniform(xs, ys []float64) float64` | TrapezoidNonUniform integrates samples ys taken at the (not necessarily equally spaced) abscissae xs using the trapezoidal rule. |
| `func TrapezoidSamples(ys []float64, h float64) float64` | TrapezoidSamples integrates uniformly spaced samples ys with spacing h using the composite trapezoidal rule. |
| `func TripleGaussLegendre(f Func3, ax, bx, ay, by, az, bz float64, n int) float64` | TripleGaussLegendre approximates the triple integral of f over the box [ax, bx] x [ay, by] x [az, bz] using the n-point Gauss-Legendre rule in each… |
| `func TripleMidpoint(f Func3, ax, bx, ay, by, az, bz float64, nx, ny, nz int) float64` | TripleMidpoint approximates the triple integral of f over the box using the composite midpoint rule with nx, ny and nz subintervals. |
| `func TripleSimpson(f Func3, ax, bx, ay, by, az, bz float64, nx, ny, nz int) float64` | TripleSimpson approximates the triple integral of f over the box using the composite Simpson rule with nx, ny and nz subintervals. |
| `func TripleTrapezoid(f Func3, ax, bx, ay, by, az, bz float64, nx, ny, nz int) float64` | TripleTrapezoid approximates the triple integral of f over the box using the composite trapezoidal rule with nx, ny and nz subintervals. |
| `func VanDerCorput(i, base int) float64` | VanDerCorput returns the i-th element (i >= 0) of the van der Corput low-discrepancy sequence in the given base, a number in [0, 1) formed by… |

### Types

| Type | What it is |
| --- | --- |
| `Box` | Box represents an axis-aligned rectangular region in n dimensions given by its lower and upper corners, which must have equal length. |
| `Func` | Func is a real-valued function of one real variable. |
| `Func2` | Func2 is a real-valued function of two real variables. |
| `Func3` | Func3 is a real-valued function of three real variables. |
| `FuncN` | FuncN is a real-valued function of a point in n-dimensional space. |
| `Interval` | Interval represents the closed one-dimensional interval [A, B]. |
| `Result` | Result bundles a computed integral value with an estimate of its absolute error and the number of function evaluations spent producing it. |
| `Rule` | Rule is a quadrature rule: a set of abscissae (Nodes) with matching Weights such that the integral of f is approximated by the weighted sum of f at… |

<details>
<summary><code>Box</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewBox(lower, upper []float64) Box` | NewBox constructs a Box from lower and upper corner slices, copying them. |
| `func (b Box) Contains(x []float64) bool` | Contains reports whether the point x lies within the closed box. |
| `func (b Box) Dim() int` | Dim reports the dimensionality of the box. |
| `func (b Box) Volume() float64` | Volume returns the product of the side lengths of the box. |

</details>

<details>
<summary><code>Interval</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (iv Interval) Contains(x float64) bool` | Contains reports whether x lies within the closed interval. |
| `func (iv Interval) Length() float64` | Length returns B - A. |
| `func (iv Interval) Mid() float64` | Mid returns the midpoint of the interval. |

</details>

<details>
<summary><code>Result</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func AdaptiveGaussKronrod(f Func, a, b, tol float64) Result` | AdaptiveGaussKronrod approximates the integral of f over [a, b] to the requested absolute tolerance by recursively subdividing the interval, applying… |
| `func AdaptiveSimpsonResult(f Func, a, b, tol float64) Result` | AdaptiveSimpsonResult is like AdaptiveSimpson but returns a Result carrying the error estimate and evaluation count. |
| `func MonteCarlo2(f Func2, ax, bx, ay, by float64, samples int, seed int64) Result` | MonteCarlo2 estimates the double integral of f over [ax, bx] x [ay, by] by plain Monte-Carlo sampling. |
| `func MonteCarlo3(f Func3, ax, bx, ay, by, az, bz float64, samples int, seed int64) Result` | MonteCarlo3 estimates the triple integral of f over the box [ax, bx] x [ay, by] x [az, bz] by plain Monte-Carlo sampling. |
| `func MonteCarloBox(f FuncN, lower, upper []float64, samples int, seed int64) Result` | MonteCarloBox estimates the integral of f over the box given by lower and upper corners using plain Monte-Carlo sampling with the given number of… |
| `func QuasiMonteCarloBox(f FuncN, lower, upper []float64, samples int) Result` | QuasiMonteCarloBox estimates the integral of f over the box using the Halton quasi-random sequence. |
| `func Romberg(f Func, a, b, tol float64, maxLevels int) Result` | Romberg approximates the integral of f over [a, b] using Romberg extrapolation, stopping when two successive diagonal estimates agree to the given… |
| `func StratifiedMonteCarloBox(f FuncN, lower, upper []float64, strata int, seed int64) Result` | StratifiedMonteCarloBox estimates the integral of f over the box using stratified sampling: each dimension is split into strata subintervals, forming… |
| `func TanhSinhResult(f Func, a, b, tol float64) Result` | TanhSinhResult is like TanhSinh but returns a Result with the error estimate. |

</details>

<details>
<summary><code>Rule</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func ClenshawCurtisRule(n int) Rule` | ClenshawCurtisRule returns the (n+1)-point Clenshaw-Curtis rule on [-1, 1]. |
| `func Fejer1Rule(n int) Rule` | Fejer1Rule returns Fejer's first rule of n points on [-1, 1]. |
| `func Fejer2Rule(n int) Rule` | Fejer2Rule returns Fejer's second rule of n points on [-1, 1]. |
| `func GaussChebyshev1Rule(n int) Rule` | GaussChebyshev1Rule returns the first-kind Gauss-Chebyshev rule of n points. |
| `func GaussChebyshev2Rule(n int) Rule` | GaussChebyshev2Rule returns the second-kind Gauss-Chebyshev rule of n points. |
| `func GaussChebyshev3Rule(n int) Rule` | GaussChebyshev3Rule returns the third-kind Gauss-Chebyshev rule of n points. |
| `func GaussChebyshev4Rule(n int) Rule` | GaussChebyshev4Rule returns the fourth-kind Gauss-Chebyshev rule of n points. |
| `func GaussGegenbauerRule(n int, lambda float64) Rule` | GaussGegenbauerRule returns the Gauss-Gegenbauer rule of n points. |
| `func GaussHermiteProbRule(n int) Rule` | GaussHermiteProbRule returns the probabilists' Gauss-Hermite rule of n points. |
| `func GaussHermiteRule(n int) Rule` | GaussHermiteRule returns the physicists' Gauss-Hermite rule of n points. |
| `func GaussJacobiRule(n int, a, b float64) Rule` | GaussJacobiRule returns the Gauss-Jacobi rule of n points with parameters a and b. |
| `func GaussLaguerreGenRule(n int, a float64) Rule` | GaussLaguerreGenRule returns the generalized Gauss-Laguerre rule of n points with parameter a. |
| `func GaussLaguerreRule(n int) Rule` | GaussLaguerreRule returns the Gauss-Laguerre rule of n points. |
| `func GaussLegendreRule(n int) Rule` | GaussLegendreRule returns the Gauss-Legendre rule of n points on [-1, 1] as a Rule value. |
| `func GaussLobattoRule(n int) Rule` | GaussLobattoRule returns the Gauss-Lobatto-Legendre rule of n points. |
| `func GaussRadauRule(n int) Rule` | GaussRadauRule returns the left Gauss-Radau-Legendre rule of n points. |
| `func NewRule(nodes, weights []float64) Rule` | NewRule builds a Rule from parallel node and weight slices. |
| `func NewtonCotesClosedRule(n int) Rule` | NewtonCotesClosedRule returns the closed (n+1)-point Newton-Cotes rule on the canonical interval [-1, 1] with equally spaced nodes. |
| `func NewtonCotesOpenRule(n int) Rule` | NewtonCotesOpenRule returns the open n-point Newton-Cotes rule on the canonical interval [-1, 1] with equally spaced interior nodes. |
| `func (r Rule) Clone() Rule` | Clone returns a deep copy of the rule. |
| `func (r Rule) Integrate(f Func) float64` | Integrate applies the rule to f, returning the weighted sum of f at the nodes. |
| `func (r Rule) IntegrateOn(f Func, a, b float64) float64` | IntegrateOn scales the canonical rule to [a, b] and applies it to f in one step. |
| `func (r Rule) Len() int` | Len reports the number of nodes in the rule. |
| `func (r Rule) Reversed() Rule` | Reversed returns a copy of the rule with nodes and weights in reverse order. |
| `func (r Rule) Scale(a, b float64) Rule` | Scale maps a rule defined on the canonical interval [-1, 1] to the interval [a, b] by the affine change of variables x = (a+b)/2 + (b-a)/2 * t,… |
| `func (r Rule) Transform(g, dg Func) Rule` | Transform returns a new rule obtained from r by the change of variables u = g(x): each node x_i becomes g(x_i) and each weight w_i is multiplied by… |
| `func (r Rule) WeightSum() float64` | WeightSum returns the sum of the rule's weights. |

</details>

### Variables

`ErrNonPositiveN`, `ErrBadParameter`, `ErrDimMismatch`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/quadrature).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
