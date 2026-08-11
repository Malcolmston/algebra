# approxtheory — implements the core tools of classical approximation theory using nothing but the Go

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/approxtheory.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/approxtheory)

Package approxtheory implements the core tools of classical approximation
theory using nothing but the Go standard library.

- Chebyshev series: computing interpolation coefficients from function samples
(a discrete cosine transform), evaluating a series with the numerically stable
Clenshaw recurrence, differentiating and integrating in coefficient space,
finding the real roots on an interval and doing ordinary arithmetic on series.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/approxtheory"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func BarycentricWeights(xs []float64) []float64` | BarycentricWeights computes the barycentric weights w_j = 1 / prod_{k!=j} (x_j - x_k) for an arbitrary set of distinct nodes. |
| `func BernsteinBasis(k, n int, t float64) float64` | BernsteinBasis evaluates the k-th Bernstein basis polynomial of degree n, b_{k,n}(t) = C(n,k) t^k (1-t)^{n-k}, at t in [0, 1]. |
| `func BernsteinBasisInterval(k, n int, x, a, b float64) float64` | BernsteinBasisInterval evaluates the k-th Bernstein basis polynomial of degree n on the interval [a, b] at x, by mapping x to [0, 1]. |
| `func BernsteinControlPoints(f func(float64) float64, n int, a, b float64) []float64` | BernsteinControlPoints returns the control values f(a + k/n (b-a)) used by the Bernstein approximation operator. |
| `func BernsteinFromMonomial(mono []float64, n int) []float64` | BernsteinFromMonomial returns the Bernstein control values of degree n for a monomial polynomial defined on [0, 1]. |
| `func BezierEval(control []float64, t float64) float64` | BezierEval evaluates a Bezier curve (Bernstein polynomial on [0,1]) with the given control values at parameter t. |
| `func Binomial(n, k int) float64` | Binomial returns the binomial coefficient C(n, k) as a float64, computed multiplicatively to avoid overflow for moderate arguments. |
| `func ChebClenshaw(c []float64, t float64) float64` | ChebClenshaw evaluates the Chebyshev series with coefficients c (index k multiplies T_k) at t in [-1, 1] using Clenshaw's recurrence. |
| `func ChebCoeffsFromValues(values []float64) []float64` | ChebCoeffsFromValues computes the Chebyshev interpolation coefficients from samples of a function taken at the n+1 Chebyshev-Gauss-Lobatto points in… |
| `func ChebGaussPoints(n int, a, b float64) []float64` | ChebGaussPoints returns the n Chebyshev-Gauss points (roots of T_n) on [a, b], ordered from a to b. |
| `func ChebPoints(n int, a, b float64) []float64` | ChebPoints returns the n+1 Chebyshev-Gauss-Lobatto points on [a, b], ordered from a to b. |
| `func ChebT(n int, x float64) float64` | ChebT evaluates the Chebyshev polynomial of the first kind T_n at x in [-1, 1] using the stable three-term recurrence. |
| `func ChebTValues(n int, x float64) []float64` | ChebTValues returns the vector [T_0(x), T_1(x), ..., T_n(x)], or nil for negative n. |
| `func ChebU(n int, x float64) float64` | ChebU evaluates the Chebyshev polynomial of the second kind U_n at x. |
| `func ChebValuesFromCoeffs(coeffs []float64) []float64` | ChebValuesFromCoeffs is the inverse of ChebCoeffsFromValues: it evaluates the Chebyshev series at the n+1 Chebyshev-Gauss-Lobatto nodes (ordered a..b… |
| `func ChebyshevBarycentricWeights(n int) []float64` | ChebyshevBarycentricWeights returns the closed-form barycentric weights for the n+1 Chebyshev-Gauss-Lobatto nodes (ordered a..b as in ChebPoints):… |
| `func ChebyshevNodes(n int, a, b float64) []float64` | ChebyshevNodes returns the n+1 Chebyshev-Gauss-Lobatto nodes on [a, b], ascending. |
| `func ChebyshevQuadrature(f func(float64) float64, n int, a, b float64) float64` | ChebyshevQuadrature integrates f over the domain of a Chebyshev fit by first forming the degree-n interpolant and integrating it exactly. |
| `func ChebyshevToMonomial(cheb []float64) []float64` | ChebyshevToMonomial converts Chebyshev coefficients on [-1, 1] into the coefficients of the equivalent polynomial in the monomial basis (ascending… |
| `func ChebyshevTruncationError(s *ChebSeries, k int) float64` | ChebyshevTruncationError estimates the uniform error of truncating a Chebyshev series after index k as the sum of the absolute values of the… |
| `func ClenshawCurtisQuadrature(f func(float64) float64, n int, a, b float64) float64` | ClenshawCurtisQuadrature approximates the integral of f over [a, b] using n+1 point Clenshaw-Curtis quadrature at the Chebyshev-Gauss-Lobatto nodes. |
| `func ClenshawCurtisWeights(n int, a, b float64) []float64` | ClenshawCurtisWeights returns the n+1 Clenshaw-Curtis quadrature weights on [a, b] associated with the Chebyshev-Gauss-Lobatto nodes (ordered a..b as… |
| `func ContinuedFractionEval(a, b []float64) float64` | ContinuedFractionEval evaluates the continued fraction |
| `func ConvergenceOrder(h1, e1, h2, e2 float64) float64` | ConvergenceOrder estimates the observed order of convergence from two step sizes and their associated errors: log(e1/e2) / log(h1/h2). |
| `func DeCasteljau(control []float64, t float64) float64` | DeCasteljau evaluates a Bezier curve with the given control values at parameter t in [0, 1] using the de Casteljau algorithm. |
| `func DiscreteMinimaxPoly(xs, ys []float64, n, maxIter int) ([]float64, float64, error)` | DiscreteMinimaxPoly computes a degree-n polynomial minimising the maximum absolute residual over the discrete data set (xs, ys) using a single-point… |
| `func EconomizeChebyshev(mono []float64, a, b, tol float64) ([]float64, float64)` | EconomizeChebyshev reduces the degree of a monomial polynomial on [a, b] by converting to the Chebyshev basis, discarding the highest-order… |
| `func EquispacedNodes(n int, a, b float64) []float64` | EquispacedNodes returns n+1 equispaced nodes spanning [a, b]. |
| `func Factorial(k int) float64` | Factorial returns k! as a float64. |
| `func FejerNodes(n int, a, b float64) []float64` | FejerNodes returns the n Fejer nodes (Chebyshev points of the first kind, the roots of T_n) on [a, b], ascending. |
| `func FejerQuadrature(f func(float64) float64, n int, a, b float64) float64` | FejerQuadrature approximates the integral of f over [a, b] using Fejer's first rule with n nodes. |
| `func FejerWeights(n int, a, b float64) []float64` | FejerWeights returns the n weights of Fejer's first quadrature rule on [a, b], associated with the roots of T_n. |
| `func GaussChebyshev2Nodes(n int, a, b float64) []float64` | GaussChebyshev2Nodes returns the n Gauss-Chebyshev nodes of the second kind, x_k = cos(k*pi/(n+1)), mapped to [a, b] in ascending order. |
| `func GaussChebyshev2Quadrature(f func(float64) float64, n int, a, b float64) float64` | GaussChebyshev2Quadrature approximates the plain integral of f over [a, b] using n-point Gauss-Chebyshev quadrature of the second kind, folding the… |
| `func GaussChebyshev2Weights(n int, a, b float64) []float64` | GaussChebyshev2Weights returns the n weights of Gauss-Chebyshev quadrature of the second kind on [a, b], w_k = pi/(n+1) sin^2(k*pi/(n+1)) scaled by… |
| `func GaussChebyshevNodes(n int, a, b float64) []float64` | GaussChebyshevNodes returns the n Gauss-Chebyshev nodes (roots of T_n) on [a, b] in ascending order. |
| `func GaussChebyshevQuadrature(f func(float64) float64, n int, a, b float64) float64` | GaussChebyshevQuadrature approximates the integral of f over [a, b] using n-point Gauss-Chebyshev quadrature of the first kind, folding the Chebyshev… |
| `func GaussChebyshevWeights(n int, a, b float64) []float64` | GaussChebyshevWeights returns the n weights of Gauss-Chebyshev quadrature of the first kind on [a, b]. |
| `func Horner(coeffs []float64, x float64) float64` | Horner is an alias for Polyval; it evaluates a monomial polynomial with Horner's method. |
| `func InterpErrorBound(nodes []float64, derivBound, a, b float64) float64` | InterpErrorBound returns the standard bound on the polynomial interpolation error, \|f(x)-p(x)\| <= (M / (n+1)!) * max\|w(x)\|, where n+1 is the… |
| `func L2Error(f, g func(float64) float64, a, b float64) float64` | L2Error returns the discrete root-integral L2 error of f-g over [a, b], estimated by the composite trapezoidal rule on a fine grid. |
| `func LagrangeBasis(xs []float64, j int, x float64) float64` | LagrangeBasis evaluates the j-th Lagrange cardinal polynomial for the nodes xs at x. |
| `func LagrangeEval(xs, ys []float64, x float64) float64` | LagrangeEval evaluates, at x, the Lagrange interpolating polynomial through (xs, ys) directly (O(n^2) per point). |
| `func LebesgueConstant(nodes []float64, a, b float64) float64` | LebesgueConstant returns the Lebesgue constant, the maximum of the Lebesgue function over [a, b], estimated on a fine grid. |
| `func LebesgueConstantChebyshevAsymptotic(n int) float64` | LebesgueConstantChebyshevAsymptotic returns the classical asymptotic estimate (2/pi) log(n+1) + 0.9625 for the Lebesgue constant of n+1 Chebyshev… |
| `func LebesgueConstantEquispacedAsymptotic(n int) float64` | LebesgueConstantEquispacedAsymptotic returns the classical asymptotic growth estimate 2^{n+1} / (e * n * ln n) of the Lebesgue constant for n+1… |
| `func LebesgueFunction(nodes []float64, x float64) float64` | LebesgueFunction evaluates the Lebesgue function lambda(x) = sum_j \|L_j(x)\| for the interpolation nodes at x, where L_j are the Lagrange cardinal… |
| `func Linspace(a, b float64, n int) []float64` | Linspace returns n equally spaced points spanning the closed interval [a, b]. |
| `func MaxAbsError(a, b []float64) float64` | MaxAbsError returns the maximum absolute difference between paired samples. |
| `func MaxError(f, g func(float64) float64, a, b float64) float64` | MaxError returns the maximum of \|f(x) - g(x)\| over [a, b], estimated on a fine grid. |
| `func MaxNodePolynomial(nodes []float64, a, b float64) float64` | MaxNodePolynomial returns the maximum of \|w(x)\| = \|prod_i (x-nodes`i`)\| over [a, b], estimated on a fine grid. |
| `func MeanAbsError(a, b []float64) float64` | MeanAbsError returns the mean absolute difference between paired samples. |
| `func MinimaxError(f func(float64) float64, coeffs []float64, a, b float64) float64` | MinimaxError returns the achieved uniform error of a monomial polynomial p against f on [a, b], sampled on a fine grid. |
| `func MonomialToChebyshev(mono []float64) []float64` | MonomialToChebyshev converts monomial coefficients of a polynomial on [-1, 1] into Chebyshev coefficients. |
| `func NevilleEval(xs, ys []float64, x float64) float64` | NevilleEval interpolates (xs, ys) at x using Neville's algorithm and returns the interpolated value. |
| `func NewtonDividedDifferences(xs, ys []float64) []float64` | NewtonDividedDifferences returns the divided-difference coefficients of the Newton interpolating polynomial through (xs, ys). |
| `func NewtonEval(coeffs, xs []float64, x float64) float64` | NewtonEval evaluates the Newton form with the given divided-difference coefficients and nodes at x using Horner-like nesting. |
| `func NodePolynomial(nodes []float64, x float64) float64` | NodePolynomial evaluates the node polynomial w(x) = prod_i (x - nodes`i`) at x. |
| `func PolyAdd(a, b []float64) []float64` | PolyAdd returns the sum of two monomial polynomials. |
| `func PolyCompose(p, q []float64) []float64` | PolyCompose returns the coefficients of the composition p(q(x)) for two monomial polynomials. |
| `func PolyDegree(coeffs []float64) int` | PolyDegree returns the degree of a monomial polynomial, ignoring trailing (high-order) coefficients that are exactly zero. |
| `func PolyDeriv(coeffs []float64) []float64` | PolyDeriv returns the derivative of a monomial polynomial. |
| `func PolyDivMod(a, b []float64) (quotient, remainder []float64, err error)` | PolyDivMod divides the monomial polynomial a by b and returns the quotient and remainder, so that a = q*b + r with deg(r) < deg(b). |
| `func PolyEqual(a, b []float64, tol float64) bool` | PolyEqual reports whether two monomial polynomials are equal within the absolute tolerance tol, ignoring trailing zero coefficients. |
| `func PolyFromRoots(roots []float64) []float64` | PolyFromRoots returns the monic monomial polynomial whose roots are exactly the supplied values, i.e. |
| `func PolyInt(coeffs []float64, constant float64) []float64` | PolyInt returns an antiderivative of a monomial polynomial whose constant term equals the supplied value. |
| `func PolyIsZero(coeffs []float64) bool` | PolyIsZero reports whether every coefficient is exactly zero. |
| `func PolyLeadingCoeff(coeffs []float64) float64` | PolyLeadingCoeff returns the coefficient of the highest-order nonzero term, or 0 for the zero polynomial. |
| `func PolyMul(a, b []float64) []float64` | PolyMul returns the product of two monomial polynomials by convolution. |
| `func PolyNormalize(coeffs []float64) []float64` | PolyNormalize returns a copy of coeffs with trailing zero coefficients removed; the zero polynomial, and an empty input, are returned as a single… |
| `func PolyRoots(coeffs []float64, a, b float64) []float64` | PolyRoots returns the real roots of a monomial polynomial that lie in [a, b], sorted ascending. |
| `func PolyScale(a []float64, s float64) []float64` | PolyScale returns the polynomial a scaled by the constant s. |
| `func PolyShift(p []float64, h float64) []float64` | PolyShift returns the coefficients of p(x + h). |
| `func PolyString(coeffs []float64) string` | PolyString renders a monomial polynomial as a human readable expression in the variable x, e.g. |
| `func PolySub(a, b []float64) []float64` | PolySub returns the difference a-b of two monomial polynomials. |
| `func Polyval(coeffs []float64, x float64) float64` | Polyval evaluates the polynomial with the given monomial coefficients (ascending order, coeffs`i` multiplies x**i) at x using Horner's method. |
| `func PolyvalSlice(coeffs []float64, xs []float64) []float64` | PolyvalSlice evaluates a monomial polynomial at every point in xs and returns the resulting slice. |
| `func RMSError(a, b []float64) float64` | RMSError returns the root-mean-square difference between paired samples. |
| `func RSquared(coeffs, xs, ys []float64) float64` | RSquared returns the coefficient of determination for the fit of a monomial polynomial to the data (xs, ys). |
| `func RelativeError(approx, exact, eps float64) float64` | RelativeError returns \|approx-exact\| / max(\|exact\|, eps) guarding against division by zero with the small floor eps. |
| `func RungeFunction(x float64) float64` | RungeFunction evaluates Runge's function 1/(1+25 x^2), the classic example where high-degree equispaced interpolation diverges. |
| `func TaylorAtan(n int) []float64` | TaylorAtan returns the first n+1 Taylor coefficients of atan about 0. |
| `func TaylorCos(n int) []float64` | TaylorCos returns the first n+1 Taylor coefficients of cos about 0. |
| `func TaylorCosh(n int) []float64` | TaylorCosh returns the first n+1 Taylor coefficients of cosh about 0. |
| `func TaylorExp(n int) []float64` | TaylorExp returns the first n+1 Taylor coefficients of exp about 0. |
| `func TaylorExpm1(n int) []float64` | TaylorExpm1 returns the first n+1 Taylor coefficients of exp(x)-1 about 0. |
| `func TaylorGeometric(n int) []float64` | TaylorGeometric returns the first n+1 Taylor coefficients of 1/(1-x), all equal to one. |
| `func TaylorLog1p(n int) []float64` | TaylorLog1p returns the first n+1 Taylor coefficients of log(1+x) about 0. |
| `func TaylorSin(n int) []float64` | TaylorSin returns the first n+1 Taylor coefficients of sin about 0. |
| `func TaylorSinh(n int) []float64` | TaylorSinh returns the first n+1 Taylor coefficients of sinh about 0. |
| `func ThieleInterp(xs, ys []float64, x float64) (float64, error)` | ThieleInterp is a convenience wrapper that builds a Thiele interpolant from (xs, ys) and evaluates it at a single point x. |

### Types

| Type | What it is |
| --- | --- |
| `Barycentric` | Barycentric holds a polynomial interpolant in the second (true) barycentric form, storing the nodes, the values at those nodes and the barycentric… |
| `BernsteinPoly` | BernsteinPoly represents a polynomial in the Bernstein (Bezier) form with control values Coeffs on the interval [A, B]. |
| `ChebSeries` | ChebSeries represents a finite Chebyshev series |
| `PadeResult` | PadeResult holds a Pade approximant P(x)/Q(x) with the numerator and denominator given as monomial coefficients in ascending order. |
| `PolyFitResult` | PolyFitResult holds the result of a polynomial fit: the monomial coefficients (ascending), the root-mean-square residual, the maximum absolute… |
| `RemezResult` | RemezResult holds the outcome of a Remez minimax approximation: the monomial coefficients of the best polynomial (ascending order), the… |
| `Thiele` | Thiele holds a rational interpolant built from Thiele's continued fraction. |

<details>
<summary><code>Barycentric</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func ChebyshevInterpolant(f func(float64) float64, n int, a, b float64) *Barycentric` | ChebyshevInterpolant builds a barycentric polynomial interpolant of f at the n+1 Chebyshev-Gauss-Lobatto nodes on [a, b], using the stable… |
| `func NewBarycentric(xs, ys []float64) (*Barycentric, error)` | NewBarycentric builds a barycentric interpolant through (xs, ys) with weights derived from the nodes. |
| `func NewBarycentricWeighted(xs, ys, weights []float64) (*Barycentric, error)` | NewBarycentricWeighted builds a barycentric interpolant with caller-supplied weights (for example the closed-form Chebyshev weights). |
| `func (b *Barycentric) Eval(x float64) float64` | Eval evaluates the barycentric interpolant at x. |
| `func (b *Barycentric) EvalSlice(xs []float64) []float64` | EvalSlice evaluates the interpolant at every point in xs. |

</details>

<details>
<summary><code>BernsteinPoly</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func BernsteinApprox(f func(float64) float64, n int, a, b float64) *BernsteinPoly` | BernsteinApprox returns the degree-n Bernstein polynomial approximation of f on [a, b], whose control values are the samples f(a + k/n (b-a)). |
| `func NewBernsteinPoly(coeffs []float64, a, b float64) *BernsteinPoly` | NewBernsteinPoly builds a Bernstein-form polynomial with the given control values on [a, b]. |
| `func (p *BernsteinPoly) Clone() *BernsteinPoly` | Clone returns a deep copy of the Bernstein polynomial. |
| `func (p *BernsteinPoly) Degree() int` | Degree returns the polynomial degree of the Bernstein form. |
| `func (p *BernsteinPoly) Derivative() *BernsteinPoly` | Derivative returns the Bernstein form of the derivative of p, correctly scaled for the domain [A, B]. |
| `func (p *BernsteinPoly) Elevate() *BernsteinPoly` | Elevate returns an equivalent Bernstein polynomial of degree n+1 (degree elevation), leaving the represented polynomial unchanged. |
| `func (p *BernsteinPoly) Eval(x float64) float64` | Eval evaluates the Bernstein polynomial at x in [A, B] using the numerically stable de Casteljau algorithm. |
| `func (p *BernsteinPoly) EvalSlice(xs []float64) []float64` | EvalSlice evaluates the Bernstein polynomial at every point in xs. |
| `func (p *BernsteinPoly) Integral() float64` | Integral returns the definite integral of the Bernstein polynomial over its whole domain [A, B]. |
| `func (p *BernsteinPoly) ToMonomial() []float64` | ToMonomial converts the Bernstein form (assumed on [0, 1]) into monomial coefficients in ascending order. |

</details>

<details>
<summary><code>ChebSeries</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func ChebFit(f func(float64) float64, n int, a, b float64) *ChebSeries` | ChebFit returns the degree-n Chebyshev interpolant of f on [a, b], sampling f at the n+1 Chebyshev-Gauss-Lobatto points and transforming to… |
| `func ChebFitGauss(f func(float64) float64, n int, a, b float64) *ChebSeries` | ChebFitGauss returns the degree n-1 Chebyshev approximation of f built from samples at the n Chebyshev-Gauss points (roots of T_n). |
| `func ChebFromMonomial(mono []float64, a, b float64) *ChebSeries` | ChebFromMonomial builds a ChebSeries on [a, b] from a monomial polynomial defined on that same interval (the polynomial is first rescaled to [-1, 1]). |
| `func ChebLeastSquares(xs, ys []float64, degree int, a, b float64) (*ChebSeries, error)` | ChebLeastSquares fits a Chebyshev series of the given degree to sampled data (xs, ys) on [a, b] by least squares in the Chebyshev basis. |
| `func ChebProduct(a, b *ChebSeries) *ChebSeries` | ChebProduct returns the Chebyshev series of the product of two series sharing the same domain. |
| `func NewChebSeries(coeffs []float64, a, b float64) *ChebSeries` | NewChebSeries builds a ChebSeries from raw coefficients on the interval [a, b]. |
| `func (s *ChebSeries) Add(o *ChebSeries) *ChebSeries` | Add returns the sum of two series. |
| `func (s *ChebSeries) Antiderivative() *ChebSeries` | Antiderivative returns a Chebyshev series whose derivative is s and whose value at A is zero (an indefinite integral pinned to the left endpoint),… |
| `func (s *ChebSeries) Clone() *ChebSeries` | Clone returns a deep copy of the series. |
| `func (s *ChebSeries) Coefficient(k int) float64` | Coefficient returns the k-th Chebyshev coefficient, or 0 when k is out of range. |
| `func (s *ChebSeries) Coefficients() []float64` | Coefficients returns a copy of the series coefficients. |
| `func (s *ChebSeries) DefiniteIntegral(x0, x1 float64) float64` | DefiniteIntegral returns the integral of the series from x0 to x1. |
| `func (s *ChebSeries) Degree() int` | Degree returns the polynomial degree of the series (len(Coeffs)-1), or 0 when empty. |
| `func (s *ChebSeries) Derivative() *ChebSeries` | Derivative returns the Chebyshev series of the derivative of s, correctly scaled for the domain [A, B]. |
| `func (s *ChebSeries) Domain() (float64, float64)` | Domain returns the interval endpoints [A, B] of the series. |
| `func (s *ChebSeries) Eval(x float64) float64` | Eval returns the value of the series at x in [A, B]. |
| `func (s *ChebSeries) EvalDerivative(x float64) float64` | EvalDerivative evaluates the derivative of the Chebyshev series at x. |
| `func (s *ChebSeries) EvalSlice(xs []float64) []float64` | EvalSlice evaluates the series at every point in xs. |
| `func (s *ChebSeries) Integral() float64` | Integral returns the integral of the series over its whole domain [A, B]. |
| `func (s *ChebSeries) L2Norm() float64` | L2Norm returns the L2 norm of the series over its domain, sqrt(int f^2 dx), computed exactly from the coefficients via the Chebyshev product. |
| `func (s *ChebSeries) Len() int` | Len returns the number of coefficients in the series. |
| `func (s *ChebSeries) MaxAbs() float64` | MaxAbs returns the maximum absolute value of the series over its domain, estimated on a fine grid. |
| `func (s *ChebSeries) MaxCoeff() float64` | MaxCoeff returns the largest absolute coefficient of the series. |
| `func (s *ChebSeries) Mul(o *ChebSeries) *ChebSeries` | Mul returns the product of two Chebyshev series sharing the same domain. |
| `func (s *ChebSeries) Remap(a, b float64) *ChebSeries` | Remap returns a new series with the same coefficients reinterpreted on the interval [a, b]; the underlying polynomial shape in the unit variable is… |
| `func (s *ChebSeries) Roots() []float64` | Roots returns the real roots of the Chebyshev series lying in [A, B], sorted ascending. |
| `func (s *ChebSeries) Scale(c float64) *ChebSeries` | Scale returns the series multiplied by the constant c. |
| `func (s *ChebSeries) Sub(o *ChebSeries) *ChebSeries` | Sub returns the difference s-o of two series sharing the same domain. |
| `func (s *ChebSeries) TailNorm(k int) float64` | TailNorm returns the sum of absolute values of the coefficients beyond index k, a convenient a-posteriori estimate of the truncation error. |
| `func (s *ChebSeries) ToMonomial() []float64` | ToMonomial converts the Chebyshev series into monomial coefficients (ascending) of the polynomial in the original variable x on [A, B]. |
| `func (s *ChebSeries) Truncate(tol float64) *ChebSeries` | Truncate returns a copy of the series with trailing coefficients whose absolute value is at most tol removed. |

</details>

<details>
<summary><code>PadeResult</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func PadeApprox(taylor []float64, m, n int) (*PadeResult, error)` | PadeApprox computes the m/n Pade approximant of the power series whose Taylor coefficients (about 0) are given in taylor, taylor`k` being the… |
| `func PadeCos(m, n int) *PadeResult` | PadeCos returns the m/n Pade approximant of cos about 0. |
| `func PadeExp(m, n int) *PadeResult` | PadeExp returns the m/n Pade approximant of exp about 0, a convenience wrapper around PadeApprox with the exponential Taylor series. |
| `func PadeLog1p(m, n int) *PadeResult` | PadeLog1p returns the m/n Pade approximant of log(1+x) about 0. |
| `func PadeSin(m, n int) *PadeResult` | PadeSin returns the m/n Pade approximant of sin about 0. |
| `func (p *PadeResult) Eval(x float64) float64` | Eval evaluates the Pade approximant at x. |

</details>

<details>
<summary><code>PolyFitResult</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func PolyFit(xs, ys []float64, degree int) (*PolyFitResult, error)` | PolyFit computes the ordinary least-squares polynomial of the given degree through the data (xs, ys) by solving the normal equations. |
| `func WeightedPolyFit(xs, ys, weights []float64, degree int) (*PolyFitResult, error)` | WeightedPolyFit computes the weighted least-squares polynomial of the given degree. |
| `func (r *PolyFitResult) Eval(x float64) float64` | Eval evaluates the fitted polynomial at x. |

</details>

<details>
<summary><code>RemezResult</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func RemezPoly(f func(float64) float64, n int, a, b float64, maxIter int, tol float64) (*RemezResult, error)` | RemezPoly computes the degree-n minimax (best uniform) polynomial approximation of f on [a, b] using the Remez exchange algorithm. |
| `func (r *RemezResult) Eval(x float64) float64` | Eval evaluates the minimax polynomial at x. |

</details>

<details>
<summary><code>Thiele</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewThiele(xs, ys []float64) (*Thiele, error)` | NewThiele builds Thiele's continued-fraction rational interpolant through the points (xs, ys). |
| `func (t *Thiele) Degree() int` | Degree returns the number of interpolation nodes, one more than the depth of the continued fraction. |
| `func (t *Thiele) Eval(x float64) float64` | Eval evaluates the Thiele continued-fraction interpolant at x. |
| `func (t *Thiele) EvalSlice(xs []float64) []float64` | EvalSlice evaluates the interpolant at every point in xs. |

</details>

### Variables

`ErrDimensionMismatch`, `ErrEmptyInput`, `ErrSingular`, `ErrZeroDivisor`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/approxtheory).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
