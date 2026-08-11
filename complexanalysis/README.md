# complexanalysis — provides numerical tools for complex analysis using only the Go standard library

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/complexanalysis.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/complexanalysis)

Package complexanalysis provides numerical tools for complex analysis using
only the Go standard library.

- Elementary and inverse functions on complex128 with well-defined principal
branches (see elementary.go): `Sqrt`, `Cbrt`, `Exp`, `Log`, `Pow`, the
trigonometric and hyperbolic families, and their inverses. - Special functions
(see special.go): the complex `Gamma` and `LogGamma` via the Lanczos
approximation, `Digamma`, `Beta`, the error functions `Erf`/`Erfc`, the
Riemann `Zeta` function via Borwein's algorithm, and the combinatorial
`Factorial`, `RisingFactorial` and `Binomial`. - Numerical contour integration
and its consequences (see contour.go): `IntegrateCircle`, `IntegrateSegment`,
`IntegratePolygon` and the general `ContourIntegral`; `Residue`,
`ResidueSimplePole` and `ResidueOrderM`; the `CauchyIntegralValue` and
`CauchyDerivative` formulas; the `WindingNumber`; and the `ArgumentPrinciple`
with `CountZeros`. - Conformal maps (see conformal.go): the `Mobius`
transformation type with composition, inversion and fixed points, together
with the `CrossRatio`, the `CayleyTransform` and the `JoukowskiMap`. - Laurent
and Taylor coefficient extraction and numeric analytic continuation (see
laurent.go): `LaurentCoefficient`, `TaylorCoefficient`, `PowerSeriesEval` and
`AnalyticContinuation`.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/complexanalysis"
```

## Usage

This is the package's own `ExampleCauchyDerivative`, so it compiles and its output is
asserted on every `go test ./complexanalysis/`.

```go
f := cmplx.Exp
	for k := 0; k <= 3; k++ {
		fmt.Println(k, fmtC(complexanalysis.CauchyDerivative(f, 0, 0, k, 1, 512)))
	}
```

```
0 1.000000+0.000000i
1 1.000000+0.000000i
2 1.000000+0.000000i
3 1.000000+0.000000i
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func Abs(z complex128) float64` | Abs returns the modulus (absolute value) \|z\|. |
| `func Acos(z complex128) complex128` | Acos returns the principal inverse cosine of z. |
| `func Acosh(z complex128) complex128` | Acosh returns the principal inverse hyperbolic cosine of z. |
| `func AnalyticContinuation(f Function, from, to complex128, order int, radius float64, n int) complex128` | AnalyticContinuation returns the value of f at the point to, reconstructed from its Taylor expansion about the point from. |
| `func AnalyticContinuationPath(f Function, path []complex128, order int, radius float64, n int) complex128` | AnalyticContinuationPath continues f numerically along the given path of centers and returns the reconstructed value at the final point. |
| `func ApproxEqual(a, b complex128, tol float64) bool` | ApproxEqual reports whether \|a-b\| <= tol. |
| `func Arg(z complex128) float64` | Arg returns the principal argument of z in the interval (-pi, pi]. |
| `func ArgumentPrinciple(f Function, center complex128, radius float64, n int) complex128` | ArgumentPrinciple returns the value of 1/(2*pi*i) * integral of f'(z)/f(z) around the positively oriented circle of the given center and radius. |
| `func Asin(z complex128) complex128` | Asin returns the principal inverse sine of z. |
| `func Asinh(z complex128) complex128` | Asinh returns the principal inverse hyperbolic sine of z. |
| `func Atan(z complex128) complex128` | Atan returns the principal inverse tangent of z. |
| `func Atanh(z complex128) complex128` | Atanh returns the principal inverse hyperbolic tangent of z. |
| `func Beta(a, b complex128) complex128` | Beta returns the Euler beta function B(a, b) = Gamma(a)Gamma(b)/Gamma(a+b), evaluated through LogGamma to remain stable for moderately large… |
| `func Binomial(z complex128, k int) complex128` | Binomial returns the generalized binomial coefficient C(z, k) = Gamma(z+1)/(Gamma(k+1)Gamma(z-k+1)) for a non-negative integer k, evaluated as a… |
| `func CauchyDerivative(f Function, z0, center complex128, k int, radius float64, n int) complex128` | CauchyDerivative evaluates the k-th derivative of f at z0 using the Cauchy integral formula for derivatives, f^(k)(z0) = k!/(2*pi*i) * integral of… |
| `func CauchyIntegralValue(f Function, z0, center complex128, radius float64, n int) complex128` | CauchyIntegralValue evaluates f at the point z0 inside a circle of the given center and radius using the Cauchy integral formula, f(z0) = 1/(2*pi*i)… |
| `func CayleyTransform(z complex128) complex128` | CayleyTransform returns the value of the Cayley transform (z-i)/(z+i), which maps the upper half-plane conformally onto the open unit disk. |
| `func Cbrt(z complex128) complex128` | Cbrt returns the principal cube root of z, exp(Log(z)/3), which is real and positive for real positive z. |
| `func Conj(z complex128) complex128` | Conj returns the complex conjugate of z. |
| `func ContourIntegral(f Function, gamma func(float64) complex128, n int) complex128` | ContourIntegral approximates the integral of f along a general contour given by the parametrization gamma over t in [0, 1]. |
| `func Cos(z complex128) complex128` | Cos returns the cosine of z. |
| `func Cosh(z complex128) complex128` | Cosh returns the hyperbolic cosine of z. |
| `func Cot(z complex128) complex128` | Cot returns the cotangent of z, cos(z)/sin(z). |
| `func Coth(z complex128) complex128` | Coth returns the hyperbolic cotangent of z, cosh(z)/sinh(z). |
| `func CountZeros(f Function, center complex128, radius float64, n int) int` | CountZeros returns the number of zeros of an analytic function f inside the circle of the given center and radius, counted with multiplicity, by… |
| `func CrossRatio(z, z1, z2, z3 complex128) complex128` | CrossRatio returns the cross-ratio (z-z1)(z2-z3) / ((z-z3)(z2-z1)) of the four points. |
| `func Csc(z complex128) complex128` | Csc returns the cosecant of z, 1/sin(z). |
| `func Digamma(z complex128) complex128` | Digamma returns the logarithmic derivative of the gamma function, psi(z) = Gamma'(z)/Gamma(z), using recurrence to raise the real part and then the… |
| `func Erf(z complex128) complex128` | Erf returns the error function of z, 2/sqrt(pi) times the integral of exp(-t^2) from 0 to z, summed from its Taylor series. |
| `func Erfc(z complex128) complex128` | Erfc returns the complementary error function 1 - Erf(z). |
| `func Exp(z complex128) complex128` | Exp returns e**z. |
| `func Factorial(z complex128) complex128` | Factorial returns z! defined as Gamma(z+1); for a non-negative integer z it equals the ordinary factorial. |
| `func Gamma(z complex128) complex128` | Gamma returns the value of the Euler gamma function at z, computed with the Lanczos approximation and the reflection formula for Re(z) < 0.5. |
| `func Im(z complex128) float64` | Im returns the imaginary part of z. |
| `func IntegrateCircle(f Function, center complex128, radius float64, n int) complex128` | IntegrateCircle approximates the contour integral of f around the positively oriented circle of the given center and radius, using the composite… |
| `func IntegratePolygon(f Function, vertices []complex128, nPerEdge int) complex128` | IntegratePolygon approximates the integral of f along the closed polygonal contour through the given vertices (the last vertex is joined back to the… |
| `func IntegrateSegment(f Function, a, b complex128, n int) complex128` | IntegrateSegment approximates the integral of f along the straight line segment from a to b using the composite Simpson rule with n subintervals (n… |
| `func InverseCayleyTransform(w complex128) complex128` | InverseCayleyTransform returns i*(1+w)/(1-w), the inverse of the Cayley transform, mapping the unit disk back to the upper half-plane. |
| `func IsZero(z complex128, tol float64) bool` | IsZero reports whether \|z\| <= tol. |
| `func JoukowskiMap(z complex128) complex128` | JoukowskiMap returns the Joukowski map (z + 1/z)/2, which maps circles about the origin to confocal ellipses and is used in classical airfoil theory. |
| `func LaurentCoefficient(f Function, center complex128, k int, radius float64, n int) complex128` | LaurentCoefficient returns the k-th Laurent coefficient c_k of f expanded about center, computed from the contour integral c_k = 1/(2*pi*i) *… |
| `func LaurentCoefficients(f Function, center complex128, lo, hi int, radius float64, n int) []complex128` | LaurentCoefficients returns the Laurent coefficients c_k of f about center for k running from lo to hi inclusive, as a slice indexed so that… |
| `func Log(z complex128) complex128` | Log returns the principal branch of the natural logarithm of z, with imaginary part (the argument) in (-pi, pi]. |
| `func LogBranch(z complex128, k int) complex128` | LogBranch returns the value of the natural logarithm of z on branch k, that is Log(z) + 2*pi*i*k. |
| `func LogGamma(z complex128) complex128` | LogGamma returns the principal branch of the logarithm of the gamma function at z. |
| `func NthRoots(z complex128, n int) []complex128` | NthRoots returns all n distinct n-th roots of z, ordered by increasing argument starting from the principal root. |
| `func Polar(r, theta float64) complex128` | Polar constructs a complex number from a modulus r and an argument theta, i.e. |
| `func Pow(z, w complex128) complex128` | Pow returns z raised to the power w using the principal branch, exp(w*Log(z)). |
| `func PowerSeriesEval(coeffs []complex128, center, z complex128) complex128` | PowerSeriesEval evaluates the power series with the given coefficients about center at the point z, returning sum_j coeffs`j` * (z-center)^j via… |
| `func Re(z complex128) float64` | Re returns the real part of z. |
| `func Reciprocal(z complex128) complex128` | Reciprocal returns 1/z. |
| `func Rect(re, im float64) complex128` | Rect constructs a complex number from its real and imaginary parts. |
| `func Residue(f Function, z0 complex128, radius float64, n int) complex128` | Residue approximates the residue of f at the isolated singularity z0 by integrating f around a small circle of the given radius and dividing by… |
| `func ResidueOrderM(f Function, z0 complex128, m int, radius float64, n int) complex128` | ResidueOrderM returns the residue of f at a pole z0 of order m using the derivative formula Res = lim_{z->z0} 1/(m-1)! d^{m-1}/dz^{m-1} [(z-z0)^m… |
| `func ResidueSimplePole(f Function, z0 complex128, eps float64) complex128` | ResidueSimplePole returns the residue of f at a simple pole z0, computed as the limit of (z-z0)f(z) as z approaches z0 via a small symmetric average… |
| `func RisingFactorial(z complex128, n int) complex128` | RisingFactorial returns the Pochhammer symbol (z)_n = z(z+1)...(z+n-1), computed as Gamma(z+n)/Gamma(z). |
| `func Sec(z complex128) complex128` | Sec returns the secant of z, 1/cos(z). |
| `func Sign(z complex128) complex128` | Sign returns z/\|z\|, the unit complex number with the same argument as z, or 0 when z is 0. |
| `func Sin(z complex128) complex128` | Sin returns the sine of z. |
| `func Sinh(z complex128) complex128` | Sinh returns the hyperbolic sine of z. |
| `func Sqrt(z complex128) complex128` | Sqrt returns the principal square root of z, the branch whose result lies in the right half-plane (or on the non-negative imaginary axis). |
| `func Tan(z complex128) complex128` | Tan returns the tangent of z. |
| `func Tanh(z complex128) complex128` | Tanh returns the hyperbolic tangent of z. |
| `func TaylorCoefficient(f Function, center complex128, k int, radius float64, n int) complex128` | TaylorCoefficient returns the k-th Taylor coefficient of an analytic function f about center, a_k = f^(k)(center)/k!, computed from the same contour… |
| `func TaylorCoefficients(f Function, center complex128, count int, radius float64, n int) []complex128` | TaylorCoefficients returns the first count Taylor coefficients (orders 0 through count-1) of an analytic function f about center. |
| `func WindingNumber(gamma func(float64) complex128, z0 complex128, n int) int` | WindingNumber returns the winding number (index) of the closed curve gamma about the point z0, computed by accumulating the change in the argument of… |
| `func Zeta(s complex128) complex128` | Zeta returns the Riemann zeta function at s via the Dirichlet eta function and analytic continuation, eta(s) = (1 - 2^(1-s)) zeta(s). |

### Types

| Type | What it is |
| --- | --- |
| `Function` | Function is a complex-valued function of a single complex variable. |
| `Mobius` | Mobius represents a Mobius (linear fractional) transformation z -> (A*z + B) / (C*z + D). |

<details>
<summary><code>Mobius</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func IdentityMobius() Mobius` | IdentityMobius returns the identity transformation z -> z. |
| `func InversionMap() Mobius` | InversionMap returns the Mobius map z -> 1/z. |
| `func MobiusFromPoints(z1, z2, z3, w1, w2, w3 complex128) Mobius` | MobiusFromPoints returns the unique Mobius transformation that sends the three distinct points z1, z2, z3 to w1, w2, w3 respectively. |
| `func NewMobius(a, b, c, d complex128) Mobius` | NewMobius constructs the Mobius transformation with the given coefficients. |
| `func RotationMap(theta float64) Mobius` | RotationMap returns the rotation z -> e^{i*theta} * z about the origin. |
| `func ScalingMap(a complex128) Mobius` | ScalingMap returns the Mobius map z -> a*z. |
| `func TranslationMap(b complex128) Mobius` | TranslationMap returns the Mobius map z -> z + b. |
| `func (m Mobius) Apply(z complex128) complex128` | Apply evaluates the transformation at z. |
| `func (m Mobius) Compose(other Mobius) Mobius` | Compose returns the transformation m(other(z)), the composition of m after other. |
| `func (m Mobius) Determinant() complex128` | Determinant returns A*D - B*C. A Mobius map is invertible exactly when this is non-zero. |
| `func (m Mobius) FixedPoints() []complex128` | FixedPoints returns the fixed points of the transformation, the solutions of C*z^2 + (D-A)*z - B = 0. |
| `func (m Mobius) Inverse() Mobius` | Inverse returns the inverse transformation z -> (D*z - B) / (-C*z + A). |
| `func (m Mobius) Normalize() Mobius` | Normalize returns an equivalent transformation whose coefficient matrix has determinant 1, dividing all coefficients by a square root of the current… |

</details>

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/complexanalysis).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
