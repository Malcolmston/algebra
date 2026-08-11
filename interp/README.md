# interp — implements interpolation and approximation of one- and multi-dimensional sampled data

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/interp.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/interp)

Package interp implements interpolation and approximation of one- and
multi-dimensional sampled data.

- Piecewise-linear interpolation (LinearInterp), together with the raw scalar
helpers Lerp, InverseLerp and LinearAt. - Nearest-neighbour (NearestInterp)
and piecewise-constant step interpolation (StepInterp). - Cubic splines with
natural, clamped and not-a-knot boundary conditions (CubicSpline), each
supporting evaluation, the first three derivatives and definite integration in
closed form. - The Akima spline (AkimaSpline), which suppresses the spurious
oscillations a global cubic spline can exhibit near sharp features. - Monotone
piecewise cubic Hermite interpolation (PCHIP), which preserves the shape and
monotonicity of the underlying data. - Bilinear and trilinear interpolation on
rectilinear grids (BilinearGrid, TrilinearGrid) with the scalar helpers
BilinearAt and TrilinearAt.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/interp"
```

## Usage

This is the package's own `ExampleBezierElevate`, so it compiles and its output is
asserted on every `go test ./interp/`.

```go
ctrl := [][]float64{{0, 0}, {1, 2}, {2, 0}}
	up := interp.BezierElevate(ctrl)
	fmt.Println(len(up))
	a := interp.BezierPoint(ctrl, 0.37)
	b := interp.BezierPoint(up, 0.37)
	fmt.Printf("%.10f %.10f\n", math.Abs(a[0]-b[0]), math.Abs(a[1]-b[1]))
```

```
4
0.0000000000 0.0000000000
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func BSplineBasis(i, degree int, knots []float64, t float64) float64` | BSplineBasis returns the value at t of the i-th B-spline basis function of the given degree over the knot vector knots, computed with the Cox-de Boor… |
| `func BarycentricEval(x, y, w []float64, xq float64) float64` | BarycentricEval evaluates the interpolating polynomial through (x`i`, y`i`) at xq using the second (true) barycentric formula with the supplied… |
| `func BarycentricWeights(x []float64) []float64` | BarycentricWeights returns the barycentric weights w`j` = 1 / prod_{k!=j} (x`j` - x`k`) for the distinct abscissae x. |
| `func BernsteinBasis(n, i int, t float64) float64` | BernsteinBasis returns the value at t of the i-th Bernstein basis polynomial of degree n, C(n, i) t^i (1-t)^{n-i}. |
| `func BezierCurve(ctrl [][]float64, ts []float64) [][]float64` | BezierCurve evaluates the vector-valued Bezier curve at every parameter in ts and returns the resulting points in order. |
| `func BezierDerivative(ctrl [][]float64, t float64) []float64` | BezierDerivative returns the derivative vector dC/dt of the Bezier curve at parameter t. |
| `func BezierElevate(ctrl [][]float64) [][]float64` | BezierElevate returns the control points of the degree-elevated Bezier curve. |
| `func BezierPoint(ctrl [][]float64, t float64) []float64` | BezierPoint evaluates the vector-valued Bezier curve with the given control points at parameter t in [0, 1] using de Casteljau's algorithm. |
| `func BilinearAt(x0, x1, y0, y1, q00, q10, q01, q11, x, y float64) float64` | BilinearAt returns the bilinear interpolation at (x, y) of the four corner values of the rectangle [x0, x1] x [y0, y1]: q00 at (x0, y0), q10 at (x1,… |
| `func ChebyshevCoefficients(f func(float64) float64, n int, a, b float64) []float64` | ChebyshevCoefficients samples f at the n Chebyshev roots on [a, b] and returns the n coefficients c of the Chebyshev series approximation f(x) ~=… |
| `func ChebyshevNodes(n int) []float64` | ChebyshevNodes returns the n roots of the degree-n Chebyshev polynomial of the first kind on [-1, 1], namely cos((2k+1)pi / (2n)) for k = 0..n-1. |
| `func ChebyshevNodesInterval(n int, a, b float64) []float64` | ChebyshevNodesInterval returns the n Chebyshev roots mapped affinely from [-1, 1] onto [a, b], in increasing order. |
| `func ChebyshevNodesSecondKind(n int) []float64` | ChebyshevNodesSecondKind returns the n Chebyshev-Gauss-Lobatto points cos(k*pi / (n-1)) for k = 0..n-1 on [-1, 1] (the extrema of T_{n-1}), returned… |
| `func ChebyshevT(n int, x float64) float64` | ChebyshevT returns the value of the Chebyshev polynomial of the first kind of degree n at x, computed by the stable three-term recurrence. |
| `func Clamp(v, lo, hi float64) float64` | Clamp constrains v to the closed interval [lo, hi]. |
| `func ClenshawEval(c []float64, y float64) float64` | ClenshawEval evaluates the Chebyshev series with coefficients c at the point y in [-1, 1] using Clenshaw recurrence, following the convention sum_{j}… |
| `func CubicHermiteSegment(p0, m0, p1, m1, t float64) float64` | CubicHermiteSegment evaluates a single cubic Hermite segment on the unit interval at parameter t in [0, 1], given endpoint values p0, p1 and endpoint… |
| `func DFTReal(y []float64) (re, im []float64)` | DFTReal returns the discrete Fourier transform of the real input y, with re`k` and im`k` holding the real and imaginary parts of X_k = sum_j y_j… |
| `func DeBoor(degree int, knots []float64, ctrl [][]float64, t float64) []float64` | DeBoor evaluates the vector-valued B-spline curve of the given degree with knot vector knots and control points ctrl at parameter t using de Boor's… |
| `func DeCasteljau(points []float64, t float64) float64` | DeCasteljau evaluates the scalar Bezier curve with the given control ordinates at parameter t in [0, 1] using de Casteljau's algorithm. |
| `func DividedDifferenceTable(x, y []float64) [][]float64` | DividedDifferenceTable returns the full lower-triangular table of divided differences. |
| `func DividedDifferences(x, y []float64) []float64` | DividedDifferences returns the Newton coefficients c of the interpolating polynomial written as c[0] + c[1]\(x-x0) + c[2]\(x-x0)(x-x1) + ... |
| `func HermiteEval(x, y, dy []float64, xq float64) float64` | HermiteEval evaluates, at xq, the osculating (Hermite) polynomial that matches both the values y and the first derivatives dy at the abscissae x. |
| `func InverseLerp(a, b, v float64) float64` | InverseLerp returns the parameter t for which Lerp(a, b, t) == v, that is (v-a)/(b-a). |
| `func LagrangeBasis(x []float64, i int, xq float64) float64` | LagrangeBasis returns the value at xq of the i-th Lagrange basis polynomial L_i built from the abscissae x. |
| `func LagrangeEval(x, y []float64, xq float64) float64` | LagrangeEval evaluates, at xq, the unique polynomial of degree < len(x) that passes through the samples (x`i`, y`i`) using the Lagrange form. |
| `func LeastSquaresLine(x, y []float64) (slope, intercept float64, err error)` | LeastSquaresLine returns the slope and intercept of the ordinary least-squares line y = slope*x + intercept fitting the samples. |
| `func Lerp(a, b, t float64) float64` | Lerp returns the linear interpolation between a and b at parameter t, that is a + t*(b-a). |
| `func LinearAt(x0, y0, x1, y1, x float64) float64` | LinearAt returns the value at x of the straight line through the points (x0, y0) and (x1, y1). |
| `func LinearInterpolate(x, y []float64, xq float64) (float64, error)` | LinearInterpolate returns the piecewise-linear interpolation of the samples (x`i`, y`i`) at the single query point xq. |
| `func NearestInterpolate(x, y []float64, xq float64) (float64, error)` | NearestInterpolate returns the nearest-neighbour interpolation of the samples (x`i`, y`i`) at the single query point xq. |
| `func NevilleEval(x, y []float64, xq float64) float64` | NevilleEval evaluates, at xq, the interpolating polynomial through the samples (x`i`, y`i`) using Neville's iterated linear interpolation. |
| `func NevilleTable(x, y []float64, xq float64) [][]float64` | NevilleTable returns the full Neville tableau evaluated at xq. |
| `func NevilleWithError(x, y []float64, xq float64) (value, errEst float64)` | NevilleWithError evaluates the interpolating polynomial at xq and also returns the magnitude of the last correction, a heuristic estimate of the… |
| `func NewtonEval(x, coeffs []float64, xq float64) float64` | NewtonEval evaluates the Newton form with nodes x and coefficients coeffs at xq using Horner's scheme. |
| `func OpenUniformKnots(numCtrl, degree int) []float64` | OpenUniformKnots returns a clamped (open uniform) knot vector on [0, 1] for a B-spline of the given degree with numCtrl control points, so that the… |
| `func PolyDerivativeCoeffs(coeffs []float64) []float64` | PolyDerivativeCoeffs returns the ascending coefficients of the derivative of the polynomial given by the ascending coefficients coeffs. |
| `func PolyFit(x, y []float64, degree int) ([]float64, error)` | PolyFit returns the ascending coefficients of the degree-`degree` least-squares polynomial that best fits the samples (x`i`, y`i`) in the ordinary… |
| `func PolyFitWeighted(x, y, w []float64, degree int) ([]float64, error)` | PolyFitWeighted returns the ascending coefficients of the degree-`degree` weighted least-squares polynomial fit, where w`i` is the weight of sample i. |
| `func PolyIntegralCoeffs(coeffs []float64, c float64) []float64` | PolyIntegralCoeffs returns the ascending coefficients of an antiderivative of the polynomial given by coeffs, using constant of integration c. |
| `func PolyVal(coeffs []float64, x float64) float64` | PolyVal evaluates the polynomial with ascending coefficients coeffs (so coeffs[0] is the constant term) at x using Horner's method. |
| `func RMSE(x, y, coeffs []float64) float64` | RMSE returns the root-mean-square error of the polynomial with ascending coefficients coeffs against the samples (x`i`, y`i`). |
| `func RSquared(x, y, coeffs []float64) float64` | RSquared returns the coefficient of determination R^2 of the polynomial with ascending coefficients coeffs against the samples (x`i`, y`i`). |
| `func RationalInterpEval(x, y []float64, xq float64) (float64, error)` | RationalInterpEval performs Bulirsch-Stoer diagonal rational function interpolation of the samples (x`i`, y`i`) at xq, returning the interpolated… |
| `func SearchInterval(xs []float64, x float64) int` | SearchInterval returns the index i such that xs`i` <= x <= xs[i+1], clamped so that i lies in [0, len(xs)-2]. |
| `func SmoothStep(edge0, edge1, x float64) float64` | SmoothStep returns the cubic Hermite smoothstep 3u^2 - 2u^3 where u is x mapped from [edge0, edge1] to [0, 1] and clamped. |
| `func SmootherStep(edge0, edge1, x float64) float64` | SmootherStep returns Perlin's quintic smoothstep 6u^5 - 15u^4 + 10u^3 where u is x mapped from [edge0, edge1] to [0, 1] and clamped. |
| `func ThieleEval(x, a []float64, xq float64) float64` | ThieleEval evaluates Thiele's continued fraction with nodes x and coefficients a (as returned by ThieleReciprocalDifferences) at xq. |
| `func ThieleReciprocalDifferences(x, y []float64) []float64` | ThieleReciprocalDifferences returns the coefficients a of Thiele's continued fraction for the samples (x`i`, y`i`): |
| `func TrigCoefficients(y []float64) (a, b []float64)` | TrigCoefficients returns the real trigonometric interpolation coefficients a and b (each of length N/2+1) for samples y taken at the equally spaced… |
| `func TrigInterpEval(y []float64, xq float64) float64` | TrigInterpEval evaluates, at xq, the trigonometric polynomial that interpolates the samples y taken at the equally spaced nodes x_j = 2*pi*j/len(y)… |
| `func TrilinearAt(x0, x1, y0, y1, z0, z1, c000, c100, c010, c110, c001, c101, c011, c111, x, y, z float64) float64` | TrilinearAt returns the trilinear interpolation at (x, y, z) of the eight corner values of the box [x0, x1] x [y0, y1] x [z0, z1]. |
| `func UniformKnots(numCtrl, degree int) []float64` | UniformKnots returns a uniform (periodic-style) knot vector 0, 1, ..., m for a B-spline of the given degree with numCtrl control points. |
| `func VandermondeSolve(x, y []float64) ([]float64, error)` | VandermondeSolve returns the ascending coefficients of the unique polynomial of degree len(x)-1 that passes exactly through the samples (x`i`, y`i`),… |

### Types

| Type | What it is |
| --- | --- |
| `AkimaSpline` | AkimaSpline is a piecewise-cubic Hermite interpolant using Akima's local slope estimates. |
| `BSpline` | BSpline is a reusable vector-valued B-spline curve. |
| `BarycentricInterpolator` | BarycentricInterpolator is a reusable interpolant that evaluates in linear time via precomputed barycentric weights. |
| `Bezier` | Bezier is a reusable vector-valued Bezier curve of arbitrary degree. |
| `BilinearGrid` | BilinearGrid is a bilinear interpolant on a rectilinear grid. |
| `ChebyshevInterpolator` | ChebyshevInterpolator approximates a function on an interval [a, b] by truncated Chebyshev series obtained from samples at the Chebyshev roots. |
| `CubicSpline` | CubicSpline is a C2-continuous piecewise-cubic interpolant. |
| `HermiteInterpolator` | HermiteInterpolator is a reusable osculating polynomial matching values and first derivatives at each node. |
| `LagrangeInterpolator` | LagrangeInterpolator is a reusable interpolant in Lagrange form built from a fixed set of samples. |
| `LinearInterp` | LinearInterp is a piecewise-linear interpolant through a set of samples. |
| `NearestInterp` | NearestInterp is a nearest-neighbour interpolant: it returns the ordinate of whichever sample abscissa is closest to the query point, breaking ties… |
| `NewtonInterpolator` | NewtonInterpolator is a reusable interpolant in Newton form that also supports incremental extension via AddPoint. |
| `PCHIP` | PCHIP is a monotone piecewise cubic Hermite interpolant using the Fritsch-Carlson slope limiter. |
| `PolynomialFit` | PolynomialFit bundles the ascending coefficients of a least-squares fit with convenient evaluation. |
| `StepInterp` | StepInterp is a piecewise-constant (zero-order-hold) interpolant. |
| `ThieleInterpolator` | ThieleInterpolator is a reusable rational interpolant in Thiele continued-fraction form. |
| `TrigInterpolator` | TrigInterpolator is a reusable trigonometric interpolant of equally spaced periodic samples. |
| `TrilinearGrid` | TrilinearGrid is a trilinear interpolant on a rectilinear grid. |

<details>
<summary><code>AkimaSpline</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewAkimaSpline(x, y []float64) (*AkimaSpline, error)` | NewAkimaSpline builds an Akima spline through the samples (x`i`, y`i`). |
| `func (s *AkimaSpline) Domain() (lo, hi float64)` | Domain returns the first and last abscissae of the samples. |
| `func (s *AkimaSpline) Eval(x float64) float64` | Eval returns the spline value at x. |
| `func (s *AkimaSpline) EvalDerivative(x float64) float64` | EvalDerivative returns the first derivative of the spline at x. |
| `func (s *AkimaSpline) EvalSecondDerivative(x float64) float64` | EvalSecondDerivative returns the second derivative of the spline at x. |
| `func (s *AkimaSpline) EvalSlice(xs []float64) []float64` | EvalSlice returns the spline evaluated at every point of xs. |
| `func (s *AkimaSpline) Integral(a, b float64) float64` | Integral returns the definite integral of the spline from a to b, with limits clamped to the sampled domain. |
| `func (s *AkimaSpline) Len() int` | Len returns the number of samples. |
| `func (s *AkimaSpline) Xs() []float64` | Xs returns a copy of the sample abscissae. |
| `func (s *AkimaSpline) Ys() []float64` | Ys returns a copy of the sample ordinates. |

</details>

<details>
<summary><code>BSpline</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewBSpline(degree int, knots []float64, ctrl [][]float64) (*BSpline, error)` | NewBSpline builds a B-spline of the given degree from the knot vector and control points. |
| `func (s *BSpline) Degree() int` | Degree reports the polynomial degree of the B-spline. |
| `func (s *BSpline) Domain() (lo, hi float64)` | Domain returns the valid parameter range [lo, hi] of the curve. |
| `func (s *BSpline) Eval(t float64) []float64` | Eval returns the point on the curve at parameter t, clamped to the domain. |

</details>

<details>
<summary><code>BarycentricInterpolator</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewBarycentricInterpolator(x, y []float64) (*BarycentricInterpolator, error)` | NewBarycentricInterpolator builds a BarycentricInterpolator from the samples (x`i`, y`i`). |
| `func (b *BarycentricInterpolator) Eval(xq float64) float64` | Eval returns the interpolating polynomial evaluated at xq. |
| `func (b *BarycentricInterpolator) EvalSlice(xs []float64) []float64` | EvalSlice evaluates the interpolant at every element of xs. |
| `func (b *BarycentricInterpolator) Weights() []float64` | Weights returns a copy of the barycentric weights. |

</details>

<details>
<summary><code>Bezier</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewBezier(ctrl [][]float64) (*Bezier, error)` | NewBezier builds a Bezier curve from the given control points. |
| `func (b *Bezier) Degree() int` | Degree reports the polynomial degree of the curve, one less than the number of control points. |
| `func (b *Bezier) Derivative(t float64) []float64` | Derivative returns the derivative vector of the curve at parameter t. |
| `func (b *Bezier) Eval(t float64) []float64` | Eval returns the point on the curve at parameter t in [0, 1]. |

</details>

<details>
<summary><code>BilinearGrid</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewBilinearGrid(x, y []float64, z [][]float64) (*BilinearGrid, error)` | NewBilinearGrid builds a bilinear grid interpolant. |
| `func (g *BilinearGrid) Dims() (nx, ny int)` | Dims returns the number of grid points along the x and y axes. |
| `func (g *BilinearGrid) Domain() (xMin, xMax, yMin, yMax float64)` | Domain returns the grid extent as (xMin, xMax, yMin, yMax). |
| `func (g *BilinearGrid) Eval(x, y float64) float64` | Eval returns the bilinearly interpolated value at (x, y). |
| `func (g *BilinearGrid) EvalClamped(x, y float64) float64` | EvalClamped returns the interpolated value at (x, y) with the query point first clamped to the grid extent, so no extrapolation occurs. |

</details>

<details>
<summary><code>ChebyshevInterpolator</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewChebyshevInterpolator(f func(float64) float64, n int, a, b float64) (*ChebyshevInterpolator, error)` | NewChebyshevInterpolator constructs a ChebyshevInterpolator of order n for f on [a, b] by sampling f at the n Chebyshev roots. |
| `func (ci *ChebyshevInterpolator) Coefficients() []float64` | Coefficients returns a copy of the Chebyshev coefficients. |
| `func (ci *ChebyshevInterpolator) Eval(xq float64) float64` | Eval returns the Chebyshev approximation at xq. |
| `func (ci *ChebyshevInterpolator) EvalSlice(xs []float64) []float64` | EvalSlice evaluates the approximation at every element of xs. |

</details>

<details>
<summary><code>CubicSpline</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewClampedCubicSpline(x, y []float64, dLeft, dRight float64) (*CubicSpline, error)` | NewClampedCubicSpline builds a cubic spline whose first derivative equals dLeft at the first node and dRight at the last node. |
| `func NewNaturalCubicSpline(x, y []float64) (*CubicSpline, error)` | NewNaturalCubicSpline builds a cubic spline with natural boundary conditions, that is with zero second derivative at both ends. |
| `func NewNotAKnotCubicSpline(x, y []float64) (*CubicSpline, error)` | NewNotAKnotCubicSpline builds a cubic spline with not-a-knot boundary conditions, requiring the third derivative to be continuous across the first… |
| `func (s *CubicSpline) Domain() (lo, hi float64)` | Domain returns the first and last abscissae of the samples. |
| `func (s *CubicSpline) Eval(x float64) float64` | Eval returns the spline value at x. |
| `func (s *CubicSpline) EvalClamped(x float64) float64` | EvalClamped returns the spline value at x with x first clamped to the sampled domain, so no extrapolation occurs. |
| `func (s *CubicSpline) EvalDerivative(x float64) float64` | EvalDerivative returns the first derivative of the spline at x. |
| `func (s *CubicSpline) EvalSecondDerivative(x float64) float64` | EvalSecondDerivative returns the second derivative of the spline at x. |
| `func (s *CubicSpline) EvalSlice(xs []float64) []float64` | EvalSlice returns the spline evaluated at every point of xs. |
| `func (s *CubicSpline) EvalThirdDerivative(x float64) float64` | EvalThirdDerivative returns the third derivative of the spline at x, constant within each interval. |
| `func (s *CubicSpline) Integral(a, b float64) float64` | Integral returns the definite integral of the spline from a to b, with the limits clamped to the sampled domain. |
| `func (s *CubicSpline) Len() int` | Len returns the number of samples. |
| `func (s *CubicSpline) Xs() []float64` | Xs returns a copy of the sample abscissae. |
| `func (s *CubicSpline) Ys() []float64` | Ys returns a copy of the sample ordinates. |

</details>

<details>
<summary><code>HermiteInterpolator</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewHermiteInterpolator(x, y, dy []float64) (*HermiteInterpolator, error)` | NewHermiteInterpolator builds a HermiteInterpolator matching the values y and first derivatives dy at the abscissae x. |
| `func (h *HermiteInterpolator) Coefficients() []float64` | Coefficients returns a copy of the Newton coefficients over the doubled node sequence. |
| `func (h *HermiteInterpolator) Eval(xq float64) float64` | Eval returns the osculating polynomial evaluated at xq. |
| `func (h *HermiteInterpolator) EvalDerivative(xq float64) float64` | EvalDerivative returns the first derivative of the osculating polynomial at xq. |

</details>

<details>
<summary><code>LagrangeInterpolator</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewLagrangeInterpolator(x, y []float64) (*LagrangeInterpolator, error)` | NewLagrangeInterpolator builds a LagrangeInterpolator from the samples (x`i`, y`i`). |
| `func (l *LagrangeInterpolator) Degree() int` | Degree reports the degree of the interpolating polynomial, which is one less than the number of samples. |
| `func (l *LagrangeInterpolator) Eval(xq float64) float64` | Eval returns the interpolating polynomial evaluated at xq. |
| `func (l *LagrangeInterpolator) EvalSlice(xs []float64) []float64` | EvalSlice evaluates the interpolant at every element of xs. |
| `func (l *LagrangeInterpolator) Len() int` | Len reports the number of samples used to build the interpolant. |

</details>

<details>
<summary><code>LinearInterp</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewLinearInterp(x, y []float64) (*LinearInterp, error)` | NewLinearInterp builds a piecewise-linear interpolant through the samples (x`i`, y`i`). |
| `func (l *LinearInterp) Domain() (lo, hi float64)` | Domain returns the first and last abscissae of the samples. |
| `func (l *LinearInterp) Eval(x float64) float64` | Eval returns the interpolated value at x. |
| `func (l *LinearInterp) EvalDerivative(x float64) float64` | EvalDerivative returns the slope of the segment containing x. |
| `func (l *LinearInterp) EvalSlice(xs []float64) []float64` | EvalSlice returns the interpolant evaluated at every point of xs. |
| `func (l *LinearInterp) Integral(a, b float64) float64` | Integral returns the definite integral from a to b, with limits clamped to the sampled domain. |
| `func (l *LinearInterp) Len() int` | Len returns the number of samples. |
| `func (l *LinearInterp) Xs() []float64` | Xs returns a copy of the sample abscissae. |
| `func (l *LinearInterp) Ys() []float64` | Ys returns a copy of the sample ordinates. |

</details>

<details>
<summary><code>NearestInterp</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewNearestInterp(x, y []float64) (*NearestInterp, error)` | NewNearestInterp builds a nearest-neighbour interpolant. |
| `func (n *NearestInterp) Domain() (lo, hi float64)` | Domain returns the first and last abscissae of the samples. |
| `func (n *NearestInterp) Eval(x float64) float64` | Eval returns the ordinate of the sample nearest to x. |
| `func (n *NearestInterp) Len() int` | Len returns the number of samples. |
| `func (n *NearestInterp) Xs() []float64` | Xs returns a copy of the sample abscissae. |
| `func (n *NearestInterp) Ys() []float64` | Ys returns a copy of the sample ordinates. |

</details>

<details>
<summary><code>NewtonInterpolator</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewNewtonInterpolator(x, y []float64) (*NewtonInterpolator, error)` | NewNewtonInterpolator builds a NewtonInterpolator from the samples (x`i`, y`i`). |
| `func (n *NewtonInterpolator) AddPoint(xi, yi float64) error` | AddPoint extends the interpolant with one additional sample (xi, yi) in O(len) time, updating the Newton coefficients in place. |
| `func (n *NewtonInterpolator) Coefficients() []float64` | Coefficients returns a copy of the Newton coefficients. |
| `func (n *NewtonInterpolator) Eval(xq float64) float64` | Eval returns the interpolating polynomial evaluated at xq. |
| `func (n *NewtonInterpolator) EvalSlice(xs []float64) []float64` | EvalSlice evaluates the interpolant at every element of xs. |
| `func (n *NewtonInterpolator) Len() int` | Len reports the number of samples currently held by the interpolant. |

</details>

<details>
<summary><code>PCHIP</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewPCHIP(x, y []float64) (*PCHIP, error)` | NewPCHIP builds a monotone piecewise cubic Hermite interpolant through the samples (x`i`, y`i`). |
| `func (s *PCHIP) Domain() (lo, hi float64)` | Domain returns the first and last abscissae of the samples. |
| `func (s *PCHIP) Eval(x float64) float64` | Eval returns the interpolated value at x. |
| `func (s *PCHIP) EvalDerivative(x float64) float64` | EvalDerivative returns the first derivative of the interpolant at x. |
| `func (s *PCHIP) EvalSecondDerivative(x float64) float64` | EvalSecondDerivative returns the second derivative of the interpolant at x. |
| `func (s *PCHIP) EvalSlice(xs []float64) []float64` | EvalSlice returns the interpolant evaluated at every point of xs. |
| `func (s *PCHIP) Integral(a, b float64) float64` | Integral returns the definite integral of the interpolant from a to b, with limits clamped to the sampled domain. |
| `func (s *PCHIP) IsMonotone() bool` | IsMonotone reports whether the underlying sample ordinates are monotone (entirely non-decreasing or entirely non-increasing), the condition under… |
| `func (s *PCHIP) Len() int` | Len returns the number of samples. |
| `func (s *PCHIP) Xs() []float64` | Xs returns a copy of the sample abscissae. |
| `func (s *PCHIP) Ys() []float64` | Ys returns a copy of the sample ordinates. |

</details>

<details>
<summary><code>PolynomialFit</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewPolynomialFit(x, y []float64, degree int) (*PolynomialFit, error)` | NewPolynomialFit builds a PolynomialFit of the given degree from the samples (x`i`, y`i`) by ordinary least squares. |
| `func (p *PolynomialFit) Coefficients() []float64` | Coefficients returns a copy of the ascending fit coefficients. |
| `func (p *PolynomialFit) Degree() int` | Degree reports the degree of the fitted polynomial. |
| `func (p *PolynomialFit) Eval(x float64) float64` | Eval evaluates the fitted polynomial at x. |

</details>

<details>
<summary><code>StepInterp</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewNextStepInterp(x, y []float64) (*StepInterp, error)` | NewNextStepInterp builds a piecewise-constant interpolant that holds the value of the nearest sample at or above the query point. |
| `func NewPreviousStepInterp(x, y []float64) (*StepInterp, error)` | NewPreviousStepInterp builds a piecewise-constant interpolant that holds the value of the nearest sample at or below the query point (a causal… |
| `func (s *StepInterp) Domain() (lo, hi float64)` | Domain returns the first and last abscissae of the samples. |
| `func (s *StepInterp) Eval(x float64) float64` | Eval returns the piecewise-constant value at x. |
| `func (s *StepInterp) Len() int` | Len returns the number of samples. |
| `func (s *StepInterp) Xs() []float64` | Xs returns a copy of the sample abscissae. |
| `func (s *StepInterp) Ys() []float64` | Ys returns a copy of the sample ordinates. |

</details>

<details>
<summary><code>ThieleInterpolator</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewThieleInterpolator(x, y []float64) (*ThieleInterpolator, error)` | NewThieleInterpolator builds a ThieleInterpolator from the samples (x`i`, y`i`). |
| `func (t *ThieleInterpolator) Eval(xq float64) float64` | Eval returns the rational interpolant evaluated at xq. |
| `func (t *ThieleInterpolator) EvalSlice(xs []float64) []float64` | EvalSlice evaluates the interpolant at every element of xs. |

</details>

<details>
<summary><code>TrigInterpolator</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewTrigInterpolator(y []float64) (*TrigInterpolator, error)` | NewTrigInterpolator builds a TrigInterpolator from samples y taken at the nodes x_j = 2*pi*j/len(y). |
| `func (t *TrigInterpolator) Eval(xq float64) float64` | Eval returns the trigonometric interpolant evaluated at xq. |

</details>

<details>
<summary><code>TrilinearGrid</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewTrilinearGrid(x, y, z []float64, v [][][]float64) (*TrilinearGrid, error)` | NewTrilinearGrid builds a trilinear grid interpolant. |
| `func (g *TrilinearGrid) Dims() (nx, ny, nz int)` | Dims returns the number of grid points along the x, y and z axes. |
| `func (g *TrilinearGrid) Domain() (xMin, xMax, yMin, yMax, zMin, zMax float64)` | Domain returns the grid extent as (xMin, xMax, yMin, yMax, zMin, zMax). |
| `func (g *TrilinearGrid) Eval(x, y, z float64) float64` | Eval returns the trilinearly interpolated value at (x, y, z). |
| `func (g *TrilinearGrid) EvalClamped(x, y, z float64) float64` | EvalClamped returns the interpolated value at (x, y, z) with the query point first clamped to the grid extent, so no extrapolation occurs. |

</details>

### Variables

`ErrTooFewPoints`, `ErrLengthMismatch`, `ErrNotSorted`, `ErrGridShape`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/interp).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
