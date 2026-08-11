# splines — implements splines and parametric curves using nothing but the Go standard library

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/splines.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/splines)

Package splines implements splines and parametric curves using nothing but the
Go standard library.

- Cubic splines with natural, clamped, not-a-knot and periodic end conditions,
together with value, derivative and integral evaluation. - Piecewise cubic
Hermite interpolation, Catmull-Rom / cardinal splines, shape-preserving
monotone interpolation (Fritsch-Carlson) and Akima splines. - Bezier curves
evaluated with de Casteljau's algorithm, plus degree elevation, subdivision,
derivatives and the Bernstein basis. - B-spline curves evaluated with de
Boor's algorithm, knot insertion (Boehm), the Cox-de Boor basis and
derivatives. - NURBS curves (rational B-splines) with control point weights. -
Arc-length computation and arc-length (unit-speed) reparameterisation. -
Tensor-product Bezier and B-spline surface evaluation.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/splines"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func BasisFuns(i int, u float64, p int, U []float64) []float64` | BasisFuns returns the p+1 non-zero B-spline basis function values at u for the span i (as returned by `FindSpan`). |
| `func Bernstein(n, i int, t float64) float64` | Bernstein returns the value of the Bernstein basis polynomial B(i,n)(t) = C(n,i) * t^i * (1-t)^(n-i). |
| `func BernsteinAll(n int, t float64) []float64` | BernsteinAll returns the full vector of degree-n Bernstein basis values at t, computed with a stable triangular recurrence (no explicit powers). |
| `func Binomial(n, k int) float64` | Binomial returns the binomial coefficient C(n, k) computed with an integer multiplicative recurrence that avoids intermediate overflow for moderate n. |
| `func CentripetalParams(pts []Vec, alpha float64) []float64` | CentripetalParams returns parameter values for pts using the exponent alpha on successive chord lengths: t[0]=0 and t`i`=t`i-1`+\|P`i`-P`i-1`\|^alpha. |
| `func ChordLengths(pts []Vec) []float64` | ChordLengths returns the cumulative chord-length parameter values of pts, starting at zero; the final entry equals the total polyline length. |
| `func ClampedKnots(numCtrl, degree int) []float64` | ClampedKnots returns an open uniform (clamped) knot vector so that the curve interpolates its first and last control points. |
| `func CubicHermite(y0, y1, m0, m1, h, t float64) float64` | CubicHermite evaluates a single cubic Hermite segment with endpoints y0, y1, endpoint tangents m0, m1, interval width h and local parameter t in… |
| `func CurveLength(c Curve) float64` | CurveLength returns the total arc length of a curve over its full parameter domain, computed by adaptive Gauss-Legendre quadrature of the speed… |
| `func CurveLengthBetween(c Curve, a, b float64) float64` | CurveLengthBetween returns the arc length of a curve between parameters a and b (a <= b assumed within the domain). |
| `func DersBasisFuns(i int, u float64, p, d int, U []float64) [][]float64` | DersBasisFuns returns the non-zero basis functions and their derivatives up to order d, evaluated at u for span i. |
| `func FindSpan(n, p int, u float64, U []float64) int` | FindSpan returns the knot span index i such that U`i` <= u < U[i+1], for a B-spline of degree p with n+1 control points (so n = len(ctrl)-1) and knot… |
| `func GaussLegendreIntegrate(f func(float64) float64, a, b float64, n int) float64` | GaussLegendreIntegrate approximates the integral of f over [a,b] with the n-point Gauss-Legendre rule. |
| `func GaussLegendreRule(n int) (nodes, weights []float64)` | GaussLegendreRule returns the nodes and weights of the n-point Gauss-Legendre quadrature rule on the interval [-1,1]. |
| `func HermiteBasis00(t float64) float64` | HermiteBasis00 is the Hermite basis function h00(t)=2t^3-3t^2+1, the weight of the start value in a unit cubic Hermite segment. |
| `func HermiteBasis01(t float64) float64` | HermiteBasis01 is the Hermite basis function h01(t)=-2t^3+3t^2, the weight of the end value. |
| `func HermiteBasis10(t float64) float64` | HermiteBasis10 is the Hermite basis function h10(t)=t^3-2t^2+t, the weight of the (scaled) start tangent. |
| `func HermiteBasis11(t float64) float64` | HermiteBasis11 is the Hermite basis function h11(t)=t^3-t^2, the weight of the (scaled) end tangent. |
| `func LinearInterpolate(x, y []float64, q float64) (float64, error)` | LinearInterpolate returns the piecewise-linear interpolant of the samples (x`i`, y`i`) evaluated at q, clamping to the end values outside the data… |
| `func NearestInterpolate(x, y []float64, q float64) (float64, error)` | NearestInterpolate returns the sample ordinate whose abscissa is closest to q (nearest-neighbour interpolation). |
| `func NormalizedChordParams(pts []Vec) []float64` | NormalizedChordParams returns chord-length parameters rescaled to the unit interval [0,1]. |
| `func OneBasisFun(p int, U []float64, i int, u float64) float64` | OneBasisFun returns the value of the single B-spline basis function N(i,p)(u) for knot vector U, implementing Algorithm A2.4. |
| `func PolylineLength(pts []Vec) float64` | PolylineLength returns the total length of the polyline through pts, the sum of the Euclidean distances between consecutive points. |
| `func SolveCyclicTridiagonal(a, b, c, d []float64, alpha, beta float64) ([]float64, error)` | SolveCyclicTridiagonal solves a tridiagonal system that additionally has the corner entries alpha (top-right, coupling the last unknown into the… |
| `func SolveTridiagonal(a, b, c, d []float64) ([]float64, error)` | SolveTridiagonal solves the tridiagonal linear system whose sub-, main- and super-diagonals are a, b and c and whose right-hand side is d, using the… |
| `func StepInterpolate(x, y []float64, q float64) (float64, error)` | StepInterpolate returns a left-continuous piecewise-constant interpolant: for q in [x`i`, x[i+1]) it returns y`i`, clamping outside the range. |
| `func UniformKnots(numCtrl, degree int) []float64` | UniformKnots returns a uniform (non-clamped) knot vector of the correct length len(U) = numCtrl + degree + 1 with unit spacing starting at zero. |
| `func VecsEqual(a, b []Vec, tol float64) bool` | VecsEqual reports whether two vector slices are element-wise equal within absolute tolerance tol. |

### Types

| Type | What it is |
| --- | --- |
| `ArcLengthParam` | ArcLengthParam provides an arc-length (unit-speed) reparameterisation of a curve. |
| `BSplineCurve` | BSplineCurve is a (non-rational) B-spline curve defined by control points, a knot vector and a degree. |
| `BSplineSurface` | BSplineSurface is a tensor-product B-spline surface with independent degrees and knot vectors in the u and v directions. |
| `BezierCurve` | BezierCurve is a polynomial Bezier curve defined by an ordered list of control points of a common dimension. |
| `BezierSurface` | BezierSurface is a tensor-product Bezier surface defined by a rectangular control net. |
| `BoundaryType` | BoundaryType selects the end conditions used when constructing a cubic spline. |
| `CatmullRomCurve` | CatmullRomCurve is a parametric Catmull-Rom spline through a sequence of points. |
| `CubicSpline` | CubicSpline is a piecewise cubic interpolant of one-dimensional data. |
| `Curve` | Curve is the common interface implemented by the parametric curve types in this package. |
| `NURBSCurve` | NURBSCurve is a Non-Uniform Rational B-Spline curve: a B-spline whose control points carry positive weights, letting it represent conic sections such… |
| `NURBSSurface` | NURBSSurface is a tensor-product rational B-spline surface, adding a positive weight to every control point of a `BSplineSurface`. |
| `PiecewiseCubic` | PiecewiseCubic is a one-dimensional interpolant made of independent cubic pieces, each expressed in a local monomial basis. |
| `Vec` | Vec is a point or vector in R^n represented by its coordinates. |

<details>
<summary><code>ArcLengthParam</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewArcLengthParam(c Curve, n int) (*ArcLengthParam, error)` | NewArcLengthParam builds an arc-length parameterisation of c using n sample segments (n >= 1) across the curve's domain. |
| `func (ap *ArcLengthParam) PointAtArcLength(s float64) Vec` | PointAtArcLength returns the point on the curve at cumulative arc length s, i.e. |
| `func (ap *ArcLengthParam) PointAtFraction(f float64) Vec` | PointAtFraction returns the point on the curve at fraction f in [0,1] of the total arc length; f=0 is the start and f=1 is the end. |
| `func (ap *ArcLengthParam) SForT(t float64) float64` | SForT returns the arc length from the start of the curve to parameter t using linear interpolation within the length table. |
| `func (ap *ArcLengthParam) TForS(s float64) float64` | TForS returns the curve parameter t at which the cumulative arc length equals s. |
| `func (ap *ArcLengthParam) TotalLength() float64` | TotalLength returns the total arc length recorded by the parameterisation. |

</details>

<details>
<summary><code>BSplineCurve</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewBSplineCurve(ctrl []Vec, knots []float64, degree int) (*BSplineCurve, error)` | NewBSplineCurve builds a B-spline curve from control points, a non-decreasing knot vector and degree p. |
| `func (b *BSplineCurve) ControlPoints() []Vec` | ControlPoints returns an independent copy of the control points. |
| `func (b *BSplineCurve) Degree() int` | Degree returns the degree of the curve. |
| `func (b *BSplineCurve) Dim() int` | Dim returns the spatial dimension of the control points. |
| `func (b *BSplineCurve) Domain() (lo, hi float64)` | Domain returns the valid parameter interval [U`p`, U[n+1]] of the curve. |
| `func (b *BSplineCurve) Eval(u float64) Vec` | Eval returns the point on the curve at parameter u using the basis functions. |
| `func (b *BSplineCurve) EvalDerivative(u float64) Vec` | EvalDerivative returns the first derivative vector of the curve at parameter u. |
| `func (b *BSplineCurve) EvalDerivatives(u float64, d int) []Vec` | EvalDerivatives returns the point on the curve together with its derivatives up to order d as a slice of length d+1 (index 0 is the point). |
| `func (b *BSplineCurve) InsertKnot(u float64, r int) (*BSplineCurve, error)` | InsertKnot returns a new curve with the knot u inserted r times using Boehm's algorithm (Algorithm A5.1). |
| `func (b *BSplineCurve) Knots() []float64` | Knots returns an independent copy of the knot vector. |

</details>

<details>
<summary><code>BSplineSurface</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewBSplineSurface(ctrl [][]Vec, knotsU, knotsV []float64, degU, degV int) (*BSplineSurface, error)` | NewBSplineSurface builds a tensor-product B-spline surface. |
| `func (s *BSplineSurface) DegreeU() int` | DegreeU returns the surface degree in the u direction. |
| `func (s *BSplineSurface) DegreeV() int` | DegreeV returns the surface degree in the v direction. |
| `func (s *BSplineSurface) Dim() int` | Dim returns the spatial dimension of the control net. |
| `func (s *BSplineSurface) Domain() (uLo, uHi, vLo, vHi float64)` | Domain returns the valid parameter rectangle [uLo,uHi] x [vLo,vHi]. |
| `func (s *BSplineSurface) Eval(u, v float64) Vec` | Eval returns the surface point at (u,v) using the tensor product of the univariate B-spline basis functions in each direction. |

</details>

<details>
<summary><code>BezierCurve</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewBezierCurve(ctrl ...Vec) (*BezierCurve, error)` | NewBezierCurve returns a Bezier curve with the given control points. |
| `func (b *BezierCurve) ControlPoints() []Vec` | ControlPoints returns an independent copy of the control points. |
| `func (b *BezierCurve) Degree() int` | Degree returns the polynomial degree of the curve (number of control points minus one). |
| `func (b *BezierCurve) Dim() int` | Dim returns the spatial dimension of the curve's control points. |
| `func (b *BezierCurve) Domain() (lo, hi float64)` | Domain returns the parameter domain of a Bezier curve, always [0,1]. |
| `func (b *BezierCurve) Elevate() *BezierCurve` | Elevate returns an equivalent Bezier curve of degree one higher. |
| `func (b *BezierCurve) Eval(t float64) Vec` | Eval returns the point on the curve at parameter t in [0,1]. |
| `func (b *BezierCurve) EvalDerivative(t float64) Vec` | EvalDerivative returns the first derivative vector of the curve at parameter t. |
| `func (b *BezierCurve) EvalDerivativeK(k int, t float64) Vec` | EvalDerivativeK returns the k-th derivative vector of the curve at parameter t (k >= 0). |
| `func (b *BezierCurve) Hodograph() *BezierCurve` | Hodograph returns the derivative curve of b as a Bezier curve of degree one lower, whose control points are Q(i) = n*(P(i+1)-P(i)). |
| `func (b *BezierCurve) Length() float64` | Length returns an arc-length estimate of the curve obtained by adaptive Gauss-Legendre quadrature of the speed \|C'(t)\| over [0,1]. |
| `func (b *BezierCurve) Polyline(n int) []Vec` | Polyline samples the curve at n+1 evenly spaced parameter values in [0,1] and returns the resulting points, a convenient coarse polygonal… |
| `func (b *BezierCurve) Reverse() *BezierCurve` | Reverse returns the curve traversed in the opposite direction (control points reversed), so Reverse().Eval(t) == Eval(1-t). |
| `func (b *BezierCurve) Split(t float64) (left, right *BezierCurve)` | Split subdivides the curve at parameter t and returns the two Bezier curves that together reproduce it: left covers the parameter range [0,t] and… |
| `func (b *BezierCurve) Tangent(t float64) Vec` | Tangent returns the unit tangent vector of the curve at parameter t. |

</details>

<details>
<summary><code>BezierSurface</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewBezierSurface(ctrl [][]Vec) (*BezierSurface, error)` | NewBezierSurface builds a tensor-product Bezier surface from a rectangular control net ctrl`i``j`, where i runs over the u direction and j over the v… |
| `func (s *BezierSurface) DegreeU() int` | DegreeU returns the surface degree in the u direction. |
| `func (s *BezierSurface) DegreeV() int` | DegreeV returns the surface degree in the v direction. |
| `func (s *BezierSurface) Dim() int` | Dim returns the spatial dimension of the control net. |
| `func (s *BezierSurface) Eval(u, v float64) Vec` | Eval returns the surface point at parameters (u,v) in [0,1]x[0,1] using the Bernstein basis in both directions. |
| `func (s *BezierSurface) EvalPartialU(u, v float64) Vec` | EvalPartialU returns the partial derivative of the surface with respect to u at (u,v). |
| `func (s *BezierSurface) EvalPartialV(u, v float64) Vec` | EvalPartialV returns the partial derivative of the surface with respect to v at (u,v). |
| `func (s *BezierSurface) Normal(u, v float64) (Vec, error)` | Normal returns the (unnormalised) surface normal at (u,v) as the cross product of the u- and v-partial derivatives. |

</details>

<details>
<summary><code>CatmullRomCurve</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewCatmullRomCurve(pts []Vec, alpha float64) (*CatmullRomCurve, error)` | NewCatmullRomCurve builds a Catmull-Rom curve through pts using the given alpha (0 uniform, 0.5 centripetal, 1 chordal). |
| `func (c *CatmullRomCurve) Domain() (lo, hi float64)` | Domain returns the parameter interval [0, n-1] of the curve. |
| `func (c *CatmullRomCurve) Eval(u float64) Vec` | Eval returns the point on the curve at parameter u in [0, n-1]. |
| `func (c *CatmullRomCurve) EvalDerivative(u float64) Vec` | EvalDerivative returns the first derivative of the curve at parameter u, approximated by a symmetric finite difference (exact to rounding for the… |
| `func (c *CatmullRomCurve) Point(k int, t float64) Vec` | Point returns the point on segment k (between control points k and k+1) at local parameter t in [0,1]. |

</details>

<details>
<summary><code>CubicSpline</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func ClampedCubicSpline(x, y []float64, dyLeft, dyRight float64) (*CubicSpline, error)` | ClampedCubicSpline builds a cubic spline whose first derivative equals dyLeft at the first knot and dyRight at the last knot. |
| `func NaturalCubicSpline(x, y []float64) (*CubicSpline, error)` | NaturalCubicSpline builds a natural cubic spline (zero second derivative at both ends) through (x`i`, y`i`). |
| `func NewCubicSpline(x, y []float64, bt BoundaryType, dyLeft, dyRight float64) (*CubicSpline, error)` | NewCubicSpline builds a cubic spline through the samples (x`i`, y`i`) with the requested boundary type. |
| `func NotAKnotCubicSpline(x, y []float64) (*CubicSpline, error)` | NotAKnotCubicSpline builds a cubic spline with not-a-knot end conditions, the default used by many numerical packages. |
| `func PeriodicCubicSpline(x, y []float64) (*CubicSpline, error)` | PeriodicCubicSpline builds a periodic cubic spline; the first and last samples must have equal ordinates (y[0] == y`len-1`). |
| `func (cs *CubicSpline) Boundary() BoundaryType` | Boundary returns the boundary type used to build the spline. |
| `func (cs *CubicSpline) Domain() (lo, hi float64)` | Domain returns the first and last knot abscissae. |
| `func (cs *CubicSpline) Eval(q float64) float64` | Eval returns the spline value at q. |
| `func (cs *CubicSpline) EvalDerivative(q float64) float64` | EvalDerivative returns the first derivative of the spline at q. |
| `func (cs *CubicSpline) EvalSecondDerivative(q float64) float64` | EvalSecondDerivative returns the second derivative of the spline at q. |
| `func (cs *CubicSpline) EvalThirdDerivative(q float64) float64` | EvalThirdDerivative returns the (piecewise constant) third derivative of the spline at q. |
| `func (cs *CubicSpline) Integrate(a, b float64) float64` | Integrate returns the definite integral of the spline from a to b. |
| `func (cs *CubicSpline) Knots() []float64` | Knots returns a copy of the spline's knot abscissae. |
| `func (cs *CubicSpline) Len() int` | Len returns the number of polynomial segments in the spline. |
| `func (cs *CubicSpline) SegmentCoeffs(i int) [4]float64` | SegmentCoeffs returns the local monomial coefficients [c0,c1,c2,c3] of segment i, so that on that segment the spline equals c0 + c1*s + c2*s^2 +… |

</details>

<details>
<summary><code>NURBSCurve</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NURBSCircle(center Vec, radius float64) (*NURBSCurve, error)` | NURBSCircle returns a full circle of the given radius centred at center, represented exactly as a quadratic NURBS curve with nine control points in… |
| `func NewNURBSCurve(ctrl []Vec, weights, knots []float64, degree int) (*NURBSCurve, error)` | NewNURBSCurve builds a NURBS curve from control points, matching positive weights, a non-decreasing knot vector and a degree. |
| `func (c *NURBSCurve) ControlPoints() []Vec` | ControlPoints returns an independent copy of the control points. |
| `func (c *NURBSCurve) Degree() int` | Degree returns the degree of the curve. |
| `func (c *NURBSCurve) Dim() int` | Dim returns the spatial dimension of the control points. |
| `func (c *NURBSCurve) Domain() (lo, hi float64)` | Domain returns the valid parameter interval of the curve. |
| `func (c *NURBSCurve) Eval(u float64) Vec` | Eval returns the point on the NURBS curve at parameter u by evaluating the weighted basis and dividing by the weight sum. |
| `func (c *NURBSCurve) EvalDerivative(u float64) Vec` | EvalDerivative returns the first derivative vector of the NURBS curve at parameter u. |
| `func (c *NURBSCurve) EvalDerivatives(u float64, d int) []Vec` | EvalDerivatives returns the point on the curve and its derivatives up to order d as a slice of length d+1, applying the rational-curve quotient rule… |
| `func (c *NURBSCurve) Knots() []float64` | Knots returns an independent copy of the knot vector. |
| `func (c *NURBSCurve) Weights() []float64` | Weights returns an independent copy of the control-point weights. |

</details>

<details>
<summary><code>NURBSSurface</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewNURBSSurface(ctrl [][]Vec, weights [][]float64, knotsU, knotsV []float64, degU, degV int) (*NURBSSurface, error)` | NewNURBSSurface builds a tensor-product NURBS surface. |
| `func (s *NURBSSurface) DegreeU() int` | DegreeU returns the surface degree in the u direction. |
| `func (s *NURBSSurface) DegreeV() int` | DegreeV returns the surface degree in the v direction. |
| `func (s *NURBSSurface) Dim() int` | Dim returns the spatial dimension of the control net. |
| `func (s *NURBSSurface) Domain() (uLo, uHi, vLo, vHi float64)` | Domain returns the valid parameter rectangle of the surface. |
| `func (s *NURBSSurface) Eval(u, v float64) Vec` | Eval returns the surface point at (u,v) using the rational tensor-product basis (weighted sum divided by the weight sum). |

</details>

<details>
<summary><code>PiecewiseCubic</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func AkimaSpline(x, y []float64) (*PiecewiseCubic, error)` | AkimaSpline builds an Akima cubic interpolant through (x`i`, y`i`). |
| `func CardinalSpline(x, y []float64, tension float64) (*PiecewiseCubic, error)` | CardinalSpline builds a cardinal spline through (x`i`, y`i`) with the given tension in [0,1]; tension 0 reproduces the Catmull-Rom spline while… |
| `func CatmullRomSpline(x, y []float64) (*PiecewiseCubic, error)` | CatmullRomSpline builds the Catmull-Rom interpolant through (x`i`, y`i`); it is the cardinal spline with zero tension. |
| `func HermiteSpline(x, y, m []float64) (*PiecewiseCubic, error)` | HermiteSpline builds a piecewise cubic Hermite interpolant from nodal values y and nodal first derivatives m at the abscissae x. |
| `func MonotoneCubicSpline(x, y []float64) (*PiecewiseCubic, error)` | MonotoneCubicSpline builds a shape-preserving monotone cubic interpolant using the Fritsch-Carlson method (the PCHIP algorithm). |
| `func (p *PiecewiseCubic) Domain() (lo, hi float64)` | Domain returns the first and last abscissae of the interpolant. |
| `func (p *PiecewiseCubic) Eval(q float64) float64` | Eval returns the interpolant value at q, extrapolating with the nearest end piece for queries outside the knot range. |
| `func (p *PiecewiseCubic) EvalDerivative(q float64) float64` | EvalDerivative returns the first derivative of the interpolant at q. |
| `func (p *PiecewiseCubic) EvalSecondDerivative(q float64) float64` | EvalSecondDerivative returns the second derivative of the interpolant at q. |
| `func (p *PiecewiseCubic) Integrate(a, b float64) float64` | Integrate returns the definite integral of the interpolant from a to b. |
| `func (p *PiecewiseCubic) Knots() []float64` | Knots returns a copy of the interpolant's abscissae. |
| `func (p *PiecewiseCubic) Len() int` | Len returns the number of cubic pieces. |
| `func (p *PiecewiseCubic) SegmentCoeffs(i int) [4]float64` | SegmentCoeffs returns the local monomial coefficients of piece i. |

</details>

<details>
<summary><code>Vec</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func CatmullRomPoint(p0, p1, p2, p3 Vec, t float64) Vec` | CatmullRomPoint evaluates the uniform Catmull-Rom segment defined by four successive control points at local parameter t in [0,1]; the returned point… |
| `func CatmullRomSegment(p0, p1, p2, p3 Vec, alpha, t float64) Vec` | CatmullRomSegment evaluates a Catmull-Rom segment between p1 and p2 using the Barry-Goldman non-uniform formulation with parameterisation exponent… |
| `func CloneVecs(pts []Vec) []Vec` | CloneVecs returns an independent deep copy of a slice of vectors. |
| `func Cross(a, b Vec) Vec` | Cross returns the 3-D cross product a x b. |
| `func DeBoor(p int, U []float64, ctrl []Vec, u float64) Vec` | DeBoor evaluates a B-spline curve of degree p at parameter u directly with de Boor's algorithm, given the control points and knot vector. |
| `func DeCasteljau(ctrl []Vec, t float64) Vec` | DeCasteljau evaluates the Bezier curve with the given control points at parameter t in [0,1] using de Casteljau's algorithm, the numerically stable… |
| `func HermiteCurvePoint(p0, p1, t0, t1 Vec, t float64) Vec` | HermiteCurvePoint evaluates a vector-valued cubic Hermite segment with endpoints p0, p1 and endpoint tangents (derivatives) t0, t1 at local parameter… |
| `func NewVec(coords ...float64) Vec` | NewVec returns a Vec containing a copy of the supplied coordinates. |
| `func SampleCurve(c Curve, n int) []Vec` | SampleCurve returns n+1 points sampled at uniformly spaced parameter values across the curve's domain (n >= 1). |
| `func VecCombine(coeff []float64, pts []Vec) (Vec, error)` | VecCombine returns the linear combination sum_i coeff`i`*pts`i`. |
| `func VecLerp(a, b Vec, t float64) Vec` | VecLerp returns the linear interpolation (1-t)*a + t*b of two vectors. |
| `func ZeroVec(n int) Vec` | ZeroVec returns the zero vector of dimension n. |
| `func (v Vec) Add(w Vec) Vec` | Add returns the component-wise sum v+w. |
| `func (v Vec) AddScaled(s float64, w Vec) Vec` | AddScaled returns v + s*w, a fused multiply-add. |
| `func (v Vec) Clone() Vec` | Clone returns an independent copy of v. |
| `func (v Vec) Dim() int` | Dim reports the number of coordinates in v. |
| `func (v Vec) Dist(w Vec) float64` | Dist returns the Euclidean distance between v and w. |
| `func (v Vec) Dist2(w Vec) float64` | Dist2 returns the squared Euclidean distance between v and w. |
| `func (v Vec) Dot(w Vec) float64` | Dot returns the Euclidean inner product of v and w. |
| `func (v Vec) Equal(w Vec, tol float64) bool` | Equal reports whether v and w are equal within absolute tolerance tol. |
| `func (v Vec) Lerp(w Vec, t float64) Vec` | Lerp returns the linear interpolation (1-t)*v + t*w. |
| `func (v Vec) Neg() Vec` | Neg returns -v. |
| `func (v Vec) Norm() float64` | Norm returns the Euclidean norm (length) of v. |
| `func (v Vec) Norm2() float64` | Norm2 returns the squared Euclidean norm of v. |
| `func (v Vec) Scale(s float64) Vec` | Scale returns v scaled component-wise by s. |
| `func (v Vec) Sub(w Vec) Vec` | Sub returns the component-wise difference v-w. |
| `func (v Vec) Unit() Vec` | Unit returns v scaled to unit length. |

</details>

### Variables

`ErrTooFewPoints`, `ErrLenMismatch`, `ErrNotIncreasing`, `ErrDim`, `ErrDegree`, `ErrKnots`, `ErrWeights`, `ErrParam`, `ErrEmpty`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/splines).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
