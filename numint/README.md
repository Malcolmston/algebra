# numint — implements numerical integration (quadrature) of functions of one and several variables

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/numint.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/numint)

Package numint implements numerical integration (quadrature) of functions of
one and several variables using only the Go standard library.

The package collects the classical fixed-order Newton-Cotes rules (trapezoid,
midpoint, Simpson 1/3, Simpson 3/8, Boole, Weddle), their composite forms, and
higher-accuracy schemes: Romberg extrapolation, adaptive Simpson and adaptive
Gauss-Kronrod, fixed Gauss-Legendre rules from 2 to 10 points plus an n-point
routine for arbitrary order, Gauss-Lobatto rules, Clenshaw-Curtis and
Gauss-Chebyshev quadrature, and a tanh-sinh (double-exponential) rule for
integrands with endpoint singularities. Sample-based helpers integrate
tabulated data, including cumulative integrals, and tensor-product routines
extend the fixed rules to two and three dimensions. A small deterministic
Monte-Carlo family covers higher-dimensional volumes.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/numint"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func AdaptiveSimpson(f Func, a, b, tol float64) float64` | AdaptiveSimpson approximates the integral of f over [a, b] using adaptive Simpson quadrature with an absolute tolerance tol. |
| `func AdaptiveSimpsonBudget(f Func, a, b, tol float64, maxDepth, maxEvals int) (value float64, converged bool)` | AdaptiveSimpsonBudget is `AdaptiveSimpson` with explicit limits on the recursion depth and the total number of integrand evaluations. |
| `func AdaptiveSimpsonMaxDepth(f Func, a, b, tol float64, maxDepth int) float64` | AdaptiveSimpsonMaxDepth approximates the integral of f over [a, b] using adaptive Simpson quadrature with an absolute tolerance tol and a bound… |
| `func AdaptiveTrapezoid(f Func, a, b, tol float64) float64` | AdaptiveTrapezoid approximates the integral of f over [a, b] by recursively bisecting subintervals until the trapezoid estimate on each meets the… |
| `func AdaptiveTrapezoidBudget(f Func, a, b, tol float64, maxDepth, maxEvals int) (value float64, converged bool)` | AdaptiveTrapezoidBudget is `AdaptiveTrapezoid` with explicit limits on the recursion depth and the total number of integrand evaluations. |
| `func ArcLength2D(p Path2, a, b float64, n int) float64` | ArcLength2D approximates the arc length of the plane curve p over the parameter interval [a, b] using n Simpson panels. |
| `func ArcLength3D(p Path3, a, b float64, n int) float64` | ArcLength3D approximates the arc length of the space curve p over the parameter interval [a, b] using n Simpson panels. |
| `func ArcLengthPolyline2D(xs, ys []float64) float64` | ArcLengthPolyline2D returns the total length of the polyline through the given vertices. |
| `func AreaBetween(ylo, yhi Func, a, b float64, n int) float64` | AreaBetween returns the area of the region enclosed between the curves y = ylo(x) and y = yhi(x) for x in [a, b]. |
| `func Boole(f Func, a, b float64, n int) float64` | Boole approximates the integral of f over [a, b] using the composite Boole rule; n is the number of panels and is rounded up to a multiple of four. |
| `func BooleRule(f Func, a, b float64) float64` | BooleRule applies the five-point Boole rule to f over [a, b]. |
| `func CentroidRegion(ylo, yhi Func, a, b float64, n int) (xbar, ybar float64)` | CentroidRegion returns the centroid (xbar, ybar) of the region between the curves y = ylo(x) and y = yhi(x) over [a, b]. |
| `func Circulation2D(field Field2, p Path2, a, b float64, n int) float64` | Circulation2D approximates the circulation of the vector field around the closed plane curve p, i.e. |
| `func ClenshawCurtis(f Func, a, b float64, n int) float64` | ClenshawCurtis approximates the integral of f over [a, b] using Clenshaw-Curtis quadrature with n+1 Chebyshev extreme points; n must be a positive… |
| `func CompositeGaussLegendre(f Func, a, b float64, panels, n int) float64` | CompositeGaussLegendre approximates the integral of f over [a, b] by splitting the interval into panels equal subintervals and applying the n-point… |
| `func CumulativeSimpson(xs, ys []float64) []float64` | CumulativeSimpson returns the cumulative integral of tabulated data (xs, ys). |
| `func CumulativeTrapezoid(xs, ys []float64) []float64` | CumulativeTrapezoid returns the cumulative integral of tabulated data (xs, ys) evaluated at each xs. |
| `func CumulativeTrapezoidUniform(ys []float64, h float64) []float64` | CumulativeTrapezoidUniform returns the cumulative integral of equally spaced samples ys with spacing h. |
| `func DoubleAverage(f Func2, ax, bx, ay, by float64) float64` | DoubleAverage returns the mean value of f over the rectangle [ax, bx] x [ay, by], i.e. |
| `func DoubleBoole(f Func2, ax, bx, ay, by float64, nx, ny int) float64` | DoubleBoole approximates the double integral of f over the rectangle [ax, bx] x [ay, by] using the composite Boole rule. |
| `func DoubleExponential(f Func, a, b float64, n int) float64` | DoubleExponential is a synonym for TanhSinh, using the alternative name for the double-exponential quadrature rule. |
| `func DoubleGaussLegendre(f Func2, ax, bx, ay, by float64, n int) float64` | DoubleGaussLegendre approximates the double integral of f over the rectangle [ax, bx] x [ay, by] using an n-by-n tensor Gauss-Legendre rule. |
| `func DoubleIntegral(f Func2, ax, bx, ay, by float64) float64` | DoubleIntegral approximates the double integral of f over the rectangle [ax, bx] x [ay, by] using a fine composite Simpson rule with sensible default… |
| `func DoubleIntegralRegion(f Func2, a, b float64, ylo, yhi Func, nx, ny int) float64` | DoubleIntegralRegion approximates the integral of f over the region { (x, y) : a <= x <= b, ylo(x) <= y <= yhi(x) } using the composite Simpson rule… |
| `func DoubleIntegralRegionGL(f Func2, a, b float64, ylo, yhi Func, n int) float64` | DoubleIntegralRegionGL approximates the integral of f over the region { (x, y) : a <= x <= b, ylo(x) <= y <= yhi(x) } using an n-point Gauss-Legendre… |
| `func DoubleMidpoint(f Func2, ax, bx, ay, by float64, nx, ny int) float64` | DoubleMidpoint approximates the double integral of f over the rectangle [ax, bx] x [ay, by] using the composite midpoint rule. |
| `func DoubleSimpson(f Func2, ax, bx, ay, by float64, nx, ny int) float64` | DoubleSimpson approximates the double integral of f over the rectangle [ax, bx] x [ay, by] using the composite Simpson rule. |
| `func DoubleTrapezoid(f Func2, ax, bx, ay, by float64, nx, ny int) float64` | DoubleTrapezoid approximates the double integral of f over the rectangle [ax, bx] x [ay, by] using the composite trapezoidal rule with nx and ny… |
| `func ExpectationNormal(f Func, mean, stddev float64, n int) float64` | ExpectationNormal approximates E[f(X)] for a normal random variable X with the given mean and standard deviation, using n-point Gauss-Hermite… |
| `func Flux2D(field Field2, p Path2, a, b float64, n int) float64` | Flux2D approximates the outward flux of the vector field across the plane curve p over [a, b], i.e. |
| `func GaussChebyshev(f Func, n int) float64` | GaussChebyshev approximates the weighted integral of f over [-1, 1] with the Chebyshev weight 1/sqrt(1-x^2), that is the integral of… |
| `func GaussHermite(f Func, n int) float64` | GaussHermite approximates the integral of exp(-x^2) * f(x) over (-inf, inf) using n-point Gauss-Hermite quadrature. |
| `func GaussHermiteNodes(n int) (nodes, weights []float64)` | GaussHermiteNodes returns the n nodes and weights of Gauss-Hermite quadrature for the weight exp(-x^2) on (-inf, inf). |
| `func GaussKronrod15(f Func, a, b float64) (value, errEst float64)` | GaussKronrod15 approximates the integral of f over [a, b] using the 15-point Gauss-Kronrod rule. |
| `func GaussKronrod15Nodes() (abscissae, kronrodWeights, gaussWeights []float64)` | GaussKronrod15Nodes returns the positive-half abscissae of the 15-point Gauss-Kronrod rule on [-1, 1], the corresponding Kronrod weights, and the… |
| `func GaussLaguerre(f Func, n int) float64` | GaussLaguerre approximates the integral of exp(-x) * f(x) over (0, inf) using n-point Gauss-Laguerre quadrature. |
| `func GaussLaguerreGen(f Func, alpha float64, n int) float64` | GaussLaguerreGen approximates the integral of x^alpha * exp(-x) * f(x) over (0, inf) using n-point generalized Gauss-Laguerre quadrature. |
| `func GaussLaguerreGenNodes(n int, alpha float64) (nodes, weights []float64)` | GaussLaguerreGenNodes returns the n nodes and weights of generalized Gauss-Laguerre quadrature for the weight x^alpha * exp(-x) on (0, inf). |
| `func GaussLaguerreNodes(n int) (nodes, weights []float64)` | GaussLaguerreNodes returns the n nodes and weights of Gauss-Laguerre quadrature for the weight exp(-x) on (0, inf). |
| `func GaussLegendre10(f Func, a, b float64) float64` | GaussLegendre10 approximates the integral of f over [a, b] using the ten-point Gauss-Legendre rule (exact for polynomials up to degree 19). |
| `func GaussLegendre2(f Func, a, b float64) float64` | GaussLegendre2 approximates the integral of f over [a, b] using the two-point Gauss-Legendre rule (exact for polynomials up to degree 3). |
| `func GaussLegendre2D(f Func2, ax, bx, ay, by float64, n int) float64` | GaussLegendre2D approximates the double integral of f over the rectangle [ax, bx] x [ay, by] using the tensor product of the n-point Gauss-Legendre… |
| `func GaussLegendre3(f Func, a, b float64) float64` | GaussLegendre3 approximates the integral of f over [a, b] using the three-point Gauss-Legendre rule (exact for polynomials up to degree 5). |
| `func GaussLegendre3D(f Func3, ax, bx, ay, by, az, bz float64, n int) float64` | GaussLegendre3D approximates the triple integral of f over the box [ax, bx] x [ay, by] x [az, bz] using the tensor product of the n-point… |
| `func GaussLegendre4(f Func, a, b float64) float64` | GaussLegendre4 approximates the integral of f over [a, b] using the four-point Gauss-Legendre rule (exact for polynomials up to degree 7). |
| `func GaussLegendre5(f Func, a, b float64) float64` | GaussLegendre5 approximates the integral of f over [a, b] using the five-point Gauss-Legendre rule (exact for polynomials up to degree 9). |
| `func GaussLegendre6(f Func, a, b float64) float64` | GaussLegendre6 approximates the integral of f over [a, b] using the six-point Gauss-Legendre rule (exact for polynomials up to degree 11). |
| `func GaussLegendre7(f Func, a, b float64) float64` | GaussLegendre7 approximates the integral of f over [a, b] using the seven-point Gauss-Legendre rule (exact for polynomials up to degree 13). |
| `func GaussLegendre8(f Func, a, b float64) float64` | GaussLegendre8 approximates the integral of f over [a, b] using the eight-point Gauss-Legendre rule (exact for polynomials up to degree 15). |
| `func GaussLegendre9(f Func, a, b float64) float64` | GaussLegendre9 approximates the integral of f over [a, b] using the nine-point Gauss-Legendre rule (exact for polynomials up to degree 17). |
| `func GaussLegendreN(f Func, a, b float64, n int) float64` | GaussLegendreN approximates the integral of f over [a, b] using the n-point Gauss-Legendre rule for any n >= 1. |
| `func GaussLegendreNodes(n int) (nodes, weights []float64)` | GaussLegendreNodes returns the abscissae and weights of the n-point Gauss-Legendre rule on the reference interval [-1, 1]. |
| `func GaussLobatto3(f Func, a, b float64) float64` | GaussLobatto3 approximates the integral of f over [a, b] using the three-point Gauss-Lobatto rule (Simpson's rule), which includes both endpoints. |
| `func GaussLobatto4(f Func, a, b float64) float64` | GaussLobatto4 approximates the integral of f over [a, b] using the four-point Gauss-Lobatto rule, which includes both endpoints and the interior… |
| `func GaussLobatto5(f Func, a, b float64) float64` | GaussLobatto5 approximates the integral of f over [a, b] using the five-point Gauss-Lobatto rule, which includes both endpoints, the centre, and the… |
| `func HaltonPoint(index int, bases []int) []float64` | HaltonPoint returns the Halton point at the given index using one base per coordinate. |
| `func HaltonPrimes(k int) []int` | HaltonPrimes returns the first k prime numbers, suitable as Halton bases for a k-dimensional sequence. |
| `func HaltonSequence(count int, bases []int) [][]float64` | HaltonSequence returns count Halton points (indices 1..count) using the given per-coordinate bases. |
| `func Integrate(f Func, a, b float64) float64` | Integrate approximates the integral of f over [a, b] to roughly full double-precision accuracy using globally adaptive Gauss-Kronrod quadrature. |
| `func Integrate2D(f Func2, ax, bx, ay, by float64) float64` | Integrate2D approximates the double integral of f over the rectangle [ax, bx] x [ay, by] to high accuracy using an 8-point Gauss-Legendre tensor… |
| `func IntegrateAnnulus(f Func2, cx, cy, rin, rout float64, nr, nt int) float64` | IntegrateAnnulus approximates the integral of f over the annulus centered at (cx, cy) with inner radius rin and outer radius rout. |
| `func IntegrateBall(f Func3, cx, cy, cz, r float64, nr, nphi, ntheta int) float64` | IntegrateBall approximates the volume integral of f over the solid ball of radius r centered at (cx, cy, cz) using spherical coordinates. |
| `func IntegrateCylinder(f Func3, cx, cy, r, z0, z1 float64, nr, ntheta, nz int) float64` | IntegrateCylinder approximates the volume integral of f over the cylinder of radius r centered on the axis x = cx, y = cy and spanning z in [z0, z1],… |
| `func IntegrateDisk(f Func2, cx, cy, r float64, nr, nt int) float64` | IntegrateDisk approximates the integral of f over the disk of radius r centered at (cx, cy), using polar coordinates with nr radial and nt angular… |
| `func IntegrateExpWeighted(f Func, rate float64, n int) float64` | IntegrateExpWeighted approximates the integral of exp(-rate*x) * f(x) over (0, inf) using n-point Gauss-Laguerre quadrature; rate must be positive. |
| `func IntegrateGaussianWeighted(f Func, n int) float64` | IntegrateGaussianWeighted approximates the integral of exp(-x^2/2) * f(x) over (-inf, inf) using n-point Gauss-Hermite quadrature. |
| `func IntegrateInfinite(f Func, n int) float64` | IntegrateInfinite approximates the integral of f over the whole real line using the substitution x = t/(1-t^2) with n-point Gauss-Legendre quadrature. |
| `func IntegrateLowerInfinite(f Func, b float64, n int) float64` | IntegrateLowerInfinite approximates the integral of f over (-inf, b] using the substitution x = b - u/(1-u) with n-point Gauss-Legendre quadrature. |
| `func IntegrateND(f FuncND, lower, upper []float64, n int) float64` | IntegrateND approximates the integral of f over the axis-aligned box defined by the lower and upper bounds using an n-point Gauss-Legendre tensor… |
| `func IntegratePolar(g Func2, r0, r1, t0, t1 float64, nr, nt int) float64` | IntegratePolar approximates the integral of g(r, theta) * r over the polar rectangle r in [r0, r1], theta in [t0, t1]; the Jacobian factor r is… |
| `func IntegrateSingularEndpoint(f Func, a, b float64, n int) float64` | IntegrateSingularEndpoint approximates the integral of f over [a, b] when f has an integrable singularity at the lower endpoint a, using the… |
| `func IntegrateSphereSurface(f Func3, cx, cy, cz, r float64, nphi, ntheta int) float64` | IntegrateSphereSurface approximates the surface integral of f over the sphere of radius r centered at (cx, cy, cz), with nphi polar and ntheta… |
| `func IntegrateTriangle(f Func2, x1, y1, x2, y2, x3, y3 float64, n int) float64` | IntegrateTriangle approximates the integral of f over the triangle with the given vertices using an n-point Gauss-Legendre rule under a Duffy… |
| `func IntegrateUpperInfinite(f Func, a float64, n int) float64` | IntegrateUpperInfinite approximates the integral of f over [a, inf) using the substitution x = a + t/(1-t) with n-point Gauss-Legendre quadrature. |
| `func LaplaceTransform(f Func, s float64, n int) float64` | LaplaceTransform approximates the Laplace transform of f evaluated at s, that is the integral of exp(-s*t) * f(t) over (0, inf), using n-point… |
| `func LineIntegralPolyline2D(f Func2, xs, ys []float64) float64` | LineIntegralPolyline2D approximates the scalar line integral of f along the polyline through the given vertices, using the midpoint of each segment. |
| `func LineIntegralScalar2D(f Func2, p Path2, a, b float64, n int) float64` | LineIntegralScalar2D approximates the scalar line integral of f along the plane curve p over [a, b], i.e. |
| `func LineIntegralScalar3D(f Func3, p Path3, a, b float64, n int) float64` | LineIntegralScalar3D approximates the scalar line integral of f along the space curve p over [a, b], i.e. |
| `func LineIntegralVector2D(field Field2, p Path2, a, b float64, n int) float64` | LineIntegralVector2D approximates the line integral of the vector field along the plane curve p over [a, b], i.e. |
| `func LineIntegralVector3D(field Field3, p Path3, a, b float64, n int) float64` | LineIntegralVector3D approximates the line integral of the vector field along the space curve p over [a, b], i.e. |
| `func Midpoint(f Func, a, b float64, n int) float64` | Midpoint approximates the integral of f over [a, b] using the composite midpoint rule with n panels. |
| `func Midpoint2D(f Func2, ax, bx, ay, by float64, nx, ny int) float64` | Midpoint2D approximates the double integral of f over the rectangle [ax, bx] x [ay, by] using the composite midpoint rule with nx panels in x and ny… |
| `func MidpointRule(f Func, a, b float64) float64` | MidpointRule applies the one-point open midpoint rule to f over [a, b]. |
| `func MidpointSamples(ys []float64, h float64) float64` | MidpointSamples integrates samples ys taken at the midpoints of n panels of width h using the composite midpoint rule. |
| `func MonteCarlo(f Func, a, b float64, n int, seed int64) float64` | MonteCarlo approximates the integral of f over [a, b] using n uniform random samples drawn from a generator seeded with seed, making the result… |
| `func MonteCarlo2D(f Func2, ax, bx, ay, by float64, n int, seed int64) float64` | MonteCarlo2D approximates the double integral of f over the rectangle [ax, bx] x [ay, by] using n uniform random samples from a generator seeded with… |
| `func MonteCarlo3D(f Func3, ax, bx, ay, by, az, bz float64, n int, seed int64) float64` | MonteCarlo3D estimates the triple integral of f over the box [ax, bx] x [ay, by] x [az, bz] by averaging n uniform random samples drawn from the… |
| `func MonteCarloAntithetic(f Func, a, b float64, n int, seed int64) float64` | MonteCarloAntithetic estimates the integral of f over [a, b] using n antithetic sample pairs, which reduces variance for monotone integrands. |
| `func MonteCarloBall(f Func3, cx, cy, cz, r float64, n int, seed int64) float64` | MonteCarloBall estimates the integral of f over the solid ball of radius r centered at (cx, cy, cz) by averaging n uniform samples over the ball. |
| `func MonteCarloDisk(f Func2, cx, cy, r float64, n int, seed int64) float64` | MonteCarloDisk estimates the integral of f over the disk of radius r centered at (cx, cy) by averaging n samples drawn uniformly over the disk. |
| `func MonteCarloImportance(f Func, sample func(rng *rand.Rand) float64, pdf Func, n int, seed int64) float64` | MonteCarloImportance estimates the integral of f using importance sampling: sample draws a point from a density proportional to pdf, and the… |
| `func MonteCarloND(f FuncND, lower, upper []float64, n int, seed int64) float64` | MonteCarloND approximates the integral of f over the axis-aligned box defined by lower and upper using n uniform random samples from a generator… |
| `func MonteCarloRegion(f Func2, a, b float64, ylo, yhi Func, n int, seed int64) float64` | MonteCarloRegion estimates the integral of f over the region { (x, y) : a <= x <= b, ylo(x) <= y <= yhi(x) } using n random samples. |
| `func MonteCarloStratified(f Func, a, b float64, strata, perStratum int, seed int64) float64` | MonteCarloStratified estimates the integral of f over [a, b] using stratified sampling with the given number of equal strata and samples per stratum. |
| `func MonteCarloVolume(pred func(x []float64) bool, lower, upper []float64, n int, seed int64) float64` | MonteCarloVolume estimates the volume of the subset of the box [lower, upper] for which pred reports true, by sampling n uniform random points. |
| `func QuasiMonteCarlo(f Func, a, b float64, n int) float64` | QuasiMonteCarlo estimates the integral of f over [a, b] using n points of the base-2 Halton (van der Corput) sequence. |
| `func QuasiMonteCarlo2D(f Func2, ax, bx, ay, by float64, n int) float64` | QuasiMonteCarlo2D estimates the double integral of f over the rectangle [ax, bx] x [ay, by] using n points of the Halton sequence with bases 2 and 3. |
| `func QuasiMonteCarloND(f FuncND, lower, upper []float64, n int) float64` | QuasiMonteCarloND estimates the integral of f over the box [lower, upper] using n Halton points whose bases are the first len(lower) primes. |
| `func RichardsonExtrapolate(coarse, fine, order float64) float64` | RichardsonExtrapolate combines a coarse and a fine estimate of the same quantity, computed with step sizes differing by a factor of two, into a… |
| `func Romberg(f Func, a, b float64, maxLevel int) float64` | Romberg approximates the integral of f over [a, b] using Romberg integration with maxLevel levels of Richardson extrapolation on the composite… |
| `func RombergTable(f Func, a, b float64, maxLevel int) (float64, [][]float64)` | RombergTable performs Romberg integration of f over [a, b] with maxLevel levels and returns the best estimate (the last diagonal entry) together with… |
| `func RombergTol(f Func, a, b, tol float64, maxLevel int) float64` | RombergTol approximates the integral of f over [a, b] using Romberg integration, stopping once two successive diagonal estimates differ by less than… |
| `func Simpson(f Func, a, b float64, n int) float64` | Simpson approximates the integral of f over [a, b] using the composite Simpson 1/3 rule; n is the number of panels and is rounded up to an even… |
| `func Simpson2D(f Func2, ax, bx, ay, by float64, nx, ny int) float64` | Simpson2D approximates the double integral of f over the rectangle [ax, bx] x [ay, by] using the composite Simpson rule with nx panels in x and ny… |
| `func Simpson38(f Func, a, b float64, n int) float64` | Simpson38 approximates the integral of f over [a, b] using the composite Simpson 3/8 rule; n is the number of panels and is rounded up to a multiple… |
| `func Simpson38Rule(f Func, a, b float64) float64` | Simpson38Rule applies the four-point Simpson 3/8 rule to f over [a, b]. |
| `func Simpson3D(f Func3, ax, bx, ay, by, az, bz float64, nx, ny, nz int) float64` | Simpson3D approximates the triple integral of f over the box [ax, bx] x [ay, by] x [az, bz] using the composite Simpson rule with the given number of… |
| `func SimpsonRule(f Func, a, b float64) float64` | SimpsonRule applies the three-point Simpson 1/3 rule to f over [a, b]. |
| `func SimpsonSamples(ys []float64, h float64) float64` | SimpsonSamples integrates equally spaced samples ys with spacing h using the composite Simpson 1/3 rule. |
| `func TanhSinh(f Func, a, b float64, n int) float64` | TanhSinh approximates the integral of f over [a, b] using the tanh-sinh (double-exponential) rule with 2n+1 abscissae per level. |
| `func TanhSinhTol(f Func, a, b, tol float64) float64` | TanhSinhTol approximates the integral of f over [a, b] using the tanh-sinh rule, refining the abscissa spacing until two successive level estimates… |
| `func Trapezoid(f Func, a, b float64, n int) float64` | Trapezoid approximates the integral of f over [a, b] using the composite trapezoid rule with n panels. |
| `func Trapezoid2D(f Func2, ax, bx, ay, by float64, nx, ny int) float64` | Trapezoid2D approximates the double integral of f over the rectangle [ax, bx] x [ay, by] using the composite trapezoid rule with nx panels in x and… |
| `func Trapezoid3D(f Func3, ax, bx, ay, by, az, bz float64, nx, ny, nz int) float64` | Trapezoid3D approximates the triple integral of f over the box [ax, bx] x [ay, by] x [az, bz] using the composite trapezoid rule with the given… |
| `func TrapezoidData(xs, ys []float64) float64` | TrapezoidData integrates tabulated data (xs, ys) with arbitrary, not necessarily uniform, spacing using the trapezoid rule. |
| `func TrapezoidRule(f Func, a, b float64) float64` | TrapezoidRule applies the two-point trapezoid rule to f over [a, b]. |
| `func TrapezoidSamples(ys []float64, h float64) float64` | TrapezoidSamples integrates equally spaced samples ys with spacing h using the composite trapezoid rule. |
| `func TripleAverage(f Func3, ax, bx, ay, by, az, bz float64) float64` | TripleAverage returns the mean value of f over the box, i.e. |
| `func TripleGaussLegendre(f Func3, ax, bx, ay, by, az, bz float64, n int) float64` | TripleGaussLegendre approximates the triple integral of f over the box using an n-by-n-by-n tensor Gauss-Legendre rule. |
| `func TripleIntegral(f Func3, ax, bx, ay, by, az, bz float64) float64` | TripleIntegral approximates the triple integral of f over the box using a composite Simpson rule with sensible default resolution. |
| `func TripleMidpoint(f Func3, ax, bx, ay, by, az, bz float64, nx, ny, nz int) float64` | TripleMidpoint approximates the triple integral of f over the box using the composite midpoint rule. |
| `func TripleSimpson(f Func3, ax, bx, ay, by, az, bz float64, nx, ny, nz int) float64` | TripleSimpson approximates the triple integral of f over the box using the composite Simpson rule. |
| `func TripleTrapezoid(f Func3, ax, bx, ay, by, az, bz float64, nx, ny, nz int) float64` | TripleTrapezoid approximates the triple integral of f over the box [ax, bx] x [ay, by] x [az, bz] using the composite trapezoidal rule. |
| `func VanDerCorput(index, base int) float64` | VanDerCorput returns the radical-inverse (van der Corput) value of index in the given base, a number in [0, 1). |
| `func Weddle(f Func, a, b float64, n int) float64` | Weddle approximates the integral of f over [a, b] using the composite Weddle rule; n is the number of panels and is rounded up to a multiple of six. |
| `func WeddleRule(f Func, a, b float64) float64` | WeddleRule applies the seven-point Weddle (Newton-Cotes degree six) rule to f over [a, b]. |

### Types

| Type | What it is |
| --- | --- |
| `Field2` | Field2 is a two-dimensional vector field (x, y) -> (u, v). |
| `Field3` | Field3 is a three-dimensional vector field (x, y, z) -> (u, v, w). |
| `Func` | Func is a real-valued function of a single real variable. |
| `Func2` | Func2 is a real-valued function of two real variables. |
| `Func3` | Func3 is a real-valued function of three real variables. |
| `FuncND` | FuncND is a real-valued function of a vector of real variables. |
| `Path2` | Path2 is a parametric plane curve t -> (x(t), y(t)). |
| `Path3` | Path3 is a parametric space curve t -> (x(t), y(t), z(t)). |
| `QuadResult` | QuadResult bundles the estimated value of an integral with an estimate of the absolute error and the number of integrand evaluations performed. |

<details>
<summary><code>QuadResult</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func AdaptiveGaussKronrod(f Func, a, b, tol float64) QuadResult` | AdaptiveGaussKronrod approximates the integral of f over [a, b] using globally adaptive 15-point Gauss-Kronrod quadrature with an absolute tolerance… |
| `func AdaptiveGaussKronrodBudget(f Func, a, b, tol float64, maxDepth, maxEvals int) QuadResult` | AdaptiveGaussKronrodBudget is `AdaptiveGaussKronrod` with explicit limits on the recursion depth and the total number of integrand evaluations. |
| `func MonteCarlo2DWithError(f Func2, ax, bx, ay, by float64, n int, seed int64) QuadResult` | MonteCarlo2DWithError estimates the double integral of f over the rectangle [ax, bx] x [ay, by] together with a standard-error estimate. |
| `func MonteCarloWithError(f Func, a, b float64, n int, seed int64) QuadResult` | MonteCarloWithError estimates the integral of f over [a, b] together with a standard-error estimate and the evaluation count. |

</details>

### Constants

`DefaultMaxDepth`, `DefaultMaxEvals`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/numint).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
