# dynamical — provides primitives for the study of discrete and continuous dynamical systems

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/dynamical.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/dynamical)

Package dynamical provides primitives for the study of discrete and continuous
dynamical systems, implemented with the Go standard library only.

- Iterated one-dimensional maps (`Map1D`) such as the logistic, tent, sine,
Gauss, cubic and circle maps, together with helpers to build parameterized
families as closures. - Iterated two-dimensional maps (`Map2D`) over `Point2D`
values, most notably the Henon map. - Orbits and trajectories: forward
iteration with optional transient discarding, and single n-th iterate
evaluation. - Fixed points and their linear stability, including closed-form
fixed points and multipliers for the logistic and tent maps. - Lyapunov
exponents estimated both from the analytic derivative and by the
trajectory-separation (renormalization) method, in one and two dimensions. -
Bifurcation-diagram sampling for parameterized one-dimensional families. -
Period detection for periodic and eventually-periodic orbits. - Cobweb
(staircase) diagram segment data for visualizing 1-D iteration. - Newton's
method in the real and complex plane, and the basin-of- attraction grid
produced by Newton iteration of a complex map. - Continuous flows in three
dimensions (`Flow3D`) integrated with a fixed-step classical Runge-Kutta (RK4)
scheme, with the Lorenz and Rossler systems provided as concrete examples.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/dynamical"
```

## Usage

This is the package's own `ExampleCircleRotation`, so it compiles and its output is
asserted on every `go test ./dynamical/`.

```go
rot := dynamical.CircleRotation(0.3)
	for _, x := range dynamical.Orbit(rot, 0.0, 4) {
		fmt.Printf("%.1f\n", x)
	}
```

```
0.0
0.3
0.6
0.9
0.2
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func CircleMap(omega, k, theta float64) float64` | CircleMap evaluates the standard circle map f(theta) = theta + omega - (k/2pi)*sin(2pi*theta), reduced modulo 1. |
| `func ClassifyRoot(roots []complex128, z complex128, tol float64) int` | ClassifyRoot returns the index in roots of the entry nearest to z, provided that nearest entry lies within tol; otherwise it returns -1. |
| `func CubicMap(a, x float64) float64` | CubicMap evaluates the cubic map f(x) = a*x - x^3 at x with parameter a. |
| `func DetectPeriod(f Map1D, x0 float64, transient, maxPeriod int, tol float64) int` | DetectPeriod attempts to determine the period of the orbit of f starting at x0. |
| `func DetectPeriod2D(f Map2D, p0 Point2D, transient, maxPeriod int, tol float64) int` | DetectPeriod2D is the two-dimensional analogue of `DetectPeriod` for a map of the plane, measuring return in Euclidean distance. |
| `func Doubling(x float64) float64` | Doubling evaluates the dyadic doubling (Bernoulli) map f(x) = 2x mod 1. |
| `func FixedPoint1D(f, df Map1D, x0 float64, maxIter int, tol float64) (root float64, iters int, converged bool)` | FixedPoint1D finds a fixed point of the map f near x0 by applying Newton's method to g(x) = f(x) - x, using the analytic derivative df of f. |
| `func GaussMap(alpha, beta, x float64) float64` | GaussMap evaluates the Gauss (mouse) map f(x) = exp(-alpha*x^2) + beta. |
| `func Henon(a, b, x, y float64) (float64, float64)` | Henon evaluates one step of the Henon map (x, y) -> (1 - a*x^2 + y, b*x) with parameters a and b. |
| `func HenonJacobian(a, b float64, p Point2D) [2][2]float64` | HenonJacobian returns the 2x2 Jacobian matrix of the Henon map at point p, stored in row-major order: {{-2*a*x, 1}, {b, 0}}. |
| `func IsFixedPoint(f Map1D, x, tol float64) bool` | IsFixedPoint reports whether x is a fixed point of f, that is whether \|f(x) - x\| <= tol. |
| `func IsPeriodicPoint(f Map1D, x0 float64, period int, tol float64) bool` | IsPeriodicPoint reports whether x0 is a periodic point of f with period exactly period, that is whether f^period(x0) returns to within tol of x0… |
| `func Logistic(r, x float64) float64` | Logistic evaluates the logistic map f(x) = r*x*(1-x) at x with parameter r. |
| `func LogisticDeriv(r, x float64) float64` | LogisticDeriv evaluates the derivative f'(x) = r*(1-2x) of the logistic map. |
| `func LogisticFixedPoints(r float64) []float64` | LogisticFixedPoints returns the fixed points of the logistic map with parameter r. |
| `func Lyapunov(f, df Map1D, x0 float64, transient, n int) float64` | Lyapunov estimates the Lyapunov exponent of the one-dimensional map f from the analytic derivative df. |
| `func Lyapunov2D(f Map2D, jac func(Point2D) [2][2]float64, p0 Point2D, transient, n int) float64` | Lyapunov2D estimates the largest Lyapunov exponent of a two-dimensional map f whose Jacobian at any point is supplied by jac (row-major 2x2). |
| `func LyapunovHenon(a, b float64, transient, n int) float64` | LyapunovHenon estimates the largest Lyapunov exponent of the Henon map with parameters a and b by `Lyapunov2D`, starting from the origin. |
| `func LyapunovSeparation(f Map1D, x0, d0 float64, transient, n int) float64` | LyapunovSeparation estimates the Lyapunov exponent of the map f without an analytic derivative, using the trajectory-separation (renormalization)… |
| `func Multiplier(df Map1D, x float64) float64` | Multiplier evaluates the derivative df at x, i.e. |
| `func NewtonBasin(f, df func(complex128) complex128, roots []complex128, xmin, xmax, ymin, ymax float64, nx, ny, maxIter int, tol float64) [][]int` | NewtonBasin computes the Newton basin-of-attraction grid for the complex map f with derivative df over the rectangle [xmin,xmax] x [ymin,ymax],… |
| `func NewtonReal(f, df Map1D, x0 float64, maxIter int, tol float64) (root float64, iters int, converged bool)` | NewtonReal applies Newton's method x -> x - f(x)/f'(x) on the real line, starting at x0, using the function f and its derivative df. |
| `func NthIterate(f Map1D, x0 float64, n int) float64` | NthIterate returns the n-th iterate f^n(x0) of the map f. |
| `func Orbit(f Map1D, x0 float64, n int) []float64` | Orbit returns the forward orbit x_0, x_1, ..., x_n of the one-dimensional map f starting from x0, that is a slice of length n+1 whose first element… |
| `func OrbitTransient(f Map1D, x0 float64, transient, n int) []float64` | OrbitTransient iterates the map f from x0 discarding the first transient points, then returns the following n+1 points of the orbit. |
| `func QuadraticMap(c, x float64) float64` | QuadraticMap evaluates the quadratic map f(x) = x^2 + c, the canonical form to which every real quadratic map is conjugate. |
| `func SineMap(r, x float64) float64` | SineMap evaluates the sine map f(x) = r*sin(pi*x) at x with parameter r. |
| `func Tent(mu, x float64) float64` | Tent evaluates the tent map at x with slope parameter mu: the value is mu*x for x < 1/2 and mu*(1-x) for x >= 1/2. |
| `func TentFixedPoints(mu float64) []float64` | TentFixedPoints returns the fixed points of the tent map with slope mu. |

### Types

| Type | What it is |
| --- | --- |
| `BifurcationPoint` | BifurcationPoint holds the sampled long-term states of a one-dimensional map at a single parameter value. |
| `CobwebSegment` | CobwebSegment is one line segment of a cobweb (staircase) diagram, joining the point (X0, Y0) to (X1, Y1). |
| `Flow3D` | Flow3D is the vector field of a continuous three-dimensional dynamical system: it maps a state s to its time derivative ds/dt = f(s). |
| `LorenzParams` | LorenzParams holds the three parameters of the Lorenz system. |
| `Map1D` | Map1D is a real function of one real variable, used as the iteration rule of a one-dimensional discrete dynamical system x_{n+1} = f(x_n). |
| `Map2D` | Map2D is a map of the plane to itself, used as the iteration rule of a two-dimensional discrete dynamical system p_{n+1} = f(p_n). |
| `NewtonResult` | NewtonResult holds the outcome of a complex Newton iteration: the point Root reached, the number of Iterations performed, and whether the iteration… |
| `Point2D` | Point2D is a point in the plane, used as the state of a two-dimensional discrete dynamical system. |
| `RosslerParams` | RosslerParams holds the three parameters of the Rossler system. |
| `Stability` | Stability classifies the linear stability of a fixed or periodic point of a one-dimensional map according to the magnitude of its multiplier (the map… |
| `Vec3` | Vec3 is a point or vector in three-dimensional space, used as the state of a continuous dynamical system. |

<details>
<summary><code>BifurcationPoint</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Bifurcation(family func(param float64) Map1D, pmin, pmax float64, steps int, x0 float64, transient, samples int) []BifurcationPoint` | Bifurcation samples the asymptotic behavior of a one-parameter family of one-dimensional maps. |
| `func LogisticBifurcation(rmin, rmax float64, steps int, x0 float64, transient, samples int) []BifurcationPoint` | LogisticBifurcation is a convenience wrapper around `Bifurcation` for the logistic family x -> r*x*(1-x), sweeping r from rmin to rmax. |

</details>

<details>
<summary><code>CobwebSegment</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Cobweb(f Map1D, x0 float64, n int) []CobwebSegment` | Cobweb returns the segments of the cobweb diagram for n iterations of the map f starting at x0. |

</details>

<details>
<summary><code>LorenzParams</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func DefaultLorenz() LorenzParams` | DefaultLorenz returns the classical chaotic Lorenz parameters sigma = 10, rho = 28, beta = 8/3. |
| `func (p LorenzParams) Field() Flow3D` | Field returns the Lorenz vector field |
| `func (p LorenzParams) LorenzFixedPoints() []Vec3` | LorenzFixedPoints returns the equilibria of the Lorenz system with parameters p. |

</details>

<details>
<summary><code>Map1D</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func CircleRotation(omega float64) Map1D` | CircleRotation returns the rigid rotation theta -> theta + omega mod 1 as a `Map1D` closure. |
| `func DoublingMap() Map1D` | DoublingMap returns the doubling map x -> 2x mod 1 as a `Map1D` closure. |
| `func LogisticMap(r float64) Map1D` | LogisticMap returns the logistic map x -> r*x*(1-x) as a `Map1D` closure for the fixed parameter r. |
| `func TentMap(mu float64) Map1D` | TentMap returns the tent map with slope parameter mu as a `Map1D` closure. |

</details>

<details>
<summary><code>Map2D</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func HenonMap(a, b float64) Map2D` | HenonMap returns the Henon map with parameters a and b as a `Map2D` closure. |

</details>

<details>
<summary><code>NewtonResult</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewtonComplex(f, df func(complex128) complex128, z0 complex128, maxIter int, tol float64) NewtonResult` | NewtonComplex applies Newton's method z -> z - f(z)/f'(z) in the complex plane, starting at z0, using the function f and its derivative df. |

</details>

<details>
<summary><code>Point2D</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func HenonFixedPoints(a, b float64) []Point2D` | HenonFixedPoints returns the real fixed points of the Henon map with parameters a and b. |
| `func NthIterate2D(f Map2D, p0 Point2D, n int) Point2D` | NthIterate2D returns the n-th iterate f^n(p0) of the two-dimensional map f. |
| `func Orbit2D(f Map2D, p0 Point2D, n int) []Point2D` | Orbit2D returns the forward orbit p_0, p_1, ..., p_n of the two-dimensional map f starting from p0, as a slice of length n+1. |
| `func OrbitTransient2D(f Map2D, p0 Point2D, transient, n int) []Point2D` | OrbitTransient2D iterates the map f from p0 discarding the first transient points, then returns the following n+1 points of the orbit. |
| `func (p Point2D) Add(q Point2D) Point2D` | Add returns the vector sum p + q. |
| `func (p Point2D) Norm() float64` | Norm returns the Euclidean norm (length) of p treated as a vector. |
| `func (p Point2D) Scale(s float64) Point2D` | Scale returns the point p with both coordinates multiplied by s. |
| `func (p Point2D) Sub(q Point2D) Point2D` | Sub returns the vector difference p - q. |

</details>

<details>
<summary><code>RosslerParams</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func DefaultRossler() RosslerParams` | DefaultRossler returns commonly used chaotic Rossler parameters a = 0.2, b = 0.2, c = 5.7. |
| `func (p RosslerParams) Field() Flow3D` | Field returns the Rossler vector field |

</details>

<details>
<summary><code>Stability</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func ClassifyStability(m, tol float64) Stability` | ClassifyStability classifies a fixed point by its multiplier m using tolerance tol to decide the borderline neutral case \|m\| == 1. |
| `func LogisticStability(r, tol float64) []Stability` | LogisticStability returns the stability class of each fixed point returned by `LogisticFixedPoints` for parameter r, using tolerance tol for the… |
| `func (s Stability) String() string` | String returns a human-readable name for the stability class. |

</details>

<details>
<summary><code>Vec3</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Integrate3D(f Flow3D, s0 Vec3, dt float64, n int) []Vec3` | Integrate3D integrates the vector field f from the initial state s0 for n steps of size dt using `RK4Step3D`, returning the trajectory as a slice of… |
| `func LorenzTrajectory(p LorenzParams, s0 Vec3, dt float64, n int) []Vec3` | LorenzTrajectory integrates the Lorenz system with parameters p from state s0 for n RK4 steps of size dt, returning the trajectory (length n+1). |
| `func RK4Step3D(f Flow3D, s Vec3, dt float64) Vec3` | RK4Step3D advances the state s by one step of size dt under the vector field f using the classical fourth-order Runge-Kutta method, returning the new… |
| `func RosslerTrajectory(p RosslerParams, s0 Vec3, dt float64, n int) []Vec3` | RosslerTrajectory integrates the Rossler system with parameters p from state s0 for n RK4 steps of size dt, returning the trajectory (length n+1). |
| `func (v Vec3) Add(w Vec3) Vec3` | Add returns the vector sum v + w. |
| `func (v Vec3) Dot(w Vec3) float64` | Dot returns the dot product v . |
| `func (v Vec3) Norm() float64` | Norm returns the Euclidean norm (length) of v. |
| `func (v Vec3) Scale(s float64) Vec3` | Scale returns v with every component multiplied by s. |
| `func (v Vec3) Sub(w Vec3) Vec3` | Sub returns the vector difference v - w. |

</details>

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/dynamical).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
