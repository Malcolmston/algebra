# optimize — implements one-dimensional root finding and function optimization for real-valued

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/optimize.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/optimize)

Package optimize implements one-dimensional root finding and function
optimization for real-valued functions of a single real variable.

- Bracketing solvers, which start from an interval [a, b] on which the
function changes sign and are guaranteed to converge: Bisection, RegulaFalsi
(with the Illinois, Pegasus and Anderson-Björck accelerations), Brent's
method (zeroin), Ridders' method and Dekker's method. - Open solvers, which
iterate from one or more starting guesses and converge quickly when close to a
simple root: Secant, Newton, Halley, Schröder (for multiple roots),
Steffensen, inverse quadratic interpolation, Muller's method and fixed-point
iteration (optionally accelerated with Aitken's Δ² process). - Polynomial
solvers, which return all real (or complex) roots of a polynomial: closed
forms for linear, quadratic and cubic equations, and the iterative
Durand-Kerner and Aberth-Ehrlich methods together with a companion-matrix
constructor. - Minimizers, which locate a local minimum of a unimodal
function: golden-section search, Brent's parabolic minimizer, ternary search
and Newton/gradient descent on the derivative, plus routines that bracket a
minimum. - Supporting numerical utilities: finite-difference derivatives, sign
and bracketing predicates, and convergence helpers.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/optimize"
```

## Usage

This is the package's own `ExampleBFGS`, so it compiles and its output is
asserted on every `go test ./optimize/`.

```go
res := optimize.BFGS(optimize.Rosenbrock, optimize.RosenbrockGrad, []float64{-1.2, 1}, 1e-10, 1000)
	fmt.Printf("x=[%.4f %.4f]\n", res.X[0], res.X[1])
```

```
x=[1.0000 1.0000]
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func Aberth(coeffs []float64, tol float64, maxIter int) []complex128` | Aberth returns all n roots (real and complex) of the degree-n polynomial given by coeffs (ascending order) using the Aberth-Ehrlich simultaneous… |
| `func AbsoluteError(a, b float64) float64` | AbsoluteError returns the absolute difference \|a - b\|. |
| `func AitkenDelta2(x0, x1, x2 float64) float64` | AitkenDelta2 applies a single step of Aitken's Δ² process to three consecutive terms x0, x1, x2 of a linearly convergent sequence and returns the… |
| `func AndersonBjorck(f Func, a, b, tol float64, maxIter int) (float64, error)` | AndersonBjorck finds a root of f in [a, b] using the Anderson-Björck variant of the method of false position, whose damping factor 1-fc/fb yields… |
| `func ArmijoStep(f MultiFunc, x, dir, grad []float64) float64` | ArmijoStep is a convenience wrapper around BacktrackingLineSearch using the conventional parameters alpha0 = 1, c = 1e-4, rho = 0.5 and 50 reductions. |
| `func BacktrackingLineSearch(f MultiFunc, x, dir, grad []float64, alpha0, c, rho float64, maxIter int) float64` | BacktrackingLineSearch returns a step length alpha along the descent direction dir from x that satisfies the Armijo (sufficient decrease) condition… |
| `func BackwardDifference(f Func, x, h float64) float64` | BackwardDifference approximates f'(x) by the backward difference (f(x)-f(x-h))/h, a first-order accurate one-sided estimate. |
| `func Bisection(f Func, a, b, tol float64, maxIter int) (float64, error)` | Bisection finds a root of f in the bracketing interval [a, b] by repeated interval halving. |
| `func Booth(x []float64) float64` | Booth is the two-dimensional convex test objective (x0 + 2*x1 - 7)^2 + (2*x0 + x1 - 5)^2, with its unique minimum value 0 at (1, 3). |
| `func BracketExpand(f Func, a, b, factor float64, maxIter int) (float64, float64, bool)` | BracketExpand geometrically expands the interval [a, b] outward (moving the endpoint with the larger \|f\| by a factor of the current width) until f… |
| `func BracketMinimum(f Func, a, b float64) (float64, float64, float64)` | BracketMinimum searches downhill from the initial points a and b to return a triple (a, b, c) with a < b < c (or c < b < a) such that f(b) is less… |
| `func Brent(f Func, a, b, tol float64, maxIter int) (float64, error)` | Brent finds a root of f in [a, b] using Brent's method (the classic zeroin algorithm), which combines the reliability of bisection with the speed of… |
| `func BrentMinimize(f Func, a, b, tol float64, maxIter int) (float64, float64, error)` | BrentMinimize returns the location and value of a minimum of f on [a, b] using Brent's method, which blends golden-section search with parabolic… |
| `func Centroid(points [][]float64) []float64` | Centroid returns the arithmetic mean of a set of points, all assumed to have the same dimension. |
| `func Clamp(x, lo, hi float64) float64` | Clamp constrains x to the closed interval [lo, hi], swapping the bounds if they are given in the wrong order. |
| `func CompanionMatrix(coeffs []float64) [][]float64` | CompanionMatrix returns the companion matrix of the polynomial given by coeffs (ascending order), a real n×n matrix whose characteristic polynomial… |
| `func Converged(prev, cur, tol float64) bool` | Converged reports whether successive iterates prev and cur satisfy a combined relative/absolute convergence test with tolerance tol. |
| `func CubicRealRoots(a, b, c, d float64) []float64` | CubicRealRoots returns the real roots of a*x^3 + b*x^2 + c*x + d = 0 in ascending order (with multiplicity), discarding any genuinely complex roots. |
| `func CubicRoots(a, b, c, d float64) []complex128` | CubicRoots returns all three roots of the cubic a*x^3 + b*x^2 + c*x + d = 0 as complex numbers, using Cardano's formula together with the… |
| `func Dekker(f Func, a, b, tol float64, maxIter int) (float64, error)` | Dekker finds a root of f in [a, b] using Dekker's method, the secant/bisection hybrid that is the historical predecessor of Brent's method. |
| `func Derivative(f Func, x, h float64) float64` | Derivative approximates f'(x) by the central difference (f(x+h)-f(x-h))/(2h), which has second-order accuracy in h. |
| `func DirectionalDerivative(f MultiFunc, x, dir []float64, h float64) float64` | DirectionalDerivative approximates the derivative of f at x along the unit- scaled direction dir using a central finite difference with step h. |
| `func DurandKerner(coeffs []float64, tol float64, maxIter int) []complex128` | DurandKerner returns all n roots (real and complex) of the degree-n polynomial given by coeffs (ascending order) using the Durand-Kerner… |
| `func FalsePosition(f Func, a, b, tol float64, maxIter int) (float64, error)` | FalsePosition is an alias for RegulaFalsi, the traditional English name for the method of false position. |
| `func FixedPoint(g Func, x0, tol float64, maxIter int) (float64, error)` | FixedPoint finds a fixed point of g (a value x with g(x) = x) by the direct iteration x_{n+1} = g(x_n) starting from x0. |
| `func FixedPointAitken(g Func, x0, tol float64, maxIter int) (float64, error)` | FixedPointAitken finds a fixed point of g using Aitken's Δ² acceleration (Steffensen's fixed-point method): it applies g twice and extrapolates,… |
| `func ForwardDifference(f Func, x, h float64) float64` | ForwardDifference approximates f'(x) by the forward difference (f(x+h)-f(x))/h, a first-order accurate one-sided estimate. |
| `func GoldenSectionMin(f Func, a, b, tol float64, maxIter int) (float64, error)` | GoldenSectionMin returns the location of a minimum of the unimodal function f on [a, b] using golden-section search, which narrows the bracket by the… |
| `func GradientDescent(df Func, x0, rate, tol float64, maxIter int) (float64, error)` | GradientDescent locates a stationary point of a function by taking fixed-rate steps against its derivative df from the starting point x0. |
| `func Halley(f, df, d2f Func, x0, tol float64, maxIter int) (float64, error)` | Halley finds a root of f using Halley's method, a third-order iteration that uses the first derivative df and second derivative d2f. |
| `func Himmelblau(x []float64) float64` | Himmelblau is the two-dimensional test objective (x0^2 + x1 - 11)^2 + (x0 + x1^2 - 7)^2, which has four equal local minima of value 0, one of which… |
| `func Illinois(f Func, a, b, tol float64, maxIter int) (float64, error)` | Illinois finds a root of f in [a, b] using the Illinois variant of the method of false position, which halves the retained (stagnant) endpoint's… |
| `func InverseQuadratic(f Func, x0, x1, x2, tol float64, maxIter int) (float64, error)` | InverseQuadratic finds a root of f by inverse quadratic interpolation through three starting points x0, x1 and x2: it fits x as a quadratic in f and… |
| `func IsBracket(f Func, a, b float64) bool` | IsBracket reports whether [a, b] strictly brackets a sign change of f. |
| `func LinearRoot(a, b float64) (float64, bool)` | LinearRoot returns the root of the linear equation a*x + b = 0. |
| `func MatCopy(a [][]float64) [][]float64` | MatCopy returns a deep copy of the matrix a. |
| `func MatIdentity(n int) [][]float64` | MatIdentity returns the n-by-n identity matrix. |
| `func MatTranspose(a [][]float64) [][]float64` | MatTranspose returns the transpose of the (rectangular) matrix a. |
| `func MatVec(a [][]float64, x []float64) []float64` | MatVec returns the matrix-vector product a*x. |
| `func Muller(f func(complex128) complex128, x0, x1, x2 complex128, tol float64, maxIter int) (complex128, error)` | Muller finds a (possibly complex) root of the complex function f using Muller's method, which fits a parabola through three starting points and takes… |
| `func Newton(f, df Func, x0, tol float64, maxIter int) (float64, error)` | Newton finds a root of f using the Newton-Raphson iteration, which requires the derivative df and a single starting guess x0. |
| `func NewtonComplex(f, df func(complex128) complex128, x0 complex128, tol float64, maxIter int) (complex128, error)` | NewtonComplex finds a complex root of the complex function f with derivative df using the Newton-Raphson iteration in the complex plane. |
| `func NewtonMinimize(df, d2f Func, x0, tol float64, maxIter int) (float64, error)` | NewtonMinimize locates a stationary point of a function by applying Newton's method to its derivative, using the first derivative df and second… |
| `func NewtonSafe(f, df Func, a, b, tol float64, maxIter int) (float64, error)` | NewtonSafe finds a root of f in the bracketing interval [a, b] using a safeguarded Newton iteration (the rtsafe algorithm): it takes a Newton step… |
| `func NumericGradient(f MultiFunc, x []float64) []float64` | NumericGradient approximates the gradient of f at x using central finite differences. |
| `func NumericGradientCentral(f MultiFunc, x []float64, h float64) []float64` | NumericGradientCentral approximates the gradient of f at x using central finite differences with step h. |
| `func NumericGradientForward(f MultiFunc, x []float64, h float64) []float64` | NumericGradientForward approximates the gradient of f at x using forward finite differences with step h. |
| `func NumericHessian(f MultiFunc, x []float64, h float64) [][]float64` | NumericHessian approximates the Hessian of f at x using central second differences with step h. |
| `func NumericJacobian(g VectorFunc, x []float64, h float64) [][]float64` | NumericJacobian approximates the Jacobian of the vector-valued function g at x using central finite differences with step h. |
| `func OppositeSign(a, b float64) bool` | OppositeSign reports whether a and b have strictly opposite signs. |
| `func OuterProduct(a, b []float64) [][]float64` | OuterProduct returns the outer product a*b^T, an len(a)-by-len(b) matrix whose (i, j) entry is a`i`*b`j`. |
| `func ParabolicMinimum(x0, x1, x2, f0, f1, f2 float64) (float64, bool)` | ParabolicMinimum fits a parabola through the three points (x0,f0), (x1,f1) and (x2,f2) and returns the abscissa of its vertex. |
| `func PartialDerivative(f MultiFunc, x []float64, i int, h float64) float64` | PartialDerivative approximates the partial derivative of f with respect to coordinate i at x by a central finite difference with step h. |
| `func Pegasus(f Func, a, b, tol float64, maxIter int) (float64, error)` | Pegasus finds a root of f in [a, b] using the Pegasus variant of the method of false position. |
| `func PolyComplexRoots(coeffs []float64) []complex128` | PolyComplexRoots returns all roots of the polynomial given by coeffs (ascending order) as complex numbers, using Durand-Kerner with default settings. |
| `func PolyDeflate(coeffs []float64, root float64) ([]float64, float64)` | PolyDeflate divides the polynomial given by coeffs (ascending order) by the linear factor (x - root) using synthetic division. |
| `func PolyDerivative(coeffs []float64) []float64` | PolyDerivative returns the coefficients (ascending order) of the derivative of the polynomial given by coeffs. |
| `func PolyEval(coeffs []float64, x float64) float64` | PolyEval evaluates the polynomial with the given ascending-order coefficients at x using Horner's scheme. |
| `func PolyEvalComplex(coeffs []float64, x complex128) complex128` | PolyEvalComplex evaluates the real-coefficient polynomial (ascending order) at a complex argument x using Horner's scheme. |
| `func PolyEvalDeriv(coeffs []float64, x float64) (float64, float64)` | PolyEvalDeriv simultaneously evaluates the polynomial (ascending-order coefficients) and its first derivative at x in a single Horner pass, returning… |
| `func PolyIntegral(coeffs []float64, constant float64) []float64` | PolyIntegral returns the coefficients (ascending order) of an antiderivative of the polynomial given by coeffs, using constant as the value of the… |
| `func PolyRoots(coeffs []float64, tol float64) []float64` | PolyRoots returns the real roots of the polynomial given by coeffs (ascending order) in ascending order. |
| `func ProjectBox(x, lo, hi []float64) []float64` | ProjectBox returns a copy of x projected onto the axis-aligned box defined by the per-coordinate bounds lo and hi. |
| `func QuadraticRealRoots(a, b, c float64) []float64` | QuadraticRealRoots returns the distinct real roots of a*x^2 + b*x + c = 0 in ascending order. |
| `func QuadraticRoots(a, b, c float64) (complex128, complex128)` | QuadraticRoots returns both roots of the quadratic a*x^2 + b*x + c = 0 as complex numbers, using a numerically stable formulation that avoids… |
| `func QuarticRoots(a, b, c, d, e float64) []complex128` | QuarticRoots returns all four roots of the quartic a*x^4 + b*x^3 + c*x^2 + d*x + e = 0 as complex numbers, computed with the Durand-Kerner iteration. |
| `func RegulaFalsi(f Func, a, b, tol float64, maxIter int) (float64, error)` | RegulaFalsi (the method of false position) finds a root of f in [a, b] by linear interpolation between the bracketing endpoints. |
| `func RelativeError(approx, exact float64) float64` | RelativeError returns the relative error of an approximation with respect to an exact value. |
| `func RichardsonDerivative(f Func, x, h float64) float64` | RichardsonDerivative approximates f'(x) by Richardson extrapolation of the central difference at step sizes h and h/2, cancelling the leading error… |
| `func Ridders(f Func, a, b, tol float64, maxIter int) (float64, error)` | Ridders finds a root of f in [a, b] using Ridders' method, which fits an exponential to the two endpoints and their midpoint to obtain a… |
| `func Rosenbrock(x []float64) float64` | Rosenbrock is the classic non-convex test objective sum_i [100*(x_{i+1} - x_i^2)^2 + (1 - x_i)^2], with its global minimum value 0 at the all-ones… |
| `func RosenbrockGrad(x []float64) []float64` | RosenbrockGrad returns the exact gradient of Rosenbrock at x. |
| `func SameSign(a, b float64) bool` | SameSign reports whether a and b are both strictly positive or both strictly negative. |
| `func Schroeder(f, df Func, x0 float64, m int, tol float64, maxIter int) (float64, error)` | Schroeder finds a root of known multiplicity m using the modified Newton (or Schröder) iteration x -= m*f/df, which restores quadratic convergence… |
| `func Secant(f Func, x0, x1, tol float64, maxIter int) (float64, error)` | Secant finds a root of f using the secant method, starting from two initial guesses x0 and x1. |
| `func SecondDerivative(f Func, x, h float64) float64` | SecondDerivative approximates f”(x) by the central difference (f(x+h)-2f(x)+f(x-h))/h^2, which has second-order accuracy in h. |
| `func Sign(x float64) float64` | Sign returns +1, -1 or 0 according to the sign of x. |
| `func SignChange(f Func, a, b float64) bool` | SignChange reports whether f takes opposite-signed values at a and b, i.e. |
| `func SolveLinearSystem(a [][]float64, b []float64) ([]float64, error)` | SolveLinearSystem solves the linear system a*x = b for x using Gaussian elimination with partial pivoting. |
| `func Sphere(x []float64) float64` | Sphere is the separable convex test objective sum_i x_i^2, with its unique minimum value 0 at the origin. |
| `func SphereGrad(x []float64) []float64` | SphereGrad returns the exact gradient 2x of Sphere at x. |
| `func Steffensen(f Func, x0, tol float64, maxIter int) (float64, error)` | Steffensen finds a root of f using Steffensen's method, a derivative-free iteration that achieves quadratic convergence by estimating the local slope… |
| `func TernarySearch(f Func, a, b, tol float64, maxIter int) (float64, error)` | TernarySearch returns the location of a minimum of the unimodal function f on [a, b] using ternary search, which discards one outer third of the… |
| `func VecAdd(a, b []float64) []float64` | VecAdd returns the element-wise sum a + b. |
| `func VecAxpy(s float64, x, y []float64) []float64` | VecAxpy returns the combination s*x + y (the classic "a x plus y"). |
| `func VecClamp(x []float64, lo, hi float64) []float64` | VecClamp returns a copy of x with every element confined to the scalar range [lo, hi]. |
| `func VecCopy(x []float64) []float64` | VecCopy returns a fresh copy of x. |
| `func VecDistance(a, b []float64) float64` | VecDistance returns the Euclidean distance between a and b. |
| `func VecDot(a, b []float64) float64` | VecDot returns the Euclidean inner product a . |
| `func VecFill(n int, v float64) []float64` | VecFill returns a new vector of length n with every element set to v. |
| `func VecInfNorm(a []float64) float64` | VecInfNorm returns the maximum-absolute-value (L-infinity) norm of a. |
| `func VecLinComb(a float64, x []float64, b float64, y []float64) []float64` | VecLinComb returns the linear combination a*x + b*y. |
| `func VecNegate(a []float64) []float64` | VecNegate returns -a. |
| `func VecNorm(a []float64) float64` | VecNorm returns the Euclidean (L2) norm \|a\|. |
| `func VecNormSquared(a []float64) float64` | VecNormSquared returns the squared Euclidean norm \|a\|^2. |
| `func VecScale(a []float64, s float64) []float64` | VecScale returns the vector s*a. |
| `func VecSub(a, b []float64) []float64` | VecSub returns the element-wise difference a - b. |
| `func VecZeros(n int) []float64` | VecZeros returns a new zero vector of length n. |
| `func WithinTolerance(a, b, tol float64) bool` | WithinTolerance reports whether a and b differ by at most tol in absolute value. |

### Types

| Type | What it is |
| --- | --- |
| `Bracket` | Bracket describes a closed interval [Lo, Hi] believed to enclose a root or a minimum. |
| `Func` | Func is a real-valued function of a single real variable, the fundamental object operated on by the solvers and minimizers in this package. |
| `GradFunc` | GradFunc returns the gradient vector (vector of first partial derivatives) of an objective evaluated at x. |
| `HessFunc` | HessFunc returns the Hessian matrix (matrix of second partial derivatives) of an objective evaluated at x, as an n-by-n row-major slice of slices. |
| `MultiFunc` | MultiFunc is a scalar-valued objective function of a real vector argument. |
| `Options` | Options bundles the common tuning parameters shared by the iterative multivariate minimizers. |
| `Result` | Result reports the outcome of a multivariate minimization: the located point X, the objective value F at that point, the number of Iterations… |
| `SAOptions` | SAOptions configures the SimulatedAnnealing minimizer. |
| `ScalarFunc` | ScalarFunc is a real-valued function of a single real variable, operated on by the one-dimensional minimizers GoldenSection and BrentParabolic. |
| `ScalarResult` | ScalarResult reports the outcome of a one-dimensional minimization: the located abscissa X, the objective value F there, the number of Iterations… |
| `VectorFunc` | VectorFunc is a vector-valued function of a vector argument, used by the Jacobian utility. |

<details>
<summary><code>Bracket</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func BracketSubdivide(f Func, a, b float64, n int) []Bracket` | BracketSubdivide divides [a, b] into n equal sub-intervals and returns every sub-interval across which f changes sign, each as a Bracket. |
| `func FindBracket(f Func, a, b float64) (Bracket, bool)` | FindBracket expands the interval [a, b] outward until it encloses a root and returns the resulting Bracket together with a flag reporting success. |
| `func (br Bracket) Contains(x float64) bool` | Contains reports whether x lies within the closed bracket, irrespective of whether Lo or Hi is the larger endpoint. |
| `func (br Bracket) Midpoint() float64` | Midpoint returns the arithmetic midpoint of the bracket. |
| `func (br Bracket) Width() float64` | Width returns the (non-negative) length of the bracket. |

</details>

<details>
<summary><code>Options</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func DefaultOptions() Options` | DefaultOptions returns an Options value populated with sensible defaults. |

</details>

<details>
<summary><code>Result</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func BFGS(f MultiFunc, grad GradFunc, x0 []float64, tol float64, maxIter int) Result` | BFGS minimizes f by the BFGS quasi-Newton method starting from x0, maintaining a dense approximation to the inverse Hessian that is refined by the… |
| `func ConjugateGradient(f MultiFunc, grad GradFunc, x0 []float64, tol float64, maxIter int) Result` | ConjugateGradient minimizes f by the nonlinear conjugate-gradient method with Fletcher-Reeves updates and Armijo backtracking line searches, starting… |
| `func CoordinateDescent(f MultiFunc, x0 []float64, window, tol float64, maxIter int) Result` | CoordinateDescent minimizes f by cyclically minimizing along each coordinate axis in turn. |
| `func GradientDescentMomentum(f MultiFunc, grad GradFunc, x0 []float64, rate, momentum, tol float64, maxIter int) Result` | GradientDescentMomentum minimizes f by gradient descent with classical (heavy-ball) momentum starting from x0. |
| `func GradientDescentNesterov(f MultiFunc, grad GradFunc, x0 []float64, rate, momentum, tol float64, maxIter int) Result` | GradientDescentNesterov minimizes f by gradient descent with Nesterov's accelerated (look-ahead) momentum starting from x0. |
| `func NelderMead(f MultiFunc, x0 []float64, step, tol float64, maxIter int) Result` | NelderMead minimizes f by the Nelder-Mead downhill simplex method starting from x0. |
| `func NewtonMultivariate(f MultiFunc, grad GradFunc, hess HessFunc, x0 []float64, tol float64, maxIter int) Result` | NewtonMultivariate minimizes f by the damped Newton method starting from x0. |
| `func SimulatedAnnealing(f MultiFunc, x0 []float64, opts SAOptions, seed int64) Result` | SimulatedAnnealing minimizes f by simulated annealing starting from x0. |

</details>

<details>
<summary><code>SAOptions</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func DefaultSAOptions() SAOptions` | DefaultSAOptions returns an SAOptions value populated with sensible defaults. |

</details>

<details>
<summary><code>ScalarResult</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func BrentParabolic(f ScalarFunc, a, b, tol float64, maxIter int) ScalarResult` | BrentParabolic locates a minimum of the function f on the interval [a, b] using Brent's method, which combines the reliability of golden-section… |
| `func GoldenSection(f ScalarFunc, a, b, tol float64, maxIter int) ScalarResult` | GoldenSection locates a minimum of the unimodal function f on the interval [a, b] by golden-section search, which shrinks the bracket by the constant… |
| `func MinimizeScalar(f ScalarFunc, a, b float64) ScalarResult` | MinimizeScalar locates a minimum of f on [a, b] using BrentParabolic with the DefaultTol and DefaultMaxIter. |

</details>

### Constants

`DefaultLearningRate`, `DefaultMaxIter`, `DefaultMaxIterations`, `DefaultStep`, `DefaultTol`, `DefaultTolerance`, `GoldenRatio`, `InvGoldenRatio`

### Variables

`ErrDimensionMismatch`, `ErrEmptyPoint`, `ErrInvalidInterval`, `ErrMaxIterations`, `ErrNoBracket`, `ErrNonFinite`, `ErrSingularMatrix`, `ErrZeroDerivative`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/optimize).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
