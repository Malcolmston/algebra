# odesolvers — provides numerical integrators for systems of first-order ordinary differential

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/odesolvers.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/odesolvers)

Package odesolvers provides numerical integrators for systems of first-order
ordinary differential equations (ODEs) of the form

implemented using only the Go standard library.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/odesolvers"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func AXPY(a float64, x, y []float64) []float64` | AXPY returns a fresh vector a*x + y (the classic "a x plus y"). |
| `func AXPYInPlace(a float64, x, y []float64) []float64` | AXPYInPlace performs y += a*x, mutating and returning y. |
| `func AdamsBashforthCoefficients(order int) ([]float64, error)` | AdamsBashforthCoefficients returns the coefficients of the explicit Adams-Bashforth method of the given order (1..5). |
| `func AdamsBashforthStep(order int, y []float64, fhist [][]float64, h float64) ([]float64, error)` | AdamsBashforthStep advances one explicit Adams-Bashforth step of the given order. |
| `func AdamsMoultonCoefficients(order int) (bNew float64, b []float64, err error)` | AdamsMoultonCoefficients returns the coefficients of the implicit Adams-Moulton method of the given order (1..5). |
| `func Add(a, b []float64) []float64` | Add returns the elementwise sum a+b. |
| `func BDFCoefficients(order int) (alpha []float64, beta float64, err error)` | BDFCoefficients returns the coefficients of the fixed-step backward- differentiation formula of the given order (1..6). |
| `func BDFStep(f Field, order int, tNew float64, hist [][]float64, h float64) ([]float64, error)` | BDFStep advances one fixed-step BDF step of the given order. |
| `func BackwardEulerStep(f Field, t float64, y []float64, h float64) ([]float64, error)` | BackwardEulerStep advances one backward-Euler step and returns the new state. |
| `func Clone(v []float64) []float64` | Clone returns a copy of v. |
| `func CountEvents(sol *Solution, g EventFunc, direction int) int` | CountEvents returns the number of zero crossings of g along sol in the given direction. |
| `func Determinant(a Matrix) float64` | Determinant returns the determinant of the square matrix a via LU factorization with partial pivoting. |
| `func Distance(a, b []float64) float64` | Distance returns the Euclidean distance between a and b. |
| `func Dot(a, b []float64) float64` | Dot returns the Euclidean inner product of a and b. |
| `func EstimateOrder(errH, errH2, errH4 float64) float64` | EstimateOrder estimates the empirical convergence order of a fixed-step integrator from three errors measured at step sizes h, h/2 and h/4: order ≈… |
| `func EulerStep(f Field, t float64, y []float64, h float64) []float64` | EulerStep advances one forward-Euler step and returns the new state. |
| `func Fill(n int, x float64) []float64` | Fill returns a newly allocated slice of length n with every element set to x. |
| `func GaussLegendre4Step(f Field, t float64, y []float64, h float64) ([]float64, error)` | GaussLegendre4Step advances one two-stage Gauss-Legendre (order 4) step. |
| `func GlobalError(sol *Solution, exact []float64) float64` | GlobalError returns the infinity-norm difference between the final state of a solution and a known exact reference state. |
| `func HarmonicEnergy(q, v []float64, omega float64) float64` | HarmonicEnergy returns the total energy 0.5*\|v\|^2 + 0.5*omega^2*\|q\|^2 of a unit-mass isotropic harmonic oscillator, a convenient conserved… |
| `func HermiteInterpolate(theta, h float64, y0, y1, f0, f1 []float64) []float64` | HermiteInterpolate evaluates the cubic Hermite interpolant on a single step. |
| `func HeunStep(f Field, t float64, y []float64, h float64) []float64` | HeunStep advances one step of Heun's method and returns the new state. |
| `func ImplicitMidpointStep(f Field, t float64, y []float64, h float64) ([]float64, error)` | ImplicitMidpointStep advances one implicit-midpoint step and returns the new state. |
| `func ImplicitRKStep(f Field, bt *ButcherTableau, t float64, y []float64, h float64) ([]float64, error)` | ImplicitRKStep advances one step of a general (possibly implicit) Runge-Kutta method from (t, y) by h and returns the new state. |
| `func IntegrateUntilEvent(f Field, g EventFunc, direction int, t0 float64, y0 []float64, tMax, h float64) (*Solution, Event, bool)` | IntegrateUntilEvent integrates y' = f(t, y) with the fixed-step RK4 method from t0 with step h and stops at the first zero crossing of g in the given… |
| `func LeapfrogStep(a Acceleration, t float64, q, v []float64, h float64) (qNew, vNew []float64)` | LeapfrogStep advances one leapfrog step. |
| `func LinearCombination(n int, coeffs []float64, vecs [][]float64) []float64` | LinearCombination returns sum_i coeffs`i` * vecs`i`. |
| `func Linspace(a, b float64, n int) []float64` | Linspace returns n evenly spaced samples over the closed interval [a, b]. |
| `func MatVec(a Matrix, x []float64) []float64` | MatVec returns the matrix-vector product a*x. |
| `func MaxAbs(v []float64) float64` | MaxAbs returns the largest absolute value among the elements of v, or 0 for an empty slice. |
| `func MaxComponentError(a, b []float64) float64` | MaxComponentError returns the maximum absolute difference between the components of a and b, a convenience for comparing trajectories in tests. |
| `func MidpointStep(f Field, t float64, y []float64, h float64) []float64` | MidpointStep advances one explicit-midpoint step and returns the new state. |
| `func Neg(v []float64) []float64` | Neg returns the elementwise negation -v. |
| `func NewtonSolve(F func(x []float64) []float64, x0 []float64, tol float64, maxIter int) ([]float64, error)` | NewtonSolve solves the square nonlinear system F(x) = 0 by Newton's method with a finite-difference Jacobian, starting from x0. |
| `func Norm1(v []float64) float64` | Norm1 returns the L1 norm (sum of absolute values) of v. |
| `func Norm2(v []float64) float64` | Norm2 returns the Euclidean (L2) norm of v. |
| `func NormInf(v []float64) float64` | NormInf returns the maximum-absolute-value (infinity) norm of v. |
| `func PositionVerletStep(a Acceleration, t float64, q, v []float64, h float64) (qNew, vNew []float64)` | PositionVerletStep advances one position-Verlet (drift-kick-drift) step of the system q” = a(t, q) and returns the updated position and velocity. |
| `func RK38Step(f Field, t float64, y []float64, h float64) []float64` | RK38Step advances one step of the fourth-order 3/8-rule and returns the new state. |
| `func RK4Step(f Field, t float64, y []float64, h float64) []float64` | RK4Step advances one classical fourth-order Runge-Kutta step and returns the new state. |
| `func RMSNorm(v []float64) float64` | RMSNorm returns the root-mean-square norm sqrt(mean(v_i^2)) of v. |
| `func RadauIIA3Step(f Field, t float64, y []float64, h float64) ([]float64, error)` | RadauIIA3Step advances one order-3 Radau IIA step and returns the new state. |
| `func RadauIIA5Step(f Field, t float64, y []float64, h float64) ([]float64, error)` | RadauIIA5Step advances one order-5 Radau IIA step and returns the new state. |
| `func RalstonStep(f Field, t float64, y []float64, h float64) []float64` | RalstonStep advances one step of Ralston's method and returns the new state. |
| `func RefineEvent(sol *Solution, g EventFunc, ta, tb float64) float64` | RefineEvent locates the crossing time of g within the bracket [ta, tb], where g(ta) and g(tb) are assumed to have opposite signs. |
| `func RichardsonExtrapolate(coarse, fine []float64, ratio float64, order int) []float64` | RichardsonExtrapolate combines two approximations coarse and fine, computed with step sizes differing by the factor ratio (>1) for a method of the… |
| `func Scale(s float64, v []float64) []float64` | Scale returns the vector s*v. |
| `func ShootingResidual(f Field, bc BoundaryFunc, t0 float64, s []float64, tEnd float64, subSteps int) []float64` | ShootingResidual returns, for a candidate initial state s at t0, the boundary residual bc(s, y(tEnd; s)) obtained by propagating s to tEnd with RK4… |
| `func SolveLinearSystem(a Matrix, b []float64) ([]float64, error)` | SolveLinearSystem solves A x = b for x by Gaussian elimination with partial pivoting. |
| `func Sub(a, b []float64) []float64` | Sub returns the elementwise difference a-b. |
| `func Trace(a Matrix) float64` | Trace returns the sum of the diagonal entries of the square matrix a. |
| `func TrapezoidalStep(f Field, t float64, y []float64, h float64) ([]float64, error)` | TrapezoidalStep advances one implicit-trapezoidal step and returns the new state. |
| `func VelocityVerletStep(a Acceleration, t float64, q, v, aCur []float64, h float64) (qNew, vNew, aNew []float64)` | VelocityVerletStep advances one velocity-Verlet (kick-drift-kick) step of the system q” = a(t, q). |
| `func WeightedRMSNorm(v, scale []float64) float64` | WeightedRMSNorm returns the scaled RMS norm sqrt(mean((v_i/scale_i)^2)), the error measure used by the adaptive integrators. |
| `func Yoshida4Step(a Acceleration, t float64, q, v []float64, h float64) (qNew, vNew []float64)` | Yoshida4Step advances one fourth-order Yoshida step by composing three velocity-Verlet substeps of sizes w1*h, w0*h and w1*h. |
| `func YoshidaCoefficients() (w1, w0 float64)` | YoshidaCoefficients returns the three composition weights (w1, w0, w1) of the fourth-order Yoshida method built from a symmetric second-order base… |
| `func Zeros(n int) []float64` | Zeros returns a newly allocated slice of n zeros. |

### Types

| Type | What it is |
| --- | --- |
| `Acceleration` | Acceleration is the force field of a separable second-order system q” = a(t, q), returning the acceleration (force per unit mass) at configuration… |
| `AdaptiveOptions` | AdaptiveOptions configures the embedded-Runge-Kutta step-size controller. |
| `BoundaryFunc` | BoundaryFunc encodes the boundary conditions of a two-point boundary-value problem. |
| `ButcherTableau` | ButcherTableau describes a Runge-Kutta method by its coefficient arrays. |
| `Event` | Event records a detected zero crossing of an `EventFunc`. |
| `EventFunc` | EventFunc is a scalar event (switching) function g(t, y). |
| `Field` | Field is the right-hand side of a first-order ODE system y' = f(t, y). |
| `Matrix` | Matrix is a dense row-major matrix used by the implicit solvers. |
| `ShootingOptions` | ShootingOptions configures the shooting solvers. |
| `Solution` | Solution stores the discrete trajectory produced by an integrator. |
| `SymplecticSolution` | SymplecticSolution stores the trajectory of a second-order system produced by a symplectic integrator: positions Q and velocities P sampled at times… |

<details>
<summary><code>AdaptiveOptions</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func DefaultAdaptiveOptions() AdaptiveOptions` | DefaultAdaptiveOptions returns a reasonable default configuration: RelTol=1e-6, AbsTol=1e-9, automatic initial step, up to 100000 steps. |

</details>

<details>
<summary><code>ButcherTableau</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func BackwardEulerTableau() *ButcherTableau` | BackwardEulerTableau returns the tableau of the backward (implicit) Euler method, order 1, L-stable. |
| `func BogackiShampineTableau() *ButcherTableau` | BogackiShampineTableau returns the embedded Bogacki-Shampine 3(2) pair used by MATLAB's ode23. |
| `func CashKarpTableau() *ButcherTableau` | CashKarpTableau returns the embedded Cash-Karp 5(4) pair. |
| `func DormandPrinceTableau() *ButcherTableau` | DormandPrinceTableau returns the embedded Dormand-Prince DOPRI5 5(4) pair used by MATLAB's ode45. |
| `func EulerTableau() *ButcherTableau` | EulerTableau returns the tableau of the forward (explicit) Euler method, order 1. |
| `func Fehlberg78Tableau() *ButcherTableau` | Fehlberg78Tableau returns the high-order Runge-Kutta-Fehlberg RKF78 7(8) embedded pair. |
| `func FehlbergTableau() *ButcherTableau` | FehlbergTableau returns the classical Runge-Kutta-Fehlberg RKF45 embedded pair. |
| `func GaussLegendre4Tableau() *ButcherTableau` | GaussLegendre4Tableau returns the two-stage, fourth-order Gauss-Legendre tableau, A-stable and symplectic. |
| `func HeunEulerTableau() *ButcherTableau` | HeunEulerTableau returns the embedded Heun-Euler 2(1) pair, the simplest adaptive Runge-Kutta method. |
| `func HeunTableau() *ButcherTableau` | HeunTableau returns the tableau of Heun's method (the explicit trapezoidal / improved Euler method), order 2. |
| `func HeunThirdOrderTableau() *ButcherTableau` | HeunThirdOrderTableau returns Heun's third-order explicit Runge-Kutta tableau. |
| `func ImplicitMidpointTableau() *ButcherTableau` | ImplicitMidpointTableau returns the tableau of the implicit midpoint method (the one-stage Gauss method), order 2, A-stable and symplectic. |
| `func KuttaThirdOrderTableau() *ButcherTableau` | KuttaThirdOrderTableau returns Kutta's classical third-order explicit Runge-Kutta tableau. |
| `func MidpointTableau() *ButcherTableau` | MidpointTableau returns the tableau of the explicit midpoint method, order 2. |
| `func RK38Tableau() *ButcherTableau` | RK38Tableau returns the "3/8-rule" fourth-order Runge-Kutta tableau, an alternative to RK4 with smaller error coefficients. |
| `func RK4Tableau() *ButcherTableau` | RK4Tableau returns the classical fourth-order Runge-Kutta tableau. |
| `func RadauIIA3Tableau() *ButcherTableau` | RadauIIA3Tableau returns the two-stage Radau IIA tableau of order 3, L-stable and stiffly accurate. |
| `func RadauIIA5Tableau() *ButcherTableau` | RadauIIA5Tableau returns the three-stage Radau IIA tableau of order 5, L-stable and stiffly accurate, the workhorse implicit method for stiff systems. |
| `func RalstonTableau() *ButcherTableau` | RalstonTableau returns the tableau of Ralston's method, the second-order explicit method with minimal truncation error. |
| `func SSPRK3Tableau() *ButcherTableau` | SSPRK3Tableau returns the strong-stability-preserving third-order explicit Runge-Kutta tableau (Shu-Osher). |
| `func TrapezoidalTableau() *ButcherTableau` | TrapezoidalTableau returns the tableau of the implicit trapezoidal rule (Crank-Nicolson), order 2, A-stable. |
| `func (bt *ButcherTableau) IsConsistent(tol float64) bool` | IsConsistent reports whether the tableau satisfies the first-order consistency conditions to within tol: the weights B sum to 1 and every row of A… |
| `func (bt *ButcherTableau) IsEmbedded() bool` | IsEmbedded reports whether the tableau carries an embedded error estimate. |
| `func (bt *ButcherTableau) IsExplicit() bool` | IsExplicit reports whether the tableau is explicit, i.e. |
| `func (bt *ButcherTableau) RowSums() []float64` | RowSums returns the sum of each row of the tableau's A matrix. |
| `func (bt *ButcherTableau) Stages() int` | Stages returns the number of stages s of the method. |
| `func (bt *ButcherTableau) Step(f Field, t float64, y []float64, h float64) []float64` | Step advances one explicit Runge-Kutta step from (t, y) by h and returns the new state. |
| `func (bt *ButcherTableau) StepEmbedded(f Field, t float64, y []float64, h float64) (ynext, err []float64)` | StepEmbedded advances one explicit step and, for an embedded tableau, also returns the local error estimate err = y_primary - y_star. |
| `func (bt *ButcherTableau) SumB() float64` | SumB returns the sum of the primary weights B, which must equal 1 for a consistent method. |

</details>

<details>
<summary><code>Event</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func FirstEvent(sol *Solution, g EventFunc, direction int) (Event, bool)` | FirstEvent returns the earliest event of g along sol in the given direction, or ok == false when none occur. |
| `func ScanEvents(sol *Solution, g EventFunc, direction int) []Event` | ScanEvents scans a computed Solution for zero crossings of the event function g. |

</details>

<details>
<summary><code>Matrix</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func CopyMatrix(a Matrix) Matrix` | CopyMatrix returns a deep copy of a. |
| `func FiniteDiffJacobian(F func(x []float64) []float64, x []float64, eps float64) Matrix` | FiniteDiffJacobian approximates the Jacobian dF/dx of the vector field F at x using first-order forward differences with step eps (a sensible default… |
| `func IdentityMatrix(n int) Matrix` | IdentityMatrix returns the n-by-n identity matrix. |
| `func MatAdd(a, b Matrix) Matrix` | MatAdd returns the elementwise sum a+b. |
| `func MatMul(a, b Matrix) Matrix` | MatMul returns the matrix product a*b. |
| `func MatScale(s float64, a Matrix) Matrix` | MatScale returns the matrix s*a. |
| `func NewMatrix(rows, cols int) Matrix` | NewMatrix returns a rows-by-cols matrix of zeros. |
| `func Transpose(a Matrix) Matrix` | Transpose returns the transpose of a. |

</details>

<details>
<summary><code>ShootingOptions</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func DefaultShootingOptions() ShootingOptions` | DefaultShootingOptions returns SubSteps=200, Tol=1e-8, MaxIter=50. |

</details>

<details>
<summary><code>Solution</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func IntegrateAdaptive(f Field, bt *ButcherTableau, t0 float64, y0 []float64, tEnd float64, opts AdaptiveOptions) (*Solution, error)` | IntegrateAdaptive integrates y' = f(t, y) from t0 to tEnd with the embedded tableau bt and automatic step-size control. |
| `func IntegrateFixed(f Field, bt *ButcherTableau, t0 float64, y0 []float64, tEnd, h float64) *Solution` | IntegrateFixed integrates y' = f(t, y) from t0 to tEnd using the explicit method described by bt with a fixed step of nominal size h. |
| `func IntegrateFixedDense(f Field, bt *ButcherTableau, t0 float64, y0 []float64, tEnd, h float64) *Solution` | IntegrateFixedDense behaves like `IntegrateFixed` but also records the stage-0 derivative at every sample so that the returned Solution supports… |
| `func IntegrateImplicit(f Field, bt *ButcherTableau, t0 float64, y0 []float64, tEnd, h float64) (*Solution, error)` | IntegrateImplicit integrates y' = f(t, y) from t0 to tEnd with a fixed step h using the implicit Runge-Kutta method bt, solving each stage system… |
| `func MultipleShooting(f Field, bc BoundaryFunc, t0, tEnd float64, nIntervals int, guess []float64, opts ShootingOptions) ([]float64, [][]float64, *Solution, error)` | MultipleShooting solves the boundary-value problem with multiple shooting. |
| `func SingleShooting(f Field, bc BoundaryFunc, t0 float64, tEnd float64, guess []float64, opts ShootingOptions) ([]float64, *Solution, error)` | SingleShooting solves the two-point boundary-value problem y' = f(t, y) on [t0, tEnd] with boundary conditions bc, using single shooting: it applies… |
| `func SolveABM(f Field, order, corrections int, t0 float64, y0 []float64, tEnd, h float64) (*Solution, error)` | SolveABM integrates with the Adams-Bashforth-Moulton predictor-corrector of the given order (1..5) and fixed step h, applying corrections iterations… |
| `func SolveAdamsBashforth(f Field, order int, t0 float64, y0 []float64, tEnd, h float64) (*Solution, error)` | SolveAdamsBashforth integrates with the explicit Adams-Bashforth method of the given order (1..5) and fixed step h, bootstrapping the required… |
| `func SolveBDF(f Field, order int, t0 float64, y0 []float64, tEnd, h float64) (*Solution, error)` | SolveBDF integrates y' = f(t, y) from t0 to tEnd with the fixed-step BDF method of the given order (1..6) and step h. |
| `func SolveBackwardEuler(f Field, t0 float64, y0 []float64, tEnd, h float64) (*Solution, error)` | SolveBackwardEuler integrates with the backward Euler method and a fixed step h. |
| `func SolveBogackiShampine(f Field, t0 float64, y0 []float64, tEnd float64, opts AdaptiveOptions) (*Solution, error)` | SolveBogackiShampine integrates adaptively with the Bogacki-Shampine 3(2) pair (ode23). |
| `func SolveCashKarp(f Field, t0 float64, y0 []float64, tEnd float64, opts AdaptiveOptions) (*Solution, error)` | SolveCashKarp integrates adaptively with the Cash-Karp 5(4) pair. |
| `func SolveDOPRI5(f Field, t0 float64, y0 []float64, tEnd float64, opts AdaptiveOptions) (*Solution, error)` | SolveDOPRI5 integrates adaptively with the Dormand-Prince 5(4) pair (ode45). |
| `func SolveEuler(f Field, t0 float64, y0 []float64, tEnd, h float64) *Solution` | SolveEuler integrates with the forward Euler method and a fixed step h. |
| `func SolveGaussLegendre4(f Field, t0 float64, y0 []float64, tEnd, h float64) (*Solution, error)` | SolveGaussLegendre4 integrates with the fourth-order Gauss-Legendre method and a fixed step h. |
| `func SolveHeun(f Field, t0 float64, y0 []float64, tEnd, h float64) *Solution` | SolveHeun integrates with Heun's method and a fixed step h. |
| `func SolveHeunEuler(f Field, t0 float64, y0 []float64, tEnd float64, opts AdaptiveOptions) (*Solution, error)` | SolveHeunEuler integrates adaptively with the Heun-Euler 2(1) pair. |
| `func SolveImplicitMidpoint(f Field, t0 float64, y0 []float64, tEnd, h float64) (*Solution, error)` | SolveImplicitMidpoint integrates with the implicit midpoint method and a fixed step h. |
| `func SolveKutta3(f Field, t0 float64, y0 []float64, tEnd, h float64) *Solution` | SolveKutta3 integrates with Kutta's third-order method and a fixed step h. |
| `func SolveMidpoint(f Field, t0 float64, y0 []float64, tEnd, h float64) *Solution` | SolveMidpoint integrates with the explicit midpoint method and a fixed step h. |
| `func SolveRK38(f Field, t0 float64, y0 []float64, tEnd, h float64) *Solution` | SolveRK38 integrates with the fourth-order 3/8-rule and a fixed step h. |
| `func SolveRK4(f Field, t0 float64, y0 []float64, tEnd, h float64) *Solution` | SolveRK4 integrates with the classical fourth-order Runge-Kutta method and a fixed step h. |
| `func SolveRKF45(f Field, t0 float64, y0 []float64, tEnd float64, opts AdaptiveOptions) (*Solution, error)` | SolveRKF45 integrates adaptively with the Runge-Kutta-Fehlberg 4(5) pair. |
| `func SolveRKF78(f Field, t0 float64, y0 []float64, tEnd float64, opts AdaptiveOptions) (*Solution, error)` | SolveRKF78 integrates adaptively with the high-order Fehlberg 7(8) pair. |
| `func SolveRadauIIA3(f Field, t0 float64, y0 []float64, tEnd, h float64) (*Solution, error)` | SolveRadauIIA3 integrates with the order-3 Radau IIA method and a fixed step h. |
| `func SolveRadauIIA5(f Field, t0 float64, y0 []float64, tEnd, h float64) (*Solution, error)` | SolveRadauIIA5 integrates with the order-5 Radau IIA method and a fixed step h. |
| `func SolveRalston(f Field, t0 float64, y0 []float64, tEnd, h float64) *Solution` | SolveRalston integrates with Ralston's method and a fixed step h. |
| `func SolveSSPRK3(f Field, t0 float64, y0 []float64, tEnd, h float64) *Solution` | SolveSSPRK3 integrates with the SSPRK3 method and a fixed step h. |
| `func SolveTrapezoidal(f Field, t0 float64, y0 []float64, tEnd, h float64) (*Solution, error)` | SolveTrapezoidal integrates with the implicit trapezoidal rule and a fixed step h. |
| `func (s *Solution) At(t float64) []float64` | At returns the interpolated state at time t. |
| `func (s *Solution) Component(j int) []float64` | Component returns the time series of state component j across all samples. |
| `func (s *Solution) Dim() int` | Dim returns the dimension of the state vectors, or 0 when the solution is empty. |
| `func (s *Solution) Final() []float64` | Final returns a copy of the last stored state, or nil when the solution holds no samples. |
| `func (s *Solution) FinalTime() float64` | FinalTime returns the last stored time, or NaN when the solution holds no samples. |
| `func (s *Solution) Initial() []float64` | Initial returns a copy of the first stored state, or nil when the solution holds no samples. |
| `func (s *Solution) InitialTime() float64` | InitialTime returns the first stored time, or NaN when the solution holds no samples. |
| `func (s *Solution) Interpolate(t float64) []float64` | Interpolate is a synonym for `Solution.At`. |
| `func (s *Solution) Len() int` | Len returns the number of stored samples. |
| `func (s *Solution) Sample(n int) ([]float64, [][]float64)` | Sample evaluates the solution at n evenly spaced times over its covered interval and returns the times together with the interpolated states. |
| `func (s *Solution) Steps() int` | Steps returns the number of integration steps, one fewer than the number of stored samples (or zero when empty). |
| `func (s *Solution) Times() []float64` | Times returns a copy of the sample-time slice. |

</details>

<details>
<summary><code>SymplecticSolution</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func SolveLeapfrog(a Acceleration, t0 float64, q0, v0 []float64, tEnd, h float64) *SymplecticSolution` | SolveLeapfrog integrates q” = a(t, q) with the leapfrog method and a fixed step h. |
| `func SolvePositionVerlet(a Acceleration, t0 float64, q0, v0 []float64, tEnd, h float64) *SymplecticSolution` | SolvePositionVerlet integrates q” = a(t, q) with the position-Verlet method and a fixed step h. |
| `func SolveVelocityVerlet(a Acceleration, t0 float64, q0, v0 []float64, tEnd, h float64) *SymplecticSolution` | SolveVelocityVerlet integrates q” = a(t, q) from t0 to tEnd with a fixed step h using the velocity-Verlet method, starting from position q0 and… |
| `func SolveYoshida4(a Acceleration, t0 float64, q0, v0 []float64, tEnd, h float64) *SymplecticSolution` | SolveYoshida4 integrates q” = a(t, q) with the fourth-order Yoshida method and a fixed step h. |
| `func (s *SymplecticSolution) Dim() int` | Dim returns the configuration-space dimension. |
| `func (s *SymplecticSolution) FinalPosition() []float64` | FinalPosition returns a copy of the last stored position, or nil when the solution holds no samples. |
| `func (s *SymplecticSolution) FinalVelocity() []float64` | FinalVelocity returns a copy of the last stored velocity, or nil when the solution holds no samples. |
| `func (s *SymplecticSolution) Len() int` | Len returns the number of stored samples. |
| `func (s *SymplecticSolution) PhaseState(i int) []float64` | PhaseState returns the combined phase-space vector [q, v] at sample i. |

</details>

### Variables

`ErrDimensionMismatch`, `ErrMaxSteps`, `ErrStepTooSmall`, `ErrSingularMatrix`, `ErrNoConvergence`, `ErrInvalidOrder`, `ErrInvalidInput`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/odesolvers).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
