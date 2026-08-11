# chaos — implements tools for the study of chaos and nonlinear dynamics using only the Go

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/chaos.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/chaos)

Package chaos implements tools for the study of chaos and nonlinear dynamics
using only the Go standard library.

- Canonical systems: the one-dimensional logistic, tent, sine, Gauss and cubic
maps (`Map1D`); the two-dimensional Henon and standard (Chirikov) maps
(`Map2D`); and the continuous Lorenz and Rossler flows (`Flow`) integrated
with fixed-step Runge-Kutta schemes.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/chaos"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func AttractorSet(f Map1D, x0 float64, transient, keep int, tol float64) []float64` | AttractorSet iterates f from x0, discards the transient and returns the distinct attractor values found among the next keep samples, sorted ascending… |
| `func Autocorrelation(xs []float64, lag int) float64` | Autocorrelation returns the sample autocorrelation of xs at the given lag. |
| `func BenettinLargestFlow(f Field, x0 Vec, h float64, transient, n int) float64` | BenettinLargestFlow estimates the largest Lyapunov exponent of the flow with field f by evolving a tangent vector along the trajectory with the… |
| `func BenettinLargestMap(F MapN, x0 Vec, transient, n int) float64` | BenettinLargestMap estimates the largest Lyapunov exponent of the map F using the Benettin single-perturbation method: a tangent vector is evolved by… |
| `func BoxCount(points []Vec, eps float64) int` | BoxCount returns the number of occupied boxes of side eps needed to cover the set of points, where each point is a Vec of common dimension. |
| `func BoxCountingDimension(points []Vec, epsMin, epsMax float64, scales int) (dim float64, logInvEps, logN []float64)` | BoxCountingDimension estimates the box-counting (capacity) dimension of a set of points by counting occupied boxes over a geometric sequence of box… |
| `func Cobweb(f Map1D, x0 float64, n int) (xs, ys []float64)` | Cobweb returns the vertices of the cobweb (staircase) diagram of f started at x0 for n steps, as an alternating sequence of points on the diagonal… |
| `func CorrelationDimension(points []Vec, rMin, rMax float64, scales int) (dim float64, logR, logC []float64)` | CorrelationDimension estimates the correlation dimension D2 by the Grassberger-Procaccia method: it evaluates the correlation sum over a geometric… |
| `func CorrelationSum(points []Vec, r float64) float64` | CorrelationSum returns the Grassberger-Procaccia correlation sum C(r): the fraction of distinct point pairs whose Euclidean distance is below r. |
| `func Det(A Mat) float64` | Det returns the determinant of a square matrix by Gaussian elimination with partial pivoting. |
| `func Det2(A Mat) float64` | Det2 returns the determinant of the 2-by-2 matrix A. |
| `func Det3(A Mat) float64` | Det3 returns the determinant of the 3-by-3 matrix A. |
| `func DetectPeriod(f Map1D, x0 float64, transient, maxPeriod int, tol float64) int` | DetectPeriod iterates f from x0, discards the transient, then determines the period of the resulting orbit up to maxPeriod, using tolerance tol. |
| `func Eigenvalues2(A Mat) [2]complex128` | Eigenvalues2 returns the two (possibly complex) eigenvalues of a 2-by-2 matrix A. |
| `func Eigenvalues3(A Mat) [3]complex128` | Eigenvalues3 returns the three (possibly complex) eigenvalues of a 3-by-3 matrix A, obtained by solving its characteristic cubic. |
| `func FeigenbaumEstimate(r1, r2, r3 float64) float64` | FeigenbaumEstimate estimates the first Feigenbaum constant delta from three consecutive period-doubling bifurcation parameters r1<r2<r3 via… |
| `func FeigenbaumFromSequence(rs []float64) []float64` | FeigenbaumFromSequence estimates the Feigenbaum delta constant from a sequence of successive bifurcation parameters, returning the ratio for each… |
| `func FitLine(xs, ys []float64) (m, b float64)` | FitLine returns the least-squares slope and intercept of the line y = m x + b through the paired samples xs, ys. |
| `func FitSlope(xs, ys []float64) float64` | FitSlope returns just the least-squares slope of ys against xs. |
| `func FixedPoint1D(f Map1D, x0, tol float64, maxIter int) (float64, error)` | FixedPoint1D refines a fixed point of the one-dimensional map f near the guess x0 by Newton's method applied to g(x)=f(x)-x, using a numerical… |
| `func GeneralizedDimension(points []Vec, q, epsMin, epsMax float64, scales int) float64` | GeneralizedDimension estimates the Renyi generalized dimension D_q of order q (q != 1) from box occupancy probabilities. |
| `func Histogram(xs []float64, lo, hi float64, bins int) (counts []int, centers []float64)` | Histogram bins the values xs into bins equal-width buckets over the range [lo, hi] and returns the bucket counts and their centres. |
| `func InformationDimension(points []Vec, epsMin, epsMax float64, scales int) (dim float64, logInvEps, entropy []float64)` | InformationDimension estimates the information (D1) dimension of a set of points using box occupancy probabilities: for each box size it forms the… |
| `func InvariantDensity(f Map1D, x0, lo, hi float64, bins, transient, n int) (centers, density []float64)` | InvariantDensity estimates the invariant density of the one-dimensional map f by iterating from x0, discarding the transient and forming a normalised… |
| `func IsHyperbolicFlow(J Mat) bool` | IsHyperbolicFlow reports whether every eigenvalue of J has non-zero real part (a hyperbolic equilibrium of a flow). |
| `func IsHyperbolicMap(J Mat) bool` | IsHyperbolicMap reports whether every eigenvalue of J has magnitude different from one (a hyperbolic fixed point of a map). |
| `func IsStable1D(f Map1D, x float64) bool` | IsStable1D reports whether the fixed point x of the map f is linearly stable, i.e. |
| `func Iterate(f Map1D, x float64, n int) float64` | Iterate applies f to x n times and returns the final value. |
| `func KaplanYorkeDimension(spectrum Vec) float64` | KaplanYorkeDimension returns the Kaplan-Yorke (Lyapunov) dimension implied by a Lyapunov spectrum. |
| `func Linspace(a, b float64, n int) []float64` | Linspace returns n evenly spaced values from a to b inclusive. |
| `func LocalMaxima(f Field, x0 Vec, h float64, transient, n, coord int) []float64` | LocalMaxima returns the successive local maxima of the coordinate coord along the trajectory of f. |
| `func LogisticAnalyticLyapunov() float64` | LogisticAnalyticLyapunov returns the exact Lyapunov exponent of the fully chaotic logistic map at r=4, which equals log 2. |
| `func LogisticFixedPoints(r float64) []float64` | LogisticFixedPoints returns the two analytic fixed points of the logistic map with parameter r: 0 and 1 - 1/r. |
| `func Logspace(a, b float64, n int) []float64` | Logspace returns n values geometrically spaced between a and b (both > 0). |
| `func Lyapunov1D(f Map1D, x0 float64, transient, n int) float64` | Lyapunov1D estimates the Lyapunov exponent of the map f using a numerical derivative, so no analytic derivative is required. |
| `func LyapunovExponentDiagram(fam Family1D, pMin, pMax float64, steps int, x0 float64, transient, n int) (params, lambdas []float64)` | LyapunovExponentDiagram returns the Lyapunov exponent of the family fam at each of steps parameter values across [pMin, pMax]. |
| `func LyapunovLog1D(f, df Map1D, x0 float64, transient, n int) float64` | LyapunovLog1D estimates the Lyapunov exponent of the one-dimensional map f from its analytic derivative df. |
| `func LyapunovSeparation1D(f Map1D, x0, d0 float64, transient, n int) float64` | LyapunovSeparation1D estimates the Lyapunov exponent of the map f by the trajectory-separation method: two orbits start a distance d0 apart and the… |
| `func Mean(xs []float64) float64` | Mean returns the arithmetic mean of xs, or NaN when empty. |
| `func Median(xs []float64) float64` | Median returns the median of xs (which is copied and sorted internally). |
| `func MetricEntropy(spectrum Vec) float64` | MetricEntropy returns the Kolmogorov-Sinai entropy estimated by Pesin's identity as the sum of the positive Lyapunov exponents. |
| `func Multiplier1D(f Map1D, x float64) float64` | Multiplier1D returns the derivative (multiplier) of the map f at the fixed point x, estimated by a central difference. |
| `func Orbit(f Map1D, x float64, n int) []float64` | Orbit returns the sequence x, f(x), ..., f^n(x) of length n+1. |
| `func OrbitAfterTransient(f Map1D, x float64, transient, n int) []float64` | OrbitAfterTransient discards the first transient iterates and then returns the next n points of the orbit as a slice of length n. |
| `func PeriodMap(fam Family1D, pMin, pMax float64, steps int, x0 float64, transient, maxPeriod int, tol float64) (params []float64, periods []int)` | PeriodMap returns, for each sampled parameter, the detected period of the attractor of the family fam, as parallel slices of parameter values and… |
| `func Recurrence(a, b Vec, tol float64) bool` | Recurrence reports whether two states are within tolerance tol in Euclidean distance, the basic predicate used in recurrence analysis. |
| `func RecurrenceRate(points []Vec, tol float64) float64` | RecurrenceRate returns the fraction of point pairs in the trajectory that lie within tol of each other (excluding the diagonal), a scalar summary of… |
| `func ReturnMap(crossings []Vec, coord int) (in, out []float64)` | ReturnMap builds a first-return (Poincare) map from the crossings of a section: it returns the pairs (crossings`i`, crossings[i+1]) as two parallel… |
| `func ShannonEntropy(counts []int) float64` | ShannonEntropy returns the Shannon entropy (in nats) of a discrete probability distribution given by counts; the counts are normalised internally. |
| `func SpectralRadius(A Mat) float64` | SpectralRadius returns the largest magnitude among the eigenvalues of a small (2-by-2 or 3-by-3) matrix. |
| `func StdDev(xs []float64) float64` | StdDev returns the population standard deviation of xs. |
| `func SuccessiveMaximaMap(maxima []float64) (in, out []float64)` | SuccessiveMaximaMap turns a sequence of successive maxima into the pairs (m`i`, m[i+1]) of the return map, as parallel slices. |
| `func SumExponents(spectrum Vec) float64` | SumExponents returns the sum of the Lyapunov exponents, which equals the average divergence of the flow (the phase-volume contraction rate). |
| `func SuperstableParameter(period, lo, hi float64, tol float64, maxIter int) (float64, error)` | SuperstableParameter finds the parameter of the logistic family for which the period-2^n superstable cycle occurs (the critical point x=1/2 is… |
| `func SuperstableSequence(count int) []float64` | SuperstableSequence returns the superstable parameters of the logistic family for periods 2^0, 2^1, ..., 2^(count-1), each found by bisection on a… |
| `func Times(h float64, n int) []float64` | Times returns the sample times 0, h, 2h, ..., n*h. |
| `func Variance(xs []float64) float64` | Variance returns the population variance of xs. |
| `func WindingNumber(omega, k, x0 float64, transient, n int) float64` | WindingNumber estimates the winding (rotation) number of the circle map with parameters omega and k, as the average increment of the lifted angle per… |
| `func WrapAngle(x float64) float64` | WrapAngle reduces an angle to the interval [0, 2pi). |
| `func WrapUnit(x float64) float64` | WrapUnit reduces x to the unit interval [0, 1). |

### Types

| Type | What it is |
| --- | --- |
| `BifurcationPoint` | BifurcationPoint holds the parameter value and the sampled attractor points recorded at that parameter in a bifurcation diagram. |
| `Family1D` | Family1D is a one-parameter family of one-dimensional maps. |
| `Field` | Field is the right-hand side of an autonomous ODE system dx/dt = f(x). |
| `Flow` | Flow bundles a vector field with the phase-space dimension it acts on. |
| `Map1D` | Map1D is a one-dimensional real map x_{n+1} = f(x_n). |
| `Map2D` | Map2D is a two-dimensional map returning the next (x, y) from the current pair. |
| `MapN` | MapN is an n-dimensional discrete map x_{k+1} = F(x_k) acting on Vec values. |
| `Mat` | Mat is a dense real matrix stored as a slice of rows; every row must have the same length, which is the number of columns. |
| `Section` | Section describes an oriented hyperplane {x : normal.(x-point)=0} used as a Poincare surface of section. |
| `StabilityKind` | StabilityKind classifies the linear stability of a fixed point or equilibrium. |
| `Vec` | Vec is a real vector represented as a slice of float64 coordinates. |

<details>
<summary><code>BifurcationPoint</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func BifurcationDiagram(fam Family1D, pMin, pMax float64, steps int, x0 float64, transient, keep int) []BifurcationPoint` | BifurcationDiagram samples the attractor of a one-parameter family over the parameter interval [pMin, pMax] at steps parameter values. |

</details>

<details>
<summary><code>Family1D</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func LogisticFamily() Family1D` | LogisticFamily returns the logistic family r -> (x -> r x(1-x)). |
| `func SineFamily() Family1D` | SineFamily returns the sine-map family a -> (x -> a sin(pi x)). |

</details>

<details>
<summary><code>Field</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func DampedPendulum(q, a, wd float64) Field` | DampedPendulum returns the field of a periodically forced, damped pendulum written as an autonomous 3-D system (theta, omega, phase) with damping q,… |
| `func DuffingField(delta, gamma, omega float64) Field` | DuffingField returns the field of the forced Duffing oscillator x” + delta x' - x + x^3 = gamma cos(omega t) as an autonomous 3-D system in (x, v,… |
| `func Lorenz(sigma, rho, beta float64) Field` | Lorenz returns the Lorenz vector field with parameters sigma, rho, beta. |
| `func LorenzStandard() Field` | LorenzStandard returns the Lorenz field with the classic chaotic parameters sigma=10, rho=28, beta=8/3. |
| `func Rossler(a, b, c float64) Field` | Rossler returns the Rossler vector field with parameters a, b, c. |
| `func RosslerStandard() Field` | RosslerStandard returns the Rossler field with the common chaotic parameters a=0.2, b=0.2, c=5.7. |
| `func VanDerPolField(mu float64) Field` | VanDerPolField returns the field of the (unforced) Van der Pol oscillator x” - mu(1-x^2)x' + x = 0 as a planar system in (x, v). |

</details>

<details>
<summary><code>Flow</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewFlow(dim int, f Field) *Flow` | NewFlow constructs a Flow of the given dimension from a vector field. |
| `func (fl *Flow) Eval(x Vec) Vec` | Eval returns the field value F(x). |

</details>

<details>
<summary><code>Map1D</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func BernoulliMap() Map1D` | BernoulliMap returns the doubling (Bernoulli) map x -> (2x) mod 1. |
| `func CircleMap(omega, k float64) Map1D` | CircleMap returns the (uncoupled) circle map theta -> theta + Omega - (K/2pi) sin(2pi theta), reduced modulo 1. |
| `func Compose(f Map1D, n int) Map1D` | Compose returns the n-fold composition f^n as a single Map1D. |
| `func CubicMap(a float64) Map1D` | CubicMap returns the cubic map x -> a*x^3 + (1-a)*x, a symmetric map with a period-doubling route to chaos. |
| `func CubicMapDeriv(a float64) Map1D` | CubicMapDeriv returns the derivative of the cubic map. |
| `func GaussMap(alpha, beta float64) Map1D` | GaussMap returns the Gaussian (mouse) map x -> exp(-alpha*x^2) + beta. |
| `func GaussMapDeriv(alpha, beta float64) Map1D` | GaussMapDeriv returns the derivative of the Gaussian map. |
| `func Logistic(r float64) Map1D` | Logistic returns the logistic map x -> r*x*(1-x). |
| `func LogisticDeriv(r float64) Map1D` | LogisticDeriv returns the derivative of the logistic map, x -> r*(1-2x). |
| `func SineMap(a float64) Map1D` | SineMap returns the sine map x -> a*sin(pi*x), a common smooth analogue of the logistic map. |
| `func SineMapDeriv(a float64) Map1D` | SineMapDeriv returns the derivative of the sine map. |
| `func Tent(mu float64) Map1D` | Tent returns the tent map with slope mu on [0,1]: 2mu*x for x<1/2 and 2mu*(1-x) otherwise. |
| `func TentDeriv(mu float64) Map1D` | TentDeriv returns the (piecewise constant) derivative of the tent map. |

</details>

<details>
<summary><code>Map2D</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func GingerbreadmanMap() Map2D` | GingerbreadmanMap returns the area-preserving Gingerbreadman map (x,y) -> (1 - y + \|x\|, x). |
| `func HenonMap(a, b float64) Map2D` | HenonMap returns the Henon map (x,y) -> (1 - a*x^2 + y, b*x). |
| `func IkedaMap(u float64) Map2D` | IkedaMap returns the Ikeda map, a dissipative planar map with a spiral attractor, using the standard parameterisation with parameter u. |
| `func StandardMap(k float64) Map2D` | StandardMap returns the Chirikov standard map on the cylinder, p' = p + K sin(theta), theta' = theta + p', with theta reduced modulo 2pi. |
| `func TinkerbellMap(a, b, c, d float64) Map2D` | TinkerbellMap returns the Tinkerbell map with parameters a, b, c, d. |

</details>

<details>
<summary><code>MapN</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Lift2D(f Map2D) MapN` | Lift2D turns a Map2D into a MapN acting on 2-vectors. |

</details>

<details>
<summary><code>Mat</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Eye(n int) Mat` | Eye returns the n-by-n identity matrix. |
| `func HenonJacobian(a, b, x, y float64) Mat` | HenonJacobian returns the Jacobian matrix of the Henon map at (x, y). |
| `func Inverse(A Mat) (Mat, error)` | Inverse returns the inverse of a square matrix, or ErrSingular if A is not invertible to working precision. |
| `func JacobianField(f Field, x Vec, eps float64) Mat` | JacobianField returns the numerical Jacobian of the vector field f at x, computed by central differences with step eps. |
| `func JacobianMap(F MapN, x Vec, eps float64) Mat` | JacobianMap returns the numerical Jacobian of the n-dimensional map F at the point x, computed by central differences with step eps. |
| `func LorenzJacobian(sigma, rho, beta float64, v Vec) Mat` | LorenzJacobian returns the Jacobian of the Lorenz field at state v. |
| `func MatFromRows(rows ...Vec) Mat` | MatFromRows builds a matrix from the given rows, copying each one. |
| `func NewMat(r, c int) Mat` | NewMat returns an r-by-c matrix of zeros. |
| `func QR(A Mat) (Q, R Mat)` | QR computes the reduced QR factorisation of A by modified Gram-Schmidt, returning an orthonormal Q with the same shape as A and an upper-triangular R… |
| `func RosslerJacobian(a, _, c float64, v Vec) Mat` | RosslerJacobian returns the Jacobian of the Rossler field at state v. |
| `func StandardMapJacobian(k, theta, _ float64) Mat` | StandardMapJacobian returns the Jacobian of the standard map at (theta, p). |
| `func (A Mat) Add(B Mat) Mat` | Add returns the elementwise sum A+B. |
| `func (A Mat) Clone() Mat` | Clone returns an independent copy of A. |
| `func (A Mat) Col(j int) Vec` | Col returns column j of A as a new vector. |
| `func (A Mat) Cols() int` | Cols returns the number of columns of A (zero for an empty matrix). |
| `func (A Mat) FrobeniusNorm() float64` | FrobeniusNorm returns the Frobenius norm of A. |
| `func (A Mat) IsSquare() bool` | IsSquare reports whether A has equal numbers of rows and columns. |
| `func (A Mat) Mul(B Mat) Mat` | Mul returns the matrix product A*B. |
| `func (A Mat) MulVec(x Vec) Vec` | MulVec returns the matrix-vector product A*x. |
| `func (A Mat) Rows() int` | Rows returns the number of rows of A. |
| `func (A Mat) Scale(s float64) Mat` | Scale returns the matrix s*A. |
| `func (A Mat) Sub(B Mat) Mat` | Sub returns the elementwise difference A-B. |
| `func (A Mat) Trace() float64` | Trace returns the sum of the diagonal entries of a square matrix. |
| `func (A Mat) Transpose() Mat` | Transpose returns the transpose of A. |

</details>

<details>
<summary><code>Section</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func PlaneSection(dim, coord int, value float64, direction int) Section` | PlaneSection returns a Section for the coordinate hyperplane x`coord`=value, keeping crossings in the given direction. |

</details>

<details>
<summary><code>StabilityKind</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func ClassifyFlow(J Mat) StabilityKind` | ClassifyFlow classifies the equilibrium of a continuous system from the eigenvalues of its Jacobian (2-by-2 or 3-by-3), by the sign of the real parts… |
| `func ClassifyMap(J Mat) StabilityKind` | ClassifyMap classifies the fixed point of a discrete system from the eigenvalues (multipliers) of its Jacobian, by their magnitude relative to the… |
| `func (k StabilityKind) String() string` | String returns a human-readable name for the stability kind. |

</details>

<details>
<summary><code>Vec</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Equilibrium(f Field, x0 Vec, tol float64, maxIter int) (Vec, error)` | Equilibrium refines an equilibrium of the vector field f (a root of f(x)=0) near x0 by Newton's method with a numerical Jacobian. |
| `func FixedPointMap(F MapN, x0 Vec, tol float64, maxIter int) (Vec, error)` | FixedPointMap refines a fixed point of the n-dimensional map F near x0 by Newton's method on G(x)=F(x)-x with a numerical Jacobian. |
| `func GramSchmidt(vs []Vec) []Vec` | GramSchmidt orthonormalises the given set of vectors in place-independent fashion, returning a new orthonormal basis for their span (with zero… |
| `func Integrate(f Field, x Vec, h float64, n int) Vec` | Integrate advances x for n steps of size h with RK4 and returns the final state. |
| `func IterateN(F MapN, v Vec, n int) Vec` | IterateN applies the map F to v n times, returning the final vector. |
| `func LorenzEquilibria(rho, beta float64) []Vec` | LorenzEquilibria returns the analytic equilibria of the Lorenz system with parameters rho and beta: the origin, and (for rho>1) the symmetric pair… |
| `func LyapunovSpectrumFlow(f Field, x0 Vec, h float64, transient, n, reortho int) Vec` | LyapunovSpectrumFlow estimates the full Lyapunov spectrum of the flow with field f by the QR method: an orthonormal frame is evolved with the… |
| `func LyapunovSpectrumMap(F MapN, x0 Vec, transient, n int) Vec` | LyapunovSpectrumMap estimates the full Lyapunov spectrum of the map F by the QR (Benettin) method: an orthonormal frame is evolved by the Jacobian… |
| `func NewVec(n int) Vec` | NewVec returns a zero vector of length n. |
| `func OrbitN(F MapN, v Vec, n int) []Vec` | OrbitN returns the sequence v, F(v), ..., F^n(v) as n+1 vectors. |
| `func OrbitNAfterTransient(F MapN, v Vec, transient, n int) []Vec` | OrbitNAfterTransient discards the first transient iterates and returns the next n points of the orbit. |
| `func PoincareSection(f Field, x0 Vec, h float64, transient, n int, s Section) ([]Vec, error)` | PoincareSection integrates the flow with field f from x0 for n steps of size h and returns the sequence of intersection points with the section,… |
| `func RandomInitialConditions(rng *rand.Rand, m, dim int, lo, hi float64) []Vec` | RandomInitialConditions returns m random points uniformly distributed in the box [lo, hi]^dim, drawn from the supplied random source for… |
| `func SolveLinear(A Mat, b Vec) (Vec, error)` | SolveLinear solves A*x = b for a square matrix A by Gaussian elimination with partial pivoting. |
| `func StepEuler(f Field, x Vec, h float64) Vec` | StepEuler advances the state x by one explicit-Euler step of size h. |
| `func StepRK2(f Field, x Vec, h float64) Vec` | StepRK2 advances the state x by one midpoint (second-order Runge-Kutta) step of size h. |
| `func StepRK4(f Field, x Vec, h float64) Vec` | StepRK4 advances the state x by one classical fourth-order Runge-Kutta step of size h. |
| `func StroboscopicMap(f Field, x0 Vec, h float64, stepsPerPeriod, transientPeriods, periods int) []Vec` | StroboscopicMap samples the trajectory of a periodically forced flow at times that are multiples of the forcing period, returning one point per… |
| `func TakensEmbedding(series []float64, m, tau int) []Vec` | TakensEmbedding reconstructs a delay-coordinate embedding of a scalar time series with embedding dimension m and delay tau, returning the sequence of… |
| `func Trajectory(f Field, x Vec, h float64, n int) []Vec` | Trajectory integrates x for n steps of size h with RK4 and returns the n+1 states x(0), x(h), ..., x(nh). |
| `func TrajectoryAfterTransient(f Field, x Vec, h float64, transient, n int) []Vec` | TrajectoryAfterTransient integrates away the first transient steps and then returns the next n+1 recorded states. |
| `func VecOf(xs ...float64) Vec` | VecOf returns a Vec containing the given components. |
| `func (v Vec) Add(w Vec) Vec` | Add returns the elementwise sum v+w. |
| `func (v Vec) AddScaled(s float64, w Vec) Vec` | AddScaled returns v + s*w, the fused multiply-add of two vectors. |
| `func (v Vec) Clone() Vec` | Clone returns an independent copy of v. |
| `func (v Vec) Dim() int` | Dim returns the number of components of v. |
| `func (v Vec) Distance(w Vec) float64` | Distance returns the Euclidean distance between v and w. |
| `func (v Vec) Dot(w Vec) float64` | Dot returns the Euclidean inner product of v and w. |
| `func (v Vec) Max() float64` | Max returns the largest component of v, or NaN if v is empty. |
| `func (v Vec) Mean() float64` | Mean returns the arithmetic mean of the components of v, or NaN if empty. |
| `func (v Vec) Min() float64` | Min returns the smallest component of v, or NaN if v is empty. |
| `func (v Vec) Norm() float64` | Norm returns the Euclidean (L2) norm of v. |
| `func (v Vec) Norm1() float64` | Norm1 returns the L1 (taxicab) norm of v. |
| `func (v Vec) NormInf() float64` | NormInf returns the maximum-absolute-value (Chebyshev) norm of v. |
| `func (v Vec) Normalize() Vec` | Normalize returns v scaled to unit Euclidean norm; a zero vector is returned unchanged. |
| `func (v Vec) Scale(s float64) Vec` | Scale returns the vector s*v. |
| `func (v Vec) Sub(w Vec) Vec` | Sub returns the elementwise difference v-w. |
| `func (v Vec) Sum() float64` | Sum returns the sum of the components of v. |

</details>

### Constants

`FeigenbaumAlpha`, `FeigenbaumDelta`, `LogisticAccumulation`

### Variables

`ErrDimensionMismatch`, `ErrEmpty`, `ErrNonSquare`, `ErrSingular`, `ErrNoConvergence`, `ErrBadParameter`, `ErrNoCrossing`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/chaos).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
