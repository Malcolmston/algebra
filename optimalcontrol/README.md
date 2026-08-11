# optimalcontrol — implements optimal control and dynamic programming algorithms using only the Go

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/optimalcontrol.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/optimalcontrol)

Package optimalcontrol implements optimal control and dynamic programming
algorithms using only the Go standard library.

The package is self-contained: it ships its own small dense linear-algebra
layer (the Matrix type with LU, Cholesky, least-squares, matrix exponential,
Faddeev–LeVerrier characteristic polynomials and Durand–Kerner
eigenvalues, symmetric Jacobi eigensolves, and Lyapunov/Sylvester solvers) and
builds the control-theoretic algorithms on top of it.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/optimalcontrol"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func CharPoly(a *Matrix) []float64` | CharPoly returns the coefficients of the characteristic polynomial det(λI − A) of a square matrix, in ascending powers of λ (index i holds the… |
| `func ClosedLoopEigenvaluesContinuous(a, b, k *Matrix) []complex128` | ClosedLoopEigenvaluesContinuous returns the eigenvalues of A − B K, the continuous closed-loop dynamics under state feedback K. |
| `func ClosedLoopEigenvaluesDiscrete(a, b, k *Matrix) []complex128` | ClosedLoopEigenvaluesDiscrete returns the eigenvalues of A − B K, the discrete closed-loop dynamics under state feedback K. |
| `func ControllabilityRank(a, b *Matrix) int` | ControllabilityRank returns the rank of the controllability matrix of (A, B). |
| `func CostateDynamicsLQ(a, q *Matrix, x, p []float64) []float64` | CostateDynamicsLQ returns the costate derivative p' = −Q x − Aᵀ p for the linear-quadratic problem. |
| `func Det(a *Matrix) (float64, error)` | Det returns the determinant of a square matrix. |
| `func Eigenvalues(a *Matrix) []complex128` | Eigenvalues returns all (generally complex) eigenvalues of a square matrix by forming its characteristic polynomial and finding its roots with the… |
| `func IsControllable(a, b *Matrix) bool` | IsControllable reports whether the pair (A, B) is controllable, i.e. |
| `func IsDetectableContinuous(a, c *Matrix) bool` | IsDetectableContinuous reports whether the pair (A, C) is detectable in continuous time (dual of stabilizability). |
| `func IsDetectableDiscrete(a, c *Matrix) bool` | IsDetectableDiscrete reports whether the pair (A, C) is detectable in discrete time. |
| `func IsObservable(a, c *Matrix) bool` | IsObservable reports whether the pair (A, C) is observable, i.e. |
| `func IsPositiveDefinite(a *Matrix) bool` | IsPositiveDefinite reports whether a symmetric matrix is positive definite. |
| `func IsPositiveSemiDefinite(a *Matrix, tol float64) bool` | IsPositiveSemiDefinite reports whether a symmetric matrix has all eigenvalues >= -tol. |
| `func IsStabilizableContinuous(a, b *Matrix) bool` | IsStabilizableContinuous reports whether the pair (A, B) is stabilizable in continuous time: every eigenvalue of A with non-negative real part is… |
| `func IsStabilizableDiscrete(a, b *Matrix) bool` | IsStabilizableDiscrete reports whether the pair (A, B) is stabilizable in discrete time: every eigenvalue of A with modulus >= 1 is controllable. |
| `func IsStableContinuous(a *Matrix, tol float64) bool` | IsStableContinuous reports whether every eigenvalue of A has real part below -tol (Hurwitz stability of x' = A x). |
| `func IsStableDiscrete(a *Matrix, tol float64) bool` | IsStableDiscrete reports whether every eigenvalue of A has modulus below 1-tol (Schur stability of x_{k+1} = A x_k). |
| `func LQRCostToGo(p *Matrix, x []float64) float64` | LQRCostToGo returns the optimal quadratic cost xᵀ P x of steering the state x to the origin under an infinite-horizon LQR with Riccati solution P. |
| `func LeastSquares(a *Matrix, b []float64) ([]float64, error)` | LeastSquares solves the overdetermined system A x = b in the least-squares sense via the normal equations AᵀA x = Aᵀb. |
| `func LinearInterp1D(grid, values []float64, x float64) float64` | LinearInterp1D linearly interpolates the tabulated function (grid, values) at x. |
| `func ObservabilityRank(a, c *Matrix) int` | ObservabilityRank returns the rank of the observability matrix of (A, C). |
| `func OptimalControlLQ(b, r *Matrix, p []float64) ([]float64, error)` | OptimalControlLQ returns the stationarity-condition control u = −R⁻¹ Bᵀ p that minimizes the linear-quadratic Hamiltonian for a given costate… |
| `func PolyRootsDK(coeffs []float64) []complex128` | PolyRootsDK finds all complex roots of a real polynomial using the Durand–Kerner (Weierstrass) iteration. |
| `func QuadraticCostDiscrete(xs [][]float64, us [][]float64, q, r, qf *Matrix) float64` | QuadraticCostDiscrete evaluates the discrete quadratic cost Σ_k (xₖᵀ Q xₖ + uₖᵀ R uₖ) + x_Nᵀ Qf x_N over a state trajectory and… |
| `func Rank(a *Matrix, tol float64) int` | Rank returns the numerical rank of a matrix using Gaussian elimination with partial pivoting and the supplied tolerance for treating pivots as zero. |
| `func SimulateDiscreteLQR(a, b *Matrix, fh *FiniteHorizonDiscrete, x0 []float64) [][]float64` | SimulateDiscreteLQR simulates the closed-loop discrete system under the time-varying gains produced by FiniteHorizonLQRDiscrete, starting from x0 and… |
| `func Solve(a *Matrix, b []float64) ([]float64, error)` | Solve solves the linear system A x = b. |
| `func SpectralAbscissa(a *Matrix) float64` | SpectralAbscissa returns the maximum real part of the eigenvalues of a matrix. |
| `func SpectralRadius(a *Matrix) float64` | SpectralRadius returns the maximum modulus of the eigenvalues of a matrix. |
| `func SymEigenvalues(a *Matrix) ([]float64, error)` | SymEigenvalues returns the eigenvalues of a symmetric matrix in ascending order. |
| `func VecAdd(a, b []float64) []float64` | VecAdd returns the element-wise sum a+b. |
| `func VecAxpy(s float64, b, a []float64) []float64` | VecAxpy returns a + s·b (the classic "axpy" operation). |
| `func VecCopy(a []float64) []float64` | VecCopy returns a copy of a. |
| `func VecDot(a, b []float64) float64` | VecDot returns the inner product aᵀb. |
| `func VecMaxAbs(a []float64) float64` | VecMaxAbs returns the largest absolute component of a. |
| `func VecNorm(a []float64) float64` | VecNorm returns the Euclidean norm of a. |
| `func VecScale(a []float64, s float64) []float64` | VecScale returns the vector a scaled by s. |
| `func VecSub(a, b []float64) []float64` | VecSub returns the element-wise difference a−b. |

### Types

| Type | What it is |
| --- | --- |
| `FiniteHorizonContinuous` | FiniteHorizonContinuous holds the sampled solution of the continuous finite-horizon Riccati differential equation on a time grid. |
| `FiniteHorizonDiscrete` | FiniteHorizonDiscrete holds the time-varying gains and cost matrices of a finite-horizon discrete LQR problem, indexed from stage 0 (initial) to N… |
| `HJBGrid1D` | HJBGrid1D describes a one-dimensional infinite-horizon optimal-control problem discretized for solution by value iteration (a semi-Lagrangian scheme… |
| `HJBResult` | HJBResult holds the value function and greedy control at each grid point. |
| `KalmanContinuousResult` | KalmanContinuousResult holds the steady-state error covariance and estimator gain of a continuous-time Kalman–Bucy filter. |
| `KalmanDiscreteResult` | KalmanDiscreteResult holds the steady-state a-priori covariance and gain of a discrete-time Kalman filter. |
| `KalmanFilter` | KalmanFilter is a recursive discrete-time Kalman filter that tracks a state estimate and its covariance. |
| `LQGResult` | LQGResult bundles the regulator and estimator designs of a linear-quadratic Gaussian controller. |
| `LQRResult` | LQRResult bundles the gain matrix, Riccati solution and closed-loop matrix of a linear-quadratic regulator design. |
| `LU` | LU holds the result of an LU factorization with partial pivoting: a combined lower/upper factor lu, the row permutation piv, and the sign of the… |
| `LinearSystem` | LinearSystem is a continuous- or discrete-time linear time-invariant model x' = A x + B u (continuous) or x_{k+1} = A x_k + B u_k (discrete), with… |
| `MDP` | MDP is a finite, discounted Markov decision process with the standard reward-maximization convention. |
| `Matrix` | Matrix is a dense real matrix stored in row-major order. |
| `PolicyIterationResult` | PolicyIterationResult holds the outcome of policy iteration. |
| `TPBVPSolution` | TPBVPSolution holds the sampled solution of a two-point boundary-value problem arising from Pontryagin's minimum principle: the state, costate and… |
| `ValueIterationResult` | ValueIterationResult holds the outcome of value iteration. |

<details>
<summary><code>FiniteHorizonContinuous</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func FiniteHorizonLQRContinuous(a, b, q, r, qf *Matrix, tFinal float64, steps int) (*FiniteHorizonContinuous, error)` | FiniteHorizonLQRContinuous integrates the matrix Riccati differential equation |
| `func (fh *FiniteHorizonContinuous) GainAt(b, r *Matrix, i int) (*Matrix, error)` | GainAt returns the continuous LQR feedback gain K(t) = R⁻¹ Bᵀ P(t) at the grid index i of a FiniteHorizonContinuous solution. |

</details>

<details>
<summary><code>FiniteHorizonDiscrete</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func FiniteHorizonLQRDiscrete(a, b, q, r, qf *Matrix, n int) (*FiniteHorizonDiscrete, error)` | FiniteHorizonLQRDiscrete solves the finite-horizon discrete LQR problem for x_{k+1} = A x_k + B u_k minimizing xₙᵀ Qf xₙ + Σ_{k<N} (xₖᵀ Q… |

</details>

<details>
<summary><code>HJBGrid1D</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (h *HJBGrid1D) ControlAt(res *HJBResult, x float64) float64` | ControlAt returns the greedy control at state x by interpolating the tabulated policy of a solved HJBResult onto the grid. |
| `func (h *HJBGrid1D) Solve(tol float64, maxIter int) *HJBResult` | Solve runs value iteration for the infinite-horizon discounted problem until the max-norm change is below tol or maxIter sweeps elapse. |
| `func (h *HJBGrid1D) SolveFiniteHorizon(stages int, terminal func(x float64) float64) *HJBResult` | SolveFiniteHorizon solves the finite-horizon problem with terminal cost g(x) by backward value iteration over the given number of stages, returning… |

</details>

<details>
<summary><code>KalmanContinuousResult</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func KalmanContinuous(a, c, w, v *Matrix) (*KalmanContinuousResult, error)` | KalmanContinuous designs the steady-state continuous-time Kalman–Bucy filter for x' = A x + w, y = C x + v with process covariance W and… |

</details>

<details>
<summary><code>KalmanDiscreteResult</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func KalmanDiscrete(a, c, w, v *Matrix) (*KalmanDiscreteResult, error)` | KalmanDiscrete designs the steady-state discrete-time Kalman filter for x_{k+1} = A x_k + w_k, y_k = C x_k + v_k with process covariance W and… |

</details>

<details>
<summary><code>KalmanFilter</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewKalmanFilter(a, b, c, q, r *Matrix, x0 []float64, p0 *Matrix) *KalmanFilter` | NewKalmanFilter constructs a Kalman filter with the given model matrices and initial estimate x0 with covariance p0. |
| `func (kf *KalmanFilter) Predict(u []float64)` | Predict advances the estimate through the process model using control input u (which may be nil when B is nil): x̂ ← A x̂ + B u, P ← A P Aᵀ +… |
| `func (kf *KalmanFilter) Update(y []float64) error` | Update corrects the prediction with measurement y using the Kalman gain: L = P Cᵀ (C P Cᵀ + R)⁻¹, x̂ ← x̂ + L (y − C x̂), P ← (I −… |

</details>

<details>
<summary><code>LQGResult</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func LQGContinuous(a, b, c, q, r, w, v *Matrix) (*LQGResult, error)` | LQGContinuous designs a continuous-time LQG controller by combining a continuous LQR (weights Q, R on the plant (A, B)) with a continuous Kalman… |
| `func LQGDiscrete(a, b, c, q, r, w, v *Matrix) (*LQGResult, error)` | LQGDiscrete designs a discrete-time LQG controller by combining a discrete LQR with a discrete Kalman filter via the separation principle. |

</details>

<details>
<summary><code>LQRResult</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func LQRContinuous(a, b, q, r *Matrix) (*LQRResult, error)` | LQRContinuous designs a continuous-time infinite-horizon LQR for the system x' = A x + B u minimizing ∫ (xᵀ Q x + uᵀ R u) dt. |
| `func LQRDiscrete(a, b, q, r *Matrix) (*LQRResult, error)` | LQRDiscrete designs a discrete-time infinite-horizon LQR for the system x_{k+1} = A x_k + B u_k minimizing Σ (xᵀ Q x + uᵀ R u). |

</details>

<details>
<summary><code>LU</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Factor(a *Matrix) (*LU, error)` | Factor computes the LU decomposition of a square matrix with partial pivoting. |
| `func (f *LU) Det() float64` | Det returns the determinant computed from the factorization. |
| `func (f *LU) Solve(b []float64) ([]float64, error)` | Solve solves A x = b for x using the factorization. |
| `func (f *LU) SolveMatrix(b *Matrix) (*Matrix, error)` | SolveMatrix solves A X = B, treating each column of B independently. |

</details>

<details>
<summary><code>LinearSystem</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewLinearSystem(a, b, c, d *Matrix) *LinearSystem` | NewLinearSystem builds a LinearSystem, filling C with the identity and D with zeros when they are nil. |
| `func (s *LinearSystem) Output(x, u []float64) []float64` | Output returns y = C x + D u for the system. |
| `func (s *LinearSystem) SimulateContinuous(x0 []float64, control func(t float64) []float64, dt float64, steps int) [][]float64` | SimulateContinuous integrates the continuous-time system with a fixed control u(t) held over each step of size dt using classical RK4, returning the… |
| `func (s *LinearSystem) SimulateDiscrete(x0 []float64, us [][]float64) [][]float64` | SimulateDiscrete simulates the discrete-time system for the given control sequence, returning the state trajectory x[0..len(us)]. |
| `func (s *LinearSystem) StepDiscrete(x, u []float64) []float64` | StepDiscrete advances the discrete-time state one step: x_{k+1} = A x + B u. |

</details>

<details>
<summary><code>MDP</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewMDP(trans []*Matrix, reward *Matrix, gamma float64) (*MDP, error)` | NewMDP constructs an MDP from per-action transition matrices, a reward matrix and a discount factor, validating the dimensions. |
| `func (m *MDP) BellmanBackup(v []float64) (newV []float64, policy []int)` | BellmanBackup applies one Bellman optimality update, returning the improved value function and the greedy policy that attains it. |
| `func (m *MDP) GaussSeidelValueIteration(tol float64, maxIter int) *ValueIterationResult` | GaussSeidelValueIteration runs value iteration using in-place (Gauss–Seidel) updates, which typically converges in fewer sweeps than the Jacobi… |
| `func (m *MDP) GreedyPolicy(v []float64) []int` | GreedyPolicy returns the greedy policy with respect to a value function. |
| `func (m *MDP) ModifiedPolicyIteration(k, maxIter int, tol float64) (*PolicyIterationResult, error)` | ModifiedPolicyIteration runs modified policy iteration, evaluating each policy with a fixed number k of Bellman expectation sweeps rather than an… |
| `func (m *MDP) PolicyEvaluationExact(policy []int) ([]float64, error)` | PolicyEvaluationExact computes the exact value function of a deterministic policy by solving the linear system (I − γ P_π) V = r_π. |
| `func (m *MDP) PolicyEvaluationIterative(policy []int, tol float64, maxIter int) []float64` | PolicyEvaluationIterative computes the value function of a policy by iterative application of the Bellman expectation operator. |
| `func (m *MDP) PolicyIteration(maxIter int) (*PolicyIterationResult, error)` | PolicyIteration runs Howard's policy iteration with exact policy evaluation, returning the optimal policy and its value function. |
| `func (m *MDP) PolicyMatrices(policy []int) (*Matrix, []float64)` | PolicyMatrices returns the transition matrix P_π and reward vector r_π induced by a deterministic policy. |
| `func (m *MDP) QValues(v []float64) *Matrix` | QValues returns the States×Actions matrix of action values Q(s, a) = R(s, a) + γ Σ_{s'} P_a(s, s') V(s') for a given value function. |
| `func (m *MDP) ValueIteration(tol float64, maxIter int) *ValueIterationResult` | ValueIteration runs value iteration until the max-norm change between sweeps falls below tol or maxIter sweeps have elapsed. |

</details>

<details>
<summary><code>Matrix</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func BlockMatrix(a, b, c, d *Matrix) *Matrix` | BlockMatrix assembles the 2×2 block matrix [[a b];[c d]]. |
| `func CAREResidual(a, b, q, r, x *Matrix) (*Matrix, error)` | CAREResidual returns Aᵀ X + X A − X B R⁻¹ Bᵀ X + Q, the residual of the continuous algebraic Riccati equation for a candidate solution X. |
| `func Cholesky(a *Matrix) (*Matrix, error)` | Cholesky computes the lower-triangular Cholesky factor L with A = L Lᵀ for a symmetric positive-definite matrix. |
| `func ColumnVector(v []float64) *Matrix` | ColumnVector builds an n×1 matrix from v. |
| `func ContinuousGain(b, r, p *Matrix) (*Matrix, error)` | ContinuousGain returns the continuous LQR feedback gain K = R⁻¹ Bᵀ P. |
| `func ControllabilityGramianContinuous(a, b *Matrix) (*Matrix, error)` | ControllabilityGramianContinuous returns the infinite-horizon controllability Gramian Wc solving A Wc + Wc Aᵀ + B Bᵀ = 0 for a stable A. |
| `func ControllabilityGramianDiscrete(a, b *Matrix) (*Matrix, error)` | ControllabilityGramianDiscrete returns the discrete controllability Gramian solving A Wc Aᵀ − Wc + B Bᵀ = 0 for a Schur-stable A. |
| `func ControllabilityMatrix(a, b *Matrix) *Matrix` | ControllabilityMatrix returns the controllability matrix [B, AB, A²B, …, Aⁿ⁻¹B] for the pair (A, B). |
| `func DAREResidual(a, b, q, r, x *Matrix) (*Matrix, error)` | DAREResidual returns the residual X − (Aᵀ X A − Aᵀ X B (R+BᵀXB)⁻¹ BᵀXA + Q) of the discrete algebraic Riccati equation for a candidate… |
| `func Diag(v []float64) *Matrix` | Diag returns a square diagonal matrix whose diagonal is v. |
| `func DiscreteGain(a, b, r, p *Matrix) (*Matrix, error)` | DiscreteGain returns the discrete LQR feedback gain K = (R + Bᵀ P B)⁻¹ Bᵀ P A for a given Riccati solution P. |
| `func DiscretizeZOH(a, b *Matrix, dt float64) (ad, bd *Matrix)` | DiscretizeZOH returns the exact zero-order-hold discretization (Ad, Bd) of the continuous pair (A, B) for sample time dt, computed from the matrix… |
| `func Eye(n int) *Matrix` | Eye returns the n×n identity matrix. |
| `func FromRows(rows [][]float64) *Matrix` | FromRows builds a matrix from a slice of equal-length rows. |
| `func HStack(a, b *Matrix) *Matrix` | HStack returns the horizontal concatenation [a b]. |
| `func HamiltonianMatrix(a, b, q, r *Matrix) (*Matrix, error)` | HamiltonianMatrix returns the 2n×2n Hamiltonian matrix |
| `func Inverse(a *Matrix) (*Matrix, error)` | Inverse returns the inverse of a square matrix. |
| `func Kron(a, b *Matrix) *Matrix` | Kron returns the Kronecker product a⊗b. |
| `func LeastSquaresMatrix(a, b *Matrix) (*Matrix, error)` | LeastSquaresMatrix solves A X = B in the least-squares sense column by column. |
| `func LyapunovContinuousResidual(a, q, x *Matrix) *Matrix` | LyapunovContinuousResidual returns Aᵀ X + X A + Q, the residual of the continuous Lyapunov equation, useful for verifying a solution. |
| `func LyapunovDiscreteResidual(a, q, x *Matrix) *Matrix` | LyapunovDiscreteResidual returns Aᵀ X A − X + Q, the residual of the discrete Lyapunov equation. |
| `func MatrixExp(a *Matrix) *Matrix` | MatrixExp returns the matrix exponential exp(A) using a scaling-and-squaring scheme with a truncated Taylor series. |
| `func MatrixPow(a *Matrix, p int) *Matrix` | MatrixPow returns A raised to the non-negative integer power p. |
| `func MatrixSign(a *Matrix, maxIter int, tol float64) (*Matrix, error)` | MatrixSign computes the matrix sign function of a square matrix with no purely imaginary eigenvalues, using the Newton iteration with determinantal… |
| `func NewMatrix(r, c int, data []float64) *Matrix` | NewMatrix builds an r×c matrix from the supplied row-major data. |
| `func ObservabilityGramianContinuous(a, c *Matrix) (*Matrix, error)` | ObservabilityGramianContinuous returns the infinite-horizon observability Gramian Wo solving Aᵀ Wo + Wo A + Cᵀ C = 0 for a stable A. |
| `func ObservabilityGramianDiscrete(a, c *Matrix) (*Matrix, error)` | ObservabilityGramianDiscrete returns the discrete observability Gramian solving Aᵀ Wo A − Wo + Cᵀ C = 0 for a Schur-stable A. |
| `func ObservabilityMatrix(a, c *Matrix) *Matrix` | ObservabilityMatrix returns the observability matrix [C; CA; CA²; …; CAⁿ⁻¹] for the pair (A, C). |
| `func Ones(r, c int) *Matrix` | Ones returns an r×c matrix whose entries are all one. |
| `func RowVector(v []float64) *Matrix` | RowVector builds a 1×n matrix from v. |
| `func SolveCARE(a, b, q, r *Matrix) (*Matrix, error)` | SolveCARE solves the continuous-time algebraic Riccati equation for its stabilizing solution. |
| `func SolveCAREKleinman(a, b, q, r, k0 *Matrix, maxIter int, tol float64) (*Matrix, error)` | SolveCAREKleinman solves the continuous-time algebraic Riccati equation by Kleinman's Newton iteration starting from a stabilizing gain k0 (so that A… |
| `func SolveCARESign(a, b, q, r *Matrix) (*Matrix, error)` | SolveCARESign solves the continuous-time algebraic Riccati equation |
| `func SolveDARE(a, b, q, r *Matrix) (*Matrix, error)` | SolveDARE solves the discrete-time algebraic Riccati equation for its stabilizing solution. |
| `func SolveDAREIter(a, b, q, r *Matrix, maxIter int, tol float64) (*Matrix, error)` | SolveDAREIter solves the discrete-time algebraic Riccati equation |
| `func SolveLyapunovContinuous(a, q *Matrix) (*Matrix, error)` | SolveLyapunovContinuous solves the continuous-time Lyapunov equation Aᵀ X + X A + Q = 0 for the symmetric matrix X. |
| `func SolveLyapunovDiscrete(a, q *Matrix) (*Matrix, error)` | SolveLyapunovDiscrete solves the discrete-time (Stein) Lyapunov equation Aᵀ X A − X + Q = 0 for the symmetric matrix X. |
| `func SolveMatrix(a, b *Matrix) (*Matrix, error)` | SolveMatrix solves A X = B for the matrix X. |
| `func SolveSylvester(a, b, c *Matrix) (*Matrix, error)` | SolveSylvester solves the Sylvester equation A X + X B = C for X, where A is m×m, B is n×n and C is m×n. |
| `func SymEigen(a *Matrix) (w []float64, v *Matrix, err error)` | SymEigen computes eigenvalues and orthonormal eigenvectors of a symmetric matrix using the cyclic Jacobi method. |
| `func Unvec(v []float64, r, c int) *Matrix` | Unvec reshapes a column-major vector into an r×c matrix. |
| `func VStack(a, b *Matrix) *Matrix` | VStack returns the vertical concatenation [a; b]. |
| `func Zeros(r, c int) *Matrix` | Zeros returns an r×c matrix of zeros. |
| `func (m *Matrix) Add(i, j int, v float64)` | Add accumulates v into the element at row i, column j. |
| `func (m *Matrix) ApproxEqual(b *Matrix, tol float64) bool` | ApproxEqual reports whether a and b have the same shape and entries within absolute tolerance tol. |
| `func (m *Matrix) At(i, j int) float64` | At returns the element at row i, column j. |
| `func (m *Matrix) Clone() *Matrix` | Clone returns a deep copy of the matrix. |
| `func (m *Matrix) Col(j int) []float64` | Col returns a copy of column j as a slice. |
| `func (m *Matrix) Cols() int` | Cols returns the number of columns. |
| `func (m *Matrix) Data() []float64` | Data returns a copy of the underlying row-major data. |
| `func (m *Matrix) Equal(b *Matrix) bool` | Equal reports whether a and b have the same shape and identical entries. |
| `func (m *Matrix) FrobeniusNorm() float64` | FrobeniusNorm returns the Frobenius norm of the matrix. |
| `func (m *Matrix) InfNorm() float64` | InfNorm returns the maximum absolute row sum of the matrix. |
| `func (m *Matrix) IsSquare() bool` | IsSquare reports whether the matrix is square. |
| `func (m *Matrix) IsSymmetric(tol float64) bool` | IsSymmetric reports whether the matrix is square and symmetric within tol. |
| `func (m *Matrix) MaxAbs() float64` | MaxAbs returns the largest absolute entry of the matrix. |
| `func (m *Matrix) Minus(b *Matrix) *Matrix` | Minus returns the difference a-b. |
| `func (m *Matrix) Mul(b *Matrix) *Matrix` | Mul returns the matrix product a·b. |
| `func (m *Matrix) MulVec(v []float64) []float64` | MulVec returns the matrix-vector product a·v. |
| `func (m *Matrix) Neg() *Matrix` | Neg returns the additive inverse of the matrix. |
| `func (m *Matrix) OneNorm() float64` | OneNorm returns the maximum absolute column sum of the matrix. |
| `func (m *Matrix) Plus(b *Matrix) *Matrix` | Plus returns the sum a+b. |
| `func (m *Matrix) Row(i int) []float64` | Row returns a copy of row i as a slice. |
| `func (m *Matrix) Rows() int` | Rows returns the number of rows. |
| `func (m *Matrix) Scale(s float64) *Matrix` | Scale returns the matrix scaled by s. |
| `func (m *Matrix) Set(i, j int, v float64)` | Set assigns v to the element at row i, column j. |
| `func (m *Matrix) SetBlock(r0, c0 int, b *Matrix)` | SetBlock copies b into the block of m with top-left corner at (r0, c0). |
| `func (m *Matrix) SetCol(j int, v []float64)` | SetCol overwrites column j with v. |
| `func (m *Matrix) SetRow(i int, v []float64)` | SetRow overwrites row i with v. |
| `func (m *Matrix) String() string` | String renders the matrix for debugging. |
| `func (m *Matrix) Submatrix(r0, r1, c0, c1 int) *Matrix` | Submatrix returns the r0..r1-1 × c0..c1-1 block of the matrix. |
| `func (m *Matrix) Symmetrize() *Matrix` | Symmetrize returns (A + Aᵀ)/2, the symmetric part of a square matrix. |
| `func (m *Matrix) Trace() float64` | Trace returns the sum of the diagonal entries of a square matrix. |
| `func (m *Matrix) Transpose() *Matrix` | Transpose returns the transpose of the matrix. |
| `func (m *Matrix) Vec() []float64` | Vec returns the column-major vectorization of the matrix (columns stacked). |

</details>

<details>
<summary><code>TPBVPSolution</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func IndirectShooting( x0 []float64, f func(x, u []float64) []float64, g func(x, p, u []float64) []float64, uOpt func(x, p []float64) []float64, termGrad func(x []float64) []float64, tFinal float64, steps, maxIter int, tol float64, ) ([]float64, *TPBVPSolution, error)` | IndirectShooting solves a fixed-final-time two-point boundary-value problem by single shooting on the initial costate. |
| `func PontryaginLQ(a, b, q, r, qf *Matrix, x0 []float64, tFinal float64, steps int) (*TPBVPSolution, error)` | PontryaginLQ solves the fixed-final-time linear-quadratic optimal-control problem |

</details>

### Variables

`ErrDim`, `ErrNotConverged`, `ErrSingular`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/optimalcontrol).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
