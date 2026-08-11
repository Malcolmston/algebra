# linprog — implements linear and convex optimization in pure Go using only the standard library

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/linprog.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/linprog)

Package linprog implements linear and convex optimization in pure Go using
only the standard library.

- `LP` holds a general linear program (minimize or maximize a linear objective
subject to linear <=, =, >= constraints with x >= 0). It can be reduced to
`StandardLP` equality/nonnegative form and solved with the two-phase primal
`Simplex` method, which uses Bland's rule to guarantee termination without
cycling.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/linprog"
```

## Usage

This is the package's own `ExampleDualCanonical`, so it compiles and its output is
asserted on every `go test ./linprog/`.

```go
c := []float64{2, 3}
	a := [][]float64{{1, 1}, {1, 3}}
	b := []float64{4, 6}
	primal := linprog.NewLP(linprog.Minimize, c, a,
		[]linprog.Relation{linprog.GreaterEqual, linprog.GreaterEqual}, b)
	dual := linprog.DualCanonical(c, a, b)
	ps := primal.Solve()
	ds := dual.Solve()
	fmt.Printf("primal %.4f dual %.4f gap %.6f\n",
		ps.Objective, ds.Objective, linprog.DualityGap(ps.Objective, ds.Objective))
```

```
primal 9.0000 dual 9.0000 gap 0.000000
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func AssignmentCost(cost [][]float64, assign []int) float64` | AssignmentCost returns the total cost of the assignment, summing cost`i`[assign`i`] over all rows. |
| `func CloneMatrix(a [][]float64) [][]float64` | CloneMatrix returns an independent deep copy of the row-major matrix A. |
| `func ComplementarySlackness(c []float64, a [][]float64, b []float64, x, y []float64, tol float64) bool` | ComplementarySlackness reports whether a primal point x and dual point y satisfy the complementary slackness conditions for the canonical pair from… |
| `func Dot(a, b []float64) float64` | Dot returns the inner product of two equal-length vectors. |
| `func DualityGap(primalObj, dualObj float64) float64` | DualityGap returns primalObj - dualObj. |
| `func Hungarian(cost [][]float64) (assign []int, total float64, err error)` | Hungarian solves the assignment problem: given an n-by-n cost matrix, it finds a permutation assignment of rows to columns that minimizes the total… |
| `func HungarianMax(value [][]float64) (assign []int, total float64, err error)` | HungarianMax solves the maximization assignment problem by negating the cost matrix, running `Hungarian`, and returning the assignment together with… |
| `func InfNorm(v []float64) float64` | InfNorm returns the infinity norm (maximum absolute value) of a vector. |
| `func IsBalanced(supply, demand []float64, tol float64) bool` | IsBalanced reports whether total supply equals total demand within tol. |
| `func IsIntegral(x []float64, tol float64) bool` | IsIntegral reports whether every entry of x is within tol of an integer. |
| `func MatTVec(a [][]float64, x []float64) []float64` | MatTVec returns the product A^T*x, i.e. |
| `func MatVec(a [][]float64, x []float64) []float64` | MatVec returns the matrix-vector product A*x where A is an m-by-n row-major matrix and x has length n. |
| `func NorthWestCorner(supply, demand []float64) [][]float64` | NorthWestCorner returns an initial basic feasible flow for a balanced transportation problem using the north-west-corner rule. |
| `func SolveLinear(a [][]float64, b []float64) ([]float64, error)` | SolveLinear solves the square linear system A*x = b for x using Gaussian elimination with partial pivoting. |
| `func SolveQPEquality(q [][]float64, c []float64, a [][]float64, b []float64) (x, lambda []float64, err error)` | SolveQPEquality solves the equality-constrained convex QP |
| `func Transportation(supply, demand []float64, cost [][]float64) (flow [][]float64, total float64, err error)` | Transportation solves the balanced transportation problem: ship from m sources with the given supplies to n sinks with the given demands at the given… |
| `func Transpose(a [][]float64) [][]float64` | Transpose returns the transpose of the m-by-n row-major matrix A as an n-by-m matrix. |
| `func WeakDualityHolds(primalObj, dualObj, tol float64) bool` | WeakDualityHolds reports whether the primal (minimizing) and dual (maximizing) objective values respect weak duality, primalObj >= dualObj, within… |

### Types

| Type | What it is |
| --- | --- |
| `KKTResidual` | KKTResidual holds the first-order (Karush-Kuhn-Tucker) optimality residuals of a candidate solution to a `QP`. |
| `LP` | LP is a general linear program over nonnegative variables: |
| `LPKKT` | LPKKT holds the four first-order optimality residuals for a canonical primal/dual LP pair. |
| `LPSolution` | LPSolution is the outcome of solving a general `LP`. |
| `QP` | QP is a convex quadratic program: |
| `QPSolution` | QPSolution is the outcome of solving a `QP` with `SolveQP`. |
| `Relation` | Relation is the comparison operator of a linear constraint row. |
| `Sense` | Sense selects whether a linear or quadratic objective is minimized or maximized. |
| `SimplexResult` | SimplexResult is the outcome of running the `Simplex` method on a `StandardLP`. |
| `StandardLP` | StandardLP is a linear program in standard equality form: |
| `Status` | Status classifies the outcome of solving an optimization problem. |

<details>
<summary><code>KKTResidual</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func ComputeKKTResidual(qp QP, x, lambda, mu []float64) KKTResidual` | ComputeKKTResidual evaluates the `KKTResidual` of the candidate primal point x with equality multipliers lambda and inequality multipliers mu for the… |
| `func (k KKTResidual) Max() float64` | Max returns the largest of the five residuals. |
| `func (k KKTResidual) Satisfied(tol float64) bool` | Satisfied reports whether every residual is within tol. |

</details>

<details>
<summary><code>LP</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func DualCanonical(c []float64, a [][]float64, b []float64) LP` | DualCanonical builds the dual of the canonical primal |
| `func NewLP(sense Sense, c []float64, a [][]float64, rel []Relation, b []float64) LP` | NewLP constructs an `LP` from its parts. |
| `func (lp LP) AddConstraint(row []float64, rel Relation, rhs float64) LP` | AddConstraint returns a copy of the program with one additional constraint row appended. |
| `func (lp LP) Clone() LP` | Clone returns an independent deep copy of the program. |
| `func (lp LP) Feasible(x []float64, tol float64) bool` | Feasible reports whether x satisfies every constraint (within tol) and is nonnegative (within tol). |
| `func (lp LP) NumConstraints() int` | NumConstraints returns the number of constraint rows. |
| `func (lp LP) NumVars() int` | NumVars returns the number of decision variables. |
| `func (lp LP) Objective(x []float64) float64` | Objective returns the raw objective value C · x, ignoring `Sense`. |
| `func (lp LP) Solve() LPSolution` | Solve solves the general `LP` by reducing it to standard form and running the two-phase `Simplex`. |
| `func (lp LP) SolveBinary(binary []bool) LPSolution` | SolveBinary solves the integer program in which the variables flagged in binary are additionally restricted to {0,1} by appending x`j` <= 1 upper… |
| `func (lp LP) SolveInteger(integer []bool) LPSolution` | SolveInteger solves the mixed-integer linear program obtained from lp by requiring the variables flagged in integer to take integer values. |
| `func (lp LP) Standard() StandardLP` | Standard reduces the general program to `StandardLP` equality form. |

</details>

<details>
<summary><code>LPKKT</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func CheckKKTLP(c []float64, a [][]float64, b []float64, x, y []float64) LPKKT` | CheckKKTLP computes the `LPKKT` residuals for the canonical primal |
| `func (k LPKKT) Max() float64` | Max returns the largest of the four residuals. |
| `func (k LPKKT) Satisfied(tol float64) bool` | Satisfied reports whether every residual is within tol. |

</details>

<details>
<summary><code>QP</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (qp QP) Feasible(x []float64, tol float64) bool` | Feasible reports whether x satisfies A x = B (within tol) and G x <= H (within tol). |
| `func (qp QP) Gradient(x []float64) []float64` | Gradient returns the objective gradient Q x + C at x. |
| `func (qp QP) IsPositiveDefinite(tol float64) bool` | IsPositiveDefinite reports whether the Hessian Q is symmetric positive definite, tested by attempting a Cholesky factorization with the diagonal… |
| `func (qp QP) IsSymmetric(tol float64) bool` | IsSymmetric reports whether the Hessian Q is square and equal to its own transpose to within tol. |
| `func (qp QP) NumVars() int` | NumVars returns the number of decision variables. |
| `func (qp QP) Objective(x []float64) float64` | Objective returns 1/2 x^T Q x + C · x. |

</details>

<details>
<summary><code>QPSolution</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func SolveQP(qp QP, x0 []float64) QPSolution` | SolveQP solves the convex `QP` with inequality (and optional equality) constraints by a primal active-set method starting from the feasible point x0. |

</details>

<details>
<summary><code>Relation</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (r Relation) String() string` | String returns "<=", "=" or ">=". |

</details>

<details>
<summary><code>Sense</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (s Sense) String() string` | String returns "Minimize" or "Maximize". |

</details>

<details>
<summary><code>SimplexResult</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Simplex(std StandardLP) SimplexResult` | Simplex solves the standard-form program by the two-phase primal simplex method. |

</details>

<details>
<summary><code>StandardLP</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (s StandardLP) Feasible(x []float64, tol float64) bool` | Feasible reports whether x satisfies A x = B (within tol) and x >= 0. |
| `func (s StandardLP) NumConstraints() int` | NumConstraints returns the number of equality rows. |
| `func (s StandardLP) NumVars() int` | NumVars returns the number of variables in the standard program. |
| `func (s StandardLP) Objective(x []float64) float64` | Objective returns C · x. |

</details>

<details>
<summary><code>Status</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (s Status) String() string` | String returns a human-readable name for the status. |

</details>

### Constants

`Eps`, `IntTol`, `MaxBranchNodes`

### Variables

`ErrDimension`, `ErrNotSquare`, `ErrSingular`, `ErrUnbalanced`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/linprog).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
