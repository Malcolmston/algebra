# numpde — provides finite-difference solvers and supporting utilities for numerically

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/numpde.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/numpde)

Package numpde provides finite-difference solvers and supporting utilities for
numerically approximating solutions to partial differential equations (PDEs)
in one and two space dimensions.

- The heat (diffusion) equation u_t = alpha * u_xx, solved with the explicit
FTCS scheme, the fully implicit BTCS scheme, and the unconditionally stable
Crank-Nicolson scheme, plus a 2D explicit solver and an
alternating-direction-implicit (ADI) solver. - The wave equation u_tt = c^2 *
u_xx, solved with the standard explicit leap-frog (central-time central-space)
scheme in 1D and 2D. - The Laplace/Poisson equation, solved with the Jacobi,
Gauss-Seidel and successive-over-relaxation (SOR) iterative methods, plus a
direct tridiagonal solver for the 1D case. - The linear advection equation u_t
+ a * u_x = 0, solved with the first-order upwind, Lax-Friedrichs and
second-order Lax-Wendroff schemes, all using periodic boundary conditions.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/numpde"
```

## Usage

This is the package's own `ExampleAdvection1DCourant`, so it compiles and its output is
asserted on every `go test ./numpde/`.

```go
fmt.Printf("%.2f\n", numpde.Advection1DCourant(2, 0.1, 0.04))
	fmt.Printf("%.2f\n", numpde.Advection1DStableDt(2, 0.1))
```

```
0.80
0.05
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func Advection1DCourant(a, dx, dt float64) float64` | Advection1DCourant returns the (signed) Courant number nu = a*dt/dx. |
| `func Advection1DLaxFriedrichs(u0 []float64, a, dx, dt float64, steps int) []float64` | Advection1DLaxFriedrichs solves u_t + a*u_x = 0 on a periodic grid with the Lax-Friedrichs scheme, taking the given number of steps. |
| `func Advection1DLaxFriedrichsStep(u []float64, nu float64) []float64` | Advection1DLaxFriedrichsStep advances one step of the Lax-Friedrichs scheme on a periodic grid with Courant number nu. |
| `func Advection1DLaxWendroff(u0 []float64, a, dx, dt float64, steps int) []float64` | Advection1DLaxWendroff solves u_t + a*u_x = 0 on a periodic grid with the second-order Lax-Wendroff scheme, taking the given number of steps. |
| `func Advection1DLaxWendroffStep(u []float64, nu float64) []float64` | Advection1DLaxWendroffStep advances one step of the second-order accurate Lax-Wendroff scheme on a periodic grid with Courant number nu. |
| `func Advection1DStableDt(a, dx float64) float64` | Advection1DStableDt returns the largest stable time step dx/\|a\| for the explicit advection schemes (the CFL limit). |
| `func Advection1DUpwind(u0 []float64, a, dx, dt float64, steps int) []float64` | Advection1DUpwind solves u_t + a*u_x = 0 on a periodic grid of spacing dx with the first-order upwind scheme, taking the given number of steps of… |
| `func Advection1DUpwindStep(u []float64, nu float64) []float64` | Advection1DUpwindStep advances one step of the first-order upwind scheme on a periodic grid with Courant number nu = a*dt/dx. |
| `func Clone(u []float64) []float64` | Clone returns a copy of the slice u. |
| `func Clone2D(u [][]float64) [][]float64` | Clone2D returns a deep copy of the matrix u. |
| `func EulerStep(f ODEFunc, t float64, y []float64, dt float64) []float64` | EulerStep advances the ODE system y' = f(t, y) by one explicit (forward) Euler step of size dt and returns the new state. |
| `func FirstDerivativeMatrix(n int, dx float64) [][]float64` | FirstDerivativeMatrix returns the dense n by n central-difference approximation of d/dx with spacing dx. |
| `func FivePointLaplacianStencil(dx, dy float64) (center, east, west, north, south float64)` | FivePointLaplacianStencil returns the coefficients of the standard second-order five-point discrete Laplacian on a grid with spacings dx and dy. |
| `func Gradient1D(u []float64, dx float64) []float64` | Gradient1D returns the discrete first derivative of u on a grid of spacing dx. |
| `func Heat1DCrankNicolson(u0 []float64, alpha, dx, dt float64, steps int) []float64` | Heat1DCrankNicolson solves u_t = alpha*u_xx with the Crank-Nicolson scheme, taking the given number of steps of size dt on a grid of spacing dx. |
| `func Heat1DCrankNicolsonStep(u []float64, r float64) []float64` | Heat1DCrankNicolsonStep advances u by one step of the second-order accurate, unconditionally stable Crank-Nicolson scheme with mesh ratio r, holding… |
| `func Heat1DExplicit(u0 []float64, alpha, dx, dt float64, steps int) []float64` | Heat1DExplicit solves u_t = alpha*u_xx on the grid implied by u0 (spacing dx) by taking the given number of explicit FTCS steps of size dt. |
| `func Heat1DExplicitStep(u []float64, r float64) []float64` | Heat1DExplicitStep advances the field u by one time step of the explicit forward-time centred-space (FTCS) scheme with mesh ratio r = alpha*dt/dx^2,… |
| `func Heat1DImplicit(u0 []float64, alpha, dx, dt float64, steps int) []float64` | Heat1DImplicit solves u_t = alpha*u_xx with the unconditionally stable BTCS scheme, taking the given number of steps of size dt on a grid of spacing… |
| `func Heat1DImplicitStep(u []float64, r float64) []float64` | Heat1DImplicitStep advances u by one step of the fully implicit backward-time centred-space (BTCS) scheme with mesh ratio r, holding the boundary… |
| `func Heat1DMeshRatio(alpha, dx, dt float64) float64` | Heat1DMeshRatio returns the dimensionless mesh (Fourier) ratio r = alpha*dt/dx^2 that governs the stability and accuracy of the finite-difference… |
| `func Heat1DStableDt(alpha, dx float64) float64` | Heat1DStableDt returns the largest time step for which the explicit FTCS scheme applied to u_t = alpha*u_xx is stable, namely dt = dx^2/(2*alpha). |
| `func Heat2DADI(u0 [][]float64, alpha, dx, dy, dt float64, steps int) [][]float64` | Heat2DADI solves u_t = alpha*(u_xx + u_yy) with the unconditionally stable Peaceman-Rachford ADI scheme, taking the given number of steps of size dt… |
| `func Heat2DADIStep(u [][]float64, alpha, dx, dy, dt float64) [][]float64` | Heat2DADIStep advances u by one Peaceman-Rachford alternating-direction implicit (ADI) step for u_t = alpha*(u_xx + u_yy). |
| `func Heat2DExplicit(u0 [][]float64, alpha, dx, dy, dt float64, steps int) [][]float64` | Heat2DExplicit solves u_t = alpha*(u_xx + u_yy) on the grid implied by u0 (spacings dx, dy) with the explicit FTCS scheme, taking the given number of… |
| `func Heat2DExplicitStep(u [][]float64, rx, ry float64) [][]float64` | Heat2DExplicitStep advances the field u (stored as u`i``j`) by one explicit FTCS step for u_t = alpha*(u_xx + u_yy), holding every boundary node… |
| `func Heat2DMeshRatio(alpha, dx, dy, dt float64) float64` | Heat2DMeshRatio returns the combined explicit stability parameter alpha*dt*(1/dx^2 + 1/dy^2) for the 2D heat equation. |
| `func Heat2DStableDt(alpha, dx, dy float64) float64` | Heat2DStableDt returns the largest stable explicit time step for u_t = alpha*(u_xx + u_yy), namely dt = 1 / (2*alpha*(1/dx^2 + 1/dy^2)). |
| `func L2Norm(u []float64) float64` | L2Norm returns the Euclidean (2-norm) length of u, sqrt(sum u_i^2). |
| `func LInfNorm(u []float64) float64` | LInfNorm returns the maximum-norm (infinity norm) of u, the largest absolute element. |
| `func Laplacian1D(u []float64, dx float64) []float64` | Laplacian1D returns the discrete second derivative of the field u on a grid of spacing dx using the three-point central stencil at interior points. |
| `func Laplacian2D(u [][]float64, dx, dy float64) [][]float64` | Laplacian2D returns the discrete Laplacian of the field u (stored as u`i``j`) using the five-point stencil with spacings dx and dy. |
| `func Linspace(a, b float64, n int) []float64` | Linspace returns n evenly spaced samples over the closed interval [a, b]. |
| `func MOLAdvection1DRHS(u []float64, a, dx float64) []float64` | MOLAdvection1DRHS returns the semi-discrete right-hand side du/dt for the linear advection equation u_t + a*u_x = 0 using the second-order central… |
| `func MOLHeat1DRHS(u []float64, alpha, dx float64) []float64` | MOLHeat1DRHS returns the semi-discrete right-hand side du/dt for the heat equation u_t = alpha*u_xx obtained by replacing u_xx with the three-point… |
| `func MOLIntegrate(f ODEFunc, y0 []float64, t0, dt float64, steps int, stepper Stepper) []float64` | MOLIntegrate advances y0 through the given number of steps of size dt using the supplied stepper, starting at time t0, and returns the final state. |
| `func MatVec(m [][]float64, x []float64) []float64` | MatVec multiplies the dense matrix m by the vector x and returns m*x. |
| `func MaxAbs(u []float64) float64` | MaxAbs returns the largest absolute value among the elements of u, or zero when u is empty. |
| `func MaxAbsDiff(a, b []float64) float64` | MaxAbsDiff returns the infinity norm of the elementwise difference a-b. |
| `func MaxAbsDiff2D(a, b [][]float64) float64` | MaxAbsDiff2D returns the infinity norm of the elementwise difference of two equally shaped matrices. |
| `func NinePointLaplacianStencil(h float64) (center, edge, corner float64)` | NinePointLaplacianStencil returns the weights of the fourth-order accurate nine-point discrete Laplacian on a square grid of spacing h. |
| `func OptimalSOROmega(nx, ny int) float64` | OptimalSOROmega returns the theoretically optimal over-relaxation factor for the SOR solution of the five-point Laplacian on an nx by ny Dirichlet… |
| `func Poisson1D(u0, f []float64, dx float64) []float64` | Poisson1D directly solves the two-point boundary-value problem u_xx = f on a uniform grid of spacing dx using the tridiagonal (Thomas) algorithm. |
| `func RK2Step(f ODEFunc, t float64, y []float64, dt float64) []float64` | RK2Step advances the ODE system by one step of the explicit midpoint method (second-order Runge-Kutta) of size dt and returns the new state. |
| `func RK4Step(f ODEFunc, t float64, y []float64, dt float64) []float64` | RK4Step advances the ODE system by one step of the classical fourth-order Runge-Kutta method of size dt and returns the new state. |
| `func RMSNorm(u []float64) float64` | RMSNorm returns the root-mean-square norm of u, sqrt(mean(u_i^2)). |
| `func SecondDerivativeMatrix(n int, dx float64) [][]float64` | SecondDerivativeMatrix returns the dense n by n central-difference approximation of the second-derivative operator d^2/dx^2 with grid spacing dx and… |
| `func ThomasSolve(a, b, c, d []float64) []float64` | ThomasSolve solves the tridiagonal linear system M x = d, where M has sub-diagonal a, main diagonal b and super-diagonal c. |
| `func Wave1DCourant(c, dx, dt float64) float64` | Wave1DCourant returns the Courant number C = c*dt/dx for the 1D wave equation. |
| `func Wave1DExplicit(u0, v0 []float64, c, dx, dt float64, steps int) []float64` | Wave1DExplicit solves u_tt = c^2*u_xx with the explicit leap-frog (central time, central space) scheme. |
| `func Wave1DStableDt(c, dx float64) float64` | Wave1DStableDt returns the largest stable explicit time step for the 1D wave equation, dt = dx/c (the Courant-Friedrichs-Lewy limit). |
| `func Wave1DStep(prev, curr []float64, c2 float64) []float64` | Wave1DStep advances one leap-frog step of the 1D wave equation given the solution at the previous two time levels. |
| `func Wave2DCourant(c, dx, dy, dt float64) float64` | Wave2DCourant returns the 2D Courant number c*dt*sqrt(1/dx^2 + 1/dy^2). |
| `func Wave2DExplicit(u0, v0 [][]float64, c, dx, dy, dt float64, steps int) [][]float64` | Wave2DExplicit solves u_tt = c^2*(u_xx + u_yy) with the explicit leap-frog scheme. |
| `func Zeros2D(nx, ny int) [][]float64` | Zeros2D allocates an nx by ny matrix with every entry set to zero. |

### Types

| Type | What it is |
| --- | --- |
| `BCType` | BCType enumerates the kinds of boundary condition understood by the solvers. |
| `Field1D` | Field1D is a scalar function of a single spatial coordinate, x -> u(x). |
| `Field2D` | Field2D is a scalar function of two spatial coordinates, (x, y) -> u(x, y). |
| `Grid1D` | Grid1D describes a uniform one-dimensional mesh of N points covering the closed interval [X0, X1]. |
| `Grid2D` | Grid2D describes a uniform two-dimensional mesh of Nx by Ny points covering the rectangle [X0, X1] x [Y0, Y1]. |
| `IterResult` | IterResult reports the outcome of an iterative linear solver. |
| `ODEFunc` | ODEFunc is the right-hand side of an autonomous or non-autonomous system of ordinary differential equations, y' = f(t, y). |
| `Stencil` | Stencil holds the offset/weight pairs of a finite-difference approximation to a derivative on a uniform grid. |
| `Stepper` | Stepper advances an ODE system by one step of size dt from state y at time t using the right-hand side f, returning the new state. |

<details>
<summary><code>BCType</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (b BCType) String() string` | String returns a human-readable name for the boundary-condition type. |

</details>

<details>
<summary><code>Grid1D</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewGrid1D(x0, x1 float64, n int) Grid1D` | NewGrid1D constructs a Grid1D over [x0, x1] with n points. |
| `func (g Grid1D) Dx() float64` | Dx returns the uniform spacing between adjacent grid points. |
| `func (g Grid1D) Len() int` | Len returns the number of grid points. |
| `func (g Grid1D) Points() []float64` | Points returns a freshly allocated slice containing every grid coordinate. |
| `func (g Grid1D) Sample(f Field1D) []float64` | Sample evaluates f at every grid point and returns the resulting values. |
| `func (g Grid1D) X(i int) float64` | X returns the coordinate of the i-th grid point. |

</details>

<details>
<summary><code>Grid2D</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewGrid2D(x0, x1, y0, y1 float64, nx, ny int) Grid2D` | NewGrid2D constructs a Grid2D over [x0,x1] x [y0,y1] with nx by ny points. |
| `func (g Grid2D) Dx() float64` | Dx returns the uniform grid spacing in the x direction. |
| `func (g Grid2D) Dy() float64` | Dy returns the uniform grid spacing in the y direction. |
| `func (g Grid2D) Sample(f Field2D) [][]float64` | Sample evaluates f at every grid node and returns the values as u`i``j`. |
| `func (g Grid2D) X(i int) float64` | X returns the x coordinate of column index i. |
| `func (g Grid2D) Y(j int) float64` | Y returns the y coordinate of row index j. |

</details>

<details>
<summary><code>IterResult</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func LaplaceGaussSeidel(u0 [][]float64, dx, dy, tol float64, maxIter int) IterResult` | LaplaceGaussSeidel solves u_xx + u_yy = 0 with the Gauss-Seidel iteration, using the border of u0 as fixed Dirichlet data. |
| `func LaplaceJacobi(u0 [][]float64, dx, dy, tol float64, maxIter int) IterResult` | LaplaceJacobi solves the Laplace equation u_xx + u_yy = 0 with the Jacobi iteration, using the border of u0 as fixed Dirichlet data. |
| `func LaplaceSOR(u0 [][]float64, dx, dy, omega, tol float64, maxIter int) IterResult` | LaplaceSOR solves u_xx + u_yy = 0 with successive over-relaxation using factor omega, using the border of u0 as fixed Dirichlet data. |
| `func PoissonGaussSeidel(u0, f [][]float64, dx, dy, tol float64, maxIter int) IterResult` | PoissonGaussSeidel solves u_xx + u_yy = f with the Gauss-Seidel iteration, which updates each node in place using the most recent neighbour values… |
| `func PoissonJacobi(u0, f [][]float64, dx, dy, tol float64, maxIter int) IterResult` | PoissonJacobi solves u_xx + u_yy = f with the Jacobi iteration. |
| `func PoissonSOR(u0, f [][]float64, dx, dy, omega, tol float64, maxIter int) IterResult` | PoissonSOR solves u_xx + u_yy = f with successive over-relaxation using the relaxation factor omega (omega == 1 reduces to Gauss-Seidel; 1 < omega <… |

</details>

<details>
<summary><code>ODEFunc</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func MOLAdvection1DField(a, dx float64) ODEFunc` | MOLAdvection1DField returns an ODEFunc that evaluates MOLAdvection1DRHS with the given advection speed and grid spacing, suitable for MOLIntegrate. |
| `func MOLHeat1DField(alpha, dx float64) ODEFunc` | MOLHeat1DField returns an ODEFunc that evaluates MOLHeat1DRHS with the given diffusivity and grid spacing, suitable for passing to MOLIntegrate. |

</details>

<details>
<summary><code>Stencil</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func BackwardFirstDerivativeStencil(dx float64) Stencil` | BackwardFirstDerivativeStencil returns the first-order accurate backward difference stencil for d/dx: weights [-1,1]/dx at offsets [-1,0]. |
| `func CentralFirstDerivativeStencil(dx float64) Stencil` | CentralFirstDerivativeStencil returns the second-order accurate central difference stencil for d/dx on a grid of spacing dx: weights [-1,0,1]/(2dx)… |
| `func ForwardFirstDerivativeStencil(dx float64) Stencil` | ForwardFirstDerivativeStencil returns the first-order accurate forward difference stencil for d/dx: weights [-1,1]/dx at offsets [0,1]. |
| `func SecondDerivativeStencil(dx float64) Stencil` | SecondDerivativeStencil returns the second-order accurate central difference stencil for d^2/dx^2: weights [1,-2,1]/dx^2 at offsets [-1,0,1]. |
| `func (s Stencil) Apply(u []float64, i int) float64` | Apply evaluates the stencil at interior index i of the slice u. |

</details>

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/numpde).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
