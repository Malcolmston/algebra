// Package numpde provides finite-difference solvers and supporting utilities
// for numerically approximating solutions to partial differential equations
// (PDEs) in one and two space dimensions.
//
// The package is organised around a handful of classic model problems:
//
//   - The heat (diffusion) equation u_t = alpha * u_xx, solved with the
//     explicit FTCS scheme, the fully implicit BTCS scheme, and the
//     unconditionally stable Crank-Nicolson scheme, plus a 2D explicit
//     solver and an alternating-direction-implicit (ADI) solver.
//   - The wave equation u_tt = c^2 * u_xx, solved with the standard explicit
//     leap-frog (central-time central-space) scheme in 1D and 2D.
//   - The Laplace/Poisson equation, solved with the Jacobi, Gauss-Seidel and
//     successive-over-relaxation (SOR) iterative methods, plus a direct
//     tridiagonal solver for the 1D case.
//   - The linear advection equation u_t + a * u_x = 0, solved with the
//     first-order upwind, Lax-Friedrichs and second-order Lax-Wendroff
//     schemes, all using periodic boundary conditions.
//
// In addition the package exposes method-of-lines (MOL) helpers that convert a
// spatial discretisation into a system of ordinary differential equations that
// can be advanced with the supplied Euler, RK2 and RK4 time steppers, as well
// as stencil builders and small dense/tridiagonal linear-algebra utilities.
//
// All routines depend only on the Go standard library and are deterministic:
// identical inputs always produce identical outputs. Grids are uniform, and
// unless stated otherwise Dirichlet boundary conditions are imposed by holding
// the first and last grid values fixed.
package numpde

import "math"

// Field1D is a scalar function of a single spatial coordinate, x -> u(x).
// It is commonly used to sample an initial condition onto a grid.
type Field1D func(x float64) float64

// Field2D is a scalar function of two spatial coordinates, (x, y) -> u(x, y).
type Field2D func(x, y float64) float64

// BCType enumerates the kinds of boundary condition understood by the solvers.
type BCType int

const (
	// Dirichlet fixes the value of the solution on the boundary.
	Dirichlet BCType = iota
	// Neumann fixes the outward normal derivative of the solution on the
	// boundary (a zero-flux/insulated boundary uses derivative zero).
	Neumann
)

// String returns a human-readable name for the boundary-condition type.
func (b BCType) String() string {
	switch b {
	case Dirichlet:
		return "Dirichlet"
	case Neumann:
		return "Neumann"
	default:
		return "Unknown"
	}
}

// Grid1D describes a uniform one-dimensional mesh of N points covering the
// closed interval [X0, X1]. The points are x_i = X0 + i*Dx for i in [0, N-1]
// with spacing Dx = (X1-X0)/(N-1).
type Grid1D struct {
	X0 float64 // left endpoint of the interval
	X1 float64 // right endpoint of the interval
	N  int     // number of grid points (must be at least 2)
}

// NewGrid1D constructs a Grid1D over [x0, x1] with n points. It panics if
// n < 2 or if x1 <= x0.
func NewGrid1D(x0, x1 float64, n int) Grid1D {
	if n < 2 {
		panic("numpde: NewGrid1D requires n >= 2")
	}
	if !(x1 > x0) {
		panic("numpde: NewGrid1D requires x1 > x0")
	}
	return Grid1D{X0: x0, X1: x1, N: n}
}

// Dx returns the uniform spacing between adjacent grid points.
func (g Grid1D) Dx() float64 { return (g.X1 - g.X0) / float64(g.N-1) }

// Len returns the number of grid points.
func (g Grid1D) Len() int { return g.N }

// X returns the coordinate of the i-th grid point.
func (g Grid1D) X(i int) float64 { return g.X0 + float64(i)*g.Dx() }

// Points returns a freshly allocated slice containing every grid coordinate.
func (g Grid1D) Points() []float64 {
	xs := make([]float64, g.N)
	dx := g.Dx()
	for i := range xs {
		xs[i] = g.X0 + float64(i)*dx
	}
	return xs
}

// Sample evaluates f at every grid point and returns the resulting values.
func (g Grid1D) Sample(f Field1D) []float64 {
	u := make([]float64, g.N)
	dx := g.Dx()
	for i := range u {
		u[i] = f(g.X0 + float64(i)*dx)
	}
	return u
}

// Grid2D describes a uniform two-dimensional mesh of Nx by Ny points covering
// the rectangle [X0, X1] x [Y0, Y1]. Field values are stored in row-major
// order as u[i][j] with i indexing x (0..Nx-1) and j indexing y (0..Ny-1).
type Grid2D struct {
	X0 float64 // left x endpoint
	X1 float64 // right x endpoint
	Y0 float64 // bottom y endpoint
	Y1 float64 // top y endpoint
	Nx int     // number of grid points in x (must be at least 2)
	Ny int     // number of grid points in y (must be at least 2)
}

// NewGrid2D constructs a Grid2D over [x0,x1] x [y0,y1] with nx by ny points.
// It panics if nx < 2, ny < 2, x1 <= x0 or y1 <= y0.
func NewGrid2D(x0, x1, y0, y1 float64, nx, ny int) Grid2D {
	if nx < 2 || ny < 2 {
		panic("numpde: NewGrid2D requires nx >= 2 and ny >= 2")
	}
	if !(x1 > x0) || !(y1 > y0) {
		panic("numpde: NewGrid2D requires x1 > x0 and y1 > y0")
	}
	return Grid2D{X0: x0, X1: x1, Y0: y0, Y1: y1, Nx: nx, Ny: ny}
}

// Dx returns the uniform grid spacing in the x direction.
func (g Grid2D) Dx() float64 { return (g.X1 - g.X0) / float64(g.Nx-1) }

// Dy returns the uniform grid spacing in the y direction.
func (g Grid2D) Dy() float64 { return (g.Y1 - g.Y0) / float64(g.Ny-1) }

// X returns the x coordinate of column index i.
func (g Grid2D) X(i int) float64 { return g.X0 + float64(i)*g.Dx() }

// Y returns the y coordinate of row index j.
func (g Grid2D) Y(j int) float64 { return g.Y0 + float64(j)*g.Dy() }

// Sample evaluates f at every grid node and returns the values as u[i][j].
func (g Grid2D) Sample(f Field2D) [][]float64 {
	u := Zeros2D(g.Nx, g.Ny)
	dx, dy := g.Dx(), g.Dy()
	for i := 0; i < g.Nx; i++ {
		x := g.X0 + float64(i)*dx
		for j := 0; j < g.Ny; j++ {
			u[i][j] = f(x, g.Y0+float64(j)*dy)
		}
	}
	return u
}

// IterResult reports the outcome of an iterative linear solver. Solution holds
// the final field, Iterations the number of sweeps performed, Residual the
// final infinity-norm residual (or update norm), and Converged whether the
// requested tolerance was met before the iteration cap.
type IterResult struct {
	Solution   [][]float64
	Iterations int
	Residual   float64
	Converged  bool
}

// numpdeAbs returns the absolute value of x. It exists so the hot loops avoid
// the function-call overhead comparison used in the standard library on some
// platforms while keeping intent explicit.
func numpdeAbs(x float64) float64 {
	if x < 0 {
		return -x
	}
	return x
}

// numpdeSqrt is a thin wrapper over math.Sqrt kept for symmetry with the other
// helpers used by the norm routines.
func numpdeSqrt(x float64) float64 { return math.Sqrt(x) }
