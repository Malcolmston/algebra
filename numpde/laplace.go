package numpde

import "math"

// The elliptic solvers in this file approximate the Poisson equation written in
// the form
//
//	u_xx + u_yy = f
//
// on a rectangular grid with Dirichlet boundary conditions. The boundary values
// are taken from the border of the supplied initial guess u0 and are held fixed
// throughout the iteration. The Laplace equation is the special case f == 0 and
// has dedicated convenience wrappers.

// OptimalSOROmega returns the theoretically optimal over-relaxation factor for
// the SOR solution of the five-point Laplacian on an nx by ny Dirichlet grid
// with equal spacing in both directions. It uses omega = 2/(1+sqrt(1-rho^2))
// where rho is the spectral radius of the corresponding Jacobi iteration,
// rho = (cos(pi/(nx-1)) + cos(pi/(ny-1)))/2. The result lies in (1, 2).
func OptimalSOROmega(nx, ny int) float64 {
	rho := 0.5 * (math.Cos(math.Pi/float64(nx-1)) + math.Cos(math.Pi/float64(ny-1)))
	return 2.0 / (1.0 + math.Sqrt(1.0-rho*rho))
}

// numpdePoissonResidual returns the infinity norm of the residual
// f - (u_xx + u_yy) evaluated with the five-point stencil over the interior of
// the grid.
func numpdePoissonResidual(u, f [][]float64, dx, dy float64) float64 {
	nx := len(u)
	ny := len(u[0])
	ix := 1.0 / (dx * dx)
	iy := 1.0 / (dy * dy)
	maxr := 0.0
	for i := 1; i < nx-1; i++ {
		for j := 1; j < ny-1; j++ {
			lap := (u[i-1][j]-2*u[i][j]+u[i+1][j])*ix +
				(u[i][j-1]-2*u[i][j]+u[i][j+1])*iy
			fij := 0.0
			if f != nil {
				fij = f[i][j]
			}
			if r := numpdeAbs(fij - lap); r > maxr {
				maxr = r
			}
		}
	}
	return maxr
}

// PoissonJacobi solves u_xx + u_yy = f with the Jacobi iteration. u0 supplies
// the initial guess and the fixed Dirichlet boundary values; f is the source
// term (pass nil for the Laplace equation). Iteration stops when the infinity
// norm of the update between sweeps falls below tol or after maxIter sweeps. The
// input matrices are not modified.
func PoissonJacobi(u0, f [][]float64, dx, dy, tol float64, maxIter int) IterResult {
	nx := len(u0)
	ny := len(u0[0])
	ix := 1.0 / (dx * dx)
	iy := 1.0 / (dy * dy)
	denom := 2*ix + 2*iy
	u := Clone2D(u0)
	next := Clone2D(u0)
	res := IterResult{Solution: u}
	for it := 1; it <= maxIter; it++ {
		diff := 0.0
		for i := 1; i < nx-1; i++ {
			for j := 1; j < ny-1; j++ {
				fij := 0.0
				if f != nil {
					fij = f[i][j]
				}
				val := ((u[i-1][j]+u[i+1][j])*ix + (u[i][j-1]+u[i][j+1])*iy - fij) / denom
				if d := numpdeAbs(val - u[i][j]); d > diff {
					diff = d
				}
				next[i][j] = val
			}
		}
		u, next = next, u
		res.Iterations = it
		if diff < tol {
			res.Converged = true
			res.Solution = u
			res.Residual = numpdePoissonResidual(u, f, dx, dy)
			return res
		}
	}
	res.Solution = u
	res.Residual = numpdePoissonResidual(u, f, dx, dy)
	return res
}

// PoissonGaussSeidel solves u_xx + u_yy = f with the Gauss-Seidel iteration,
// which updates each node in place using the most recent neighbour values and
// typically converges about twice as fast as Jacobi. Arguments and termination
// mirror PoissonJacobi. The input matrices are not modified.
func PoissonGaussSeidel(u0, f [][]float64, dx, dy, tol float64, maxIter int) IterResult {
	return numpdeSOR(u0, f, dx, dy, 1.0, tol, maxIter)
}

// PoissonSOR solves u_xx + u_yy = f with successive over-relaxation using the
// relaxation factor omega (omega == 1 reduces to Gauss-Seidel; 1 < omega < 2
// accelerates convergence). Arguments and termination mirror PoissonJacobi. The
// input matrices are not modified.
func PoissonSOR(u0, f [][]float64, dx, dy, omega, tol float64, maxIter int) IterResult {
	return numpdeSOR(u0, f, dx, dy, omega, tol, maxIter)
}

// numpdeSOR is the shared Gauss-Seidel/SOR kernel.
func numpdeSOR(u0, f [][]float64, dx, dy, omega, tol float64, maxIter int) IterResult {
	nx := len(u0)
	ny := len(u0[0])
	ix := 1.0 / (dx * dx)
	iy := 1.0 / (dy * dy)
	denom := 2*ix + 2*iy
	u := Clone2D(u0)
	res := IterResult{Solution: u}
	for it := 1; it <= maxIter; it++ {
		diff := 0.0
		for i := 1; i < nx-1; i++ {
			for j := 1; j < ny-1; j++ {
				fij := 0.0
				if f != nil {
					fij = f[i][j]
				}
				gs := ((u[i-1][j]+u[i+1][j])*ix + (u[i][j-1]+u[i][j+1])*iy - fij) / denom
				val := (1-omega)*u[i][j] + omega*gs
				if d := numpdeAbs(val - u[i][j]); d > diff {
					diff = d
				}
				u[i][j] = val
			}
		}
		res.Iterations = it
		if diff < tol {
			res.Converged = true
			res.Residual = numpdePoissonResidual(u, f, dx, dy)
			return res
		}
	}
	res.Residual = numpdePoissonResidual(u, f, dx, dy)
	return res
}

// LaplaceJacobi solves the Laplace equation u_xx + u_yy = 0 with the Jacobi
// iteration, using the border of u0 as fixed Dirichlet data. It is a
// convenience wrapper around PoissonJacobi with a zero source term.
func LaplaceJacobi(u0 [][]float64, dx, dy, tol float64, maxIter int) IterResult {
	return PoissonJacobi(u0, nil, dx, dy, tol, maxIter)
}

// LaplaceGaussSeidel solves u_xx + u_yy = 0 with the Gauss-Seidel iteration,
// using the border of u0 as fixed Dirichlet data.
func LaplaceGaussSeidel(u0 [][]float64, dx, dy, tol float64, maxIter int) IterResult {
	return PoissonGaussSeidel(u0, nil, dx, dy, tol, maxIter)
}

// LaplaceSOR solves u_xx + u_yy = 0 with successive over-relaxation using factor
// omega, using the border of u0 as fixed Dirichlet data.
func LaplaceSOR(u0 [][]float64, dx, dy, omega, tol float64, maxIter int) IterResult {
	return PoissonSOR(u0, nil, dx, dy, omega, tol, maxIter)
}

// Poisson1D directly solves the two-point boundary-value problem u_xx = f on a
// uniform grid of spacing dx using the tridiagonal (Thomas) algorithm. The
// boundary values u0[0] and u0[last] are imposed as Dirichlet data; f holds the
// source term at every node (interior entries are used). The exact discrete
// solution is returned in a new slice. It panics if len(u0) != len(f) or the
// grid has fewer than three points.
func Poisson1D(u0, f []float64, dx float64) []float64 {
	n := len(u0)
	if len(f) != n {
		panic("numpde: Poisson1D u0 and f length mismatch")
	}
	if n < 3 {
		panic("numpde: Poisson1D requires at least three grid points")
	}
	a := make([]float64, n)
	b := make([]float64, n)
	c := make([]float64, n)
	d := make([]float64, n)
	h2 := dx * dx
	b[0], d[0] = 1, u0[0]
	b[n-1], d[n-1] = 1, u0[n-1]
	for i := 1; i < n-1; i++ {
		a[i] = 1
		b[i] = -2
		c[i] = 1
		d[i] = h2 * f[i]
	}
	return ThomasSolve(a, b, c, d)
}
