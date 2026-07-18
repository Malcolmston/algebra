package numpde

import (
	"math"
	"testing"
)

const tinyTol = 1e-12

func maxDiffInterior(got []float64, exact func(x float64) float64, x0, dx float64) float64 {
	m := 0.0
	for i := 1; i < len(got)-1; i++ {
		d := math.Abs(got[i] - exact(x0+float64(i)*dx))
		if d > m {
			m = d
		}
	}
	return m
}

// ---------------------------------------------------------------------------
// Grid / helper utilities
// ---------------------------------------------------------------------------

func TestLinspace(t *testing.T) {
	xs := Linspace(0, 1, 5)
	want := []float64{0, 0.25, 0.5, 0.75, 1}
	for i := range want {
		if math.Abs(xs[i]-want[i]) > tinyTol {
			t.Fatalf("Linspace[%d]=%v want %v", i, xs[i], want[i])
		}
	}
	if xs[len(xs)-1] != 1 {
		t.Fatalf("Linspace endpoint not exact: %v", xs[len(xs)-1])
	}
}

func TestGrid1D(t *testing.T) {
	g := NewGrid1D(0, 2, 5)
	if math.Abs(g.Dx()-0.5) > tinyTol {
		t.Fatalf("Dx=%v want 0.5", g.Dx())
	}
	if g.Len() != 5 {
		t.Fatalf("Len=%d want 5", g.Len())
	}
	if math.Abs(g.X(2)-1) > tinyTol {
		t.Fatalf("X(2)=%v want 1", g.X(2))
	}
	u := g.Sample(func(x float64) float64 { return x * x })
	if math.Abs(u[4]-4) > tinyTol {
		t.Fatalf("Sample x^2 at x=2 = %v want 4", u[4])
	}
}

func TestGrid2D(t *testing.T) {
	g := NewGrid2D(0, 1, 0, 2, 3, 5)
	if math.Abs(g.Dx()-0.5) > tinyTol || math.Abs(g.Dy()-0.5) > tinyTol {
		t.Fatalf("Dx=%v Dy=%v", g.Dx(), g.Dy())
	}
	u := g.Sample(func(x, y float64) float64 { return x + y })
	if math.Abs(u[2][4]-3) > tinyTol {
		t.Fatalf("Sample x+y at (1,2)=%v want 3", u[2][4])
	}
}

func TestNorms(t *testing.T) {
	u := []float64{3, -4}
	if math.Abs(L2Norm(u)-5) > tinyTol {
		t.Fatalf("L2Norm=%v want 5", L2Norm(u))
	}
	if math.Abs(MaxAbs(u)-4) > tinyTol {
		t.Fatalf("MaxAbs=%v want 4", MaxAbs(u))
	}
	if math.Abs(LInfNorm(u)-4) > tinyTol {
		t.Fatalf("LInfNorm=%v want 4", LInfNorm(u))
	}
	if math.Abs(RMSNorm(u)-5/math.Sqrt2) > tinyTol {
		t.Fatalf("RMSNorm=%v", RMSNorm(u))
	}
}

func TestThomasSolve(t *testing.T) {
	// Solve the 3x3 system:
	//   2x0 -  x1        = 1
	//  -x0 + 2x1 - x2    = 0
	//         -x1 + 2x2  = 1
	// Exact solution x = (1, 1, 1).
	a := []float64{0, -1, -1}
	b := []float64{2, 2, 2}
	c := []float64{-1, -1, 0}
	d := []float64{1, 0, 1}
	x := ThomasSolve(a, b, c, d)
	for i := range x {
		if math.Abs(x[i]-1) > 1e-12 {
			t.Fatalf("ThomasSolve x[%d]=%v want 1", i, x[i])
		}
	}
}

func TestMatVec(t *testing.T) {
	m := [][]float64{{1, 2}, {3, 4}}
	y := MatVec(m, []float64{1, 1})
	if y[0] != 3 || y[1] != 7 {
		t.Fatalf("MatVec=%v want [3 7]", y)
	}
}

// ---------------------------------------------------------------------------
// Stencils
// ---------------------------------------------------------------------------

func TestStencilsOnPolynomials(t *testing.T) {
	dx := 0.1
	// u = x^2 sampled around x=1: values at 0.9,1.0,1.1
	u := []float64{0.81, 1.0, 1.21}
	// second derivative of x^2 is 2
	s2 := SecondDerivativeStencil(dx)
	if got := s2.Apply(u, 1); math.Abs(got-2) > 1e-9 {
		t.Fatalf("second derivative stencil=%v want 2", got)
	}
	// central first derivative of x^2 at x=1 is 2
	s1 := CentralFirstDerivativeStencil(dx)
	if got := s1.Apply(u, 1); math.Abs(got-2) > 1e-9 {
		t.Fatalf("central first derivative=%v want 2", got)
	}
}

func TestFivePointLaplacianStencil(t *testing.T) {
	c, e, w, n, s := FivePointLaplacianStencil(0.5, 0.5)
	// sum of weights must be zero (annihilates constants)
	if math.Abs(c+e+w+n+s) > 1e-12 {
		t.Fatalf("stencil weights do not sum to zero: %v", c+e+w+n+s)
	}
}

func TestLaplacian2DOnQuadratic(t *testing.T) {
	// u = x^2 - y^2 is harmonic; discrete 5-point Laplacian must be ~0.
	g := NewGrid2D(0, 1, 0, 1, 11, 11)
	u := g.Sample(func(x, y float64) float64 { return x*x - y*y })
	lap := Laplacian2D(u, g.Dx(), g.Dy())
	for i := 1; i < 10; i++ {
		for j := 1; j < 10; j++ {
			if math.Abs(lap[i][j]) > 1e-9 {
				t.Fatalf("Laplacian of harmonic at (%d,%d)=%v", i, j, lap[i][j])
			}
		}
	}
}

// ---------------------------------------------------------------------------
// 1D heat equation: exact u = sin(pi x) exp(-alpha pi^2 t) on [0,1]
// ---------------------------------------------------------------------------

func heat1DExact(x, t, alpha float64) float64 {
	return math.Sin(math.Pi*x) * math.Exp(-alpha*math.Pi*math.Pi*t)
}

func TestHeat1DExplicit(t *testing.T) {
	alpha := 1.0
	n := 41
	g := NewGrid1D(0, 1, n)
	dx := g.Dx()
	dt := 0.4 * dx * dx / alpha // mesh ratio 0.4 (stable)
	steps := 200
	tf := float64(steps) * dt
	u0 := g.Sample(func(x float64) float64 { return heat1DExact(x, 0, alpha) })
	u := Heat1DExplicit(u0, alpha, dx, dt, steps)
	err := maxDiffInterior(u, func(x float64) float64 { return heat1DExact(x, tf, alpha) }, 0, dx)
	if err > 2e-3 {
		t.Fatalf("Heat1DExplicit max error %v exceeds 2e-3", err)
	}
}

func TestHeat1DImplicit(t *testing.T) {
	alpha := 1.0
	n := 41
	g := NewGrid1D(0, 1, n)
	dx := g.Dx()
	dt := 2 * dx * dx / alpha // large step (unconditionally stable)
	steps := 20
	tf := float64(steps) * dt
	u0 := g.Sample(func(x float64) float64 { return heat1DExact(x, 0, alpha) })
	u := Heat1DImplicit(u0, alpha, dx, dt, steps)
	err := maxDiffInterior(u, func(x float64) float64 { return heat1DExact(x, tf, alpha) }, 0, dx)
	if err > 5e-3 {
		t.Fatalf("Heat1DImplicit max error %v exceeds 5e-3", err)
	}
}

func TestHeat1DCrankNicolson(t *testing.T) {
	alpha := 1.0
	n := 41
	g := NewGrid1D(0, 1, n)
	dx := g.Dx()
	dt := dx // large step; CN is second order and unconditionally stable
	steps := 50
	tf := float64(steps) * dt
	u0 := g.Sample(func(x float64) float64 { return heat1DExact(x, 0, alpha) })
	u := Heat1DCrankNicolson(u0, alpha, dx, dt, steps)
	err := maxDiffInterior(u, func(x float64) float64 { return heat1DExact(x, tf, alpha) }, 0, dx)
	if err > 5e-3 {
		t.Fatalf("Heat1DCrankNicolson max error %v exceeds 5e-3", err)
	}
}

func TestHeat1DExplicitUnstablePanics(t *testing.T) {
	defer func() {
		if recover() == nil {
			t.Fatal("expected panic for unstable mesh ratio")
		}
	}()
	u0 := make([]float64, 10)
	Heat1DExplicit(u0, 1, 0.1, 1, 1) // r = 100, unstable
}

// ---------------------------------------------------------------------------
// 2D heat equation: exact u = sin(pi x) sin(pi y) exp(-2 alpha pi^2 t)
// ---------------------------------------------------------------------------

func heat2DExact(x, y, t, alpha float64) float64 {
	return math.Sin(math.Pi*x) * math.Sin(math.Pi*y) * math.Exp(-2*alpha*math.Pi*math.Pi*t)
}

func maxDiff2DInterior(got [][]float64, g Grid2D, exact func(x, y float64) float64) float64 {
	m := 0.0
	for i := 1; i < g.Nx-1; i++ {
		for j := 1; j < g.Ny-1; j++ {
			d := math.Abs(got[i][j] - exact(g.X(i), g.Y(j)))
			if d > m {
				m = d
			}
		}
	}
	return m
}

func TestHeat2DExplicit(t *testing.T) {
	alpha := 1.0
	g := NewGrid2D(0, 1, 0, 1, 21, 21)
	dx, dy := g.Dx(), g.Dy()
	dt := 0.2 * Heat2DStableDt(alpha, dx, dy)
	steps := 100
	tf := float64(steps) * dt
	u0 := g.Sample(func(x, y float64) float64 { return heat2DExact(x, y, 0, alpha) })
	u := Heat2DExplicit(u0, alpha, dx, dy, dt, steps)
	err := maxDiff2DInterior(u, g, func(x, y float64) float64 { return heat2DExact(x, y, tf, alpha) })
	if err > 5e-3 {
		t.Fatalf("Heat2DExplicit max error %v exceeds 5e-3", err)
	}
}

func TestHeat2DADI(t *testing.T) {
	alpha := 1.0
	g := NewGrid2D(0, 1, 0, 1, 21, 21)
	dx, dy := g.Dx(), g.Dy()
	dt := 5 * Heat2DStableDt(alpha, dx, dy) // large step, ADI is stable
	steps := 40
	tf := float64(steps) * dt
	u0 := g.Sample(func(x, y float64) float64 { return heat2DExact(x, y, 0, alpha) })
	u := Heat2DADI(u0, alpha, dx, dy, dt, steps)
	err := maxDiff2DInterior(u, g, func(x, y float64) float64 { return heat2DExact(x, y, tf, alpha) })
	if err > 1e-2 {
		t.Fatalf("Heat2DADI max error %v exceeds 1e-2", err)
	}
}

// ---------------------------------------------------------------------------
// Wave equation: exact u = sin(pi x) cos(pi c t) on [0,1]
// ---------------------------------------------------------------------------

func TestWave1DExplicit(t *testing.T) {
	c := 1.0
	g := NewGrid1D(0, 1, 81)
	dx := g.Dx()
	dt := 0.5 * dx / c // Courant 0.5
	steps := 100
	tf := float64(steps) * dt
	u0 := g.Sample(func(x float64) float64 { return math.Sin(math.Pi * x) })
	v0 := make([]float64, g.N) // zero initial velocity
	u := Wave1DExplicit(u0, v0, c, dx, dt, steps)
	exact := func(x float64) float64 { return math.Sin(math.Pi*x) * math.Cos(math.Pi*c*tf) }
	err := maxDiffInterior(u, exact, 0, dx)
	if err > 5e-3 {
		t.Fatalf("Wave1DExplicit max error %v exceeds 5e-3", err)
	}
}

func TestWave2DExplicit(t *testing.T) {
	c := 1.0
	g := NewGrid2D(0, 1, 0, 1, 41, 41)
	dx, dy := g.Dx(), g.Dy()
	dt := 0.4 / (c * math.Sqrt(1/(dx*dx)+1/(dy*dy)))
	steps := 60
	tf := float64(steps) * dt
	freq := c * math.Pi * math.Sqrt(2)
	u0 := g.Sample(func(x, y float64) float64 { return math.Sin(math.Pi*x) * math.Sin(math.Pi*y) })
	v0 := Zeros2D(g.Nx, g.Ny)
	u := Wave2DExplicit(u0, v0, c, dx, dy, dt, steps)
	exact := func(x, y float64) float64 {
		return math.Sin(math.Pi*x) * math.Sin(math.Pi*y) * math.Cos(freq*tf)
	}
	err := maxDiff2DInterior(u, g, exact)
	if err > 2e-2 {
		t.Fatalf("Wave2DExplicit max error %v exceeds 2e-2", err)
	}
}

// ---------------------------------------------------------------------------
// Laplace / Poisson
// ---------------------------------------------------------------------------

func setupHarmonicBoundary(g Grid2D, f func(x, y float64) float64) [][]float64 {
	u := Zeros2D(g.Nx, g.Ny)
	for i := 0; i < g.Nx; i++ {
		for j := 0; j < g.Ny; j++ {
			if i == 0 || i == g.Nx-1 || j == 0 || j == g.Ny-1 {
				u[i][j] = f(g.X(i), g.Y(j))
			}
		}
	}
	return u
}

func TestLaplaceGaussSeidelHarmonic(t *testing.T) {
	g := NewGrid2D(0, 1, 0, 1, 21, 21)
	exact := func(x, y float64) float64 { return x*x - y*y }
	u0 := setupHarmonicBoundary(g, exact)
	res := LaplaceGaussSeidel(u0, g.Dx(), g.Dy(), 1e-10, 20000)
	if !res.Converged {
		t.Fatalf("Gauss-Seidel did not converge (iters=%d res=%v)", res.Iterations, res.Residual)
	}
	err := maxDiff2DInterior(res.Solution, g, exact)
	if err > 1e-4 {
		t.Fatalf("Laplace GS error %v exceeds 1e-4", err)
	}
}

func TestLaplaceSORHarmonic(t *testing.T) {
	g := NewGrid2D(0, 1, 0, 1, 31, 31)
	exact := func(x, y float64) float64 { return x*x - y*y }
	u0 := setupHarmonicBoundary(g, exact)
	omega := OptimalSOROmega(g.Nx, g.Ny)
	if omega <= 1 || omega >= 2 {
		t.Fatalf("OptimalSOROmega out of range: %v", omega)
	}
	res := LaplaceSOR(u0, g.Dx(), g.Dy(), omega, 1e-11, 20000)
	if !res.Converged {
		t.Fatalf("SOR did not converge (iters=%d)", res.Iterations)
	}
	err := maxDiff2DInterior(res.Solution, g, exact)
	if err > 1e-4 {
		t.Fatalf("Laplace SOR error %v exceeds 1e-4", err)
	}
}

func TestLaplaceJacobiHarmonic(t *testing.T) {
	g := NewGrid2D(0, 1, 0, 1, 17, 17)
	exact := func(x, y float64) float64 { return x*x - y*y }
	u0 := setupHarmonicBoundary(g, exact)
	res := LaplaceJacobi(u0, g.Dx(), g.Dy(), 1e-10, 50000)
	if !res.Converged {
		t.Fatalf("Jacobi did not converge (iters=%d)", res.Iterations)
	}
	err := maxDiff2DInterior(res.Solution, g, exact)
	if err > 1e-3 {
		t.Fatalf("Laplace Jacobi error %v exceeds 1e-3", err)
	}
}

func TestPoissonSOR(t *testing.T) {
	// -laplacian u = 2 pi^2 sin sin, i.e. laplacian u = -2 pi^2 sin sin.
	g := NewGrid2D(0, 1, 0, 1, 41, 41)
	exact := func(x, y float64) float64 { return math.Sin(math.Pi*x) * math.Sin(math.Pi*y) }
	f := g.Sample(func(x, y float64) float64 {
		return -2 * math.Pi * math.Pi * math.Sin(math.Pi*x) * math.Sin(math.Pi*y)
	})
	u0 := Zeros2D(g.Nx, g.Ny) // zero boundaries
	omega := OptimalSOROmega(g.Nx, g.Ny)
	res := PoissonSOR(u0, f, g.Dx(), g.Dy(), omega, 1e-10, 20000)
	if !res.Converged {
		t.Fatalf("Poisson SOR did not converge (iters=%d)", res.Iterations)
	}
	err := maxDiff2DInterior(res.Solution, g, exact)
	if err > 3e-3 {
		t.Fatalf("Poisson SOR error %v exceeds 3e-3", err)
	}
}

func TestPoisson1D(t *testing.T) {
	// u_xx = f with u = sin(pi x), f = -pi^2 sin(pi x), u(0)=u(1)=0.
	g := NewGrid1D(0, 1, 101)
	dx := g.Dx()
	f := g.Sample(func(x float64) float64 { return -math.Pi * math.Pi * math.Sin(math.Pi*x) })
	u0 := make([]float64, g.N) // zero boundaries
	u := Poisson1D(u0, f, dx)
	err := maxDiffInterior(u, func(x float64) float64 { return math.Sin(math.Pi * x) }, 0, dx)
	if err > 2e-3 {
		t.Fatalf("Poisson1D error %v exceeds 2e-3", err)
	}
}

// ---------------------------------------------------------------------------
// Advection
// ---------------------------------------------------------------------------

func TestAdvectionUpwindExactShift(t *testing.T) {
	// At Courant number 1 the upwind scheme is an exact one-cell shift.
	n := 40
	L := 1.0
	dx := L / float64(n)
	a := 1.0
	dt := dx // nu = 1
	u0 := make([]float64, n)
	for i := 0; i < n; i++ {
		u0[i] = math.Sin(2 * math.Pi * float64(i) * dx)
	}
	steps := 10
	u := Advection1DUpwind(u0, a, dx, dt, steps)
	// exact: shift by steps cells
	for i := 0; i < n; i++ {
		want := math.Sin(2 * math.Pi * float64(numpdeWrap(i-steps, n)) * dx)
		if math.Abs(u[i]-want) > 1e-10 {
			t.Fatalf("upwind nu=1 shift error at %d: got %v want %v", i, u[i], want)
		}
	}
}

func TestAdvectionLaxWendroffAccuracy(t *testing.T) {
	n := 128
	L := 1.0
	dx := L / float64(n)
	a := 1.0
	dt := 0.5 * dx // nu = 0.5
	steps := 128   // t = 0.5*128*dx = 0.5 -> shift 0.5
	tf := float64(steps) * dt
	u0 := make([]float64, n)
	for i := 0; i < n; i++ {
		u0[i] = math.Sin(2 * math.Pi * float64(i) * dx)
	}
	u := Advection1DLaxWendroff(u0, a, dx, dt, steps)
	maxErr := 0.0
	for i := 0; i < n; i++ {
		want := math.Sin(2 * math.Pi * (float64(i)*dx - a*tf))
		if d := math.Abs(u[i] - want); d > maxErr {
			maxErr = d
		}
	}
	if maxErr > 3e-2 {
		t.Fatalf("Lax-Wendroff error %v exceeds 3e-2", maxErr)
	}
}

func TestAdvectionLaxFriedrichsStable(t *testing.T) {
	// Lax-Friedrichs is stable and bounded; verify no blow-up and sane range.
	n := 64
	dx := 1.0 / float64(n)
	u0 := make([]float64, n)
	for i := 0; i < n; i++ {
		u0[i] = math.Sin(2 * math.Pi * float64(i) * dx)
	}
	u := Advection1DLaxFriedrichs(u0, 1, dx, 0.8*dx, 50)
	if MaxAbs(u) > 1.0001 {
		t.Fatalf("Lax-Friedrichs amplified beyond initial amplitude: %v", MaxAbs(u))
	}
}

// ---------------------------------------------------------------------------
// Method of lines
// ---------------------------------------------------------------------------

func TestMOLHeatRK4(t *testing.T) {
	alpha := 1.0
	g := NewGrid1D(0, 1, 41)
	dx := g.Dx()
	dt := 0.5 * dx * dx / alpha // within RK4 stability for this operator
	steps := 100
	tf := float64(steps) * dt
	u0 := g.Sample(func(x float64) float64 { return heat1DExact(x, 0, alpha) })
	f := MOLHeat1DField(alpha, dx)
	u := MOLIntegrate(f, u0, 0, dt, steps, RK4Step)
	err := maxDiffInterior(u, func(x float64) float64 { return heat1DExact(x, tf, alpha) }, 0, dx)
	if err > 2e-3 {
		t.Fatalf("MOL heat RK4 error %v exceeds 2e-3", err)
	}
}

func TestSteppersOnLinearODE(t *testing.T) {
	// y' = -y, y(0)=1, exact y(T)=exp(-T).
	f := func(_ float64, y []float64) []float64 { return []float64{-y[0]} }
	T := 1.0
	steps := 100
	dt := T / float64(steps)
	exact := math.Exp(-T)
	yr := MOLIntegrate(f, []float64{1}, 0, dt, steps, RK4Step)
	if math.Abs(yr[0]-exact) > 1e-8 {
		t.Fatalf("RK4 on y'=-y: %v want %v", yr[0], exact)
	}
	ye := MOLIntegrate(f, []float64{1}, 0, dt, steps, EulerStep)
	if math.Abs(ye[0]-exact) > 1e-2 {
		t.Fatalf("Euler on y'=-y too far: %v want ~%v", ye[0], exact)
	}
	y2 := MOLIntegrate(f, []float64{1}, 0, dt, steps, RK2Step)
	if math.Abs(y2[0]-exact) > 1e-4 {
		t.Fatalf("RK2 on y'=-y: %v want %v", y2[0], exact)
	}
}

// ---------------------------------------------------------------------------
// Benchmark for the heaviest routine
// ---------------------------------------------------------------------------

func BenchmarkPoissonSOR(b *testing.B) {
	g := NewGrid2D(0, 1, 0, 1, 65, 65)
	f := g.Sample(func(x, y float64) float64 {
		return -2 * math.Pi * math.Pi * math.Sin(math.Pi*x) * math.Sin(math.Pi*y)
	})
	omega := OptimalSOROmega(g.Nx, g.Ny)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		u0 := Zeros2D(g.Nx, g.Ny)
		_ = PoissonSOR(u0, f, g.Dx(), g.Dy(), omega, 1e-8, 5000)
	}
}
