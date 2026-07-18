package optimize

import (
	"math"
	"testing"
)

// optmzApprox reports whether a and b agree to within tol.
func optmzApprox(a, b, tol float64) bool { return math.Abs(a-b) <= tol }

// optmzVecApprox reports whether every component of a and b agrees to within
// tol.
func optmzVecApprox(a, b []float64, tol float64) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if math.Abs(a[i]-b[i]) > tol {
			return false
		}
	}
	return true
}

func TestVectorUtilities(t *testing.T) {
	a := []float64{1, 2, 3}
	b := []float64{4, 5, 6}
	if got := VecDot(a, b); !optmzApprox(got, 32, 1e-12) {
		t.Errorf("VecDot = %v, want 32", got)
	}
	if got := VecNorm([]float64{3, 4}); !optmzApprox(got, 5, 1e-12) {
		t.Errorf("VecNorm = %v, want 5", got)
	}
	if got := VecNormSquared([]float64{3, 4}); !optmzApprox(got, 25, 1e-12) {
		t.Errorf("VecNormSquared = %v, want 25", got)
	}
	if got := VecInfNorm([]float64{-7, 2, 3}); !optmzApprox(got, 7, 1e-12) {
		t.Errorf("VecInfNorm = %v, want 7", got)
	}
	if got := VecAdd(a, b); !optmzVecApprox(got, []float64{5, 7, 9}, 1e-12) {
		t.Errorf("VecAdd = %v", got)
	}
	if got := VecSub(b, a); !optmzVecApprox(got, []float64{3, 3, 3}, 1e-12) {
		t.Errorf("VecSub = %v", got)
	}
	if got := VecScale(a, 2); !optmzVecApprox(got, []float64{2, 4, 6}, 1e-12) {
		t.Errorf("VecScale = %v", got)
	}
	if got := VecAxpy(2, a, b); !optmzVecApprox(got, []float64{6, 9, 12}, 1e-12) {
		t.Errorf("VecAxpy = %v", got)
	}
	if got := VecLinComb(2, a, 3, b); !optmzVecApprox(got, []float64{14, 19, 24}, 1e-12) {
		t.Errorf("VecLinComb = %v", got)
	}
	if got := VecDistance([]float64{0, 0}, []float64{3, 4}); !optmzApprox(got, 5, 1e-12) {
		t.Errorf("VecDistance = %v, want 5", got)
	}
	if got := VecClamp([]float64{-2, 0.5, 9}, 0, 1); !optmzVecApprox(got, []float64{0, 0.5, 1}, 1e-12) {
		t.Errorf("VecClamp = %v", got)
	}
	if got := ProjectBox([]float64{5, -5}, []float64{0, 0}, []float64{1, 1}); !optmzVecApprox(got, []float64{1, 0}, 1e-12) {
		t.Errorf("ProjectBox = %v", got)
	}
	if got := Centroid([][]float64{{0, 0}, {2, 4}}); !optmzVecApprox(got, []float64{1, 2}, 1e-12) {
		t.Errorf("Centroid = %v", got)
	}
}

func TestMatrixUtilities(t *testing.T) {
	id := MatIdentity(3)
	x := []float64{1, 2, 3}
	if got := MatVec(id, x); !optmzVecApprox(got, x, 1e-12) {
		t.Errorf("MatVec(I, x) = %v, want %v", got, x)
	}
	a := [][]float64{{1, 2}, {3, 4}}
	if got := MatTranspose(a); got[0][1] != 3 || got[1][0] != 2 {
		t.Errorf("MatTranspose = %v", got)
	}
	op := OuterProduct([]float64{1, 2}, []float64{3, 4})
	if op[0][0] != 3 || op[1][1] != 8 {
		t.Errorf("OuterProduct = %v", op)
	}
}

func TestSolveLinearSystem(t *testing.T) {
	// 2x + y = 5 ; x - y = 1  => x = 2, y = 1.
	a := [][]float64{{2, 1}, {1, -1}}
	b := []float64{5, 1}
	sol, err := SolveLinearSystem(a, b)
	if err != nil {
		t.Fatalf("SolveLinearSystem error: %v", err)
	}
	if !optmzVecApprox(sol, []float64{2, 1}, 1e-10) {
		t.Errorf("SolveLinearSystem = %v, want [2 1]", sol)
	}
	// Singular matrix must be reported.
	if _, err := SolveLinearSystem([][]float64{{1, 1}, {2, 2}}, []float64{1, 2}); err != ErrSingularMatrix {
		t.Errorf("expected ErrSingularMatrix, got %v", err)
	}
}

func TestNumericDerivatives(t *testing.T) {
	x := []float64{1.5, -2, 3}
	// Gradient of Sphere is 2x.
	want := SphereGrad(x)
	if got := NumericGradientCentral(Sphere, x, 1e-6); !optmzVecApprox(got, want, 1e-4) {
		t.Errorf("NumericGradientCentral = %v, want %v", got, want)
	}
	if got := NumericGradientForward(Sphere, x, 1e-6); !optmzVecApprox(got, want, 1e-3) {
		t.Errorf("NumericGradientForward = %v, want %v", got, want)
	}
	if got := PartialDerivative(Sphere, x, 0, 1e-6); !optmzApprox(got, 3, 1e-4) {
		t.Errorf("PartialDerivative = %v, want 3", got)
	}
	// Hessian of Sphere is 2I.
	h := NumericHessian(Sphere, x, 1e-4)
	for i := range h {
		for j := range h[i] {
			want := 0.0
			if i == j {
				want = 2
			}
			if !optmzApprox(h[i][j], want, 1e-3) {
				t.Errorf("NumericHessian[%d][%d] = %v, want %v", i, j, h[i][j], want)
			}
		}
	}
	// Jacobian of g(x) = [x0^2, x0*x1] at (2,3) is [[4,0],[3,2]].
	g := func(v []float64) []float64 { return []float64{v[0] * v[0], v[0] * v[1]} }
	jac := NumericJacobian(g, []float64{2, 3}, 1e-6)
	if !optmzApprox(jac[0][0], 4, 1e-3) || !optmzApprox(jac[1][0], 3, 1e-3) || !optmzApprox(jac[1][1], 2, 1e-3) {
		t.Errorf("NumericJacobian = %v", jac)
	}
	// Directional derivative of Sphere at (1,1) along (1,0) is 2.
	if got := DirectionalDerivative(Sphere, []float64{1, 1}, []float64{1, 0}, 1e-6); !optmzApprox(got, 2, 1e-4) {
		t.Errorf("DirectionalDerivative = %v, want 2", got)
	}
}

func TestScalarMinimizers(t *testing.T) {
	// (x-2)^2 has its minimum at x = 2.
	quad := func(x float64) float64 { return (x - 2) * (x - 2) }
	gs := GoldenSection(quad, -5, 10, 1e-10, 200)
	if !gs.Converged || !optmzApprox(gs.X, 2, 1e-4) {
		t.Errorf("GoldenSection = %+v, want X~2", gs)
	}
	bp := BrentParabolic(quad, -5, 10, 1e-10, 200)
	if !bp.Converged || !optmzApprox(bp.X, 2, 1e-6) {
		t.Errorf("BrentParabolic = %+v, want X~2", bp)
	}
	// cos has its minimum at x = pi on [0, 2pi], value -1.
	c := BrentParabolic(math.Cos, 0, 2*math.Pi, 1e-10, 200)
	if !optmzApprox(c.X, math.Pi, 1e-5) || !optmzApprox(c.F, -1, 1e-8) {
		t.Errorf("BrentParabolic(cos) = %+v, want X~pi F~-1", c)
	}
	ms := MinimizeScalar(quad, -5, 10)
	if !optmzApprox(ms.X, 2, 1e-5) {
		t.Errorf("MinimizeScalar = %+v, want X~2", ms)
	}
}

func TestLineSearchArmijo(t *testing.T) {
	x := []float64{2, 2}
	g := SphereGrad(x)  // (4, 4)
	dir := VecNegate(g) // steepest descent
	a := ArmijoStep(Sphere, x, dir, g)
	if a <= 0 || a > 1 {
		t.Fatalf("ArmijoStep returned %v", a)
	}
	// The accepted step must actually reduce the objective.
	if Sphere(VecAxpy(a, dir, x)) >= Sphere(x) {
		t.Errorf("Armijo step did not decrease the objective")
	}
}

func TestFirstOrderMinimizers(t *testing.T) {
	x0 := []float64{3, -4, 1.5}
	want := []float64{0, 0, 0}

	gm := GradientDescentMomentum(Sphere, SphereGrad, x0, 0.1, 0.9, 1e-9, 5000)
	if !optmzVecApprox(gm.X, want, 1e-3) {
		t.Errorf("GradientDescentMomentum = %v, want origin", gm.X)
	}
	gn := GradientDescentNesterov(Sphere, SphereGrad, x0, 0.1, 0.9, 1e-9, 5000)
	if !optmzVecApprox(gn.X, want, 1e-3) {
		t.Errorf("GradientDescentNesterov = %v, want origin", gn.X)
	}
	// Momentum descent with a numeric gradient (grad == nil) must also work.
	gnil := GradientDescentMomentum(Sphere, nil, x0, 0.1, 0.9, 1e-9, 5000)
	if !optmzVecApprox(gnil.X, want, 1e-3) {
		t.Errorf("GradientDescentMomentum(numeric) = %v, want origin", gnil.X)
	}
}

func TestCoordinateDescent(t *testing.T) {
	// Separable quadratic with minimum at (1, -2, 3).
	target := []float64{1, -2, 3}
	f := func(x []float64) float64 {
		var s float64
		for i := range x {
			d := x[i] - target[i]
			s += d * d
		}
		return s
	}
	res := CoordinateDescent(f, []float64{0, 0, 0}, 5, 1e-9, 500)
	if !optmzVecApprox(res.X, target, 1e-5) {
		t.Errorf("CoordinateDescent = %v, want %v", res.X, target)
	}
}

func TestConjugateGradient(t *testing.T) {
	// Quadratic bowl (x-a)^T(x-a) with a = (2, -1).
	a := []float64{2, -1}
	f := func(x []float64) float64 { return (x[0]-a[0])*(x[0]-a[0]) + (x[1]-a[1])*(x[1]-a[1]) }
	res := ConjugateGradient(f, nil, []float64{0, 0}, 1e-10, 200)
	if !optmzVecApprox(res.X, a, 1e-4) {
		t.Errorf("ConjugateGradient = %v, want %v", res.X, a)
	}
}

func TestNewtonMultivariate(t *testing.T) {
	// Booth is a convex quadratic with minimum (1, 3); Newton solves it fast.
	res := NewtonMultivariate(Booth, nil, nil, []float64{0, 0}, 1e-10, 100)
	if !optmzVecApprox(res.X, []float64{1, 3}, 1e-5) {
		t.Errorf("NewtonMultivariate(Booth) = %v, want [1 3]", res.X)
	}
	if !optmzApprox(res.F, 0, 1e-8) {
		t.Errorf("NewtonMultivariate F = %v, want 0", res.F)
	}
}

func TestBFGS(t *testing.T) {
	// BFGS on Rosenbrock from the classic hard start.
	res := BFGS(Rosenbrock, RosenbrockGrad, []float64{-1.2, 1}, 1e-8, 1000)
	if !optmzVecApprox(res.X, []float64{1, 1}, 1e-3) {
		t.Errorf("BFGS(Rosenbrock) = %v, want [1 1]", res.X)
	}
	if !optmzApprox(res.F, 0, 1e-6) {
		t.Errorf("BFGS F = %v, want 0", res.F)
	}
	// With a numeric gradient too.
	rn := BFGS(Sphere, nil, []float64{5, -3, 2}, 1e-9, 500)
	if !optmzVecApprox(rn.X, []float64{0, 0, 0}, 1e-4) {
		t.Errorf("BFGS(Sphere, numeric) = %v, want origin", rn.X)
	}
}

func TestNelderMead(t *testing.T) {
	// Sphere from an offset start.
	rs := NelderMead(Sphere, []float64{2, -3, 1}, 0.1, 1e-10, 2000)
	if !optmzVecApprox(rs.X, []float64{0, 0, 0}, 1e-4) {
		t.Errorf("NelderMead(Sphere) = %v, want origin", rs.X)
	}
	// Booth: minimum at (1, 3).
	rb := NelderMead(Booth, []float64{0, 0}, 0.1, 1e-12, 2000)
	if !optmzVecApprox(rb.X, []float64{1, 3}, 1e-4) {
		t.Errorf("NelderMead(Booth) = %v, want [1 3]", rb.X)
	}
	// Rosenbrock: minimum at (1, 1).
	rr := NelderMead(Rosenbrock, []float64{-1.2, 1}, 0.1, 1e-12, 5000)
	if !optmzVecApprox(rr.X, []float64{1, 1}, 1e-3) {
		t.Errorf("NelderMead(Rosenbrock) = %v, want [1 1]", rr.X)
	}
}

func TestSimulatedAnnealing(t *testing.T) {
	opts := DefaultSAOptions()
	opts.StepSize = 0.5
	opts.MaxIter = 20000
	r1 := SimulatedAnnealing(Sphere, []float64{5, -5}, opts, 42)
	if VecNorm(r1.X) > 0.2 {
		t.Errorf("SimulatedAnnealing did not approach origin: %v", r1.X)
	}
	// Determinism: identical seed => identical result.
	r2 := SimulatedAnnealing(Sphere, []float64{5, -5}, opts, 42)
	if !optmzVecApprox(r1.X, r2.X, 0) || r1.F != r2.F {
		t.Errorf("SimulatedAnnealing not deterministic: %v vs %v", r1, r2)
	}
}

func TestTestObjectives(t *testing.T) {
	if !optmzApprox(Sphere([]float64{0, 0, 0}), 0, 0) {
		t.Errorf("Sphere origin != 0")
	}
	if !optmzApprox(Rosenbrock([]float64{1, 1, 1}), 0, 1e-12) {
		t.Errorf("Rosenbrock(ones) != 0")
	}
	if g := RosenbrockGrad([]float64{1, 1, 1}); VecNorm(g) > 1e-12 {
		t.Errorf("RosenbrockGrad(ones) = %v, want 0", g)
	}
	if !optmzApprox(Booth([]float64{1, 3}), 0, 1e-12) {
		t.Errorf("Booth(1,3) != 0")
	}
	if !optmzApprox(Himmelblau([]float64{3, 2}), 0, 1e-12) {
		t.Errorf("Himmelblau(3,2) != 0")
	}
}

// BenchmarkNelderMead exercises the heaviest routine on Rosenbrock.
func BenchmarkNelderMead(b *testing.B) {
	for i := 0; i < b.N; i++ {
		NelderMead(Rosenbrock, []float64{-1.2, 1, -1, 1}, 0.1, 1e-10, 2000)
	}
}
