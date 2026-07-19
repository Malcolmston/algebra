package optimalcontrol

import (
	"fmt"
	"math"
	"math/cmplx"
	"testing"
)

func approx(a, b, tol float64) bool { return math.Abs(a-b) <= tol }

func approxMat(a, b *Matrix, tol float64) bool { return a.ApproxEqual(b, tol) }

// ---- Linear algebra ----------------------------------------------------------

func TestSolveInverseDet(t *testing.T) {
	a := FromRows([][]float64{{2, 1}, {1, 3}})
	x, err := Solve(a, []float64{3, 5})
	if err != nil {
		t.Fatal(err)
	}
	if !approx(x[0], 0.8, 1e-12) || !approx(x[1], 1.4, 1e-12) {
		t.Errorf("Solve = %v want [0.8 1.4]", x)
	}
	d, _ := Det(a)
	if !approx(d, 5, 1e-12) {
		t.Errorf("Det = %v want 5", d)
	}
	inv, err := Inverse(a)
	if err != nil {
		t.Fatal(err)
	}
	prod := a.Mul(inv)
	if !approxMat(prod, Eye(2), 1e-12) {
		t.Errorf("A A^-1 != I:\n%v", prod)
	}
}

func TestCholesky(t *testing.T) {
	a := FromRows([][]float64{{4, 2}, {2, 3}})
	l, err := Cholesky(a)
	if err != nil {
		t.Fatal(err)
	}
	if !approxMat(l.Mul(l.Transpose()), a, 1e-12) {
		t.Errorf("L Lᵀ != A")
	}
	if !IsPositiveDefinite(a) {
		t.Errorf("expected SPD")
	}
	if IsPositiveDefinite(FromRows([][]float64{{1, 2}, {2, 1}})) {
		t.Errorf("indefinite matrix reported SPD")
	}
}

func TestMatrixExp(t *testing.T) {
	// exp([[0,1],[0,0]]) = [[1,1],[0,1]].
	n := FromRows([][]float64{{0, 1}, {0, 0}})
	if !approxMat(MatrixExp(n), FromRows([][]float64{{1, 1}, {0, 1}}), 1e-12) {
		t.Errorf("nilpotent exp wrong")
	}
	// exp of a rotation generator by pi/2.
	th := math.Pi / 2
	rot := FromRows([][]float64{{0, -th}, {th, 0}})
	got := MatrixExp(rot)
	want := FromRows([][]float64{{0, -1}, {1, 0}})
	if !approxMat(got, want, 1e-9) {
		t.Errorf("rotation exp = %v", got)
	}
}

func TestEigenvalues(t *testing.T) {
	a := FromRows([][]float64{{2, 0}, {0, 3}})
	ev := Eigenvalues(a)
	got := []float64{real(ev[0]), real(ev[1])}
	sortFloats(got)
	if !approx(got[0], 2, 1e-9) || !approx(got[1], 3, 1e-9) {
		t.Errorf("eigs = %v want [2 3]", got)
	}
	// Purely imaginary pair.
	rot := FromRows([][]float64{{0, -1}, {1, 0}})
	evr := Eigenvalues(rot)
	for _, z := range evr {
		if !approx(real(z), 0, 1e-9) || !approx(math.Abs(imag(z)), 1, 1e-9) {
			t.Errorf("rotation eig = %v want ±i", z)
		}
	}
	if !approx(SpectralRadius(rot), 1, 1e-9) {
		t.Errorf("spectral radius wrong")
	}
	if !approx(SpectralAbscissa(a), 3, 1e-9) {
		t.Errorf("spectral abscissa wrong")
	}
}

func TestCharPoly(t *testing.T) {
	// A = [[0,1],[-2,-3]] has char poly λ² + 3λ + 2.
	a := FromRows([][]float64{{0, 1}, {-2, -3}})
	c := CharPoly(a)
	want := []float64{2, 3, 1}
	for i := range want {
		if !approx(c[i], want[i], 1e-12) {
			t.Errorf("charpoly[%d]=%v want %v", i, c[i], want[i])
		}
	}
}

func TestSymEigen(t *testing.T) {
	a := FromRows([][]float64{{2, 1}, {1, 2}})
	w, err := SymEigenvalues(a)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(w[0], 1, 1e-9) || !approx(w[1], 3, 1e-9) {
		t.Errorf("sym eigs = %v want [1 3]", w)
	}
	vals, vecs, _ := SymEigen(a)
	// Check A V = V diag(vals).
	av := a.Mul(vecs)
	vd := vecs.Mul(Diag(vals))
	if !approxMat(av, vd, 1e-8) {
		t.Errorf("eigenvectors wrong")
	}
}

func TestStability(t *testing.T) {
	stable := FromRows([][]float64{{-1, 0}, {0, -2}})
	if !IsStableContinuous(stable, 1e-9) {
		t.Errorf("expected Hurwitz")
	}
	if IsStableContinuous(FromRows([][]float64{{1, 0}, {0, -2}}), 1e-9) {
		t.Errorf("unstable reported stable")
	}
	d := FromRows([][]float64{{0.5, 0}, {0, -0.4}})
	if !IsStableDiscrete(d, 1e-9) {
		t.Errorf("expected Schur stable")
	}
	if IsStableDiscrete(FromRows([][]float64{{1.2, 0}, {0, 0.1}}), 1e-9) {
		t.Errorf("unstable discrete reported stable")
	}
}

// ---- Lyapunov / Gramians -----------------------------------------------------

func TestLyapunovContinuous(t *testing.T) {
	a := FromRows([][]float64{{-1, 0}, {0, -2}})
	q := Eye(2)
	x, err := SolveLyapunovContinuous(a, q)
	if err != nil {
		t.Fatal(err)
	}
	want := Diag([]float64{0.5, 0.25})
	if !approxMat(x, want, 1e-10) {
		t.Errorf("lyap X = %v want diag(0.5,0.25)", x)
	}
	res := LyapunovContinuousResidual(a, q, x)
	if res.MaxAbs() > 1e-10 {
		t.Errorf("residual too big: %v", res.MaxAbs())
	}
}

func TestLyapunovDiscrete(t *testing.T) {
	a := FromRows([][]float64{{0.5}})
	q := FromRows([][]float64{{1}})
	x, err := SolveLyapunovDiscrete(a, q)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(x.At(0, 0), 4.0/3.0, 1e-12) {
		t.Errorf("discrete lyap = %v want 1.3333", x.At(0, 0))
	}
}

func TestSylvester(t *testing.T) {
	a := FromRows([][]float64{{2, 1}, {0, 3}})
	b := FromRows([][]float64{{-5, 0}, {1, -6}})
	c := FromRows([][]float64{{1, 1}, {2, 0}})
	x, err := SolveSylvester(a, b, c)
	if err != nil {
		t.Fatal(err)
	}
	res := a.Mul(x).Plus(x.Mul(b)).Minus(c)
	if res.MaxAbs() > 1e-10 {
		t.Errorf("sylvester residual %v", res.MaxAbs())
	}
}

// ---- Structure ---------------------------------------------------------------

func TestControllabilityObservability(t *testing.T) {
	a := FromRows([][]float64{{0, 1}, {0, 0}})
	b := FromRows([][]float64{{0}, {1}})
	if !IsControllable(a, b) {
		t.Errorf("double integrator should be controllable")
	}
	if ControllabilityRank(a, b) != 2 {
		t.Errorf("ctrb rank wrong")
	}
	c := FromRows([][]float64{{1, 0}})
	if !IsObservable(a, c) {
		t.Errorf("should be observable")
	}
	// Uncontrollable pair.
	au := Eye(2)
	bu := FromRows([][]float64{{1}, {0}})
	if IsControllable(au, bu) {
		t.Errorf("expected uncontrollable")
	}
	if !IsStabilizableContinuous(a, b) {
		t.Errorf("double integrator stabilizable")
	}
	if !IsDetectableContinuous(a, c) {
		t.Errorf("expected detectable")
	}
}

// ---- Riccati / LQR -----------------------------------------------------------

func TestCAREDoubleIntegrator(t *testing.T) {
	a := FromRows([][]float64{{0, 1}, {0, 0}})
	b := FromRows([][]float64{{0}, {1}})
	q := Eye(2)
	r := Eye(1)
	p, err := SolveCARE(a, b, q, r)
	if err != nil {
		t.Fatal(err)
	}
	s3 := math.Sqrt(3)
	want := FromRows([][]float64{{s3, 1}, {1, s3}})
	if !approxMat(p, want, 1e-6) {
		t.Errorf("CARE P = %v want [[√3,1],[1,√3]]", p)
	}
	res, _ := CAREResidual(a, b, q, r, p)
	if res.MaxAbs() > 1e-6 {
		t.Errorf("CARE residual %v", res.MaxAbs())
	}
	lqr, err := LQRContinuous(a, b, q, r)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(lqr.K.At(0, 0), 1, 1e-6) || !approx(lqr.K.At(0, 1), s3, 1e-6) {
		t.Errorf("LQR gain = %v want [1, √3]", lqr.K)
	}
	if !IsStableContinuous(lqr.ClosedLoop, 1e-6) {
		t.Errorf("closed loop unstable")
	}
}

func TestCAREScalar(t *testing.T) {
	// a=1 (unstable), b=1, q=1, r=1 -> x^2 - 2x - 1 = 0 -> x = 1+√2.
	a := FromRows([][]float64{{1}})
	b := FromRows([][]float64{{1}})
	q := FromRows([][]float64{{1}})
	r := FromRows([][]float64{{1}})
	p, err := SolveCARE(a, b, q, r)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(p.At(0, 0), 1+math.Sqrt2, 1e-8) {
		t.Errorf("scalar CARE = %v want %v", p.At(0, 0), 1+math.Sqrt2)
	}
}

func TestCAREKleinman(t *testing.T) {
	a := FromRows([][]float64{{0, 1}, {0, 0}})
	b := FromRows([][]float64{{0}, {1}})
	q := Eye(2)
	r := Eye(1)
	// Stabilizing initial gain placing poles at -1,-1: A-BK = [[0,1],[-k1,-k2]].
	k0 := FromRows([][]float64{{1, 2}})
	p, err := SolveCAREKleinman(a, b, q, r, k0, 50, 1e-12)
	if err != nil && err != ErrNotConverged {
		t.Fatal(err)
	}
	s3 := math.Sqrt(3)
	want := FromRows([][]float64{{s3, 1}, {1, s3}})
	if !approxMat(p, want, 1e-6) {
		t.Errorf("Kleinman P = %v", p)
	}
}

func TestDAREScalar(t *testing.T) {
	// a=b=q=r=1 -> x^2 - x - 1 = 0 -> golden ratio.
	a := FromRows([][]float64{{1}})
	b := FromRows([][]float64{{1}})
	q := FromRows([][]float64{{1}})
	r := FromRows([][]float64{{1}})
	p, err := SolveDARE(a, b, q, r)
	if err != nil {
		t.Fatal(err)
	}
	phi := (1 + math.Sqrt(5)) / 2
	if !approx(p.At(0, 0), phi, 1e-8) {
		t.Errorf("DARE = %v want golden ratio %v", p.At(0, 0), phi)
	}
	lqr, err := LQRDiscrete(a, b, q, r)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(lqr.K.At(0, 0), 1/phi, 1e-6) {
		t.Errorf("DARE gain = %v want %v", lqr.K.At(0, 0), 1/phi)
	}
	res, _ := DAREResidual(a, b, q, r, p)
	if res.MaxAbs() > 1e-7 {
		t.Errorf("DARE residual %v", res.MaxAbs())
	}
}

func TestDAREMatrix(t *testing.T) {
	a := FromRows([][]float64{{1, 0.1}, {0, 1}})
	b := FromRows([][]float64{{0}, {0.1}})
	q := Eye(2)
	r := FromRows([][]float64{{1}})
	p, err := SolveDARE(a, b, q, r)
	if err != nil {
		t.Fatal(err)
	}
	res, _ := DAREResidual(a, b, q, r, p)
	if res.MaxAbs() > 1e-6 {
		t.Errorf("DARE matrix residual %v", res.MaxAbs())
	}
	lqr, _ := LQRDiscrete(a, b, q, r)
	if !IsStableDiscrete(lqr.ClosedLoop, 1e-6) {
		t.Errorf("discrete closed loop unstable")
	}
}

// ---- Finite horizon ----------------------------------------------------------

func TestFiniteHorizonDiscreteConverges(t *testing.T) {
	a := FromRows([][]float64{{1, 0.1}, {0, 1}})
	b := FromRows([][]float64{{0}, {0.1}})
	q := Eye(2)
	r := FromRows([][]float64{{1}})
	// Long horizon should approach the infinite-horizon DARE solution.
	fh, err := FiniteHorizonLQRDiscrete(a, b, q, r, Eye(2), 300)
	if err != nil {
		t.Fatal(err)
	}
	pInf, _ := SolveDARE(a, b, q, r)
	if !approxMat(fh.P[0], pInf, 1e-4) {
		t.Errorf("finite horizon did not approach DARE:\n%v\n%v", fh.P[0], pInf)
	}
	// Simulated trajectory should drive the state toward the origin.
	traj := SimulateDiscreteLQR(a, b, fh, []float64{1, 0})
	if VecNorm(traj[len(traj)-1]) >= VecNorm(traj[0]) {
		t.Errorf("state did not decrease")
	}
}

func TestFiniteHorizonContinuousConverges(t *testing.T) {
	a := FromRows([][]float64{{0, 1}, {0, 0}})
	b := FromRows([][]float64{{0}, {1}})
	q := Eye(2)
	r := Eye(1)
	fh, err := FiniteHorizonLQRContinuous(a, b, q, r, Eye(2), 5, 2000)
	if err != nil {
		t.Fatal(err)
	}
	pInf, _ := SolveCARE(a, b, q, r)
	if !approxMat(fh.P[0], pInf, 1e-3) {
		t.Errorf("continuous finite horizon P(0) = %v want %v", fh.P[0], pInf)
	}
}

// ---- MDP ---------------------------------------------------------------------

func buildTestMDP() *MDP {
	// Action 0: stay; Action 1: jump to state 0.
	stay := FromRows([][]float64{{1, 0}, {0, 1}})
	jump := FromRows([][]float64{{1, 0}, {1, 0}})
	reward := FromRows([][]float64{{1, 1}, {0, -0.1}})
	m, _ := NewMDP([]*Matrix{stay, jump}, reward, 0.5)
	return m
}

func TestMDPValueIteration(t *testing.T) {
	m := buildTestMDP()
	res := m.ValueIteration(1e-12, 1000)
	if !res.Converged {
		t.Errorf("value iteration did not converge")
	}
	if !approx(res.Value[0], 2, 1e-6) || !approx(res.Value[1], 0.9, 1e-6) {
		t.Errorf("VI value = %v want [2 0.9]", res.Value)
	}
	if res.Policy[1] != 1 {
		t.Errorf("policy[1] = %d want 1", res.Policy[1])
	}
}

func TestMDPPolicyIteration(t *testing.T) {
	m := buildTestMDP()
	res, err := m.PolicyIteration(100)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(res.Value[0], 2, 1e-9) || !approx(res.Value[1], 0.9, 1e-9) {
		t.Errorf("PI value = %v want [2 0.9]", res.Value)
	}
	if res.Policy[1] != 1 {
		t.Errorf("PI policy[1] = %d want 1", res.Policy[1])
	}
	// Consistency across the three DP methods.
	gs := m.GaussSeidelValueIteration(1e-12, 1000)
	mpi, _ := m.ModifiedPolicyIteration(20, 100, 1e-12)
	for s := 0; s < m.States; s++ {
		if !approx(gs.Value[s], res.Value[s], 1e-5) {
			t.Errorf("Gauss-Seidel mismatch at %d", s)
		}
		if !approx(mpi.Value[s], res.Value[s], 1e-5) {
			t.Errorf("MPI mismatch at %d", s)
		}
	}
}

func TestMDPPolicyEvaluation(t *testing.T) {
	m := buildTestMDP()
	policy := []int{0, 1}
	exact, err := m.PolicyEvaluationExact(policy)
	if err != nil {
		t.Fatal(err)
	}
	iter := m.PolicyEvaluationIterative(policy, 1e-12, 10000)
	for s := range exact {
		if !approx(exact[s], iter[s], 1e-5) {
			t.Errorf("policy eval mismatch at %d: %v vs %v", s, exact[s], iter[s])
		}
	}
	if !approx(exact[0], 2, 1e-9) {
		t.Errorf("exact eval V(0) = %v want 2", exact[0])
	}
}

// ---- HJB ---------------------------------------------------------------------

func TestHJBGrid1D(t *testing.T) {
	grid := make([]float64, 41)
	for i := range grid {
		grid[i] = -2 + 0.1*float64(i)
	}
	h := &HJBGrid1D{
		Grid:        grid,
		Controls:    []float64{-1, 0, 1},
		Dynamics:    func(x, u float64) float64 { return u },
		RunningCost: func(x, u float64) float64 { return x*x + 0.01*u*u },
		Dt:          0.1,
		Rho:         1.0,
	}
	res := h.Solve(1e-10, 5000)
	if !res.Converged {
		t.Errorf("HJB did not converge")
	}
	// Value at the origin (index 20) should be near zero and smallest.
	zeroIdx := 20
	if math.Abs(grid[zeroIdx]) > 1e-9 {
		t.Fatalf("grid[20] should be 0, got %v", grid[zeroIdx])
	}
	for i := range grid {
		if res.Value[i] < res.Value[zeroIdx]-1e-9 {
			t.Errorf("value at %v below value at origin", grid[i])
		}
	}
	// Control should push toward the origin.
	if h.ControlAt(res, 1.5) > 0 {
		t.Errorf("control at x=1.5 should be <= 0")
	}
	if h.ControlAt(res, -1.5) < 0 {
		t.Errorf("control at x=-1.5 should be >= 0")
	}
}

func TestLinearInterp1D(t *testing.T) {
	g := []float64{0, 1, 2}
	v := []float64{0, 10, 20}
	if !approx(LinearInterp1D(g, v, 0.5), 5, 1e-12) {
		t.Errorf("interp wrong")
	}
	if !approx(LinearInterp1D(g, v, -1), 0, 1e-12) {
		t.Errorf("clamp low wrong")
	}
	if !approx(LinearInterp1D(g, v, 5), 20, 1e-12) {
		t.Errorf("clamp high wrong")
	}
}

// ---- Pontryagin --------------------------------------------------------------

func TestPontryaginLQ(t *testing.T) {
	a := FromRows([][]float64{{0, 1}, {0, 0}})
	b := FromRows([][]float64{{0}, {1}})
	q := Eye(2)
	r := Eye(1)
	qf := Eye(2)
	x0 := []float64{1, 0.5}
	sol, err := PontryaginLQ(a, b, q, r, qf, x0, 2.0, 100)
	if err != nil {
		t.Fatal(err)
	}
	// In LQ, p(0) = P(0) x0 where P is the Riccati solution.
	fh, _ := FiniteHorizonLQRContinuous(a, b, q, r, qf, 2.0, 4000)
	pStart := fh.P[0].MulVec(x0)
	for i := range pStart {
		if !approx(sol.P[0][i], pStart[i], 5e-3) {
			t.Errorf("costate p(0)[%d] = %v want %v", i, sol.P[0][i], pStart[i])
		}
	}
	// Terminal transversality p(T) = Qf x(T).
	last := len(sol.Times) - 1
	pT := qf.MulVec(sol.X[last])
	for i := range pT {
		if !approx(sol.P[last][i], pT[i], 1e-6) {
			t.Errorf("transversality violated: %v vs %v", sol.P[last][i], pT[i])
		}
	}
}

func TestIndirectShootingLQ(t *testing.T) {
	// Solve a scalar LQ problem by indirect shooting and compare to PontryaginLQ.
	a := FromRows([][]float64{{0.5}})
	b := FromRows([][]float64{{1}})
	q := FromRows([][]float64{{1}})
	r := FromRows([][]float64{{1}})
	qf := FromRows([][]float64{{2}})
	x0 := []float64{1}
	T := 1.5
	ref, _ := PontryaginLQ(a, b, q, r, qf, x0, T, 200)
	f := func(x, u []float64) []float64 {
		return []float64{a.At(0, 0)*x[0] + b.At(0, 0)*u[0]}
	}
	g := func(x, p, u []float64) []float64 {
		return []float64{-q.At(0, 0)*x[0] - a.At(0, 0)*p[0]}
	}
	uOpt := func(x, p []float64) []float64 {
		return []float64{-p[0] / r.At(0, 0)}
	}
	term := func(x []float64) []float64 {
		return []float64{qf.At(0, 0) * x[0]}
	}
	_, sol, err := IndirectShooting(x0, f, g, uOpt, term, T, 400, 50, 1e-12)
	if err != nil {
		t.Fatal(err)
	}
	last := len(sol.Times) - 1
	if !approx(sol.X[last][0], ref.X[len(ref.X)-1][0], 1e-4) {
		t.Errorf("indirect shooting x(T)=%v want %v", sol.X[last][0], ref.X[len(ref.X)-1][0])
	}
}

// ---- Kalman / LQG ------------------------------------------------------------

func TestKalmanContinuous(t *testing.T) {
	a := FromRows([][]float64{{-1}})
	c := FromRows([][]float64{{1}})
	w := FromRows([][]float64{{1}})
	v := FromRows([][]float64{{1}})
	kf, err := KalmanContinuous(a, c, w, v)
	if err != nil {
		t.Fatal(err)
	}
	want := math.Sqrt2 - 1
	if !approx(kf.P.At(0, 0), want, 1e-8) {
		t.Errorf("Kalman P = %v want %v", kf.P.At(0, 0), want)
	}
	if !approx(kf.L.At(0, 0), want, 1e-8) {
		t.Errorf("Kalman L = %v want %v", kf.L.At(0, 0), want)
	}
}

func TestKalmanDiscreteFilter(t *testing.T) {
	a := FromRows([][]float64{{1}})
	c := FromRows([][]float64{{1}})
	w := FromRows([][]float64{{0.01}})
	v := FromRows([][]float64{{0.1}})
	kd, err := KalmanDiscrete(a, c, w, v)
	if err != nil {
		t.Fatal(err)
	}
	// P^2 - W P - W V = 0 -> P = (W + sqrt(W^2+4WV))/2.
	W, V := 0.01, 0.1
	want := (W + math.Sqrt(W*W+4*W*V)) / 2
	if !approx(kd.P.At(0, 0), want, 1e-6) {
		t.Errorf("discrete Kalman P = %v want %v", kd.P.At(0, 0), want)
	}
	// Recursive filter a-priori covariance should approach the same value.
	filt := NewKalmanFilter(a, nil, c, w, v, []float64{0}, FromRows([][]float64{{1}}))
	for i := 0; i < 500; i++ {
		filt.Predict(nil)
		_ = filt.Update([]float64{0})
	}
	filt.Predict(nil)
	if !approx(filt.P.At(0, 0), want, 1e-5) {
		t.Errorf("recursive filter P = %v want %v", filt.P.At(0, 0), want)
	}
}

func TestLQG(t *testing.T) {
	a := FromRows([][]float64{{0, 1}, {0, 0}})
	b := FromRows([][]float64{{0}, {1}})
	c := FromRows([][]float64{{1, 0}})
	q := Eye(2)
	r := Eye(1)
	w := Eye(2)
	v := Eye(1)
	lqg, err := LQGContinuous(a, b, c, q, r, w, v)
	if err != nil {
		t.Fatal(err)
	}
	// Regulator and estimator must both be stabilizing.
	if !IsStableContinuous(a.Minus(b.Mul(lqg.K)), 1e-6) {
		t.Errorf("regulator unstable")
	}
	if !IsStableContinuous(a.Minus(lqg.L.Mul(c)), 1e-6) {
		t.Errorf("estimator unstable")
	}
}

// ---- Simulation / utilities --------------------------------------------------

func TestDiscretizeZOH(t *testing.T) {
	a := FromRows([][]float64{{0, 1}, {0, 0}})
	b := FromRows([][]float64{{0}, {1}})
	ad, bd := DiscretizeZOH(a, b, 0.5)
	wantAd := FromRows([][]float64{{1, 0.5}, {0, 1}})
	wantBd := FromRows([][]float64{{0.125}, {0.5}})
	if !approxMat(ad, wantAd, 1e-9) {
		t.Errorf("Ad = %v", ad)
	}
	if !approxMat(bd, wantBd, 1e-9) {
		t.Errorf("Bd = %v", bd)
	}
}

func TestSimulateDiscrete(t *testing.T) {
	a := FromRows([][]float64{{0.5}})
	b := FromRows([][]float64{{1}})
	sys := NewLinearSystem(a, b, nil, nil)
	traj := sys.SimulateDiscrete([]float64{1}, [][]float64{{0}, {0}})
	if !approx(traj[1][0], 0.5, 1e-12) || !approx(traj[2][0], 0.25, 1e-12) {
		t.Errorf("sim = %v", traj)
	}
}

func TestVectorOps(t *testing.T) {
	if !approx(VecDot([]float64{1, 2}, []float64{3, 4}), 11, 1e-12) {
		t.Errorf("dot wrong")
	}
	if !approx(VecNorm([]float64{3, 4}), 5, 1e-12) {
		t.Errorf("norm wrong")
	}
	s := VecAxpy(2, []float64{1, 1}, []float64{1, 2})
	if !approx(s[0], 3, 1e-12) || !approx(s[1], 4, 1e-12) {
		t.Errorf("axpy wrong: %v", s)
	}
}

func TestKron(t *testing.T) {
	a := FromRows([][]float64{{1, 2}, {3, 4}})
	i2 := Eye(2)
	k := Kron(a, i2)
	if k.Rows() != 4 || k.Cols() != 4 {
		t.Fatalf("kron dims %dx%d", k.Rows(), k.Cols())
	}
	if !approx(k.At(0, 0), 1, 1e-12) || !approx(k.At(0, 2), 2, 1e-12) {
		t.Errorf("kron entries wrong")
	}
}

// ---- Examples ----------------------------------------------------------------

func ExampleLQRContinuous() {
	// Double integrator with unit weights: known gain [1, √3].
	a := FromRows([][]float64{{0, 1}, {0, 0}})
	b := FromRows([][]float64{{0}, {1}})
	q := Eye(2)
	r := Eye(1)
	res, _ := LQRContinuous(a, b, q, r)
	fmt.Printf("K = [%.4f %.4f]\n", res.K.At(0, 0), res.K.At(0, 1))
	// Output: K = [1.0000 1.7321]
}

func ExampleMDP_ValueIteration() {
	stay := FromRows([][]float64{{1, 0}, {0, 1}})
	jump := FromRows([][]float64{{1, 0}, {1, 0}})
	reward := FromRows([][]float64{{1, 1}, {0, -0.1}})
	m, _ := NewMDP([]*Matrix{stay, jump}, reward, 0.5)
	res := m.ValueIteration(1e-12, 1000)
	fmt.Printf("V = [%.2f %.2f], best action in state 1 = %d\n",
		res.Value[0], res.Value[1], res.Policy[1])
	// Output: V = [2.00 0.90], best action in state 1 = 1
}

func ExampleSolveDARE() {
	a := FromRows([][]float64{{1}})
	b := FromRows([][]float64{{1}})
	q := FromRows([][]float64{{1}})
	r := FromRows([][]float64{{1}})
	p, _ := SolveDARE(a, b, q, r)
	fmt.Printf("P = %.4f\n", p.At(0, 0)) // golden ratio
	// Output: P = 1.6180
}

var _ = cmplx.Abs
