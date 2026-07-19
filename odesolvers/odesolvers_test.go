package odesolvers

import (
	"fmt"
	"math"
	"testing"
)

const eulerE = math.E // exp(1)

// expField is the scalar test problem y' = y with solution y(t) = e^t.
func expField(t float64, y []float64) []float64 { return []float64{y[0]} }

// decayField is y' = -y with solution y(t) = e^{-t}.
func decayField(t float64, y []float64) []float64 { return []float64{-y[0]} }

// harmonic is the system y0' = y1, y1' = -y0, solution [cos t, -sin t] from
// [1, 0].
func harmonic(t float64, y []float64) []float64 { return []float64{y[1], -y[0]} }

func approx(t *testing.T, got, want, tol float64, name string) {
	t.Helper()
	if math.Abs(got-want) > tol {
		t.Errorf("%s: got %.12g, want %.12g (|diff|=%.3e > tol %.3e)", name, got, want, math.Abs(got-want), tol)
	}
}

func TestVectorOps(t *testing.T) {
	a := []float64{1, 2, 3}
	b := []float64{4, 5, 6}
	if got := Dot(a, b); got != 32 {
		t.Errorf("Dot = %v, want 32", got)
	}
	if got := Add(a, b); got[0] != 5 || got[2] != 9 {
		t.Errorf("Add = %v", got)
	}
	if got := Sub(b, a); got[0] != 3 || got[2] != 3 {
		t.Errorf("Sub = %v", got)
	}
	approx(t, Norm2([]float64{3, 4}), 5, 1e-12, "Norm2")
	approx(t, Norm1([]float64{-1, 2, -3}), 6, 1e-12, "Norm1")
	approx(t, NormInf([]float64{-1, 5, -3}), 5, 1e-12, "NormInf")
	approx(t, RMSNorm([]float64{3, 4}), 5/math.Sqrt2, 1e-12, "RMSNorm")
	if got := AXPY(2, a, b); got[0] != 6 || got[2] != 12 {
		t.Errorf("AXPY = %v", got)
	}
	ls := Linspace(0, 1, 5)
	approx(t, ls[0], 0, 1e-15, "Linspace0")
	approx(t, ls[4], 1, 1e-15, "Linspace4")
	approx(t, ls[2], 0.5, 1e-15, "Linspace2")
	lc := LinearCombination(3, []float64{1, -1}, [][]float64{a, b})
	if lc[0] != -3 || lc[1] != -3 || lc[2] != -3 {
		t.Errorf("LinearCombination = %v", lc)
	}
}

func TestLinearAlgebra(t *testing.T) {
	A := Matrix{{2, 1, 1}, {1, 3, 2}, {1, 0, 0}}
	b := []float64{4, 5, 6}
	x, err := SolveLinearSystem(A, b)
	if err != nil {
		t.Fatalf("SolveLinearSystem: %v", err)
	}
	// Verify by residual.
	r := Sub(MatVec(A, x), b)
	if NormInf(r) > 1e-10 {
		t.Errorf("residual too large: %v", r)
	}
	// Determinant of a known matrix.
	D := Matrix{{1, 2}, {3, 4}}
	approx(t, Determinant(D), -2, 1e-12, "Determinant")
	// Singular system.
	if _, err := SolveLinearSystem(Matrix{{1, 1}, {1, 1}}, []float64{1, 2}); err == nil {
		t.Errorf("expected singular-matrix error")
	}
	// Matrix product and transpose.
	P := MatMul(Matrix{{1, 2}, {3, 4}}, IdentityMatrix(2))
	if P[0][0] != 1 || P[1][1] != 4 {
		t.Errorf("MatMul identity = %v", P)
	}
	Tr := Transpose(Matrix{{1, 2, 3}, {4, 5, 6}})
	if len(Tr) != 3 || Tr[2][1] != 6 {
		t.Errorf("Transpose = %v", Tr)
	}
	approx(t, Trace(Matrix{{1, 0}, {0, 5}}), 6, 1e-12, "Trace")
}

func TestNewtonSolve(t *testing.T) {
	// Solve x^2 - 2 = 0, y - 3 = 0.
	F := func(v []float64) []float64 {
		return []float64{v[0]*v[0] - 2, v[1] - 3}
	}
	sol, err := NewtonSolve(F, []float64{1, 1}, 1e-12, 50)
	if err != nil {
		t.Fatalf("NewtonSolve: %v", err)
	}
	approx(t, sol[0], math.Sqrt2, 1e-8, "Newton x")
	approx(t, sol[1], 3, 1e-8, "Newton y")
}

func TestExplicitFixed(t *testing.T) {
	cases := []struct {
		name  string
		solve func(Field, float64, []float64, float64, float64) *Solution
		h     float64
		tol   float64
	}{
		{"Euler", SolveEuler, 1e-4, 2e-4},
		{"Midpoint", SolveMidpoint, 1e-2, 1e-4},
		{"Heun", SolveHeun, 1e-2, 1e-4},
		{"Ralston", SolveRalston, 1e-2, 1e-4},
		{"SSPRK3", SolveSSPRK3, 1e-2, 1e-6},
		{"Kutta3", SolveKutta3, 1e-2, 1e-6},
		{"RK4", SolveRK4, 1e-1, 1e-5},
		{"RK38", SolveRK38, 1e-1, 1e-5},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			sol := c.solve(expField, 0, []float64{1}, 1, c.h)
			approx(t, sol.Final()[0], eulerE, c.tol, c.name)
			if sol.Steps() < 1 {
				t.Errorf("no steps taken")
			}
		})
	}
}

func TestExplicitStepFuncs(t *testing.T) {
	// One RK4 step of y'=y from y=1 with h=0.5 matches the Taylor value well.
	y := RK4Step(expField, 0, []float64{1}, 0.5)
	approx(t, y[0], math.Exp(0.5), 1e-3, "RK4Step")
	// Euler step is exact to first order.
	ye := EulerStep(expField, 0, []float64{1}, 0.1)
	approx(t, ye[0], 1.1, 1e-15, "EulerStep")
	// Midpoint step.
	ym := MidpointStep(expField, 0, []float64{1}, 0.1)
	approx(t, ym[0], 1+0.1+0.005, 1e-15, "MidpointStep")
}

func TestAdaptive(t *testing.T) {
	cases := []struct {
		name  string
		solve func(Field, float64, []float64, float64, AdaptiveOptions) (*Solution, error)
		tol   float64
	}{
		{"HeunEuler", SolveHeunEuler, 1e-3},
		{"BogackiShampine", SolveBogackiShampine, 1e-4},
		{"RKF45", SolveRKF45, 1e-5},
		{"CashKarp", SolveCashKarp, 1e-5},
		{"DOPRI5", SolveDOPRI5, 1e-5},
		{"RKF78", SolveRKF78, 1e-6},
	}
	opts := DefaultAdaptiveOptions()
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			sol, err := c.solve(expField, 0, []float64{1}, 1, opts)
			if err != nil {
				t.Fatalf("%s: %v", c.name, err)
			}
			approx(t, sol.Final()[0], eulerE, c.tol, c.name)
			if sol.Accepted == 0 {
				t.Errorf("%s: no accepted steps recorded", c.name)
			}
		})
	}
}

func TestAdaptiveBackward(t *testing.T) {
	// Integrate y'=y from t=1 back to t=0; from y(1)=e we should recover y(0)=1.
	opts := DefaultAdaptiveOptions()
	sol, err := SolveDOPRI5(expField, 1, []float64{math.E}, 0, opts)
	if err != nil {
		t.Fatalf("backward DOPRI5: %v", err)
	}
	approx(t, sol.Final()[0], 1, 1e-5, "backward integration")
	if sol.FinalTime() != 0 {
		t.Errorf("final time = %v, want 0", sol.FinalTime())
	}
}

func TestHarmonicSystem(t *testing.T) {
	opts := DefaultAdaptiveOptions()
	opts.RelTol = 1e-9
	opts.AbsTol = 1e-12
	sol, err := SolveDOPRI5(harmonic, 0, []float64{1, 0}, 2*math.Pi, opts)
	if err != nil {
		t.Fatalf("harmonic: %v", err)
	}
	fin := sol.Final()
	approx(t, fin[0], 1, 1e-6, "cos(2pi)")
	approx(t, fin[1], 0, 1e-6, "-sin(2pi)")
	// Dense output midway should give cos/-sin at t=pi/2.
	mid := sol.At(math.Pi / 2)
	approx(t, mid[0], 0, 1e-3, "dense cos(pi/2)")
	approx(t, mid[1], -1, 1e-3, "dense -sin(pi/2)")
}

func TestImplicit(t *testing.T) {
	cases := []struct {
		name  string
		solve func(Field, float64, []float64, float64, float64) (*Solution, error)
		h     float64
		tol   float64
	}{
		{"BackwardEuler", SolveBackwardEuler, 1e-3, 5e-3},
		{"Trapezoidal", SolveTrapezoidal, 1e-2, 1e-4},
		{"ImplicitMidpoint", SolveImplicitMidpoint, 1e-2, 1e-4},
		{"GaussLegendre4", SolveGaussLegendre4, 5e-2, 1e-6},
		{"RadauIIA3", SolveRadauIIA3, 1e-2, 1e-5},
		{"RadauIIA5", SolveRadauIIA5, 1e-1, 1e-7},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			sol, err := c.solve(expField, 0, []float64{1}, 1, c.h)
			if err != nil {
				t.Fatalf("%s: %v", c.name, err)
			}
			approx(t, sol.Final()[0], eulerE, c.tol, c.name)
		})
	}
}

func TestStiffDecay(t *testing.T) {
	// A moderately stiff scalar problem y' = -50 y; backward Euler is stable at
	// a step where forward Euler would blow up.
	stiff := func(t float64, y []float64) []float64 { return []float64{-50 * y[0]} }
	sol, err := SolveBackwardEuler(stiff, 0, []float64{1}, 1, 0.1)
	if err != nil {
		t.Fatalf("stiff: %v", err)
	}
	// Solution must decay monotonically toward 0 and stay bounded.
	prev := math.Inf(1)
	for _, y := range sol.Y {
		if y[0] < 0 || y[0] > prev+1e-12 {
			t.Errorf("backward Euler not stable/monotone: %v", sol.Component(0))
			break
		}
		prev = y[0]
	}
	approx(t, sol.Final()[0], math.Exp(-50), 1e-2, "stiff final")
}

func TestBDF(t *testing.T) {
	for order := 1; order <= 6; order++ {
		sol, err := SolveBDF(decayField, order, 0, []float64{1}, 1, 0.005)
		if err != nil {
			t.Fatalf("BDF%d: %v", order, err)
		}
		approx(t, sol.Final()[0], math.Exp(-1), 5e-3, fmt.Sprintf("BDF%d", order))
	}
	// Coefficient sanity: BDF2 alpha and beta.
	alpha, beta, err := BDFCoefficients(2)
	if err != nil {
		t.Fatal(err)
	}
	approx(t, alpha[0], 1, 1e-15, "BDF2 a0")
	approx(t, alpha[1], -4.0/3.0, 1e-15, "BDF2 a1")
	approx(t, beta, 2.0/3.0, 1e-15, "BDF2 beta")
	if _, _, err := BDFCoefficients(7); err == nil {
		t.Errorf("BDF7 should be invalid")
	}
}

func TestMultistep(t *testing.T) {
	// Adams-Bashforth orders.
	for order := 1; order <= 5; order++ {
		sol, err := SolveAdamsBashforth(expField, order, 0, []float64{1}, 1, 0.005)
		if err != nil {
			t.Fatalf("AB%d: %v", order, err)
		}
		approx(t, sol.Final()[0], eulerE, 1e-2, fmt.Sprintf("AB%d", order))
	}
	// Predictor-corrector is markedly more accurate.
	sol, err := SolveABM(expField, 4, 2, 0, []float64{1}, 1, 0.01)
	if err != nil {
		t.Fatalf("ABM4: %v", err)
	}
	approx(t, sol.Final()[0], eulerE, 1e-6, "ABM4")
	// Coefficient checks.
	ab, _ := AdamsBashforthCoefficients(2)
	approx(t, ab[0], 1.5, 1e-15, "AB2 b0")
	approx(t, ab[1], -0.5, 1e-15, "AB2 b1")
	bNew, bOld, _ := AdamsMoultonCoefficients(2)
	approx(t, bNew, 0.5, 1e-15, "AM2 new")
	approx(t, bOld[0], 0.5, 1e-15, "AM2 old")
}

func TestSymplecticEnergy(t *testing.T) {
	acc := func(t float64, q []float64) []float64 { return []float64{-q[0]} }
	q0, v0 := []float64{1}, []float64{0}
	e0 := HarmonicEnergy(q0, v0, 1)
	periods := 10.0
	tEnd := periods * 2 * math.Pi

	methods := []struct {
		name  string
		solve func(Acceleration, float64, []float64, []float64, float64, float64) *SymplecticSolution
		drift float64
	}{
		{"VelocityVerlet", SolveVelocityVerlet, 1e-3},
		{"PositionVerlet", SolvePositionVerlet, 1e-3},
		{"Leapfrog", SolveLeapfrog, 1e-3},
		{"Yoshida4", SolveYoshida4, 1e-6},
	}
	for _, m := range methods {
		t.Run(m.name, func(t *testing.T) {
			sol := m.solve(acc, 0, q0, v0, tEnd, 0.02)
			ef := HarmonicEnergy(sol.FinalPosition(), sol.FinalVelocity(), 1)
			if d := math.Abs(ef - e0); d > m.drift {
				t.Errorf("%s energy drift %.3e exceeds %.3e", m.name, d, m.drift)
			}
			// After a whole number of periods the position returns near q0.
			approx(t, sol.FinalPosition()[0], 1, 1e-2, m.name+" position")
		})
	}
	// Yoshida composition weights sum to 1 (w1 + w0 + w1).
	w1, w0 := YoshidaCoefficients()
	approx(t, 2*w1+w0, 1, 1e-12, "Yoshida weight sum")
}

func TestEvents(t *testing.T) {
	opts := AdaptiveOptions{RelTol: 1e-9, AbsTol: 1e-12, Dense: true}
	sol, err := IntegrateAdaptive(harmonic, DormandPrinceTableau(), 0, []float64{1, 0}, 10, opts)
	if err != nil {
		t.Fatalf("integrate: %v", err)
	}
	g := func(t float64, y []float64) float64 { return y[0] } // cos t
	all := ScanEvents(sol, g, 0)
	// Zeros of cos in (0,10): pi/2, 3pi/2, 5pi/2 ... up to ~3.14 periods.
	if len(all) < 3 {
		t.Fatalf("expected several crossings, got %d", len(all))
	}
	approx(t, all[0].T, math.Pi/2, 1e-5, "first crossing")
	approx(t, all[1].T, 3*math.Pi/2, 1e-5, "second crossing")
	// Directional filter: first downward crossing is at pi/2 (cos goes + to -).
	first, ok := FirstEvent(sol, g, -1)
	if !ok {
		t.Fatal("expected a downward crossing")
	}
	approx(t, first.T, math.Pi/2, 1e-5, "first downward")
	if first.Direction != -1 {
		t.Errorf("direction = %d, want -1", first.Direction)
	}
	// IntegrateUntilEvent stops at the first crossing.
	_, ev, found := IntegrateUntilEvent(harmonic, g, -1, 0, []float64{1, 0}, 10, 0.01)
	if !found {
		t.Fatal("IntegrateUntilEvent found no event")
	}
	approx(t, ev.T, math.Pi/2, 1e-3, "until-event time")
}

func TestBVPSingleShooting(t *testing.T) {
	// y'' = -y, y(0) = 0, y(pi/2) = 1  =>  y = sin, y'(0) = 1.
	bc := func(ya, yb []float64) []float64 {
		return []float64{ya[0], yb[0] - 1}
	}
	s0, sol, err := SingleShooting(harmonic, bc, 0, math.Pi/2, []float64{0, 0}, DefaultShootingOptions())
	if err != nil {
		t.Fatalf("SingleShooting: %v", err)
	}
	approx(t, s0[0], 0, 1e-6, "y(0)")
	approx(t, s0[1], 1, 1e-6, "y'(0)")
	// Trajectory should match sin at the midpoint.
	approx(t, sol.At(math.Pi / 4)[0], math.Sin(math.Pi/4), 1e-4, "sin(pi/4)")
}

func TestBVPMultipleShooting(t *testing.T) {
	bc := func(ya, yb []float64) []float64 {
		return []float64{ya[0], yb[0] - 1}
	}
	nodes, states, _, err := MultipleShooting(harmonic, bc, 0, math.Pi/2, 4, []float64{0, 0}, DefaultShootingOptions())
	if err != nil {
		t.Fatalf("MultipleShooting: %v", err)
	}
	if len(nodes) != 5 {
		t.Fatalf("expected 5 nodes, got %d", len(nodes))
	}
	approx(t, states[0][1], 1, 1e-5, "y'(0)")
	// Each node value should equal sin(node time).
	for k, tk := range nodes {
		approx(t, states[k][0], math.Sin(tk), 1e-4, fmt.Sprintf("node %d", k))
	}
}

func TestDenseOutputHermite(t *testing.T) {
	// Cubic Hermite must reproduce endpoint values and derivatives exactly.
	y0 := []float64{2}
	y1 := []float64{5}
	f0 := []float64{1}
	f1 := []float64{-1}
	h := 2.0
	approx(t, HermiteInterpolate(0, h, y0, y1, f0, f1)[0], 2, 1e-14, "Hermite theta=0")
	approx(t, HermiteInterpolate(1, h, y0, y1, f0, f1)[0], 5, 1e-14, "Hermite theta=1")
	// A fixed-step dense RK4 solution interpolates y'=y accurately off-node.
	sol := IntegrateFixedDense(expField, RK4Tableau(), 0, []float64{1}, 1, 0.1)
	approx(t, sol.At(0.55)[0], math.Exp(0.55), 1e-5, "dense RK4 interp")
}

func TestRichardsonAndOrder(t *testing.T) {
	// Two Euler solves at h and h/2 extrapolate to second order.
	coarse := SolveEuler(expField, 0, []float64{1}, 1, 0.02).Final()
	fine := SolveEuler(expField, 0, []float64{1}, 1, 0.01).Final()
	ex := RichardsonExtrapolate(coarse, fine, 2, 1)
	if math.Abs(ex[0]-eulerE) >= math.Abs(fine[0]-eulerE) {
		t.Errorf("Richardson did not improve accuracy: ex=%.8f fine=%.8f", ex[0], fine[0])
	}
	// Empirical order of RK4 near 4.
	e := func(h float64) float64 {
		return math.Abs(SolveRK4(expField, 0, []float64{1}, 1, h).Final()[0] - eulerE)
	}
	ord := EstimateOrder(e(0.2), e(0.1), e(0.05))
	if math.Abs(ord-4) > 0.5 {
		t.Errorf("estimated RK4 order = %.3f, want ~4", ord)
	}
}

func TestTableauConsistency(t *testing.T) {
	tabs := []*ButcherTableau{
		EulerTableau(), MidpointTableau(), HeunTableau(), RalstonTableau(),
		SSPRK3Tableau(), RK4Tableau(), RK38Tableau(), KuttaThirdOrderTableau(),
		HeunThirdOrderTableau(), HeunEulerTableau(), BogackiShampineTableau(),
		FehlbergTableau(), CashKarpTableau(), DormandPrinceTableau(),
		Fehlberg78Tableau(), BackwardEulerTableau(), TrapezoidalTableau(),
		ImplicitMidpointTableau(), GaussLegendre4Tableau(), RadauIIA3Tableau(),
		RadauIIA5Tableau(),
	}
	for _, bt := range tabs {
		if !bt.IsConsistent(1e-10) {
			t.Errorf("%s is not consistent (sumB=%.15g, rowSums=%v, C=%v)",
				bt.Name, bt.SumB(), bt.RowSums(), bt.C)
		}
	}
	// Explicitness classification.
	if !RK4Tableau().IsExplicit() {
		t.Errorf("RK4 should be explicit")
	}
	if BackwardEulerTableau().IsExplicit() {
		t.Errorf("backward Euler should be implicit")
	}
	if !DormandPrinceTableau().IsEmbedded() {
		t.Errorf("DOPRI5 should be embedded")
	}
}

// ExampleSolveRK4 integrates the logistic equation y' = y (1 - y) with the
// classical fourth-order Runge-Kutta method and prints the value at t = 5,
// which approaches the carrying capacity 1.
func ExampleSolveRK4() {
	logistic := func(t float64, y []float64) []float64 {
		return []float64{y[0] * (1 - y[0])}
	}
	sol := SolveRK4(logistic, 0, []float64{0.1}, 5, 0.1)
	fmt.Printf("y(5) = %.5f\n", sol.Final()[0])
	// Output: y(5) = 0.94283
}

// ExampleSolveDOPRI5 shows adaptive integration of the exponential-growth
// problem y' = y and recovery of Euler's number at t = 1.
func ExampleSolveDOPRI5() {
	f := func(t float64, y []float64) []float64 { return []float64{y[0]} }
	sol, _ := SolveDOPRI5(f, 0, []float64{1}, 1, DefaultAdaptiveOptions())
	fmt.Printf("e ~= %.6f\n", sol.Final()[0])
	// Output: e ~= 2.718282
}
