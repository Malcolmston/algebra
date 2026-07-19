package chaos

import (
	"fmt"
	"math"
	"math/cmplx"
	"math/rand"
	"testing"
)

func approx(t *testing.T, got, want, tol float64, name string) {
	t.Helper()
	if math.Abs(got-want) > tol {
		t.Errorf("%s: got %g, want %g (tol %g)", name, got, want, tol)
	}
}

// --- linear algebra ---

func TestVecOps(t *testing.T) {
	a := VecOf(3, 4)
	approx(t, a.Norm(), 5, 1e-12, "norm")
	approx(t, a.Norm1(), 7, 1e-12, "norm1")
	approx(t, a.NormInf(), 4, 1e-12, "norminf")
	b := VecOf(1, 2)
	approx(t, a.Dot(b), 11, 1e-12, "dot")
	approx(t, a.Sub(b).Norm(), math.Hypot(2, 2), 1e-12, "distance")
	approx(t, a.AddScaled(2, b).Norm(), VecOf(5, 8).Norm(), 1e-12, "addscaled")
	u := a.Normalize()
	approx(t, u.Norm(), 1, 1e-12, "normalize")
	approx(t, a.Mean(), 3.5, 1e-12, "mean")
	approx(t, a.Max(), 4, 1e-12, "max")
	approx(t, a.Min(), 3, 1e-12, "min")
}

func TestMatBasics(t *testing.T) {
	A := MatFromRows(VecOf(1, 2), VecOf(3, 4))
	approx(t, A.Trace(), 5, 1e-12, "trace")
	approx(t, Det2(A), -2, 1e-12, "det2")
	approx(t, Det(A), -2, 1e-12, "det")
	I := Eye(2)
	P := A.Mul(I)
	approx(t, P.Sub(A).FrobeniusNorm(), 0, 1e-12, "mul identity")
	AT := A.Transpose()
	approx(t, AT[0][1], 3, 1e-12, "transpose")
}

func TestDet3(t *testing.T) {
	A := MatFromRows(VecOf(2, -3, 1), VecOf(2, 0, -1), VecOf(1, 4, 5))
	approx(t, Det3(A), 49, 1e-9, "det3")
	approx(t, Det(A), 49, 1e-9, "det gauss")
}

func TestSolveAndInverse(t *testing.T) {
	A := MatFromRows(VecOf(2, 1), VecOf(1, 3))
	b := VecOf(3, 5)
	x, err := SolveLinear(A, b)
	if err != nil {
		t.Fatal(err)
	}
	// 2x+y=3, x+3y=5 => x=0.8, y=1.4
	approx(t, x[0], 0.8, 1e-12, "solve x")
	approx(t, x[1], 1.4, 1e-12, "solve y")
	inv, err := Inverse(A)
	if err != nil {
		t.Fatal(err)
	}
	prod := A.Mul(inv)
	approx(t, prod.Sub(Eye(2)).FrobeniusNorm(), 0, 1e-12, "A*inv=I")
	if _, err := SolveLinear(MatFromRows(VecOf(1, 2), VecOf(2, 4)), VecOf(1, 1)); err != ErrSingular {
		t.Errorf("expected ErrSingular, got %v", err)
	}
}

func TestQR(t *testing.T) {
	A := MatFromRows(VecOf(12, -51, 4), VecOf(6, 167, -68), VecOf(-4, 24, -41))
	Q, R := QR(A)
	// Q orthonormal: Q^T Q = I.
	QtQ := Q.Transpose().Mul(Q)
	approx(t, QtQ.Sub(Eye(3)).FrobeniusNorm(), 0, 1e-10, "Q orthonormal")
	// Q R = A.
	approx(t, Q.Mul(R).Sub(A).FrobeniusNorm(), 0, 1e-9, "QR=A")
}

func TestEigenvalues(t *testing.T) {
	// Diagonal.
	D := MatFromRows(VecOf(2, 0), VecOf(0, 3))
	e := Eigenvalues2(D)
	vals := []float64{real(e[0]), real(e[1])}
	if !(nearEither(vals, 2) && nearEither(vals, 3)) {
		t.Errorf("eig2 diagonal got %v", e)
	}
	// Rotation-like with complex eigenvalues [[0,-1],[1,0]] -> +/- i.
	R := MatFromRows(VecOf(0, -1), VecOf(1, 0))
	er := Eigenvalues2(R)
	approx(t, math.Abs(imag(er[0])), 1, 1e-12, "rot imag")
	approx(t, real(er[0]), 0, 1e-12, "rot real")
	// 3x3 diagonal.
	D3 := MatFromRows(VecOf(1, 0, 0), VecOf(0, 5, 0), VecOf(0, 0, -2))
	e3 := Eigenvalues3(D3)
	sr := SpectralRadius(D3)
	approx(t, sr, 5, 1e-9, "spectral radius")
	sum := real(e3[0]) + real(e3[1]) + real(e3[2])
	approx(t, sum, 4, 1e-9, "eig3 trace")
}

func nearEither(vals []float64, want float64) bool {
	for _, v := range vals {
		if math.Abs(v-want) < 1e-9 {
			return true
		}
	}
	return false
}

func TestSolveCubic(t *testing.T) {
	// (x-1)(x-2)(x-3) = x^3 -6x^2 +11x -6.
	r := solveCubic(1, -6, 11, -6)
	got := []float64{real(r[0]), real(r[1]), real(r[2])}
	for _, want := range []float64{1, 2, 3} {
		if !nearEither(got, want) {
			t.Errorf("cubic missing root %v in %v", want, got)
		}
	}
	for _, root := range r {
		approx(t, imag(root), 0, 1e-9, "cubic imag")
	}
}

// --- 1D maps ---

func TestLogisticFixedPoints(t *testing.T) {
	fp := LogisticFixedPoints(4)
	approx(t, fp[0], 0, 1e-12, "fp0")
	approx(t, fp[1], 0.75, 1e-12, "fp1")
	// Refine numerically.
	r, err := FixedPoint1D(Logistic(3.2), 0.6, 1e-12, 100)
	if err != nil {
		t.Fatal(err)
	}
	approx(t, r, 1-1/3.2, 1e-9, "fixedpoint refine")
	approx(t, Multiplier1D(Logistic(3.2), r), 3.2*(1-2*(1-1/3.2)), 1e-6, "multiplier")
}

func TestMapStability(t *testing.T) {
	// Logistic fixed point 1-1/r stable for 1<r<3.
	f := Logistic(2.5)
	if !IsStable1D(f, 1-1/2.5) {
		t.Error("expected stable fixed point at r=2.5")
	}
	if IsStable1D(Logistic(3.5), 1-1/3.5) {
		t.Error("expected unstable fixed point at r=3.5")
	}
}

func TestOrbitAndCompose(t *testing.T) {
	f := Logistic(2)
	// Fixed point of r=2 is 0.5; starting there stays.
	o := Orbit(f, 0.5, 5)
	for _, x := range o {
		approx(t, x, 0.5, 1e-12, "orbit fp")
	}
	g := Compose(f, 3)
	approx(t, g(0.5), Iterate(f, 0.5, 3), 1e-12, "compose")
	xs, ys := Cobweb(f, 0.3, 4)
	if len(xs) != 9 || len(ys) != 9 {
		t.Errorf("cobweb length got %d,%d", len(xs), len(ys))
	}
}

// --- Lyapunov ---

func TestLyapunov1D(t *testing.T) {
	cases := []struct {
		name string
		f    Map1D
		x0   float64
		want float64
		tol  float64
	}{
		{"logistic r=4", Logistic(4), 0.123, math.Ln2, 0.02},
		{"tent mu=1", Tent(1.0), 0.11, math.Ln2, 1e-3},
	}
	for _, c := range cases {
		got := Lyapunov1D(c.f, c.x0, 2000, 40000)
		approx(t, got, c.want, c.tol, "lyap "+c.name)
	}
	// Analytic derivative form on Bernoulli map: derivative 2 everywhere.
	lam := LyapunovLog1D(BernoulliMap(), func(float64) float64 { return 2 }, 0.1234567, 1000, 10000)
	approx(t, lam, math.Ln2, 1e-9, "bernoulli analytic")
	// Separation method should agree for logistic r=4.
	sep := LyapunovSeparation1D(Logistic(4), 0.2, 1e-8, 2000, 40000)
	approx(t, sep, math.Ln2, 0.03, "separation")
}

func TestHenonSpectrum(t *testing.T) {
	H := Lift2D(HenonMap(1.4, 0.3))
	sp := LyapunovSpectrumMap(H, VecOf(0, 0), 2000, 40000)
	if sp[0] < sp[1] {
		t.Error("spectrum not sorted descending")
	}
	approx(t, sp[0], 0.419, 0.01, "henon largest")
	// Sum of exponents equals log|det J| = log(b) = log(0.3).
	approx(t, sp.Sum(), math.Log(0.3), 0.02, "henon sum")
	// Largest agrees with Benettin single-vector method.
	bl := BenettinLargestMap(H, VecOf(0, 0), 2000, 40000)
	approx(t, bl, sp[0], 0.01, "benettin vs qr")
	ky := KaplanYorkeDimension(sp)
	approx(t, ky, 1.258, 0.01, "henon KY")
}

func TestLorenzSpectrum(t *testing.T) {
	f := LorenzStandard()
	sp := LyapunovSpectrumFlow(f, VecOf(1, 1, 1), 0.01, 3000, 30000, 1)
	// Sum equals trace of Jacobian = -(sigma+1+beta).
	wantSum := -(10.0 + 1.0 + 8.0/3.0)
	approx(t, sp.Sum(), wantSum, 0.3, "lorenz sum")
	approx(t, sp[0], 0.9, 0.06, "lorenz largest")
	// Middle exponent close to zero (flow direction).
	approx(t, sp[1], 0, 0.05, "lorenz zero exponent")
	ky := KaplanYorkeDimension(sp)
	approx(t, ky, 2.06, 0.05, "lorenz KY")
	bl := BenettinLargestFlow(f, VecOf(1, 1, 1), 0.01, 3000, 20000)
	approx(t, bl, sp[0], 0.05, "lorenz benettin")
}

func TestKaplanYorke(t *testing.T) {
	// Two exponents 0.42, -1.62: D = 1 + 0.42/1.62.
	ky := KaplanYorkeDimension(VecOf(0.42, -1.62))
	approx(t, ky, 1+0.42/1.62, 1e-12, "KY explicit")
	// All negative -> 0.
	approx(t, KaplanYorkeDimension(VecOf(-0.1, -0.2)), 0, 1e-12, "KY all neg")
	// All non-negative -> full dimension.
	approx(t, KaplanYorkeDimension(VecOf(0.1, 0.2)), 2, 1e-12, "KY all pos")
	approx(t, MetricEntropy(VecOf(0.9, 0, -14)), 0.9, 1e-12, "pesin entropy")
}

// --- flows / integrators ---

func TestIntegratorExponential(t *testing.T) {
	// dx/dt = x, x(0)=1 -> x(1)=e.
	f := Field(func(v Vec) Vec { return Vec{v[0]} })
	got := Integrate(f, VecOf(1), 0.001, 1000)
	approx(t, got[0], math.E, 1e-4, "RK4 exp")
	// Euler is less accurate but in the right direction.
	e := StepEuler(f, VecOf(1), 0.1)
	approx(t, e[0], 1.1, 1e-12, "euler step")
	m := StepRK2(f, VecOf(1), 0.1)
	approx(t, m[0], 1.105, 1e-12, "rk2 step")
}

func TestHarmonicOscillator(t *testing.T) {
	// x''=-x -> energy conserved. State (x, v).
	f := Field(func(u Vec) Vec { return Vec{u[1], -u[0]} })
	traj := Trajectory(f, VecOf(1, 0), 0.001, 6283) // ~ one period 2pi
	last := traj[len(traj)-1]
	approx(t, last[0], 1, 1e-3, "oscillator returns")
	approx(t, last[1], 0, 1e-3, "oscillator velocity")
	ts := Times(0.001, 6283)
	approx(t, ts[len(ts)-1], 6.283, 1e-6, "times")
}

func TestLorenzEquilibriaAndClassify(t *testing.T) {
	eq := LorenzEquilibria(28, 8.0/3.0)
	if len(eq) != 3 {
		t.Fatalf("expected 3 equilibria, got %d", len(eq))
	}
	// Origin is a saddle for rho>1.
	J0 := LorenzJacobian(10, 28, 8.0/3.0, eq[0])
	if k := ClassifyFlow(J0); k != Saddle {
		t.Errorf("origin should be saddle, got %v", k)
	}
	// Verify equilibrium refinement lands on a genuine root.
	root, err := Equilibrium(LorenzStandard(), VecOf(7, 7, 26), 1e-10, 100)
	if err != nil {
		t.Fatal(err)
	}
	approx(t, LorenzStandard()(root).Norm(), 0, 1e-8, "equilibrium residual")
}

func TestClassifyFlowKinds(t *testing.T) {
	// Stable node: diag(-1,-2).
	if k := ClassifyFlow(MatFromRows(VecOf(-1, 0), VecOf(0, -2))); k != StableNode {
		t.Errorf("stable node got %v", k)
	}
	// Unstable node.
	if k := ClassifyFlow(MatFromRows(VecOf(1, 0), VecOf(0, 2))); k != UnstableNode {
		t.Errorf("unstable node got %v", k)
	}
	// Stable focus: [[-0.1,-1],[1,-0.1]].
	if k := ClassifyFlow(MatFromRows(VecOf(-0.1, -1), VecOf(1, -0.1))); k != StableFocus {
		t.Errorf("stable focus got %v", k)
	}
	// Center: [[0,-1],[1,0]].
	if k := ClassifyFlow(MatFromRows(VecOf(0, -1), VecOf(1, 0))); k != Center {
		t.Errorf("center got %v", k)
	}
	if !IsHyperbolicFlow(MatFromRows(VecOf(-1, 0), VecOf(0, -2))) {
		t.Error("expected hyperbolic")
	}
	if k := StableNode.String(); k != "stable node" {
		t.Errorf("string got %q", k)
	}
}

func TestClassifyMap(t *testing.T) {
	// Henon fixed point is a saddle.
	a, b := 1.4, 0.3
	// Fixed point x* solves x = 1 - a x^2 + b x.
	F := Lift2D(HenonMap(a, b))
	fp, err := FixedPointMap(F, VecOf(0.6, 0.2), 1e-12, 100)
	if err != nil {
		t.Fatal(err)
	}
	J := HenonJacobian(a, b, fp[0], fp[1])
	if k := ClassifyMap(J); k != Saddle {
		t.Errorf("henon fixed point should be saddle, got %v", k)
	}
	// Contraction inside unit circle -> stable node.
	if k := ClassifyMap(MatFromRows(VecOf(0.5, 0), VecOf(0, 0.3))); k != StableNode {
		t.Errorf("stable node map got %v", k)
	}
}

// --- 2D maps ---

func TestStandardMapAreaPreserving(t *testing.T) {
	for _, th := range []float64{0.1, 1.0, 2.5} {
		J := StandardMapJacobian(0.9, th, 0)
		approx(t, Det2(J), 1, 1e-12, "standard map det")
	}
	// One orbit stays bounded in p for small K on an invariant curve.
	f := StandardMap(0.5)
	th, p := 0.1, 0.2
	for i := 0; i < 100; i++ {
		th, p = f(th, p)
		if th < 0 || th >= 2*math.Pi+1e-9 {
			t.Errorf("theta out of range: %g", th)
		}
	}
}

func TestHenonJacobianNumeric(t *testing.T) {
	a, b := 1.4, 0.3
	F := Lift2D(HenonMap(a, b))
	num := JacobianMap(F, VecOf(0.3, 0.1), 1e-6)
	ana := HenonJacobian(a, b, 0.3, 0.1)
	approx(t, num.Sub(ana).FrobeniusNorm(), 0, 1e-4, "henon jacobian")
}

// --- bifurcation / Feigenbaum ---

func TestDetectPeriod(t *testing.T) {
	// r=3.2 -> period 2.
	if p := DetectPeriod(Logistic(3.2), 0.4, 5000, 16, 1e-6); p != 2 {
		t.Errorf("expected period 2, got %d", p)
	}
	// r=3.5 -> period 4.
	if p := DetectPeriod(Logistic(3.5), 0.4, 5000, 16, 1e-6); p != 4 {
		t.Errorf("expected period 4, got %d", p)
	}
	// r=2.8 -> fixed point (period 1).
	if p := DetectPeriod(Logistic(2.8), 0.4, 5000, 16, 1e-8); p != 1 {
		t.Errorf("expected period 1, got %d", p)
	}
}

func TestBifurcationDiagram(t *testing.T) {
	d := BifurcationDiagram(LogisticFamily(), 2.5, 4.0, 50, 0.3, 500, 64)
	if len(d) != 50 {
		t.Fatalf("expected 50 points, got %d", len(d))
	}
	// At r=2.9 the attractor is a single fixed point.
	set := AttractorSet(Logistic(2.9), 0.3, 2000, 200, 1e-6)
	if len(set) != 1 {
		t.Errorf("expected single attractor point, got %d", len(set))
	}
	approx(t, set[0], 1-1/2.9, 1e-4, "attractor value")
}

func TestFeigenbaum(t *testing.T) {
	seq := SuperstableSequence(6)
	// Known superstable parameters.
	wants := []float64{2.0, 3.236068, 3.498562, 3.554641, 3.566667, 3.569244}
	for i, w := range wants {
		approx(t, seq[i], w, 1e-4, fmt.Sprintf("superstable s%d", i))
	}
	ratios := FeigenbaumFromSequence(seq)
	// Later ratios approach the Feigenbaum constant.
	last := ratios[len(ratios)-1]
	approx(t, last, FeigenbaumDelta, 0.05, "feigenbaum ratio")
	approx(t, FeigenbaumEstimate(1, 2, 2.2), 5, 1e-12, "feigenbaum estimate")
}

func TestLyapunovDiagram(t *testing.T) {
	_, lam := LyapunovExponentDiagram(LogisticFamily(), 2.8, 4.0, 40, 0.3, 1000, 2000)
	// At r=4 (last) exponent positive; at r=2.8 (first) negative.
	if lam[0] >= 0 {
		t.Errorf("expected negative exponent at r=2.8, got %g", lam[0])
	}
	if lam[len(lam)-1] <= 0 {
		t.Errorf("expected positive exponent near r=4, got %g", lam[len(lam)-1])
	}
}

// --- Poincare ---

func TestPoincareSection(t *testing.T) {
	f := LorenzStandard()
	sec := PlaneSection(3, 2, 27, 1) // z = 27, upward crossings
	cr, err := PoincareSection(f, VecOf(1, 1, 1), 0.01, 3000, 20000, sec)
	if err != nil {
		t.Fatal(err)
	}
	if len(cr) < 50 {
		t.Errorf("too few crossings: %d", len(cr))
	}
	// Every crossing lies on the section plane.
	for _, p := range cr {
		approx(t, p[2], 27, 1e-2, "crossing on plane")
	}
	in, out := ReturnMap(cr, 0)
	if len(in) != len(cr)-1 || len(out) != len(cr)-1 {
		t.Error("return map length mismatch")
	}
	// No crossing case.
	if _, err := PoincareSection(f, VecOf(1, 1, 1), 0.01, 0, 10, PlaneSection(3, 2, 1000, 1)); err != ErrNoCrossing {
		t.Errorf("expected ErrNoCrossing, got %v", err)
	}
}

func TestLorenzReturnMap(t *testing.T) {
	f := LorenzStandard()
	maxima := LocalMaxima(f, VecOf(1, 1, 1), 0.005, 4000, 40000, 2)
	if len(maxima) < 50 {
		t.Fatalf("too few maxima: %d", len(maxima))
	}
	in, out := SuccessiveMaximaMap(maxima)
	if len(in) != len(maxima)-1 {
		t.Error("maxima map length")
	}
	// Classic Lorenz map: z-maxima lie above ~30.
	for _, m := range maxima {
		if m < 25 {
			t.Errorf("unexpected small maximum %g", m)
		}
	}
	_ = out
}

// --- dimensions ---

func TestFitLine(t *testing.T) {
	xs := []float64{0, 1, 2, 3, 4}
	ys := []float64{1, 3, 5, 7, 9}
	m, b := FitLine(xs, ys)
	approx(t, m, 2, 1e-12, "slope")
	approx(t, b, 1, 1e-12, "intercept")
	approx(t, FitSlope(xs, ys), 2, 1e-12, "fitslope")
}

func TestBoxCount(t *testing.T) {
	pts := []Vec{{0, 0}, {0.4, 0.4}, {0.6, 0.6}}
	// eps=0.5: first two share box (0,0); third in box (1,1) -> 2 boxes.
	if n := BoxCount(pts, 0.5); n != 2 {
		t.Errorf("boxcount got %d want 2", n)
	}
}

func TestCorrelationSum(t *testing.T) {
	pts := []Vec{{0}, {1}, {2}, {10}}
	// Pairs with distance < 1.5: (0,1),(1,2) = 2 of 6.
	approx(t, CorrelationSum(pts, 1.5), 2.0/6.0, 1e-12, "corr sum")
}

func TestDimensionsUniform(t *testing.T) {
	rng := rand.New(rand.NewSource(42))
	// Points filling a 2D square -> dimension ~2.
	pts := RandomInitialConditions(rng, 1500, 2, 0, 1)
	d2, _, _ := CorrelationDimension(pts, 0.03, 0.25, 12)
	approx(t, d2, 2, 0.25, "corr dim 2d")
	// Points on a line segment (x, x) -> dimension ~1.
	line := make([]Vec, 1200)
	for i := range line {
		t := rng.Float64()
		line[i] = Vec{t, t}
	}
	d1, _, _ := CorrelationDimension(line, 0.02, 0.2, 12)
	approx(t, d1, 1, 0.2, "corr dim line")
	// Generalized dimension order 0 ~ box-counting on the square.
	dq0, _, _ := BoxCountingDimension(pts, 0.02, 0.2, 10)
	gd := GeneralizedDimension(pts, 0, 0.02, 0.2, 10)
	approx(t, gd, dq0, 1e-9, "generalized q0 == box dim")
}

func TestTakensEmbedding(t *testing.T) {
	series := []float64{1, 2, 3, 4, 5}
	emb := TakensEmbedding(series, 2, 1)
	if len(emb) != 4 {
		t.Fatalf("expected 4 vectors, got %d", len(emb))
	}
	approx(t, emb[0][0], 1, 1e-12, "embed[0][0]")
	approx(t, emb[0][1], 2, 1e-12, "embed[0][1]")
}

// --- util ---

func TestUtil(t *testing.T) {
	ls := Linspace(0, 1, 5)
	approx(t, ls[0], 0, 1e-12, "linspace0")
	approx(t, ls[4], 1, 1e-12, "linspace4")
	approx(t, ls[2], 0.5, 1e-12, "linspace mid")
	lg := Logspace(1, 100, 3)
	approx(t, lg[1], 10, 1e-9, "logspace mid")
	xs := []float64{2, 4, 4, 4, 5, 5, 7, 9}
	approx(t, Mean(xs), 5, 1e-12, "mean")
	approx(t, StdDev(xs), 2, 1e-12, "stddev")
	approx(t, Median([]float64{3, 1, 2}), 2, 1e-12, "median")
	approx(t, Autocorrelation(xs, 0), 1, 1e-12, "autocorr lag0")
	approx(t, WrapUnit(1.25), 0.25, 1e-12, "wrapunit")
	approx(t, WrapAngle(-0.5), 2*math.Pi-0.5, 1e-12, "wrapangle")
	approx(t, LogisticAnalyticLyapunov(), math.Ln2, 1e-12, "logistic analytic")
}

func TestWindingNumber(t *testing.T) {
	// With K=0 the circle map is a pure rotation by omega.
	approx(t, WindingNumber(0.3, 0, 0.1, 100, 1000), 0.3, 1e-9, "winding K=0")
}

func TestInvariantDensityLogistic(t *testing.T) {
	// Logistic r=4 invariant density is 1/(pi sqrt(x(1-x))).
	centers, dens := InvariantDensity(Logistic(4), 0.2, 0, 1, 20, 5000, 400000)
	// Check the density near x=0.5 (minimum of theoretical density = 2/pi).
	for i, c := range centers {
		if c > 0.45 && c < 0.55 {
			theory := 1 / (math.Pi * math.Sqrt(c*(1-c)))
			approx(t, dens[i], theory, 0.15, "invariant density mid")
		}
	}
}

func TestHistogramEntropy(t *testing.T) {
	counts, _ := Histogram([]float64{0.1, 0.2, 0.6, 0.9}, 0, 1, 2)
	if counts[0] != 2 || counts[1] != 2 {
		t.Errorf("histogram got %v", counts)
	}
	approx(t, ShannonEntropy(counts), math.Ln2, 1e-12, "entropy uniform")
}

func TestSpectralRadiusComplex(t *testing.T) {
	R := MatFromRows(VecOf(0, -1), VecOf(1, 0))
	e := Eigenvalues2(R)
	approx(t, cmplx.Abs(e[0]), 1, 1e-12, "complex eig magnitude")
	approx(t, SpectralRadius(R), 1, 1e-12, "spectral radius complex")
}

// --- runnable examples ---

func ExampleLyapunov1D() {
	// The tent map with unit slope has Lyapunov exponent log 2.
	lam := Lyapunov1D(Tent(1.0), 0.1234, 1000, 100000)
	fmt.Printf("%.3f\n", lam)
	// Output: 0.693
}

func ExampleKaplanYorkeDimension() {
	// Kaplan-Yorke dimension from a two-exponent spectrum.
	d := KaplanYorkeDimension(Vec{0.42, -1.62})
	fmt.Printf("%.4f\n", d)
	// Output: 1.2593
}

func ExampleLogisticFixedPoints() {
	fp := LogisticFixedPoints(4)
	fmt.Printf("%.2f %.2f\n", fp[0], fp[1])
	// Output: 0.00 0.75
}
