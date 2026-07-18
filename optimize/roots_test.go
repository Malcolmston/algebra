package optimize

import (
	"math"
	"sort"
	"testing"
)

const testTol = 1e-9

func approx(a, b, tol float64) bool { return math.Abs(a-b) <= tol }

// sqrt2Func is f(x) = x^2 - 2, whose positive root is sqrt(2).
func sqrt2Func(x float64) float64 { return x*x - 2 }

// bracketSolver is the common signature of the always-bracketed solvers.
type bracketSolver func(f Func, a, b, tol float64, maxIter int) (float64, error)

func TestBracketingSolvers(t *testing.T) {
	solvers := map[string]bracketSolver{
		"Bisection":      Bisection,
		"RegulaFalsi":    RegulaFalsi,
		"FalsePosition":  FalsePosition,
		"Illinois":       Illinois,
		"Pegasus":        Pegasus,
		"AndersonBjorck": AndersonBjorck,
		"Brent":          Brent,
		"Ridders":        Ridders,
		"Dekker":         Dekker,
	}
	want := math.Sqrt2
	for name, solve := range solvers {
		got, err := solve(sqrt2Func, 1, 2, 1e-12, 200)
		if err != nil {
			t.Fatalf("%s: unexpected error: %v", name, err)
		}
		if !approx(got, want, 1e-8) {
			t.Errorf("%s: got %v, want %v", name, got, want)
		}
	}
}

func TestBracketingSolversNoBracket(t *testing.T) {
	// f(x) = x^2 + 1 has no real root; [1,2] cannot bracket one.
	f := func(x float64) float64 { return x*x + 1 }
	if _, err := Bisection(f, 1, 2, 1e-12, 100); err != ErrNoBracket {
		t.Errorf("expected ErrNoBracket, got %v", err)
	}
	if _, err := Brent(f, 1, 2, 1e-12, 100); err != ErrNoBracket {
		t.Errorf("expected ErrNoBracket, got %v", err)
	}
}

func TestSecant(t *testing.T) {
	got, err := Secant(sqrt2Func, 1, 2, 1e-12, 100)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(got, math.Sqrt2, 1e-9) {
		t.Errorf("Secant: got %v, want %v", got, math.Sqrt2)
	}
}

func TestNewton(t *testing.T) {
	df := func(x float64) float64 { return 2 * x }
	got, err := Newton(sqrt2Func, df, 2, 1e-12, 100)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(got, math.Sqrt2, 1e-12) {
		t.Errorf("Newton: got %v, want %v", got, math.Sqrt2)
	}
}

func TestNewtonSafe(t *testing.T) {
	df := func(x float64) float64 { return 2 * x }
	got, err := NewtonSafe(sqrt2Func, df, 1, 2, 1e-12, 100)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(got, math.Sqrt2, 1e-10) {
		t.Errorf("NewtonSafe: got %v, want %v", got, math.Sqrt2)
	}
}

func TestHalley(t *testing.T) {
	// f(x) = x^3 - 2, root = cbrt(2).
	f := func(x float64) float64 { return x*x*x - 2 }
	df := func(x float64) float64 { return 3 * x * x }
	d2f := func(x float64) float64 { return 6 * x }
	want := math.Cbrt(2)
	got, err := Halley(f, df, d2f, 1.5, 1e-12, 100)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(got, want, 1e-12) {
		t.Errorf("Halley: got %v, want %v", got, want)
	}
}

func TestSchroeder(t *testing.T) {
	// f(x) = (x-1)^2 has a double root at x = 1.
	f := func(x float64) float64 { return (x - 1) * (x - 1) }
	df := func(x float64) float64 { return 2 * (x - 1) }
	got, err := Schroeder(f, df, 2, 2, 1e-12, 100)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(got, 1, 1e-10) {
		t.Errorf("Schroeder: got %v, want 1", got)
	}
}

func TestSteffensen(t *testing.T) {
	got, err := Steffensen(sqrt2Func, 2, 1e-12, 100)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(got, math.Sqrt2, 1e-9) {
		t.Errorf("Steffensen: got %v, want %v", got, math.Sqrt2)
	}
}

func TestInverseQuadratic(t *testing.T) {
	got, err := InverseQuadratic(sqrt2Func, 1, 1.5, 2, 1e-12, 100)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(got, math.Sqrt2, 1e-9) {
		t.Errorf("InverseQuadratic: got %v, want %v", got, math.Sqrt2)
	}
}

func TestFixedPoint(t *testing.T) {
	// x = cos(x) converges to the Dottie number.
	want := 0.7390851332151607
	got, err := FixedPoint(math.Cos, 1, 1e-12, 500)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(got, want, 1e-8) {
		t.Errorf("FixedPoint: got %v, want %v", got, want)
	}
}

func TestFixedPointAitken(t *testing.T) {
	want := 0.7390851332151607
	got, err := FixedPointAitken(math.Cos, 1, 1e-12, 100)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(got, want, 1e-10) {
		t.Errorf("FixedPointAitken: got %v, want %v", got, want)
	}
}

func TestAitkenDelta2(t *testing.T) {
	// Geometric sequence x_n = 1 - r^n -> 1, r = 0.5.
	x0, x1, x2 := 0.0, 0.5, 0.75
	got := AitkenDelta2(x0, x1, x2)
	if !approx(got, 1, 1e-12) {
		t.Errorf("AitkenDelta2: got %v, want 1", got)
	}
}

func TestPolyEval(t *testing.T) {
	// p(x) = 1 - 2x + x^3, ascending coeffs {1,-2,0,1}.
	coeffs := []float64{1, -2, 0, 1}
	if got := PolyEval(coeffs, 2); !approx(got, 5, 1e-12) {
		t.Errorf("PolyEval: got %v, want 5", got)
	}
	p, dp := PolyEvalDeriv(coeffs, 2)
	if !approx(p, 5, 1e-12) || !approx(dp, 10, 1e-12) {
		t.Errorf("PolyEvalDeriv: got (%v,%v), want (5,10)", p, dp)
	}
}

func TestPolyDerivativeIntegral(t *testing.T) {
	coeffs := []float64{1, -2, 0, 1} // 1 - 2x + x^3
	d := PolyDerivative(coeffs)      // -2 + 3x^2 -> {-2,0,3}
	want := []float64{-2, 0, 3}
	if len(d) != len(want) {
		t.Fatalf("PolyDerivative length: got %v", d)
	}
	for i := range want {
		if !approx(d[i], want[i], 1e-12) {
			t.Errorf("PolyDerivative[%d]: got %v, want %v", i, d[i], want[i])
		}
	}
	// Integral of derivative should recover original up to the constant term.
	in := PolyIntegral(d, coeffs[0])
	for i := range coeffs {
		if !approx(in[i], coeffs[i], 1e-12) {
			t.Errorf("PolyIntegral[%d]: got %v, want %v", i, in[i], coeffs[i])
		}
	}
}

func TestPolyDeflate(t *testing.T) {
	// (x-2)(x-3) = 6 - 5x + x^2, ascending {6,-5,1}. Deflate by root 2.
	q, rem := PolyDeflate([]float64{6, -5, 1}, 2)
	if !approx(rem, 0, 1e-12) {
		t.Errorf("PolyDeflate remainder: got %v, want 0", rem)
	}
	want := []float64{-3, 1} // x - 3
	if len(q) != len(want) {
		t.Fatalf("PolyDeflate quotient length: got %v", q)
	}
	for i := range want {
		if !approx(q[i], want[i], 1e-12) {
			t.Errorf("PolyDeflate quotient[%d]: got %v, want %v", i, q[i], want[i])
		}
	}
}

func TestLinearRoot(t *testing.T) {
	r, ok := LinearRoot(2, -6)
	if !ok || !approx(r, 3, 1e-12) {
		t.Errorf("LinearRoot: got %v, %v; want 3, true", r, ok)
	}
	if _, ok := LinearRoot(0, 1); ok {
		t.Errorf("LinearRoot(0,1): expected ok=false")
	}
}

func TestQuadraticRealRoots(t *testing.T) {
	// x^2 - 3x + 2 -> {1, 2}.
	got := QuadraticRealRoots(1, -3, 2)
	want := []float64{1, 2}
	if len(got) != 2 || !approx(got[0], want[0], 1e-12) || !approx(got[1], want[1], 1e-12) {
		t.Errorf("QuadraticRealRoots: got %v, want %v", got, want)
	}
	// Complex roots: x^2 + 1 -> none real.
	if r := QuadraticRealRoots(1, 0, 1); len(r) != 0 {
		t.Errorf("QuadraticRealRoots(x^2+1): got %v, want empty", r)
	}
}

func TestQuadraticRoots(t *testing.T) {
	// x^2 + 1 -> ±i.
	r1, r2 := QuadraticRoots(1, 0, 1)
	if !approx(real(r1), 0, 1e-12) || !approx(math.Abs(imag(r1)), 1, 1e-12) {
		t.Errorf("QuadraticRoots imag: got %v, %v", r1, r2)
	}
}

func TestCubicRealRoots(t *testing.T) {
	// x^3 - 6x^2 + 11x - 6 -> {1, 2, 3}.
	got := CubicRealRoots(1, -6, 11, -6)
	want := []float64{1, 2, 3}
	if len(got) != 3 {
		t.Fatalf("CubicRealRoots: got %v, want %v", got, want)
	}
	for i := range want {
		if !approx(got[i], want[i], 1e-6) {
			t.Errorf("CubicRealRoots[%d]: got %v, want %v", i, got[i], want[i])
		}
	}
}

func TestCubicOneRealRoot(t *testing.T) {
	// x^3 + x + 1 has a single real root near -0.6823278.
	got := CubicRealRoots(1, 0, 1, 1)
	if len(got) != 1 {
		t.Fatalf("CubicRealRoots: got %v, want 1 real root", got)
	}
	if !approx(got[0], -0.6823278038280193, 1e-9) {
		t.Errorf("CubicRealRoots: got %v", got[0])
	}
}

func realParts(roots []complex128, tol float64) []float64 {
	var out []float64
	for _, r := range roots {
		if math.Abs(imag(r)) <= tol {
			out = append(out, real(r))
		}
	}
	sort.Float64s(out)
	return out
}

func TestQuarticRoots(t *testing.T) {
	// (x-1)(x-2)(x-3)(x-4) = x^4 -10x^3 +35x^2 -50x +24.
	roots := QuarticRoots(1, -10, 35, -50, 24)
	got := realParts(roots, 1e-6)
	want := []float64{1, 2, 3, 4}
	if len(got) != 4 {
		t.Fatalf("QuarticRoots: got %v real parts, want %v", got, want)
	}
	for i := range want {
		if !approx(got[i], want[i], 1e-6) {
			t.Errorf("QuarticRoots[%d]: got %v, want %v", i, got[i], want[i])
		}
	}
}

func TestDurandKernerAberth(t *testing.T) {
	// x^2 - 2 -> ±sqrt(2), ascending coeffs {-2,0,1}.
	coeffs := []float64{-2, 0, 1}
	for _, tc := range []struct {
		name  string
		solve func([]float64, float64, int) []complex128
	}{
		{"DurandKerner", DurandKerner},
		{"Aberth", Aberth},
	} {
		got := realParts(tc.solve(coeffs, 1e-14, 500), 1e-8)
		if len(got) != 2 || !approx(got[0], -math.Sqrt2, 1e-8) || !approx(got[1], math.Sqrt2, 1e-8) {
			t.Errorf("%s: got %v, want ±sqrt2", tc.name, got)
		}
	}
}

func TestPolyRoots(t *testing.T) {
	// (x+1)(x-2)(x-4) = x^3 -5x^2 +2x +8, ascending {8,2,-5,1}.
	got := PolyRoots([]float64{8, 2, -5, 1}, 1e-7)
	want := []float64{-1, 2, 4}
	if len(got) != 3 {
		t.Fatalf("PolyRoots: got %v, want %v", got, want)
	}
	for i := range want {
		if !approx(got[i], want[i], 1e-7) {
			t.Errorf("PolyRoots[%d]: got %v, want %v", i, got[i], want[i])
		}
	}
}

func TestCompanionMatrix(t *testing.T) {
	// x^2 - 3x + 2, ascending {2,-3,1}. Companion last column {-2, 3}.
	m := CompanionMatrix([]float64{2, -3, 1})
	if len(m) != 2 || len(m[0]) != 2 {
		t.Fatalf("CompanionMatrix: wrong shape %v", m)
	}
	if !approx(m[0][1], -2, 1e-12) || !approx(m[1][1], 3, 1e-12) || !approx(m[1][0], 1, 1e-12) {
		t.Errorf("CompanionMatrix: got %v", m)
	}
}

func TestMuller(t *testing.T) {
	// f(z) = z^2 + 1 -> ±i.
	f := func(z complex128) complex128 { return z*z + 1 }
	got, err := Muller(f, complex(1, 0), complex(0.5, 0.5), complex(0, 0.5), 1e-12, 100)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(math.Abs(imag(got)), 1, 1e-8) || !approx(real(got), 0, 1e-8) {
		t.Errorf("Muller: got %v, want ±i", got)
	}
}

func TestNewtonComplex(t *testing.T) {
	// z^2 + 1 -> i from a complex start.
	f := func(z complex128) complex128 { return z*z + 1 }
	df := func(z complex128) complex128 { return 2 * z }
	got, err := NewtonComplex(f, df, complex(0.1, 1), 1e-12, 100)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(real(got), 0, 1e-9) || !approx(imag(got), 1, 1e-9) {
		t.Errorf("NewtonComplex: got %v, want i", got)
	}
}

func TestGoldenSectionMin(t *testing.T) {
	// (x-2)^2 has its minimum at x = 2.
	f := func(x float64) float64 { return (x - 2) * (x - 2) }
	got, err := GoldenSectionMin(f, 0, 5, 1e-10, 200)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(got, 2, 1e-6) {
		t.Errorf("GoldenSectionMin: got %v, want 2", got)
	}
}

func TestBrentMinimize(t *testing.T) {
	f := func(x float64) float64 { return (x - 2) * (x - 2) }
	xmin, fmin, err := BrentMinimize(f, 0, 5, 1e-10, 200)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(xmin, 2, 1e-6) || !approx(fmin, 0, 1e-10) {
		t.Errorf("BrentMinimize: got (%v,%v), want (2,0)", xmin, fmin)
	}
}

func TestTernarySearch(t *testing.T) {
	f := func(x float64) float64 { return (x - 2) * (x - 2) }
	got, err := TernarySearch(f, -1, 6, 1e-9, 500)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(got, 2, 1e-6) {
		t.Errorf("TernarySearch: got %v, want 2", got)
	}
}

func TestParabolicMinimum(t *testing.T) {
	// (x-2)^2 sampled at x = 1, 2.5, 4.
	f := func(x float64) float64 { return (x - 2) * (x - 2) }
	got, ok := ParabolicMinimum(1, 2.5, 4, f(1), f(2.5), f(4))
	if !ok || !approx(got, 2, 1e-9) {
		t.Errorf("ParabolicMinimum: got %v, %v; want 2, true", got, ok)
	}
}

func TestBracketMinimum(t *testing.T) {
	f := func(x float64) float64 { return (x - 2) * (x - 2) }
	a, b, c := BracketMinimum(f, 0, 1)
	if !(f(b) < f(a) && f(b) < f(c)) {
		t.Errorf("BracketMinimum: f(b) not the lowest: a=%v b=%v c=%v", a, b, c)
	}
	// Confirmed minimum lies inside the bracket.
	lo, hi := a, c
	if lo > hi {
		lo, hi = hi, lo
	}
	if !(2 >= lo && 2 <= hi) {
		t.Errorf("BracketMinimum: minimum 2 not in [%v,%v]", lo, hi)
	}
}

func TestNewtonMinimize(t *testing.T) {
	// Minimize (x-2)^2: df = 2(x-2), d2f = 2.
	df := func(x float64) float64 { return 2 * (x - 2) }
	d2f := func(x float64) float64 { return 2 }
	got, err := NewtonMinimize(df, d2f, 0, 1e-12, 100)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(got, 2, 1e-10) {
		t.Errorf("NewtonMinimize: got %v, want 2", got)
	}
}

func TestGradientDescent(t *testing.T) {
	df := func(x float64) float64 { return 2 * (x - 2) }
	got, err := GradientDescent(df, 0, 0.4, 1e-12, 1000)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(got, 2, 1e-6) {
		t.Errorf("GradientDescent: got %v, want 2", got)
	}
}

func TestDerivatives(t *testing.T) {
	// d/dx sin(x) = cos(x); at x = 1, cos(1).
	x := 1.0
	wantD := math.Cos(1)
	if got := Derivative(math.Sin, x, 1e-5); !approx(got, wantD, 1e-8) {
		t.Errorf("Derivative: got %v, want %v", got, wantD)
	}
	if got := RichardsonDerivative(math.Sin, x, 1e-2); !approx(got, wantD, 1e-8) {
		t.Errorf("RichardsonDerivative: got %v, want %v", got, wantD)
	}
	if got := ForwardDifference(math.Sin, x, 1e-6); !approx(got, wantD, 1e-4) {
		t.Errorf("ForwardDifference: got %v, want %v", got, wantD)
	}
	if got := BackwardDifference(math.Sin, x, 1e-6); !approx(got, wantD, 1e-4) {
		t.Errorf("BackwardDifference: got %v, want %v", got, wantD)
	}
	// d2/dx2 sin(x) = -sin(x).
	wantD2 := -math.Sin(1)
	if got := SecondDerivative(math.Sin, x, 1e-4); !approx(got, wantD2, 1e-4) {
		t.Errorf("SecondDerivative: got %v, want %v", got, wantD2)
	}
}

func TestBracketingHelpers(t *testing.T) {
	if !SignChange(sqrt2Func, 1, 2) {
		t.Error("SignChange(1,2) should be true")
	}
	if SignChange(sqrt2Func, 1.5, 1.6) {
		t.Error("SignChange(1.5,1.6) should be false")
	}
	if !SameSign(3, 5) || SameSign(3, -5) {
		t.Error("SameSign wrong")
	}
	if !OppositeSign(3, -5) || OppositeSign(3, 5) {
		t.Error("OppositeSign wrong")
	}
	if !IsBracket(sqrt2Func, 0, 2) {
		t.Error("IsBracket(0,2) should be true")
	}
}

func TestBracketExpandAndSubdivide(t *testing.T) {
	lo, hi, ok := BracketExpand(sqrt2Func, 1.3, 1.4, 1.6, 60)
	if !ok || !SignChange(sqrt2Func, lo, hi) {
		t.Errorf("BracketExpand failed: lo=%v hi=%v ok=%v", lo, hi, ok)
	}
	// x^3 - x = x(x-1)(x+1) has roots at -1, 0, 1.
	f := func(x float64) float64 { return x*x*x - x }
	brs := BracketSubdivide(f, -2, 2, 40)
	if len(brs) < 3 {
		t.Fatalf("BracketSubdivide: got %d brackets, want >=3", len(brs))
	}
	// Each bracket should refine to a known root.
	roots := map[float64]bool{}
	for _, br := range brs {
		r, err := Brent(f, br.Lo, br.Hi, 1e-12, 100)
		if err != nil {
			t.Fatal(err)
		}
		roots[math.Round(r)] = true
	}
	for _, want := range []float64{-1, 0, 1} {
		if !roots[want] {
			t.Errorf("BracketSubdivide/Brent missed root %v (found %v)", want, roots)
		}
	}
}

func TestFindBracket(t *testing.T) {
	br, ok := FindBracket(sqrt2Func, 1.3, 1.35)
	if !ok || !br.Contains(math.Sqrt2) {
		t.Errorf("FindBracket: br=%v ok=%v", br, ok)
	}
}

func TestBracketType(t *testing.T) {
	br := Bracket{Lo: 1, Hi: 4}
	if !approx(br.Width(), 3, 1e-12) {
		t.Errorf("Width: got %v", br.Width())
	}
	if !approx(br.Midpoint(), 2.5, 1e-12) {
		t.Errorf("Midpoint: got %v", br.Midpoint())
	}
	if !br.Contains(2) || br.Contains(5) {
		t.Errorf("Contains wrong")
	}
}

func TestUtilities(t *testing.T) {
	if Sign(-3) != -1 || Sign(3) != 1 || Sign(0) != 0 {
		t.Error("Sign wrong")
	}
	if Clamp(5, 0, 3) != 3 || Clamp(-1, 0, 3) != 0 || Clamp(2, 0, 3) != 2 {
		t.Error("Clamp wrong")
	}
	if !approx(RelativeError(1.01, 1), 0.01, 1e-12) {
		t.Error("RelativeError wrong")
	}
	if !approx(AbsoluteError(2, 5), 3, 1e-12) {
		t.Error("AbsoluteError wrong")
	}
	if !WithinTolerance(1, 1.0001, 1e-3) || WithinTolerance(1, 2, 1e-3) {
		t.Error("WithinTolerance wrong")
	}
	if !Converged(1.0, 1.0+1e-12, 1e-9) {
		t.Error("Converged wrong")
	}
}

func TestErrorConstantsAndDefaults(t *testing.T) {
	if DefaultTolerance <= 0 || DefaultMaxIterations <= 0 {
		t.Error("defaults should be positive")
	}
	if ErrInvalidInterval == nil {
		t.Error("ErrInvalidInterval should be defined")
	}
}

// BenchmarkDurandKerner exercises the heaviest routine: simultaneous
// root finding on a degree-four polynomial.
func BenchmarkDurandKerner(b *testing.B) {
	coeffs := []float64{24, -50, 35, -10, 1} // (x-1)(x-2)(x-3)(x-4)
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		_ = DurandKerner(coeffs, 1e-14, 500)
	}
}
