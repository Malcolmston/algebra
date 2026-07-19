package splines

import (
	"fmt"
	"math"
	"testing"
)

const tol = 1e-9

func close(a, b float64) bool { return math.Abs(a-b) <= tol }

func closeT(a, b, t float64) bool { return math.Abs(a-b) <= t }

func vecClose(a, b Vec, t float64) bool { return a.Equal(b, t) }

// --- Vec ----------------------------------------------------------------

func TestVecArithmetic(t *testing.T) {
	a := NewVec(1, 2, 3)
	b := NewVec(4, 5, 6)
	if got := a.Add(b); !vecClose(got, NewVec(5, 7, 9), tol) {
		t.Errorf("Add=%v", got)
	}
	if got := b.Sub(a); !vecClose(got, NewVec(3, 3, 3), tol) {
		t.Errorf("Sub=%v", got)
	}
	if got := a.Scale(2); !vecClose(got, NewVec(2, 4, 6), tol) {
		t.Errorf("Scale=%v", got)
	}
	if got := a.Dot(b); !close(got, 32) {
		t.Errorf("Dot=%v", got)
	}
	if got := NewVec(3, 4).Norm(); !close(got, 5) {
		t.Errorf("Norm=%v", got)
	}
	if got := a.AddScaled(2, b); !vecClose(got, NewVec(9, 12, 15), tol) {
		t.Errorf("AddScaled=%v", got)
	}
	if got := NewVec(0, 0).Lerp(NewVec(10, 20), 0.25); !vecClose(got, NewVec(2.5, 5), tol) {
		t.Errorf("Lerp=%v", got)
	}
	if got := Cross(NewVec(1, 0, 0), NewVec(0, 1, 0)); !vecClose(got, NewVec(0, 0, 1), tol) {
		t.Errorf("Cross=%v", got)
	}
}

// --- Tridiagonal solvers ------------------------------------------------

func TestSolveTridiagonal(t *testing.T) {
	// System:
	// [2 1 0][x]   [3]
	// [1 2 1][ ] = [4]
	// [0 1 2][ ]   [3]
	// Solution x = (1,1,1).
	a := []float64{0, 1, 1}
	b := []float64{2, 2, 2}
	c := []float64{1, 1, 0}
	d := []float64{3, 4, 3}
	x, err := SolveTridiagonal(a, b, c, d)
	if err != nil {
		t.Fatal(err)
	}
	for i, v := range x {
		if !close(v, 1) {
			t.Errorf("x[%d]=%v want 1", i, v)
		}
	}
}

func TestSolveCyclicTridiagonal(t *testing.T) {
	// Circulant-like symmetric system with corners; verify by residual.
	n := 5
	a := make([]float64, n)
	b := make([]float64, n)
	c := make([]float64, n)
	d := make([]float64, n)
	for i := 0; i < n; i++ {
		a[i], b[i], c[i], d[i] = 1, 4, 1, float64(i+1)
	}
	alpha, beta := 1.0, 1.0
	x, err := SolveCyclicTridiagonal(a, b, c, d, alpha, beta)
	if err != nil {
		t.Fatal(err)
	}
	// Residual check for the full cyclic matrix.
	for i := 0; i < n; i++ {
		lo := (i - 1 + n) % n
		hi := (i + 1) % n
		res := b[i]*x[i] + a[i]*x[lo]*btoi(i != 0) + c[i]*x[hi]*btoi(i != n-1)
		if i == 0 {
			res += alpha * x[n-1]
		}
		if i == n-1 {
			res += beta * x[0]
		}
		if !closeT(res, d[i], 1e-9) {
			t.Errorf("row %d residual=%v want %v", i, res, d[i])
		}
	}
}

func btoi(b bool) float64 {
	if b {
		return 1
	}
	return 0
}

// --- Cubic splines ------------------------------------------------------

func refCubic(x float64) float64  { return 2 - x + 0.5*x*x + 0.3*x*x*x }
func refCubicD(x float64) float64 { return -1 + x + 0.9*x*x }

func TestNotAKnotReproducesCubic(t *testing.T) {
	x := []float64{0, 1, 2, 3, 4, 5}
	y := make([]float64, len(x))
	for i := range x {
		y[i] = refCubic(x[i])
	}
	cs, err := NotAKnotCubicSpline(x, y)
	if err != nil {
		t.Fatal(err)
	}
	for _, q := range []float64{0.3, 1.5, 2.7, 3.9, 4.8} {
		if !closeT(cs.Eval(q), refCubic(q), 1e-8) {
			t.Errorf("Eval(%v)=%v want %v", q, cs.Eval(q), refCubic(q))
		}
		if !closeT(cs.EvalDerivative(q), refCubicD(q), 1e-8) {
			t.Errorf("Eval'(%v)=%v want %v", q, cs.EvalDerivative(q), refCubicD(q))
		}
	}
}

func TestClampedReproducesCubic(t *testing.T) {
	x := []float64{0, 1, 2, 3, 4, 5}
	y := make([]float64, len(x))
	for i := range x {
		y[i] = refCubic(x[i])
	}
	cs, err := ClampedCubicSpline(x, y, refCubicD(0), refCubicD(5))
	if err != nil {
		t.Fatal(err)
	}
	for _, q := range []float64{0.5, 2.2, 4.6} {
		if !closeT(cs.Eval(q), refCubic(q), 1e-8) {
			t.Errorf("Eval(%v)=%v want %v", q, cs.Eval(q), refCubic(q))
		}
	}
	// Endpoint derivatives must match the clamp.
	if !closeT(cs.EvalDerivative(0), refCubicD(0), 1e-9) {
		t.Errorf("left slope=%v", cs.EvalDerivative(0))
	}
	if !closeT(cs.EvalDerivative(5), refCubicD(5), 1e-9) {
		t.Errorf("right slope=%v", cs.EvalDerivative(5))
	}
}

func TestNaturalSecondDerivZeroAndInterpolates(t *testing.T) {
	x := []float64{0, 1, 2, 3, 4}
	y := []float64{0, 1, 0, 1, 0}
	cs, err := NaturalCubicSpline(x, y)
	if err != nil {
		t.Fatal(err)
	}
	if !close(cs.EvalSecondDerivative(0), 0) || !close(cs.EvalSecondDerivative(4), 0) {
		t.Errorf("natural end 2nd derivs not zero: %v %v",
			cs.EvalSecondDerivative(0), cs.EvalSecondDerivative(4))
	}
	for i := range x {
		if !close(cs.Eval(x[i]), y[i]) {
			t.Errorf("interp fail at node %d: %v want %v", i, cs.Eval(x[i]), y[i])
		}
	}
}

func TestPeriodicContinuity(t *testing.T) {
	x := []float64{0, 1, 2, 3, 4}
	y := []float64{0, 1, 0, -1, 0}
	cs, err := PeriodicCubicSpline(x, y)
	if err != nil {
		t.Fatal(err)
	}
	// First and second derivatives at the two ends must agree.
	if !closeT(cs.EvalDerivative(0), cs.EvalDerivative(4), 1e-8) {
		t.Errorf("periodic slope mismatch: %v vs %v", cs.EvalDerivative(0), cs.EvalDerivative(4))
	}
	if !closeT(cs.EvalSecondDerivative(0), cs.EvalSecondDerivative(4), 1e-8) {
		t.Errorf("periodic curvature mismatch: %v vs %v",
			cs.EvalSecondDerivative(0), cs.EvalSecondDerivative(4))
	}
}

func TestCubicIntegrate(t *testing.T) {
	// Integrate the reproduced cubic and compare with the analytic integral.
	x := []float64{0, 1, 2, 3, 4, 5}
	y := make([]float64, len(x))
	for i := range x {
		y[i] = refCubic(x[i])
	}
	cs, _ := NotAKnotCubicSpline(x, y)
	// analytic: int_0^5 (2 - x + 0.5x^2 + 0.3x^3) dx
	F := func(x float64) float64 { return 2*x - x*x/2 + 0.5*x*x*x/3 + 0.3*x*x*x*x/4 }
	want := F(5) - F(0)
	if got := cs.Integrate(0, 5); !closeT(got, want, 1e-7) {
		t.Errorf("Integrate=%v want %v", got, want)
	}
	if got := cs.Integrate(4, 4); !close(got, 0) {
		t.Errorf("empty integral=%v", got)
	}
	if got := cs.Integrate(5, 0); !closeT(got, -want, 1e-7) {
		t.Errorf("reversed integral=%v", got)
	}
}

// --- Hermite / Catmull-Rom / monotone / Akima ---------------------------

func TestHermiteInterpolates(t *testing.T) {
	x := []float64{0, 1, 3}
	y := []float64{0, 2, 1}
	m := []float64{1, 0, -1}
	h, err := HermiteSpline(x, y, m)
	if err != nil {
		t.Fatal(err)
	}
	for i := range x {
		if !close(h.Eval(x[i]), y[i]) {
			t.Errorf("value at node %d", i)
		}
		if !close(h.EvalDerivative(x[i]), m[i]) {
			t.Errorf("tangent at node %d = %v want %v", i, h.EvalDerivative(x[i]), m[i])
		}
	}
}

func TestHermiteBasisPartition(t *testing.T) {
	for _, tt := range []float64{0, 0.25, 0.5, 0.75, 1} {
		s := HermiteBasis00(tt) + HermiteBasis01(tt)
		if !close(s, 1) {
			t.Errorf("h00+h01 at %v = %v", tt, s)
		}
	}
	if !close(HermiteBasis00(0), 1) || !close(HermiteBasis01(1), 1) {
		t.Errorf("endpoint interpolation basis wrong")
	}
}

func TestCatmullRomInterpolates(t *testing.T) {
	x := []float64{0, 1, 2, 3}
	y := []float64{1, 3, 2, 5}
	cr, err := CatmullRomSpline(x, y)
	if err != nil {
		t.Fatal(err)
	}
	for i := range x {
		if !close(cr.Eval(x[i]), y[i]) {
			t.Errorf("CR node %d = %v want %v", i, cr.Eval(x[i]), y[i])
		}
	}
}

func TestMonotonePreservesMonotonicity(t *testing.T) {
	x := []float64{0, 1, 2, 3, 4, 5}
	y := []float64{0, 0, 0, 1, 1, 1} // monotone non-decreasing with a flat then a step
	m, err := MonotoneCubicSpline(x, y)
	if err != nil {
		t.Fatal(err)
	}
	prev := m.Eval(0)
	for q := 0.0; q <= 5.0; q += 0.05 {
		v := m.Eval(q)
		if v < prev-1e-12 {
			t.Errorf("monotonicity violated at %v: %v < %v", q, v, prev)
		}
		if v < -1e-12 || v > 1+1e-12 {
			t.Errorf("overshoot at %v: %v", q, v)
		}
		prev = v
	}
	for i := range x {
		if !close(m.Eval(x[i]), y[i]) {
			t.Errorf("monotone interp at node %d", i)
		}
	}
}

func TestAkimaInterpolatesAndValue(t *testing.T) {
	x := []float64{0, 1, 2, 3, 4}
	y := []float64{0, 0, 0, 1, 1}
	ak, err := AkimaSpline(x, y)
	if err != nil {
		t.Fatal(err)
	}
	for i := range x {
		if !close(ak.Eval(x[i]), y[i]) {
			t.Errorf("akima interp at node %d", i)
		}
	}
	if !closeT(ak.Eval(2.5), 0.3125, 1e-9) {
		t.Errorf("akima Eval(2.5)=%v want 0.3125", ak.Eval(2.5))
	}
}

// --- Bernstein / Bezier -------------------------------------------------

func TestBinomial(t *testing.T) {
	cases := []struct {
		n, k int
		want float64
	}{
		{0, 0, 1}, {5, 0, 1}, {5, 5, 1}, {5, 2, 10}, {6, 3, 20}, {10, 4, 210}, {5, 6, 0},
	}
	for _, c := range cases {
		if got := Binomial(c.n, c.k); got != c.want {
			t.Errorf("Binomial(%d,%d)=%v want %v", c.n, c.k, got, c.want)
		}
	}
}

func TestBernsteinPartitionOfUnity(t *testing.T) {
	for _, u := range []float64{0, 0.2, 0.5, 0.9, 1} {
		all := BernsteinAll(4, u)
		var s float64
		for i, v := range all {
			s += v
			if !close(v, Bernstein(4, i, u)) {
				t.Errorf("BernsteinAll mismatch at i=%d u=%v", i, u)
			}
		}
		if !close(s, 1) {
			t.Errorf("Bernstein sum at %v = %v", u, s)
		}
	}
}

func TestBezierEvalDerivativeElevateSplit(t *testing.T) {
	b, err := NewBezierCurve(NewVec(0, 0), NewVec(1, 2), NewVec(3, 2), NewVec(4, 0))
	if err != nil {
		t.Fatal(err)
	}
	if !vecClose(b.Eval(0), NewVec(0, 0), tol) || !vecClose(b.Eval(1), NewVec(4, 0), tol) {
		t.Errorf("bezier endpoints wrong")
	}
	if !vecClose(b.Eval(0.5), NewVec(2, 1.5), tol) {
		t.Errorf("bezier mid=%v", b.Eval(0.5))
	}
	if !vecClose(b.EvalDerivative(0.5), NewVec(4.5, 0), tol) {
		t.Errorf("bezier deriv=%v", b.EvalDerivative(0.5))
	}
	// Degree elevation preserves the curve.
	e := b.Elevate()
	if e.Degree() != 4 {
		t.Errorf("elevated degree=%d", e.Degree())
	}
	for _, u := range []float64{0, 0.3, 0.6, 1} {
		if !vecClose(e.Eval(u), b.Eval(u), 1e-12) {
			t.Errorf("elevation changed curve at %v", u)
		}
	}
	// Subdivision reproduces the curve.
	l, r := b.Split(0.4)
	for _, u := range []float64{0, 0.5, 1} {
		if !vecClose(l.Eval(u), b.Eval(0.4*u), 1e-12) {
			t.Errorf("split-left mismatch at %v", u)
		}
		if !vecClose(r.Eval(u), b.Eval(0.4+0.6*u), 1e-12) {
			t.Errorf("split-right mismatch at %v", u)
		}
	}
	// Reverse.
	rev := b.Reverse()
	for _, u := range []float64{0, 0.3, 1} {
		if !vecClose(rev.Eval(u), b.Eval(1-u), 1e-12) {
			t.Errorf("reverse mismatch at %v", u)
		}
	}
	// DeCasteljau free function agrees with method.
	if !vecClose(DeCasteljau(b.ControlPoints(), 0.7), b.Eval(0.7), 1e-12) {
		t.Errorf("DeCasteljau disagreement")
	}
}

func TestBezierLineLength(t *testing.T) {
	b, _ := NewBezierCurve(NewVec(0, 0), NewVec(1, 1), NewVec(2, 2), NewVec(3, 3))
	want := math.Sqrt(18)
	if got := b.Length(); !closeT(got, want, 1e-7) {
		t.Errorf("length=%v want %v", got, want)
	}
}

// --- B-spline -----------------------------------------------------------

func TestBSplineBasisPartitionOfUnity(t *testing.T) {
	U := []float64{0, 0, 0, 1, 2, 3, 3, 3}
	p := 2
	n := len(U) - p - 2 // last control index
	for _, u := range []float64{0, 0.5, 1, 1.7, 2.5, 3} {
		span := FindSpan(n, p, u, U)
		N := BasisFuns(span, u, p, U)
		var s float64
		for _, v := range N {
			s += v
		}
		if !close(s, 1) {
			t.Errorf("basis sum at %v = %v", u, s)
		}
		// OneBasisFun must agree with the corresponding non-zero basis value.
		for j := 0; j <= p; j++ {
			i := span - p + j
			if !closeT(OneBasisFun(p, U, i, u), N[j], 1e-12) {
				t.Errorf("OneBasisFun mismatch u=%v i=%d", u, i)
			}
		}
	}
}

func TestBSplineEndpointsAndDeBoor(t *testing.T) {
	ctrl := []Vec{NewVec(0, 0), NewVec(1, 2), NewVec(3, 2), NewVec(4, 0), NewVec(5, 2)}
	U := ClampedKnots(len(ctrl), 3)
	bs, err := NewBSplineCurve(ctrl, U, 3)
	if err != nil {
		t.Fatal(err)
	}
	if !vecClose(bs.Eval(0), ctrl[0], tol) || !vecClose(bs.Eval(1), ctrl[len(ctrl)-1], tol) {
		t.Errorf("clamped B-spline does not interpolate ends")
	}
	for _, u := range []float64{0, 0.25, 0.5, 0.75, 1} {
		if !vecClose(DeBoor(3, U, ctrl, u), bs.Eval(u), 1e-12) {
			t.Errorf("DeBoor vs Eval mismatch at %v", u)
		}
	}
}

func TestBSplineReducesToBezier(t *testing.T) {
	ctrl := []Vec{NewVec(0, 0), NewVec(1, 2), NewVec(3, 2), NewVec(4, 0)}
	U := []float64{0, 0, 0, 0, 1, 1, 1, 1}
	bs, _ := NewBSplineCurve(ctrl, U, 3)
	bz, _ := NewBezierCurve(ctrl...)
	for _, u := range []float64{0, 0.2, 0.5, 0.8, 1} {
		if !vecClose(bs.Eval(u), bz.Eval(u), 1e-12) {
			t.Errorf("B-spline != Bezier at %v", u)
		}
	}
}

func TestKnotInsertionInvariance(t *testing.T) {
	ctrl := []Vec{NewVec(0, 0), NewVec(1, 2), NewVec(3, 2), NewVec(4, 0), NewVec(5, 2)}
	U := ClampedKnots(len(ctrl), 3)
	bs, _ := NewBSplineCurve(ctrl, U, 3)
	bs2, err := bs.InsertKnot(0.5, 1)
	if err != nil {
		t.Fatal(err)
	}
	if len(bs2.ControlPoints()) != len(ctrl)+1 {
		t.Errorf("knot insertion did not add a control point")
	}
	for _, u := range []float64{0, 0.25, 0.5, 0.75, 1} {
		if !vecClose(bs2.Eval(u), bs.Eval(u), 1e-11) {
			t.Errorf("knot insertion changed curve at %v", u)
		}
	}
}

func TestBSplineDerivativeFiniteDiff(t *testing.T) {
	ctrl := []Vec{NewVec(0, 0), NewVec(1, 2), NewVec(3, 2), NewVec(4, 0), NewVec(5, 2)}
	U := ClampedKnots(len(ctrl), 3)
	bs, _ := NewBSplineCurve(ctrl, U, 3)
	for _, u := range []float64{0.2, 0.5, 0.8} {
		h := 1e-6
		fd := bs.Eval(u + h).Sub(bs.Eval(u - h)).Scale(1 / (2 * h))
		if !vecClose(bs.EvalDerivative(u), fd, 1e-4) {
			t.Errorf("deriv mismatch at %v: %v vs %v", u, bs.EvalDerivative(u), fd)
		}
	}
}

// --- NURBS --------------------------------------------------------------

func TestNURBSCircleExact(t *testing.T) {
	c, err := NURBSCircle(NewVec(0, 0), 2)
	if err != nil {
		t.Fatal(err)
	}
	for _, u := range []float64{0, 0.1, 0.125, 0.25, 0.4, 0.5, 0.7, 0.875, 1} {
		p := c.Eval(u)
		if !closeT(math.Hypot(p[0], p[1]), 2, 1e-12) {
			t.Errorf("circle radius at u=%v = %v", u, math.Hypot(p[0], p[1]))
		}
	}
}

func TestNURBSEqualsWeightedBSpline(t *testing.T) {
	ctrl := []Vec{NewVec(0, 0), NewVec(1, 2), NewVec(3, 2), NewVec(4, 0)}
	U := []float64{0, 0, 0, 0, 1, 1, 1, 1}
	w := []float64{1, 1, 1, 1}
	nb, err := NewNURBSCurve(ctrl, w, U, 3)
	if err != nil {
		t.Fatal(err)
	}
	bs, _ := NewBSplineCurve(ctrl, U, 3)
	for _, u := range []float64{0, 0.3, 0.6, 1} {
		if !vecClose(nb.Eval(u), bs.Eval(u), 1e-12) {
			t.Errorf("NURBS(w=1) != B-spline at %v", u)
		}
	}
	// Derivative sanity via finite difference.
	u := 0.5
	h := 1e-6
	fd := nb.Eval(u + h).Sub(nb.Eval(u - h)).Scale(1 / (2 * h))
	if !vecClose(nb.EvalDerivative(u), fd, 1e-4) {
		t.Errorf("NURBS derivative mismatch: %v vs %v", nb.EvalDerivative(u), fd)
	}
}

// --- Arc length ---------------------------------------------------------

func TestGaussLegendreIntegrate(t *testing.T) {
	// int_0^1 x^2 = 1/3, int_0^pi sin = 2.
	if got := GaussLegendreIntegrate(func(x float64) float64 { return x * x }, 0, 1, 5); !close(got, 1.0/3) {
		t.Errorf("int x^2 = %v", got)
	}
	if got := GaussLegendreIntegrate(math.Sin, 0, math.Pi, 8); !closeT(got, 2, 1e-9) {
		t.Errorf("int sin = %v", got)
	}
}

func TestArcLengthParamUnitSpeed(t *testing.T) {
	b, _ := NewBezierCurve(NewVec(0, 0), NewVec(1, 2), NewVec(3, 2), NewVec(4, 0))
	ap, err := NewArcLengthParam(b, 200)
	if err != nil {
		t.Fatal(err)
	}
	total := ap.TotalLength()
	// Points at equal arc-length fractions must be equally spaced in arc length.
	var prev Vec
	step := total / 10
	for k := 0; k <= 10; k++ {
		p := ap.PointAtArcLength(float64(k) * step)
		if k > 0 {
			// The chord underestimates arc length but should be close for small steps.
			d := p.Dist(prev)
			if d > step+1e-6 {
				t.Errorf("chord %v exceeds arc step %v", d, step)
			}
		}
		prev = p
	}
	// Round trip: TForS(SForT(t)) == t.
	for _, tt := range []float64{0.1, 0.37, 0.62, 0.9} {
		s := ap.SForT(tt)
		back := ap.TForS(s)
		if !closeT(back, tt, 1e-4) {
			t.Errorf("arc-length round trip: %v -> %v -> %v", tt, s, back)
		}
	}
}

func TestCurveLengthLine(t *testing.T) {
	b, _ := NewBezierCurve(NewVec(0, 0), NewVec(3, 4))
	if got := CurveLength(b); !closeT(got, 5, 1e-9) {
		t.Errorf("line length=%v want 5", got)
	}
}

// --- Surfaces -----------------------------------------------------------

func bezierTestNet() [][]Vec {
	return [][]Vec{
		{NewVec(0, 0, 0), NewVec(1, 0, 1), NewVec(2, 0, 0)},
		{NewVec(0, 1, 1), NewVec(1, 1, 2), NewVec(2, 1, 1)},
		{NewVec(0, 2, 0), NewVec(1, 2, 1), NewVec(2, 2, 0)},
	}
}

func TestBezierSurface(t *testing.T) {
	s, err := NewBezierSurface(bezierTestNet())
	if err != nil {
		t.Fatal(err)
	}
	// Corners interpolate the corner control points.
	if !vecClose(s.Eval(0, 0), NewVec(0, 0, 0), tol) {
		t.Errorf("corner (0,0)=%v", s.Eval(0, 0))
	}
	if !vecClose(s.Eval(1, 1), NewVec(2, 2, 0), tol) {
		t.Errorf("corner (1,1)=%v", s.Eval(1, 1))
	}
	if !vecClose(s.Eval(0.5, 0.5), NewVec(1, 1, 1), tol) {
		t.Errorf("center=%v", s.Eval(0.5, 0.5))
	}
	// Partial derivative via finite difference.
	h := 1e-6
	fd := s.Eval(0.5+h, 0.5).Sub(s.Eval(0.5-h, 0.5)).Scale(1 / (2 * h))
	if !vecClose(s.EvalPartialU(0.5, 0.5), fd, 1e-4) {
		t.Errorf("partial-u mismatch: %v vs %v", s.EvalPartialU(0.5, 0.5), fd)
	}
	n, err := s.Normal(0.5, 0.5)
	if err != nil {
		t.Fatal(err)
	}
	// Normal must be perpendicular to both partials.
	if !close(n.Dot(s.EvalPartialU(0.5, 0.5)), 0) || !close(n.Dot(s.EvalPartialV(0.5, 0.5)), 0) {
		t.Errorf("normal not orthogonal to tangent plane")
	}
}

func TestBSplineSurfaceReducesToBezier(t *testing.T) {
	net := bezierTestNet()
	U := []float64{0, 0, 0, 1, 1, 1}
	V := []float64{0, 0, 0, 1, 1, 1}
	bss, err := NewBSplineSurface(net, U, V, 2, 2)
	if err != nil {
		t.Fatal(err)
	}
	bez, _ := NewBezierSurface(net)
	for _, u := range []float64{0, 0.3, 1} {
		for _, v := range []float64{0, 0.5, 1} {
			if !vecClose(bss.Eval(u, v), bez.Eval(u, v), 1e-12) {
				t.Errorf("B-spline surface != Bezier at (%v,%v)", u, v)
			}
		}
	}
}

func TestNURBSSurfaceReducesToBSpline(t *testing.T) {
	net := bezierTestNet()
	U := []float64{0, 0, 0, 1, 1, 1}
	V := []float64{0, 0, 0, 1, 1, 1}
	w := [][]float64{{1, 1, 1}, {1, 1, 1}, {1, 1, 1}}
	ns, err := NewNURBSSurface(net, w, U, V, 2, 2)
	if err != nil {
		t.Fatal(err)
	}
	bez, _ := NewBezierSurface(net)
	for _, u := range []float64{0.2, 0.5, 0.9} {
		for _, v := range []float64{0.1, 0.6} {
			if !vecClose(ns.Eval(u, v), bez.Eval(u, v), 1e-12) {
				t.Errorf("NURBS(w=1) surface != Bezier at (%v,%v)", u, v)
			}
		}
	}
}

// --- Util ---------------------------------------------------------------

func TestLinearAndStepInterpolate(t *testing.T) {
	x := []float64{0, 1, 2}
	y := []float64{0, 10, 30}
	if got, _ := LinearInterpolate(x, y, 0.5); !close(got, 5) {
		t.Errorf("linear=%v", got)
	}
	if got, _ := LinearInterpolate(x, y, 1.5); !close(got, 20) {
		t.Errorf("linear=%v", got)
	}
	if got, _ := LinearInterpolate(x, y, -1); !close(got, 0) {
		t.Errorf("clamp left=%v", got)
	}
	if got, _ := StepInterpolate(x, y, 1.9); !close(got, 10) {
		t.Errorf("step=%v", got)
	}
	if got, _ := NearestInterpolate(x, y, 1.4); !close(got, 10) {
		t.Errorf("nearest=%v", got)
	}
	if got, _ := NearestInterpolate(x, y, 1.6); !close(got, 30) {
		t.Errorf("nearest=%v", got)
	}
}

func TestChordAndCentripetalParams(t *testing.T) {
	pts := []Vec{NewVec(0, 0), NewVec(3, 4), NewVec(3, 4+5)}
	cl := ChordLengths(pts)
	if !close(cl[1], 5) || !close(cl[2], 10) {
		t.Errorf("chord lengths=%v", cl)
	}
	if got := PolylineLength(pts); !close(got, 10) {
		t.Errorf("polyline length=%v", got)
	}
	nc := NormalizedChordParams(pts)
	if !close(nc[0], 0) || !close(nc[2], 1) || !close(nc[1], 0.5) {
		t.Errorf("normalized params=%v", nc)
	}
}

func TestCatmullRomCurveCentripetal(t *testing.T) {
	pts := []Vec{NewVec(0, 0), NewVec(1, 1), NewVec(2, 0), NewVec(3, 1)}
	c, err := NewCatmullRomCurve(pts, 0.5)
	if err != nil {
		t.Fatal(err)
	}
	// Interpolates control points at integer parameters.
	for i := range pts {
		if !vecClose(c.Eval(float64(i)), pts[i], 1e-9) {
			t.Errorf("centripetal CR node %d = %v want %v", i, c.Eval(float64(i)), pts[i])
		}
	}
	// CatmullRomPoint uniform interpolates endpoints of a segment.
	p := CatmullRomPoint(pts[0], pts[1], pts[2], pts[3], 0)
	if !vecClose(p, pts[1], 1e-12) {
		t.Errorf("CatmullRomPoint t=0 = %v want %v", p, pts[1])
	}
}

// --- Error paths --------------------------------------------------------

func TestConstructorErrors(t *testing.T) {
	if _, err := NaturalCubicSpline([]float64{0}, []float64{1}); err == nil {
		t.Error("expected too-few-points error")
	}
	if _, err := NaturalCubicSpline([]float64{0, 0}, []float64{1, 2}); err == nil {
		t.Error("expected non-increasing error")
	}
	if _, err := NewBSplineCurve([]Vec{NewVec(0)}, []float64{0, 1}, 3); err == nil {
		t.Error("expected too-few-points error for B-spline")
	}
	if _, err := NURBSCircle(NewVec(0), 1); err == nil {
		t.Error("expected dimension error for circle")
	}
	if _, err := PeriodicCubicSpline([]float64{0, 1, 2}, []float64{0, 1, 5}); err == nil {
		t.Error("expected periodic endpoint mismatch error")
	}
}

// --- Example ------------------------------------------------------------

func ExampleBezierCurve() {
	b, _ := NewBezierCurve(
		NewVec(0, 0), NewVec(1, 2), NewVec(3, 2), NewVec(4, 0),
	)
	p := b.Eval(0.5)
	fmt.Printf("point at t=0.5: (%.1f, %.1f)\n", p[0], p[1])
	fmt.Printf("degree after elevation: %d\n", b.Elevate().Degree())
	// Output:
	// point at t=0.5: (2.0, 1.5)
	// degree after elevation: 4
}

func ExampleNURBSCircle() {
	c, _ := NURBSCircle(NewVec(0, 0), 1)
	p := c.Eval(0.25) // quarter turn
	fmt.Printf("(%.3f, %.3f)\n", p[0], p[1])
	// Output:
	// (0.000, 1.000)
}
