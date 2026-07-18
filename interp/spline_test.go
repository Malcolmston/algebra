package interp

import (
	"math"
	"testing"
)

func interpClose(t *testing.T, name string, got, want, tol float64) {
	t.Helper()
	if math.IsNaN(got) || math.Abs(got-want) > tol {
		t.Errorf("%s = %v, want %v (tol %v)", name, got, want, tol)
	}
}

// ---------------------------------------------------------------------------
// Scalar helpers
// ---------------------------------------------------------------------------

func TestScalarHelpers(t *testing.T) {
	interpClose(t, "Lerp mid", Lerp(2, 4, 0.5), 3, 1e-15)
	interpClose(t, "Lerp extrap", Lerp(2, 4, 2), 6, 1e-15)
	interpClose(t, "InverseLerp", InverseLerp(2, 4, 3), 0.5, 1e-15)
	interpClose(t, "InverseLerp deg", InverseLerp(5, 5, 5), 0, 1e-15)
	interpClose(t, "LinearAt", LinearAt(0, 0, 2, 10, 1), 5, 1e-15)
	interpClose(t, "Clamp lo", Clamp(-1, 0, 1), 0, 1e-15)
	interpClose(t, "Clamp hi", Clamp(3, 0, 1), 1, 1e-15)
	interpClose(t, "Clamp in", Clamp(0.4, 0, 1), 0.4, 1e-15)

	interpClose(t, "SmoothStep lo", SmoothStep(0, 1, -1), 0, 1e-15)
	interpClose(t, "SmoothStep hi", SmoothStep(0, 1, 2), 1, 1e-15)
	interpClose(t, "SmoothStep mid", SmoothStep(0, 1, 0.5), 0.5, 1e-15)
	interpClose(t, "SmootherStep mid", SmootherStep(0, 1, 0.5), 0.5, 1e-15)
	interpClose(t, "SmootherStep lo", SmootherStep(0, 1, 0), 0, 1e-15)
}

func TestBilinearTrilinearAt(t *testing.T) {
	// f = 1 + 2x + 3y + 4xy on the unit square, reproduced exactly.
	f := func(x, y float64) float64 { return 1 + 2*x + 3*y + 4*x*y }
	got := BilinearAt(0, 1, 0, 1, f(0, 0), f(1, 0), f(0, 1), f(1, 1), 0.3, 0.7)
	interpClose(t, "BilinearAt", got, f(0.3, 0.7), 1e-13)

	g := func(x, y, z float64) float64 { return x + 2*y + 3*z + x*y*z }
	c := func(x, y, z float64) float64 { return g(x, y, z) }
	gt := TrilinearAt(0, 1, 0, 1, 0, 1,
		c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(1, 1, 0),
		c(0, 0, 1), c(1, 0, 1), c(0, 1, 1), c(1, 1, 1),
		0.2, 0.6, 0.9)
	interpClose(t, "TrilinearAt", gt, g(0.2, 0.6, 0.9), 1e-13)
}

// ---------------------------------------------------------------------------
// LinearInterp
// ---------------------------------------------------------------------------

func TestLinearInterp(t *testing.T) {
	x := []float64{0, 1, 2, 4}
	y := []float64{1, 3, 5, 9} // exactly y = 2x + 1
	l, err := NewLinearInterp(x, y)
	if err != nil {
		t.Fatal(err)
	}
	for _, xq := range []float64{0, 0.5, 1.5, 3, 4} {
		interpClose(t, "LinearInterp", l.Eval(xq), 2*xq+1, 1e-13)
	}
	interpClose(t, "LinearInterp extrap", l.Eval(5), 11, 1e-13)
	interpClose(t, "LinearInterp deriv", l.EvalDerivative(1.5), 2, 1e-13)
	// integral of 2x+1 from 0 to 4 = x^2+x = 16+4 = 20
	interpClose(t, "LinearInterp integral", l.Integral(0, 4), 20, 1e-12)
	interpClose(t, "LinearInterp integral rev", l.Integral(4, 0), -20, 1e-12)
	if l.Len() != 4 {
		t.Errorf("Len = %d", l.Len())
	}
	lo, hi := l.Domain()
	if lo != 0 || hi != 4 {
		t.Errorf("Domain = %v,%v", lo, hi)
	}
	v, err := LinearInterpolate(x, y, 2.5)
	if err != nil {
		t.Fatal(err)
	}
	interpClose(t, "LinearInterpolate", v, 6, 1e-13)
}

func TestLinearInterpErrors(t *testing.T) {
	if _, err := NewLinearInterp([]float64{0}, []float64{0}); err != ErrTooFewPoints {
		t.Errorf("want ErrTooFewPoints, got %v", err)
	}
	if _, err := NewLinearInterp([]float64{0, 1}, []float64{0}); err != ErrLengthMismatch {
		t.Errorf("want ErrLengthMismatch, got %v", err)
	}
	if _, err := NewLinearInterp([]float64{0, 0}, []float64{0, 1}); err != ErrNotSorted {
		t.Errorf("want ErrNotSorted, got %v", err)
	}
}

// ---------------------------------------------------------------------------
// Nearest and step
// ---------------------------------------------------------------------------

func TestNearestInterp(t *testing.T) {
	x := []float64{0, 1, 2, 3}
	y := []float64{10, 20, 30, 40}
	n, _ := NewNearestInterp(x, y)
	interpClose(t, "nearest below", n.Eval(-5), 10, 0)
	interpClose(t, "nearest above", n.Eval(9), 40, 0)
	interpClose(t, "nearest 0.4", n.Eval(0.4), 10, 0)
	interpClose(t, "nearest 0.6", n.Eval(0.6), 20, 0)
	interpClose(t, "nearest tie", n.Eval(1.5), 20, 0) // tie -> lower index
	v, _ := NearestInterpolate(x, y, 2.9)
	interpClose(t, "NearestInterpolate", v, 40, 0)
}

func TestStepInterp(t *testing.T) {
	x := []float64{0, 1, 2, 3}
	y := []float64{10, 20, 30, 40}
	prev, _ := NewPreviousStepInterp(x, y)
	interpClose(t, "prev 0.5", prev.Eval(0.5), 10, 0)
	interpClose(t, "prev 1.0", prev.Eval(1.0), 20, 0)
	interpClose(t, "prev 2.9", prev.Eval(2.9), 30, 0)
	interpClose(t, "prev below", prev.Eval(-1), 10, 0)

	next, _ := NewNextStepInterp(x, y)
	interpClose(t, "next 0.5", next.Eval(0.5), 20, 0)
	interpClose(t, "next 1.0", next.Eval(1.0), 20, 0)
	interpClose(t, "next 2.1", next.Eval(2.1), 40, 0)
	interpClose(t, "next above", next.Eval(9), 40, 0)
}

// ---------------------------------------------------------------------------
// Cubic splines
// ---------------------------------------------------------------------------

func TestNaturalCubicReproducesLine(t *testing.T) {
	x := []float64{0, 1, 2, 3, 5}
	y := make([]float64, len(x))
	for i, xi := range x {
		y[i] = 3*xi - 2
	}
	s, err := NewNaturalCubicSpline(x, y)
	if err != nil {
		t.Fatal(err)
	}
	for _, xq := range []float64{0, 0.3, 1.7, 2.5, 4.9} {
		interpClose(t, "natural line", s.Eval(xq), 3*xq-2, 1e-12)
		interpClose(t, "natural line d2", s.EvalSecondDerivative(xq), 0, 1e-10)
	}
	// integral of 3x-2 over [0,5] = 3/2*25 - 10 = 27.5
	interpClose(t, "natural integral", s.Integral(0, 5), 27.5, 1e-11)
}

func TestNotAKnotReproducesCubic(t *testing.T) {
	// A not-a-knot spline reproduces any cubic exactly.
	poly := func(x float64) float64 { return 2*x*x*x - 3*x*x + x - 5 }
	dpoly := func(x float64) float64 { return 6*x*x - 6*x + 1 }
	d2poly := func(x float64) float64 { return 12*x - 6 }
	x := []float64{-2, -1, 0, 1, 2, 3, 4}
	y := make([]float64, len(x))
	for i, xi := range x {
		y[i] = poly(xi)
	}
	s, err := NewNotAKnotCubicSpline(x, y)
	if err != nil {
		t.Fatal(err)
	}
	for _, xq := range []float64{-1.5, -0.3, 0.5, 1.9, 3.4} {
		interpClose(t, "notaknot val", s.Eval(xq), poly(xq), 1e-9)
		interpClose(t, "notaknot d1", s.EvalDerivative(xq), dpoly(xq), 1e-8)
		interpClose(t, "notaknot d2", s.EvalSecondDerivative(xq), d2poly(xq), 1e-7)
		interpClose(t, "notaknot d3", s.EvalThirdDerivative(xq), 12, 1e-6)
	}
	// integral of the cubic from -2 to 4: F = x^4/2 - x^3 + x^2/2 - 5x
	F := func(x float64) float64 { return x*x*x*x/2 - x*x*x + x*x/2 - 5*x }
	interpClose(t, "notaknot integral", s.Integral(-2, 4), F(4)-F(-2), 1e-7)
}

func TestClampedReproducesCubic(t *testing.T) {
	poly := func(x float64) float64 { return x*x*x - 2*x }
	dpoly := func(x float64) float64 { return 3*x*x - 2 }
	x := []float64{0, 1, 2, 3, 4}
	y := make([]float64, len(x))
	for i, xi := range x {
		y[i] = poly(xi)
	}
	s, err := NewClampedCubicSpline(x, y, dpoly(0), dpoly(4))
	if err != nil {
		t.Fatal(err)
	}
	for _, xq := range []float64{0.25, 1.5, 2.7, 3.9} {
		interpClose(t, "clamped cubic", s.Eval(xq), poly(xq), 1e-9)
		interpClose(t, "clamped deriv", s.EvalDerivative(xq), dpoly(xq), 1e-8)
	}
	// endpoint slopes enforced exactly
	interpClose(t, "clamped left slope", s.EvalDerivative(0), dpoly(0), 1e-10)
	interpClose(t, "clamped right slope", s.EvalDerivative(4), dpoly(4), 1e-10)
}

func TestCubicNodesInterpolated(t *testing.T) {
	x := []float64{0, 0.5, 1.3, 2.0, 3.1, 4.0}
	y := []float64{0, 1, -1, 2, 0.5, 3}
	s, _ := NewNaturalCubicSpline(x, y)
	for i := range x {
		interpClose(t, "node", s.Eval(x[i]), y[i], 1e-10)
	}
	// C2 continuity: second derivative continuous at interior knots.
	for i := 1; i < len(x)-1; i++ {
		left := s.EvalSecondDerivative(x[i] - 1e-6)
		right := s.EvalSecondDerivative(x[i] + 1e-6)
		interpClose(t, "C2", left, right, 1e-3)
	}
}

func TestNotAKnotNeedsFour(t *testing.T) {
	if _, err := NewNotAKnotCubicSpline([]float64{0, 1, 2}, []float64{0, 1, 4}); err != ErrTooFewPoints {
		t.Errorf("want ErrTooFewPoints, got %v", err)
	}
}

// ---------------------------------------------------------------------------
// Akima
// ---------------------------------------------------------------------------

func TestAkimaReproducesLineAndNodes(t *testing.T) {
	x := []float64{0, 1, 2, 3, 4, 5}
	y := make([]float64, len(x))
	for i, xi := range x {
		y[i] = -2*xi + 7
	}
	s, err := NewAkimaSpline(x, y)
	if err != nil {
		t.Fatal(err)
	}
	for _, xq := range []float64{0.2, 1.5, 3.3, 4.8} {
		interpClose(t, "akima line", s.Eval(xq), -2*xq+7, 1e-12)
		interpClose(t, "akima line deriv", s.EvalDerivative(xq), -2, 1e-11)
	}
	// node interpolation for arbitrary data
	yd := []float64{0, 2, 1, 3, 2, 5}
	s2, _ := NewAkimaSpline(x, yd)
	for i := range x {
		interpClose(t, "akima node", s2.Eval(x[i]), yd[i], 1e-11)
	}
}

func TestAkimaFlatRegion(t *testing.T) {
	// Akima's signature property: a flat segment stays flat without overshoot.
	x := []float64{0, 1, 2, 3, 4, 5, 6, 7}
	y := []float64{0, 0, 0, 1, 1, 1, 1, 1}
	s, _ := NewAkimaSpline(x, y)
	// Within the leading flat region the interpolant stays at 0.
	for _, xq := range []float64{0.5, 1.5} {
		interpClose(t, "akima flat", s.Eval(xq), 0, 1e-12)
	}
	// Within the trailing flat region it stays at 1.
	for _, xq := range []float64{4.5, 5.5, 6.5} {
		interpClose(t, "akima flat hi", s.Eval(xq), 1, 1e-12)
	}
}

// ---------------------------------------------------------------------------
// PCHIP
// ---------------------------------------------------------------------------

func TestPCHIPReproducesLine(t *testing.T) {
	x := []float64{0, 1, 2, 3, 4}
	y := make([]float64, len(x))
	for i, xi := range x {
		y[i] = 4*xi + 1
	}
	s, err := NewPCHIP(x, y)
	if err != nil {
		t.Fatal(err)
	}
	for _, xq := range []float64{0.3, 1.6, 2.9, 3.5} {
		interpClose(t, "pchip line", s.Eval(xq), 4*xq+1, 1e-12)
		interpClose(t, "pchip line deriv", s.EvalDerivative(xq), 4, 1e-11)
	}
}

func TestPCHIPMonotone(t *testing.T) {
	x := []float64{0, 1, 2, 3, 4, 5}
	y := []float64{0, 0.1, 0.2, 5, 5.1, 5.2} // monotone increasing, sharp step
	s, _ := NewPCHIP(x, y)
	if !s.IsMonotone() {
		t.Error("IsMonotone should be true")
	}
	// No overshoot: sampled values are monotone and within global bounds.
	prev := math.Inf(-1)
	for i := 0; i <= 500; i++ {
		xq := 5 * float64(i) / 500
		v := s.Eval(xq)
		if v < -1e-12 || v > 5.2+1e-12 {
			t.Fatalf("overshoot at %v: %v", xq, v)
		}
		if v < prev-1e-12 {
			t.Fatalf("non-monotone at %v: %v < %v", xq, v, prev)
		}
		prev = v
	}
	// nodes interpolated exactly
	for i := range x {
		interpClose(t, "pchip node", s.Eval(x[i]), y[i], 1e-11)
	}
}

func TestPCHIPNonMonotone(t *testing.T) {
	x := []float64{0, 1, 2, 3}
	y := []float64{0, 1, 0, 1}
	s, _ := NewPCHIP(x, y)
	if s.IsMonotone() {
		t.Error("IsMonotone should be false")
	}
	// Local extrema at the interior peaks: slope forced to zero there.
	interpClose(t, "pchip peak slope", s.EvalDerivative(1), 0, 1e-12)
}

func TestPCHIPTwoPoints(t *testing.T) {
	s, err := NewPCHIP([]float64{0, 2}, []float64{1, 5})
	if err != nil {
		t.Fatal(err)
	}
	interpClose(t, "pchip 2pt", s.Eval(1), 3, 1e-12) // linear
}

// ---------------------------------------------------------------------------
// Grid interpolation
// ---------------------------------------------------------------------------

func TestBilinearGrid(t *testing.T) {
	x := []float64{0, 1, 2}
	y := []float64{0, 2, 4}
	f := func(x, y float64) float64 { return 1 + 2*x + 3*y + 4*x*y }
	z := make([][]float64, len(x))
	for i := range x {
		z[i] = make([]float64, len(y))
		for j := range y {
			z[i][j] = f(x[i], y[j])
		}
	}
	g, err := NewBilinearGrid(x, y, z)
	if err != nil {
		t.Fatal(err)
	}
	for _, p := range [][2]float64{{0.3, 0.7}, {1.5, 3.2}, {0.9, 1.1}} {
		interpClose(t, "bilinear", g.Eval(p[0], p[1]), f(p[0], p[1]), 1e-11)
	}
	// clamped query at the corner stays put
	interpClose(t, "bilinear clamped", g.EvalClamped(9, 9), f(2, 4), 1e-11)
	nx, ny := g.Dims()
	if nx != 3 || ny != 3 {
		t.Errorf("Dims = %d,%d", nx, ny)
	}
}

func TestBilinearGridErrors(t *testing.T) {
	if _, err := NewBilinearGrid([]float64{0}, []float64{0, 1}, [][]float64{{0}}); err != ErrTooFewPoints {
		t.Errorf("want ErrTooFewPoints, got %v", err)
	}
	bad := [][]float64{{0, 1}} // only 1 row but 2 x-values
	if _, err := NewBilinearGrid([]float64{0, 1}, []float64{0, 1}, bad); err != ErrGridShape {
		t.Errorf("want ErrGridShape, got %v", err)
	}
}

func TestTrilinearGrid(t *testing.T) {
	x := []float64{0, 1, 2}
	y := []float64{0, 1, 2}
	z := []float64{0, 1, 2}
	f := func(x, y, z float64) float64 { return 2 + x + 3*y + 5*z + x*y + y*z + x*z + x*y*z }
	v := make([][][]float64, len(x))
	for i := range x {
		v[i] = make([][]float64, len(y))
		for j := range y {
			v[i][j] = make([]float64, len(z))
			for k := range z {
				v[i][j][k] = f(x[i], y[j], z[k])
			}
		}
	}
	g, err := NewTrilinearGrid(x, y, z, v)
	if err != nil {
		t.Fatal(err)
	}
	for _, p := range [][3]float64{{0.3, 0.7, 0.2}, {1.5, 0.4, 1.9}, {0.9, 1.1, 0.5}} {
		interpClose(t, "trilinear", g.Eval(p[0], p[1], p[2]), f(p[0], p[1], p[2]), 1e-11)
	}
	interpClose(t, "trilinear clamped", g.EvalClamped(-1, -1, -1), f(0, 0, 0), 1e-11)
	nx, ny, nz := g.Dims()
	if nx != 3 || ny != 3 || nz != 3 {
		t.Errorf("Dims = %d,%d,%d", nx, ny, nz)
	}
}

// ---------------------------------------------------------------------------
// Search helper
// ---------------------------------------------------------------------------

func TestSearchInterval(t *testing.T) {
	xs := []float64{0, 1, 2, 3, 4}
	cases := []struct {
		x    float64
		want int
	}{
		{-1, 0}, {0, 0}, {0.5, 0}, {1, 1}, {2.5, 2}, {3.9, 3}, {4, 3}, {10, 3},
	}
	for _, c := range cases {
		if got := SearchInterval(xs, c.x); got != c.want {
			t.Errorf("SearchInterval(%v) = %d, want %d", c.x, got, c.want)
		}
	}
}

// ---------------------------------------------------------------------------
// Benchmark: building and evaluating a natural cubic spline (the heaviest
// routine, combining a tridiagonal solve with many polynomial evaluations).
// ---------------------------------------------------------------------------

func BenchmarkNaturalCubicSpline(b *testing.B) {
	const n = 1024
	x := make([]float64, n)
	y := make([]float64, n)
	for i := 0; i < n; i++ {
		x[i] = float64(i) * 0.1
		y[i] = math.Sin(x[i]) + 0.3*math.Cos(3*x[i])
	}
	b.ResetTimer()
	var acc float64
	for i := 0; i < b.N; i++ {
		s, err := NewNaturalCubicSpline(x, y)
		if err != nil {
			b.Fatal(err)
		}
		for q := 0; q < 4096; q++ {
			acc += s.Eval(float64(q) * 0.025)
		}
	}
	_ = acc
}
