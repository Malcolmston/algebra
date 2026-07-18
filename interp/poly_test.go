package interp

import (
	"math"
	"testing"
)

const polyTol = 1e-9

func polyClose(a, b, tol float64) bool { return math.Abs(a-b) <= tol }

// reference cubic used across several exactness tests: p(x) = 2 - x + 3x^2.
func polyRef(x float64) float64 { return 2 - x + 3*x*x }

func TestLagrangeExact(t *testing.T) {
	x := []float64{0, 1, 2, 3}
	y := make([]float64, len(x))
	for i, xi := range x {
		y[i] = polyRef(xi)
	}
	for _, xq := range []float64{-1, 0.5, 1.7, 2.5, 4} {
		got := LagrangeEval(x, y, xq)
		if !polyClose(got, polyRef(xq), polyTol) {
			t.Errorf("LagrangeEval(%v)=%v want %v", xq, got, polyRef(xq))
		}
	}
	// LagrangeBasis property: L_i(x_j) = delta_ij.
	for i := range x {
		for j := range x {
			got := LagrangeBasis(x, i, x[j])
			want := 0.0
			if i == j {
				want = 1
			}
			if !polyClose(got, want, polyTol) {
				t.Errorf("LagrangeBasis i=%d at x[%d]=%v want %v", i, j, got, want)
			}
		}
	}
	li, err := NewLagrangeInterpolator(x, y)
	if err != nil {
		t.Fatal(err)
	}
	if li.Degree() != 3 || li.Len() != 4 {
		t.Errorf("degree/len = %d/%d", li.Degree(), li.Len())
	}
	if !polyClose(li.Eval(0.5), polyRef(0.5), polyTol) {
		t.Errorf("interp eval mismatch")
	}
}

func TestNewtonDividedDifferences(t *testing.T) {
	x := []float64{0, 1, 2}
	y := []float64{1, 4, 9} // (x+1)^2
	c := DividedDifferences(x, y)
	want := []float64{1, 3, 1}
	for i := range want {
		if !polyClose(c[i], want[i], polyTol) {
			t.Errorf("coeff[%d]=%v want %v", i, c[i], want[i])
		}
	}
	tab := DividedDifferenceTable(x, y)
	if !polyClose(tab[0][2], 1, polyTol) || !polyClose(tab[1][1], 5, polyTol) {
		t.Errorf("table wrong: %v", tab)
	}
	for _, xq := range []float64{-1, 0.5, 1.5, 3} {
		got := NewtonEval(x, c, xq)
		w := (xq + 1) * (xq + 1)
		if !polyClose(got, w, polyTol) {
			t.Errorf("NewtonEval(%v)=%v want %v", xq, got, w)
		}
	}
	ni, err := NewNewtonInterpolator(x, y)
	if err != nil {
		t.Fatal(err)
	}
	if err := ni.AddPoint(3, 16); err != nil { // (3+1)^2 = 16
		t.Fatal(err)
	}
	if ni.Len() != 4 {
		t.Errorf("len after add = %d", ni.Len())
	}
	if !polyClose(ni.Eval(2.5), (2.5+1)*(2.5+1), polyTol) {
		t.Errorf("after AddPoint eval mismatch: %v", ni.Eval(2.5))
	}
}

func TestBarycentricMatchesLagrange(t *testing.T) {
	x := []float64{-1, 0, 2, 3, 5}
	y := make([]float64, len(x))
	for i, xi := range x {
		y[i] = math.Sin(xi) + xi*xi
	}
	bi, err := NewBarycentricInterpolator(x, y)
	if err != nil {
		t.Fatal(err)
	}
	for _, xq := range []float64{-0.5, 1.2, 2.5, 4} {
		if !polyClose(bi.Eval(xq), LagrangeEval(x, y, xq), 1e-8) {
			t.Errorf("barycentric != lagrange at %v", xq)
		}
	}
	// exact recovery at nodes
	for i := range x {
		if !polyClose(bi.Eval(x[i]), y[i], polyTol) {
			t.Errorf("node recovery failed at %d", i)
		}
	}
}

func TestHermiteCubic(t *testing.T) {
	// Hermite matching x^3 and derivative 3x^2 at {0,1} recovers x^3.
	x := []float64{0, 1}
	y := []float64{0, 1}
	dy := []float64{0, 3}
	hi, err := NewHermiteInterpolator(x, y, dy)
	if err != nil {
		t.Fatal(err)
	}
	for _, xq := range []float64{0.25, 0.5, 0.75} {
		if !polyClose(hi.Eval(xq), xq*xq*xq, polyTol) {
			t.Errorf("hermite eval at %v = %v want %v", xq, hi.Eval(xq), xq*xq*xq)
		}
		if !polyClose(hi.EvalDerivative(xq), 3*xq*xq, 1e-8) {
			t.Errorf("hermite deriv at %v = %v want %v", xq, hi.EvalDerivative(xq), 3*xq*xq)
		}
	}
	if !polyClose(HermiteEval(x, y, dy, 0.5), 0.125, polyTol) {
		t.Errorf("HermiteEval mismatch")
	}
	// CubicHermiteSegment endpoints and midpoint.
	if !polyClose(CubicHermiteSegment(1, 0, 2, 0, 0), 1, polyTol) ||
		!polyClose(CubicHermiteSegment(1, 0, 2, 0, 1), 2, polyTol) ||
		!polyClose(CubicHermiteSegment(1, 0, 2, 0, 0.5), 1.5, polyTol) {
		t.Errorf("CubicHermiteSegment wrong")
	}
}

func TestNeville(t *testing.T) {
	x := []float64{0, 1, 2, 3}
	y := make([]float64, len(x))
	for i, xi := range x {
		y[i] = polyRef(xi)
	}
	for _, xq := range []float64{0.5, 1.5, 2.5} {
		if !polyClose(NevilleEval(x, y, xq), polyRef(xq), polyTol) {
			t.Errorf("NevilleEval mismatch at %v", xq)
		}
		v, e := NevilleWithError(x, y, xq)
		if !polyClose(v, polyRef(xq), polyTol) {
			t.Errorf("NevilleWithError value mismatch at %v: %v", xq, v)
		}
		if e < 0 {
			t.Errorf("error estimate negative")
		}
	}
	tab := NevilleTable(x, y, 0.5)
	if len(tab) != 4 || len(tab[3]) != 1 {
		t.Fatalf("bad table shape")
	}
	if !polyClose(tab[3][0], polyRef(0.5), polyTol) {
		t.Errorf("NevilleTable final = %v want %v", tab[3][0], polyRef(0.5))
	}
}

func TestChebyshev(t *testing.T) {
	if !polyClose(ChebyshevT(0, 0.3), 1, polyTol) ||
		!polyClose(ChebyshevT(1, 0.3), 0.3, polyTol) ||
		!polyClose(ChebyshevT(2, 0.5), -0.5, polyTol) ||
		!polyClose(ChebyshevT(3, 0.5), -1, polyTol) {
		t.Errorf("ChebyshevT values wrong")
	}
	// Nodes lie in [-1,1] and are increasing.
	nodes := ChebyshevNodes(5)
	for i := 1; i < len(nodes); i++ {
		if nodes[i] <= nodes[i-1] {
			t.Errorf("nodes not increasing")
		}
	}
	if nodes[0] < -1 || nodes[len(nodes)-1] > 1 {
		t.Errorf("nodes out of range")
	}
	sk := ChebyshevNodesSecondKind(3)
	want := []float64{-1, 0, 1}
	for i := range want {
		if !polyClose(sk[i], want[i], polyTol) {
			t.Errorf("second-kind node %d = %v want %v", i, sk[i], want[i])
		}
	}
	// A cubic is reproduced exactly by 4-point Chebyshev interpolation.
	f := func(x float64) float64 { return x*x*x - 2*x + 1 }
	ci, err := NewChebyshevInterpolator(f, 4, -1, 1)
	if err != nil {
		t.Fatal(err)
	}
	for _, xq := range []float64{-0.9, -0.3, 0.2, 0.8} {
		if !polyClose(ci.Eval(xq), f(xq), 1e-8) {
			t.Errorf("cheb interp at %v = %v want %v", xq, ci.Eval(xq), f(xq))
		}
	}
	// Smooth function on a general interval.
	g := func(x float64) float64 { return math.Exp(x) }
	cg, _ := NewChebyshevInterpolator(g, 20, 0, 2)
	for _, xq := range []float64{0.1, 0.7, 1.3, 1.9} {
		if !polyClose(cg.Eval(xq), g(xq), 1e-9) {
			t.Errorf("cheb exp at %v = %v want %v", xq, cg.Eval(xq), g(xq))
		}
	}
	if len(ci.Coefficients()) != 4 {
		t.Errorf("coeff length")
	}
}

func TestBezier(t *testing.T) {
	if !polyClose(BernsteinBasis(2, 1, 0.5), 0.5, polyTol) {
		t.Errorf("BernsteinBasis wrong")
	}
	if !polyClose(DeCasteljau([]float64{0, 1}, 0.3), 0.3, polyTol) {
		t.Errorf("DeCasteljau linear wrong")
	}
	ctrl := [][]float64{{0, 0}, {1, 2}, {2, 0}}
	p := BezierPoint(ctrl, 0.5)
	if !polyClose(p[0], 1, polyTol) || !polyClose(p[1], 1, polyTol) {
		t.Errorf("BezierPoint = %v want [1 1]", p)
	}
	d := BezierDerivative(ctrl, 0.5)
	if !polyClose(d[0], 2, polyTol) || !polyClose(d[1], 0, polyTol) {
		t.Errorf("BezierDerivative = %v want [2 0]", d)
	}
	// Degree elevation preserves the curve.
	el := BezierElevate(ctrl)
	if len(el) != 4 {
		t.Fatalf("elevate len = %d", len(el))
	}
	for _, tt := range []float64{0.1, 0.4, 0.9} {
		a := BezierPoint(ctrl, tt)
		b := BezierPoint(el, tt)
		if !polyClose(a[0], b[0], polyTol) || !polyClose(a[1], b[1], polyTol) {
			t.Errorf("elevation changed curve at %v: %v vs %v", tt, a, b)
		}
	}
	bz, err := NewBezier(ctrl)
	if err != nil {
		t.Fatal(err)
	}
	if bz.Degree() != 2 {
		t.Errorf("bezier degree")
	}
	q := bz.Eval(0.5)
	if !polyClose(q[0], 1, polyTol) {
		t.Errorf("bezier eval")
	}
}

func TestBSpline(t *testing.T) {
	// Clamped cubic single-segment (Bezier equivalent): interpolates ends.
	ctrl := [][]float64{{0, 0}, {1, 3}, {3, 3}, {4, 0}}
	knots := OpenUniformKnots(4, 3)
	if len(knots) != 8 {
		t.Fatalf("knot len = %d", len(knots))
	}
	bs, err := NewBSpline(3, knots, ctrl)
	if err != nil {
		t.Fatal(err)
	}
	lo, hi := bs.Domain()
	start := bs.Eval(lo)
	end := bs.Eval(hi)
	if !polyClose(start[0], 0, 1e-8) || !polyClose(start[1], 0, 1e-8) {
		t.Errorf("start = %v want [0 0]", start)
	}
	if !polyClose(end[0], 4, 1e-8) || !polyClose(end[1], 0, 1e-8) {
		t.Errorf("end = %v want [4 0]", end)
	}
	// Clamped cubic with these knots equals the Bezier curve of same controls.
	for _, tt := range []float64{0.2, 0.5, 0.85} {
		bp := BezierPoint(ctrl, tt)
		sp := bs.Eval(tt)
		if !polyClose(bp[0], sp[0], 1e-8) || !polyClose(bp[1], sp[1], 1e-8) {
			t.Errorf("bspline != bezier at %v: %v vs %v", tt, sp, bp)
		}
	}
	// Partition of unity of the basis functions at an interior parameter.
	deg := 2
	kc := 5
	kn := OpenUniformKnots(kc, deg)
	for _, tt := range []float64{0.1, 0.37, 0.63, 0.9} {
		var sum float64
		for i := 0; i < kc; i++ {
			sum += BSplineBasis(i, deg, kn, tt)
		}
		if !polyClose(sum, 1, 1e-9) {
			t.Errorf("partition of unity at %v = %v", tt, sum)
		}
	}
	// Uniform knots produce a valid domain.
	uk := UniformKnots(5, 2)
	if len(uk) != 8 || uk[2] != 2 || uk[5] != 5 {
		t.Errorf("UniformKnots wrong: %v", uk)
	}
}

func TestThieleRational(t *testing.T) {
	// f(x) = (1 + 2x)/(1 + x), a degree (1,1) rational, exact from 3 points.
	f := func(x float64) float64 { return (1 + 2*x) / (1 + x) }
	x := []float64{0, 1, 2}
	y := make([]float64, len(x))
	for i, xi := range x {
		y[i] = f(xi)
	}
	ti, err := NewThieleInterpolator(x, y)
	if err != nil {
		t.Fatal(err)
	}
	for _, xq := range []float64{0.5, 1.5, 3, 5} {
		if !polyClose(ti.Eval(xq), f(xq), 1e-9) {
			t.Errorf("Thiele at %v = %v want %v", xq, ti.Eval(xq), f(xq))
		}
	}
	// Bulirsch-Stoer rational interpolation of a non-rational smooth function
	// agrees with the (unique) Thiele rational interpolant through the same
	// points, and reproduces the samples at the nodes.
	gx := []float64{0, 0.5, 1, 1.5, 2}
	gy := make([]float64, len(gx))
	for i, xi := range gx {
		gy[i] = math.Exp(xi)
	}
	gt, err := NewThieleInterpolator(gx, gy)
	if err != nil {
		t.Fatal(err)
	}
	for _, xq := range []float64{0.2, 0.8, 1.3, 1.8} {
		v, err := RationalInterpEval(gx, gy, xq)
		if err != nil {
			t.Fatal(err)
		}
		if !polyClose(v, gt.Eval(xq), 1e-7) {
			t.Errorf("RationalInterpEval at %v = %v want %v (Thiele)", xq, v, gt.Eval(xq))
		}
	}
	// Node query returns node value.
	v, err := RationalInterpEval(gx, gy, 1)
	if err != nil {
		t.Fatal(err)
	}
	if !polyClose(v, math.Exp(1), 1e-9) {
		t.Errorf("rational node value wrong: %v", v)
	}
}

func TestLeastSquaresFit(t *testing.T) {
	// Exact quadratic 1 + 2x + 3x^2 over-determined recovers coefficients.
	coeffsTrue := []float64{1, 2, 3}
	x := []float64{-2, -1, 0, 1, 2, 3}
	y := make([]float64, len(x))
	for i, xi := range x {
		y[i] = PolyVal(coeffsTrue, xi)
	}
	c, err := PolyFit(x, y, 2)
	if err != nil {
		t.Fatal(err)
	}
	for i := range coeffsTrue {
		if !polyClose(c[i], coeffsTrue[i], 1e-6) {
			t.Errorf("PolyFit coeff[%d]=%v want %v", i, c[i], coeffsTrue[i])
		}
	}
	if r2 := RSquared(x, y, c); !polyClose(r2, 1, 1e-9) {
		t.Errorf("RSquared = %v want 1", r2)
	}
	if e := RMSE(x, y, c); e > 1e-6 {
		t.Errorf("RMSE = %v want ~0", e)
	}
	// Weighted fit reduces to ordinary with unit weights.
	w := make([]float64, len(x))
	for i := range w {
		w[i] = 1
	}
	cw, err := PolyFitWeighted(x, y, w, 2)
	if err != nil {
		t.Fatal(err)
	}
	for i := range c {
		if !polyClose(cw[i], c[i], 1e-9) {
			t.Errorf("weighted != ordinary")
		}
	}
	// Vandermonde exact interpolation.
	vx := []float64{0, 1, 2}
	vy := []float64{1, 4, 9}
	vc, err := VandermondeSolve(vx, vy)
	if err != nil {
		t.Fatal(err)
	}
	want := []float64{1, 2, 1} // (x+1)^2
	for i := range want {
		if !polyClose(vc[i], want[i], 1e-9) {
			t.Errorf("Vandermonde coeff[%d]=%v want %v", i, vc[i], want[i])
		}
	}
	// Line fit.
	lx := []float64{0, 1, 2, 3}
	ly := []float64{1, 3, 5, 7} // 1 + 2x
	slope, intercept, err := LeastSquaresLine(lx, ly)
	if err != nil {
		t.Fatal(err)
	}
	if !polyClose(slope, 2, 1e-9) || !polyClose(intercept, 1, 1e-9) {
		t.Errorf("line fit = %v x + %v", slope, intercept)
	}
	// Polynomial calculus helpers.
	dc := PolyDerivativeCoeffs([]float64{1, 2, 3}) // d/dx (1+2x+3x^2) = 2+6x
	if !polyClose(dc[0], 2, polyTol) || !polyClose(dc[1], 6, polyTol) {
		t.Errorf("PolyDerivativeCoeffs = %v", dc)
	}
	ic := PolyIntegralCoeffs([]float64{2, 6}, 1) // integral = 1 + 2x + 3x^2
	if !polyClose(ic[0], 1, polyTol) || !polyClose(ic[1], 2, polyTol) || !polyClose(ic[2], 3, polyTol) {
		t.Errorf("PolyIntegralCoeffs = %v", ic)
	}
	pf, err := NewPolynomialFit(x, y, 2)
	if err != nil {
		t.Fatal(err)
	}
	if pf.Degree() != 2 || !polyClose(pf.Eval(1.5), PolyVal(coeffsTrue, 1.5), 1e-6) {
		t.Errorf("PolynomialFit eval wrong")
	}
}

func TestTrigInterp(t *testing.T) {
	// Band-limited f with frequencies below Nyquist is reproduced exactly.
	f := func(x float64) float64 { return 1 + math.Cos(x) + math.Sin(2*x) }
	N := 8
	y := make([]float64, N)
	for j := 0; j < N; j++ {
		y[j] = f(2 * math.Pi * float64(j) / float64(N))
	}
	ti, err := NewTrigInterpolator(y)
	if err != nil {
		t.Fatal(err)
	}
	for _, xq := range []float64{0.3, 1.1, 2.7, 4.5, 5.9} {
		if !polyClose(ti.Eval(xq), f(xq), 1e-9) {
			t.Errorf("trig interp at %v = %v want %v", xq, ti.Eval(xq), f(xq))
		}
	}
	// Interpolation at the sample nodes.
	for j := 0; j < N; j++ {
		xq := 2 * math.Pi * float64(j) / float64(N)
		if !polyClose(TrigInterpEval(y, xq), y[j], 1e-9) {
			t.Errorf("trig node %d mismatch", j)
		}
	}
	a, b := TrigCoefficients(y)
	if !polyClose(a[0]/2, 1, 1e-9) || !polyClose(a[1], 1, 1e-9) || !polyClose(b[2], 1, 1e-9) {
		t.Errorf("TrigCoefficients wrong: a=%v b=%v", a, b)
	}
	// DFT of a constant signal concentrates all energy in bin 0.
	re, im := DFTReal([]float64{2, 2, 2, 2})
	if !polyClose(re[0], 8, 1e-9) || !polyClose(im[0], 0, 1e-9) {
		t.Errorf("DFT bin 0 = %v+%vi", re[0], im[0])
	}
	for k := 1; k < 4; k++ {
		if !polyClose(re[k], 0, 1e-9) || !polyClose(im[k], 0, 1e-9) {
			t.Errorf("DFT bin %d nonzero", k)
		}
	}
}

func TestErrorPaths(t *testing.T) {
	if _, err := NewLagrangeInterpolator([]float64{1, 2}, []float64{1}); err == nil {
		t.Errorf("expected length mismatch error")
	}
	if _, err := NewNewtonInterpolator([]float64{1, 1}, []float64{1, 2}); err == nil {
		t.Errorf("expected distinct error")
	}
	if _, err := PolyFit([]float64{1, 2}, []float64{1, 2}, 5); err == nil {
		t.Errorf("expected too-few error")
	}
	if _, err := NewBSpline(0, []float64{0, 1}, [][]float64{{0}, {1}}); err == nil {
		t.Errorf("expected degree error")
	}
}

func BenchmarkNewChebyshevInterpolator(b *testing.B) {
	f := func(x float64) float64 { return math.Exp(-x) * math.Sin(3*x) }
	b.ReportAllocs()
	var sink float64
	for i := 0; i < b.N; i++ {
		ci, _ := NewChebyshevInterpolator(f, 64, 0, 5)
		sink += ci.Eval(2.5)
	}
	_ = sink
}
