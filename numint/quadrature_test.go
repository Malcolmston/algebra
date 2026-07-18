package numint

import (
	"math"
	"testing"
)

// approx reports whether got is within tol of want.
func approx(got, want, tol float64) bool {
	return math.Abs(got-want) <= tol
}

// Reference integrands with closed-form integrals used across the tests.
var (
	// f1(x) = x^2, integral on [0,1] = 1/3.
	f1 = func(x float64) float64 { return x * x }
	// fexp(x) = e^x, integral on [0,1] = e-1.
	fexp = func(x float64) float64 { return math.Exp(x) }
	// fsin(x) = sin(x), integral on [0,pi] = 2.
	fsin = func(x float64) float64 { return math.Sin(x) }
	// fpi(x) = 4/(1+x^2), integral on [0,1] = pi.
	fpi = func(x float64) float64 { return 4 / (1 + x*x) }
)

func TestSinglePanelRules(t *testing.T) {
	// Each rule is exact for polynomials up to its degree; test on x^2 over
	// [0,1] except the trapezoid/midpoint which are not exact for x^2.
	cases := []struct {
		name string
		got  float64
		want float64
		tol  float64
	}{
		{"SimpsonRule x^2", SimpsonRule(f1, 0, 1), 1.0 / 3.0, 1e-14},
		{"Simpson38Rule x^2", Simpson38Rule(f1, 0, 1), 1.0 / 3.0, 1e-14},
		{"BooleRule x^2", BooleRule(f1, 0, 1), 1.0 / 3.0, 1e-14},
		{"WeddleRule x^2", WeddleRule(f1, 0, 1), 1.0 / 3.0, 1e-13},
		{"TrapezoidRule linear", TrapezoidRule(func(x float64) float64 { return 2*x + 1 }, 0, 1), 2.0, 1e-14},
		{"MidpointRule linear", MidpointRule(func(x float64) float64 { return 2*x + 1 }, 0, 1), 2.0, 1e-14},
	}
	for _, c := range cases {
		if !approx(c.got, c.want, c.tol) {
			t.Errorf("%s = %.15g, want %.15g", c.name, c.got, c.want)
		}
	}
}

func TestCompositeRules(t *testing.T) {
	e1 := math.E - 1
	cases := []struct {
		name string
		got  float64
		want float64
		tol  float64
	}{
		{"Trapezoid e^x", Trapezoid(fexp, 0, 1, 2000), e1, 1e-6},
		{"Midpoint e^x", Midpoint(fexp, 0, 1, 2000), e1, 1e-6},
		{"Simpson e^x", Simpson(fexp, 0, 1, 100), e1, 1e-10},
		{"Simpson38 e^x", Simpson38(fexp, 0, 1, 99), e1, 1e-8},
		{"Boole e^x", Boole(fexp, 0, 1, 100), e1, 1e-12},
		{"Weddle e^x", Weddle(fexp, 0, 1, 96), e1, 1e-12},
		{"Simpson sin", Simpson(fsin, 0, math.Pi, 100), 2, 1e-7},
	}
	for _, c := range cases {
		if !approx(c.got, c.want, c.tol) {
			t.Errorf("%s = %.15g, want %.15g", c.name, c.got, c.want)
		}
	}
}

func TestSampleBased(t *testing.T) {
	n := 401
	xs := make([]float64, n)
	ys := make([]float64, n)
	h := 1.0 / float64(n-1)
	for i := 0; i < n; i++ {
		xs[i] = float64(i) * h
		ys[i] = math.Exp(xs[i])
	}
	e1 := math.E - 1
	if got := TrapezoidSamples(ys, h); !approx(got, e1, 1e-5) {
		t.Errorf("TrapezoidSamples = %.15g want %.15g", got, e1)
	}
	if got := SimpsonSamples(ys, h); !approx(got, e1, 1e-10) {
		t.Errorf("SimpsonSamples = %.15g want %.15g", got, e1)
	}
	if got := TrapezoidData(xs, ys); !approx(got, e1, 1e-5) {
		t.Errorf("TrapezoidData = %.15g want %.15g", got, e1)
	}

	// Midpoint samples: integrate x^2 on [0,1] via 1000 midpoints.
	m := 1000
	mh := 1.0 / float64(m)
	mid := make([]float64, m)
	for i := 0; i < m; i++ {
		xm := (float64(i) + 0.5) * mh
		mid[i] = xm * xm
	}
	if got := MidpointSamples(mid, mh); !approx(got, 1.0/3.0, 1e-6) {
		t.Errorf("MidpointSamples = %.15g want %.15g", got, 1.0/3.0)
	}
}

func TestCumulative(t *testing.T) {
	n := 201
	xs := make([]float64, n)
	ys := make([]float64, n)
	h := math.Pi / float64(n-1)
	for i := 0; i < n; i++ {
		xs[i] = float64(i) * h
		ys[i] = math.Sin(xs[i])
	}
	ct := CumulativeTrapezoid(xs, ys)
	if ct[0] != 0 {
		t.Errorf("CumulativeTrapezoid[0] = %g, want 0", ct[0])
	}
	// Full integral of sin over [0,pi] = 2.
	if !approx(ct[n-1], 2, 1e-3) {
		t.Errorf("CumulativeTrapezoid total = %.15g, want 2", ct[n-1])
	}
	// Cumulative to pi/2 (index n/2) should be about 1.
	if !approx(ct[(n-1)/2], 1, 1e-3) {
		t.Errorf("CumulativeTrapezoid mid = %.15g, want 1", ct[(n-1)/2])
	}

	cu := CumulativeTrapezoidUniform(ys, h)
	if !approx(cu[n-1], 2, 1e-3) {
		t.Errorf("CumulativeTrapezoidUniform total = %.15g, want 2", cu[n-1])
	}

	cs := CumulativeSimpson(xs, ys)
	if !approx(cs[n-1], 2, 1e-6) {
		t.Errorf("CumulativeSimpson total = %.15g, want 2", cs[n-1])
	}
}

func TestRomberg(t *testing.T) {
	if got := Romberg(fpi, 0, 1, 8); !approx(got, math.Pi, 1e-10) {
		t.Errorf("Romberg = %.15g, want %.15g", got, math.Pi)
	}
	if got := RombergTol(fexp, 0, 1, 1e-11, 12); !approx(got, math.E-1, 1e-10) {
		t.Errorf("RombergTol = %.15g, want %.15g", got, math.E-1)
	}
	v, table := RombergTable(fsin, 0, math.Pi, 8)
	if !approx(v, 2, 1e-9) {
		t.Errorf("RombergTable value = %.15g, want 2", v)
	}
	if len(table) != 8 {
		t.Errorf("RombergTable rows = %d, want 8", len(table))
	}
}

func TestAdaptive(t *testing.T) {
	if got := AdaptiveTrapezoid(fpi, 0, 1, 1e-9); !approx(got, math.Pi, 1e-7) {
		t.Errorf("AdaptiveTrapezoid = %.15g, want %.15g", got, math.Pi)
	}
	if got := AdaptiveSimpson(fpi, 0, 1, 1e-12); !approx(got, math.Pi, 1e-10) {
		t.Errorf("AdaptiveSimpson = %.15g, want %.15g", got, math.Pi)
	}
	// Integrand with a sharp peak: 1/(1+ (10(x-0.3))^2 ).
	peak := func(x float64) float64 { return 1 / (1 + 100*(x-0.3)*(x-0.3)) }
	want := 0.1 * (math.Atan(10*(1-0.3)) - math.Atan(10*(0-0.3)))
	res := AdaptiveGaussKronrod(peak, 0, 1, 1e-11)
	if !approx(res.Value, want, 1e-9) {
		t.Errorf("AdaptiveGaussKronrod = %.15g, want %.15g", res.Value, want)
	}
	if res.Evals <= 0 {
		t.Errorf("AdaptiveGaussKronrod Evals = %d, want > 0", res.Evals)
	}
}

func TestGaussLegendreFixed(t *testing.T) {
	// n-point Gauss-Legendre is exact for polynomials of degree 2n-1.
	// Use e^x on [0,1]; higher order should be markedly more accurate.
	e1 := math.E - 1
	cases := []struct {
		got float64
		tol float64
	}{
		{GaussLegendre2(fexp, 0, 1), 1e-3},
		{GaussLegendre3(fexp, 0, 1), 1e-6},
		{GaussLegendre4(fexp, 0, 1), 1e-8},
		{GaussLegendre5(fexp, 0, 1), 1e-10},
		{GaussLegendre6(fexp, 0, 1), 1e-12},
		{GaussLegendre7(fexp, 0, 1), 1e-13},
		{GaussLegendre8(fexp, 0, 1), 1e-14},
		{GaussLegendre9(fexp, 0, 1), 1e-14},
		{GaussLegendre10(fexp, 0, 1), 1e-14},
	}
	for i, c := range cases {
		if !approx(c.got, e1, c.tol) {
			t.Errorf("GaussLegendre%d = %.15g, want %.15g", i+2, c.got, e1)
		}
	}

	// Exactness check: GaussLegendre3 integrates a degree-5 polynomial
	// exactly. p(x) = x^5, integral on [0,2] = 2^6/6 = 64/6.
	p5 := func(x float64) float64 { return x * x * x * x * x }
	if got := GaussLegendre3(p5, 0, 2); !approx(got, 64.0/6.0, 1e-12) {
		t.Errorf("GaussLegendre3 x^5 = %.15g, want %.15g", got, 64.0/6.0)
	}
}

func TestGaussLegendreN(t *testing.T) {
	e1 := math.E - 1
	if got := GaussLegendreN(fexp, 0, 1, 5); !approx(got, e1, 1e-10) {
		t.Errorf("GaussLegendreN(5) = %.15g, want %.15g", got, e1)
	}
	// High order computed via Newton iteration (n=16 not tabulated).
	if got := GaussLegendreN(fexp, 0, 1, 16); !approx(got, e1, 1e-14) {
		t.Errorf("GaussLegendreN(16) = %.15g, want %.15g", got, e1)
	}
	// Nodes/weights on [-1,1] must sum to 2 and be symmetric.
	nodes, weights := GaussLegendreNodes(16)
	sw := 0.0
	for _, w := range weights {
		sw += w
	}
	if !approx(sw, 2, 1e-13) {
		t.Errorf("sum of GaussLegendreNodes(16) weights = %.15g, want 2", sw)
	}
	if len(nodes) != 16 {
		t.Errorf("GaussLegendreNodes(16) len = %d, want 16", len(nodes))
	}
	if got := CompositeGaussLegendre(fpi, 0, 1, 4, 5); !approx(got, math.Pi, 1e-10) {
		t.Errorf("CompositeGaussLegendre = %.15g, want %.15g", got, math.Pi)
	}
}

func TestGaussLobatto(t *testing.T) {
	e1 := math.E - 1
	if got := GaussLobatto3(fexp, 0, 1); !approx(got, e1, 1e-3) {
		t.Errorf("GaussLobatto3 = %.15g, want %.15g", got, e1)
	}
	if got := GaussLobatto4(fexp, 0, 1); !approx(got, e1, 1e-5) {
		t.Errorf("GaussLobatto4 = %.15g, want %.15g", got, e1)
	}
	if got := GaussLobatto5(fexp, 0, 1); !approx(got, e1, 1e-7) {
		t.Errorf("GaussLobatto5 = %.15g, want %.15g", got, e1)
	}
	// GaussLobatto5 is exact up to degree 7; test x^6 on [0,1] = 1/7.
	p6 := func(x float64) float64 { return math.Pow(x, 6) }
	if got := GaussLobatto5(p6, 0, 1); !approx(got, 1.0/7.0, 1e-13) {
		t.Errorf("GaussLobatto5 x^6 = %.15g, want %.15g", got, 1.0/7.0)
	}
}

func TestGaussKronrod15(t *testing.T) {
	v, e := GaussKronrod15(fpi, 0, 1)
	if !approx(v, math.Pi, 1e-9) {
		t.Errorf("GaussKronrod15 = %.15g, want %.15g", v, math.Pi)
	}
	if e < 0 {
		t.Errorf("GaussKronrod15 error estimate = %.3g, want >= 0", e)
	}
	// The error estimate should bound the true error (loosely).
	if math.Abs(v-math.Pi) > e+1e-12 {
		t.Errorf("true error %.3g exceeds estimate %.3g", math.Abs(v-math.Pi), e)
	}
	abs, wk, wg := GaussKronrod15Nodes()
	if len(abs) != 8 || len(wk) != 8 || len(wg) != 4 {
		t.Errorf("GaussKronrod15Nodes lengths = %d,%d,%d", len(abs), len(wk), len(wg))
	}
}

func TestClenshawCurtisAndChebyshev(t *testing.T) {
	// Clenshaw-Curtis of e^x on [-1,1] = 2*sinh(1).
	want := 2 * math.Sinh(1)
	if got := ClenshawCurtis(fexp, -1, 1, 32); !approx(got, want, 1e-12) {
		t.Errorf("ClenshawCurtis = %.15g, want %.15g", got, want)
	}
	// Gauss-Chebyshev: integral of 1/sqrt(1-x^2) over [-1,1] = pi.
	one := func(float64) float64 { return 1 }
	if got := GaussChebyshev(one, 20); !approx(got, math.Pi, 1e-12) {
		t.Errorf("GaussChebyshev(1) = %.15g, want %.15g", got, math.Pi)
	}
	// integral of x^2/sqrt(1-x^2) over [-1,1] = pi/2.
	if got := GaussChebyshev(f1, 20); !approx(got, math.Pi/2, 1e-12) {
		t.Errorf("GaussChebyshev(x^2) = %.15g, want %.15g", got, math.Pi/2)
	}
}

func TestTanhSinh(t *testing.T) {
	// Smooth integrand.
	if got := TanhSinh(fexp, 0, 1, 100); !approx(got, math.E-1, 1e-12) {
		t.Errorf("TanhSinh smooth = %.15g, want %.15g", got, math.E-1)
	}
	// Endpoint singularity: integral of 1/sqrt(x) over [0,1] = 2.
	sing := func(x float64) float64 { return 1 / math.Sqrt(x) }
	if got := TanhSinh(sing, 0, 1, 200); !approx(got, 2, 1e-6) {
		t.Errorf("TanhSinh singular = %.15g, want 2", got)
	}
	// Both endpoints singular: integral of 1/sqrt(1-x^2) over [-1,1] = pi.
	sing2 := func(x float64) float64 { return 1 / math.Sqrt(1-x*x) }
	if got := TanhSinh(sing2, -1, 1, 200); !approx(got, math.Pi, 1e-6) {
		t.Errorf("TanhSinh double-singular = %.15g, want %.15g", got, math.Pi)
	}
	if got := TanhSinhTol(fexp, 0, 1, 1e-12); !approx(got, math.E-1, 1e-10) {
		t.Errorf("TanhSinhTol = %.15g, want %.15g", got, math.E-1)
	}
	if got := DoubleExponential(fexp, 0, 1, 100); !approx(got, math.E-1, 1e-12) {
		t.Errorf("DoubleExponential = %.15g, want %.15g", got, math.E-1)
	}
}

func TestHighLevelAndRichardson(t *testing.T) {
	if got := Integrate(fpi, 0, 1); !approx(got, math.Pi, 1e-11) {
		t.Errorf("Integrate = %.15g, want %.15g", got, math.Pi)
	}
	// Trapezoid on x^4 has O(h^2) error; Richardson (order 2) removes it.
	f4 := func(x float64) float64 { return x * x * x * x }
	coarse := Trapezoid(f4, 0, 1, 50)
	fine := Trapezoid(f4, 0, 1, 100)
	ex := RichardsonExtrapolate(coarse, fine, 2)
	if !approx(ex, 0.2, 1e-6) {
		t.Errorf("RichardsonExtrapolate = %.15g, want 0.2", ex)
	}
}

func TestMultiDimensional(t *testing.T) {
	// integral of x*y over unit square = 1/4.
	xy := func(x, y float64) float64 { return x * y }
	if got := Trapezoid2D(xy, 0, 1, 0, 1, 200, 200); !approx(got, 0.25, 1e-5) {
		t.Errorf("Trapezoid2D = %.15g, want 0.25", got)
	}
	if got := Simpson2D(xy, 0, 1, 0, 1, 20, 20); !approx(got, 0.25, 1e-12) {
		t.Errorf("Simpson2D = %.15g, want 0.25", got)
	}
	if got := Midpoint2D(xy, 0, 1, 0, 1, 200, 200); !approx(got, 0.25, 1e-5) {
		t.Errorf("Midpoint2D = %.15g, want 0.25", got)
	}
	if got := GaussLegendre2D(xy, 0, 1, 0, 1, 4); !approx(got, 0.25, 1e-13) {
		t.Errorf("GaussLegendre2D = %.15g, want 0.25", got)
	}
	if got := Integrate2D(xy, 0, 1, 0, 1); !approx(got, 0.25, 1e-12) {
		t.Errorf("Integrate2D = %.15g, want 0.25", got)
	}

	// integral of x*y*z over unit cube = 1/8.
	xyz := func(x, y, z float64) float64 { return x * y * z }
	if got := Trapezoid3D(xyz, 0, 1, 0, 1, 0, 1, 60, 60, 60); !approx(got, 0.125, 1e-4) {
		t.Errorf("Trapezoid3D = %.15g, want 0.125", got)
	}
	if got := Simpson3D(xyz, 0, 1, 0, 1, 0, 1, 10, 10, 10); !approx(got, 0.125, 1e-12) {
		t.Errorf("Simpson3D = %.15g, want 0.125", got)
	}
	if got := GaussLegendre3D(xyz, 0, 1, 0, 1, 0, 1, 4); !approx(got, 0.125, 1e-13) {
		t.Errorf("GaussLegendre3D = %.15g, want 0.125", got)
	}

	// integral of x0*x1*x2*x3 over unit 4-cube = 1/16.
	prod := func(x []float64) float64 {
		p := 1.0
		for _, v := range x {
			p *= v
		}
		return p
	}
	low := []float64{0, 0, 0, 0}
	up := []float64{1, 1, 1, 1}
	if got := IntegrateND(prod, low, up, 4); !approx(got, 1.0/16.0, 1e-12) {
		t.Errorf("IntegrateND = %.15g, want %.15g", got, 1.0/16.0)
	}
}

func TestMonteCarloDeterministic(t *testing.T) {
	// Deterministic for a fixed seed: two calls agree exactly.
	a := MonteCarlo(f1, 0, 1, 100000, 42)
	b := MonteCarlo(f1, 0, 1, 100000, 42)
	if a != b {
		t.Errorf("MonteCarlo not deterministic: %v != %v", a, b)
	}
	if !approx(a, 1.0/3.0, 5e-3) {
		t.Errorf("MonteCarlo = %.6g, want ~0.3333", a)
	}
	xy := func(x, y float64) float64 { return x * y }
	if got := MonteCarlo2D(xy, 0, 1, 0, 1, 200000, 7); !approx(got, 0.25, 5e-3) {
		t.Errorf("MonteCarlo2D = %.6g, want ~0.25", got)
	}
	prod := func(x []float64) float64 { return x[0] * x[1] * x[2] }
	low := []float64{0, 0, 0}
	up := []float64{1, 1, 1}
	if got := MonteCarloND(prod, low, up, 300000, 11); !approx(got, 0.125, 5e-3) {
		t.Errorf("MonteCarloND = %.6g, want ~0.125", got)
	}
}

// BenchmarkAdaptiveGaussKronrod exercises the heaviest routine: globally
// adaptive Gauss-Kronrod on a peaked integrand.
func BenchmarkAdaptiveGaussKronrod(b *testing.B) {
	peak := func(x float64) float64 { return 1 / (1 + 1000*(x-0.4)*(x-0.4)) }
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = AdaptiveGaussKronrod(peak, 0, 1, 1e-10)
	}
}
