package quadrature

import (
	"fmt"
	"math"
	"testing"
)

const tol = 1e-9

func close(got, want, eps float64) bool {
	if math.IsNaN(got) || math.IsInf(got, 0) {
		return false
	}
	return math.Abs(got-want) <= eps
}

func assertClose(t *testing.T, name string, got, want, eps float64) {
	t.Helper()
	if !close(got, want, eps) {
		t.Errorf("%s = %.15g, want %.15g (tol %g, diff %g)", name, got, want, eps, math.Abs(got-want))
	}
}

// sumWeights is a small helper used by several weight-sum checks.
func sumWeights(w []float64) float64 {
	var s float64
	for _, x := range w {
		s += x
	}
	return s
}

func TestGaussLegendre(t *testing.T) {
	// Known 2-point nodes are +-1/sqrt(3), weights 1.
	n2, w2 := GaussLegendre(2)
	assertClose(t, "GL2 node-", n2[0], -1/math.Sqrt(3), 1e-14)
	assertClose(t, "GL2 node+", n2[1], 1/math.Sqrt(3), 1e-14)
	assertClose(t, "GL2 w0", w2[0], 1, 1e-14)
	assertClose(t, "GL2 w1", w2[1], 1, 1e-14)

	// 3-point nodes are 0, +-sqrt(3/5); weights 8/9, 5/9.
	n3, w3 := GaussLegendre(3)
	assertClose(t, "GL3 node0", n3[0], -math.Sqrt(3.0/5), 1e-14)
	assertClose(t, "GL3 node1", n3[1], 0, 1e-14)
	assertClose(t, "GL3 w mid", w3[1], 8.0/9, 1e-14)
	assertClose(t, "GL3 w side", w3[0], 5.0/9, 1e-14)

	// Weight sum equals interval length 2 for many n.
	for _, n := range []int{1, 4, 7, 12, 25} {
		_, w := GaussLegendre(n)
		assertClose(t, fmt.Sprintf("GL%d wsum", n), sumWeights(w), 2, 1e-12)
	}

	// Exactness: n-point rule integrates degree 2n-1 exactly.
	poly := func(x float64) float64 { return 3*x*x*x*x*x - 2*x*x*x + x - 4 }
	// exact integral over [-1,1] of the above = -8 (odd terms vanish, -4*2)
	assertClose(t, "GL3 poly", IntegrateGaussLegendre(poly, -1, 1, 3), -8, 1e-12)
	// e^x over [0,1] = e-1
	assertClose(t, "GL exp", IntegrateGaussLegendre(math.Exp, 0, 1, 8), math.E-1, 1e-12)
}

func TestGaussLegendreExactness(t *testing.T) {
	for n := 1; n <= 10; n++ {
		// x^(2n-1) integrates to 0 over [-1,1]; x^(2n-2) to 2/(2n-1).
		deg := 2*n - 2
		f := func(x float64) float64 { return math.Pow(x, float64(deg)) }
		want := 2.0 / float64(deg+1)
		got := IntegrateGaussLegendre(f, -1, 1, n)
		assertClose(t, fmt.Sprintf("GL%d x^%d", n, deg), got, want, 1e-10)
	}
}

func TestGaussChebyshev(t *testing.T) {
	// First kind: int 1/sqrt(1-x^2) = pi; int x^2/sqrt = pi/2.
	assertClose(t, "GC1 const", IntegrateGaussChebyshev1(func(x float64) float64 { return 1 }, 6), math.Pi, 1e-13)
	assertClose(t, "GC1 x^2", IntegrateGaussChebyshev1(func(x float64) float64 { return x * x }, 6), math.Pi/2, 1e-13)
	// x^4/sqrt(1-x^2) = 3pi/8
	assertClose(t, "GC1 x^4", IntegrateGaussChebyshev1(func(x float64) float64 { return x * x * x * x }, 6), 3*math.Pi/8, 1e-13)
	// Second kind: int sqrt(1-x^2) = pi/2; int x^2 sqrt = pi/8.
	assertClose(t, "GC2 const", IntegrateGaussChebyshev2(func(x float64) float64 { return 1 }, 6), math.Pi/2, 1e-13)
	assertClose(t, "GC2 x^2", IntegrateGaussChebyshev2(func(x float64) float64 { return x * x }, 6), math.Pi/8, 1e-13)
	// Third/fourth kind masses.
	_, w3 := GaussChebyshev3(7)
	assertClose(t, "GC3 mass", sumWeights(w3), math.Pi, 1e-12)
	_, w4 := GaussChebyshev4(7)
	assertClose(t, "GC4 mass", sumWeights(w4), math.Pi, 1e-12)
}

func TestGaussHermite(t *testing.T) {
	sp := math.Sqrt(math.Pi)
	// int e^{-x^2} = sqrt(pi); int x^2 e^{-x^2} = sqrt(pi)/2; int x^4 = 3 sqrt(pi)/4.
	assertClose(t, "GH const", IntegrateGaussHermite(func(x float64) float64 { return 1 }, 6), sp, 1e-12)
	assertClose(t, "GH x^2", IntegrateGaussHermite(func(x float64) float64 { return x * x }, 6), sp/2, 1e-12)
	assertClose(t, "GH x^4", IntegrateGaussHermite(func(x float64) float64 { return x * x * x * x }, 6), 3*sp/4, 1e-12)
	// Probabilists': int e^{-x^2/2} = sqrt(2pi); E[X^2]=1 for N(0,1).
	assertClose(t, "GHprob mass", IntegrateGaussHermiteProb(func(x float64) float64 { return 1 }, 6), math.Sqrt(2*math.Pi), 1e-12)
	assertClose(t, "E[X^2]", ExpectationGaussHermite(func(x float64) float64 { return x * x }, 6), 1, 1e-12)
	// E[X^4] = 3 for standard normal.
	assertClose(t, "E[X^4]", ExpectationGaussHermite(func(x float64) float64 { return x * x * x * x }, 6), 3, 1e-12)
	// Folded weights recover an un-weighted integral: int e^{-x^2} cos x dx
	// over the real line equals sqrt(pi) * e^{-1/4}.
	got := IntegrateGaussHermiteFunc(func(x float64) float64 { return math.Exp(-x*x) * math.Cos(x) }, 20)
	assertClose(t, "GH folded", got, math.Sqrt(math.Pi)*math.Exp(-0.25), 1e-9)
}

func TestGaussLaguerre(t *testing.T) {
	// int x^k e^{-x} = k! over [0,inf).
	facts := []float64{1, 1, 2, 6, 24, 120}
	for k, want := range facts {
		f := func(x float64) float64 { return math.Pow(x, float64(k)) }
		got := IntegrateGaussLaguerre(f, 6)
		assertClose(t, fmt.Sprintf("Lag x^%d", k), got, want, 1e-9)
	}
	// Generalized: int x^a e^{-x} = Gamma(a+1); with a=1.5, k=0 -> Gamma(2.5).
	got := IntegrateGaussLaguerreGen(func(x float64) float64 { return 1 }, 6, 1.5)
	assertClose(t, "genLag mass", got, math.Gamma(2.5), 1e-10)
	// Exponential expectation E[X]=1.
	assertClose(t, "Exp E[X]", ExpectationExponential(func(x float64) float64 { return x }, 8), 1, 1e-9)
}

func TestGaussJacobi(t *testing.T) {
	// Jacobi with a=b=0 reduces to Legendre; mass = 2.
	_, w := GaussJacobi(5, 0, 0)
	assertClose(t, "Jac00 mass", sumWeights(w), 2, 1e-12)
	// Mass for general a,b equals 2^{a+b+1} B(a+1,b+1).
	a, b := 1.5, 0.5
	_, wj := GaussJacobi(6, a, b)
	wantMass := math.Pow(2, a+b+1) * math.Gamma(a+1) * math.Gamma(b+1) / math.Gamma(a+b+2)
	assertClose(t, "Jac mass", sumWeights(wj), wantMass, 1e-11)
	// Gegenbauer lambda=1 -> weight (1-x^2), mass = pi/2.
	_, wg := GaussGegenbauer(6, 1)
	assertClose(t, "Geg mass", sumWeights(wg), math.Pi/2, 1e-11)
}

func TestGaussLobattoRadau(t *testing.T) {
	// Lobatto n=3: nodes -1,0,1; weights 1/3,4/3,1/3.
	ln, lw := GaussLobatto(3)
	assertClose(t, "Lob node-", ln[0], -1, 1e-13)
	assertClose(t, "Lob node0", ln[1], 0, 1e-13)
	assertClose(t, "Lob node+", ln[2], 1, 1e-13)
	assertClose(t, "Lob w0", lw[0], 1.0/3, 1e-13)
	assertClose(t, "Lob wm", lw[1], 4.0/3, 1e-13)
	// Exactness of n-point Lobatto: degree 2n-3.
	// n=4 -> degree 5 exact. int x^4 over [-1,1] = 2/5.
	assertClose(t, "Lob x^4", IntegrateGaussLobatto(func(x float64) float64 { return x * x * x * x }, -1, 1, 4), 0.4, 1e-12)

	// Radau n=2: nodes -1,1/3; weights 1/2,3/2.
	rn, rw := GaussRadau(2)
	assertClose(t, "Rad node-", rn[0], -1, 1e-13)
	assertClose(t, "Rad node+", rn[1], 1.0/3, 1e-13)
	assertClose(t, "Rad w0", rw[0], 0.5, 1e-13)
	assertClose(t, "Rad w1", rw[1], 1.5, 1e-13)
	// Radau n=3 known nodes -1, (1-sqrt6)/5, (1+sqrt6)/5.
	rn3, _ := GaussRadau(3)
	assertClose(t, "Rad3 node1", rn3[1], (1-math.Sqrt(6))/5, 1e-12)
	assertClose(t, "Rad3 node2", rn3[2], (1+math.Sqrt(6))/5, 1e-12)
	// Right Radau fixes +1.
	rr, _ := GaussRadauRight(3)
	assertClose(t, "RadR last", rr[2], 1, 1e-13)
	// Weight sums equal 2.
	assertClose(t, "Lob wsum", sumWeights(lw), 2, 1e-13)
	assertClose(t, "Rad wsum", sumWeights(rw), 2, 1e-13)
}

func TestNewtonCotes(t *testing.T) {
	cube := func(x float64) float64 { return x * x * x }
	// Simpson and beyond integrate cubics exactly: int_0^2 x^3 = 4.
	assertClose(t, "Simpson cube", Simpson(cube, 0, 2), 4, 1e-12)
	assertClose(t, "Simpson38 cube", Simpson38(cube, 0, 2), 4, 1e-12)
	assertClose(t, "Boole cube", Boole(cube, 0, 2), 4, 1e-12)
	// Boole integrates degree-5 exactly: int_0^1 x^5 = 1/6.
	assertClose(t, "Boole x^5", Boole(func(x float64) float64 { return math.Pow(x, 5) }, 0, 1), 1.0/6, 1e-12)
	// Trapezoid/midpoint exact for linear.
	lin := func(x float64) float64 { return 2*x + 1 }
	assertClose(t, "Trap lin", Trapezoid(lin, 1, 3), 10, 1e-12)
	assertClose(t, "Mid lin", Midpoint(lin, 1, 3), 10, 1e-12)

	// Composite rules converge on int_0^pi sin = 2.
	assertClose(t, "compSimpson", CompositeSimpson(math.Sin, 0, math.Pi, 64), 2, 1e-7)
	assertClose(t, "compBoole", CompositeBoole(math.Sin, 0, math.Pi, 32), 2, 1e-8)
	assertClose(t, "compTrap", CompositeTrapezoid(math.Sin, 0, math.Pi, 400), 2, 1e-4)
	assertClose(t, "compMid", CompositeMidpoint(math.Sin, 0, math.Pi, 400), 2, 1e-4)
	assertClose(t, "compS38", CompositeSimpson38(math.Sin, 0, math.Pi, 60), 2, 1e-6)

	// Exact Newton-Cotes weights: closed n=2 is 1/3,4/3,1/3.
	w := NewtonCotesClosed(2)
	assertClose(t, "NC2 w0", w[0], 1.0/3, 1e-14)
	assertClose(t, "NC2 w1", w[1], 4.0/3, 1e-14)
	assertClose(t, "NC2 w2", w[2], 1.0/3, 1e-14)
	// Closed n=4 (Boole) weights sum to 4.
	assertClose(t, "NC4 sum", sumWeights(NewtonCotesClosed(4)), 4, 1e-13)
	// Open n=2 weights are 3/2,3/2.
	ow := NewtonCotesOpen(2)
	assertClose(t, "NCopen2 w0", ow[0], 1.5, 1e-14)
	// Single-panel NC of order 5 integrates degree 5 exactly.
	assertClose(t, "NC5 poly", IntegrateNewtonCotesClosed(func(x float64) float64 { return math.Pow(x, 5) }, 0, 1, 5), 1.0/6, 1e-12)
}

func TestSampleIntegration(t *testing.T) {
	// Samples of x^2 at 0..4, h=1: exact integral 64/3.
	ys := []float64{0, 1, 4, 9, 16}
	assertClose(t, "SimpsonSamples", SimpsonSamples(ys, 1), 64.0/3, 1e-12)
	assertClose(t, "BooleSamples", BooleSamples(ys, 1), 64.0/3, 1e-12)
	assertClose(t, "IntegrateSamples", IntegrateSamples(ys, 1), 64.0/3, 1e-12)
	assertClose(t, "TrapSamples", TrapezoidSamples(ys, 1), 22, 1e-12) // trapezoid overestimates

	// Non-uniform trapezoid: linear function integrates exactly.
	xs := []float64{0, 0.3, 1.1, 2.0}
	yv := make([]float64, len(xs))
	for i, x := range xs {
		yv[i] = 3*x + 1
	}
	assertClose(t, "TrapNonUniform", TrapezoidNonUniform(xs, yv), 3*2.0*2.0/2+2.0, 1e-12)

	// Cumulative trapezoid: last entry equals total.
	cum := CumulativeTrapezoid(xs, yv)
	assertClose(t, "cumTrap last", cum[len(cum)-1], TrapezoidNonUniform(xs, yv), 1e-14)

	// Midpoint samples of a constant.
	assertClose(t, "MidpointSamples", MidpointSamples([]float64{2, 2, 2}, 0.5), 3, 1e-14)
}

func TestRomberg(t *testing.T) {
	assertClose(t, "Romberg sin", Romberg(math.Sin, 0, math.Pi, 1e-12, 20).Value, 2, 1e-11)
	assertClose(t, "Romberg exp", Romberg(math.Exp, 0, 1, 1e-12, 20).Value, math.E-1, 1e-11)
	assertClose(t, "RombergValue", RombergValue(math.Cos, 0, math.Pi/2, 12), 1, 1e-10)
	// Richardson extrapolation of the trapezoid rule (order 2) on sin.
	c := CompositeTrapezoid(math.Sin, 0, math.Pi, 4)
	f := CompositeTrapezoid(math.Sin, 0, math.Pi, 8)
	rich := RichardsonExtrapolate(c, f, 2, 2)
	assertClose(t, "Richardson", rich, CompositeSimpson(math.Sin, 0, math.Pi, 8), 1e-12)
}

func TestAdaptive(t *testing.T) {
	// Smooth integrand.
	assertClose(t, "adaptSimpson", AdaptiveSimpson(math.Sin, 0, math.Pi, 1e-10), 2, 1e-9)
	assertClose(t, "adaptTrap", AdaptiveTrapezoid(math.Sin, 0, math.Pi, 1e-8), 2, 1e-6)
	assertClose(t, "adaptLobatto", AdaptiveLobatto(math.Sin, 0, math.Pi, 1e-10), 2, 1e-9)
	// Peaked integrand: 1/(1+25x^2) over [-1,1] = (2/5) arctan(5).
	runge := func(x float64) float64 { return 1 / (1 + 25*x*x) }
	want := 2.0 / 5 * math.Atan(5)
	assertClose(t, "adaptSimpson runge", AdaptiveSimpson(runge, -1, 1, 1e-10), want, 1e-8)
	res := AdaptiveSimpsonResult(runge, -1, 1, 1e-11)
	if !res.Success || res.Evals == 0 {
		t.Errorf("AdaptiveSimpsonResult reported no work: %+v", res)
	}
}

func TestGaussKronrod(t *testing.T) {
	// Exact for polynomials up to degree 2*15-1: use degree 10.
	f := func(x float64) float64 { return math.Pow(x, 10) }
	v, e := GaussKronrod15Eval(f, 0, 1)
	assertClose(t, "K15 x^10", v, 1.0/11, 1e-12)
	if e > 1e-9 {
		t.Errorf("K15 error estimate too large: %g", e)
	}
	// Adaptive on a harder integrand: int_0^1 sqrt(x) = 2/3.
	r := AdaptiveGaussKronrod(math.Sqrt, 0, 1, 1e-10)
	assertClose(t, "AGK sqrt", r.Value, 2.0/3, 1e-8)
	// int_0^pi sin = 2.
	r2 := AdaptiveGaussKronrod(math.Sin, 0, math.Pi, 1e-12)
	assertClose(t, "AGK sin", r2.Value, 2, 1e-11)
	// Embedded Gauss rule weights sum to 2.
	_, _, gw := GaussKronrod15()
	assertClose(t, "GK gauss wsum", sumWeights(gw), 2, 1e-12)
}

func TestTanhSinh(t *testing.T) {
	// Smooth: int_0^1 e^x = e-1.
	assertClose(t, "TS exp", TanhSinh(math.Exp, 0, 1, 1e-12), math.E-1, 1e-11)
	// Endpoint singularity: int_0^1 1/sqrt(x) = 2.
	assertClose(t, "TS 1/sqrt", TanhSinh(func(x float64) float64 { return 1 / math.Sqrt(x) }, 0, 1, 1e-12), 2, 1e-7)
	// int_0^1 ln x = -1.
	assertClose(t, "TS log", TanhSinh(math.Log, 0, 1, 1e-12), -1, 1e-9)
	// int_0^1 x^{-1/2}(1-x)^{-1/2} = pi (both endpoints singular).
	beta := func(x float64) float64 { return 1 / math.Sqrt(x*(1-x)) }
	assertClose(t, "TS beta", TanhSinh(beta, 0, 1, 1e-12), math.Pi, 1e-6)
	// Semi-infinite: int_0^inf e^{-x} = 1; int_0^inf e^{-x} x = 1.
	assertClose(t, "ExpSinh", ExpSinh(func(x float64) float64 { return math.Exp(-x) }, 0, 1e-11), 1, 1e-9)
	assertClose(t, "ExpSinh x", ExpSinh(func(x float64) float64 { return x * math.Exp(-x) }, 0, 1e-11), 1, 1e-8)
	// Whole line: int e^{-x^2} = sqrt(pi).
	assertClose(t, "SinhSinh", SinhSinh(func(x float64) float64 { return math.Exp(-x * x) }, 1e-11), math.Sqrt(math.Pi), 1e-9)
	// 1/(1+x^2) over the whole line = pi.
	assertClose(t, "IntInfinite", IntegrateInfinite(func(x float64) float64 { return 1 / (1 + x*x) }, 1e-10), math.Pi, 1e-6)
	// Left semi-infinite: int_{-inf}^0 e^{x} = 1.
	assertClose(t, "ExpSinhLeft", ExpSinhLeft(func(x float64) float64 { return math.Exp(x) }, 0, 1e-11), 1, 1e-9)
	// DoubleExponential alias matches.
	assertClose(t, "DE alias", DoubleExponential(math.Exp, 0, 1, 1e-12), TanhSinh(math.Exp, 0, 1, 1e-12), 1e-15)
}

func TestClenshawCurtisFejer(t *testing.T) {
	// Weight sums equal 2 on [-1,1].
	_, ccw := ClenshawCurtisWeights(8)
	assertClose(t, "CC wsum", sumWeights(ccw), 2, 1e-13)
	_, f1w := Fejer1Weights(8)
	assertClose(t, "F1 wsum", sumWeights(f1w), 2, 1e-13)
	_, f2w := Fejer2Weights(8)
	assertClose(t, "F2 wsum", sumWeights(f2w), 2, 1e-13)

	// Polynomial exactness / convergence.
	x4 := func(x float64) float64 { return x * x * x * x }
	assertClose(t, "CC x^4", ClenshawCurtis(x4, -1, 1, 8), 0.4, 1e-12)
	assertClose(t, "F1 x^4", IntegrateFejer1(x4, -1, 1, 8), 0.4, 1e-12)
	assertClose(t, "F2 x^4", IntegrateFejer2(x4, -1, 1, 8), 0.4, 1e-12)
	// Smooth transcendental.
	assertClose(t, "CC exp", ClenshawCurtis(math.Exp, 0, 1, 24), math.E-1, 1e-12)
	assertClose(t, "CC cos", ClenshawCurtis(math.Cos, 0, math.Pi, 32), 0, 1e-12)
}

func TestOrthogonalPolynomials(t *testing.T) {
	// ChebyshevT(n, cos t) = cos(n t).
	for _, th := range []float64{0.3, 1.0, 2.5} {
		assertClose(t, "T5", ChebyshevT(5, math.Cos(th)), math.Cos(5*th), 1e-12)
	}
	// Legendre matches JacobiP(n,0,0,.) and known P2, P3.
	assertClose(t, "P2", LegendreP(2, 0.7), 0.5*(3*0.49-1), 1e-13)
	assertClose(t, "P3 vs Jacobi", JacobiP(3, 0, 0, 0.4), LegendreP(3, 0.4), 1e-13)
	// Hermite physicists H3(x) = 8x^3-12x.
	assertClose(t, "H3", HermiteH(3, 1.3), 8*math.Pow(1.3, 3)-12*1.3, 1e-12)
	// Probabilists He3(x) = x^3-3x.
	assertClose(t, "He3", HermiteHe(3, 1.3), math.Pow(1.3, 3)-3*1.3, 1e-12)
	// Laguerre L2(x) = 1 - 2x + x^2/2.
	assertClose(t, "L2", LaguerreL(2, 0.9), 1-2*0.9+0.9*0.9/2, 1e-13)
	// ChebyshevU(n, cos t) sin t = sin((n+1) t).
	th := 0.7
	assertClose(t, "U4", ChebyshevU(4, math.Cos(th))*math.Sin(th), math.Sin(5*th), 1e-12)
}

func TestMultiDimensional(t *testing.T) {
	// int_[0,1]^2 x*y = 1/4.
	xy := func(x, y float64) float64 { return x * y }
	assertClose(t, "GL2D", DoubleGaussLegendre(xy, 0, 1, 0, 1, 4), 0.25, 1e-12)
	assertClose(t, "Simpson2D", DoubleSimpson(xy, 0, 1, 0, 1, 8, 8), 0.25, 1e-12)
	assertClose(t, "Boole2D", DoubleBoole(xy, 0, 1, 0, 1, 8, 8), 0.25, 1e-12)
	assertClose(t, "Trap2D", DoubleTrapezoid(xy, 0, 1, 0, 1, 200, 200), 0.25, 1e-4)
	assertClose(t, "Mid2D", DoubleMidpoint(xy, 0, 1, 0, 1, 200, 200), 0.25, 1e-4)
	assertClose(t, "DoubleAvg", DoubleAverage(func(x, y float64) float64 { return x + y }, 0, 2, 0, 2, 4), 2, 1e-12)

	// int_[0,1]^2 x^2 y^2 = 1/9.
	x2y2 := func(x, y float64) float64 { return x * x * y * y }
	assertClose(t, "GL2D x2y2", DoubleGaussLegendre(x2y2, 0, 1, 0, 1, 5), 1.0/9, 1e-12)

	// Triple integral of 1 over unit cube = 1.
	assertClose(t, "GL3D vol", TripleGaussLegendre(func(x, y, z float64) float64 { return 1 }, 0, 1, 0, 1, 0, 1, 3), 1, 1e-12)
	assertClose(t, "Simpson3D", TripleSimpson(func(x, y, z float64) float64 { return x * y * z }, 0, 1, 0, 1, 0, 1, 4, 4, 4), 0.125, 1e-12)
	assertClose(t, "Trap3D", TripleTrapezoid(func(x, y, z float64) float64 { return 1 }, 0, 1, 0, 1, 0, 1, 4, 4, 4), 1, 1e-12)
	assertClose(t, "Mid3D", TripleMidpoint(func(x, y, z float64) float64 { return 1 }, 0, 1, 0, 1, 0, 1, 4, 4, 4), 1, 1e-12)

	// Tensor Gauss over a box via FuncN.
	got := IntegrateBoxGaussLegendre(func(p []float64) float64 { return p[0] * p[1] * p[2] }, []float64{0, 0, 0}, []float64{1, 1, 1}, 3)
	assertClose(t, "boxGL", got, 0.125, 1e-12)

	// Product rule directly.
	rx := GaussLegendreRule(4).Scale(0, 1)
	ry := GaussLegendreRule(4).Scale(0, 1)
	assertClose(t, "product2", IntegrateProduct2(xy, rx, ry), 0.25, 1e-12)
}

func TestMonteCarlo(t *testing.T) {
	xy := func(p []float64) float64 { return p[0] * p[1] }
	mc := MonteCarloBox(xy, []float64{0, 0}, []float64{1, 1}, 200000, 12345)
	assertClose(t, "MC xy", mc.Value, 0.25, 5e-3)
	if mc.AbsErr <= 0 {
		t.Errorf("MonteCarloBox stderr should be positive, got %g", mc.AbsErr)
	}
	// Quasi Monte Carlo converges more tightly for smooth integrands.
	qmc := QuasiMonteCarloBox(xy, []float64{0, 0}, []float64{1, 1}, 20000)
	assertClose(t, "QMC xy", qmc.Value, 0.25, 2e-3)
	// Stratified.
	sm := StratifiedMonteCarloBox(xy, []float64{0, 0}, []float64{1, 1}, 100, 7)
	assertClose(t, "stratified", sm.Value, 0.25, 5e-3)
	// MonteCarlo2/3 convenience wrappers.
	mc2 := MonteCarlo2(func(x, y float64) float64 { return x + y }, 0, 1, 0, 1, 100000, 99)
	assertClose(t, "MC2", mc2.Value, 1, 5e-3)
	mc3 := MonteCarlo3(func(x, y, z float64) float64 { return 1 }, 0, 1, 0, 1, 0, 1, 50000, 3)
	assertClose(t, "MC3 vol", mc3.Value, 1, 1e-9)
}

func TestQuasiRandom(t *testing.T) {
	// Van der Corput base 2: 1->0.5, 2->0.25, 3->0.75.
	assertClose(t, "vdc1", VanDerCorput(1, 2), 0.5, 1e-15)
	assertClose(t, "vdc2", VanDerCorput(2, 2), 0.25, 1e-15)
	assertClose(t, "vdc3", VanDerCorput(3, 2), 0.75, 1e-15)
	// Base 3: 1 -> 1/3.
	assertClose(t, "vdc base3", VanDerCorput(1, 3), 1.0/3, 1e-15)
	// Halton points lie in the unit cube.
	seq := HaltonSequence(50, 3)
	for _, p := range seq {
		for _, c := range p {
			if c < 0 || c >= 1 {
				t.Errorf("Halton coordinate out of range: %g", c)
			}
		}
	}
}

func TestRuleMethods(t *testing.T) {
	r := GaussLegendreRule(5)
	assertClose(t, "rule wsum", r.WeightSum(), 2, 1e-13)
	if r.Len() != 5 {
		t.Errorf("Len = %d, want 5", r.Len())
	}
	// Scale to [0,1] and integrate x^2 -> 1/3.
	scaled := r.Scale(0, 1)
	assertClose(t, "scaled x^2", scaled.Integrate(func(x float64) float64 { return x * x }), 1.0/3, 1e-12)
	assertClose(t, "IntegrateOn", r.IntegrateOn(func(x float64) float64 { return x * x }, 0, 1), 1.0/3, 1e-12)
	// Reversed keeps the integral.
	assertClose(t, "reversed", r.Reversed().Integrate(func(x float64) float64 { return x * x }), r.Integrate(func(x float64) float64 { return x * x }), 1e-15)
	// Transform: substitution u=(x+1)/2 maps [-1,1] rule to [0,1].
	tr := r.Transform(func(x float64) float64 { return (x + 1) / 2 }, func(x float64) float64 { return 0.5 })
	assertClose(t, "transform x^2", tr.Integrate(func(x float64) float64 { return x * x }), 1.0/3, 1e-12)
	// Rescale helper.
	rn, rw := Rescale(r.Nodes, r.Weights, 0, 1)
	assertClose(t, "rescale wsum", sumWeights(rw), 1, 1e-13)
	if len(rn) != 5 {
		t.Errorf("Rescale returned %d nodes", len(rn))
	}
}

func TestGeometryTypes(t *testing.T) {
	iv := Interval{A: 1, B: 4}
	assertClose(t, "iv len", iv.Length(), 3, 0)
	assertClose(t, "iv mid", iv.Mid(), 2.5, 0)
	if !iv.Contains(2) || iv.Contains(5) {
		t.Errorf("Interval.Contains wrong")
	}
	b := NewBox([]float64{0, 0}, []float64{2, 3})
	if b.Dim() != 2 {
		t.Errorf("Box.Dim = %d", b.Dim())
	}
	assertClose(t, "box vol", b.Volume(), 6, 0)
	if !b.Contains([]float64{1, 1}) || b.Contains([]float64{3, 1}) {
		t.Errorf("Box.Contains wrong")
	}
}

func TestGolubWelschDirect(t *testing.T) {
	// Feeding the Legendre recurrence to GolubWelsch reproduces GaussLegendre.
	alpha, beta := RecurrenceLegendre(6)
	n1, w1 := GolubWelsch(alpha, beta)
	n2, w2 := GaussLegendre(6)
	for i := range n1 {
		assertClose(t, "GW node", n1[i], n2[i], 1e-13)
		assertClose(t, "GW weight", w1[i], w2[i], 1e-13)
	}
}

// ExampleAdaptiveSimpson demonstrates integrating a smooth function.
func ExampleAdaptiveSimpson() {
	// Integrate sin over [0, pi]; the exact value is 2.
	got := AdaptiveSimpson(math.Sin, 0, math.Pi, 1e-10)
	fmt.Printf("%.10f\n", got)
	// Output: 2.0000000000
}

// ExampleGaussLegendre shows generating nodes and weights.
func ExampleGaussLegendre() {
	nodes, weights := GaussLegendre(2)
	fmt.Printf("nodes: %.6f, %.6f\n", nodes[0], nodes[1])
	fmt.Printf("weights: %.6f, %.6f\n", weights[0], weights[1])
	// Output:
	// nodes: -0.577350, 0.577350
	// weights: 1.000000, 1.000000
}

// ExampleTanhSinh integrates a function with an endpoint singularity.
func ExampleTanhSinh() {
	// Integrate 1/sqrt(x) over [0, 1]; the exact value is 2.
	got := TanhSinh(func(x float64) float64 { return 1 / math.Sqrt(x) }, 0, 1, 1e-9)
	fmt.Printf("%.5f\n", got)
	// Output: 2.00000
}
