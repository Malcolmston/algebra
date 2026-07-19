package approxtheory

import (
	"fmt"
	"math"
	"testing"
)

// close reports whether a and b agree to within an absolute tolerance.
func close(a, b, tol float64) bool {
	if math.IsNaN(a) || math.IsNaN(b) {
		return false
	}
	return math.Abs(a-b) <= tol
}

// sliceClose reports whether two slices agree elementwise within tol.
func sliceClose(a, b []float64, tol float64) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if !close(a[i], b[i], tol) {
			return false
		}
	}
	return true
}

func TestPolyval(t *testing.T) {
	tests := []struct {
		coeffs  []float64
		x, want float64
	}{
		{[]float64{1, 2, 3}, 2, 1 + 2*2 + 3*4}, // 17
		{[]float64{0}, 5, 0},
		{[]float64{-1, 0, 0, 4}, 1, 3},
		{nil, 3, 0},
	}
	for _, tt := range tests {
		if got := Polyval(tt.coeffs, tt.x); !close(got, tt.want, 1e-12) {
			t.Errorf("Polyval(%v,%v)=%v want %v", tt.coeffs, tt.x, got, tt.want)
		}
	}
}

func TestPolyArithmetic(t *testing.T) {
	a := []float64{1, 2, 3}
	b := []float64{0, 1}
	if got := PolyAdd(a, b); !sliceClose(got, []float64{1, 3, 3}, 1e-12) {
		t.Errorf("PolyAdd=%v", got)
	}
	if got := PolySub(a, b); !sliceClose(got, []float64{1, 1, 3}, 1e-12) {
		t.Errorf("PolySub=%v", got)
	}
	if got := PolyMul([]float64{1, 1}, []float64{-1, 1}); !sliceClose(got, []float64{-1, 0, 1}, 1e-12) {
		t.Errorf("PolyMul=%v", got)
	}
	if got := PolyDeriv([]float64{5, 4, 3}); !sliceClose(got, []float64{4, 6}, 1e-12) {
		t.Errorf("PolyDeriv=%v", got)
	}
	if got := PolyInt([]float64{4, 6}, 5); !sliceClose(got, []float64{5, 4, 3}, 1e-12) {
		t.Errorf("PolyInt=%v", got)
	}
	if got := PolyFromRoots([]float64{1, -1}); !sliceClose(got, []float64{-1, 0, 1}, 1e-12) {
		t.Errorf("PolyFromRoots=%v", got)
	}
	if PolyDegree([]float64{1, 2, 0, 0}) != 1 {
		t.Errorf("PolyDegree wrong")
	}
	if !close(PolyLeadingCoeff([]float64{1, 2, 3, 0}), 3, 0) {
		t.Errorf("PolyLeadingCoeff wrong")
	}
	if s := PolyString([]float64{1, -2, 3}); s != "1 - 2*x + 3*x^2" {
		t.Errorf("PolyString=%q", s)
	}
	if !PolyEqual([]float64{1, 2}, []float64{1, 2, 0}, 1e-12) {
		t.Errorf("PolyEqual wrong")
	}
}

func TestChebT(t *testing.T) {
	tests := []struct {
		n    int
		x    float64
		want float64
	}{
		{0, 0.5, 1},
		{1, 0.5, 0.5},
		{2, 0.5, 2*0.25 - 1},      // -0.5
		{3, 0.5, 4*0.125 - 3*0.5}, // -1
		{5, 1, 1},
		{4, -1, 1},
	}
	for _, tt := range tests {
		if got := ChebT(tt.n, tt.x); !close(got, tt.want, 1e-12) {
			t.Errorf("ChebT(%d,%v)=%v want %v", tt.n, tt.x, got, tt.want)
		}
	}
	// T_n(cos t) = cos(n t)
	for _, x := range []float64{-0.9, -0.3, 0.2, 0.75} {
		th := math.Acos(x)
		for n := 0; n <= 6; n++ {
			if !close(ChebT(n, x), math.Cos(float64(n)*th), 1e-12) {
				t.Errorf("ChebT identity failed n=%d x=%v", n, x)
			}
		}
	}
}

func TestChebFitAndEval(t *testing.T) {
	fns := []struct {
		name string
		f    func(float64) float64
		a, b float64
		tol  float64
	}{
		{"exp", math.Exp, -1, 1, 1e-12},
		{"sin", math.Sin, -3, 3, 1e-10},
		{"runge", RungeFunction, -1, 1, 1e-3},
	}
	for _, fn := range fns {
		s := ChebFit(fn.f, 40, fn.a, fn.b)
		if err := MaxError(fn.f, s.Eval, fn.a, fn.b); err > fn.tol {
			t.Errorf("%s: ChebFit error %g > %g", fn.name, err, fn.tol)
		}
	}
	// exact reproduction of a quadratic
	q := func(x float64) float64 { return 3 + 2*x + x*x }
	s := ChebFit(q, 2, -2, 5)
	if !close(s.Eval(4), q(4), 1e-9) {
		t.Errorf("quadratic reproduction failed: %v", s.Eval(4))
	}
	// coefficient of x^2 in [-1,1] Chebyshev basis: x^2 = (T0+T2)/2
	sc := ChebFit(func(x float64) float64 { return x * x }, 2, -1, 1)
	if !sliceClose(sc.Coeffs, []float64{0.5, 0, 0.5}, 1e-12) {
		t.Errorf("x^2 cheb coeffs=%v", sc.Coeffs)
	}
}

func TestChebCalculus(t *testing.T) {
	s := ChebFit(math.Exp, 25, -1, 1)
	d := s.Derivative()
	if err := MaxError(math.Exp, d.Eval, -1, 1); err > 1e-10 {
		t.Errorf("derivative error %g", err)
	}
	// integral of exp over [-1,1] = e - 1/e
	if !close(s.Integral(), math.Exp(1)-math.Exp(-1), 1e-12) {
		t.Errorf("Integral=%v", s.Integral())
	}
	// definite integral 0..1 = e - 1
	if !close(s.DefiniteIntegral(0, 1), math.Exp(1)-1, 1e-12) {
		t.Errorf("DefiniteIntegral=%v", s.DefiniteIntegral(0, 1))
	}
	// antiderivative pinned at A
	F := s.Antiderivative()
	if !close(F.Eval(-1), 0, 1e-12) {
		t.Errorf("Antiderivative(A)=%v", F.Eval(-1))
	}
	if !close(F.Eval(0.5)-F.Eval(-0.5), s.DefiniteIntegral(-0.5, 0.5), 1e-12) {
		t.Errorf("antiderivative/definite mismatch")
	}
}

func TestChebRoots(t *testing.T) {
	// cos on [0,pi] has a single root at pi/2
	cs := ChebFit(math.Cos, 20, 0, math.Pi)
	roots := cs.Roots()
	if len(roots) != 1 || !close(roots[0], math.Pi/2, 1e-9) {
		t.Errorf("cos roots=%v", roots)
	}
	// polynomial (x-0.2)(x+0.5)(x-0.9) on [-1,1]
	p := func(x float64) float64 { return (x - 0.2) * (x + 0.5) * (x - 0.9) }
	ps := ChebFit(p, 3, -1, 1)
	got := ps.Roots()
	want := []float64{-0.5, 0.2, 0.9}
	if !sliceClose(got, want, 1e-8) {
		t.Errorf("cubic roots=%v want %v", got, want)
	}
}

func TestChebProductAndConvert(t *testing.T) {
	a := NewChebSeries([]float64{0, 1}, -1, 1) // x
	b := NewChebSeries([]float64{0, 1}, -1, 1) // x
	prod := ChebProduct(a, b)                  // x^2 = (T0+T2)/2
	if !close(prod.Eval(0.7), 0.49, 1e-12) {
		t.Errorf("ChebProduct eval=%v", prod.Eval(0.7))
	}
	// round trip monomial<->chebyshev
	mono := []float64{1, -2, 0, 3}
	cheb := MonomialToChebyshev(mono)
	back := ChebyshevToMonomial(cheb)
	if !PolyEqual(mono, back, 1e-10) {
		t.Errorf("mono round trip: %v -> %v", mono, back)
	}
	// ToMonomial of a fitted quadratic on a shifted interval
	q := func(x float64) float64 { return 3 + 2*x + x*x }
	s := ChebFit(q, 2, -2, 5)
	if !sliceClose(s.ToMonomial(), []float64{3, 2, 1}, 1e-9) {
		t.Errorf("ToMonomial=%v", s.ToMonomial())
	}
}

func TestRemez(t *testing.T) {
	tests := []struct {
		name    string
		f       func(float64) float64
		deg     int
		want    []float64
		wantErr float64
		tol     float64
	}{
		// best degree-1 approx of x^2 on [-1,1] is 1/2, error 1/2
		{"x^2", func(x float64) float64 { return x * x }, 1, []float64{0.5, 0}, 0.5, 1e-6},
		// best degree-2 approx of x^3 is (3/4)x, error 1/4
		{"x^3", func(x float64) float64 { return x * x * x }, 2, []float64{0, 0.75, 0}, 0.25, 1e-6},
	}
	for _, tt := range tests {
		r, err := RemezPoly(tt.f, tt.deg, -1, 1, 100, 1e-4)
		if err != nil {
			t.Fatalf("%s: %v", tt.name, err)
		}
		if !sliceClose(r.Coeffs, tt.want, tt.tol) {
			t.Errorf("%s coeffs=%v want %v", tt.name, r.Coeffs, tt.want)
		}
		if !close(r.Error, tt.wantErr, tt.tol) {
			t.Errorf("%s error=%v want %v", tt.name, r.Error, tt.wantErr)
		}
		if !r.Converged {
			t.Errorf("%s did not converge", tt.name)
		}
	}
	// Remez must beat (or match) the Chebyshev interpolant in the sup norm.
	f := math.Exp
	r, _ := RemezPoly(f, 5, -1, 1, 100, 1e-4)
	remezErr := MinimaxError(f, r.Coeffs, -1, 1)
	cheb := ChebFit(f, 5, -1, 1)
	chebErr := MaxError(f, cheb.Eval, -1, 1)
	if remezErr > chebErr*1.001 {
		t.Errorf("Remez error %g not <= Chebyshev error %g", remezErr, chebErr)
	}
	// known minimax error for exp on [-1,1], degree 5 ~ 4.52e-5
	if !close(r.Error, 4.52e-5, 5e-7) {
		t.Errorf("exp deg5 minimax error=%v", r.Error)
	}
}

func TestPade(t *testing.T) {
	// [1/1] Pade of exp = (1+x/2)/(1-x/2)
	p, err := PadeApprox(TaylorExp(2), 1, 1)
	if err != nil {
		t.Fatal(err)
	}
	if !sliceClose(p.Num, []float64{1, 0.5}, 1e-12) || !sliceClose(p.Den, []float64{1, -0.5}, 1e-12) {
		t.Errorf("pade exp [1/1] num=%v den=%v", p.Num, p.Den)
	}
	// [2/2] Pade of exp is far better than the degree-4 Taylor polynomial.
	p22 := PadeExp(2, 2)
	taylorErr := math.Abs(math.Exp(1) - Polyval(TaylorExp(4), 1))
	padeErr := math.Abs(math.Exp(1) - p22.Eval(1))
	if padeErr >= taylorErr {
		t.Errorf("pade error %g not < taylor error %g", padeErr, taylorErr)
	}
	// Pade of a rational function reproduces it exactly.
	// 1/(1+x) has Taylor 1,-1,1,-1,...  [0/1] Pade = 1/(1+x)
	tay := []float64{1, -1, 1, -1, 1}
	pr, _ := PadeApprox(tay, 0, 1)
	if !close(pr.Eval(0.3), 1/1.3, 1e-12) {
		t.Errorf("pade rational eval=%v", pr.Eval(0.3))
	}
}

func TestPolyFit(t *testing.T) {
	xs := []float64{0, 1, 2, 3, 4, 5}
	ys := make([]float64, len(xs))
	for i, x := range xs {
		ys[i] = 2 + 3*x - x*x
	}
	res, err := PolyFit(xs, ys, 2)
	if err != nil {
		t.Fatal(err)
	}
	if !sliceClose(res.Coeffs, []float64{2, 3, -1}, 1e-6) {
		t.Errorf("PolyFit coeffs=%v", res.Coeffs)
	}
	if !close(res.R2, 1, 1e-9) {
		t.Errorf("R2=%v", res.R2)
	}
	// noisy linear data: slope/intercept recovery
	lx := []float64{0, 1, 2, 3}
	ly := []float64{1.1, 1.9, 3.2, 3.8}
	lf, _ := PolyFit(lx, ly, 1)
	if !close(lf.Coeffs[1], 0.94, 0.1) {
		t.Errorf("linear slope=%v", lf.Coeffs[1])
	}
	// discrete minimax on exact quadratic returns it with zero error
	dm, de, _ := DiscreteMinimaxPoly(xs, ys, 2, 50)
	if !sliceClose(dm, []float64{2, 3, -1}, 1e-6) || de > 1e-6 {
		t.Errorf("DiscreteMinimax coeffs=%v err=%v", dm, de)
	}
	// Chebyshev least squares
	cs, _ := ChebLeastSquares(xs, ys, 2, 0, 5)
	if !close(cs.Eval(2.5), 2+3*2.5-2.5*2.5, 1e-6) {
		t.Errorf("ChebLeastSquares eval=%v", cs.Eval(2.5))
	}
}

func TestThiele(t *testing.T) {
	// interpolate a rational function; Thiele reproduces it well
	f := func(x float64) float64 { return (1 + x) / (1 + x*x) }
	xs := []float64{-2, -1, 0, 1, 2, 3}
	ys := make([]float64, len(xs))
	for i, x := range xs {
		ys[i] = f(x)
	}
	th, err := NewThiele(xs, ys)
	if err != nil {
		t.Fatal(err)
	}
	// reproduces data exactly
	for i, x := range xs {
		if !close(th.Eval(x), ys[i], 1e-9) {
			t.Errorf("Thiele node %v: %v want %v", x, th.Eval(x), ys[i])
		}
	}
	// interpolates between nodes accurately for a rational target
	if !close(th.Eval(0.5), f(0.5), 1e-9) {
		t.Errorf("Thiele mid=%v want %v", th.Eval(0.5), f(0.5))
	}
	v, _ := ThieleInterp(xs, ys, 1.5)
	if !close(v, f(1.5), 1e-9) {
		t.Errorf("ThieleInterp=%v want %v", v, f(1.5))
	}
}

func TestBernstein(t *testing.T) {
	// basis partition of unity
	for _, tval := range []float64{0, 0.3, 0.7, 1} {
		var sum float64
		for k := 0; k <= 5; k++ {
			sum += BernsteinBasis(k, 5, tval)
		}
		if !close(sum, 1, 1e-12) {
			t.Errorf("partition of unity at %v = %v", tval, sum)
		}
	}
	// Bernstein reproduces linear functions exactly
	lin := func(x float64) float64 { return 2 + 3*x }
	bp := BernsteinApprox(lin, 4, 0, 2)
	if !close(bp.Eval(1.3), lin(1.3), 1e-9) {
		t.Errorf("Bernstein linear=%v want %v", bp.Eval(1.3), lin(1.3))
	}
	// Bernstein converges to a continuous function
	be := BernsteinApprox(math.Sin, 200, 0, math.Pi)
	if err := MaxError(math.Sin, be.Eval, 0, math.Pi); err > 0.05 {
		t.Errorf("Bernstein sin error=%v", err)
	}
	// control points from monomial reproduce the polynomial
	ctrl := BernsteinFromMonomial([]float64{1, 2, 3}, 2)
	poly := NewBernsteinPoly(ctrl, 0, 1)
	if !close(poly.Eval(0.7), 1+2*0.7+3*0.49, 1e-9) {
		t.Errorf("BernsteinFromMonomial eval=%v", poly.Eval(0.7))
	}
	// ToMonomial inverts BernsteinFromMonomial
	if !sliceClose(poly.ToMonomial(), []float64{1, 2, 3}, 1e-9) {
		t.Errorf("Bernstein ToMonomial=%v", poly.ToMonomial())
	}
	// degree elevation preserves the polynomial
	el := poly.Elevate()
	if !close(el.Eval(0.4), poly.Eval(0.4), 1e-12) {
		t.Errorf("Elevate changed values")
	}
	// derivative of Bernstein form
	d := poly.Derivative()
	if !close(d.Eval(0.5), 2+6*0.5, 1e-9) { // d/dx(1+2x+3x^2)=2+6x
		t.Errorf("Bernstein derivative=%v", d.Eval(0.5))
	}
	// integral of 1+2x+3x^2 over [0,1] = 1+1+1 = 3
	if !close(poly.Integral(), 3, 1e-9) {
		t.Errorf("Bernstein integral=%v", poly.Integral())
	}
}

func TestQuadrature(t *testing.T) {
	want := math.Exp(1) - math.Exp(-1)
	if got := ClenshawCurtisQuadrature(math.Exp, 24, -1, 1); !close(got, want, 1e-10) {
		t.Errorf("ClenshawCurtis=%v want %v", got, want)
	}
	if got := FejerQuadrature(math.Exp, 24, -1, 1); !close(got, want, 1e-8) {
		t.Errorf("Fejer=%v want %v", got, want)
	}
	if got := ChebyshevQuadrature(math.Exp, 24, -1, 1); !close(got, want, 1e-10) {
		t.Errorf("ChebyshevQuadrature=%v want %v", got, want)
	}
	// Gauss-Chebyshev rules converge more slowly; use a looser tolerance.
	if got := GaussChebyshevQuadrature(math.Exp, 60, -1, 1); !close(got, want, 1e-3) {
		t.Errorf("GaussChebyshev=%v want %v", got, want)
	}
	if got := GaussChebyshev2Quadrature(math.Exp, 60, -1, 1); !close(got, want, 1e-3) {
		t.Errorf("GaussChebyshev2=%v want %v", got, want)
	}
	// Clenshaw-Curtis weights sum to the interval length.
	w := ClenshawCurtisWeights(10, 0, 4)
	var sum float64
	for _, wi := range w {
		sum += wi
	}
	if !close(sum, 4, 1e-12) {
		t.Errorf("CC weights sum=%v", sum)
	}
}

func TestBarycentricAndInterp(t *testing.T) {
	f := math.Exp
	bc := ChebyshevInterpolant(f, 30, -1, 1)
	if err := MaxError(f, bc.Eval, -1, 1); err > 1e-12 {
		t.Errorf("Chebyshev barycentric error=%v", err)
	}
	// barycentric reproduces node values
	xs := []float64{0, 1, 2, 4}
	ys := []float64{1, 3, 2, 5}
	b, _ := NewBarycentric(xs, ys)
	for i, x := range xs {
		if !close(b.Eval(x), ys[i], 1e-9) {
			t.Errorf("barycentric node reproduction failed")
		}
	}
	// Lagrange, Newton and barycentric agree off-node
	x0 := 1.7
	lg := LagrangeEval(xs, ys, x0)
	dd := NewtonDividedDifferences(xs, ys)
	nv := NewtonEval(dd, xs, x0)
	if !close(lg, b.Eval(x0), 1e-9) || !close(lg, nv, 1e-9) {
		t.Errorf("interp forms disagree: lag=%v bary=%v newton=%v", lg, b.Eval(x0), nv)
	}
	// Neville agrees too
	if !close(NevilleEval(xs, ys, x0), lg, 1e-9) {
		t.Errorf("Neville disagrees")
	}
}

func TestErrorsAndLebesgue(t *testing.T) {
	// Lebesgue constant of Chebyshev nodes is close to the asymptotic estimate.
	for _, n := range []int{10, 20, 40} {
		lc := LebesgueConstant(ChebPoints(n, -1, 1), -1, 1)
		asymp := LebesgueConstantChebyshevAsymptotic(n)
		if math.Abs(lc-asymp) > 0.2 {
			t.Errorf("n=%d Lebesgue %v vs asymptotic %v", n, lc, asymp)
		}
	}
	// Chebyshev nodes give a smaller Lebesgue constant than equispaced nodes.
	if LebesgueConstant(ChebPoints(16, -1, 1), -1, 1) >= LebesgueConstant(EquispacedNodes(16, -1, 1), -1, 1) {
		t.Errorf("Chebyshev Lebesgue not smaller than equispaced")
	}
	// node polynomial and error bound
	nodes := ChebPoints(4, -1, 1)
	if MaxNodePolynomial(nodes, -1, 1) <= 0 {
		t.Errorf("MaxNodePolynomial nonpositive")
	}
	// error norms
	f := math.Sin
	g := func(x float64) float64 { return x - x*x*x/6 } // 3rd order Taylor
	if L2Error(f, g, -0.5, 0.5) <= 0 {
		t.Errorf("L2Error nonpositive")
	}
	a := []float64{1, 2, 3}
	b := []float64{1, 2, 5}
	if !close(MaxAbsError(a, b), 2, 1e-12) {
		t.Errorf("MaxAbsError=%v", MaxAbsError(a, b))
	}
	if !close(Factorial(5), 120, 0) {
		t.Errorf("Factorial=%v", Factorial(5))
	}
	if !close(Binomial(6, 2), 15, 0) {
		t.Errorf("Binomial=%v", Binomial(6, 2))
	}
}

func TestTaylorSeries(t *testing.T) {
	if !sliceClose(TaylorExp(3), []float64{1, 1, 0.5, 1.0 / 6}, 1e-12) {
		t.Errorf("TaylorExp wrong")
	}
	if !sliceClose(TaylorSin(5), []float64{0, 1, 0, -1.0 / 6, 0, 1.0 / 120}, 1e-12) {
		t.Errorf("TaylorSin wrong")
	}
	if !sliceClose(TaylorCos(4), []float64{1, 0, -0.5, 0, 1.0 / 24}, 1e-12) {
		t.Errorf("TaylorCos wrong")
	}
	if !sliceClose(TaylorLog1p(3), []float64{0, 1, -0.5, 1.0 / 3}, 1e-12) {
		t.Errorf("TaylorLog1p wrong")
	}
	if !sliceClose(TaylorAtan(5), []float64{0, 1, 0, -1.0 / 3, 0, 1.0 / 5}, 1e-12) {
		t.Errorf("TaylorAtan wrong")
	}
	// verify a Taylor series against the true function at a small argument
	if !close(Polyval(TaylorExp(12), 0.5), math.Exp(0.5), 1e-9) {
		t.Errorf("TaylorExp eval mismatch")
	}
}

func TestContinuedFractionAndCompose(t *testing.T) {
	// 1 + 1/(1 + 1/(1 + 1/1)) evaluates the CF [1;1,1,1]
	got := ContinuedFractionEval([]float64{1, 1, 1, 1}, []float64{0, 1, 1, 1})
	if !close(got, 1+1.0/(1+1.0/(1+1.0/1)), 1e-12) {
		t.Errorf("ContinuedFractionEval=%v", got)
	}
	// PolyCompose: p(q(x)) with p=x^2, q=x+1 -> (x+1)^2
	c := PolyCompose([]float64{0, 0, 1}, []float64{1, 1})
	if !sliceClose(c, []float64{1, 2, 1}, 1e-12) {
		t.Errorf("PolyCompose=%v", c)
	}
	// PolyShift: p(x+1) for p=x^2
	s := PolyShift([]float64{0, 0, 1}, 1)
	if !sliceClose(s, []float64{1, 2, 1}, 1e-12) {
		t.Errorf("PolyShift=%v", s)
	}
}

func TestEconomize(t *testing.T) {
	// x^3 on [-1,1] = (3x + T3)/4 in Chebyshev terms; dropping the T3 term
	// (magnitude 0.25) gives (3/4)x with error 0.25.
	mono := []float64{0, 0, 0, 1}
	reduced, dropped := EconomizeChebyshev(mono, -1, 1, 0.3)
	if !close(dropped, 0.25, 1e-9) {
		t.Errorf("economize dropped=%v want 0.25", dropped)
	}
	if !sliceClose(reduced, []float64{0, 0.75}, 1e-9) {
		t.Errorf("economize reduced=%v", reduced)
	}
}

// Example demonstrates fitting exp with a Chebyshev series and evaluating it.
func ExampleChebFit() {
	s := ChebFit(math.Exp, 15, -1, 1)
	fmt.Printf("%.10f\n", s.Eval(0.5))
	// Output: 1.6487212707
}

// Example shows a [1/1] Pade approximant of the exponential function.
func ExamplePadeApprox() {
	p, _ := PadeApprox(TaylorExp(2), 1, 1)
	fmt.Println(p.Num, p.Den)
	// Output: [1 0.5] [1 -0.5]
}

// Example computes the best degree-1 minimax polynomial of x^2 on [-1,1].
func ExampleRemezPoly() {
	r, _ := RemezPoly(func(x float64) float64 { return x * x }, 1, -1, 1, 100, 1e-6)
	fmt.Printf("%.4f %.4f\n", r.Coeffs[0], r.Error)
	// Output: 0.5000 0.5000
}

// Example integrates exp over [-1,1] with Clenshaw-Curtis quadrature.
func ExampleClenshawCurtisQuadrature() {
	got := ClenshawCurtisQuadrature(math.Exp, 20, -1, 1)
	fmt.Printf("%.10f\n", got)
	// Output: 2.3504023873
}
