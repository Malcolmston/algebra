package diffalgebra

import (
	"fmt"
	"math"
	"math/big"
	"testing"
)

const tol = 1e-7

func approx(a, b float64) bool { return math.Abs(a-b) <= 1e-6*(1+math.Abs(b)) }

func ri(n int64) *big.Rat { return big.NewRat(n, 1) }

// ---------------------------------------------------------------------------
// Polynomials
// ---------------------------------------------------------------------------

func TestPolyArithmetic(t *testing.T) {
	a := PolyFromInts(1, 2, 3) // 3x^2+2x+1
	b := PolyFromInts(0, 1)    // x
	tests := []struct {
		name string
		got  Poly
		want Poly
	}{
		{"add", a.Add(b), PolyFromInts(1, 3, 3)},
		{"sub", a.Sub(b), PolyFromInts(1, 1, 3)},
		{"mul", a.Mul(b), PolyFromInts(0, 1, 2, 3)},
		{"neg", a.Neg(), PolyFromInts(-1, -2, -3)},
		{"scalar", a.ScalarMul(ri(2)), PolyFromInts(2, 4, 6)},
		{"pow", b.Pow(3), PolyFromInts(0, 0, 0, 1)},
		{"deriv", a.Derivative(), PolyFromInts(2, 6)},
		{"integral", PolyFromInts(2, 6).Integral(), PolyFromInts(0, 2, 3)},
	}
	for _, tc := range tests {
		if !tc.got.Equal(tc.want) {
			t.Errorf("%s: got %s want %s", tc.name, tc.got, tc.want)
		}
	}
}

func TestPolyDivMod(t *testing.T) {
	tests := []struct {
		num, den, q, r []int64
	}{
		{[]int64{-1, 0, 1}, []int64{-1, 1}, []int64{1, 1}, nil},       // (x^2-1)/(x-1)=x+1
		{[]int64{0, 0, 1}, []int64{1, 1}, []int64{-1, 1}, []int64{1}}, // x^2/(x+1)=x-1 r 1
	}
	for i, tc := range tests {
		num := PolyFromInts(tc.num...)
		den := PolyFromInts(tc.den...)
		q, r, err := num.DivMod(den)
		if err != nil {
			t.Fatalf("case %d: %v", i, err)
		}
		if !q.Equal(PolyFromInts(tc.q...)) {
			t.Errorf("case %d: quotient got %s", i, q)
		}
		if !r.Equal(PolyFromInts(tc.r...)) {
			t.Errorf("case %d: remainder got %s", i, r)
		}
		// verify num == q*den + r
		if !q.Mul(den).Add(r).Equal(num) {
			t.Errorf("case %d: reconstruction failed", i)
		}
	}
}

func TestPolyGCD(t *testing.T) {
	// (x^2-1) and (x^2-2x+1) share (x-1)
	a := PolyFromInts(-1, 0, 1)
	b := PolyFromInts(1, -2, 1)
	g := a.GCD(b)
	if !g.Equal(PolyFromInts(-1, 1)) {
		t.Errorf("gcd got %s want x-1", g)
	}
	// extended gcd identity s*a+t*b = g
	gg, s, tt := a.ExtendedGCD(b)
	if !s.Mul(a).Add(tt.Mul(b)).Equal(gg) {
		t.Errorf("extended gcd identity failed")
	}
}

func TestPolyResultantDiscriminant(t *testing.T) {
	// disc(x^2+bx+c) = b^2-4c ; for x^2-5x+6 -> 25-24=1
	p := PolyFromInts(6, -5, 1)
	if p.Discriminant().Cmp(ri(1)) != 0 {
		t.Errorf("disc got %s want 1", p.Discriminant().RatString())
	}
	// resultant of coprime polynomials is nonzero, of sharing ones is zero
	if p.Resultant(PolyFromInts(-6, 1)).Sign() == 0 {
		// x-6 not a root of p, so nonzero
		t.Errorf("expected nonzero resultant")
	}
	if PolyFromInts(-2, 1).Resultant(PolyFromInts(-2, 1)).Sign() != 0 {
		t.Errorf("expected zero resultant for equal factors")
	}
}

func TestPolyRationalRoots(t *testing.T) {
	// 2x^2 - x - 1 = (2x+1)(x-1): roots 1 and -1/2
	p := PolyFromInts(-1, -1, 2)
	roots := p.RationalRoots()
	found := map[string]bool{}
	for _, r := range roots {
		found[r.RatString()] = true
	}
	if !found["1"] || !found["-1/2"] {
		t.Errorf("roots got %v", roots)
	}
}

func TestSquareFree(t *testing.T) {
	// (x-1)^2 (x+2) : squarefree part (x-1)(x+2)=x^2+x-2
	p := PolyFromInts(1, -2, 1).Mul(PolyFromInts(2, 1)) // (x-1)^2*(x+2)
	sf := p.SquareFreePart()
	want := PolyFromInts(-2, 1, 1) // x^2+x-2
	if !sf.Equal(want) {
		t.Errorf("squarefree part got %s want %s", sf, want)
	}
	facs := p.SquareFreeFactorization()
	// expect (x-1)^2 and (x+2)^1
	mults := map[int]int{}
	for _, f := range facs {
		mults[f.Mult] += f.Factor.Degree()
	}
	if mults[2] != 1 || mults[1] != 1 {
		t.Errorf("squarefree factorization structure wrong: %+v", facs)
	}
}

// ---------------------------------------------------------------------------
// Rational functions
// ---------------------------------------------------------------------------

func TestRatFuncArithmetic(t *testing.T) {
	f, _ := NewRatFunc(OnePoly(), XPoly())            // 1/x
	g, _ := NewRatFunc(OnePoly(), PolyFromInts(1, 1)) // 1/(x+1)
	sum := f.Add(g)                                   // (2x+1)/(x(x+1))
	if v := sum.EvalFloat(1); !approx(v, 1.5) {
		t.Errorf("sum eval got %v want 1.5", v)
	}
	prod := f.Mul(g)
	if v := prod.EvalFloat(1); !approx(v, 0.5) {
		t.Errorf("prod eval got %v want 0.5", v)
	}
	// derivative of 1/x is -1/x^2
	d := f.Derivative()
	if v := d.EvalFloat(2); !approx(v, -0.25) {
		t.Errorf("deriv got %v want -0.25", v)
	}
}

func TestLogDerivative(t *testing.T) {
	// f = x^2 ; f'/f = 2/x
	f := RatFuncFromPoly(XPoly().Pow(2))
	ld, err := f.LogDerivative()
	if err != nil {
		t.Fatal(err)
	}
	want, _ := NewRatFunc(ConstPolyInt(2), XPoly())
	if !ld.Equal(want) {
		t.Errorf("log derivative got %s want 2/x", ld)
	}
}

func TestPartialFractions(t *testing.T) {
	// 1/(x^2 (x+1)) = -1/x + 1/x^2 + 1/(x+1): three terms grouped by the
	// square-free factors x^2 and (x+1).
	den := XPoly().Pow(2).Mul(PolyFromInts(1, 1))
	f, _ := NewRatFunc(OnePoly(), den)
	poly, terms := f.PartialFractions()
	if !poly.IsZero() {
		t.Errorf("expected no polynomial part, got %s", poly)
	}
	if len(terms) != 3 {
		t.Fatalf("expected 3 terms got %d", len(terms))
	}
	// reconstruct and compare
	recon := RatFuncFromPoly(poly)
	for _, tm := range terms {
		den := tm.Factor.Pow(tm.Power)
		part, _ := NewRatFunc(tm.Numerator, den)
		recon = recon.Add(part)
	}
	if !recon.Equal(f) {
		t.Errorf("partial fraction reconstruction failed: %s", recon)
	}
}

// ---------------------------------------------------------------------------
// Derivation and operators
// ---------------------------------------------------------------------------

func TestDerivationLeibniz(t *testing.T) {
	d := StandardDerivation()
	f, _ := NewRatFunc(OnePoly(), XPoly())
	g := RatFuncFromPoly(PolyFromInts(0, 0, 1))
	direct, leib := d.Leibniz(f, g)
	if !direct.Equal(leib) {
		t.Errorf("leibniz mismatch: %s vs %s", direct, leib)
	}
	if !d.IsConstant(ConstRatFunc(ri(3))) {
		t.Errorf("constant not recognised")
	}
}

func TestOperatorAlgebra(t *testing.T) {
	// D*x = x*D + 1
	Dx := DOperator().Mul(NewOperator(XRatFunc()))
	want := NewOperator(OneRatFunc(), XRatFunc()) // 1 + x*D
	if !Dx.Equal(want) {
		t.Errorf("D*x got %s want %s", Dx, want)
	}
	// apply (D*x) to x^2 => D(x^3)=3x^2
	got := Dx.ApplyPoly(XPoly().Pow(2))
	if !got.Equal(RatFuncFromPoly(PolyFromInts(0, 0, 3))) {
		t.Errorf("(D*x)(x^2) got %s want 3x^2", got)
	}
	// order and coeffs
	if Dx.Order() != 1 {
		t.Errorf("order got %d", Dx.Order())
	}
	// D^2 applied to x^3 = 6x
	got2 := DOperator().Pow(2).ApplyPoly(XPoly().Pow(3))
	if !got2.Equal(RatFuncFromPoly(PolyFromInts(0, 6))) {
		t.Errorf("D^2(x^3) got %s want 6x", got2)
	}
}

func TestIndicialPolynomial(t *testing.T) {
	// Euler operator x^2 D^2 + x D - 1 has indicial poly s^2 - 1 (roots +-1)
	op := OperatorFromPolys(ConstPolyInt(-1)).
		Add(NewOperator(ZeroRatFunc(), XRatFunc())).
		Add(NewOperator(ZeroRatFunc(), ZeroRatFunc(), RatFuncFromPoly(XPoly().Pow(2))))
	ind, ok := op.IndicialPolynomial()
	if !ok {
		t.Fatal("indicial polynomial not computed")
	}
	roots := ind.RationalRoots()
	found := map[string]bool{}
	for _, r := range roots {
		found[r.RatString()] = true
	}
	if !found["1"] || !found["-1"] {
		t.Errorf("indicial roots got %v (poly %s)", roots, ind)
	}
}

// ---------------------------------------------------------------------------
// Wronskians
// ---------------------------------------------------------------------------

func TestWronskian(t *testing.T) {
	tests := []struct {
		name string
		fns  []Poly
		want Poly
	}{
		{"1,x,x^2", []Poly{OnePoly(), XPoly(), XPoly().Pow(2)}, ConstPolyInt(2)},
		{"1,x", []Poly{OnePoly(), XPoly()}, OnePoly()},
		{"x,x^2", []Poly{XPoly(), XPoly().Pow(2)}, XPoly().Pow(2)},
	}
	for _, tc := range tests {
		w, err := WronskianPoly(tc.fns)
		if err != nil {
			t.Fatalf("%s: %v", tc.name, err)
		}
		if !w.Equal(tc.want) {
			t.Errorf("%s: wronskian got %s want %s", tc.name, w, tc.want)
		}
	}
	// dependent set {x, 2x} has zero Wronskian
	if LinearlyIndependentPoly([]Poly{XPoly(), XPoly().ScalarMul(ri(2))}) {
		t.Errorf("dependent set reported independent")
	}
	if !LinearlyIndependentPoly([]Poly{OnePoly(), XPoly()}) {
		t.Errorf("independent set reported dependent")
	}
}

// ---------------------------------------------------------------------------
// Rational integration
// ---------------------------------------------------------------------------

func TestIntegrateRational(t *testing.T) {
	tests := []struct {
		name    string
		num     []int64
		den     []int64
		check   float64 // sample point
		want    float64 // expected antiderivative value (abs-log convention)
		numLogs int
	}{
		{"1/x^2", []int64{1}, []int64{0, 0, 1}, 2, -0.5, 0},
		{"2x/(x^2+1)", []int64{0, 2}, []int64{1, 0, 1}, 1, math.Log(2), 1},
		{"1/(x^2-1)", []int64{1}, []int64{-1, 0, 1}, 3, 0.5*math.Log(2) - 0.5*math.Log(4), 2},
		{"1/x", []int64{1}, []int64{0, 1}, math.E, 1, 1},
		{"x^3", []int64{0, 0, 0, 1}, []int64{1}, 2, 4, 0}, // integral = x^4/4 => 4
	}
	for _, tc := range tests {
		f, err := NewRatFunc(PolyFromInts(tc.num...), PolyFromInts(tc.den...))
		if err != nil {
			t.Fatalf("%s: %v", tc.name, err)
		}
		res, err := IntegrateRational(f)
		if err != nil {
			t.Fatalf("%s: %v", tc.name, err)
		}
		if len(res.Logs) != tc.numLogs {
			t.Errorf("%s: got %d logs want %d (%s)", tc.name, len(res.Logs), tc.numLogs, res)
		}
		if v := res.EvalFloat(tc.check); !approx(v, tc.want) {
			t.Errorf("%s: eval got %v want %v", tc.name, v, tc.want)
		}
		// verify derivative of the computed integral matches f (rational part only,
		// away from log arguments): d/dx(rational + logs) = f.
		checkDerivativeMatches(t, tc.name, res, f, tc.check)
	}
}

// checkDerivativeMatches differentiates the closed-form integral numerically and
// compares to f at the sample point.
func checkDerivativeMatches(t *testing.T, name string, res RationalIntegral, f RatFunc, x float64) {
	h := 1e-6
	d := (res.EvalFloat(x+h) - res.EvalFloat(x-h)) / (2 * h)
	want := f.EvalFloat(x)
	if math.Abs(d-want) > 1e-3*(1+math.Abs(want)) {
		t.Errorf("%s: numeric derivative %v != f %v", name, d, want)
	}
}

func TestHermiteReduce(t *testing.T) {
	// 1/(x-1)^2 : rational part -1/(x-1), no log remainder
	f, _ := NewRatFunc(OnePoly(), PolyFromInts(1, -2, 1))
	rational, remaining := HermiteReduce(f)
	if !remaining.IsZero() {
		t.Errorf("expected no log remainder, got %s", remaining)
	}
	want, _ := NewRatFunc(ConstPolyInt(-1), PolyFromInts(-1, 1))
	if !rational.Equal(want) {
		t.Errorf("hermite rational got %s want -1/(x-1)", rational)
	}
}

// ---------------------------------------------------------------------------
// Risch exponential heuristic
// ---------------------------------------------------------------------------

func TestRischExp(t *testing.T) {
	tests := []struct {
		name   string
		f, g   Poly
		wantOK bool
		wantR  Poly
	}{
		{"2x e^{x^2}", PolyFromInts(0, 2), PolyFromInts(0, 0, 1), true, OnePoly()},
		{"e^x", OnePoly(), PolyFromInts(0, 1), true, OnePoly()},
		{"x e^x", PolyFromInts(0, 1), PolyFromInts(0, 1), true, PolyFromInts(-1, 1)}, // (x-1)e^x
		{"e^{x^2}", OnePoly(), PolyFromInts(0, 0, 1), false, ZeroPoly()},             // not elementary
	}
	for _, tc := range tests {
		R, ok := RischExpIntegrate(tc.f, tc.g)
		if ok != tc.wantOK {
			t.Errorf("%s: ok got %v want %v", tc.name, ok, tc.wantOK)
			continue
		}
		if ok {
			// verify R' + g' R = f
			lhs := R.Derivative().Add(tc.g.Derivative().Mul(R))
			if !lhs.Equal(tc.f) {
				t.Errorf("%s: R'+g'R got %s want %s", tc.name, lhs, tc.f)
			}
			if !R.Equal(tc.wantR) {
				t.Errorf("%s: R got %s want %s", tc.name, R, tc.wantR)
			}
		}
	}
}

// ---------------------------------------------------------------------------
// Constant-coefficient ODEs
// ---------------------------------------------------------------------------

func TestConstantODE(t *testing.T) {
	// y'' + y = 0 -> cos, sin ; y(0)=1, y'(0)=0 -> cos(x)
	c, sol, err := SolveODEIVP([]float64{1, 0, 1}, 0, []float64{1, 0}, 7, 1e-6)
	if err != nil {
		t.Fatal(err)
	}
	for _, x := range []float64{0, 0.5, 1.0, 2.0} {
		if v := sol.Evaluate(c, x); !approx(v, math.Cos(x)) {
			t.Errorf("cos ivp at %v got %v want %v", x, v, math.Cos(x))
		}
	}
	// y'' - 3y' + 2y = 0 -> e^x, e^2x ; y(0)=0, y'(0)=1 -> e^{2x}-e^x
	c2, sol2, err := SolveODEIVP([]float64{2, -3, 1}, 0, []float64{0, 1}, 7, 1e-6)
	if err != nil {
		t.Fatal(err)
	}
	for _, x := range []float64{0, 0.5, 1.0} {
		want := math.Exp(2*x) - math.Exp(x)
		if v := sol2.Evaluate(c2, x); !approx(v, want) {
			t.Errorf("exp ivp at %v got %v want %v", x, v, want)
		}
	}
	// repeated root: y'' - 2y' + y = 0 -> e^x, x e^x
	sol3, err := SolveLinearConstantODE([]float64{1, -2, 1}, 7, 1e-4)
	if err != nil {
		t.Fatal(err)
	}
	if sol3.Dimension() != 2 {
		t.Errorf("repeated root dimension got %d want 2", sol3.Dimension())
	}
}

func TestConstantODESolutionSatisfies(t *testing.T) {
	// Verify each basis solution actually satisfies y''+y=0 numerically.
	sol, _ := SolveLinearConstantODE([]float64{1, 0, 1}, 3, 1e-6)
	for _, term := range sol.Basis() {
		x := 0.7
		ypp := term.nthDerivValue(2, x)
		y := term.Eval(x)
		if math.Abs(ypp+y) > 1e-6 {
			t.Errorf("basis term %s does not satisfy y''+y=0: %v", term, ypp+y)
		}
	}
}

// ---------------------------------------------------------------------------
// Recurrences
// ---------------------------------------------------------------------------

func TestRecurrenceFibonacci(t *testing.T) {
	// a_{n+2} = a_{n+1} + a_n : coeffs [-1,-1,1]
	coeffs := []float64{-1, -1, 1}
	c, sol, err := SolveRecurrenceIVP(coeffs, []float64{0, 1}, 11, 1e-6)
	if err != nil {
		t.Fatal(err)
	}
	fib := []float64{0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55}
	for n, want := range fib {
		if v := sol.Evaluate(c, n); math.Abs(v-want) > 1e-4 {
			t.Errorf("fib(%d) got %v want %v", n, v, want)
		}
	}
	// direct iteration cross-check
	vals, err := RecurrenceValues(coeffs, []float64{0, 1}, 11)
	if err != nil {
		t.Fatal(err)
	}
	for n, want := range fib {
		if math.Abs(vals[n]-want) > 1e-9 {
			t.Errorf("iterated fib(%d) got %v want %v", n, vals[n], want)
		}
	}
}

func TestRecurrenceGeometric(t *testing.T) {
	// a_{n+1} = 2 a_n : coeffs [-2, 1], a_0 = 3 -> 3*2^n
	c, sol, err := SolveRecurrenceIVP([]float64{-2, 1}, []float64{3}, 5, 1e-6)
	if err != nil {
		t.Fatal(err)
	}
	for n := 0; n < 6; n++ {
		want := 3 * math.Pow(2, float64(n))
		if v := sol.Evaluate(c, n); math.Abs(v-want) > 1e-6 {
			t.Errorf("geom(%d) got %v want %v", n, v, want)
		}
	}
}

// ---------------------------------------------------------------------------
// Variation of parameters
// ---------------------------------------------------------------------------

func TestVariationOfParameters(t *testing.T) {
	// Fundamental system {1, x} for y'' = g. Integrands: u1' = -x g, u2' = g.
	ys := []RatFunc{OneRatFunc(), XRatFunc()}
	g := RatFuncFromPoly(PolyFromInts(0, 0, 1)) // g = x^2
	integ, err := VariationOfParametersIntegrands(ys, g)
	if err != nil {
		t.Fatal(err)
	}
	// u1' = -x * x^2 = -x^3 ; u2' = x^2
	if !integ[0].Equal(RatFuncFromPoly(PolyFromInts(0, 0, 0, -1))) {
		t.Errorf("u1' got %s want -x^3", integ[0])
	}
	if !integ[1].Equal(RatFuncFromPoly(PolyFromInts(0, 0, 1))) {
		t.Errorf("u2' got %s want x^2", integ[1])
	}
	// Particular solution y_p = u1*1 + u2*x with u1 = -x^4/4, u2 = x^3/3
	// => y_p = -x^4/4 + x^4/3 = x^4/12
	yp, err := VariationOfParameters(ys, g)
	if err != nil {
		t.Fatal(err)
	}
	want, _ := NewRatFunc(PolyFromInts(0, 0, 0, 0, 1), ConstPolyInt(12))
	if !yp.Equal(want) {
		t.Errorf("particular solution got %s want x^4/12", yp)
	}
	// check y_p'' = x^2
	if !yp.Derivative().Derivative().Equal(g) {
		t.Errorf("y_p'' got %s want x^2", yp.Derivative().Derivative())
	}
}

// ---------------------------------------------------------------------------
// Kovacic
// ---------------------------------------------------------------------------

func TestKovacic(t *testing.T) {
	// y'' = (2/x^2) y has solution x^2 (and 1/x).
	r, _ := NewRatFunc(ConstPolyInt(2), PolyFromInts(0, 0, 1))
	res, err := Kovacic(r)
	if err != nil {
		t.Fatal(err)
	}
	if !res.Found {
		t.Fatal("expected Liouvillian solution")
	}
	// Verify the constructed solution satisfies y'' = r y numerically.
	for _, x := range []float64{1.3, 2.0, 3.5} {
		h := 1e-4
		ypp := (res.EvalFloat(x+h) - 2*res.EvalFloat(x) + res.EvalFloat(x-h)) / (h * h)
		rhs := r.EvalFloat(x) * res.EvalFloat(x)
		if math.Abs(ypp-rhs) > 1e-2*(1+math.Abs(rhs)) {
			t.Errorf("kovacic solution fails ODE at %v: y''=%v r*y=%v", x, ypp, rhs)
		}
	}
}

func TestKovacicReduceNormalForm(t *testing.T) {
	// y'' + 0*y' + 0*y => r = 0 (a2=1). Just check reduction runs.
	r, err := ReduceToNormalForm(OnePoly(), ZeroPoly(), ZeroPoly())
	if err != nil {
		t.Fatal(err)
	}
	if !r.IsZero() {
		t.Errorf("normal form r got %s want 0", r)
	}
	// Euler: x^2 y'' + x y' - y = 0 reduces and Kovacic finds a solution.
	r2, res, err := KovacicSolveSecondOrder(XPoly().Pow(2), XPoly(), ConstPolyInt(-1))
	if err != nil {
		t.Fatal(err)
	}
	_ = r2
	if !res.Found {
		t.Errorf("expected Kovacic solution for Euler equation")
	}
}

// ---------------------------------------------------------------------------
// Root finding
// ---------------------------------------------------------------------------

func TestComplexRoots(t *testing.T) {
	// x^2 + 1 -> +-i
	p := PolyFromInts(1, 0, 1)
	clusters := p.PolyComplexRootClusters(5, 1e-6)
	if len(clusters) != 2 {
		t.Fatalf("got %d clusters", len(clusters))
	}
	for _, cl := range clusters {
		if math.Abs(math.Abs(imag(cl.Value))-1) > 1e-6 || math.Abs(real(cl.Value)) > 1e-6 {
			t.Errorf("root got %v want +-i", cl.Value)
		}
	}
	// (x-2)^3 -> triple root 2
	q := PolyFromInts(-2, 1).Pow(3)
	cl := q.PolyComplexRootClusters(5, 1e-3)
	if len(cl) != 1 || cl[0].Mult != 3 {
		t.Errorf("triple root clustering got %+v", cl)
	}
}

// ---------------------------------------------------------------------------
// Example
// ---------------------------------------------------------------------------

func ExampleIntegrateRational() {
	// Integrate 1/(x^2 - 1) over Q.
	f, _ := NewRatFunc(OnePoly(), PolyFromInts(-1, 0, 1))
	res, _ := IntegrateRational(f)
	fmt.Println(res)
	// Output: 1/2*log(x - 1) + -1/2*log(x + 1)
}

func ExampleSolveLinearRecurrence() {
	// Fibonacci recurrence a_{n+2} = a_{n+1} + a_n with a_0=0, a_1=1.
	coeffs := []float64{-1, -1, 1}
	c, sol, _ := SolveRecurrenceIVP(coeffs, []float64{0, 1}, 1, 1e-6)
	fmt.Printf("a_10 = %.0f\n", sol.Evaluate(c, 10))
	// Output: a_10 = 55
}

func ExampleRischExpIntegrate() {
	// Integrate 2x * e^(x^2): result is 1 * e^(x^2).
	R, ok := RischExpIntegrate(PolyFromInts(0, 2), PolyFromInts(0, 0, 1))
	fmt.Println(R, ok)
	// Output: 1 true
}
