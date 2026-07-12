package algebra

import (
	"fmt"
	"math"
	"math/big"
	"math/cmplx"
	"testing"
)

// numDeriv cross-checks a symbolic derivative against a central finite
// difference at the given points.
func numDeriv(t *testing.T, expr string, pts []float64) {
	t.Helper()
	x := Sym("x")
	e := mustParse(t, expr)
	d := Diff(e, x)
	const h = 1e-6
	for _, p := range pts {
		fp, e1 := Eval(e, map[string]float64{"x": p + h})
		fm, e2 := Eval(e, map[string]float64{"x": p - h})
		da, e3 := Eval(d, map[string]float64{"x": p})
		if e1 != nil || e2 != nil || e3 != nil {
			t.Fatalf("%s eval error at %g", expr, p)
		}
		approx := (fp - fm) / (2 * h)
		if math.Abs(approx-da) > 1e-4*(1+math.Abs(approx)) {
			t.Errorf("d/dx %s at %g: symbolic %g vs numeric %g", expr, p, da, approx)
		}
	}
}

func TestTrigExactValues(t *testing.T) {
	cases := []struct {
		got  Expr
		want Expr
	}{
		{Sin(piMul(big.NewRat(1, 6))), Rat(1, 2)},
		{Cos(piMul(big.NewRat(1, 3))), Rat(1, 2)},
		{Sin(Pi), Int(0)},
		{Cos(Pi), Int(-1)},
		{Sin(piMul(big.NewRat(1, 2))), Int(1)},
		{Cos(piMul(big.NewRat(1, 2))), Int(0)},
		{Tan(piMul(big.NewRat(1, 4))), Int(1)},
		{Sin(piMul(big.NewRat(-1, 6))), Rat(-1, 2)},
		{Cos(Mul(Int(2), Pi)), Int(1)},
	}
	for i, c := range cases {
		if !c.got.Equal(c.want) {
			t.Errorf("case %d: got %s, want %s", i, c.got, c.want)
		}
	}
	// sqrt-valued angles verified numerically.
	v, _ := Evalf(Cos(piMul(big.NewRat(1, 4))))
	if math.Abs(v-math.Sqrt2/2) > 1e-12 {
		t.Errorf("cos(pi/4) = %g", v)
	}
}

func TestTrigIdentities(t *testing.T) {
	if got := Simplify(MustParse("sin(x)^2 + cos(x)^2")); !got.Equal(Int(1)) {
		t.Errorf("sin^2+cos^2 = %s, want 1", got)
	}
	if got := Simplify(MustParse("3*sin(y)^2 + 3*cos(y)^2")); !got.Equal(Int(3)) {
		t.Errorf("3sin^2+3cos^2 = %s, want 3", got)
	}
	if got := Simplify(MustParse("cos(x)^2 - sin(x)^2")); !got.Equal(Cos(Mul(Int(2), Sym("x")))) {
		t.Errorf("cos^2-sin^2 = %s, want cos(2*x)", got)
	}
	if got := Simplify(MustParse("sin(x)/cos(x)")); !got.Equal(Tan(Sym("x"))) {
		t.Errorf("sin/cos = %s, want tan(x)", got)
	}
	// sin(x)*cos(x) -> sin(2x)/2, checked numerically.
	e := Simplify(MustParse("sin(x)*cos(x)"))
	for _, p := range []float64{0.3, 1.1, 2.4} {
		a, _ := Eval(e, map[string]float64{"x": p})
		if math.Abs(a-math.Sin(p)*math.Cos(p)) > 1e-12 {
			t.Errorf("sin*cos rewrite wrong at %g: %g", p, a)
		}
	}
}

func TestLogExpRules(t *testing.T) {
	x := Sym("x")
	if got := Simplify(Log(Exp(x))); !got.Equal(x) {
		t.Errorf("log(exp(x)) = %s", got)
	}
	if got := Simplify(Exp(Log(x))); !got.Equal(x) {
		t.Errorf("exp(log(x)) = %s", got)
	}
}

func TestInverseTrigHyperbolicDiff(t *testing.T) {
	inUnit := []float64{-0.4, 0.2, 0.6} // domain (-1,1)
	numDeriv(t, "asin(x)", inUnit)
	numDeriv(t, "acos(x)", inUnit)
	numDeriv(t, "atanh(x)", inUnit)
	allReal := []float64{0.3, 1.2, 2.5}
	numDeriv(t, "atan(x)", allReal)
	numDeriv(t, "acot(x)", allReal)
	numDeriv(t, "sinh(x)", allReal)
	numDeriv(t, "cosh(x)", allReal)
	numDeriv(t, "tanh(x)", allReal)
	numDeriv(t, "coth(x)", []float64{0.6, 1.2, 2.5})
	numDeriv(t, "sech(x)", allReal)
	numDeriv(t, "csch(x)", []float64{0.6, 1.2, 2.5})
	numDeriv(t, "asinh(x)", allReal)
	gt1 := []float64{1.4, 2.2, 3.1} // domain > 1
	numDeriv(t, "acosh(x)", gt1)
	numDeriv(t, "asec(x)", gt1)
	numDeriv(t, "acsc(x)", gt1)
	numDeriv(t, "sec(x)", []float64{0.3, 1.0, 2.4})
	numDeriv(t, "csc(x)", []float64{0.6, 1.2, 2.5})
	numDeriv(t, "cot(x)", []float64{0.6, 1.2, 2.5})
	numDeriv(t, "abs(x)", []float64{0.7, 1.5, 2.3})
	numDeriv(t, "erf(x)", allReal)
	numDeriv(t, "erfc(x)", allReal)
	numDeriv(t, "gamma(x)", []float64{1.3, 2.1, 3.4})
}

func TestSpecialFunctions(t *testing.T) {
	cases := []struct {
		got  Expr
		want Expr
	}{
		{Abs(Int(-3)), Int(3)},
		{Abs(Rat(-3, 2)), Rat(3, 2)},
		{Sign(Int(-2)), Int(-1)},
		{Sign(Int(0)), Int(0)},
		{Floor(Rat(7, 2)), Int(3)},
		{Ceil(Rat(7, 2)), Int(4)},
		{Floor(Rat(-7, 2)), Int(-4)},
		{Factorial(Int(5)), Int(120)},
		{Factorial(Int(0)), Int(1)},
		{Gamma(Int(5)), Int(24)},
		{Gamma(Rat(1, 2)), Sqrt(Pi)},
	}
	for i, c := range cases {
		if !c.got.Equal(c.want) {
			t.Errorf("case %d: got %s, want %s", i, c.got, c.want)
		}
	}
	// Beta(2,3) = 1/12.
	if v, _ := Evalf(Beta(Int(2), Int(3))); math.Abs(v-1.0/12) > 1e-12 {
		t.Errorf("Beta(2,3) = %g, want 1/12", v)
	}
	// erf/gamma numeric.
	if v, _ := Evalf(Erf(Flt(1))); math.Abs(v-math.Erf(1)) > 1e-12 {
		t.Errorf("erf(1) = %g", v)
	}
	if v, _ := Evalf(Gamma(Flt(4.5))); math.Abs(v-math.Gamma(4.5)) > 1e-9 {
		t.Errorf("gamma(4.5) = %g", v)
	}
}

func TestComplexArithmetic(t *testing.T) {
	if got := Pow(I, Int(2)); !got.Equal(Int(-1)) {
		t.Errorf("I^2 = %s, want -1", got)
	}
	if got := Pow(I, Int(3)); !got.Equal(neg(I)) {
		t.Errorf("I^3 = %s, want -I", got)
	}
	if got := Exp(Mul(I, Pi)); !got.Equal(Int(-1)) {
		t.Errorf("exp(I*pi) = %s, want -1", got)
	}
	if got := Exp(Mul(I, piMul(big.NewRat(1, 2)))); !got.Equal(I) {
		t.Errorf("exp(I*pi/2) = %s, want I", got)
	}
	// (2+3I)*(1-I) = 5 + I.
	prod := Simplify(Expand(Mul(Add(Int(2), Mul(Int(3), I)), Add(Int(1), neg(I)))))
	if !prod.Equal(Add(Int(5), I)) {
		t.Errorf("(2+3I)(1-I) = %s, want 5 + I", prod)
	}
	z := Add(Int(3), Mul(Int(4), I))
	if !Re(z).Equal(Int(3)) || !Im(z).Equal(Int(4)) {
		t.Errorf("Re/Im wrong: %s %s", Re(z), Im(z))
	}
	if !Conjugate(z).Equal(Add(Int(3), Mul(Int(-4), I))) {
		t.Errorf("Conjugate = %s", Conjugate(z))
	}
	if !Abs(z).Equal(Int(5)) {
		t.Errorf("Abs(3+4I) = %s, want 5", Abs(z))
	}
}

func TestEvalc(t *testing.T) {
	z, err := Evalc(MustParse("sqrt(-1)"))
	if err != nil || cmplx.Abs(z-1i) > 1e-12 {
		t.Errorf("sqrt(-1) = %v, %v", z, err)
	}
	z, _ = Evalc(Exp(Mul(I, piMul(big.NewRat(1, 3)))))
	want := cmplx.Exp(1i * math.Pi / 3)
	if cmplx.Abs(z-want) > 1e-12 {
		t.Errorf("exp(I*pi/3) = %v, want %v", z, want)
	}
}

// TestIntegrateExtended differentiates each antiderivative and checks it
// recovers the integrand numerically at sample points.
func TestIntegrateExtended(t *testing.T) {
	x := Sym("x")
	cases := []struct {
		expr string
		pts  []float64
	}{
		{"1/(1+x^2)", []float64{0.3, 1.2, 2.4}},
		{"1/sqrt(1-x^2)", []float64{-0.4, 0.2, 0.6}},
		{"tan(x)", []float64{0.3, 1.0}},
		{"cot(x)", []float64{0.6, 1.2}},
		{"sinh(x)", []float64{0.3, 1.4}},
		{"cosh(x)", []float64{0.3, 1.4}},
		{"tanh(x)", []float64{0.3, 1.4}},
		{"x*exp(x)", []float64{0.3, 1.4, 2.1}},
		{"x*sin(x)", []float64{0.3, 1.4, 2.1}},
		{"x^2*exp(x)", []float64{0.3, 1.4}},
		{"x^2*cos(x)", []float64{0.3, 1.4}},
		{"1/((x-3)*(x-5))", []float64{0.3, 1.4, 6.2}},
		{"1/(x^2+4)", []float64{0.3, 1.4, 2.1}},
		{"1/sqrt(x^2+1)", []float64{0.3, 1.4}},
		{"(x^2+1)/(x-2)", []float64{0.3, 3.4, 5.1}},
	}
	for _, c := range cases {
		f := mustParse(t, c.expr)
		F := Integrate(f, x)
		if _, un := F.(*integral); un {
			t.Errorf("Integrate(%s) unresolved", c.expr)
			continue
		}
		back := Diff(F, x)
		for _, p := range c.pts {
			a, e1 := Eval(f, map[string]float64{"x": p})
			b, e2 := Eval(back, map[string]float64{"x": p})
			if e1 != nil || e2 != nil {
				t.Fatalf("%s eval error", c.expr)
			}
			if math.Abs(a-b) > 1e-5*(1+math.Abs(a)) {
				t.Errorf("d/dx Integrate(%s) at %g: %g vs %g", c.expr, p, b, a)
			}
		}
	}
}

func TestSolveCubicQuartic(t *testing.T) {
	x := Sym("x")
	for _, s := range []string{
		"x^3 + 1", "x^3 - 6*x^2 + 11*x - 6", "x^4 - 1",
		"x^3 - 2", "x^4 - 5*x^2 + 4", "x^2 + 1", "2*x^2 + 2*x + 5",
	} {
		roots, err := Solve(mustParse(t, s), x)
		if err != nil {
			t.Fatalf("Solve(%s) error: %v", s, err)
		}
		if len(roots) == 0 {
			t.Errorf("Solve(%s) returned no roots", s)
		}
		for _, r := range roots {
			res, err := Evalc(Subs(mustParse(t, s), x, r))
			if err != nil {
				t.Fatalf("residual eval %s: %v", r, err)
			}
			if cmplx.Abs(res) > 1e-6 {
				t.Errorf("Solve(%s): root %s gives residual %v", s, r, res)
			}
		}
	}
	// Complex quadratic gives exact conjugate pair.
	roots, _ := Solve(mustParse(t, "x^2 + 1"), x)
	if len(roots) != 2 {
		t.Fatalf("x^2+1 should have 2 roots, got %v", roots)
	}
}

func TestSolveSystem(t *testing.T) {
	x, y, z := Sym("x"), Sym("y"), Sym("z")
	sol, err := SolveSystem([]Expr{MustParse("x + y - 3"), MustParse("x - y - 1")}, []Expr{x, y})
	if err != nil {
		t.Fatalf("SolveSystem error: %v", err)
	}
	if !sol[0].Equal(Int(2)) || !sol[1].Equal(Int(1)) {
		t.Errorf("system solution = %v, want [2 1]", sol)
	}
	// 3x3 system.
	sol, err = SolveSystem([]Expr{
		MustParse("x + y + z - 6"),
		MustParse("2*x - y + z - 3"),
		MustParse("x + 2*y - z - 2"),
	}, []Expr{x, y, z})
	if err != nil {
		t.Fatalf("3x3 error: %v", err)
	}
	if !sol[0].Equal(Int(1)) || !sol[1].Equal(Int(2)) || !sol[2].Equal(Int(3)) {
		t.Errorf("3x3 solution = %v, want [1 2 3]", sol)
	}
	// Singular system errors.
	if _, err := SolveSystem([]Expr{MustParse("x + y - 1"), MustParse("2*x + 2*y - 2")}, []Expr{x, y}); err == nil {
		t.Error("singular system should error")
	}
}

func TestSeries(t *testing.T) {
	x := Sym("x")
	// exp: 1 + x + x^2/2 + x^3/6 + x^4/24.
	got := Series(MustParse("exp(x)"), x, Int(0), 5)
	want := Add(Int(1), x, Mul(Rat(1, 2), Pow(x, Int(2))),
		Mul(Rat(1, 6), Pow(x, Int(3))), Mul(Rat(1, 24), Pow(x, Int(4))))
	if !got.Equal(want) {
		t.Errorf("Series(exp) = %s", got)
	}
	// 1/(1-x): 1 + x + x^2 + x^3.
	got = Series(MustParse("1/(1-x)"), x, Int(0), 4)
	if !got.Equal(Add(Int(1), x, Pow(x, Int(2)), Pow(x, Int(3)))) {
		t.Errorf("Series(1/(1-x)) = %s", got)
	}
	// sin matches numerically.
	s := Series(MustParse("sin(x)"), x, Int(0), 8)
	for _, p := range []float64{0.1, 0.5, 1.0} {
		v, _ := Eval(s, map[string]float64{"x": p})
		if math.Abs(v-math.Sin(p)) > 1e-4 {
			t.Errorf("sin series at %g: %g vs %g", p, v, math.Sin(p))
		}
	}
}

func TestLimit(t *testing.T) {
	x := Sym("x")
	cases := []struct {
		expr string
		to   Expr
		want Expr
	}{
		{"sin(x)/x", Int(0), Int(1)},
		{"(1-cos(x))/x^2", Int(0), Rat(1, 2)},
		{"(x^2 - 1)/(x - 1)", Int(1), Int(2)},
		{"(x^2+1)/(2*x^2+3)", Inf, Rat(1, 2)},
		{"(3*x^2+2)/(x^3+1)", Inf, Int(0)},
	}
	for _, c := range cases {
		got := Limit(mustParse(t, c.expr), x, c.to)
		if !got.Equal(c.want) {
			t.Errorf("Limit(%s -> %s) = %s, want %s", c.expr, c.to, got, c.want)
		}
	}
}

func TestSummationProduct(t *testing.T) {
	n := Sym("n")
	k := Sym("k")
	// sum_{k=1}^{n} k = n(n+1)/2.
	got := Summation(k, k, Int(1), n)
	for _, nv := range []int64{1, 5, 10} {
		v, _ := Eval(got, map[string]float64{"n": float64(nv)})
		want := float64(nv*(nv+1)) / 2
		if math.Abs(v-want) > 1e-9 {
			t.Errorf("sum k at n=%d: %g want %g", nv, v, want)
		}
	}
	// geometric sum_{k=0}^{n} 2^k = 2^(n+1) - 1.
	g := Summation(Pow(Int(2), k), k, Int(0), n)
	for _, nv := range []int64{0, 3, 6} {
		v, _ := Eval(g, map[string]float64{"n": float64(nv)})
		want := math.Pow(2, float64(nv+1)) - 1
		if math.Abs(v-want) > 1e-6 {
			t.Errorf("geo sum at n=%d: %g want %g", nv, v, want)
		}
	}
	// finite numeric sum.
	if s := Summation(k, k, Int(1), Int(100)); !s.Equal(Int(5050)) {
		t.Errorf("sum 1..100 = %s, want 5050", s)
	}
	// product k = n!.
	if p := Product(k, k, Int(1), n); !p.Equal(Factorial(n)) {
		t.Errorf("prod k = %s, want n!", p)
	}
	if p := Product(k, k, Int(1), Int(5)); !p.Equal(Int(120)) {
		t.Errorf("prod 1..5 = %s, want 120", p)
	}
}

func TestParseNewFunctions(t *testing.T) {
	for _, s := range []string{
		"sec(x)", "asin(x)", "sinh(x)", "atan2(y, x)", "beta(a, b)",
		"gamma(x)", "erf(x)", "floor(x)", "conjugate(x)",
	} {
		if _, err := Parse(s); err != nil {
			t.Errorf("Parse(%q) error: %v", s, err)
		}
	}
	// Postfix factorial.
	if got := MustParse("5!"); !got.Equal(Int(120)) {
		t.Errorf("5! = %s", got)
	}
	if got := MustParse("(x+1)!"); !got.Equal(Factorial(Add(Sym("x"), Int(1)))) {
		t.Errorf("(x+1)! = %s", got)
	}
	// Imaginary unit constant parses.
	if got := MustParse("I^2"); !got.Equal(Int(-1)) {
		t.Errorf("I^2 = %s", got)
	}
}

// --- examples --------------------------------------------------------------

func ExampleSeries() {
	x := Sym("x")
	fmt.Println(Series(MustParse("exp(x)"), x, Int(0), 5))
	// Output: 1/24*x^4 + 1/6*x^3 + 1/2*x^2 + x + 1
}

func ExampleLimit() {
	x := Sym("x")
	fmt.Println(Limit(MustParse("sin(x)/x"), x, Int(0)))
	// Output: 1
}

func ExampleSummation() {
	k, n := Sym("k"), Sym("n")
	fmt.Println(Summation(k, k, Int(1), n))
	// Output: 1/2*n^2 + 1/2*n
}

func ExampleSolveSystem() {
	x, y := Sym("x"), Sym("y")
	sol, _ := SolveSystem([]Expr{MustParse("x + y - 3"), MustParse("x - y - 1")}, []Expr{x, y})
	fmt.Println(sol[0], sol[1])
	// Output: 2 1
}

func ExampleConjugate() {
	z := Add(Int(2), Mul(Int(3), I))
	fmt.Println(Conjugate(z))
	// Output: 2 - 3*I
}

func ExampleIntegrate_byParts() {
	x := Sym("x")
	fmt.Println(Integrate(MustParse("x*exp(x)"), x))
	// Output: x*exp(x) - exp(x)
}
