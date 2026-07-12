package algebra

import (
	"fmt"
	"math"
	"math/big"
	"testing"
)

func bigFromInt64(v int64) *big.Int { return big.NewInt(v) }

func mustParse(t *testing.T, s string) Expr {
	t.Helper()
	e, err := Parse(s)
	if err != nil {
		t.Fatalf("Parse(%q) error: %v", s, err)
	}
	return e
}

// --- parsing ---------------------------------------------------------------

func TestParseRoundTrip(t *testing.T) {
	cases := []struct{ in, want string }{
		{"x^2 + 2*x + 1", "x^2 + 2*x + 1"},
		{"2*x + 3*x", "5*x"},
		{"x*x", "x^2"},
		{"(x+1)*(x-1)", "(x - 1)*(x + 1)"},
		{"2x", "2*x"},
		{"2(x+1)", "2*(x + 1)"},
		{"-x^2", "-x^2"},
		{"sin(x) + cos(x)", "cos(x) + sin(x)"},
		{"3/4", "3/4"},
		{"a/b", "a*b^(-1)"},
		{"2^10", "1024"},
	}
	for _, c := range cases {
		got := mustParse(t, c.in).String()
		if got != c.want {
			t.Errorf("Parse(%q).String() = %q, want %q", c.in, got, c.want)
		}
	}
}

func TestParseReparseStable(t *testing.T) {
	inputs := []string{
		"x^2 + 2*x + 1", "sin(x^2)*cos(x)", "(a+b)^3", "1/2*x + 3", "exp(x) - log(x)",
	}
	for _, s := range inputs {
		e1 := mustParse(t, s)
		e2 := mustParse(t, e1.String())
		if !e1.Equal(e2) {
			t.Errorf("reparse of %q not stable: %q vs %q", s, e1, e2)
		}
	}
}

func TestParseErrors(t *testing.T) {
	for _, s := range []string{"x +", "(x", "x)", "2 @ 3"} {
		if _, err := Parse(s); err == nil {
			t.Errorf("Parse(%q) expected error, got nil", s)
		}
	}
}

// --- simplification identities --------------------------------------------

func TestIdentities(t *testing.T) {
	x := Sym("x")
	cases := []struct {
		got  Expr
		want Expr
	}{
		{Add(x, Int(0)), x},
		{Mul(x, Int(1)), x},
		{Mul(x, Int(0)), Int(0)},
		{Pow(x, Int(1)), x},
		{Pow(x, Int(0)), Int(1)},
		{Add(x, x), Mul(Int(2), x)},
		{Mul(x, x), Pow(x, Int(2))},
		{Add(Int(2), Int(3)), Int(5)},
		{Mul(Rat(1, 2), Int(2)), Int(1)},
		{Add(x, Mul(Int(2), x)), Mul(Int(3), x)},
		{Add(Mul(Int(2), x), Mul(Int(-2), x)), Int(0)},
		{Pow(Pow(x, Int(2)), Int(3)), Pow(x, Int(6))},
	}
	for i, c := range cases {
		if !c.got.Equal(c.want) {
			t.Errorf("case %d: got %s, want %s", i, c.got, c.want)
		}
	}
}

func TestSimplifyFoldsFunctions(t *testing.T) {
	cases := []struct {
		got  Expr
		want Expr
	}{
		{Sin(Int(0)), Int(0)},
		{Cos(Int(0)), Int(1)},
		{Exp(Int(0)), Int(1)},
		{Log(Int(1)), Int(0)},
		{Log(E), Int(1)},
		{Sqrt(Int(9)), Int(3)},
		{Sqrt(Int(8)), Mul(Int(2), Sqrt(Int(2)))},
	}
	for i, c := range cases {
		if !c.got.Equal(c.want) {
			t.Errorf("case %d: got %s, want %s", i, c.got, c.want)
		}
	}
}

// --- differentiation -------------------------------------------------------

func TestDiffRules(t *testing.T) {
	x := Sym("x")
	cases := []struct {
		in   string
		want string
	}{
		{"x^2", "2*x"},
		{"x^3", "3*x^2"},
		{"3*x + 5", "3"},
		{"x*sin(x)", "x*cos(x) + sin(x)"}, // product rule
		{"sin(x^2)", "2*x*cos(x^2)"},      // chain rule
		{"sin(x)", "cos(x)"},
		{"cos(x)", "-sin(x)"},
		{"exp(x)", "exp(x)"},
		{"log(x)", "x^(-1)"},
		{"exp(x^2)", "2*x*exp(x^2)"},
		{"1/x", "-x^(-2)"}, // quotient/power rule
	}
	for _, c := range cases {
		got := Diff(mustParse(t, c.in), x).String()
		if got != c.want {
			t.Errorf("d/dx %s = %s, want %s", c.in, got, c.want)
		}
	}
}

// TestDiffNumeric cross-checks symbolic derivatives against a finite-difference
// approximation at several sample points.
func TestDiffNumeric(t *testing.T) {
	x := Sym("x")
	exprs := []string{
		"x^3 + 2*x", "sin(x)*cos(x)", "exp(x)*x", "x^2*log(x)",
		"tan(x)", "sqrt(x)", "1/(x^2 + 1)", "sin(x^2)",
	}
	pts := []float64{0.3, 0.7, 1.2, 2.1}
	const h = 1e-6
	for _, s := range exprs {
		e := mustParse(t, s)
		d := Diff(e, x)
		for _, p := range pts {
			fp, err1 := Eval(e, map[string]float64{"x": p + h})
			fm, err2 := Eval(e, map[string]float64{"x": p - h})
			da, err3 := Eval(d, map[string]float64{"x": p})
			if err1 != nil || err2 != nil || err3 != nil {
				t.Fatalf("%s eval error at %g", s, p)
			}
			approx := (fp - fm) / (2 * h)
			if math.Abs(approx-da) > 1e-4*(1+math.Abs(approx)) {
				t.Errorf("d/dx %s at %g: symbolic %g vs numeric %g", s, p, da, approx)
			}
		}
	}
}

// --- expansion -------------------------------------------------------------

func TestExpand(t *testing.T) {
	cases := []struct{ in, want string }{
		{"(x+1)^2", "x^2 + 2*x + 1"},
		{"(x+1)^3", "x^3 + 3*x^2 + 3*x + 1"},
		{"(x+1)*(x-1)", "x^2 - 1"},
		{"(x+2)*(x+3)", "x^2 + 5*x + 6"},
		{"(a+b)^2", "a^2 + 2*a*b + b^2"},
		{"2*(x+3)", "2*x + 6"},
	}
	for _, c := range cases {
		got := Expand(mustParse(t, c.in)).String()
		if got != c.want {
			t.Errorf("Expand(%s) = %s, want %s", c.in, got, c.want)
		}
	}
}

// --- substitution ----------------------------------------------------------

func TestSubs(t *testing.T) {
	x, y := Sym("x"), Sym("y")
	got := Subs(mustParse(t, "x^2 + 1"), x, Int(3))
	if !got.Equal(Int(10)) {
		t.Errorf("subs x=3 in x^2+1 = %s, want 10", got)
	}
	got2 := Subs(mustParse(t, "y^2"), y, mustParse(t, "x+1")).Expand()
	if got2.String() != "x^2 + 2*x + 1" {
		t.Errorf("subs y=x+1 in y^2 = %s", got2)
	}
}

// --- numeric evaluation ----------------------------------------------------

func TestEval(t *testing.T) {
	cases := []struct {
		in  string
		env map[string]float64
		val float64
	}{
		{"x^2 + 2*x + 1", map[string]float64{"x": 3}, 16},
		{"sin(pi)", nil, math.Sin(math.Pi)},
		{"exp(0)", nil, 1},
		{"log(E)", nil, 1},
		{"sqrt(2)", nil, math.Sqrt2},
		{"2*x + y", map[string]float64{"x": 1, "y": 5}, 7},
	}
	for _, c := range cases {
		v, err := Eval(mustParse(t, c.in), c.env)
		if err != nil {
			t.Errorf("Eval(%s) error: %v", c.in, err)
			continue
		}
		if math.Abs(v-c.val) > 1e-9 {
			t.Errorf("Eval(%s) = %g, want %g", c.in, v, c.val)
		}
	}
	if _, err := Eval(Sym("z"), nil); err == nil {
		t.Error("Eval of unbound symbol should error")
	}
}

// --- integration -----------------------------------------------------------

func TestIntegrate(t *testing.T) {
	x := Sym("x")
	cases := []struct{ in, want string }{
		{"x^2", "1/3*x^3"},
		{"x", "1/2*x^2"},
		{"1/x", "log(x)"},
		{"exp(x)", "exp(x)"},
		{"sin(x)", "-cos(x)"},
		{"cos(x)", "sin(x)"},
		{"3", "3*x"},
		{"cos(2*x)", "1/2*sin(2*x)"},
		{"2*x + 1", "x^2 + x"},
	}
	for _, c := range cases {
		got := Integrate(mustParse(t, c.in), x).String()
		if got != c.want {
			t.Errorf("Integrate(%s) = %s, want %s", c.in, got, c.want)
		}
	}
	// Unhandled integrand returns an unevaluated Integral node.
	un := Integrate(mustParse(t, "x*sin(x)"), x)
	if _, ok := un.(*integral); !ok {
		t.Errorf("expected unevaluated integral, got %s", un)
	}
}

// TestIntegrateInverseOfDiff verifies that differentiating an antiderivative
// recovers the integrand numerically for the covered cases.
func TestIntegrateInverseOfDiff(t *testing.T) {
	x := Sym("x")
	for _, s := range []string{"x^2", "x^3 + x", "exp(x)", "sin(x)", "cos(2*x)"} {
		f := mustParse(t, s)
		F := Integrate(f, x)
		back := Simplify(Diff(F, x))
		for _, p := range []float64{0.4, 1.3, 2.2} {
			a, _ := Eval(f, map[string]float64{"x": p})
			b, _ := Eval(back, map[string]float64{"x": p})
			if math.Abs(a-b) > 1e-6*(1+math.Abs(a)) {
				t.Errorf("d/dx int %s at %g: %g vs %g", s, p, b, a)
			}
		}
	}
}

// --- solving ---------------------------------------------------------------

func solve(t *testing.T, s string) []Expr {
	t.Helper()
	r, err := Solve(mustParse(t, s), Sym("x"))
	if err != nil {
		t.Fatalf("Solve(%s) error: %v", s, err)
	}
	return r
}

func TestSolveLinear(t *testing.T) {
	r := solve(t, "2*x + 4")
	if len(r) != 1 || !r[0].Equal(Int(-2)) {
		t.Errorf("Solve(2x+4) = %v, want [-2]", r)
	}
	r = solve(t, "3*x - 9")
	if len(r) != 1 || !r[0].Equal(Int(3)) {
		t.Errorf("Solve(3x-9) = %v, want [3]", r)
	}
}

func TestSolveQuadratic(t *testing.T) {
	// Two distinct rational roots.
	r := solve(t, "x^2 - 5*x + 6")
	if len(r) != 2 || !r[0].Equal(Int(2)) || !r[1].Equal(Int(3)) {
		t.Errorf("Solve(x^2-5x+6) = %v, want [2 3]", r)
	}
	// Repeated root.
	r = solve(t, "x^2 - 4*x + 4")
	if len(r) != 1 || !r[0].Equal(Int(2)) {
		t.Errorf("Solve(x^2-4x+4) = %v, want [2]", r)
	}
	// Irrational roots via sqrt.
	r = solve(t, "x^2 - 2")
	if len(r) != 2 {
		t.Fatalf("Solve(x^2-2) = %v, want 2 roots", r)
	}
	// Verify each root satisfies the equation numerically.
	for _, root := range r {
		v, err := Evalf(root)
		if err != nil {
			t.Fatalf("eval root %s: %v", root, err)
		}
		if math.Abs(v*v-2) > 1e-9 {
			t.Errorf("root %s = %g does not satisfy x^2=2", root, v)
		}
	}
}

func TestSolveErrors(t *testing.T) {
	if _, err := Solve(mustParse(t, "x^3 + 1"), Sym("x")); err == nil {
		t.Error("cubic should be unsupported")
	}
	if _, err := Solve(mustParse(t, "sin(x)"), Sym("x")); err == nil {
		t.Error("non-polynomial should error")
	}
}

// --- bonus: factor / collect ----------------------------------------------

func TestFactorCollect(t *testing.T) {
	x := Sym("x")
	f := Factor(mustParse(t, "x^2 - 1"), x)
	// Re-expanding the factorization must recover the polynomial.
	if !Expand(f).Equal(mustParse(t, "x^2 - 1")) {
		t.Errorf("Factor(x^2-1) = %s does not expand back", f)
	}
	c := Collect(mustParse(t, "x + 3*x^2 + 2 + x"), x)
	if c.String() != "3*x^2 + 2*x + 2" {
		t.Errorf("Collect = %s, want 3*x^2 + 2*x + 2", c)
	}
}

// --- equality --------------------------------------------------------------

func TestEqual(t *testing.T) {
	if !mustParse(t, "x + y").Equal(mustParse(t, "y + x")) {
		t.Error("x+y should equal y+x")
	}
	if mustParse(t, "x + y").Equal(mustParse(t, "x + z")) {
		t.Error("x+y should not equal x+z")
	}
	if !Int(2).Equal(Add(Int(1), Int(1))) {
		t.Error("2 should equal 1+1")
	}
}

// --- fluent builder API ----------------------------------------------------

func TestBuilderMethods(t *testing.T) {
	x := Sym("x")
	// x^2 + 2*x + 1 assembled fluently.
	e := x.Pow(Int(2)).Add(x.Mul(Int(2)), Int(1))
	if e.String() != "x^2 + 2*x + 1" {
		t.Errorf("builder built %s", e)
	}
	if !e.Diff(x).Equal(mustParse(t, "2*x + 2")) {
		t.Errorf("e.Diff = %s", e.Diff(x))
	}
	if !x.Add(x).Simplify().Equal(Mul(Int(2), x)) {
		t.Error("x.Add(x) should simplify to 2*x")
	}
	if !mustParse(t, "(x+1)^2").Expand().Equal(mustParse(t, "x^2 + 2*x + 1")) {
		t.Error("Expand via method failed")
	}
	if !mustParse(t, "y^2").Subs(Sym("y"), Int(2)).Equal(Int(4)) {
		t.Error("Subs via method failed")
	}
}

func TestConstructorsAndPrinting(t *testing.T) {
	if Pi.String() != "pi" || E.String() != "E" {
		t.Errorf("constant printing: %s %s", Pi, E)
	}
	if IntBig(bigFromInt64(42)).String() != "42" {
		t.Error("IntBig failed")
	}
	if Rat(2, 4).String() != "1/2" {
		t.Errorf("Rat(2,4) = %s", Rat(2, 4))
	}
	if Flt(1.5).String() != "1.5" {
		t.Errorf("Flt = %s", Flt(1.5))
	}
	v, err := Evalf(mustParse(t, "2^10"))
	if err != nil || v != 1024 {
		t.Errorf("Evalf(2^10) = %g, %v", v, err)
	}
}

func TestRatPanicsOnZero(t *testing.T) {
	defer func() {
		if recover() == nil {
			t.Error("Rat(1,0) should panic")
		}
	}()
	Rat(1, 0)
}

func TestMustParsePanics(t *testing.T) {
	defer func() {
		if recover() == nil {
			t.Error("MustParse of garbage should panic")
		}
	}()
	MustParse("(((")
}

// --- godoc examples --------------------------------------------------------

func ExampleParse() {
	e, _ := Parse("x^2 + 2*x + 1")
	fmt.Println(e)
	// Output: x^2 + 2*x + 1
}

func ExampleDiff() {
	x := Sym("x")
	fmt.Println(Diff(MustParse("x^2"), x))
	fmt.Println(Diff(MustParse("sin(x^2)"), x))
	// Output:
	// 2*x
	// 2*x*cos(x^2)
}

func ExampleExpand() {
	fmt.Println(Expand(MustParse("(x + 1)^2")))
	// Output: x^2 + 2*x + 1
}

func ExampleSolve() {
	roots, _ := Solve(MustParse("x^2 - 5*x + 6"), Sym("x"))
	fmt.Println(roots[0], roots[1])
	// Output: 2 3
}

func ExampleSimplify() {
	x := Sym("x")
	fmt.Println(Simplify(Add(x, x, Int(3), Int(4))))
	// Output: 2*x + 7
}
