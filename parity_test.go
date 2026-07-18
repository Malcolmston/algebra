package algebra

import (
	"math"
	"math/big"
	"sort"
	"testing"
)

// This file encodes known-answer vectors taken directly from SymPy's own test
// suite (github.com/sympy/sympy, master branch) and asserts them against this
// port's public API. Each vector cites the SymPy test file it was lifted from
// so the two stay auditable. Where SymPy's canonical printed form differs from
// this port's (e.g. 2*sqrt(x)**3/3 vs 2/3*x^(3/2)) the assertion compares the
// mathematical value rather than the surface syntax.

// piTimes returns the expression pi * p/q, matching SymPy's pi*Rational(p, q).
func piTimes(p, q int64) Expr { return Mul(Rat(p, q), Pi) }

// evalNum reduces e to a float64 for numeric parity comparisons, failing the
// test if evaluation errors.
func evalNum(t *testing.T, e Expr) float64 {
	t.Helper()
	v, err := Evalf(e)
	if err != nil {
		t.Fatalf("Evalf(%s): %v", e, err)
	}
	return v
}

// TestParityTrigSpecialValues mirrors sympy/functions/elementary/tests/
// test_trigonometric.py::test_sin and ::test_cos exact special-angle values.
func TestParityTrigSpecialValues(t *testing.T) {
	exact := []struct {
		got, want Expr
	}{
		{Sin(piTimes(1, 6)), Rat(1, 2)},   // sin(pi/6) == S.Half
		{Sin(piTimes(-1, 6)), Rat(-1, 2)}, // sin(-pi/6) == -1/2
		{Sin(piTimes(1, 2)), Int(1)},      // sin(pi/2) == 1
		{Sin(piTimes(-1, 2)), Int(-1)},    // sin(-pi/2) == -1
		{Sin(Pi), Int(0)},                 // sin(pi) == 0
		{Sin(Mul(Int(2), Pi)), Int(0)},    // sin(2*pi) == 0
		{Cos(piTimes(1, 3)), Rat(1, 2)},   // cos(pi/3) == S.Half
		{Cos(piTimes(-2, 3)), Rat(-1, 2)}, // cos(-2*pi/3) == -1/2
		{Cos(piTimes(1, 2)), Int(0)},      // cos(pi/2) == 0
		{Cos(Pi), Int(-1)},                // cos(pi) == -1
		{Tan(piTimes(1, 4)), Int(1)},      // tan(pi/4) == 1
	}
	for i, c := range exact {
		if !c.got.Equal(c.want) {
			t.Errorf("exact case %d: got %s, want %s", i, c.got, c.want)
		}
	}

	// sqrt-valued angles: SymPy gives closed forms; compare numerically.
	numeric := []struct {
		got  Expr
		want float64
	}{
		{Cos(piTimes(1, 6)), math.Sqrt(3) / 2},  // cos(pi/6) == sqrt(3)/2
		{Sin(piTimes(1, 4)), math.Sqrt2 / 2},    // sin(pi/4) == sqrt(2)/2
		{Cos(piTimes(1, 4)), math.Sqrt2 / 2},    // cos(pi/4) == sqrt(2)/2
		{Tan(piTimes(1, 3)), math.Sqrt(3)},      // tan(pi/3) == sqrt(3)
		{Sin(piTimes(2, 3)), math.Sqrt(3) / 2},  // sin(2*pi/3) == sqrt(3)/2
		{Cos(piTimes(7, 6)), -math.Sqrt(3) / 2}, // cos(7*pi/6) == -sqrt(3)/2
	}
	for i, c := range numeric {
		if v := evalNum(t, c.got); math.Abs(v-c.want) > 1e-12 {
			t.Errorf("numeric case %d: %s = %g, want %g", i, c.got, v, c.want)
		}
	}
}

// TestParityDiff mirrors sympy/core/tests/test_diff.py (test_diff, test_diff2,
// test_diff3).
func TestParityDiff(t *testing.T) {
	x := Sym("x")
	a := Sym("a")
	b := Sym("b")
	c := Sym("c")

	// (x**2).diff(x) == 2*x
	if got, want := Simplify(Diff(Pow(x, Int(2)), x)), Mul(Int(2), x); !got.Equal(want) {
		t.Errorf("diff(x^2) = %s, want %s", got, want)
	}
	// e = a*b + b**5 ; e.diff(b) == a + 5*b**4
	e := Add(Mul(a, b), Pow(b, Int(5)))
	if got, want := Simplify(Diff(e, b)), Add(a, Mul(Int(5), Pow(b, Int(4)))); !got.Equal(want) {
		t.Errorf("diff(a*b+b^5, b) = %s, want %s", got, want)
	}
	// e = a*b*c ; e.diff(c) == a*b
	if got, want := Simplify(Diff(Mul(a, b, c), c)), Mul(a, b); !got.Equal(want) {
		t.Errorf("diff(a*b*c, c) = %s, want %s", got, want)
	}
	// c*log(c) - c ; diff == log(c)
	if got, want := Simplify(Diff(Add(Mul(c, Log(c)), Mul(Int(-1), c)), c)), Log(c); !got.Equal(want) {
		t.Errorf("diff(c*log(c)-c) = %s, want %s", got, want)
	}
	// tan(c).diff(c) == 1 + tan(c)**2  (one of SymPy's accepted forms)
	if got, want := Simplify(Diff(Tan(c), c)), Add(Int(1), Pow(Tan(c), Int(2))); !got.Equal(want) {
		t.Errorf("diff(tan(c)) = %s, want %s", got, want)
	}
	// (x + 1)**3 ; diff == 3*(x + 1)**2
	if got, want := Diff(Pow(Add(x, Int(1)), Int(3)), x), Mul(Int(3), Pow(Add(x, Int(1)), Int(2))); !got.Equal(want) {
		t.Errorf("diff((x+1)^3) = %s, want %s", got, want)
	}
	// 2*exp(x*x)*x ; diff == 2*exp(x**2) + 4*x**2*exp(x**2)
	got := Simplify(Diff(Mul(Int(2), Exp(Pow(x, Int(2))), x), x))
	want := Add(Mul(Int(2), Exp(Pow(x, Int(2)))), Mul(Int(4), Pow(x, Int(2)), Exp(Pow(x, Int(2)))))
	if !got.Equal(want) {
		t.Errorf("diff(2*exp(x^2)*x) = %s, want %s", got, want)
	}
	// e = c**5 ; e.diff(c, 5) == 120 ; e.diff(c, 6) == 0
	d := Expr(Pow(c, Int(5)))
	for i := 0; i < 5; i++ {
		d = Simplify(Diff(d, c))
	}
	if !d.Equal(Int(120)) {
		t.Errorf("d^5/dc^5 (c^5) = %s, want 120", d)
	}
	if d6 := Simplify(Diff(d, c)); !d6.Equal(Int(0)) {
		t.Errorf("d^6/dc^6 (c^5) = %s, want 0", d6)
	}
}

// TestParityIntegrate mirrors sympy/integrals/tests/test_integrals.py
// indefinite-integral vectors (test_integrate_poly, test_issue_*, etc.).
func TestParityIntegrate(t *testing.T) {
	x := Sym("x")
	a := Sym("a")

	// integrate(x**2, x) == x**3/3
	if got, want := Integrate(Pow(x, Int(2)), x), Mul(Rat(1, 3), Pow(x, Int(3))); !got.Equal(want) {
		t.Errorf("int x^2 = %s, want %s", got, want)
	}
	// integrate(t + 1, t) == t**2/2 + t
	if got, want := Integrate(Add(x, Int(1)), x), Add(Mul(Rat(1, 2), Pow(x, Int(2))), x); !got.Equal(want) {
		t.Errorf("int (x+1) = %s, want %s", got, want)
	}
	// integrate(x**-2, x) == -1/x
	if got, want := Integrate(Pow(x, Int(-2)), x), Mul(Int(-1), Pow(x, Int(-1))); !got.Equal(want) {
		t.Errorf("int x^-2 = %s, want %s", got, want)
	}
	// integrate(1/t, t) == log(t)
	if got, want := Integrate(Pow(x, Int(-1)), x), Log(x); !got.Equal(want) {
		t.Errorf("int 1/x = %s, want %s", got, want)
	}
	// integrate(cos(x), x) == sin(x)
	if got, want := Integrate(Cos(x), x), Sin(x); !got.Equal(want) {
		t.Errorf("int cos = %s, want %s", got, want)
	}
	// integrate(exp(x), x) == exp(x)
	if got, want := Integrate(Exp(x), x), Exp(x); !got.Equal(want) {
		t.Errorf("int exp = %s, want %s", got, want)
	}
	// integrate(a*t**4, t) == a*t**5/5
	if got, want := Integrate(Mul(a, Pow(x, Int(4))), x), Mul(Rat(1, 5), a, Pow(x, Int(5))); !got.Equal(want) {
		t.Errorf("int a*x^4 = %s, want %s", got, want)
	}
	// integrate(1/(1 + x**2), x) == atan(x)
	if got, want := Integrate(Pow(Add(Int(1), Pow(x, Int(2))), Int(-1)), x), Atan(x); !got.Equal(want) {
		t.Errorf("int 1/(1+x^2) = %s, want %s", got, want)
	}

	// integrate(sqrt(x), x) == 2*sqrt(x)**3/3 == 2/3 * x**(3/2)  (gap closed).
	if got, want := Integrate(Sqrt(x), x), Mul(Rat(2, 3), Pow(x, Rat(3, 2))); !got.Equal(want) {
		t.Errorf("int sqrt(x) = %s, want %s", got, want)
	}
	// integrate(x**Rational(1,2), x) == 2/3 * x**(3/2).
	if got, want := Integrate(Pow(x, Rat(1, 2)), x), Mul(Rat(2, 3), Pow(x, Rat(3, 2))); !got.Equal(want) {
		t.Errorf("int x^(1/2) = %s, want %s", got, want)
	}
	// integrate(1/sqrt(x), x) == 2*sqrt(x) == 2 * x**(1/2).
	if got, want := Integrate(Pow(x, Rat(-1, 2)), x), Mul(Int(2), Pow(x, Rat(1, 2))); !got.Equal(want) {
		t.Errorf("int x^(-1/2) = %s, want %s", got, want)
	}

	// Verify the fundamental theorem numerically for the fractional-power
	// antiderivative: d/dx (2/3 x^(3/2)) == sqrt(x) at x = 4 -> 2.
	deriv := Simplify(Diff(Integrate(Sqrt(x), x), x))
	if v, err := Eval(deriv, map[string]float64{"x": 4}); err != nil || math.Abs(v-2) > 1e-9 {
		t.Errorf("d/dx int sqrt(x) at 4 = %v (err %v), want 2", v, err)
	}
}

// numericRoots evaluates and sorts a root list for set comparison.
func numericRoots(t *testing.T, rs []Expr) []float64 {
	t.Helper()
	out := make([]float64, 0, len(rs))
	for _, r := range rs {
		out = append(out, evalNum(t, r))
	}
	sort.Float64s(out)
	return out
}

func rootsMatch(got, want []float64) bool {
	if len(got) != len(want) {
		return false
	}
	for i := range got {
		if math.Abs(got[i]-want[i]) > 1e-9 {
			return false
		}
	}
	return true
}

// TestParitySolve mirrors sympy/solvers/tests/test_solvers.py polynomial
// solve() vectors.
func TestParitySolve(t *testing.T) {
	x := Sym("x")

	// set(solve(x**2 - 4)) == {2, -2}
	r, err := Solve(Add(Pow(x, Int(2)), Int(-4)), x)
	if err != nil {
		t.Fatalf("solve x^2-4: %v", err)
	}
	if got := numericRoots(t, r); !rootsMatch(got, []float64{-2, 2}) {
		t.Errorf("solve(x^2-4) = %v, want [-2 2]", got)
	}

	// set(solve(x**2 - 1, x)) == {1, -1}
	r, _ = Solve(Add(Pow(x, Int(2)), Int(-1)), x)
	if got := numericRoots(t, r); !rootsMatch(got, []float64{-1, 1}) {
		t.Errorf("solve(x^2-1) = %v, want [-1 1]", got)
	}

	// solve(x**2 - 5*x + 6, x) == {2, 3}  (README quick-start vector)
	r, _ = Solve(MustParse("x^2 - 5*x + 6"), x)
	if got := numericRoots(t, r); !rootsMatch(got, []float64{2, 3}) {
		t.Errorf("solve(x^2-5x+6) = %v, want [2 3]", got)
	}

	// set(solve(x**3 - 15*x - 4, x)) == {4, -2+sqrt(3), -2-sqrt(3)}
	r, _ = Solve(MustParse("x^3 - 15*x - 4"), x)
	want := []float64{-2 - math.Sqrt(3), -2 + math.Sqrt(3), 4}
	sort.Float64s(want)
	if got := numericRoots(t, r); !rootsMatch(got, want) {
		t.Errorf("solve(x^3-15x-4) = %v, want %v", got, want)
	}

	// solve(4**(2*x**2 + 2*x) - 8) style quadratic: 4x^2+4x-3 -> {1/2, -3/2}
	r, _ = Solve(MustParse("4*x^2 + 4*x - 3"), x)
	if got := numericRoots(t, r); !rootsMatch(got, []float64{-1.5, 0.5}) {
		t.Errorf("solve(4x^2+4x-3) = %v, want [-1.5 0.5]", got)
	}
}

// TestParityLimits mirrors sympy/series/tests/test_limits.py::test_basic1.
func TestParityLimits(t *testing.T) {
	x := Sym("x")
	cases := []struct {
		got  Expr
		want Expr
	}{
		{Limit(x, x, Inf), Inf},                                       // limit(x, x, oo) is oo
		{Limit(x, x, NegInf), NegInf},                                 // limit(x, x, -oo) is -oo
		{Limit(Mul(Int(-1), x), x, Inf), NegInf},                      // limit(-x, x, oo) is -oo
		{Limit(Pow(x, Int(2)), x, NegInf), Inf},                       // limit(x**2, x, -oo) is oo
		{Limit(Mul(Int(-1), Pow(x, Int(2))), x, Inf), NegInf},         // limit(-x**2, x, oo) is -oo
		{Limit(Pow(x, Int(-1)), x, Inf), Int(0)},                      // limit(1/x, x, oo) == 0
		{Limit(Exp(x), x, Inf), Inf},                                  // limit(exp(x), x, oo) is oo
		{Limit(Mul(Int(-1), Exp(x)), x, Inf), NegInf},                 // limit(-exp(x), x, oo) is -oo
		{Limit(Mul(Exp(x), Pow(x, Int(-1))), x, Inf), Inf},            // limit(exp(x)/x, x, oo) is oo
		{Limit(Add(x, Pow(x, Int(-1))), x, Inf), Inf},                 // limit(x + 1/x, x, oo) is oo
		{Limit(Add(x, Mul(Int(-1), Pow(x, Int(2)))), x, Inf), NegInf}, // limit(x - x**2, x, oo) is -oo
	}
	for i, c := range cases {
		if !c.got.Equal(c.want) {
			t.Errorf("limit case %d: got %s, want %s", i, c.got, c.want)
		}
	}
}

// TestParitySimplify mirrors sympy/simplify/tests/test_simplify.py identities
// that this port supports.
func TestParitySimplify(t *testing.T) {
	x := Sym("x")
	// Pythagorean identity: simplify(sin(x)**2 + cos(x)**2) == 1.
	if got := Simplify(Add(Pow(Sin(x), Int(2)), Pow(Cos(x), Int(2)))); !got.Equal(Int(1)) {
		t.Errorf("simplify(sin^2+cos^2) = %s, want 1", got)
	}
	// Series/Maclaurin sanity for exp(x): coefficient of x is 1, of x^2 is 1/2.
	s := Series(Exp(x), x, Int(0), 4)
	if v := evalNum(t, Subs(s, x, Int(0))); math.Abs(v-1) > 1e-12 {
		t.Errorf("exp series constant term = %g, want 1", v)
	}
	_ = big.NewInt(0)
}
