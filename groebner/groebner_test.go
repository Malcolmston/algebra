package groebner

import (
	"fmt"
	"math/big"
	"math/cmplx"
	"sort"
	"testing"
)

func r(a, b int64) *big.Rat { return big.NewRat(a, b) }

func TestMonomialOps(t *testing.T) {
	a := NewMonomial(1, 2, 0)
	b := NewMonomial(0, 1, 3)
	if got := a.Mul(b); !got.Equal(NewMonomial(1, 3, 3)) {
		t.Errorf("Mul = %v", got)
	}
	if got := a.LCM(b); !got.Equal(NewMonomial(1, 2, 3)) {
		t.Errorf("LCM = %v", got)
	}
	if got := a.GCD(b); !got.Equal(NewMonomial(0, 1, 0)) {
		t.Errorf("GCD = %v", got)
	}
	if a.Coprime(b) {
		t.Error("a and b share x2, not coprime")
	}
	if !NewMonomial(1, 0, 0).Coprime(NewMonomial(0, 1, 0)) {
		t.Error("x and y should be coprime")
	}
	if a.Degree() != 3 {
		t.Errorf("Degree = %d", a.Degree())
	}
	q, ok := NewMonomial(1, 1, 0).Div(NewMonomial(2, 3, 0))
	if !ok || !q.Equal(NewMonomial(1, 2, 0)) {
		t.Errorf("Div = %v %v", q, ok)
	}
	if _, ok := NewMonomial(2, 0, 0).Div(NewMonomial(1, 0, 0)); ok {
		t.Error("division should fail")
	}
}

func TestMonomialOrders(t *testing.T) {
	// Cox, Little, O'Shea examples.
	a := NewMonomial(1, 2, 0) // x y^2
	b := NewMonomial(0, 3, 1) // y^3 z
	// lex: compare first coord: 1 vs 0 -> a > b
	if CompareLex(a, b) <= 0 {
		t.Error("lex: a should exceed b")
	}
	// grlex: deg a =3, deg b=4 -> b > a
	if CompareGrlex(a, b) >= 0 {
		t.Error("grlex: b should exceed a")
	}
	// grevlex classic: x^2 z vs x y^2 (deg 3 each); grevlex ranks x^2z > xy^2? Actually
	// grevlex(x^2z, xy^2): last differing var z: x^2z has z=1, xy^2 has z=0; smaller
	// exponent in last var is larger, so xy^2 > x^2z.
	m1 := NewMonomial(2, 0, 1)
	m2 := NewMonomial(1, 2, 0)
	if CompareGrevlex(m1, m2) >= 0 {
		t.Error("grevlex: xy^2 should exceed x^2z")
	}
}

func TestPolyArithmetic(t *testing.T) {
	x := Var(2, 0)
	y := Var(2, 1)
	// (x+y)^2 = x^2 + 2xy + y^2
	p := x.Add(y).Pow(2)
	want := NewPoly(2,
		NewTerm(r(1, 1), NewMonomial(2, 0)),
		NewTerm(r(2, 1), NewMonomial(1, 1)),
		NewTerm(r(1, 1), NewMonomial(0, 2)),
	)
	if !p.Equal(want) {
		t.Errorf("(x+y)^2 = %v, want %v", p, want)
	}
	// (x+y)(x-y) = x^2 - y^2
	q := x.Add(y).Mul(x.Sub(y))
	wantq := NewPoly(2,
		NewTerm(r(1, 1), NewMonomial(2, 0)),
		NewTerm(r(-1, 1), NewMonomial(0, 2)),
	)
	if !q.Equal(wantq) {
		t.Errorf("(x+y)(x-y) = %v", q)
	}
	if p.TotalDegree() != 2 {
		t.Errorf("degree = %d", p.TotalDegree())
	}
	// derivative d/dx (x^2 + 2xy + y^2) = 2x + 2y
	d := p.Derivative(0)
	wantd := x.ScalarMul(r(2, 1)).Add(y.ScalarMul(r(2, 1)))
	if !d.Equal(wantd) {
		t.Errorf("derivative = %v", d)
	}
}

func TestEval(t *testing.T) {
	x := Var(2, 0)
	y := Var(2, 1)
	p := x.Mul(x).Add(y.Mul(y)) // x^2+y^2
	got := p.Eval([]*big.Rat{r(3, 1), r(4, 1)})
	if got.Cmp(r(25, 1)) != 0 {
		t.Errorf("eval = %v, want 25", got)
	}
	c := p.EvalComplex([]complex128{complex(3, 0), complex(4, 0)})
	if cmplx.Abs(c-complex(25, 0)) > 1e-12 {
		t.Errorf("evalComplex = %v", c)
	}
}

func TestDivisionAndSPoly(t *testing.T) {
	x := Var(2, 0)
	y := Var(2, 1)
	// Divide x^2 y + x y^2 + y^2 by [xy - 1, y^2 - 1] (CLO example) in lex.
	f := x.Mul(x).Mul(y).Add(x.Mul(y).Mul(y)).Add(y.Mul(y))
	g1 := x.Mul(y).Sub(One(2))
	g2 := y.Mul(y).Sub(One(2))
	res := MultivariateDivide(f, []Poly{g1, g2}, LexOrder())
	// Reconstruct f = q1 g1 + q2 g2 + r
	recon := res.Quotients[0].Mul(g1).Add(res.Quotients[1].Mul(g2)).Add(res.Remainder)
	if !recon.Equal(f) {
		t.Errorf("division reconstruction failed: %v != %v", recon, f)
	}
	// S-poly leading terms must cancel.
	s := SPolynomial(g1, g2, LexOrder())
	// S(xy-1, y^2-1) = x*(y^2-1)... let's just check it is in the ideal remainder 0 later.
	_ = s
}

func TestGroebnerBasisMembership(t *testing.T) {
	x := Var(2, 0)
	y := Var(2, 1)
	f := x.Mul(x).Add(y.Mul(y)).Sub(One(2)) // circle
	g := x.Sub(y)                           // line
	gens := []Poly{f, g}
	basis := ReducedGroebnerBasis(gens, LexOrder())
	if !IsGroebnerBasis(basis, LexOrder()) {
		t.Error("reduced basis is not a Gröbner basis")
	}
	// f and g are in the ideal.
	if !InIdeal(f, gens, LexOrder()) || !InIdeal(g, gens, LexOrder()) {
		t.Error("generators must be members")
	}
	// x^2+y^2-1 combined... a random ideal element:
	comb := f.Mul(x).Add(g.Mul(y.Add(One(2))))
	if !InIdeal(comb, gens, LexOrder()) {
		t.Error("combination must be a member")
	}
	// 1 is not in the ideal (nontrivial variety exists).
	if InIdeal(One(2), gens, LexOrder()) {
		t.Error("1 should not be in the ideal")
	}
	// Expected reduced GB under lex: {x - y, y^2 - 1/2}
	want := []Poly{
		x.Sub(y),
		y.Mul(y).Sub(Constant(2, r(1, 2))),
	}
	if len(basis) != len(want) {
		t.Fatalf("basis size %d, want %d: %v", len(basis), len(want), basis)
	}
	for i := range want {
		if !basis[i].Equal(want[i]) {
			t.Errorf("basis[%d] = %v, want %v", i, basis[i], want[i])
		}
	}
}

func TestUnitIdeal(t *testing.T) {
	x := Var(1, 0)
	// (x, x+1) = whole ring.
	id := NewIdeal(LexOrder(), x, x.Add(One(1)))
	if !id.IsUnit() {
		t.Error("(x, x+1) should be the unit ideal")
	}
	if !id.Contains(One(1)) {
		t.Error("unit ideal contains 1")
	}
}

func TestIdealIntersection(t *testing.T) {
	x := Var(2, 0)
	y := Var(2, 1)
	ix := NewIdeal(LexOrder(), x)
	iy := NewIdeal(LexOrder(), y)
	inter := ix.Intersect(iy)
	// (x) ∩ (y) = (xy)
	gb := inter.GroebnerBasis()
	if len(gb) != 1 || !gb[0].Equal(x.Mul(y)) {
		t.Errorf("(x)∩(y) = %v, want (xy)", gb)
	}
}

func TestIdealQuotient(t *testing.T) {
	x := Var(2, 0)
	y := Var(2, 1)
	// (xy) : (x) = (y)
	id := NewIdeal(LexOrder(), x.Mul(y))
	q := id.QuotientElem(x)
	gb := q.GroebnerBasis()
	if len(gb) != 1 || !gb[0].Equal(y) {
		t.Errorf("(xy):(x) = %v, want (y)", gb)
	}
}

func TestIdealSumProduct(t *testing.T) {
	x := Var(2, 0)
	y := Var(2, 1)
	ix := NewIdeal(LexOrder(), x)
	iy := NewIdeal(LexOrder(), y)
	sum := ix.Sum(iy)
	if !sum.Contains(x) || !sum.Contains(y) || !sum.Contains(x.Add(y)) {
		t.Error("sum should contain x, y, x+y")
	}
	prod := ix.Product(iy)
	gb := prod.GroebnerBasis()
	if len(gb) != 1 || !gb[0].Equal(x.Mul(y)) {
		t.Errorf("(x)(y) = %v, want (xy)", gb)
	}
}

func TestElimination(t *testing.T) {
	// System: x = t, y = t^2 -> eliminate t to get y - x^2? Use vars (t,x,y).
	// t is var0, x var1, y var2. gens: x - t, y - t^2. Eliminate t (first var).
	t0 := Var(3, 0)
	x := Var(3, 1)
	y := Var(3, 2)
	g1 := x.Sub(t0)
	g2 := y.Sub(t0.Mul(t0))
	elim := EliminationIdeal([]Poly{g1, g2}, 1, LexOrder())
	// Expect y - x^2 (up to sign / scaling).
	if len(elim) != 1 {
		t.Fatalf("elimination ideal size %d: %v", len(elim), elim)
	}
	want := y.Sub(x.Mul(x))
	if !elim[0].Equal(want) && !elim[0].Equal(want.Neg()) {
		t.Errorf("elimination = %v, want ±(y - x^2)", elim[0])
	}
}

func TestZeroDimensional(t *testing.T) {
	x := Var(2, 0)
	y := Var(2, 1)
	f := x.Mul(x).Sub(One(2)) // x^2-1
	g := y.Mul(y).Sub(One(2)) // y^2-1
	id := NewIdeal(GrevlexOrder(), f, g)
	if !id.IsZeroDimensional() {
		t.Error("should be zero-dimensional")
	}
	dim, ok := id.VectorSpaceDimension()
	if !ok || dim != 4 {
		t.Errorf("vector space dimension = %d, want 4", dim)
	}
	// Positive-dimensional example: single equation in 2 vars.
	id2 := NewIdeal(GrevlexOrder(), f)
	if id2.IsZeroDimensional() {
		t.Error("single equation in 2 vars is not zero-dimensional")
	}
}

func TestSolveUnivariate(t *testing.T) {
	// x^2 - 5x + 6 = 0 -> roots 2, 3.
	roots := SolveUnivariate([]complex128{complex(6, 0), complex(-5, 0), complex(1, 0)}, 1)
	if len(roots) != 2 {
		t.Fatalf("got %d roots", len(roots))
	}
	vals := []float64{real(roots[0]), real(roots[1])}
	sort.Float64s(vals)
	if abs(vals[0]-2) > 1e-6 || abs(vals[1]-3) > 1e-6 {
		t.Errorf("roots = %v, want {2,3}", vals)
	}
	// x^2 + 1 = 0 -> ±i.
	roots2 := SolveUnivariate([]complex128{complex(1, 0), complex(0, 0), complex(1, 0)}, 3)
	for _, z := range roots2 {
		if abs(cmplx.Abs(z)-1) > 1e-6 || abs(real(z)) > 1e-6 {
			t.Errorf("root %v not ±i", z)
		}
	}
}

func TestSolveZeroDimensional(t *testing.T) {
	x := Var(2, 0)
	y := Var(2, 1)
	// Circle ∩ line: x^2+y^2=1, x=y -> (±1/√2, ±1/√2).
	f := x.Mul(x).Add(y.Mul(y)).Sub(One(2))
	g := x.Sub(y)
	sols, err := SolveZeroDimensional([]Poly{f, g}, 42, 1e-9)
	if err != nil {
		t.Fatal(err)
	}
	if len(sols) != 2 {
		t.Fatalf("got %d solutions, want 2", len(sols))
	}
	inv := 0.7071067811865476
	for _, s := range sols {
		if abs(real(s[0])-real(s[1])) > 1e-6 {
			t.Errorf("x != y in solution %v", s)
		}
		if abs(abs(real(s[0]))-inv) > 1e-6 {
			t.Errorf("|x| = %v, want %v", real(s[0]), inv)
		}
		if abs(imag(s[0])) > 1e-6 {
			t.Errorf("solution not real: %v", s)
		}
	}

	// Four corners: x^2-1, y^2-1 -> 4 real solutions.
	f2 := x.Mul(x).Sub(One(2))
	g2 := y.Mul(y).Sub(One(2))
	sols2, err := SolveZeroDimensional([]Poly{f2, g2}, 7, 1e-9)
	if err != nil {
		t.Fatal(err)
	}
	if len(sols2) != 4 {
		t.Fatalf("got %d solutions, want 4", len(sols2))
	}
	reals := RealSolutions(sols2, 1e-6)
	if len(reals) != 4 {
		t.Errorf("expected 4 real solutions, got %d", len(reals))
	}
	for _, s := range reals {
		if abs(abs(s[0])-1) > 1e-6 || abs(abs(s[1])-1) > 1e-6 {
			t.Errorf("corner solution off unit: %v", s)
		}
	}
}

func TestSolveNoSolutions(t *testing.T) {
	x := Var(1, 0)
	// (x^2+1, ...) unit? No. Use (x, x-1) unit ideal -> no solutions.
	sols, err := SolveZeroDimensional([]Poly{x, x.Sub(One(1))}, 1, 1e-9)
	if err != nil {
		t.Fatal(err)
	}
	if len(sols) != 0 {
		t.Errorf("unit ideal should have no solutions, got %v", sols)
	}
}

func TestReducedBasisCanonical(t *testing.T) {
	x := Var(2, 0)
	y := Var(2, 1)
	gens := []Poly{
		x.Mul(x).Add(y.Mul(y)).Sub(One(2)),
		x.Sub(y),
	}
	// Different generator orderings give the same reduced basis.
	b1 := ReducedGroebnerBasis(gens, LexOrder())
	b2 := ReducedGroebnerBasis([]Poly{gens[1], gens[0]}, LexOrder())
	if len(b1) != len(b2) {
		t.Fatalf("basis sizes differ: %d vs %d", len(b1), len(b2))
	}
	for i := range b1 {
		if !b1[i].Equal(b2[i]) {
			t.Errorf("reduced basis not canonical: %v vs %v", b1[i], b2[i])
		}
	}
}

func TestRingHelpers(t *testing.T) {
	ring := NewRing(3, GrevlexOrder())
	if ring.Nvars() != 3 {
		t.Errorf("nvars = %d", ring.Nvars())
	}
	v := ring.Vars()
	p := v[0].Add(v[1]).Add(v[2])
	if p.TotalDegree() != 1 {
		t.Errorf("degree = %d", p.TotalDegree())
	}
	id := ring.Ideal(v[0], v[1], v[2])
	if !id.Contains(p) {
		t.Error("ideal should contain x+y+z")
	}
}

func abs(x float64) float64 {
	if x < 0 {
		return -x
	}
	return x
}

func ExampleReducedGroebnerBasis() {
	// Compute the reduced Gröbner basis of <x^2 + y^2 - 1, x - y> under lex.
	x := Var(2, 0)
	y := Var(2, 1)
	f := x.Mul(x).Add(y.Mul(y)).Sub(One(2))
	g := x.Sub(y)
	basis := ReducedGroebnerBasis([]Poly{f, g}, LexOrder())
	for _, p := range basis {
		fmt.Println(p)
	}
	// Output:
	// x1 - x2
	// x2^2 - 1/2
}

func ExampleInIdeal() {
	x := Var(2, 0)
	y := Var(2, 1)
	gens := []Poly{x.Mul(x).Sub(y), x.Mul(y).Sub(x)}
	// x^2*y - y^2 = y*(x^2 - y) is in the ideal.
	f := x.Mul(x).Mul(y).Sub(y.Mul(y))
	fmt.Println(InIdeal(f, gens, GrevlexOrder()))
	// Output:
	// true
}
