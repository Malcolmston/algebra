package rootfind

import (
	"fmt"
	"math"
	"math/cmplx"
	"sort"
	"testing"
)

// --- helpers ---------------------------------------------------------------

func approx(a, b, tol float64) bool { return math.Abs(a-b) <= tol }

func approxC(a, b complex128, tol float64) bool { return cmplx.Abs(a-b) <= tol }

// matchReal checks that got and want contain the same real values (order
// independent) within tol.
func matchReal(got, want []float64, tol float64) bool {
	if len(got) != len(want) {
		return false
	}
	g := append([]float64(nil), got...)
	w := append([]float64(nil), want...)
	sort.Float64s(g)
	sort.Float64s(w)
	for i := range g {
		if !approx(g[i], w[i], tol) {
			return false
		}
	}
	return true
}

// matchComplex checks that got contains each wanted complex value within tol,
// order independent.
func matchComplex(got, want []complex128, tol float64) bool {
	if len(got) != len(want) {
		return false
	}
	used := make([]bool, len(got))
	for _, w := range want {
		found := false
		for i, g := range got {
			if !used[i] && approxC(g, w, tol) {
				used[i] = true
				found = true
				break
			}
		}
		if !found {
			return false
		}
	}
	return true
}

// --- polynomial arithmetic -------------------------------------------------

func TestPolyEvalAndDeriv(t *testing.T) {
	// p(x) = 1 + 2x + 3x^2
	p := NewPoly(1, 2, 3)
	tests := []struct {
		x, val, d1, d2 float64
	}{
		{0, 1, 2, 6},
		{1, 6, 8, 6},
		{2, 17, 14, 6},
		{-1, 2, -4, 6},
	}
	for _, tc := range tests {
		if v := p.Eval(tc.x); !approx(v, tc.val, 1e-12) {
			t.Errorf("Eval(%v)=%v want %v", tc.x, v, tc.val)
		}
		v, d := p.EvalDeriv(tc.x)
		if !approx(v, tc.val, 1e-12) || !approx(d, tc.d1, 1e-12) {
			t.Errorf("EvalDeriv(%v)=(%v,%v) want (%v,%v)", tc.x, v, d, tc.val, tc.d1)
		}
		v2, d1b, d2 := p.EvalDeriv2(tc.x)
		if !approx(v2, tc.val, 1e-12) || !approx(d1b, tc.d1, 1e-12) || !approx(d2, tc.d2, 1e-12) {
			t.Errorf("EvalDeriv2(%v)=(%v,%v,%v) want (%v,%v,%v)", tc.x, v2, d1b, d2, tc.val, tc.d1, tc.d2)
		}
	}
}

func TestPolyDegreeTrim(t *testing.T) {
	tests := []struct {
		p    Poly
		deg  int
		zero bool
	}{
		{NewPoly(0, 0, 0), -1, true},
		{NewPoly(5), 0, false},
		{NewPoly(1, 2, 0, 0), 1, false},
		{NewPoly(0, 0, 7), 2, false},
	}
	for _, tc := range tests {
		if d := tc.p.Degree(); d != tc.deg {
			t.Errorf("Degree(%v)=%d want %d", tc.p, d, tc.deg)
		}
		if z := tc.p.IsZero(); z != tc.zero {
			t.Errorf("IsZero(%v)=%v want %v", tc.p, z, tc.zero)
		}
	}
}

func TestPolyArithmetic(t *testing.T) {
	a := NewPoly(1, 2, 3) // 3x^2+2x+1
	b := NewPoly(0, 1, 0, 4)
	if s := a.Add(b); !s.Equal(NewPoly(1, 3, 3, 4), 1e-12) {
		t.Errorf("Add=%v", s)
	}
	if d := b.Sub(a); !d.Equal(NewPoly(-1, -1, -3, 4), 1e-12) {
		t.Errorf("Sub=%v", d)
	}
	// (x+1)(x+2) = x^2+3x+2
	m := NewPoly(1, 1).Mul(NewPoly(2, 1))
	if !m.Equal(NewPoly(2, 3, 1), 1e-12) {
		t.Errorf("Mul=%v", m)
	}
	if sc := a.Scale(2); !sc.Equal(NewPoly(2, 4, 6), 1e-12) {
		t.Errorf("Scale=%v", sc)
	}
	if p := NewPoly(1, 1).Pow(3); !p.Equal(NewPoly(1, 3, 3, 1), 1e-12) {
		t.Errorf("Pow=%v", p)
	}
}

func TestPolyDivMod(t *testing.T) {
	// (x^3 - 6x^2 + 11x - 6) / (x - 1) = x^2 - 5x + 6, r=0
	p := FromRoots(1, 2, 3)
	q, r, err := p.DivMod(NewPoly(-1, 1))
	if err != nil {
		t.Fatal(err)
	}
	if !q.Equal(NewPoly(6, -5, 1), 1e-9) {
		t.Errorf("quo=%v", q)
	}
	if r.Degree() >= 0 && math.Abs(r.Eval(0)) > 1e-9 {
		t.Errorf("rem=%v want 0", r)
	}
	// reconstruct
	recon := q.Mul(NewPoly(-1, 1)).Add(r)
	if !recon.Equal(p, 1e-9) {
		t.Errorf("q*d+r=%v want %v", recon, p)
	}
}

func TestPolyGCD(t *testing.T) {
	// gcd((x-1)^2(x-2), (x-1)(x-3)) = (x-1)
	a := FromRoots(1, 1, 2)
	b := FromRoots(1, 3)
	g := a.GCD(b)
	if !approx(g.Eval(1), 0, 1e-7) {
		t.Errorf("gcd(1)=%v want 0; gcd=%v", g.Eval(1), g)
	}
	if g.Degree() != 1 {
		t.Errorf("gcd degree=%d want 1 (%v)", g.Degree(), g)
	}
}

func TestPolyComposeReverse(t *testing.T) {
	p := NewPoly(1, 0, 1) // x^2+1
	// p(x+1) = (x+1)^2+1 = x^2+2x+2
	c := p.Compose(NewPoly(1, 1))
	if !c.Equal(NewPoly(2, 2, 1), 1e-12) {
		t.Errorf("compose=%v", c)
	}
	// reverse of 1+2x+3x^2 is 3+2x+x^2
	r := NewPoly(1, 2, 3).Reverse()
	if !r.Equal(NewPoly(3, 2, 1), 1e-12) {
		t.Errorf("reverse=%v", r)
	}
}

func TestFromRootsVieta(t *testing.T) {
	p := FromRoots(1, 2, 3)
	if !p.Equal(NewPoly(-6, 11, -6, 1), 1e-12) {
		t.Errorf("FromRoots=%v", p)
	}
	if !approx(p.SumOfRoots(), 6, 1e-12) {
		t.Errorf("SumOfRoots=%v want 6", p.SumOfRoots())
	}
	if !approx(p.ProductOfRoots(), 6, 1e-12) {
		t.Errorf("ProductOfRoots=%v want 6", p.ProductOfRoots())
	}
}

func TestPolyString(t *testing.T) {
	tests := []struct {
		p    Poly
		want string
	}{
		{NewPoly(0), "0"},
		{NewPoly(-6, 11, -6, 1), "x^3 - 6x^2 + 11x - 6"},
		{NewPoly(5), "5"},
		{NewPoly(0, -1), "-x"},
	}
	for _, tc := range tests {
		if s := tc.p.String(); s != tc.want {
			t.Errorf("String(%v)=%q want %q", []float64(tc.p), s, tc.want)
		}
	}
}

// --- scalar solvers --------------------------------------------------------

func TestScalarSolvers(t *testing.T) {
	// root of cos(x) - x is the Dottie number ~0.7390851332151607
	const want = 0.7390851332151607
	f := func(x float64) float64 { return math.Cos(x) - x }
	df := func(x float64) float64 { return -math.Sin(x) - 1 }
	d2f := func(x float64) float64 { return -math.Cos(x) }

	type solver struct {
		name string
		run  func() (Result, error)
	}
	solvers := []solver{
		{"Bisection", func() (Result, error) { return Bisection(f, 0, 1, 1e-14, 200) }},
		{"FalsePosition", func() (Result, error) { return FalsePosition(f, 0, 1, 1e-14, 200) }},
		{"Illinois", func() (Result, error) { return Illinois(f, 0, 1, 1e-14, 200) }},
		{"Brent", func() (Result, error) { return Brent(f, 0, 1, 1e-14, 200) }},
		{"Ridders", func() (Result, error) { return Ridders(f, 0, 1, 1e-14, 200) }},
		{"Secant", func() (Result, error) { return Secant(f, 0, 1, 1e-14, 200) }},
		{"Steffensen", func() (Result, error) { return Steffensen(f, 0.5, 1e-14, 200) }},
		{"Newton", func() (Result, error) { return Newton(f, df, 0.5, 1e-14, 200) }},
		{"Halley", func() (Result, error) { return Halley(f, df, d2f, 0.5, 1e-14, 200) }},
	}
	for _, s := range solvers {
		res, err := s.run()
		if err != nil {
			t.Errorf("%s error: %v", s.name, err)
			continue
		}
		if !res.Converged {
			t.Errorf("%s did not converge", s.name)
		}
		if !approx(res.Root, want, 1e-9) {
			t.Errorf("%s root=%.15f want %.15f", s.name, res.Root, want)
		}
	}
}

func TestFixedPoint(t *testing.T) {
	// x = cos(x) fixed point == Dottie number
	res, err := FixedPoint(math.Cos, 0.5, 1e-12, 500)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(res.Root, 0.7390851332151607, 1e-8) {
		t.Errorf("FixedPoint=%v", res.Root)
	}
}

func TestNoBracket(t *testing.T) {
	f := func(x float64) float64 { return x*x + 1 } // no real root
	if _, err := Bisection(f, -1, 1, 1e-12, 50); err != ErrNoBracket {
		t.Errorf("want ErrNoBracket, got %v", err)
	}
	if _, err := Brent(f, -1, 1, 1e-12, 50); err != ErrNoBracket {
		t.Errorf("want ErrNoBracket, got %v", err)
	}
}

func TestSignChangeAndBrackets(t *testing.T) {
	// sin over [0, 10] has roots at 0, pi, 2pi, 3pi
	f := math.Sin
	br := FindBrackets(f, 0.5, 10, 100)
	if len(br) != 3 {
		t.Errorf("FindBrackets found %d want 3", len(br))
	}
	if !SignChange(f, 3, 3.5) {
		t.Errorf("expected sign change around pi")
	}
	lo, hi, err := BracketOutward(func(x float64) float64 { return x - 5 }, 0, 1, 1.6, 100)
	if err != nil || SignChange(func(x float64) float64 { return x - 5 }, lo, hi) == false {
		t.Errorf("BracketOutward failed lo=%v hi=%v err=%v", lo, hi, err)
	}
}

// --- polynomial global solvers ---------------------------------------------

func TestGlobalSolversRealRoots(t *testing.T) {
	p := FromRoots(-2, 1, 3, 5) // quartic
	want := []complex128{
		complex(-2, 0), complex(1, 0), complex(3, 0), complex(5, 0),
	}
	dk, _, err := DurandKerner(p.ToComplex(), 1e-14, 300)
	if err != nil || !matchComplex(dk, want, 1e-7) {
		t.Errorf("DurandKerner=%v err=%v", dk, err)
	}
	ab, _, err := AberthEhrlich(p.ToComplex(), 1e-14, 300)
	if err != nil || !matchComplex(ab, want, 1e-7) {
		t.Errorf("AberthEhrlich=%v err=%v", ab, err)
	}
	comp, err := CompanionEigenvalues(p)
	if err != nil || !matchComplex(comp, want, 1e-6) {
		t.Errorf("Companion=%v err=%v", comp, err)
	}
	bs, err := BairstowRoots(p, 1e-13, 400)
	if err != nil || !matchComplex(bs, want, 1e-6) {
		t.Errorf("Bairstow=%v err=%v", bs, err)
	}
	lg, err := LaguerreRoots(p.ToComplex(), 1e-14, 300)
	if err != nil || !matchComplex(lg, want, 1e-7) {
		t.Errorf("Laguerre=%v err=%v", lg, err)
	}
}

func TestGlobalSolversComplexRoots(t *testing.T) {
	// (x^2+1)(x^2+4): roots +-i, +-2i
	p := NewPoly(4, 0, 5, 0, 1)
	want := []complex128{
		complex(0, 1), complex(0, -1), complex(0, 2), complex(0, -2),
	}
	comp, err := CompanionEigenvalues(p)
	if err != nil || !matchComplex(comp, want, 1e-7) {
		t.Errorf("companion=%v err=%v", comp, err)
	}
	bs, err := BairstowRoots(p, 1e-13, 400)
	if err != nil || !matchComplex(bs, want, 1e-7) {
		t.Errorf("bairstow=%v err=%v", bs, err)
	}
	ab, _, err := AberthEhrlich(p.ToComplex(), 1e-14, 400)
	if err != nil || !matchComplex(ab, want, 1e-7) {
		t.Errorf("aberth=%v err=%v", ab, err)
	}
}

func TestComplexCoeffPolynomial(t *testing.T) {
	// (x - (1+2i))(x - (3-i)) with complex coefficients
	c := CFromRoots(complex(1, 2), complex(3, -1))
	want := []complex128{complex(1, 2), complex(3, -1)}
	roots, err := CPolyRoots(c)
	if err != nil || !matchComplex(roots, want, 1e-9) {
		t.Errorf("CPolyRoots=%v err=%v", roots, err)
	}
}

func TestQuadraticRoots(t *testing.T) {
	tests := []struct {
		a, b, c float64
		want    []complex128
	}{
		{1, -3, 2, []complex128{complex(1, 0), complex(2, 0)}},
		{1, 0, 1, []complex128{complex(0, 1), complex(0, -1)}},
		{1, -2, 1, []complex128{complex(1, 0), complex(1, 0)}},
	}
	for _, tc := range tests {
		r1, r2 := QuadraticRoots(tc.a, tc.b, tc.c)
		if !matchComplex([]complex128{r1, r2}, tc.want, 1e-12) {
			t.Errorf("QuadraticRoots(%v,%v,%v)=%v,%v want %v", tc.a, tc.b, tc.c, r1, r2, tc.want)
		}
	}
}

func TestSolveCubic(t *testing.T) {
	// x^3 - 6x^2 + 11x - 6 -> 1,2,3
	roots, err := SolveCubicReal(1, -6, 11, -6)
	if err != nil {
		t.Fatal(err)
	}
	if !matchReal(roots, []float64{1, 2, 3}, 1e-9) {
		t.Errorf("SolveCubicReal=%v", roots)
	}
	// x^3 - 1 -> one real root 1, two complex
	all, _ := SolveCubic(1, 0, 0, -1)
	want := []complex128{
		complex(1, 0),
		complex(-0.5, math.Sqrt(3)/2),
		complex(-0.5, -math.Sqrt(3)/2),
	}
	if !matchComplex(all, want, 1e-9) {
		t.Errorf("SolveCubic(x^3-1)=%v", all)
	}
	if DiscriminantCubic(1, -6, 11, -6) <= 0 {
		t.Errorf("discriminant should be positive for 3 distinct real roots")
	}
}

// --- Sturm / isolation / counting ------------------------------------------

func TestSturmCounting(t *testing.T) {
	tests := []struct {
		p     Poly
		total int
	}{
		{FromRoots(1, 2, 3), 3},
		{NewPoly(1, 0, 1), 0},       // x^2+1
		{NewPoly(4, 0, 5, 0, 1), 0}, // (x^2+1)(x^2+4)
		{FromRoots(-2, 1, 3, 5), 4}, // quartic
		{FromRoots(1, 1, 2), 2},     // repeated root -> 2 distinct
		{NewPoly(-1, 0, 1), 2},      // x^2-1
	}
	for _, tc := range tests {
		if n := CountRealRoots(tc.p); n != tc.total {
			t.Errorf("CountRealRoots(%v)=%d want %d", tc.p, n, tc.total)
		}
	}
}

func TestCountRealRootsInInterval(t *testing.T) {
	p := FromRoots(-2, 1, 3, 5)
	tests := []struct {
		a, b float64
		n    int
	}{
		{0, 4, 2},  // 1 and 3
		{-3, 0, 1}, // -2
		{4, 6, 1},  // 5
		{-3, 6, 4}, // all
		{1.5, 2.5, 0},
	}
	for _, tc := range tests {
		if n := CountRealRootsInInterval(p, tc.a, tc.b); n != tc.n {
			t.Errorf("count in (%v,%v]=%d want %d", tc.a, tc.b, n, tc.n)
		}
	}
}

func TestIsolateAndRefine(t *testing.T) {
	p := FromRoots(-2, 1, 3, 5)
	iv := IsolateRoots(p)
	if len(iv) != 4 {
		t.Fatalf("IsolateRoots got %d intervals want 4", len(iv))
	}
	roots := SturmRealRoots(p, 1e-12)
	if !matchReal(roots, []float64{-2, 1, 3, 5}, 1e-9) {
		t.Errorf("SturmRealRoots=%v", roots)
	}
}

func TestSturmEvenMultiplicity(t *testing.T) {
	// (x-2)^2 (x+1): double root at 2 (even multiplicity, no sign change)
	p := FromRoots(2, 2, -1)
	roots := SturmRealRoots(p, 1e-10)
	if !matchReal(roots, []float64{-1, 2}, 1e-5) {
		t.Errorf("SturmRealRoots(even mult)=%v want {-1,2}", roots)
	}
	if CountRealRoots(p) != 2 {
		t.Errorf("distinct count=%d want 2", CountRealRoots(p))
	}
}

// --- Descartes / Budan-Fourier ---------------------------------------------

func TestDescartes(t *testing.T) {
	tests := []struct {
		p        Poly
		pos, neg int
	}{
		{FromRoots(1, 2, 3), 3, 0},
		{FromRoots(-1, -2, -3), 0, 3},
		{FromRoots(-2, 1, 3, 5), 3, 1},
		{NewPoly(1, 0, 1), 0, 0}, // x^2+1
	}
	for _, tc := range tests {
		pos, neg := DescartesRuleOfSigns(tc.p)
		if pos != tc.pos || neg != tc.neg {
			t.Errorf("Descartes(%v)=(%d,%d) want (%d,%d)", tc.p, pos, neg, tc.pos, tc.neg)
		}
	}
}

func TestSignVariations(t *testing.T) {
	if v := SignVariations([]float64{1, -1, 1, -1}); v != 3 {
		t.Errorf("SignVariations=%d want 3", v)
	}
	if v := SignVariations([]float64{1, 0, 0, 1}); v != 0 {
		t.Errorf("SignVariations with zeros=%d want 0", v)
	}
	if v := SignVariations([]float64{-2, 0, 3, -1}); v != 2 {
		t.Errorf("SignVariations=%d want 2", v)
	}
}

func TestBudanFourier(t *testing.T) {
	// three positive roots 1,2,3 -> exactly 3 roots in (0,4]
	p := FromRoots(1, 2, 3)
	if c := BudanFourierCount(p, 0, 4); c != 3 {
		t.Errorf("BudanFourierCount=%d want 3", c)
	}
	if c := BudanFourierCount(p, 0, 1.5); c != 1 {
		t.Errorf("BudanFourierCount(0,1.5)=%d want 1", c)
	}
}

// --- bounds ----------------------------------------------------------------

func TestRootBounds(t *testing.T) {
	p := FromRoots(-2, 1, 3, 5) // max |root| = 5
	bounds := map[string]float64{
		"Cauchy":   CauchyBound(p),
		"Lagrange": LagrangeBound(p),
		"Fujiwara": FujiwaraBound(p),
		"Kojima":   KojimaBound(p),
	}
	for name, b := range bounds {
		if b < 5 {
			t.Errorf("%s bound=%v must be >= 5 (max root modulus)", name, b)
		}
	}
	// all roots must lie within Cauchy bound in modulus
	roots, _ := PolyRoots(p)
	cb := CauchyBound(p)
	for _, r := range roots {
		if cmplx.Abs(r) > cb+1e-9 {
			t.Errorf("root %v exceeds Cauchy bound %v", r, cb)
		}
	}
}

func TestLowerBound(t *testing.T) {
	// roots 2,4,8: smallest modulus 2
	p := FromRoots(2, 4, 8)
	lo := LowerRootBound(p)
	if lo <= 0 || lo > 2 {
		t.Errorf("LowerRootBound=%v want in (0,2]", lo)
	}
}

// --- multiplicity / squarefree ---------------------------------------------

func TestMultiplicity(t *testing.T) {
	// (x-1)^3 (x-2)
	p := FromRoots(1, 1, 1, 2)
	if m := Multiplicity(p, 1, 1e-6); m != 3 {
		t.Errorf("Multiplicity(1)=%d want 3", m)
	}
	if m := Multiplicity(p, 2, 1e-6); m != 1 {
		t.Errorf("Multiplicity(2)=%d want 1", m)
	}
	if m := Multiplicity(p, 5, 1e-6); m != 0 {
		t.Errorf("Multiplicity(5)=%d want 0", m)
	}
}

func TestSquareFree(t *testing.T) {
	// (x-1)^2 (x-3)^3 -> squarefree part (x-1)(x-3), roots 1 and 3 simple
	p := FromRoots(1, 1, 3, 3, 3)
	sf, err := SquareFree(p)
	if err != nil {
		t.Fatal(err)
	}
	if sf.Degree() != 2 {
		t.Errorf("squarefree degree=%d want 2 (%v)", sf.Degree(), sf)
	}
	if !approx(sf.Eval(1), 0, 1e-6) || !approx(sf.Eval(3), 0, 1e-6) {
		t.Errorf("squarefree should vanish at 1 and 3: %v", sf)
	}
}

func TestSquareFreeFactorization(t *testing.T) {
	// (x-1)^1 (x-2)^2 (x-3)^3
	p := FromRoots(1).Mul(FromRoots(2, 2)).Mul(FromRoots(3, 3, 3))
	facs, err := SquareFreeFactorization(p)
	if err != nil {
		t.Fatal(err)
	}
	got := map[int]float64{} // multiplicity -> root
	for _, f := range facs {
		if f.Factor.Degree() != 1 {
			t.Errorf("factor of unexpected degree: %v", f.Factor)
			continue
		}
		got[f.Multiplicity] = -f.Factor.Coeff(0) / f.Factor.Coeff(1)
	}
	want := map[int]float64{1: 1, 2: 2, 3: 3}
	for mult, root := range want {
		if !approx(got[mult], root, 1e-6) {
			t.Errorf("multiplicity %d root=%v want %v", mult, got[mult], root)
		}
	}
}

func TestRealRootsWithMultiplicity(t *testing.T) {
	p := FromRoots(-1, 2, 2)
	rm := RealRootsWithMultiplicity(p, 1e-9)
	if len(rm) != 2 {
		t.Fatalf("got %d distinct roots want 2: %v", len(rm), rm)
	}
	total := TotalRealRoots(p, 1e-9)
	if total != 3 {
		t.Errorf("TotalRealRoots=%d want 3", total)
	}
}

func TestGroupComplexRoots(t *testing.T) {
	roots := []complex128{
		complex(1, 0), complex(1+1e-9, 0), complex(1-1e-9, 0),
		complex(5, 0),
	}
	groups := GroupComplexRoots(roots, 1e-6)
	if len(groups) != 2 {
		t.Fatalf("got %d groups want 2", len(groups))
	}
	// find the triple cluster
	for _, g := range groups {
		if approxC(g.Root, complex(1, 0), 1e-6) && g.Multiplicity != 3 {
			t.Errorf("cluster at 1 multiplicity=%d want 3", g.Multiplicity)
		}
	}
}

// --- deflation -------------------------------------------------------------

func TestDeflation(t *testing.T) {
	p := FromRoots(1, 2, 3)
	q, rem := p.DeflateReal(2)
	if math.Abs(rem) > 1e-9 {
		t.Errorf("DeflateReal remainder=%v want 0", rem)
	}
	// q should have roots 1 and 3
	if !approx(q.Eval(1), 0, 1e-9) || !approx(q.Eval(3), 0, 1e-9) {
		t.Errorf("deflated poly wrong: %v", q)
	}
	// complex deflation
	c := CFromRoots(complex(0, 1), complex(0, -1))
	cq, crem := c.Deflate(complex(0, 1))
	if cmplx.Abs(crem) > 1e-9 {
		t.Errorf("CPoly.Deflate remainder=%v", crem)
	}
	if cmplx.Abs(cq.Eval(complex(0, -1))) > 1e-9 {
		t.Errorf("complex deflation wrong")
	}
	// deflate list
	dr := DeflateRoots(p, []float64{1, 2})
	if dr.Degree() != 1 || !approx(dr.Eval(3), 0, 1e-9) {
		t.Errorf("DeflateRoots=%v", dr)
	}
}

// --- companion matrix structure --------------------------------------------

func TestCompanionMatrix(t *testing.T) {
	p := NewPoly(-6, 11, -6, 1) // monic x^3-6x^2+11x-6
	m, err := CompanionMatrix(p)
	if err != nil {
		t.Fatal(err)
	}
	if len(m) != 3 {
		t.Fatalf("companion size=%d want 3", len(m))
	}
	// first row = -a2,-a1,-a0 = 6,-11,6
	wantRow := []float64{6, -11, 6}
	for j := 0; j < 3; j++ {
		if !approx(m[0][j], wantRow[j], 1e-12) {
			t.Errorf("companion[0][%d]=%v want %v", j, m[0][j], wantRow[j])
		}
	}
	if m[1][0] != 1 || m[2][1] != 1 {
		t.Errorf("subdiagonal ones missing")
	}
}

// --- complex Newton / Halley -----------------------------------------------

func TestComplexNewtonHalley(t *testing.T) {
	c := CFromRoots(complex(2, 3), complex(-1, 1))
	r, _, err := NewtonComplex(c, complex(2.1, 3.1), 1e-14, 100)
	if err != nil || !approxC(r, complex(2, 3), 1e-9) {
		t.Errorf("NewtonComplex=%v err=%v", r, err)
	}
	r2, _, err := HalleyComplex(c, complex(-1.1, 0.9), 1e-14, 100)
	if err != nil || !approxC(r2, complex(-1, 1), 1e-9) {
		t.Errorf("HalleyComplex=%v err=%v", r2, err)
	}
}

func TestHornerFree(t *testing.T) {
	if !approx(Horner([]float64{1, 2, 3}, 2), 17, 1e-12) {
		t.Errorf("Horner=%v", Horner([]float64{1, 2, 3}, 2))
	}
	if !approxC(HornerComplex([]complex128{1, 0, 1}, complex(0, 1)), 0, 1e-12) {
		t.Errorf("HornerComplex(i) of x^2+1 should be 0")
	}
}

// --- residual / separation -------------------------------------------------

func TestSeparateAndResidual(t *testing.T) {
	p := FromRoots(-2, 1, 3, 5).Mul(NewPoly(1, 0, 1)) // add +-i
	roots, err := PolyRoots(p)
	if err != nil {
		t.Fatal(err)
	}
	reals, comps := SeparateRoots(roots, 1e-6)
	if !matchReal(reals, []float64{-2, 1, 3, 5}, 1e-6) {
		t.Errorf("real roots=%v", reals)
	}
	if len(comps) != 2 {
		t.Errorf("complex roots=%v want 2", comps)
	}
	if res := MaxResidual(p.ToComplex(), roots); res > 1e-6 {
		t.Errorf("MaxResidual=%v too large", res)
	}
}

// --- example ---------------------------------------------------------------

func ExamplePolyRoots() {
	// x^3 - 6x^2 + 11x - 6 = (x-1)(x-2)(x-3)
	p := NewPoly(-6, 11, -6, 1)
	roots, _ := PolyRoots(p)
	for _, r := range roots {
		fmt.Printf("%.4f\n", real(r))
	}
	// Output:
	// 1.0000
	// 2.0000
	// 3.0000
}

func ExampleBrent() {
	// Solve cos(x) = x on [0, 1].
	f := func(x float64) float64 { return math.Cos(x) - x }
	res, _ := Brent(f, 0, 1, 1e-12, 100)
	fmt.Printf("%.10f\n", res.Root)
	// Output:
	// 0.7390851332
}

func ExampleSturmRealRoots() {
	// (x+2)(x-1)(x-3): three distinct real roots.
	p := FromRoots(-2, 1, 3)
	for _, r := range SturmRealRoots(p, 1e-12) {
		fmt.Printf("%.4f\n", r)
	}
	// Output:
	// -2.0000
	// 1.0000
	// 3.0000
}

func ExampleDescartesRuleOfSigns() {
	// x^2 - 1 has one sign change -> one positive, one negative root.
	p := NewPoly(-1, 0, 1)
	pos, neg := DescartesRuleOfSigns(p)
	fmt.Println(pos, neg)
	// Output:
	// 1 1
}
