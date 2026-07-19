package exterioralgebra

import (
	"fmt"
	"math"
	"testing"
)

const tol = 1e-9

// ---------------------------------------------------------------------------
// Blade / bit helpers
// ---------------------------------------------------------------------------

func TestPopcountAndMasks(t *testing.T) {
	cases := []struct {
		mask uint
		pop  int
		idx  []int
	}{
		{0, 0, nil},
		{0b1, 1, []int{0}},
		{0b1011, 3, []int{0, 1, 3}},
		{0b10000, 1, []int{4}},
	}
	for _, c := range cases {
		if got := Popcount(c.mask); got != c.pop {
			t.Errorf("Popcount(%b)=%d want %d", c.mask, got, c.pop)
		}
		got := MaskToIndices(c.mask)
		if len(got) != len(c.idx) {
			t.Errorf("MaskToIndices(%b)=%v want %v", c.mask, got, c.idx)
			continue
		}
		for i := range got {
			if got[i] != c.idx[i] {
				t.Errorf("MaskToIndices(%b)=%v want %v", c.mask, got, c.idx)
			}
		}
	}
	if FullMask(4) != 0b1111 {
		t.Errorf("FullMask(4)=%b", FullMask(4))
	}
}

func TestIndicesToMask(t *testing.T) {
	cases := []struct {
		idx  []int
		mask uint
		sign int
		ok   bool
	}{
		{[]int{0, 1}, 0b11, 1, true},
		{[]int{1, 0}, 0b11, -1, true},
		{[]int{2, 0, 1}, 0b111, 1, true},  // even permutation
		{[]int{0, 2, 1}, 0b111, -1, true}, // odd permutation
		{[]int{1, 1}, 0, 0, false},        // repeated
		{[]int{0, 9}, 0, 0, false},        // out of range (n=4)
	}
	for _, c := range cases {
		m, s, ok := IndicesToMask(4, c.idx...)
		if ok != c.ok || (ok && (m != c.mask || s != c.sign)) {
			t.Errorf("IndicesToMask(%v)=(%b,%d,%v) want (%b,%d,%v)",
				c.idx, m, s, ok, c.mask, c.sign, c.ok)
		}
	}
}

// ---------------------------------------------------------------------------
// Construction & inspection
// ---------------------------------------------------------------------------

func TestConstructionAndCoeff(t *testing.T) {
	f := FromVector([]float64{1, 2, 3})
	if f.Dim() != 3 {
		t.Fatalf("Dim=%d", f.Dim())
	}
	if got := f.Coeff(1); got != 2 {
		t.Errorf("Coeff(1)=%v want 2", got)
	}
	// grade-2 blade with swapped indices flips sign
	b, err := BasisBlade(4, 2, 1)
	if err != nil {
		t.Fatal(err)
	}
	if got := b.Coeff(1, 2); got != -1 {
		t.Errorf("Coeff(1,2)=%v want -1", got)
	}
	if got := b.Coeff(2, 1); got != 1 {
		t.Errorf("Coeff(2,1)=%v want 1", got)
	}
	// repeated index blade is zero
	z, _ := BasisBlade(4, 1, 1)
	if !z.IsZero() {
		t.Errorf("e1∧e1 should be zero, got %v", z)
	}
}

func TestGradeInspection(t *testing.T) {
	e01, _ := BasisBlade(3, 0, 1)
	mixed := Scalar(3, 5).Add(e01)
	if mixed.IsHomogeneous() {
		t.Error("scalar+bivector should be inhomogeneous")
	}
	if _, ok := mixed.Grade(); ok {
		t.Error("Grade should report not-ok for mixed form")
	}
	if got := mixed.GradeProject(0).ScalarPart(); got != 5 {
		t.Errorf("scalar part=%v want 5", got)
	}
	if k, ok := e01.Grade(); !ok || k != 2 {
		t.Errorf("Grade(e01)=(%d,%v)", k, ok)
	}
	gs := mixed.Grades()
	if len(gs) != 2 || gs[0] != 0 || gs[1] != 2 {
		t.Errorf("Grades=%v want [0 2]", gs)
	}
	if mixed.MaxGrade() != 2 || mixed.MinGrade() != 0 {
		t.Errorf("Max/Min grade = %d/%d", mixed.MaxGrade(), mixed.MinGrade())
	}
}

// ---------------------------------------------------------------------------
// Linear structure
// ---------------------------------------------------------------------------

func TestLinearOps(t *testing.T) {
	a := FromVector([]float64{1, 2, 3})
	b := FromVector([]float64{4, 5, 6})
	sum := a.Add(b)
	if !sum.Equal(FromVector([]float64{5, 7, 9})) {
		t.Errorf("Add=%v", sum)
	}
	diff := b.Sub(a)
	if !diff.Equal(FromVector([]float64{3, 3, 3})) {
		t.Errorf("Sub=%v", diff)
	}
	if !a.Scale(2).Equal(FromVector([]float64{2, 4, 6})) {
		t.Errorf("Scale=%v", a.Scale(2))
	}
	if !a.Neg().Equal(a.Scale(-1)) {
		t.Error("Neg != Scale(-1)")
	}
	lc := LinearCombination([]float64{2, -1}, []*Form{a, b})
	if !lc.Equal(FromVector([]float64{-2, -1, 0})) {
		t.Errorf("LinearCombination=%v", lc)
	}
}

// ---------------------------------------------------------------------------
// Wedge product
// ---------------------------------------------------------------------------

func TestWedge(t *testing.T) {
	e0 := Basis1(3, 0)
	e1 := Basis1(3, 1)
	e2 := Basis1(3, 2)

	// antisymmetry: e0∧e1 = -(e1∧e0)
	if !Wedge(e0, e1).Equal(Wedge(e1, e0).Neg()) {
		t.Error("wedge not antisymmetric")
	}
	// self wedge vanishes
	if !Wedge(e0, e0).IsZero() {
		t.Error("e0∧e0 should vanish")
	}
	// triple product is the volume form
	vol := WedgeAll(e0, e1, e2)
	if !vol.Equal(VolumeForm(3)) {
		t.Errorf("e0∧e1∧e2=%v want volume form", vol)
	}
	// graded commutativity for a vector (p=1) and bivector (q=2): sign (-1)^2=+1
	e12 := Wedge(e1, e2)
	if !Wedge(e0, e12).Equal(Wedge(e12, e0)) {
		t.Error("vector∧bivector should commute")
	}
	// distributivity
	x := FromVector([]float64{1, 1, 0})
	y := FromVector([]float64{0, 1, 1})
	lhs := Wedge(e0, x.Add(y))
	rhs := Wedge(e0, x).Add(Wedge(e0, y))
	if !lhs.EqualTol(rhs, tol) {
		t.Error("wedge not distributive")
	}
}

func TestWedgeDeterminant(t *testing.T) {
	// (a e0 + b e1) ∧ (c e0 + d e1) = (ad - bc) e0∧e1
	a, b, c, d := 2.0, 3.0, 5.0, 7.0
	u := FromVector([]float64{a, b})
	v := FromVector([]float64{c, d})
	w := u.Wedge(v)
	if got := w.Coeff(0, 1); math.Abs(got-(a*d-b*c)) > tol {
		t.Errorf("2x2 wedge det=%v want %v", got, a*d-b*c)
	}
	det, err := Determinant([]float64{1, 2, 3}, []float64{0, 1, 4}, []float64{5, 6, 0})
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(det-1) > tol { // det of that matrix is 1
		t.Errorf("Determinant=%v want 1", det)
	}
}

func TestReverseAndInvolution(t *testing.T) {
	// reversion sign on grades 0..3 is +,+,-,-
	for k, want := range map[int]float64{0: 1, 1: 1, 2: -1, 3: -1} {
		idx := make([]int, k)
		for i := range idx {
			idx[i] = i
		}
		blade, _ := BasisBlade(3, idx...)
		if blade.IsZero() && k <= 3 {
			blade = VolumeForm(3) // k==3 top form
		}
		rev := blade.Reverse()
		if k <= 3 && k != 3 {
			if math.Abs(rev.CoeffMask(blade.Masks()[0])-want) > tol {
				t.Errorf("Reverse grade %d sign=%v want %v", k, rev.CoeffMask(blade.Masks()[0]), want)
			}
		}
	}
	e01, _ := BasisBlade(3, 0, 1)
	if !e01.GradeInvolution().Equal(e01) {
		t.Error("grade involution should fix even blade")
	}
	e0 := Basis1(3, 0)
	if !e0.GradeInvolution().Equal(e0.Neg()) {
		t.Error("grade involution should negate odd blade")
	}
}

// ---------------------------------------------------------------------------
// Interior product & contraction
// ---------------------------------------------------------------------------

func TestInteriorProduct(t *testing.T) {
	e0 := Basis1(3, 0)
	e012 := VolumeForm(3)
	// ι_{e0}(e0∧e1∧e2) = e1∧e2
	got := InteriorProduct(e0, e012)
	want, _ := BasisBlade(3, 1, 2)
	if !got.Equal(want) {
		t.Errorf("ι_e0 vol=%v want %v", got, want)
	}
	// ι_{e1}(e0∧e1) = -e0
	e01, _ := BasisBlade(3, 0, 1)
	got = InteriorProduct(Basis1(3, 1), e01)
	if !got.Equal(Basis1(3, 0).Neg()) {
		t.Errorf("ι_e1(e0∧e1)=%v want -e0", got)
	}
	// antiderivation property: ι_v(α∧β) = (ι_v α)∧β + (-1)^|α| α∧(ι_v β)
	v := FromVector([]float64{1, 2, 3})
	alpha := Basis1(3, 0)                         // grade 1
	beta, _ := BasisBlade(3, 1, 2)                // grade 2
	lhs := InteriorProduct(v, alpha.Wedge(beta))  //
	rhs := InteriorProduct(v, alpha).Wedge(beta). //
							Sub(alpha.Wedge(InteriorProduct(v, beta)))
	if !lhs.EqualTol(rhs, tol) {
		t.Errorf("antiderivation failed: %v vs %v", lhs, rhs)
	}
}

func TestContractionsAndInner(t *testing.T) {
	// left contraction adjointness: <a⌋b, c> = <b, a∧c>
	a := FromVector([]float64{1, 0, 2})
	b, _ := BasisBlade(3, 0, 1)
	b = b.Add(VolumeForm(3)) // e0∧e1 + e0∧e1∧e2
	c := Basis1(3, 2)
	lhs := InnerProduct(a.LeftContract(b), c)
	rhs := InnerProduct(b, a.Wedge(c))
	if math.Abs(lhs-rhs) > tol {
		t.Errorf("left-contract adjoint: %v vs %v", lhs, rhs)
	}
	// Euclidean inner product & norm on orthonormal blades
	if got := InnerProduct(VolumeForm(3), VolumeForm(3)); got != 1 {
		t.Errorf("<vol,vol>=%v", got)
	}
	u := FromVector([]float64{3, 4, 0})
	if math.Abs(u.Norm()-5) > tol {
		t.Errorf("Norm=%v want 5", u.Norm())
	}
	// angle between orthogonal vectors is pi/2
	if got := Angle(Basis1(3, 0), Basis1(3, 1)); math.Abs(got-math.Pi/2) > tol {
		t.Errorf("Angle=%v want pi/2", got)
	}
}

// ---------------------------------------------------------------------------
// Hodge star
// ---------------------------------------------------------------------------

func TestHodgeStarBasic(t *testing.T) {
	// R^3: *e0 = e1∧e2, *e1 = e2∧e0 = -e0∧e2, *e2 = e0∧e1
	e12, _ := BasisBlade(3, 1, 2)
	if !HodgeStar(Basis1(3, 0)).Equal(e12) {
		t.Errorf("*e0=%v", HodgeStar(Basis1(3, 0)))
	}
	e02, _ := BasisBlade(3, 0, 2)
	if !HodgeStar(Basis1(3, 1)).Equal(e02.Neg()) {
		t.Errorf("*e1=%v want -e0∧e2", HodgeStar(Basis1(3, 1)))
	}
	// *1 = vol and *vol = 1
	if !HodgeStar(One(3)).Equal(VolumeForm(3)) {
		t.Error("*1 should be volume form")
	}
	if !HodgeStar(VolumeForm(3)).Equal(One(3)) {
		t.Error("*vol should be 1")
	}
}

func TestHodgeStarSquared(t *testing.T) {
	// ** = (-1)^{k(n-k)} id, tested across dimensions and grades.
	for n := 1; n <= 5; n++ {
		full := FullMask(n)
		for mask := uint(0); mask <= full; mask++ {
			if mask & ^full != 0 {
				continue
			}
			k := Popcount(mask)
			f := New(n)
			f.SetMask(mask, 1)
			ss := HodgeStar(HodgeStar(f))
			want := f.Scale(float64(HodgeSignSquared(k, n)))
			if !ss.EqualTol(want, tol) {
				t.Errorf("n=%d mask=%b: **=%v want %v", n, mask, ss, want)
			}
			// inverse star round-trips
			if !InverseHodgeStar(HodgeStar(f)).EqualTol(f, tol) {
				t.Errorf("n=%d mask=%b inverse star failed", n, mask)
			}
		}
	}
}

func TestHodgeStarMetric(t *testing.T) {
	// Euclidean signature agrees with HodgeStar.
	f := FromVector([]float64{1, 2, 3})
	got, err := HodgeStarMetric(f, EuclideanSignature(3))
	if err != nil {
		t.Fatal(err)
	}
	if !got.Equal(HodgeStar(f)) {
		t.Error("metric star with +1 signature should equal HodgeStar")
	}
	// Minkowski (n=4): *dt where signature = (-1,+1,+1,+1).
	sig := LorentzSignature(4)
	star, err := HodgeStarMetric(Basis1(4, 0), sig)
	if err != nil {
		t.Fatal(err)
	}
	// *e0 = reorderSign(e0, e1e2e3)*g00 * e1∧e2∧e3 = 1 * (-1) * e123
	want, _ := BasisBlade(4, 1, 2, 3)
	if !star.Equal(want.Neg()) {
		t.Errorf("Minkowski *e0=%v want -e1∧e2∧e3", star)
	}
	if _, err := HodgeStarMetric(f, []int{1, 2, 1}); err != ErrGrade {
		t.Errorf("bad signature should error, got %v", err)
	}
}

func TestCrossProduct(t *testing.T) {
	got, err := CrossProduct([]float64{1, 0, 0}, []float64{0, 1, 0})
	if err != nil {
		t.Fatal(err)
	}
	want := []float64{0, 0, 1}
	for i := range want {
		if math.Abs(got[i]-want[i]) > tol {
			t.Errorf("cross=%v want %v", got, want)
		}
	}
}

// ---------------------------------------------------------------------------
// Polynomials
// ---------------------------------------------------------------------------

func TestPolyArithmetic(t *testing.T) {
	x := Var(2, 0)
	y := Var(2, 1)
	// (x+y)^2 = x^2 + 2xy + y^2
	p := x.Add(y).Pow(2)
	want := x.Pow(2).Add(x.Mul(y).Scale(2)).Add(y.Pow(2))
	if !p.Equal(want) {
		t.Errorf("(x+y)^2=%s want %s", p, want)
	}
	if got := p.Eval([]float64{2, 3}); math.Abs(got-25) > tol {
		t.Errorf("eval (2+3)^2=%v want 25", got)
	}
	// partial derivatives
	if !p.Partial(0).Equal(x.Scale(2).Add(y.Scale(2))) {
		t.Errorf("d/dx (x+y)^2 = %s", p.Partial(0))
	}
	if p.Degree() != 2 {
		t.Errorf("degree=%d", p.Degree())
	}
}

func TestPolyCompose(t *testing.T) {
	// p(x0,x1) = x0^2 + x1 ; substitute x0=t, x1=t^2 -> t^2 + t^2 = 2t^2? no: t^2 + t^2
	n := 2
	p := Var(n, 0).Pow(2).Add(Var(n, 1))
	t2 := Var(1, 0).Pow(2)
	subs := []*Poly{Var(1, 0), t2}
	c, err := p.Compose(subs)
	if err != nil {
		t.Fatal(err)
	}
	// x0^2 -> t^2, x1 -> t^2, sum = 2 t^2
	want := Var(1, 0).Pow(2).Scale(2)
	if !c.Equal(want) {
		t.Errorf("compose=%s want %s", c, want)
	}
	// consistency with evaluation: c(t) == p(t, t^2)
	for _, tv := range []float64{-1, 0, 2.5} {
		got := c.Eval([]float64{tv})
		exp := p.Eval([]float64{tv, tv * tv})
		if math.Abs(got-exp) > tol {
			t.Errorf("compose eval at %v: %v want %v", tv, got, exp)
		}
	}
}

// ---------------------------------------------------------------------------
// Differential forms: d, d^2, pullback, Laplacian
// ---------------------------------------------------------------------------

func TestExteriorDerivativeAndClosure(t *testing.T) {
	n := 3
	// f = x0^2 x1 + x2  -> df = 2 x0 x1 dx0 + x0^2 dx1 + dx2
	f := Var(n, 0).Pow(2).Mul(Var(n, 1)).Add(Var(n, 2))
	df := PConst(f).ExteriorDerivative()
	want := NewPForm(n)
	c0 := Var(n, 0).Mul(Var(n, 1)).Scale(2)
	want = want.Add(mustBlade(t, n, c0, 0))
	want = want.Add(mustBlade(t, n, Var(n, 0).Pow(2), 1))
	want = want.Add(mustBlade(t, n, ConstPoly(n, 1), 2))
	if !df.Equal(want) {
		t.Errorf("df=%s want %s", df, want)
	}
	// d^2 = 0 on a range of random-ish polynomial forms
	forms := []*PForm{
		PConst(f),
		mustBlade(t, n, Var(n, 0).Mul(Var(n, 2)), 1),
		mustBlade(t, n, Var(n, 1).Pow(3), 0).Add(mustBlade(t, n, Var(n, 0), 2)),
	}
	for i, w := range forms {
		if !w.ExteriorDerivative().ExteriorDerivative().IsZero() {
			t.Errorf("d^2 != 0 on form %d", i)
		}
	}
}

func TestPullbackCommutesWithD(t *testing.T) {
	// omega on R^2, map phi: R^2 -> R^2, phi(y0,y1) = (y0^2 - y1, y0 y1)
	n := 2
	omega := mustBlade(t, n, Var(n, 0).Mul(Var(n, 1)), 0) // (x0 x1) dx0
	phi := []*Poly{
		Var(2, 0).Pow(2).Sub(Var(2, 1)),
		Var(2, 0).Mul(Var(2, 1)),
	}
	// d(phi* omega) should equal phi*(d omega)
	pbThenD, err := omega.Pullback(phi)
	if err != nil {
		t.Fatal(err)
	}
	lhs := pbThenD.ExteriorDerivative()
	rhs, err := omega.ExteriorDerivative().Pullback(phi)
	if err != nil {
		t.Fatal(err)
	}
	if !lhs.EqualTol(rhs, tol) {
		t.Errorf("pullback/d mismatch:\n d(phi*w)=%s\n phi*(dw)=%s", lhs, rhs)
	}
}

func TestPullbackWedgeHomomorphism(t *testing.T) {
	// phi*(a ∧ b) = phi*a ∧ phi*b
	n := 3
	a := Dx(n, 0).ScalePoly(Var(n, 1))
	b := Dx(n, 2)
	phi := []*Poly{
		Var(2, 0),
		Var(2, 1).Pow(2),
		Var(2, 0).Mul(Var(2, 1)),
	}
	left, err := a.Wedge(b).Pullback(phi)
	if err != nil {
		t.Fatal(err)
	}
	pa, _ := a.Pullback(phi)
	pb, _ := b.Pullback(phi)
	right := pa.Wedge(pb)
	if !left.EqualTol(right, tol) {
		t.Errorf("pullback not a wedge homomorphism:\n %s\n %s", left, right)
	}
}

func TestHodgeLaplacianScalar(t *testing.T) {
	// For a 0-form f, Delta f = -sum d^2 f/dxi^2.
	n := 3
	f := Var(n, 0).Pow(3).Add(Var(n, 1).Pow(2).Mul(Var(n, 2)))
	lap := PConst(f).HodgeLaplacian()
	// analytic scalar Laplacian
	sl := NewPoly(n)
	for i := 0; i < n; i++ {
		sl = sl.Add(f.Partial(i).Partial(i))
	}
	want := PConst(sl.Neg())
	if !lap.EqualTol(want, tol) {
		t.Errorf("HodgeLaplacian=%s want %s", lap, want)
	}
}

func TestPFormEval(t *testing.T) {
	n := 2
	w := mustBlade(t, n, Var(n, 0).Mul(Var(n, 1)), 0).Add(mustBlade(t, n, Var(n, 0), 1))
	f := w.Eval([]float64{2, 3})
	if got := f.Coeff(0); math.Abs(got-6) > tol {
		t.Errorf("eval dx0 coeff=%v want 6", got)
	}
	if got := f.Coeff(1); math.Abs(got-2) > tol {
		t.Errorf("eval dx1 coeff=%v want 2", got)
	}
}

// ---------------------------------------------------------------------------
// Numerical exterior derivative
// ---------------------------------------------------------------------------

func TestNumGradient(t *testing.T) {
	f := func(x []float64) float64 { return x[0]*x[0] + x[1]*x[2] }
	g := NumGradient(f, []float64{1, 2, 3}, 1e-5)
	// grad = (2 x0, x2, x1) = (2, 3, 2)
	want := []float64{2, 3, 2}
	for i := range want {
		if math.Abs(g.Coeff(i)-want[i]) > 1e-4 {
			t.Errorf("NumGradient[%d]=%v want %v", i, g.Coeff(i), want[i])
		}
	}
}

func TestNumExteriorDerivativeMatchesPoly(t *testing.T) {
	// Build a 1-form field a0(x) dx0 + a1(x) dx1 on R^2 and compare numeric d
	// against the exact polynomial exterior derivative at a point.
	n := 2
	a0 := Var(n, 0).Mul(Var(n, 1)) // x0 x1
	a1 := Var(n, 0).Pow(2)         // x0^2
	pform := mustBlade(t, n, a0, 0).Add(mustBlade(t, n, a1, 1))
	field := func(x []float64) *Form {
		return pform.Eval(x)
	}
	x := []float64{1.3, -0.7}
	numeric := NumExteriorDerivative(field, x, 1e-5)
	exact := pform.ExteriorDerivative().Eval(x)
	if !numeric.EqualTol(exact, 1e-4) {
		t.Errorf("numeric d=%v exact d=%v", numeric, exact)
	}
}

// ---------------------------------------------------------------------------
// Extras: basis enumeration, grade helpers, vector calculus
// ---------------------------------------------------------------------------

func TestBinomialAndBasis(t *testing.T) {
	if Binomial(5, 2) != 10 || Binomial(4, 0) != 1 || Binomial(3, 4) != 0 {
		t.Errorf("Binomial wrong: %d %d %d", Binomial(5, 2), Binomial(4, 0), Binomial(3, 4))
	}
	if AlgebraDimension(4) != 16 {
		t.Errorf("AlgebraDimension(4)=%d", AlgebraDimension(4))
	}
	b := Basis(4, 2)
	if len(b) != GradeDimension(4, 2) || len(b) != 6 {
		t.Fatalf("Basis(4,2) has %d blades want 6", len(b))
	}
	// every returned blade is a unit grade-2 form
	for _, f := range b {
		if k, ok := f.Grade(); !ok || k != 2 || math.Abs(f.Norm()-1) > tol {
			t.Errorf("bad basis blade %v", f)
		}
	}
	if len(AllBasisBlades(3)) != 8 {
		t.Errorf("AllBasisBlades(3)=%d want 8", len(AllBasisBlades(3)))
	}
}

func TestEvenOddScalarVector(t *testing.T) {
	f := Scalar(3, 2).Add(Basis1(3, 0))
	e01, _ := BasisBlade(3, 0, 1)
	f = f.Add(e01)
	if !f.EvenPart().Equal(Scalar(3, 2).Add(e01)) {
		t.Errorf("EvenPart=%v", f.EvenPart())
	}
	if !f.OddPart().Equal(Basis1(3, 0)) {
		t.Errorf("OddPart=%v", f.OddPart())
	}
	if !Basis1(3, 1).IsVector() || Scalar(3, 4).IsVector() {
		t.Error("IsVector wrong")
	}
	if !Scalar(3, 4).IsScalar() || Basis1(3, 0).IsScalar() {
		t.Error("IsScalar wrong")
	}
}

func TestPolyVectorCalculusDeRham(t *testing.T) {
	// curl(grad f) = 0 for a scalar field on R^3.
	f := Var(3, 0).Pow(2).Mul(Var(3, 1)).Add(Var(3, 2).Pow(3)).Add(Var(3, 0).Mul(Var(3, 2)))
	grad := PolyGradient(f)
	curl, err := PolyCurl(grad)
	if err != nil {
		t.Fatal(err)
	}
	for i, c := range curl {
		if !c.IsZero() {
			t.Errorf("curl(grad f)[%d]=%s want 0", i, c)
		}
	}
	// div(curl F) = 0 for a vector field on R^3.
	field := []*Poly{
		Var(3, 1).Mul(Var(3, 2)),
		Var(3, 0).Pow(2),
		Var(3, 0).Mul(Var(3, 1)).Add(Var(3, 2).Pow(2)),
	}
	cf, err := PolyCurl(field)
	if err != nil {
		t.Fatal(err)
	}
	div, err := PolyDivergence(cf)
	if err != nil {
		t.Fatal(err)
	}
	if !div.IsZero() {
		t.Errorf("div(curl F)=%s want 0", div)
	}
	// Laplacian consistency: Laplacian(x0^2+x1^2+x2^2) = 6.
	sq := Var(3, 0).Pow(2).Add(Var(3, 1).Pow(2)).Add(Var(3, 2).Pow(2))
	if got := PolyLaplacian(sq); !got.Equal(ConstPoly(3, 6)) {
		t.Errorf("Laplacian=%s want 6", got)
	}
}

func TestPFormInteriorConstAndClosed(t *testing.T) {
	n := 3
	// dx0∧dx1 contracted with e0 gives dx1.
	w := Dx(n, 0).Wedge(Dx(n, 1))
	got := w.InteriorConst([]float64{1, 0, 0})
	if !got.Equal(Dx(n, 1)) {
		t.Errorf("interior=%s want dx1", got)
	}
	// A constant-coefficient top form is closed, and any exact form is closed.
	if !VolumeFormP(n).IsClosed() {
		t.Error("constant volume form should be closed")
	}
	f := Var(n, 0).Mul(Var(n, 1).Pow(2))
	if !PConst(f).ExteriorDerivative().IsClosed() {
		t.Error("df should be closed (d^2=0)")
	}
}

// VolumeFormP is a small test helper building the constant top-degree form.
func VolumeFormP(n int) *PForm {
	idx := make([]int, n)
	for i := range idx {
		idx[i] = i
	}
	b, _ := PBasisBlade(n, ConstPoly(n, 1), idx...)
	return b
}

func mustBlade(t *testing.T, n int, coeff *Poly, idx ...int) *PForm {
	t.Helper()
	w, err := PBasisBlade(n, coeff, idx...)
	if err != nil {
		t.Fatal(err)
	}
	return w
}

// ---------------------------------------------------------------------------
// Runnable examples
// ---------------------------------------------------------------------------

func ExampleHodgeStar() {
	// In R^3 the Hodge dual of the basis 1-form e0 is the area element e1∧e2.
	fmt.Println(HodgeStar(Basis1(3, 0)))
	// Output: 1 e1∧e2
}

func ExampleForm_Wedge() {
	u := FromVector([]float64{1, 0, 0})
	v := FromVector([]float64{0, 1, 0})
	fmt.Println(u.Wedge(v))
	// Output: 1 e0∧e1
}

func ExamplePForm_ExteriorDerivative() {
	// d(x0 x1) = x1 dx0 + x0 dx1 on R^2.
	n := 2
	f := Var(n, 0).Mul(Var(n, 1))
	fmt.Println(PConst(f).ExteriorDerivative())
	// Output: (x1) dx0 + (x0) dx1
}
