package matrix

import (
	"errors"
	"math"
	"testing"

	"github.com/malcolmston/algebra"
)

// numEqual reports whether the numeric value of e is within tol of want.
func numEqual(t *testing.T, e algebra.Expr, want, tol float64) {
	t.Helper()
	got, err := algebra.Evalf(e)
	if err != nil {
		t.Fatalf("Evalf(%s): %v", e.String(), err)
	}
	if math.Abs(got-want) > tol {
		t.Fatalf("got %v (%s), want %v", got, e.String(), want)
	}
}

func TestDet3x3(t *testing.T) {
	// Known determinant: det = 1*(5*9-6*8) - 2*(4*9-6*7) + 3*(4*8-5*7)
	//                        = 1*(-3) - 2*(-6) + 3*(-3) = -3 +12 -9 = 0.
	m := FromInts([][]int64{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}})
	d, err := m.Det()
	if err != nil {
		t.Fatal(err)
	}
	if !d.Equal(algebra.Int(0)) {
		t.Fatalf("det = %s, want 0", d.String())
	}

	// A matrix with a known non-zero determinant.
	m2 := FromInts([][]int64{{2, 1, 1}, {1, 3, 2}, {1, 0, 0}})
	// det = 1*(1*2-1*3) - 0 + 0 (expand along last row) = 1*(2-3) = -1.
	d2, _ := m2.Det()
	if !d2.Equal(algebra.Int(-1)) {
		t.Fatalf("det2 = %s, want -1", d2.String())
	}
}

func TestDetLaplaceLarge(t *testing.T) {
	// 4x4 upper-triangular: determinant is product of the diagonal = 24.
	m := FromInts([][]int64{
		{2, 5, 7, 1},
		{0, 3, 4, 9},
		{0, 0, 1, 6},
		{0, 0, 0, 4},
	})
	d, _ := m.Det()
	if !d.Equal(algebra.Int(24)) {
		t.Fatalf("det = %s, want 24", d.String())
	}
}

func TestInverseTimesSelfIsIdentity(t *testing.T) {
	m := FromInts([][]int64{{4, 7}, {2, 6}})
	inv, err := m.Inverse()
	if err != nil {
		t.Fatal(err)
	}
	prod, err := m.Mul(inv)
	if err != nil {
		t.Fatal(err)
	}
	if !prod.Equal(Identity(2)) {
		t.Fatalf("A*A^-1 =\n%s\nwant identity", prod.String())
	}
	// And the other order.
	prod2, _ := inv.Mul(m)
	if !prod2.Equal(Identity(2)) {
		t.Fatalf("A^-1*A =\n%s\nwant identity", prod2.String())
	}
}

func TestInverse3x3(t *testing.T) {
	m := FromInts([][]int64{{1, 2, 3}, {0, 1, 4}, {5, 6, 0}})
	inv, err := m.Inverse()
	if err != nil {
		t.Fatal(err)
	}
	prod, _ := m.Mul(inv)
	if !prod.Equal(Identity(3)) {
		t.Fatalf("A*A^-1 =\n%s\nwant identity", prod.String())
	}
}

func TestSingularInverse(t *testing.T) {
	m := FromInts([][]int64{{1, 2}, {2, 4}})
	if _, err := m.Inverse(); !errors.Is(err, ErrSingular) {
		t.Fatalf("err = %v, want ErrSingular", err)
	}
}

func TestSymbolicInverse(t *testing.T) {
	a := algebra.Sym("a")
	b := algebra.Sym("b")
	c := algebra.Sym("c")
	d := algebra.Sym("d")
	m := FromExpr([][]algebra.Expr{{a, b}, {c, d}})
	inv, err := m.Inverse()
	if err != nil {
		t.Fatal(err)
	}
	// The inverse of [[a,b],[c,d]] is adj/det = [[d,-b],[-c,a]]/(a*d-b*c). The
	// parent simplifier does not combine rational functions, so we verify the
	// inverse against this exact closed form rather than expecting A*A^-1 to
	// collapse to the literal identity.
	det := algebra.Add(algebra.Mul(a, d), algebra.Mul(algebra.Int(-1), b, c))
	invDet := algebra.Pow(det, algebra.Int(-1))
	want := FromExpr([][]algebra.Expr{
		{algebra.Mul(d, invDet), algebra.Mul(algebra.Int(-1), b, invDet)},
		{algebra.Mul(algebra.Int(-1), c, invDet), algebra.Mul(a, invDet)},
	})
	if !inv.Equal(want) {
		t.Fatalf("symbolic inverse =\n%s\nwant\n%s", inv.String(), want.String())
	}
	// Numeric substitution sanity check: with a,b,c,d = 1,2,3,4 the product
	// A*A^-1 must evaluate to the identity.
	prod, _ := m.Mul(inv)
	env := map[string]float64{"a": 1, "b": 2, "c": 3, "d": 4}
	for i := 0; i < 2; i++ {
		for j := 0; j < 2; j++ {
			v, err := algebra.Eval(prod.At(i, j), env)
			if err != nil {
				t.Fatalf("eval entry (%d,%d): %v", i, j, err)
			}
			want := 0.0
			if i == j {
				want = 1.0
			}
			if math.Abs(v-want) > 1e-12 {
				t.Fatalf("A*A^-1 entry (%d,%d) = %v, want %v", i, j, v, want)
			}
		}
	}
}

func TestSolveUnique(t *testing.T) {
	// 2x + y = 5 ; x - y = 1  ->  x=2, y=1.
	a := FromInts([][]int64{{2, 1}, {1, -1}})
	b := VectorFromInts(5, 1)
	x, err := Solve(a, b)
	if err != nil {
		t.Fatal(err)
	}
	if !x.At(0).Equal(algebra.Int(2)) || !x.At(1).Equal(algebra.Int(1)) {
		t.Fatalf("solution = %s, want [2, 1]", x.String())
	}
}

func TestSolveRationalSolution(t *testing.T) {
	// x + y = 1 ; x - y = 0 -> x=y=1/2.
	a := FromInts([][]int64{{1, 1}, {1, -1}})
	b := VectorFromInts(1, 0)
	x, _ := Solve(a, b)
	if !x.At(0).Equal(algebra.Rat(1, 2)) || !x.At(1).Equal(algebra.Rat(1, 2)) {
		t.Fatalf("solution = %s, want [1/2, 1/2]", x.String())
	}
}

func TestSolveInconsistent(t *testing.T) {
	a := FromInts([][]int64{{1, 1}, {1, 1}})
	b := VectorFromInts(1, 2)
	if _, err := Solve(a, b); !errors.Is(err, ErrInconsistent) {
		t.Fatalf("err = %v, want ErrInconsistent", err)
	}
}

func TestSolveUnderdetermined(t *testing.T) {
	a := FromInts([][]int64{{1, 1}, {2, 2}})
	b := VectorFromInts(1, 2)
	if _, err := Solve(a, b); !errors.Is(err, ErrUnderdetermined) {
		t.Fatalf("err = %v, want ErrUnderdetermined", err)
	}
}

func TestRankAndRREF(t *testing.T) {
	m := FromInts([][]int64{{1, 2, 3}, {2, 4, 6}, {1, 1, 1}})
	if got := m.Rank(); got != 2 {
		t.Fatalf("rank = %d, want 2", got)
	}
	if got := Identity(3).Rank(); got != 3 {
		t.Fatalf("identity rank = %d, want 3", got)
	}
	rref, _ := m.RREF()
	// Third row should be all zeros in RREF.
	for j := 0; j < 3; j++ {
		if !rref.At(2, j).Equal(algebra.Int(0)) {
			t.Fatalf("rref row 2 not zero:\n%s", rref.String())
		}
	}
}

func TestArithmetic(t *testing.T) {
	a := FromInts([][]int64{{1, 2}, {3, 4}})
	b := FromInts([][]int64{{5, 6}, {7, 8}})
	sum, _ := a.Add(b)
	if !sum.Equal(FromInts([][]int64{{6, 8}, {10, 12}})) {
		t.Fatalf("sum wrong:\n%s", sum.String())
	}
	diff, _ := b.Sub(a)
	if !diff.Equal(FromInts([][]int64{{4, 4}, {4, 4}})) {
		t.Fatalf("diff wrong:\n%s", diff.String())
	}
	prod, _ := a.Mul(b)
	// [[1*5+2*7, 1*6+2*8],[3*5+4*7,3*6+4*8]] = [[19,22],[43,50]].
	if !prod.Equal(FromInts([][]int64{{19, 22}, {43, 50}})) {
		t.Fatalf("prod wrong:\n%s", prod.String())
	}
	if !a.Neg().Equal(FromInts([][]int64{{-1, -2}, {-3, -4}})) {
		t.Fatalf("neg wrong")
	}
	if !a.ScalarMul(algebra.Int(2)).Equal(FromInts([][]int64{{2, 4}, {6, 8}})) {
		t.Fatalf("scalarmul wrong")
	}
}

func TestMulVec(t *testing.T) {
	a := FromInts([][]int64{{1, 2}, {3, 4}})
	v := VectorFromInts(1, 1)
	got, _ := a.MulVec(v)
	if !got.Equal(VectorFromInts(3, 7)) {
		t.Fatalf("MulVec = %s, want [3, 7]", got.String())
	}
}

func TestTransposeAndPow(t *testing.T) {
	a := FromInts([][]int64{{1, 2, 3}, {4, 5, 6}})
	tr := a.Transpose()
	if tr.Rows() != 3 || tr.Cols() != 2 || !tr.At(2, 1).Equal(algebra.Int(6)) {
		t.Fatalf("transpose wrong:\n%s", tr.String())
	}
	m := FromInts([][]int64{{2, 0}, {0, 3}})
	p3, _ := m.Pow(3)
	if !p3.Equal(FromInts([][]int64{{8, 0}, {0, 27}})) {
		t.Fatalf("pow wrong:\n%s", p3.String())
	}
	p0, _ := m.Pow(0)
	if !p0.Equal(Identity(2)) {
		t.Fatalf("pow0 not identity")
	}
}

func TestKron(t *testing.T) {
	a := FromInts([][]int64{{1, 0}, {0, 1}})
	b := FromInts([][]int64{{1, 2}, {3, 4}})
	k := a.Kron(b)
	if k.Rows() != 4 || k.Cols() != 4 {
		t.Fatalf("kron shape %dx%d", k.Rows(), k.Cols())
	}
	// Identity ⊗ B is block-diagonal with B on the diagonal.
	want := FromInts([][]int64{
		{1, 2, 0, 0},
		{3, 4, 0, 0},
		{0, 0, 1, 2},
		{0, 0, 3, 4},
	})
	if !k.Equal(want) {
		t.Fatalf("kron wrong:\n%s", k.String())
	}
}

func TestTraceMinorCofactorAdjugate(t *testing.T) {
	m := FromInts([][]int64{{1, 2, 3}, {4, 5, 6}, {7, 8, 10}})
	tr, _ := m.Trace()
	if !tr.Equal(algebra.Int(16)) {
		t.Fatalf("trace = %s, want 16", tr.String())
	}
	// Minor(0,0) = det[[5,6],[8,10]] = 50-48 = 2.
	minor, _ := m.Minor(0, 0)
	if !minor.Equal(algebra.Int(2)) {
		t.Fatalf("minor = %s, want 2", minor.String())
	}
	// Cofactor(0,1) = -det[[4,6],[7,10]] = -(40-42) = 2.
	cof, _ := m.Cofactor(0, 1)
	if !cof.Equal(algebra.Int(2)) {
		t.Fatalf("cofactor = %s, want 2", cof.String())
	}
	// A * adj(A) = det(A) * I.
	adj, _ := m.Adjugate()
	prod, _ := m.Mul(adj)
	det, _ := m.Det()
	want := Identity(3).ScalarMul(det)
	if !prod.Equal(want) {
		t.Fatalf("A*adj(A) =\n%s\nwant det*I", prod.String())
	}
}

func TestVectorOps(t *testing.T) {
	v := VectorFromInts(3, 4)
	d, _ := v.Dot(v)
	if !d.Equal(algebra.Int(25)) {
		t.Fatalf("dot = %s, want 25", d.String())
	}
	if !v.Norm().Equal(algebra.Int(5)) {
		t.Fatalf("norm = %s, want 5", v.Norm().String())
	}
	unit, _ := v.Normalize()
	if !unit.At(0).Equal(algebra.Rat(3, 5)) || !unit.At(1).Equal(algebra.Rat(4, 5)) {
		t.Fatalf("normalize = %s, want [3/5, 4/5]", unit.String())
	}
	// Non-perfect-square norm stays symbolic: |[1,1]| = sqrt(2).
	w := VectorFromInts(1, 1)
	numEqual(t, w.Norm(), math.Sqrt2, 1e-12)
}

func TestCross(t *testing.T) {
	e1 := VectorFromInts(1, 0, 0)
	e2 := VectorFromInts(0, 1, 0)
	e3 := VectorFromInts(0, 0, 1)
	c, err := e1.Cross(e2)
	if err != nil {
		t.Fatal(err)
	}
	if !c.Equal(e3) {
		t.Fatalf("e1 x e2 = %s, want [0, 0, 1]", c.String())
	}
	// e2 x e3 = e1, e3 x e1 = e2.
	c23, _ := e2.Cross(e3)
	if !c23.Equal(e1) {
		t.Fatalf("e2 x e3 = %s, want e1", c23.String())
	}
	c31, _ := e3.Cross(e1)
	if !c31.Equal(e2) {
		t.Fatalf("e3 x e1 = %s, want e2", c31.String())
	}
}

func TestAngle(t *testing.T) {
	a := VectorFromInts(1, 0)
	b := VectorFromInts(0, 1)
	ang, err := a.Angle(b)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(ang-math.Pi/2) > 1e-12 {
		t.Fatalf("angle = %v, want pi/2", ang)
	}
}

func TestCharPoly(t *testing.T) {
	m := FromInts([][]int64{{2, 0}, {0, 3}})
	p, err := m.CharPoly()
	if err != nil {
		t.Fatal(err)
	}
	// det(A - λI) = (2-λ)(3-λ) = λ^2 - 5λ + 6.
	want := algebra.Add(
		algebra.Pow(Lambda, algebra.Int(2)),
		algebra.Mul(algebra.Int(-5), Lambda),
		algebra.Int(6),
	)
	if !simp(p).Equal(simp(want)) {
		t.Fatalf("charpoly = %s, want lambda^2 - 5*lambda + 6", p.String())
	}
}

func TestEigenvalues2x2Diagonal(t *testing.T) {
	m := FromInts([][]int64{{2, 0}, {0, 3}})
	ev, err := m.Eigenvalues()
	if err != nil {
		t.Fatal(err)
	}
	if len(ev) != 2 {
		t.Fatalf("got %d eigenvalues, want 2", len(ev))
	}
	got := map[string]bool{ev[0].String(): true, ev[1].String(): true}
	if !got["2"] || !got["3"] {
		t.Fatalf("eigenvalues = %v, want {2, 3}", []string{ev[0].String(), ev[1].String()})
	}
}

func TestEigenvalues2x2Symbolic(t *testing.T) {
	// [[0,1],[1,0]] has eigenvalues ±1.
	m := FromInts([][]int64{{0, 1}, {1, 0}})
	ev, _ := m.Eigenvalues()
	if len(ev) != 2 {
		t.Fatalf("got %d eigenvalues, want 2", len(ev))
	}
	vals := map[string]bool{}
	for _, e := range ev {
		vals[e.String()] = true
	}
	if !vals["1"] || !vals["-1"] {
		t.Fatalf("eigenvalues = %v, want {1, -1}", []string{ev[0].String(), ev[1].String()})
	}
}

func TestEigenvalues3x3Rational(t *testing.T) {
	// Lower-triangular: eigenvalues are the diagonal {2, 3, 4}.
	m := FromInts([][]int64{{2, 0, 0}, {1, 3, 0}, {5, 6, 4}})
	ev, err := m.Eigenvalues()
	if err != nil {
		t.Fatal(err)
	}
	if len(ev) != 3 {
		t.Fatalf("got %d eigenvalues, want 3", len(ev))
	}
	vals := map[string]bool{}
	for _, e := range ev {
		vals[e.String()] = true
	}
	for _, w := range []string{"2", "3", "4"} {
		if !vals[w] {
			t.Fatalf("eigenvalues = %v, want {2,3,4}", vals)
		}
	}
}

func TestEigenvalues3x3Irrational(t *testing.T) {
	// Tridiagonal [[2,-1,0],[-1,2,-1],[0,-1,2]] has eigenvalues 2, 2±sqrt(2).
	m := FromInts([][]int64{{2, -1, 0}, {-1, 2, -1}, {0, -1, 2}})
	ev, err := m.Eigenvalues()
	if err != nil {
		t.Fatal(err)
	}
	if len(ev) != 3 {
		t.Fatalf("got %d eigenvalues, want 3", len(ev))
	}
	want := []float64{2 - math.Sqrt2, 2, 2 + math.Sqrt2}
	for _, w := range want {
		found := false
		for _, e := range ev {
			if v, err := algebra.Evalf(e); err == nil && math.Abs(v-w) < 1e-9 {
				found = true
				break
			}
		}
		if !found {
			t.Fatalf("missing eigenvalue %v in %v", w, exprStrings(ev))
		}
	}
}

func exprStrings(es []algebra.Expr) []string {
	out := make([]string, len(es))
	for i, e := range es {
		out[i] = e.String()
	}
	return out
}

func TestEigenvalues3x3NumericOneReal(t *testing.T) {
	// Companion matrix of x^3 - x - 1 (plastic number): one real root ~1.3247,
	// two complex. No rational root, so the numeric fallback returns the single
	// real eigenvalue as a float.
	m := FromInts([][]int64{{0, 0, 1}, {1, 0, 1}, {0, 1, 0}})
	ev, err := m.Eigenvalues()
	if err != nil {
		t.Fatal(err)
	}
	if len(ev) != 1 {
		t.Fatalf("got %d real eigenvalues, want 1: %v", len(ev), exprStrings(ev))
	}
	numEqual(t, ev[0], 1.3247179572, 1e-6)
}

func TestEigenvalues3x3NumericThreeReal(t *testing.T) {
	// Companion matrix of x^3 - 3x - 1: three distinct irrational real roots.
	m := FromInts([][]int64{{0, 0, 1}, {1, 0, 3}, {0, 1, 0}})
	ev, err := m.Eigenvalues()
	if err != nil {
		t.Fatal(err)
	}
	if len(ev) != 3 {
		t.Fatalf("got %d eigenvalues, want 3: %v", len(ev), exprStrings(ev))
	}
	for _, w := range []float64{1.8793852416, -0.3472963553, -1.5320888862} {
		found := false
		for _, e := range ev {
			if v, err := algebra.Evalf(e); err == nil && math.Abs(v-w) < 1e-6 {
				found = true
			}
		}
		if !found {
			t.Fatalf("missing eigenvalue %v in %v", w, exprStrings(ev))
		}
	}
}

func TestAccessorsAndConversions(t *testing.T) {
	m := FromInts([][]int64{{1, 2}, {3, 4}})
	if !m.Row(0).Equal(VectorFromInts(1, 2)) {
		t.Fatal("Row wrong")
	}
	if !m.Col(1).Equal(VectorFromInts(2, 4)) {
		t.Fatal("Col wrong")
	}
	m.Set(0, 0, algebra.Int(9))
	if !m.At(0, 0).Equal(algebra.Int(9)) {
		t.Fatal("Set/At wrong")
	}
	v := VectorFromInts(5, 6)
	if !v.ColMatrix().Equal(FromInts([][]int64{{5}, {6}})) {
		t.Fatal("ColMatrix wrong")
	}
	if !v.RowMatrix().Equal(FromInts([][]int64{{5, 6}})) {
		t.Fatal("RowMatrix wrong")
	}
	if got := New(2, 2).String(); got == "" {
		t.Fatal("empty String")
	}
}

func TestVectorArithmetic(t *testing.T) {
	a := VectorFromInts(1, 2, 3)
	b := VectorFromInts(4, 5, 6)
	sum, _ := a.Add(b)
	if !sum.Equal(VectorFromInts(5, 7, 9)) {
		t.Fatal("vector Add wrong")
	}
	diff, _ := b.Sub(a)
	if !diff.Equal(VectorFromInts(3, 3, 3)) {
		t.Fatal("vector Sub wrong")
	}
	if !a.Neg().Equal(VectorFromInts(-1, -2, -3)) {
		t.Fatal("vector Neg wrong")
	}
	if !a.ScalarMul(algebra.Int(2)).Equal(VectorFromInts(2, 4, 6)) {
		t.Fatal("vector ScalarMul wrong")
	}
}

func TestEigenvaluesUnsupported(t *testing.T) {
	m := Identity(4)
	if _, err := m.Eigenvalues(); !errors.Is(err, ErrUnsupported) {
		t.Fatalf("err = %v, want ErrUnsupported", err)
	}
}

func TestDimensionErrors(t *testing.T) {
	a := FromInts([][]int64{{1, 2, 3}})
	b := FromInts([][]int64{{1, 2}})
	if _, err := a.Add(b); !errors.Is(err, ErrDimension) {
		t.Fatalf("add err = %v, want ErrDimension", err)
	}
	if _, err := a.Mul(b); !errors.Is(err, ErrDimension) {
		t.Fatalf("mul err = %v, want ErrDimension", err)
	}
	if _, err := a.Det(); !errors.Is(err, ErrNotSquare) {
		t.Fatalf("det err = %v, want ErrNotSquare", err)
	}
}

func TestConstructors(t *testing.T) {
	if !Zeros(2, 2).Equal(New(2, 2)) {
		t.Fatal("Zeros != New")
	}
	if !Ones(1, 3).Equal(FromInts([][]int64{{1, 1, 1}})) {
		t.Fatal("Ones wrong")
	}
	d := Diag(algebra.Int(5), algebra.Int(7))
	if !d.Equal(FromInts([][]int64{{5, 0}, {0, 7}})) {
		t.Fatal("Diag wrong")
	}
	f := FromFloats([][]float64{{1.5, 2.5}})
	if got, _ := f.Floats(); got[0][0] != 1.5 || got[0][1] != 2.5 {
		t.Fatal("FromFloats/Floats round-trip wrong")
	}
}
