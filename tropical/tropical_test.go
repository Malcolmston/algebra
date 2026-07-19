package tropical

import (
	"fmt"
	"math"
	"testing"
)

const tol = 1e-9

func inf() float64  { return math.Inf(1) }
func ninf() float64 { return math.Inf(-1) }

func approx(a, b float64) bool { return closeScalar(a, b, tol) }

// ---------------------------------------------------------------------------
// Scalars and semirings
// ---------------------------------------------------------------------------

func TestScalarOps(t *testing.T) {
	tests := []struct {
		name string
		got  float64
		want float64
	}{
		{"minplus add", MinPlusAdd(3, 5), 3},
		{"maxplus add", MaxPlusAdd(3, 5), 5},
		{"minplus mul", MinPlusMul(3, 5), 8},
		{"maxplus mul", MaxPlusMul(3, 5), 8},
		{"minplus mul zero", MinPlusMul(inf(), 5), inf()},
		{"maxplus mul zero", MaxPlusMul(ninf(), 5), ninf()},
		{"minplus div", MinPlusDiv(8, 3), 5},
		{"maxplus div", MaxPlusDiv(8, 3), 5},
		{"minplus pow", MinPlusPow(2, 3), 6},
		{"maxplus pow", MaxPlusPow(2, 3), 6},
		{"minplus pow0", MinPlusPow(7, 0), 0},
		{"minplus star pos", MinPlusStar(2), 0},
		{"minplus star neg", MinPlusStar(-2), ninf()},
		{"maxplus star neg", MaxPlusStar(-2), 0},
		{"maxplus star pos", MaxPlusStar(2), inf()},
		{"minplus sum", MinPlusSum(4, 1, 9), 1},
		{"maxplus sum", MaxPlusSum(4, 1, 9), 9},
		{"minplus prod", MinPlusProd(4, 1, 9), 14},
	}
	for _, tc := range tests {
		if !approx(tc.got, tc.want) {
			t.Errorf("%s: got %v want %v", tc.name, tc.got, tc.want)
		}
	}
}

func TestSemiringBasics(t *testing.T) {
	mp := MinPlusSemiring()
	xp := MaxPlusSemiring()
	if mp.Zero() != inf() || xp.Zero() != ninf() {
		t.Fatalf("zeros wrong: %v %v", mp.Zero(), xp.Zero())
	}
	if mp.One() != 0 || xp.One() != 0 {
		t.Fatalf("ones wrong")
	}
	if !mp.IsMinPlus() || !xp.IsMaxPlus() {
		t.Fatalf("kind predicates wrong")
	}
	if mp.Dual().Kind() != MaxPlus || xp.Dual().Kind() != MinPlus {
		t.Fatalf("dual wrong")
	}
	if !mp.AtLeastAsGood(2, 5) || mp.AtLeastAsGood(5, 2) {
		t.Fatalf("minplus AtLeastAsGood wrong")
	}
	if !xp.AtLeastAsGood(5, 2) || xp.AtLeastAsGood(2, 5) {
		t.Fatalf("maxplus AtLeastAsGood wrong")
	}
	if MinPlus.String() != "min-plus" || MaxPlus.String() != "max-plus" {
		t.Fatalf("kind string wrong")
	}
}

// ---------------------------------------------------------------------------
// Vectors
// ---------------------------------------------------------------------------

func TestVector(t *testing.T) {
	v := MinPlusVector([]float64{3, 1, 4})
	w := MinPlusVector([]float64{2, 5, 1})
	if !v.Add(w).Equal(MinPlusVector([]float64{2, 1, 1})) {
		t.Errorf("vector add wrong: %v", v.Add(w))
	}
	if got := v.Dot(w); !approx(got, 5) { // min(3+2,1+5,4+1)=min(5,6,5)=5
		t.Errorf("dot got %v want 5", got)
	}
	if !approx(v.Sum(), 1) {
		t.Errorf("sum got %v want 1", v.Sum())
	}
	if !approx(v.Prod(), 8) {
		t.Errorf("prod got %v want 8", v.Prod())
	}
	if v.ScalarMul(10).At(0) != 13 {
		t.Errorf("scalar mul wrong")
	}
	if v.Argmin() != 1 || v.Argmax() != 2 {
		t.Errorf("argmin/argmax wrong: %d %d", v.Argmin(), v.Argmax())
	}
	u := UnitVector(MaxPlusSemiring(), 3, 1)
	if u.At(1) != 0 || u.At(0) != ninf() {
		t.Errorf("unit vector wrong: %v", u)
	}
}

// ---------------------------------------------------------------------------
// Matrix arithmetic
// ---------------------------------------------------------------------------

func TestMatMul(t *testing.T) {
	A := MinPlusMatrix([][]float64{{1, 2}, {3, 4}})
	B := MinPlusMatrix([][]float64{{5, 6}, {7, 8}})
	C, err := A.Mul(B)
	if err != nil {
		t.Fatal(err)
	}
	// (0,0) = min(1+5, 2+7)=min(6,9)=6; (0,1)=min(1+6,2+8)=7
	// (1,0) = min(3+5,4+7)=8;           (1,1)=min(3+6,4+8)=9
	want := MinPlusMatrix([][]float64{{6, 7}, {8, 9}})
	if !C.EqualTol(want, tol) {
		t.Errorf("matmul got %v want %v", C, want)
	}
	// identity behaves as unit
	I := MinPlusIdentity(2)
	AI, _ := A.Mul(I)
	if !AI.EqualTol(A, tol) {
		t.Errorf("A*I != A: %v", AI)
	}
}

func TestMatPow(t *testing.T) {
	A := MaxPlusMatrix([][]float64{{0, 2}, {1, 0}})
	P, err := A.Pow(3)
	if err != nil {
		t.Fatal(err)
	}
	// verify against repeated multiplication
	M := A.Clone()
	M2, _ := M.Mul(A)
	M3, _ := M2.Mul(A)
	if !P.EqualTol(M3, tol) {
		t.Errorf("Pow(3) got %v want %v", P, M3)
	}
	P0, _ := A.Pow(0)
	if !P0.IsIdentity(tol) {
		t.Errorf("Pow(0) not identity")
	}
}

func TestTransposeTrace(t *testing.T) {
	A := MinPlusMatrix([][]float64{{1, 2, 3}, {4, 5, 6}})
	AT := A.Transpose()
	if AT.Rows() != 3 || AT.Cols() != 2 || AT.At(2, 1) != 6 {
		t.Errorf("transpose wrong: %v", AT)
	}
	S := MinPlusMatrix([][]float64{{7, 2}, {3, 4}})
	tr, _ := S.Trace()
	if !approx(tr, 4) { // min(7,4)
		t.Errorf("trace got %v want 4", tr)
	}
}

// ---------------------------------------------------------------------------
// Closure / shortest paths
// ---------------------------------------------------------------------------

func TestShortestPaths(t *testing.T) {
	z := inf()
	A := MinPlusMatrix([][]float64{
		{0, 3, z, 7},
		{8, 0, 2, z},
		{5, z, 0, 1},
		{2, z, z, 0},
	})
	sp, err := A.ShortestPaths()
	if err != nil {
		t.Fatal(err)
	}
	want := MinPlusMatrix([][]float64{
		{0, 3, 5, 6},
		{5, 0, 2, 3},
		{3, 6, 0, 1},
		{2, 5, 7, 0},
	})
	if !sp.EqualTol(want, tol) {
		t.Errorf("shortest paths got\n%v\nwant\n%v", sp, want)
	}
}

func TestNegativeCycleDetection(t *testing.T) {
	// A 2-cycle 0->1->0 with total weight -1 is a negative cycle.
	A := MinPlusMatrix([][]float64{{inf(), 1}, {-2, inf()}})
	if _, err := A.Star(); err != ErrDivergent {
		t.Errorf("expected ErrDivergent, got %v", err)
	}
	bad, _ := A.HasBadCycle()
	if !bad {
		t.Errorf("HasBadCycle should be true")
	}
	// A non-negative version converges.
	B := MinPlusMatrix([][]float64{{inf(), 1}, {2, inf()}})
	if _, err := B.Star(); err != nil {
		t.Errorf("unexpected error: %v", err)
	}
}

func TestStarSeriesMatchesStar(t *testing.T) {
	A := MinPlusMatrix([][]float64{
		{inf(), 3, inf()},
		{inf(), inf(), 2},
		{1, inf(), inf()},
	})
	star, err := A.Star()
	if err != nil {
		t.Fatal(err)
	}
	series, _ := A.StarSeries(A.Rows() - 1)
	if !star.EqualTol(series, tol) {
		t.Errorf("star != series:\n%v\n%v", star, series)
	}
}

// ---------------------------------------------------------------------------
// Eigenvalues (Karp)
// ---------------------------------------------------------------------------

func TestCycleMean(t *testing.T) {
	M := MaxPlusMatrix([][]float64{
		{ninf(), 3},
		{2, ninf()},
	})
	lam, ok, err := M.MaxCycleMean()
	if err != nil || !ok || !approx(lam, 2.5) {
		t.Errorf("max cycle mean got %v ok=%v err=%v want 2.5", lam, ok, err)
	}
	// Min cycle mean of a matrix with self-loops.
	N := MinPlusMatrix([][]float64{
		{4, 1},
		{1, 5},
	})
	// cycles: self loop 0 (mean 4), self loop 1 (mean 5), 0-1-0 (1+1)/2=1.
	mn, ok2, _ := N.MinCycleMean()
	if !ok2 || !approx(mn, 1) {
		t.Errorf("min cycle mean got %v want 1", mn)
	}
	// DAG has no cycle.
	D := MaxPlusMatrix([][]float64{{ninf(), 2}, {ninf(), ninf()}})
	if _, ok3, _ := D.MaxCycleMean(); ok3 {
		t.Errorf("DAG should report no cycle")
	}
}

func TestEigenvector(t *testing.T) {
	M := MaxPlusMatrix([][]float64{
		{ninf(), 3},
		{2, ninf()},
	})
	v, lam, ok, err := M.Eigenvector()
	if err != nil || !ok {
		t.Fatalf("eigenvector failed: ok=%v err=%v", ok, err)
	}
	Av, _ := M.MulVec(v)
	lv := v.ScalarMul(lam)
	if !Av.EqualTol(lv, 1e-7) {
		t.Errorf("A*v=%v but lam*v=%v", Av, lv)
	}
}

// ---------------------------------------------------------------------------
// Determinant / permanent / assignment
// ---------------------------------------------------------------------------

func TestAssignment(t *testing.T) {
	C := MinPlusMatrix([][]float64{
		{4, 2, 8},
		{4, 3, 7},
		{3, 1, 6},
	})
	val, perm, err := C.OptimalAssignment()
	if err != nil {
		t.Fatal(err)
	}
	if !approx(val, 12) {
		t.Errorf("assignment value got %v want 12", val)
	}
	if C.AssignmentValue(perm) != val {
		t.Errorf("AssignmentValue disagrees with OptimalAssignment")
	}
	brute, _ := C.PermanentBrute()
	if !approx(brute, val) {
		t.Errorf("brute %v != hungarian %v", brute, val)
	}
	// Max-plus permanent is the max-weight matching.
	Cx := MaxPlusMatrix([][]float64{
		{4, 2, 8},
		{4, 3, 7},
		{3, 1, 6},
	})
	mx, _ := Cx.Permanent()
	mb, _ := Cx.PermanentBrute()
	if !approx(mx, mb) {
		t.Errorf("maxplus permanent %v != brute %v", mx, mb)
	}
}

func TestSingular(t *testing.T) {
	// Two permutations attain the same optimum -> singular.
	A := MinPlusMatrix([][]float64{
		{0, 0},
		{0, 0},
	})
	s, err := A.IsTropicallySingular(tol)
	if err != nil || !s {
		t.Errorf("expected singular, got %v err=%v", s, err)
	}
	// Unique optimum -> non-singular.
	B := MinPlusMatrix([][]float64{
		{1, 5},
		{5, 1},
	})
	s2, _ := B.IsTropicallySingular(tol)
	if s2 {
		t.Errorf("expected non-singular")
	}
}

func TestCofactor(t *testing.T) {
	A := MinPlusMatrix([][]float64{
		{1, 2, 3},
		{4, 5, 6},
		{7, 8, 9},
	})
	// Minor deleting row0,col0 = [[5,6],[8,9]]; min assignment = min(5+9,6+8)=14.
	c, err := A.Cofactor(0, 0)
	if err != nil || !approx(c, 14) {
		t.Errorf("cofactor got %v want 14 err=%v", c, err)
	}
}

// ---------------------------------------------------------------------------
// Polynomials, Newton polygon, roots
// ---------------------------------------------------------------------------

func TestPolyEvalMul(t *testing.T) {
	p := MinPlusPoly([]float64{3, 0, 0}) // min(3, x, 2x)
	tests := []struct {
		x, want float64
	}{
		{-2, -4}, // min(3,-2,-4)
		{0, 0},   // min(3,0,0)
		{5, 3},   // min(3,5,10)
	}
	for _, tc := range tests {
		if got := p.Eval(tc.x); !approx(got, tc.want) {
			t.Errorf("Eval(%v)=%v want %v", tc.x, got, tc.want)
		}
	}
	// (x (+) 0)(x (+) 3) as convolution.
	f1 := MinPlusPoly([]float64{0, 0}) // min(x,0)
	f2 := MinPlusPoly([]float64{3, 0}) // min(x,3)
	prod := f1.Mul(f2)
	// product = min(2x, x+3, x, 3) = min(2x, x, 3) with coeffs [3,0,0]
	if !prod.Equal(MinPlusPoly([]float64{3, 0, 0})) {
		t.Errorf("poly mul got coeffs %v", prod.Coeffs())
	}
}

func TestPolyRoots(t *testing.T) {
	tests := []struct {
		name  string
		poly  Poly
		roots []Root
	}{
		{
			"minplus double",
			MinPlusPoly([]float64{0, 0, 0}), // min(0,x,2x)
			[]Root{{Value: 0, Multiplicity: 2}},
		},
		{
			"minplus two simple",
			MinPlusPoly([]float64{3, 0, 0}), // min(3,x,2x)
			[]Root{{Value: 0, Multiplicity: 1}, {Value: 3, Multiplicity: 1}},
		},
		{
			"maxplus double",
			MaxPlusPoly([]float64{3, 0, 0}), // max(3,x,2x)
			[]Root{{Value: 1.5, Multiplicity: 2}},
		},
	}
	for _, tc := range tests {
		got := tc.poly.Roots()
		if len(got) != len(tc.roots) {
			t.Errorf("%s: got %d roots %v want %d", tc.name, len(got), got, len(tc.roots))
			continue
		}
		for i := range got {
			if !approx(got[i].Value, tc.roots[i].Value) || got[i].Multiplicity != tc.roots[i].Multiplicity {
				t.Errorf("%s: root %d got %+v want %+v", tc.name, i, got[i], tc.roots[i])
			}
		}
		// Each returned root must make the polynomial tropically vanish.
		for _, r := range got {
			if !tc.poly.IsRoot(r.Value, 1e-7) {
				t.Errorf("%s: %v is not a tropical root", tc.name, r.Value)
			}
		}
	}
}

func TestFromRoots(t *testing.T) {
	roots := []float64{-1, 2, 5}
	p := MinPlusFromRoots(roots)
	got := p.RootValues()
	want := []float64{-1, 2, 5}
	if len(got) != len(want) {
		t.Fatalf("got %v want %v", got, want)
	}
	for i := range want {
		if !approx(got[i], want[i]) {
			t.Errorf("root %d got %v want %v", i, got[i], want[i])
		}
	}
	if !p.IsMonic() {
		t.Errorf("FromRoots should be monic")
	}
}

// ---------------------------------------------------------------------------
// Linear systems
// ---------------------------------------------------------------------------

func TestSolveExact(t *testing.T) {
	A := MaxPlusMatrix([][]float64{
		{1, 2},
		{3, 0},
	})
	b := MaxPlusVector([]float64{5, 6})
	x, ok, err := A.SolveExact(b, tol)
	if err != nil {
		t.Fatal(err)
	}
	Ax, _ := A.MulVec(x)
	if ok && !Ax.EqualTol(b, tol) {
		t.Errorf("claims exact but A*x=%v != b=%v", Ax, b)
	}
	// Greatest subsolution always satisfies A*x <= b.
	for i := 0; i < b.Len(); i++ {
		if Ax.At(i) > b.At(i)+tol {
			t.Errorf("subsolution violated at %d: %v > %v", i, Ax.At(i), b.At(i))
		}
	}
}

func TestSolveAffine(t *testing.T) {
	// x = A x (+) b, min-plus. A strictly upper triangular so star converges.
	A := MinPlusMatrix([][]float64{
		{inf(), 2},
		{inf(), inf()},
	})
	b := MinPlusVector([]float64{0, 5})
	x, err := A.SolveAffine(b)
	if err != nil {
		t.Fatal(err)
	}
	// Check the fixed point x == A x (+) b.
	Ax, _ := A.MulVec(x)
	fp := Ax.Add(b)
	if !x.EqualTol(fp, tol) {
		t.Errorf("not a fixed point: x=%v Ax+b=%v", x, fp)
	}
	if !x.EqualTol(MinPlusVector([]float64{0, 5}), tol) {
		t.Errorf("affine solution got %v want [0 5]", x)
	}
}

func TestResidual(t *testing.T) {
	A := MaxPlusMatrix([][]float64{{1, 2}, {3, 4}})
	B := MaxPlusMatrix([][]float64{{5, 6}, {7, 8}})
	X, err := A.LeftResidual(B)
	if err != nil {
		t.Fatal(err)
	}
	// A (*) X must be <= B everywhere (greatest subsolution property).
	AX, _ := A.Mul(X)
	for i := 0; i < B.Rows(); i++ {
		for j := 0; j < B.Cols(); j++ {
			if AX.At(i, j) > B.At(i, j)+tol {
				t.Errorf("A*X > B at (%d,%d): %v > %v", i, j, AX.At(i, j), B.At(i, j))
			}
		}
	}
}

// ---------------------------------------------------------------------------
// Structural helpers
// ---------------------------------------------------------------------------

func TestKroneckerAndStack(t *testing.T) {
	A := MinPlusMatrix([][]float64{{1, 2}})
	B := MinPlusMatrix([][]float64{{0, 10}})
	K, err := A.KroneckerProduct(B)
	if err != nil {
		t.Fatal(err)
	}
	// row = [1+0,1+10,2+0,2+10] = [1,11,2,12]
	if !K.EqualTol(MinPlusMatrix([][]float64{{1, 11, 2, 12}}), tol) {
		t.Errorf("kronecker got %v", K)
	}
	H, _ := A.HStack(B)
	if H.Cols() != 4 || H.At(0, 2) != 0 {
		t.Errorf("hstack wrong: %v", H)
	}
	D, _ := A.Transpose().DirectSum(B.Transpose())
	if D.Rows() != 4 || D.Cols() != 2 || !D.sr.IsZero(D.At(0, 1)) {
		t.Errorf("direct sum wrong: %v", D)
	}
}

func TestDual(t *testing.T) {
	A := MinPlusMatrix([][]float64{{1, inf()}, {-2, 3}})
	D := A.Dual()
	if !D.Semiring().IsMaxPlus() {
		t.Errorf("dual semiring wrong")
	}
	if D.At(0, 0) != -1 || D.At(0, 1) != ninf() || D.At(1, 0) != 2 {
		t.Errorf("dual entries wrong: %v", D)
	}
}

func TestSingleSource(t *testing.T) {
	z := inf()
	A := MinPlusMatrix([][]float64{
		{0, 3, z},
		{z, 0, 2},
		{z, z, 0},
	})
	d, err := A.SingleSourceShortestPaths(0)
	if err != nil {
		t.Fatal(err)
	}
	if !d.EqualTol(MinPlusVector([]float64{0, 3, 5}), tol) {
		t.Errorf("single source got %v want [0 3 5]", d)
	}
}

// ---------------------------------------------------------------------------
// Example
// ---------------------------------------------------------------------------

func ExampleMatrix_ShortestPaths() {
	inf := MinPlusZero()
	// Weighted digraph as a min-plus adjacency matrix; inf means "no edge".
	A := MinPlusMatrix([][]float64{
		{0, 3, inf},
		{inf, 0, 2},
		{1, inf, 0},
	})
	sp, _ := A.ShortestPaths()
	// Shortest distance from node 0 to node 2 is 0->1->2 = 5.
	fmt.Println(sp.At(0, 2))
	// Output: 5
}

func ExamplePoly_Roots() {
	// Tropical polynomial min(3, x, 2x) over the min-plus semiring.
	p := MinPlusPoly([]float64{3, 0, 0})
	for _, r := range p.Roots() {
		fmt.Printf("root %.0f (mult %d)\n", r.Value, r.Multiplicity)
	}
	// Output:
	// root 0 (mult 1)
	// root 3 (mult 1)
}
