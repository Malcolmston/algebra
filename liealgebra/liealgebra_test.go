package liealgebra

import (
	"fmt"
	"math"
	"math/cmplx"
	"testing"
)

const tol = 1e-9

// ---------------------------------------------------------------------------
// Matrix arithmetic
// ---------------------------------------------------------------------------

func TestMatrixBasics(t *testing.T) {
	a, _ := NewMatrixFromRows([][]float64{{1, 2}, {3, 4}})
	b, _ := NewMatrixFromRows([][]float64{{0, 1}, {1, 0}})
	sum, _ := a.Add(b)
	if !sum.Equal(mustM(t, [][]float64{{1, 3}, {4, 4}})) {
		t.Fatalf("Add: %v", sum.Data)
	}
	prod, _ := a.Mul(b)
	if !prod.Equal(mustM(t, [][]float64{{2, 1}, {4, 3}})) {
		t.Fatalf("Mul: %v", prod.Data)
	}
	tr, _ := a.Trace()
	if tr != 5 {
		t.Fatalf("Trace=%v", tr)
	}
	if d, _ := Det(a); math.Abs(d-(-2)) > tol {
		t.Fatalf("Det=%v", d)
	}
	inv, err := Inverse(a)
	if err != nil {
		t.Fatal(err)
	}
	id, _ := a.Mul(inv)
	if !id.ApproxEqual(IdentityMatrix(2), tol) {
		t.Fatalf("A*inv != I: %v", id.Data)
	}
}

func TestMatrixPredicates(t *testing.T) {
	sym, _ := NewMatrixFromRows([][]float64{{1, 2}, {2, 3}})
	if !sym.IsSymmetric(tol) {
		t.Fatal("expected symmetric")
	}
	skew, _ := NewMatrixFromRows([][]float64{{0, 2}, {-2, 0}})
	if !skew.IsAntisymmetric(tol) {
		t.Fatal("expected antisymmetric")
	}
	if !IdentityMatrix(3).IsDiagonal(tol) {
		t.Fatal("identity should be diagonal")
	}
}

func TestSolveAndLstsq(t *testing.T) {
	a, _ := NewMatrixFromRows([][]float64{{2, 1}, {1, 3}})
	x, err := Solve(a, []float64{3, 5})
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(x[0]-0.8) > tol || math.Abs(x[1]-1.4) > tol {
		t.Fatalf("Solve x=%v", x)
	}
	// Overdetermined exact system.
	m, _ := NewMatrixFromRows([][]float64{{1, 0}, {0, 1}, {1, 1}})
	sol, err := SolveLeastSquares(m, []float64{1, 2, 3})
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(sol[0]-1) > tol || math.Abs(sol[1]-2) > tol {
		t.Fatalf("lstsq=%v", sol)
	}
}

// ---------------------------------------------------------------------------
// Complex matrices / Pauli
// ---------------------------------------------------------------------------

func TestPauliAlgebra(t *testing.T) {
	sx, sy, sz := PauliMatrices()
	for _, s := range []*CMatrix{sx, sy, sz} {
		if !s.IsHermitian(tol) {
			t.Fatal("Pauli not Hermitian")
		}
		if !s.IsUnitary(tol) {
			t.Fatal("Pauli not unitary")
		}
		if !s.IsTraceless(tol) {
			t.Fatal("Pauli not traceless")
		}
	}
	// [sx,sy] = 2i sz.
	br, _ := CBracket(sx, sy)
	want := sz.Scale(complex(0, 2))
	if !br.ApproxEqual(want, tol) {
		t.Fatalf("[sx,sy] wrong: %v", br.Data)
	}
	// sx^2 = I.
	sq, _ := sx.Mul(sx)
	if !sq.ApproxEqual(IdentityCMatrix(2), tol) {
		t.Fatalf("sx^2=%v", sq.Data)
	}
}

// ---------------------------------------------------------------------------
// Brackets / Jacobi / structure constants
// ---------------------------------------------------------------------------

func TestSL2Relations(t *testing.T) {
	e, f, h := SL2Generators()
	// [H,E]=2E, [H,F]=-2F, [E,F]=H.
	he, _ := Bracket(h, e)
	if !he.ApproxEqual(e.Scale(2), tol) {
		t.Fatalf("[H,E]=%v", he.Data)
	}
	hf, _ := Bracket(h, f)
	if !hf.ApproxEqual(f.Scale(-2), tol) {
		t.Fatalf("[H,F]=%v", hf.Data)
	}
	ef, _ := Bracket(e, f)
	if !ef.ApproxEqual(h, tol) {
		t.Fatalf("[E,F]=%v", ef.Data)
	}
	if !SatisfiesJacobi(e, f, h, tol) {
		t.Fatal("Jacobi failed for sl2")
	}
}

func TestSO3Structure(t *testing.T) {
	lx, ly, lz := SO3Generators()
	// [Lx,Ly]=Lz.
	br, _ := Bracket(lx, ly)
	if !br.ApproxEqual(lz, tol) {
		t.Fatalf("[Lx,Ly]=%v", br.Data)
	}
	basis := []*Matrix{lx, ly, lz}
	c, err := StructureConstants(basis)
	if err != nil {
		t.Fatal(err)
	}
	// Antisymmetric Levi-Civita: c[2][0][1]=+1, c[2][1][0]=-1.
	if math.Abs(c[2][0][1]-1) > tol || math.Abs(c[2][1][0]+1) > tol {
		t.Fatalf("structure constants wrong: %v %v", c[2][0][1], c[2][1][0])
	}
	if r := JacobiResidualConstants(c); r > tol {
		t.Fatalf("Jacobi residual (constants) = %v", r)
	}
	if !IsClosedUnderBracket(basis, tol) {
		t.Fatal("so3 basis not closed")
	}
}

func TestKillingForms(t *testing.T) {
	// so(3): Killing form = -2 I (negative definite, compact).
	lx, ly, lz := SO3Generators()
	k, err := KillingForm([]*Matrix{lx, ly, lz})
	if err != nil {
		t.Fatal(err)
	}
	if !k.ApproxEqual(IdentityMatrix(3).Scale(-2), tol) {
		t.Fatalf("so3 Killing=%v", k.Data)
	}
	semis, _ := IsSemisimple([]*Matrix{lx, ly, lz}, tol)
	if !semis {
		t.Fatal("so3 should be semisimple")
	}
	// sl(2): K(E,F)=4, K(H,H)=8.
	e, f, h := SL2Generators()
	k2, _ := KillingForm([]*Matrix{e, f, h})
	if math.Abs(k2.At(0, 1)-4) > tol || math.Abs(k2.At(2, 2)-8) > tol {
		t.Fatalf("sl2 Killing=%v", k2.Data)
	}
	// Heisenberg algebra is nilpotent -> Killing form identically zero.
	x, y, z := HeisenbergGenerators()
	k3, _ := KillingForm([]*Matrix{x, y, z})
	if k3.MaxAbs() > tol {
		t.Fatalf("Heisenberg Killing should vanish: %v", k3.Data)
	}
	if IsAbelian([]*Matrix{x, y, z}, tol) {
		t.Fatal("Heisenberg is not abelian")
	}
}

// ---------------------------------------------------------------------------
// Exponential map
// ---------------------------------------------------------------------------

func TestMatExp(t *testing.T) {
	// exp(0)=I.
	z := ZeroMatrix(3)
	e0, _ := MatExp(z)
	if !e0.ApproxEqual(IdentityMatrix(3), tol) {
		t.Fatal("exp(0)!=I")
	}
	// Rotation: exp(theta*Lz) rotates by theta.
	_, _, lz := SO3Generators()
	theta := math.Pi / 2
	r, _ := MatExp(lz.Scale(theta))
	want, _ := NewMatrixFromRows([][]float64{{0, -1, 0}, {1, 0, 0}, {0, 0, 1}})
	if !r.ApproxEqual(want, 1e-9) {
		t.Fatalf("rotation exp=%v", r.Data)
	}
	// Diagonal check: exp(diag(1,2)) = diag(e,e^2).
	d := DiagMatrix([]float64{1, 2})
	ed, _ := MatExp(d)
	if math.Abs(ed.At(0, 0)-math.E) > 1e-9 || math.Abs(ed.At(1, 1)-math.E*math.E) > 1e-9 {
		t.Fatalf("exp(diag)=%v", ed.Data)
	}
}

func TestCMatExpUnitary(t *testing.T) {
	// exp(i pi/2 sx) = i sx.
	sx := PauliX()
	arg := sx.Scale(complex(0, math.Pi/2))
	u, _ := CMatExp(arg)
	want := sx.Scale(complex(0, 1))
	if !u.ApproxEqual(want, 1e-9) {
		t.Fatalf("cmat exp=%v", u.Data)
	}
	if !u.IsUnitary(1e-9) {
		t.Fatal("exp of anti-Hermitian should be unitary")
	}
}

func TestBCH(t *testing.T) {
	// For commuting matrices BCH reduces to X+Y.
	x := DiagMatrix([]float64{1, 2})
	y := DiagMatrix([]float64{3, -1})
	b, _ := BCHApprox(x, y)
	sum, _ := x.Add(y)
	if !b.ApproxEqual(sum, tol) {
		t.Fatalf("BCH commuting=%v", b.Data)
	}
	// Ad_{exp X}(Y) = exp(X) Y exp(-X); compare series to direct computation.
	e, f, _ := SL2Generators()
	xs := e.Scale(0.1)
	series, _ := ExpBracketSeries(xs, f, 25)
	ex, _ := MatExp(xs)
	exn, _ := MatExp(xs.Scale(-1))
	tmp, _ := ex.Mul(f)
	direct, _ := tmp.Mul(exn)
	if !series.ApproxEqual(direct, 1e-9) {
		t.Fatalf("Ad series != direct: %v vs %v", series.Data, direct.Data)
	}
}

// ---------------------------------------------------------------------------
// Spin representations
// ---------------------------------------------------------------------------

func TestSpinMatrices(t *testing.T) {
	// Spin-1/2 equals sigma/2.
	jx, jy, jz, err := SpinMatrices(0.5)
	if err != nil {
		t.Fatal(err)
	}
	sx, sy, sz := SU2SpinMatrices()
	if !jx.ApproxEqual(sx, tol) || !jy.ApproxEqual(sy, tol) || !jz.ApproxEqual(sz, tol) {
		t.Fatal("spin-1/2 mismatch")
	}
	// General relation [Jx,Jy]=i Jz for spin 1.
	jx1, jy1, jz1, _ := SpinMatrices(1)
	br, _ := CBracket(jx1, jy1)
	want := jz1.Scale(complex(0, 1))
	if !br.ApproxEqual(want, tol) {
		t.Fatalf("[Jx,Jy]=%v", br.Data)
	}
	// J^2 = j(j+1) I for spin 1 (=2).
	xx, _ := jx1.Mul(jx1)
	yy, _ := jy1.Mul(jy1)
	zz, _ := jz1.Mul(jz1)
	s1, _ := xx.Add(yy)
	j2, _ := s1.Add(zz)
	if !j2.ApproxEqual(IdentityCMatrix(3).Scale(2), tol) {
		t.Fatalf("J^2=%v", j2.Data)
	}
}

func TestGellMann(t *testing.T) {
	gm := GellMannMatrices()
	if len(gm) != 8 {
		t.Fatal("expected 8 Gell-Mann matrices")
	}
	for i, g := range gm {
		if !g.IsHermitian(tol) {
			t.Fatalf("lambda%d not Hermitian", i+1)
		}
		if !g.IsTraceless(tol) {
			t.Fatalf("lambda%d not traceless", i+1)
		}
	}
	// tr(lambda_i lambda_j) = 2 delta_ij.
	for i := 0; i < 8; i++ {
		for j := 0; j < 8; j++ {
			tr, _ := CTraceForm(gm[i], gm[j])
			want := complex(0, 0)
			if i == j {
				want = 2
			}
			if cmplx.Abs(tr-want) > tol {
				t.Fatalf("tr(l%d l%d)=%v", i+1, j+1, tr)
			}
		}
	}
}

// ---------------------------------------------------------------------------
// Root systems
// ---------------------------------------------------------------------------

func TestCartanMatricesAndCounts(t *testing.T) {
	tests := []struct {
		family      string
		rank        int
		numPos      int
		dim         int
		coxeter     int
		dualCoxeter int
		weyl        string
	}{
		{"A", 1, 1, 3, 2, 2, "2"},
		{"A", 2, 3, 8, 3, 3, "6"},
		{"A", 3, 6, 15, 4, 4, "24"},
		{"A", 4, 10, 24, 5, 5, "120"},
		{"B", 2, 4, 10, 4, 3, "8"},
		{"B", 3, 9, 21, 6, 5, "48"},
		{"C", 3, 9, 21, 6, 4, "48"},
		{"D", 4, 12, 28, 6, 6, "192"},
		{"D", 5, 20, 45, 8, 8, "1920"},
		{"G", 2, 6, 14, 6, 4, "12"},
		{"F", 4, 24, 52, 12, 9, "1152"},
		{"E", 6, 36, 78, 12, 12, "51840"},
		{"E", 7, 63, 133, 18, 18, "2903040"},
		{"E", 8, 120, 248, 30, 30, "696729600"},
	}
	for _, tc := range tests {
		name := fmt.Sprintf("%s%d", tc.family, tc.rank)
		np, err := NumPositiveRoots(tc.family, tc.rank)
		if err != nil {
			t.Fatalf("%s: %v", name, err)
		}
		if np != tc.numPos {
			t.Errorf("%s numPos=%d want %d", name, np, tc.numPos)
		}
		dim, _ := LieAlgebraDimension(tc.family, tc.rank)
		if dim != tc.dim {
			t.Errorf("%s dim=%d want %d", name, dim, tc.dim)
		}
		h, _ := CoxeterNumber(tc.family, tc.rank)
		if h != tc.coxeter {
			t.Errorf("%s coxeter=%d want %d", name, h, tc.coxeter)
		}
		hv, _ := DualCoxeterNumber(tc.family, tc.rank)
		if hv != tc.dualCoxeter {
			t.Errorf("%s dualCoxeter=%d want %d", name, hv, tc.dualCoxeter)
		}
		hvr, _ := DualCoxeterFromRoots(tc.family, tc.rank)
		if hvr != tc.dualCoxeter {
			t.Errorf("%s dualCoxeterFromRoots=%d want %d", name, hvr, tc.dualCoxeter)
		}
		w, _ := WeylGroupOrder(tc.family, tc.rank)
		if w.String() != tc.weyl {
			t.Errorf("%s |W|=%s want %s", name, w.String(), tc.weyl)
		}
		// number of roots = h * rank.
		if 2*np != h*tc.rank {
			t.Errorf("%s: 2*numPos=%d but h*rank=%d", name, 2*np, h*tc.rank)
		}
	}
}

func TestCartanMatrixKnown(t *testing.T) {
	a2, _ := CartanMatrix("A", 2)
	if !a2.Equal(mustM(t, [][]float64{{2, -1}, {-1, 2}})) {
		t.Fatalf("A2 Cartan=%v", a2.Data)
	}
	g2, _ := CartanMatrix("G", 2)
	if !g2.Equal(mustM(t, [][]float64{{2, -1}, {-3, 2}})) {
		t.Fatalf("G2 Cartan=%v", g2.Data)
	}
	// Determinant of Cartan matrix = order of center.
	det, _ := CartanMatrixDeterminant("A", 3)
	if math.Abs(det-4) > tol {
		t.Fatalf("det A3=%v", det)
	}
	det8, _ := CartanMatrixDeterminant("E", 8)
	if math.Abs(det8-1) > 1e-6 {
		t.Fatalf("det E8=%v", det8)
	}
}

func TestCartanFromCoordinateRoots(t *testing.T) {
	for _, tc := range []struct {
		f string
		r int
	}{{"A", 3}, {"B", 4}, {"C", 4}, {"D", 5}} {
		sr, err := SimpleRoots(tc.f, tc.r)
		if err != nil {
			t.Fatal(err)
		}
		fromRoots, err := CartanMatrixFromRoots(sr)
		if err != nil {
			t.Fatal(err)
		}
		cm, _ := CartanMatrix(tc.f, tc.r)
		if !fromRoots.ApproxEqual(cm, tol) {
			t.Fatalf("%s%d: coordinate Cartan mismatch %v vs %v", tc.f, tc.r, fromRoots.Data, cm.Data)
		}
		// Coordinate positive-root count matches the label count.
		pr, _ := PositiveRoots(tc.f, tc.r)
		np, _ := NumPositiveRoots(tc.f, tc.r)
		if len(pr) != np {
			t.Fatalf("%s%d: coord posroots=%d labels=%d", tc.f, tc.r, len(pr), np)
		}
	}
}

func TestWeylReflectionAndOrbit(t *testing.T) {
	// Reflect e1 across root e1-e2 -> e2.
	root := []float64{1, -1, 0}
	v := []float64{1, 0, 0}
	r, _ := WeylReflection(root, v)
	if !VecEqual(r, []float64{0, 1, 0}, tol) {
		t.Fatalf("reflection=%v", r)
	}
	// Reflection matrix is orthogonal and an involution.
	m, _ := WeylReflectionMatrix(root)
	mm, _ := m.Mul(m)
	if !mm.ApproxEqual(IdentityMatrix(3), tol) {
		t.Fatal("reflection not involutive")
	}
	// Weyl orbit of e1-e2 under A2 simple reflections has 6 elements (all roots).
	sr, _ := SimpleRoots("A", 2)
	orbit, _ := WeylOrbit(sr, sr[0], tol)
	if len(orbit) != 6 {
		t.Fatalf("A2 orbit size=%d want 6", len(orbit))
	}
}

func TestFundamentalWeightsAndRho(t *testing.T) {
	// For A2, <omega_i, alpha_j^vee> = delta_ij.
	fw, _ := FundamentalWeightsRootBasis("A", 2)
	cm, _ := CartanMatrix("A", 2)
	n := 2
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			// <omega_i, alpha_j^vee> = sum_k fw[i][k] * A[k][j].
			s := 0.0
			for k := 0; k < n; k++ {
				s += fw[i][k] * cm.At(k, j)
			}
			want := 0.0
			if i == j {
				want = 1
			}
			if math.Abs(s-want) > tol {
				t.Fatalf("pairing[%d][%d]=%v", i, j, s)
			}
		}
	}
	// rho = half sum of positive roots (check via inner products with coroots =1).
	rho, _ := WeylVectorRootBasis("A", 2)
	for j := 0; j < n; j++ {
		s := 0.0
		for k := 0; k < n; k++ {
			s += rho[k] * cm.At(k, j)
		}
		if math.Abs(s-1) > tol {
			t.Fatalf("<rho,alpha_%d^vee>=%v want 1", j, s)
		}
	}
}

func TestWeylDimension(t *testing.T) {
	tests := []struct {
		family string
		rank   int
		weight []int
		dim    int
	}{
		{"A", 1, []int{1}, 2},    // spin 1/2
		{"A", 1, []int{2}, 3},    // adjoint
		{"A", 2, []int{1, 0}, 3}, // fundamental of su(3)
		{"A", 2, []int{0, 1}, 3}, // anti-fundamental
		{"A", 2, []int{1, 1}, 8}, // adjoint
		{"A", 3, []int{1, 0, 0}, 4},
		{"B", 2, []int{1, 0}, 5}, // vector of so(5)
		{"B", 2, []int{0, 1}, 4}, // spinor of so(5)
		{"G", 2, []int{1, 0}, 7},
		{"G", 2, []int{0, 1}, 14},
		{"D", 4, []int{1, 0, 0, 0}, 8},
		{"E", 8, []int{0, 0, 0, 0, 0, 0, 0, 1}, 248},
	}
	for _, tc := range tests {
		d, err := WeylDimension(tc.family, tc.rank, tc.weight)
		if err != nil {
			t.Fatalf("%s%d %v: %v", tc.family, tc.rank, tc.weight, err)
		}
		if d != tc.dim {
			t.Errorf("%s%d %v dim=%d want %d", tc.family, tc.rank, tc.weight, d, tc.dim)
		}
	}
	// Adjoint dimension equals algebra dimension.
	adjDim, _ := WeylDimension("A", 2, []int{1, 1})
	algDim, _ := LieAlgebraDimension("A", 2)
	if adjDim != algDim {
		t.Fatalf("adjoint dim %d != alg dim %d", adjDim, algDim)
	}
}

func TestCasimir(t *testing.T) {
	// su(2): (lambda,lambda+2rho) with long-root^2=2 gives 2 j(j+1) for a=2j.
	c, _ := CasimirEigenvalue("A", 1, []int{2}) // adjoint, j=1 -> 2*1*2 = 4
	if math.Abs(c-4) > tol {
		t.Fatalf("Casimir A1 adjoint=%v", c)
	}
	// Physics-convention helper.
	if math.Abs(CasimirSU2(0.5)-0.75) > tol {
		t.Fatalf("CasimirSU2(1/2)=%v", CasimirSU2(0.5))
	}
	if d, _ := DimensionSU2(1.5); d != 4 {
		t.Fatalf("DimensionSU2(3/2)=%d", d)
	}
	// A2 adjoint Casimir in these units = 2*h_dual*... just check positivity & symmetry.
	c1, _ := CasimirEigenvalue("A", 2, []int{1, 0})
	c2, _ := CasimirEigenvalue("A", 2, []int{0, 1})
	if math.Abs(c1-c2) > tol {
		t.Fatalf("fundamental/anti-fundamental Casimir differ: %v %v", c1, c2)
	}
}

func TestDimensionHelpers(t *testing.T) {
	if d, _ := DimensionSLn(3); d != 8 {
		t.Fatalf("sl3 dim=%d", d)
	}
	if d, _ := DimensionSOn(5); d != 10 {
		t.Fatalf("so5 dim=%d", d)
	}
	if d, _ := DimensionSPn(3); d != 21 {
		t.Fatalf("sp6 dim=%d", d)
	}
	name, _ := LieAlgebraName("B", 3)
	if name != "so(7)" {
		t.Fatalf("name=%s", name)
	}
}

func TestDynkinDiagram(t *testing.T) {
	adj, _ := DynkinDiagramAdjacency("A", 3)
	want, _ := NewMatrixFromRows([][]float64{{0, 1, 0}, {1, 0, 1}, {0, 1, 0}})
	if !adj.Equal(want) {
		t.Fatalf("A3 adjacency=%v", adj.Data)
	}
	sl, _ := IsSimplyLaced("D", 4)
	if !sl {
		t.Fatal("D4 should be simply laced")
	}
	sl2, _ := IsSimplyLaced("B", 3)
	if sl2 {
		t.Fatal("B3 should not be simply laced")
	}
	bond, _ := DynkinBondMatrix("G", 2)
	if bond.At(0, 1) != 3 || bond.At(1, 0) != 3 {
		t.Fatalf("G2 bond=%v", bond.Data)
	}
}

func TestErrors(t *testing.T) {
	if _, err := CartanMatrix("Z", 3); err != ErrType {
		t.Fatalf("expected ErrType, got %v", err)
	}
	if _, err := CartanMatrix("E", 5); err != ErrRange {
		t.Fatalf("expected ErrRange, got %v", err)
	}
	a, _ := NewMatrixFromRows([][]float64{{1, 2}})
	b, _ := NewMatrixFromRows([][]float64{{1, 2}})
	if _, err := a.Mul(b); err != ErrDim {
		t.Fatalf("expected ErrDim, got %v", err)
	}
	sing := ZeroMatrix(2)
	if _, err := Inverse(sing); err != ErrSingular {
		t.Fatalf("expected ErrSingular, got %v", err)
	}
}

func TestSO3VectorRoundTrip(t *testing.T) {
	w := []float64{0.3, -1.2, 2.5}
	m, _ := SO3FromVector(w)
	if !m.IsAntisymmetric(tol) {
		t.Fatal("hat map not antisymmetric")
	}
	back, _ := SO3ToVector(m)
	if !VecEqual(back, w, tol) {
		t.Fatalf("round trip=%v", back)
	}
	// [w]_x u = w x u.
	u := []float64{1, 0, 0}
	mu, _ := m.MatVec(u)
	// w x u = (w2*u3-w3*u2, w3*u1-w1*u3, w1*u2-w2*u1)
	cross := []float64{w[1]*u[2] - w[2]*u[1], w[2]*u[0] - w[0]*u[2], w[0]*u[1] - w[1]*u[0]}
	if !VecEqual(mu, cross, tol) {
		t.Fatalf("cross product mismatch %v vs %v", mu, cross)
	}
}

// ---------------------------------------------------------------------------
// Helpers and example
// ---------------------------------------------------------------------------

func mustM(t *testing.T, rows [][]float64) *Matrix {
	t.Helper()
	m, err := NewMatrixFromRows(rows)
	if err != nil {
		t.Fatal(err)
	}
	return m
}

// ExampleWeylDimension prints the dimension of the adjoint representation of
// su(3), whose highest weight has Dynkin labels (1,1).
func ExampleWeylDimension() {
	dim, _ := WeylDimension("A", 2, []int{1, 1})
	fmt.Println(dim)
	// Output: 8
}

// ExampleBracket demonstrates the sl(2) commutation relation [E,F]=H.
func ExampleBracket() {
	e, f, h := SL2Generators()
	br, _ := Bracket(e, f)
	fmt.Println(br.Equal(h))
	// Output: true
}

// ExampleWeylGroupOrder prints the order of the Weyl group of E8.
func ExampleWeylGroupOrder() {
	w, _ := WeylGroupOrder("E", 8)
	fmt.Println(w.String())
	// Output: 696729600
}
