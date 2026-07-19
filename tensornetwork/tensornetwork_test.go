package tensornetwork

import (
	"fmt"
	"math"
	"testing"
)

const tol = 1e-9

func approx(a, b, t float64) bool { return math.Abs(a-b) <= t }

func mustTensor(t *testing.T, data []float64, shape ...int) *Tensor {
	t.Helper()
	tn, err := NewWithData(data, shape...)
	if err != nil {
		t.Fatalf("NewWithData: %v", err)
	}
	return tn
}

func mustMatrix(t *testing.T, r, c int, data []float64) *Matrix {
	t.Helper()
	m, err := NewMatrixData(r, c, data)
	if err != nil {
		t.Fatalf("NewMatrixData: %v", err)
	}
	return m
}

func TestTensorBasics(t *testing.T) {
	tn := mustTensor(t, []float64{1, 2, 3, 4, 5, 6, 7, 8}, 2, 2, 2)
	if tn.Rank() != 3 {
		t.Errorf("rank = %d, want 3", tn.Rank())
	}
	if tn.Size() != 8 {
		t.Errorf("size = %d, want 8", tn.Size())
	}
	if got := tn.At(1, 1, 0); got != 7 {
		t.Errorf("At(1,1,0) = %v, want 7", got)
	}
	tn.Set(99, 0, 0, 0)
	if got := tn.At(0, 0, 0); got != 99 {
		t.Errorf("Set/At = %v, want 99", got)
	}
	mi := tn.MultiIndex(5)
	if mi[0] != 1 || mi[1] != 0 || mi[2] != 1 {
		t.Errorf("MultiIndex(5) = %v, want [1 0 1]", mi)
	}
}

func TestReshapePermute(t *testing.T) {
	tn := mustTensor(t, []float64{1, 2, 3, 4, 5, 6}, 2, 3)
	r, err := tn.Reshape(3, 2)
	if err != nil {
		t.Fatal(err)
	}
	if r.At(2, 1) != 6 {
		t.Errorf("reshape At(2,1) = %v, want 6", r.At(2, 1))
	}
	r2, err := tn.Reshape(-1)
	if err != nil || r2.Rank() != 1 || r2.Size() != 6 {
		t.Errorf("reshape -1 failed: %v shape=%v", err, r2.Shape())
	}
	tp := tn.Transpose()
	want := []float64{1, 4, 2, 5, 3, 6}
	if !tp.Equal(mustTensor(t, want, 3, 2), tol) {
		t.Errorf("transpose = %v, want %v", tp.Data(), want)
	}
	p, err := tn.Permute(1, 0)
	if err != nil || !p.Equal(tp, tol) {
		t.Errorf("permute(1,0) != transpose")
	}
}

func TestElementwiseReductions(t *testing.T) {
	a := mustTensor(t, []float64{1, 2, 3, 4}, 2, 2)
	b := mustTensor(t, []float64{5, 6, 7, 8}, 2, 2)
	sum, _ := a.Add(b)
	if !sum.Equal(mustTensor(t, []float64{6, 8, 10, 12}, 2, 2), tol) {
		t.Errorf("Add wrong: %v", sum.Data())
	}
	if got := a.Sum(); got != 10 {
		t.Errorf("Sum = %v want 10", got)
	}
	if got := a.Norm(); !approx(got, math.Sqrt(30), tol) {
		t.Errorf("Norm = %v want sqrt(30)", got)
	}
	if got := a.Max(); got != 4 {
		t.Errorf("Max = %v want 4", got)
	}
	sa, _ := a.SumAxis(0)
	if !sa.Equal(FromVector([]float64{4, 6}), tol) {
		t.Errorf("SumAxis(0) = %v want [4 6]", sa.Data())
	}
	sa1, _ := a.SumAxis(1)
	if !sa1.Equal(FromVector([]float64{3, 7}), tol) {
		t.Errorf("SumAxis(1) = %v want [3 7]", sa1.Data())
	}
	if got := a.NormP(1); got != 10 {
		t.Errorf("NormP(1) = %v want 10", got)
	}
}

func TestUnfoldFold(t *testing.T) {
	tn := mustTensor(t, []float64{1, 2, 3, 4, 5, 6, 7, 8}, 2, 2, 2)
	tests := []struct {
		mode int
		want []float64
	}{
		{0, []float64{1, 2, 3, 4, 5, 6, 7, 8}},
		{1, []float64{1, 2, 5, 6, 3, 4, 7, 8}},
		{2, []float64{1, 3, 5, 7, 2, 4, 6, 8}},
	}
	for _, tc := range tests {
		m, err := tn.Unfold(tc.mode)
		if err != nil {
			t.Fatal(err)
		}
		for i, w := range tc.want {
			if m.data[i] != w {
				t.Errorf("mode %d unfold[%d] = %v want %v", tc.mode, i, m.data[i], w)
			}
		}
		back, err := Fold(m, tc.mode, tn.Shape())
		if err != nil {
			t.Fatal(err)
		}
		if !back.Equal(tn, tol) {
			t.Errorf("mode %d fold roundtrip failed: %v", tc.mode, back.Data())
		}
	}
}

func TestKroneckerKhatriRao(t *testing.T) {
	a := mustMatrix(t, 2, 2, []float64{1, 2, 3, 4})
	b := mustMatrix(t, 2, 2, []float64{0, 1, 1, 0})
	k := KroneckerMatrix(a, b)
	want := []float64{
		0, 1, 0, 2,
		1, 0, 2, 0,
		0, 3, 0, 4,
		3, 0, 4, 0,
	}
	for i, w := range want {
		if k.data[i] != w {
			t.Errorf("kron[%d] = %v want %v", i, k.data[i], w)
		}
	}
	// Khatri-Rao of two 2x2 matrices.
	c := mustMatrix(t, 2, 2, []float64{1, 2, 3, 4})
	d := mustMatrix(t, 2, 2, []float64{5, 6, 7, 8})
	kr, err := KhatriRaoMatrix(c, d)
	if err != nil {
		t.Fatal(err)
	}
	// column 0: c[:,0]=[1,3] kron d[:,0]=[5,7] -> [5,7,15,21]
	wantKR := []float64{
		5, 12,
		7, 16,
		15, 24,
		21, 32,
	}
	for i, w := range wantKR {
		if !approx(kr.data[i], w, tol) {
			t.Errorf("khatrirao[%d] = %v want %v", i, kr.data[i], w)
		}
	}
}

func TestModeProduct(t *testing.T) {
	tn := mustTensor(t, []float64{1, 2, 3, 4, 5, 6, 7, 8}, 2, 2, 2)
	id := IdentityMatrix(2)
	mp, err := ModeProduct(tn, id, 0)
	if err != nil {
		t.Fatal(err)
	}
	if !mp.Equal(tn, tol) {
		t.Errorf("mode product with identity changed tensor")
	}
	// Doubling matrix on mode 1.
	dbl := mustMatrix(t, 2, 2, []float64{2, 0, 0, 2})
	mp2, _ := ModeProduct(tn, dbl, 1)
	if !mp2.Equal(tn.Scale(2), tol) {
		t.Errorf("mode product doubling failed: %v", mp2.Data())
	}
}

func TestTensorDotEinsumMatmul(t *testing.T) {
	a := mustTensor(t, []float64{1, 2, 3, 4, 5, 6}, 2, 3)
	b := mustTensor(t, []float64{1, 0, 0, 1, 1, 1}, 3, 2)
	// TensorDot contracting a axis1 with b axis0 == matmul.
	td, err := TensorDot(a, b, []int{1}, []int{0})
	if err != nil {
		t.Fatal(err)
	}
	want := mustTensor(t, []float64{4, 5, 10, 11}, 2, 2)
	if !td.Equal(want, tol) {
		t.Errorf("TensorDot matmul = %v want %v", td.Data(), want.Data())
	}
	es, err := Einsum("ij,jk->ik", a, b)
	if err != nil {
		t.Fatal(err)
	}
	if !es.Equal(want, tol) {
		t.Errorf("Einsum matmul = %v want %v", es.Data(), want.Data())
	}
}

func TestEinsumTraceTranspose(t *testing.T) {
	m := mustTensor(t, []float64{1, 2, 3, 4}, 2, 2)
	tr, err := Einsum("ii->", m)
	if err != nil {
		t.Fatal(err)
	}
	if tr.Rank() != 0 || tr.AtFlat(0) != 5 {
		t.Errorf("Einsum trace = %v want 5", tr.Data())
	}
	tp, err := Einsum("ij->ji", m)
	if err != nil {
		t.Fatal(err)
	}
	if !tp.Equal(mustTensor(t, []float64{1, 3, 2, 4}, 2, 2), tol) {
		t.Errorf("Einsum transpose = %v", tp.Data())
	}
	// Outer via einsum.
	v := FromVector([]float64{1, 2})
	w := FromVector([]float64{3, 4})
	o, err := Einsum("i,j->ij", v, w)
	if err != nil {
		t.Fatal(err)
	}
	if !o.Equal(mustTensor(t, []float64{3, 4, 6, 8}, 2, 2), tol) {
		t.Errorf("Einsum outer = %v", o.Data())
	}
}

func TestSVDReconstruct(t *testing.T) {
	m := mustMatrix(t, 2, 2, []float64{1, 2, 3, 4})
	U, s, V, err := m.SVD()
	if err != nil {
		t.Fatal(err)
	}
	if !approx(s[0], 5.4649857, 1e-6) || !approx(s[1], 0.36596619, 1e-6) {
		t.Errorf("singular values = %v want [5.4649857 0.36596619]", s)
	}
	// Reconstruct U diag(s) V^T.
	rec := reconstructSVD(U, s, V)
	if !rec.Equal(m, 1e-9) {
		t.Errorf("SVD reconstruct = %v want %v", rec.data, m.data)
	}
	// Non-square tall matrix.
	m2 := mustMatrix(t, 3, 2, []float64{1, 2, 3, 4, 5, 6})
	U2, s2, V2, _ := m2.SVD()
	rec2 := reconstructSVD(U2, s2, V2)
	if !rec2.Equal(m2, 1e-8) {
		t.Errorf("SVD tall reconstruct failed")
	}
}

func reconstructSVD(U *Matrix, s []float64, V *Matrix) *Matrix {
	sv := DiagMatrix(s)
	us := MatMul(U, sv)
	return MatMul(us, V.Transpose())
}

func TestQR(t *testing.T) {
	m := mustMatrix(t, 3, 2, []float64{1, 2, 3, 4, 5, 6})
	Q, R, err := m.QR()
	if err != nil {
		t.Fatal(err)
	}
	rec := MatMul(Q, R)
	if !rec.Equal(m, 1e-9) {
		t.Errorf("QR reconstruct failed: %v", rec.data)
	}
	// Q columns orthonormal: Q^T Q = I.
	qtq := MatMul(Q.Transpose(), Q)
	if !qtq.Equal(IdentityMatrix(2), 1e-9) {
		t.Errorf("Q not orthonormal: %v", qtq.data)
	}
}

func TestEigSym(t *testing.T) {
	m := mustMatrix(t, 2, 2, []float64{2, 1, 1, 2})
	vals, vecs, err := m.EigSym()
	if err != nil {
		t.Fatal(err)
	}
	if !approx(vals[0], 3, 1e-9) || !approx(vals[1], 1, 1e-9) {
		t.Errorf("eigenvalues = %v want [3 1]", vals)
	}
	// A V = V diag(vals).
	av := MatMul(m, vecs)
	vd := MatMul(vecs, DiagMatrix(vals))
	if !av.Equal(vd, 1e-8) {
		t.Errorf("eigenvector relation failed")
	}
}

func TestPinvSolve(t *testing.T) {
	m := mustMatrix(t, 2, 2, []float64{4, 0, 0, 2})
	p := m.Pinv(0)
	if !approx(p.At(0, 0), 0.25, tol) || !approx(p.At(1, 1), 0.5, tol) {
		t.Errorf("pinv diag = %v", p.data)
	}
	x, err := m.Solve([]float64{8, 6})
	if err != nil {
		t.Fatal(err)
	}
	if !approx(x[0], 2, 1e-9) || !approx(x[1], 3, 1e-9) {
		t.Errorf("Solve = %v want [2 3]", x)
	}
}

func TestMatrixRankCond(t *testing.T) {
	// Rank-1 matrix.
	m := mustMatrix(t, 2, 2, []float64{1, 2, 2, 4})
	if r := m.Rank(1e-9); r != 1 {
		t.Errorf("rank = %d want 1", r)
	}
	id := IdentityMatrix(3)
	if c := id.Cond(); !approx(c, 1, 1e-9) {
		t.Errorf("cond(I) = %v want 1", c)
	}
}

func TestHOSVDExact(t *testing.T) {
	tn := RandTensor(7, 3, 4, 2)
	full := HOSVDMust(t, tn, []int{3, 4, 2})
	re, _ := full.RelError(tn)
	if re > 1e-9 {
		t.Errorf("full HOSVD relerror = %v, want ~0", re)
	}
	// Core dimensions match requested ranks.
	if !shapeEqual(full.Core.Shape(), []int{3, 4, 2}) {
		t.Errorf("core shape = %v", full.Core.Shape())
	}
}

func HOSVDMust(t *testing.T, tn *Tensor, ranks []int) *TuckerDecomposition {
	t.Helper()
	d, err := HOSVD(tn, ranks)
	if err != nil {
		t.Fatal(err)
	}
	return d
}

func TestHOSVDTruncateLowRank(t *testing.T) {
	// Build an exactly multilinear-rank-(1,1,1) tensor: outer product of vectors.
	a := FromVector([]float64{1, 2, 3})
	b := FromVector([]float64{1, -1})
	c := FromVector([]float64{2, 1})
	ab := Outer(a, b)
	tn := Outer(ab, c) // shape (3,2,2), multilinear rank (1,1,1)
	d, err := HOSVD(tn, []int{1, 1, 1})
	if err != nil {
		t.Fatal(err)
	}
	re, _ := d.RelError(tn)
	if re > 1e-9 {
		t.Errorf("rank-1 HOSVD relerror = %v want ~0", re)
	}
}

func TestHOOIImproves(t *testing.T) {
	a := FromVector([]float64{1, 2, 3})
	b := FromVector([]float64{1, -1, 2})
	c := FromVector([]float64{2, 1})
	tn := Outer(Outer(a, b), c)
	d, err := HOOI(tn, []int{1, 1, 1}, 20)
	if err != nil {
		t.Fatal(err)
	}
	re, _ := d.RelError(tn)
	if re > 1e-8 {
		t.Errorf("HOOI relerror = %v want ~0", re)
	}
}

func TestCPALSLowRank(t *testing.T) {
	// Construct a rank-2 tensor from known factors and recover it.
	A := mustMatrix(t, 3, 2, []float64{1, 0.5, 0, 1, 2, -1})
	B := mustMatrix(t, 3, 2, []float64{1, 1, 0.5, 0, -1, 2})
	C := mustMatrix(t, 2, 2, []float64{1, -1, 2, 1})
	cp0 := &CPDecomposition{Factors: []*Matrix{A, B, C}, Weights: []float64{1, 1}, Rank: 2, Shape: []int{3, 3, 2}}
	tn := cp0.Reconstruct()
	opts := DefaultCPOptions()
	opts.MaxIter = 1000
	opts.Seed = 3
	cp, err := CPALS(tn, 2, opts)
	if err != nil {
		t.Fatal(err)
	}
	re, _ := cp.RelError(tn)
	if re > 1e-4 {
		t.Errorf("CP-ALS relerror = %v want < 1e-4", re)
	}
}

func TestTTSVDExact(t *testing.T) {
	tn := RandTensor(11, 2, 3, 2, 2)
	tt, err := TTSVD(tn, 0)
	if err != nil {
		t.Fatal(err)
	}
	re, _ := tt.RelError(tn)
	if re > 1e-9 {
		t.Errorf("TTSVD exact relerror = %v want ~0", re)
	}
	if tt.Ranks[0] != 1 || tt.Ranks[len(tt.Ranks)-1] != 1 {
		t.Errorf("boundary ranks = %v", tt.Ranks)
	}
}

func TestTTSVDRankCap(t *testing.T) {
	tn := RandTensor(5, 3, 3, 3)
	tt, err := TTSVDRank(tn, []int{2, 2})
	if err != nil {
		t.Fatal(err)
	}
	if tt.MaxRank() > 2 {
		t.Errorf("max rank = %d want <= 2", tt.MaxRank())
	}
	if tt.Compression() <= 0 {
		t.Errorf("compression = %v", tt.Compression())
	}
}

func TestTTReconstructMatchesReshape(t *testing.T) {
	// Exact TT of a small tensor reconstructs entrywise.
	tn := mustTensor(t, []float64{1, 2, 3, 4, 5, 6, 7, 8}, 2, 2, 2)
	tt, err := TTSVD(tn, 0)
	if err != nil {
		t.Fatal(err)
	}
	rec := tt.Reconstruct()
	if !rec.Equal(tn, 1e-9) {
		t.Errorf("TT reconstruct = %v want %v", rec.Data(), tn.Data())
	}
}

func TestNetworkContractMatchesEinsum(t *testing.T) {
	a := mustTensor(t, []float64{1, 2, 3, 4, 5, 6}, 2, 3)
	b := mustTensor(t, []float64{1, 0, 0, 1, 1, 1}, 3, 2)
	nw := NewNetwork()
	if err := nw.Add(a, "i", "j"); err != nil {
		t.Fatal(err)
	}
	if err := nw.Add(b, "j", "k"); err != nil {
		t.Fatal(err)
	}
	res, labels, err := nw.Contract()
	if err != nil {
		t.Fatal(err)
	}
	want, _ := Einsum("ij,jk->ik", a, b)
	if !res.Equal(want, tol) {
		t.Errorf("network contract = %v want %v (labels %v)", res.Data(), want.Data(), labels)
	}
}

func TestNetworkOptimalOrder(t *testing.T) {
	// Chain of three matrices; optimal and greedy should give same result.
	a := RandTensor(1, 2, 5)
	b := RandTensor(2, 5, 4)
	c := RandTensor(3, 4, 3)
	build := func() *Network {
		nw := NewNetwork()
		nw.Add(a, "i", "j")
		nw.Add(b, "j", "k")
		nw.Add(c, "k", "l")
		return nw
	}
	nw := build()
	order, cost, err := nw.OptimalOrder()
	if err != nil {
		t.Fatal(err)
	}
	if cost <= 0 {
		t.Errorf("optimal cost = %v", cost)
	}
	res, _, err := nw.ContractOptimal()
	if err != nil {
		t.Fatal(err)
	}
	// Compare to einsum ijk chain.
	want, _ := Einsum("ij,jk,kl->il", a, b, c)
	if !res.Equal(want, 1e-8) {
		t.Errorf("optimal contract mismatch")
	}
	gc, _ := nw.ContractionCost(order)
	if !approx(gc, cost, 1e-6) {
		t.Errorf("ContractionCost(order) = %v want %v", gc, cost)
	}
}

func TestConcatenateStack(t *testing.T) {
	a := mustTensor(t, []float64{1, 2, 3, 4}, 2, 2)
	b := mustTensor(t, []float64{5, 6, 7, 8}, 2, 2)
	c, err := Concatenate(0, a, b)
	if err != nil {
		t.Fatal(err)
	}
	if !shapeEqual(c.Shape(), []int{4, 2}) || c.At(3, 1) != 8 {
		t.Errorf("concatenate shape=%v", c.Shape())
	}
	s, err := Stack(0, a, b)
	if err != nil {
		t.Fatal(err)
	}
	if !shapeEqual(s.Shape(), []int{2, 2, 2}) || s.At(1, 0, 0) != 5 {
		t.Errorf("stack shape=%v", s.Shape())
	}
}

func TestOuterKronecker(t *testing.T) {
	a := FromVector([]float64{1, 2})
	b := FromVector([]float64{3, 4})
	o := Outer(a, b)
	if !o.Equal(mustTensor(t, []float64{3, 4, 6, 8}, 2, 2), tol) {
		t.Errorf("Outer = %v", o.Data())
	}
	ka := mustTensor(t, []float64{1, 2}, 2)
	kb := mustTensor(t, []float64{10, 20}, 2)
	k, err := Kronecker(ka, kb)
	if err != nil {
		t.Fatal(err)
	}
	if !k.Equal(mustTensor(t, []float64{10, 20, 20, 40}, 4), tol) {
		t.Errorf("Kronecker = %v", k.Data())
	}
}

func TestMultilinearRank(t *testing.T) {
	a := FromVector([]float64{1, 2, 3})
	b := FromVector([]float64{1, -1})
	tn := Outer(Outer(a, b), FromVector([]float64{2, 5}))
	mr := MultilinearRank(tn, 1e-9)
	for i, r := range mr {
		if r != 1 {
			t.Errorf("multilinear rank[%d] = %d want 1", i, r)
		}
	}
}

func ExampleEinsum() {
	a, _ := NewWithData([]float64{1, 2, 3, 4, 5, 6}, 2, 3)
	b, _ := NewWithData([]float64{1, 0, 0, 1, 1, 1}, 3, 2)
	c, _ := Einsum("ij,jk->ik", a, b)
	fmt.Println(c.Data())
	// Output: [4 5 10 11]
}

func ExampleTensor_Unfold() {
	tn, _ := NewWithData([]float64{1, 2, 3, 4, 5, 6, 7, 8}, 2, 2, 2)
	m, _ := tn.Unfold(0)
	fmt.Println(m.Row(0))
	fmt.Println(m.Row(1))
	// Output:
	// [1 2 3 4]
	// [5 6 7 8]
}

func ExampleCPALS() {
	tn, _ := NewWithData([]float64{1, 0, 0, 1, 0, 1, 1, 0}, 2, 2, 2)
	opts := DefaultCPOptions()
	opts.MaxIter = 500
	cp, _ := CPALS(tn, 2, opts)
	re, _ := cp.RelError(tn)
	fmt.Println(re < 1e-3)
	// Output: true
}
