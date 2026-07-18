package matrix

import (
	"errors"
	"math"
	"testing"

	"github.com/malcolmston/algebra"
)

// matrixFloatAt evaluates entry (i,j) of m to a float64, failing the test on
// any evaluation error.
func matrixFloatAt(t *testing.T, m *Matrix, i, j int) float64 {
	t.Helper()
	f, err := algebra.Evalf(m.At(i, j))
	if err != nil {
		t.Fatalf("Evalf(%s): %v", m.At(i, j).String(), err)
	}
	return f
}

// matrixApproxEqual fails the test unless every entry of got is within tol of
// the corresponding entry of want.
func matrixApproxEqual(t *testing.T, got, want *Matrix, tol float64) {
	t.Helper()
	if got.Rows() != want.Rows() || got.Cols() != want.Cols() {
		t.Fatalf("shape mismatch: got %dx%d, want %dx%d",
			got.Rows(), got.Cols(), want.Rows(), want.Cols())
	}
	for i := 0; i < got.Rows(); i++ {
		for j := 0; j < got.Cols(); j++ {
			g := matrixFloatAt(t, got, i, j)
			w := matrixFloatAt(t, want, i, j)
			if math.Abs(g-w) > tol {
				t.Fatalf("entry (%d,%d) = %v, want %v", i, j, g, w)
			}
		}
	}
}

func TestLUReconstruct(t *testing.T) {
	cases := [][][]float64{
		{{4, 3}, {6, 3}},
		{{2, 1, 1}, {1, 3, 2}, {1, 0, 0}},
		{{1, 2, 3}, {4, 5, 6}, {7, 8, 10}},
		{{10, -7, 0}, {-3, 2, 6}, {5, -1, 5}},
	}
	for idx, c := range cases {
		a := FromFloats(c)
		l, u, p, sign, err := a.LU()
		if err != nil {
			t.Fatalf("case %d: LU error: %v", idx, err)
		}
		if sign != 1 && sign != -1 {
			t.Fatalf("case %d: sign = %d, want ±1", idx, sign)
		}
		// P*A == L*U.
		pa, err := p.Mul(a)
		if err != nil {
			t.Fatalf("case %d: P*A: %v", idx, err)
		}
		lu, err := l.Mul(u)
		if err != nil {
			t.Fatalf("case %d: L*U: %v", idx, err)
		}
		matrixApproxEqual(t, pa, lu, 1e-9)

		// L is unit lower triangular; U is upper triangular.
		n := a.Rows()
		for i := 0; i < n; i++ {
			for j := 0; j < n; j++ {
				lij := matrixFloatAt(t, l, i, j)
				uij := matrixFloatAt(t, u, i, j)
				if j > i && math.Abs(lij) > 1e-12 {
					t.Fatalf("case %d: L[%d][%d] = %v, want 0", idx, i, j, lij)
				}
				if i == j && math.Abs(lij-1) > 1e-12 {
					t.Fatalf("case %d: L[%d][%d] = %v, want 1", idx, i, j, lij)
				}
				if j < i && math.Abs(uij) > 1e-12 {
					t.Fatalf("case %d: U[%d][%d] = %v, want 0", idx, i, j, uij)
				}
			}
		}
	}
}

func TestLUNotSquare(t *testing.T) {
	a := FromFloats([][]float64{{1, 2, 3}, {4, 5, 6}})
	if _, _, _, _, err := a.LU(); !errors.Is(err, ErrNotSquare) {
		t.Fatalf("LU on non-square: err = %v, want ErrNotSquare", err)
	}
}

func TestLUSymbolicUnsupported(t *testing.T) {
	a := FromExpr([][]algebra.Expr{
		{algebra.Sym("a"), algebra.Int(1)},
		{algebra.Int(2), algebra.Int(3)},
	})
	if _, _, _, _, err := a.LU(); !errors.Is(err, ErrUnsupported) {
		t.Fatalf("LU on symbolic: err = %v, want ErrUnsupported", err)
	}
}

func TestDetLUKnown(t *testing.T) {
	cases := []struct {
		in   [][]float64
		want float64
	}{
		{[][]float64{{4, 3}, {6, 3}}, -6},
		{[][]float64{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, 0},
		{[][]float64{{2, 1, 1}, {1, 3, 2}, {1, 0, 0}}, -1},
		{[][]float64{{2, 0, 0, 0}, {0, 3, 0, 0}, {0, 0, 4, 0}, {0, 0, 0, 5}}, 120},
		{[][]float64{{6, 1, 1}, {4, -2, 5}, {2, 8, 7}}, -306},
	}
	for idx, c := range cases {
		a := FromFloats(c.in)
		d, err := a.DetLU()
		if err != nil {
			t.Fatalf("case %d: DetLU error: %v", idx, err)
		}
		got, err := algebra.Evalf(d)
		if err != nil {
			t.Fatalf("case %d: Evalf: %v", idx, err)
		}
		if math.Abs(got-c.want) > 1e-9 {
			t.Fatalf("case %d: DetLU = %v, want %v", idx, got, c.want)
		}
	}
}

func TestDetLUMatchesCofactor(t *testing.T) {
	a := FromInts([][]int64{
		{3, 1, 0, 2},
		{1, 4, 2, 1},
		{0, 2, 5, 3},
		{2, 1, 3, 6},
	})
	exact, err := a.Det()
	if err != nil {
		t.Fatal(err)
	}
	fast, err := a.DetLU()
	if err != nil {
		t.Fatal(err)
	}
	ex, _ := algebra.Evalf(exact)
	fa, _ := algebra.Evalf(fast)
	if math.Abs(ex-fa) > 1e-9 {
		t.Fatalf("DetLU = %v, cofactor Det = %v", fa, ex)
	}
}

func TestDetLUNotSquare(t *testing.T) {
	a := FromFloats([][]float64{{1, 2, 3}, {4, 5, 6}})
	if _, err := a.DetLU(); !errors.Is(err, ErrNotSquare) {
		t.Fatalf("DetLU non-square: err = %v, want ErrNotSquare", err)
	}
}

func TestSolveLUKnown(t *testing.T) {
	cases := []struct {
		a    [][]float64
		b    []float64
		want []float64
	}{
		{[][]float64{{2, 1}, {1, -1}}, []float64{5, 1}, []float64{2, 1}},
		{[][]float64{{3, 2, -1}, {2, -2, 4}, {-1, 0.5, -1}}, []float64{1, -2, 0}, []float64{1, -2, -2}},
		{[][]float64{{10, -7, 0}, {-3, 2, 6}, {5, -1, 5}}, []float64{7, 4, 6}, []float64{0, -1, 1}},
	}
	for idx, c := range cases {
		a := FromFloats(c.a)
		b := NewVector(floatsToExpr(c.b)...)
		x, err := a.SolveLU(b)
		if err != nil {
			t.Fatalf("case %d: SolveLU error: %v", idx, err)
		}
		if x.Len() != len(c.want) {
			t.Fatalf("case %d: len = %d, want %d", idx, x.Len(), len(c.want))
		}
		for i := range c.want {
			got, _ := algebra.Evalf(x.At(i))
			if math.Abs(got-c.want[i]) > 1e-9 {
				t.Fatalf("case %d: x[%d] = %v, want %v", idx, i, got, c.want[i])
			}
		}
	}
}

func TestSolveLUSingular(t *testing.T) {
	a := FromFloats([][]float64{{1, 2}, {2, 4}})
	b := NewVector(algebra.Flt(1), algebra.Flt(2))
	if _, err := a.SolveLU(b); !errors.Is(err, ErrSingular) {
		t.Fatalf("SolveLU singular: err = %v, want ErrSingular", err)
	}
}

func TestSolveLUDimension(t *testing.T) {
	a := FromFloats([][]float64{{1, 0}, {0, 1}})
	b := NewVector(algebra.Flt(1), algebra.Flt(2), algebra.Flt(3))
	if _, err := a.SolveLU(b); !errors.Is(err, ErrDimension) {
		t.Fatalf("SolveLU bad dim: err = %v, want ErrDimension", err)
	}
}

func TestQRReconstruct(t *testing.T) {
	cases := [][][]float64{
		{{12, -51, 4}, {6, 167, -68}, {-4, 24, -41}},
		{{1, 2}, {3, 4}, {5, 6}},
		{{2, 0}, {0, 3}},
		{{1, 1, 0}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}},
	}
	for idx, c := range cases {
		a := FromFloats(c)
		q, r, err := a.QR()
		if err != nil {
			t.Fatalf("case %d: QR error: %v", idx, err)
		}
		mRows, nCols := a.Rows(), a.Cols()
		if q.Rows() != mRows || q.Cols() != mRows {
			t.Fatalf("case %d: Q is %dx%d, want %dx%d", idx, q.Rows(), q.Cols(), mRows, mRows)
		}
		if r.Rows() != mRows || r.Cols() != nCols {
			t.Fatalf("case %d: R is %dx%d, want %dx%d", idx, r.Rows(), r.Cols(), mRows, nCols)
		}
		// A == Q*R.
		qr, err := q.Mul(r)
		if err != nil {
			t.Fatalf("case %d: Q*R: %v", idx, err)
		}
		matrixApproxEqual(t, qr, a, 1e-9)

		// Q orthogonal: Qᵀ*Q == I.
		qtq, err := q.Transpose().Mul(q)
		if err != nil {
			t.Fatalf("case %d: Qᵀ*Q: %v", idx, err)
		}
		matrixApproxEqual(t, qtq, Identity(mRows), 1e-9)

		// R upper triangular.
		for i := 0; i < mRows; i++ {
			for j := 0; j < nCols && j < i; j++ {
				if v := matrixFloatAt(t, r, i, j); math.Abs(v) > 1e-9 {
					t.Fatalf("case %d: R[%d][%d] = %v, want 0", idx, i, j, v)
				}
			}
		}
	}
}

func TestQRWideUnsupported(t *testing.T) {
	a := FromFloats([][]float64{{1, 2, 3}})
	if _, _, err := a.QR(); !errors.Is(err, ErrDimension) {
		t.Fatalf("QR wide: err = %v, want ErrDimension", err)
	}
}

func TestCholeskyReconstruct(t *testing.T) {
	cases := [][][]float64{
		{{4, 2}, {2, 3}},
		{{25, 15, -5}, {15, 18, 0}, {-5, 0, 11}},
		{{2, -1, 0}, {-1, 2, -1}, {0, -1, 2}},
	}
	for idx, c := range cases {
		a := FromFloats(c)
		l, err := a.Cholesky()
		if err != nil {
			t.Fatalf("case %d: Cholesky error: %v", idx, err)
		}
		// L*Lᵀ == A.
		llt, err := l.Mul(l.Transpose())
		if err != nil {
			t.Fatalf("case %d: L*Lᵀ: %v", idx, err)
		}
		matrixApproxEqual(t, llt, a, 1e-9)
		// L lower triangular.
		n := a.Rows()
		for i := 0; i < n; i++ {
			for j := i + 1; j < n; j++ {
				if v := matrixFloatAt(t, l, i, j); math.Abs(v) > 1e-12 {
					t.Fatalf("case %d: L[%d][%d] = %v, want 0", idx, i, j, v)
				}
			}
		}
	}
}

func TestCholeskyNotPositiveDefinite(t *testing.T) {
	// Symmetric but indefinite.
	a := FromFloats([][]float64{{1, 2}, {2, 1}})
	if _, err := a.Cholesky(); !errors.Is(err, ErrNotPositiveDefinite) {
		t.Fatalf("indefinite: err = %v, want ErrNotPositiveDefinite", err)
	}
	// Asymmetric.
	b := FromFloats([][]float64{{4, 2}, {1, 3}})
	if _, err := b.Cholesky(); !errors.Is(err, ErrNotPositiveDefinite) {
		t.Fatalf("asymmetric: err = %v, want ErrNotPositiveDefinite", err)
	}
}

func TestCholeskyNotSquare(t *testing.T) {
	a := FromFloats([][]float64{{1, 0, 0}, {0, 1, 0}})
	if _, err := a.Cholesky(); !errors.Is(err, ErrNotSquare) {
		t.Fatalf("Cholesky non-square: err = %v, want ErrNotSquare", err)
	}
}

// floatsToExpr converts a []float64 to a []algebra.Expr of Flt literals.
func floatsToExpr(vs []float64) []algebra.Expr {
	out := make([]algebra.Expr, len(vs))
	for i, v := range vs {
		out[i] = algebra.Flt(v)
	}
	return out
}

// matrixBenchMatrix builds a deterministic, diagonally dominant n×n matrix that
// is nonsingular (and symmetric positive definite) for benchmarking.
func matrixBenchMatrix(n int) *Matrix {
	vals := make([][]float64, n)
	for i := 0; i < n; i++ {
		vals[i] = make([]float64, n)
		for j := 0; j < n; j++ {
			vals[i][j] = 1.0 / float64(i+j+1)
		}
		vals[i][i] += float64(n)
	}
	return FromFloats(vals)
}

func BenchmarkDetLU(b *testing.B) {
	m := matrixBenchMatrix(16)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if _, err := m.DetLU(); err != nil {
			b.Fatal(err)
		}
	}
}

func BenchmarkSolveLU(b *testing.B) {
	m := matrixBenchMatrix(16)
	rhs := make([]algebra.Expr, 16)
	for i := range rhs {
		rhs[i] = algebra.Flt(float64(i) + 1)
	}
	vec := NewVector(rhs...)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if _, err := m.SolveLU(vec); err != nil {
			b.Fatal(err)
		}
	}
}

func BenchmarkQR(b *testing.B) {
	m := matrixBenchMatrix(16)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if _, _, err := m.QR(); err != nil {
			b.Fatal(err)
		}
	}
}

func BenchmarkCholesky(b *testing.B) {
	m := matrixBenchMatrix(16)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if _, err := m.Cholesky(); err != nil {
			b.Fatal(err)
		}
	}
}
