package matrix

import (
	"errors"
	"math"
	"testing"

	"github.com/malcolmston/algebra"
)

// lstsqFloatsOf evaluates every component of v to float64 for comparison.
func lstsqFloatsOf(t *testing.T, v *Vector) []float64 {
	t.Helper()
	out := make([]float64, v.Len())
	for i := 0; i < v.Len(); i++ {
		f, err := algebra.Evalf(v.At(i))
		if err != nil {
			t.Fatalf("Evalf component %d: %v", i, err)
		}
		out[i] = f
	}
	return out
}

func lstsqClose(a, b []float64, tol float64) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if math.Abs(a[i]-b[i]) > tol {
			return false
		}
	}
	return true
}

func TestLeastSquaresKnown(t *testing.T) {
	const tol = 1e-9
	cases := []struct {
		name string
		a    [][]float64
		b    []float64
		want []float64
	}{
		{
			// Square, exact: diag(2,3) x = (4,9).
			name: "square",
			a:    [][]float64{{2, 0}, {0, 3}},
			b:    []float64{4, 9},
			want: []float64{2, 3},
		},
		{
			// Overdetermined straight-line fit y = c0 + c1*x through
			// (1,6),(2,5),(3,7),(4,10); normal equations give c0=3.5, c1=1.4.
			name: "overdetermined-line",
			a:    [][]float64{{1, 1}, {1, 2}, {1, 3}, {1, 4}},
			b:    []float64{6, 5, 7, 10},
			want: []float64{3.5, 1.4},
		},
		{
			// Underdetermined x1+x2 = 2; minimum-norm solution is (1,1).
			name: "underdetermined-simple",
			a:    [][]float64{{1, 1}},
			b:    []float64{2},
			want: []float64{1, 1},
		},
		{
			// Underdetermined: x1+2x2 = 3, x3 = 5; minimum-norm solution
			// is (0.6, 1.2, 5).
			name: "underdetermined-two-rows",
			a:    [][]float64{{1, 2, 0}, {0, 0, 1}},
			b:    []float64{3, 5},
			want: []float64{0.6, 1.2, 5},
		},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			a := FromFloats(tc.a)
			b := NewVector(lstsqFloatsToExpr(tc.b)...)
			x, err := LeastSquares(a, b)
			if err != nil {
				t.Fatalf("LeastSquares: %v", err)
			}
			got := lstsqFloatsOf(t, x)
			if !lstsqClose(got, tc.want, tol) {
				t.Fatalf("x = %v, want %v", got, tc.want)
			}
		})
	}
}

// TestLeastSquaresResidualOrthogonal checks the defining property of the
// least-squares solution for an overdetermined system: the residual A·x − b is
// orthogonal to every column of A.
func TestLeastSquaresResidualOrthogonal(t *testing.T) {
	a := FromFloats([][]float64{
		{1, 0, 1},
		{0, 1, 1},
		{1, 1, 0},
		{2, 0, 1},
	})
	b := VectorFromInts(1, 2, 3, 4)
	x, err := LeastSquares(a, b)
	if err != nil {
		t.Fatalf("LeastSquares: %v", err)
	}
	af, _ := a.Floats()
	xf := lstsqFloatsOf(t, x)
	m, n := a.Rows(), a.Cols()
	// residual r = A x - b
	r := make([]float64, m)
	for i := 0; i < m; i++ {
		s := 0.0
		for j := 0; j < n; j++ {
			s += af[i][j] * xf[j]
		}
		bi, _ := algebra.Evalf(b.At(i))
		r[i] = s - bi
	}
	// Aᵀ r must be ~0.
	for j := 0; j < n; j++ {
		s := 0.0
		for i := 0; i < m; i++ {
			s += af[i][j] * r[i]
		}
		if math.Abs(s) > 1e-9 {
			t.Fatalf("column %d not orthogonal to residual: %v", j, s)
		}
	}
}

func TestLeastSquaresDimension(t *testing.T) {
	a := FromFloats([][]float64{{1, 2}, {3, 4}})
	b := VectorFromInts(1, 2, 3)
	if _, err := LeastSquares(a, b); !errors.Is(err, ErrDimension) {
		t.Fatalf("err = %v, want ErrDimension", err)
	}
}

func TestLeastSquaresSymbolicUnsupported(t *testing.T) {
	a := FromExpr([][]algebra.Expr{
		{algebra.Sym("a"), algebra.Int(1)},
		{algebra.Int(2), algebra.Int(3)},
	})
	b := VectorFromInts(1, 2)
	if _, err := LeastSquares(a, b); !errors.Is(err, ErrUnsupported) {
		t.Fatalf("err = %v, want ErrUnsupported", err)
	}
}

func TestLeastSquaresSingular(t *testing.T) {
	// Rank-deficient overdetermined system (columns are identical).
	a := FromFloats([][]float64{{1, 1}, {2, 2}, {3, 3}})
	b := VectorFromInts(1, 2, 3)
	if _, err := LeastSquares(a, b); !errors.Is(err, ErrSingular) {
		t.Fatalf("err = %v, want ErrSingular", err)
	}
}

func TestNullspaceExactKnown(t *testing.T) {
	// A = [[1,2,3],[4,5,6]] has rank 2, nullity 1. Row-reduce:
	// RREF = [[1,0,-1],[0,1,2]]; free var x3 -> basis vector [1,-2,1].
	a := FromInts([][]int64{{1, 2, 3}, {4, 5, 6}})
	basis, err := a.NullspaceExact()
	if err != nil {
		t.Fatalf("NullspaceExact: %v", err)
	}
	if len(basis) != 1 {
		t.Fatalf("nullity = %d, want 1", len(basis))
	}
	want := NewVector(algebra.Int(1), algebra.Int(-2), algebra.Int(1))
	if !basis[0].Equal(want) {
		t.Fatalf("basis[0] = %v, want %v", basis[0], want)
	}
	// A·v must be zero for every basis vector.
	for _, v := range basis {
		prod, _ := a.MulVec(v)
		for i := 0; i < prod.Len(); i++ {
			if !isZeroExpr(prod.At(i)) {
				t.Fatalf("A·basis not zero: %v", prod)
			}
		}
	}
}

func TestNullspaceExactFullRank(t *testing.T) {
	a := FromInts([][]int64{{1, 0}, {0, 1}})
	basis, err := a.NullspaceExact()
	if err != nil {
		t.Fatalf("NullspaceExact: %v", err)
	}
	if len(basis) != 0 {
		t.Fatalf("nullity = %d, want 0", len(basis))
	}
}

func TestNullspaceExactSymbolic(t *testing.T) {
	// [[1, a]] has nullspace spanned by [-a, 1] (free column is index 1).
	a := FromExpr([][]algebra.Expr{{algebra.Int(1), algebra.Sym("a")}})
	basis, err := a.NullspaceExact()
	if err != nil {
		t.Fatalf("NullspaceExact: %v", err)
	}
	if len(basis) != 1 {
		t.Fatalf("nullity = %d, want 1", len(basis))
	}
	want := NewVector(algebra.Mul(algebra.Int(-1), algebra.Sym("a")), algebra.Int(1))
	if !basis[0].Equal(want) {
		t.Fatalf("basis[0] = %v, want %v", basis[0], want)
	}
}

func TestColumnSpaceExactKnown(t *testing.T) {
	// A = [[1,2,3],[4,5,6]] has pivot columns 0 and 1; the basis is the
	// original columns [1,4] and [2,5].
	a := FromInts([][]int64{{1, 2, 3}, {4, 5, 6}})
	basis, err := a.ColumnSpaceExact()
	if err != nil {
		t.Fatalf("ColumnSpaceExact: %v", err)
	}
	if len(basis) != 2 {
		t.Fatalf("dim = %d, want 2", len(basis))
	}
	want0 := NewVector(algebra.Int(1), algebra.Int(4))
	want1 := NewVector(algebra.Int(2), algebra.Int(5))
	if !basis[0].Equal(want0) || !basis[1].Equal(want1) {
		t.Fatalf("basis = %v, %v; want %v, %v", basis[0], basis[1], want0, want1)
	}
}

func TestRowSpaceExactKnown(t *testing.T) {
	// A = [[1,2,3],[4,5,6]] -> RREF [[1,0,-1],[0,1,2]].
	a := FromInts([][]int64{{1, 2, 3}, {4, 5, 6}})
	basis, err := a.RowSpaceExact()
	if err != nil {
		t.Fatalf("RowSpaceExact: %v", err)
	}
	if len(basis) != 2 {
		t.Fatalf("dim = %d, want 2", len(basis))
	}
	want0 := NewVector(algebra.Int(1), algebra.Int(0), algebra.Int(-1))
	want1 := NewVector(algebra.Int(0), algebra.Int(1), algebra.Int(2))
	if !basis[0].Equal(want0) || !basis[1].Equal(want1) {
		t.Fatalf("basis = %v, %v; want %v, %v", basis[0], basis[1], want0, want1)
	}
}

func TestRowSpaceExactRankDeficient(t *testing.T) {
	// Second row is twice the first; row space is one-dimensional.
	a := FromInts([][]int64{{1, 2}, {2, 4}})
	basis, err := a.RowSpaceExact()
	if err != nil {
		t.Fatalf("RowSpaceExact: %v", err)
	}
	if len(basis) != 1 {
		t.Fatalf("dim = %d, want 1", len(basis))
	}
	want := NewVector(algebra.Int(1), algebra.Int(2))
	if !basis[0].Equal(want) {
		t.Fatalf("basis[0] = %v, want %v", basis[0], want)
	}
}

// lstsqFloatsToExpr converts a float slice into algebra.Expr components for NewVector.
func lstsqFloatsToExpr(vals []float64) []algebra.Expr {
	out := make([]algebra.Expr, len(vals))
	for i, v := range vals {
		out[i] = algebra.Flt(v)
	}
	return out
}

func BenchmarkLeastSquaresOverdetermined(b *testing.B) {
	const m, n = 64, 16
	vals := make([][]float64, m)
	for i := 0; i < m; i++ {
		vals[i] = make([]float64, n)
		for j := 0; j < n; j++ {
			vals[i][j] = math.Sin(float64(i*7+j*3+1)) + float64(j)
		}
	}
	a := FromFloats(vals)
	rhs := make([]int64, m)
	for i := 0; i < m; i++ {
		rhs[i] = int64(i%5 + 1)
	}
	bv := VectorFromInts(rhs...)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if _, err := LeastSquares(a, bv); err != nil {
			b.Fatalf("LeastSquares: %v", err)
		}
	}
}

func BenchmarkLeastSquaresUnderdetermined(b *testing.B) {
	const m, n = 16, 64
	vals := make([][]float64, m)
	for i := 0; i < m; i++ {
		vals[i] = make([]float64, n)
		for j := 0; j < n; j++ {
			vals[i][j] = math.Cos(float64(i*5+j*2+1)) + float64(i)
		}
	}
	a := FromFloats(vals)
	rhs := make([]int64, m)
	for i := 0; i < m; i++ {
		rhs[i] = int64(i%3 + 1)
	}
	bv := VectorFromInts(rhs...)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if _, err := LeastSquares(a, bv); err != nil {
			b.Fatalf("LeastSquares: %v", err)
		}
	}
}
