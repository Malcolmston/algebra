package matrix

import (
	"errors"
	"math"
	"testing"

	"github.com/malcolmston/algebra"
)

// eigTol is the absolute tolerance for numeric eigenvalue comparisons.
const eigTol = 1e-7

func eigFloatsOf(t *testing.T, m *Matrix) [][]float64 {
	t.Helper()
	f, err := m.Floats()
	if err != nil {
		t.Fatalf("Floats: %v", err)
	}
	return f
}

func eigApproxSlice(a, b []float64, tol float64) bool {
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

// --- EigSym / EigSymValues ---------------------------------------------------

func TestEigSymValuesKnown(t *testing.T) {
	cases := []struct {
		name string
		in   [][]float64
		want []float64
	}{
		{
			name: "diagonal",
			in:   [][]float64{{2, 0, 0}, {0, 3, 0}, {0, 0, 1}},
			want: []float64{1, 2, 3},
		},
		{
			name: "2x2 symmetric",
			in:   [][]float64{{2, 1}, {1, 2}},
			want: []float64{1, 3},
		},
		{
			name: "all-ones-plus-diag", // [[2,1,1],[1,2,1],[1,1,2]] -> 1,1,4
			in:   [][]float64{{2, 1, 1}, {1, 2, 1}, {1, 1, 2}},
			want: []float64{1, 1, 4},
		},
		{
			name: "4x4 symmetric", // diag(4)+off; eigenvalues of J-form
			in:   [][]float64{{2, 1, 0, 0}, {1, 2, 1, 0}, {0, 1, 2, 1}, {0, 0, 1, 2}},
			want: []float64{
				2 - 2*math.Cos(math.Pi*1/5),
				2 - 2*math.Cos(math.Pi*2/5),
				2 - 2*math.Cos(math.Pi*3/5),
				2 - 2*math.Cos(math.Pi*4/5),
			},
		},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			m := FromFloats(tc.in)
			got, err := m.EigSymValues()
			if err != nil {
				t.Fatalf("EigSymValues: %v", err)
			}
			if !eigApproxSlice(got, tc.want, eigTol) {
				t.Fatalf("values = %v, want %v", got, tc.want)
			}
		})
	}
}

func TestEigSymVectorsResidual(t *testing.T) {
	inputs := [][][]float64{
		{{2, 1}, {1, 2}},
		{{2, 1, 1}, {1, 2, 1}, {1, 1, 2}},
		{{4, 1, -2}, {1, 3, 0}, {-2, 0, 5}},
	}
	for _, in := range inputs {
		m := FromFloats(in)
		vals, vecs, err := m.EigSym()
		if err != nil {
			t.Fatalf("EigSym: %v", err)
		}
		n := len(in)
		if vecs.Rows() != n || vecs.Cols() != n {
			t.Fatalf("vectors shape = %dx%d, want %dx%d", vecs.Rows(), vecs.Cols(), n, n)
		}
		// Ascending order check.
		for i := 1; i < n; i++ {
			if vals[i]+eigTol < vals[i-1] {
				t.Fatalf("values not ascending: %v", vals)
			}
		}
		a := in
		vf := eigFloatsOf(t, vecs)
		// Residual: A·v_i - λ_i·v_i ~ 0 for each column i.
		for c := 0; c < n; c++ {
			for r := 0; r < n; r++ {
				var av float64
				for k := 0; k < n; k++ {
					av += a[r][k] * vf[k][c]
				}
				if math.Abs(av-vals[c]*vf[r][c]) > 1e-6 {
					t.Fatalf("residual too large at col %d row %d: %g", c, r, av-vals[c]*vf[r][c])
				}
			}
			// Orthonormal columns: unit norm.
			var norm float64
			for r := 0; r < n; r++ {
				norm += vf[r][c] * vf[r][c]
			}
			if math.Abs(norm-1) > 1e-6 {
				t.Fatalf("column %d not unit norm: %g", c, norm)
			}
		}
	}
}

func TestEigSymErrors(t *testing.T) {
	// Non-square.
	if _, _, err := FromFloats([][]float64{{1, 2, 3}, {4, 5, 6}}).EigSym(); !errors.Is(err, ErrNotSquare) {
		t.Fatalf("non-square: got %v, want ErrNotSquare", err)
	}
	if _, err := FromFloats([][]float64{{1, 2, 3}, {4, 5, 6}}).EigSymValues(); !errors.Is(err, ErrNotSquare) {
		t.Fatalf("non-square values: got %v, want ErrNotSquare", err)
	}
	// Not symmetric.
	if _, _, err := FromFloats([][]float64{{0, 1}, {-1, 0}}).EigSym(); !errors.Is(err, ErrUnsupported) {
		t.Fatalf("asymmetric: got %v, want ErrUnsupported", err)
	}
	// Symbolic entries.
	sym := FromExpr([][]algebra.Expr{
		{algebra.Sym("x"), algebra.Int(1)},
		{algebra.Int(1), algebra.Int(2)},
	})
	if _, _, err := sym.EigSym(); !errors.Is(err, ErrUnsupported) {
		t.Fatalf("symbolic: got %v, want ErrUnsupported", err)
	}
}

// --- EigenvaluesNumeric ------------------------------------------------------

func complexApproxUnordered(got, want []complex128, tol float64) bool {
	if len(got) != len(want) {
		return false
	}
	used := make([]bool, len(want))
	for _, g := range got {
		found := false
		for j, w := range want {
			if used[j] {
				continue
			}
			if math.Abs(real(g)-real(w)) <= tol && math.Abs(imag(g)-imag(w)) <= tol {
				used[j] = true
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

func TestEigenvaluesNumericKnown(t *testing.T) {
	cases := []struct {
		name string
		in   [][]float64
		want []complex128 // in the deterministic sorted order the function returns
	}{
		{
			name: "diagonal",
			in:   [][]float64{{2, 0, 0}, {0, 3, 0}, {0, 0, 4}},
			want: []complex128{2, 3, 4},
		},
		{
			name: "upper-triangular",
			in:   [][]float64{{1, 2, 3}, {0, 4, 5}, {0, 0, 6}},
			want: []complex128{1, 4, 6},
		},
		{
			name: "companion-real", // roots 1,2,3
			in:   [][]float64{{0, 1, 0}, {0, 0, 1}, {6, -11, 6}},
			want: []complex128{1, 2, 3},
		},
		{
			name: "pure-imaginary", // x^2+1 -> +-i
			in:   [][]float64{{0, 1}, {-1, 0}},
			want: []complex128{complex(0, -1), complex(0, 1)},
		},
		{
			name: "complex-block", // 1+-i
			in:   [][]float64{{1, -1}, {1, 1}},
			want: []complex128{complex(1, -1), complex(1, 1)},
		},
		{
			name: "mixed-real-complex", // (x-2)(x^2+1) -> -i, i, 2
			in:   [][]float64{{0, 1, 0}, {0, 0, 1}, {2, -1, 2}},
			want: []complex128{complex(0, -1), complex(0, 1), complex(2, 0)},
		},
		{
			name: "4x4-companion", // roots 1,2,3,4
			in: [][]float64{
				{0, 1, 0, 0},
				{0, 0, 1, 0},
				{0, 0, 0, 1},
				{-24, 50, -35, 10},
			},
			want: []complex128{1, 2, 3, 4},
		},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			m := FromFloats(tc.in)
			got, err := m.EigenvaluesNumeric()
			if err != nil {
				t.Fatalf("EigenvaluesNumeric: %v", err)
			}
			if len(got) != len(tc.want) {
				t.Fatalf("len = %d, want %d (%v)", len(got), len(tc.want), got)
			}
			// Value correctness (as a set).
			if !complexApproxUnordered(got, tc.want, 1e-6) {
				t.Fatalf("values = %v, want %v", got, tc.want)
			}
			// Deterministic order: sorted by real then imaginary part.
			for i := 1; i < len(got); i++ {
				if real(got[i]) < real(got[i-1])-1e-9 {
					t.Fatalf("not sorted by real part: %v", got)
				}
				if math.Abs(real(got[i])-real(got[i-1])) <= 1e-9 && imag(got[i]) < imag(got[i-1])-1e-9 {
					t.Fatalf("conjugate pair not ordered by imag: %v", got)
				}
			}
			// Determinism: repeated call yields identical results.
			got2, _ := m.EigenvaluesNumeric()
			for i := range got {
				if got[i] != got2[i] {
					t.Fatalf("non-deterministic result: %v vs %v", got, got2)
				}
			}
		})
	}
}

func TestEigenvaluesNumericErrors(t *testing.T) {
	if _, err := FromFloats([][]float64{{1, 2, 3}, {4, 5, 6}}).EigenvaluesNumeric(); !errors.Is(err, ErrNotSquare) {
		t.Fatalf("non-square: got %v, want ErrNotSquare", err)
	}
	sym := FromExpr([][]algebra.Expr{
		{algebra.Sym("y"), algebra.Int(1)},
		{algebra.Int(0), algebra.Int(2)},
	})
	if _, err := sym.EigenvaluesNumeric(); !errors.Is(err, ErrUnsupported) {
		t.Fatalf("symbolic: got %v, want ErrUnsupported", err)
	}
}

// --- benchmarks --------------------------------------------------------------

// benchSymMatrix builds a deterministic n×n symmetric matrix.
func benchSymMatrix(n int) *Matrix {
	rows := make([][]float64, n)
	for i := 0; i < n; i++ {
		rows[i] = make([]float64, n)
	}
	for i := 0; i < n; i++ {
		for j := i; j < n; j++ {
			v := math.Sin(float64(i*31+j*7+1)) + float64(i+j)*0.01
			rows[i][j] = v
			rows[j][i] = v
		}
	}
	return FromFloats(rows)
}

// benchGeneralMatrix builds a deterministic n×n general (non-symmetric) matrix.
func benchGeneralMatrix(n int) *Matrix {
	rows := make([][]float64, n)
	for i := 0; i < n; i++ {
		rows[i] = make([]float64, n)
		for j := 0; j < n; j++ {
			rows[i][j] = math.Cos(float64(i*17+j*5+3)) + float64(i-j)*0.02
		}
	}
	return FromFloats(rows)
}

func BenchmarkEigSym(b *testing.B) {
	m := benchSymMatrix(24)
	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if _, _, err := m.EigSym(); err != nil {
			b.Fatal(err)
		}
	}
}

func BenchmarkEigSymValues(b *testing.B) {
	m := benchSymMatrix(24)
	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if _, err := m.EigSymValues(); err != nil {
			b.Fatal(err)
		}
	}
}

func BenchmarkEigenvaluesNumeric(b *testing.B) {
	m := benchGeneralMatrix(24)
	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if _, err := m.EigenvaluesNumeric(); err != nil {
			b.Fatal(err)
		}
	}
}
