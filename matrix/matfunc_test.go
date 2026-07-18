package matrix

import (
	"errors"
	"math"
	"testing"

	"github.com/malcolmston/algebra"
)

// expTol is the absolute tolerance for numeric matrix-exponential comparisons.
const expTol = 1e-9

// matExpFloats evaluates every entry of m to a float64, failing the test on a
// non-numeric entry.
func matExpFloats(t *testing.T, m *Matrix) [][]float64 {
	t.Helper()
	f, err := m.Floats()
	if err != nil {
		t.Fatalf("Floats: %v", err)
	}
	return f
}

// matExpEqual reports whether the numeric matrix m matches want within expTol.
func matExpEqual(t *testing.T, m *Matrix, want [][]float64) {
	t.Helper()
	got := matExpFloats(t, m)
	if len(got) != len(want) {
		t.Fatalf("shape: got %d rows, want %d", len(got), len(want))
	}
	for i := range want {
		if len(got[i]) != len(want[i]) {
			t.Fatalf("row %d: got %d cols, want %d", i, len(got[i]), len(want[i]))
		}
		for j := range want[i] {
			if math.Abs(got[i][j]-want[i][j]) > expTol {
				t.Fatalf("entry (%d,%d): got %.12g, want %.12g", i, j, got[i][j], want[i][j])
			}
		}
	}
}

func TestExpKnownAnswers(t *testing.T) {
	e := math.E
	c1, s1 := math.Cos(1), math.Sin(1)
	tests := []struct {
		name string
		in   [][]float64
		want [][]float64
	}{
		{
			name: "zero-is-identity",
			in:   [][]float64{{0, 0}, {0, 0}},
			want: [][]float64{{1, 0}, {0, 1}},
		},
		{
			name: "identity-gives-e-on-diagonal",
			in:   [][]float64{{1, 0}, {0, 1}},
			want: [][]float64{{e, 0}, {0, e}},
		},
		{
			name: "diagonal",
			in:   [][]float64{{1, 0}, {0, 2}},
			want: [][]float64{{e, 0}, {0, e * e}},
		},
		{
			name: "nilpotent-jordan-block",
			in:   [][]float64{{0, 1}, {0, 0}},
			want: [][]float64{{1, 1}, {0, 1}},
		},
		{
			name: "shifted-jordan-block",
			// [[1,1],[0,1]] has exp e*[[1,1],[0,1]].
			in:   [][]float64{{1, 1}, {0, 1}},
			want: [][]float64{{e, e}, {0, e}},
		},
		{
			name: "rotation-generator",
			// exp([[0,-1],[1,0]]) is the rotation by 1 radian.
			in:   [][]float64{{0, -1}, {1, 0}},
			want: [][]float64{{c1, -s1}, {s1, c1}},
		},
		{
			name: "large-norm-needs-scaling",
			// Diagonal with a large entry stresses scaling-and-squaring.
			in:   [][]float64{{5, 0}, {0, -3}},
			want: [][]float64{{math.Exp(5), 0}, {0, math.Exp(-3)}},
		},
		{
			name: "one-by-one",
			in:   [][]float64{{2}},
			want: [][]float64{{e * e}},
		},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			got, err := FromFloats(tc.in).Exp()
			if err != nil {
				t.Fatalf("Exp: %v", err)
			}
			matExpEqual(t, got, tc.want)
		})
	}
}

func TestExpZeroExactIdentity(t *testing.T) {
	got, err := FromFloats([][]float64{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}).Exp()
	if err != nil {
		t.Fatalf("Exp: %v", err)
	}
	f := matExpFloats(t, got)
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			want := 0.0
			if i == j {
				want = 1.0
			}
			if f[i][j] != want {
				t.Fatalf("entry (%d,%d): got %.17g, want exactly %.17g", i, j, f[i][j], want)
			}
		}
	}
}

// TestExpNonDiagonalizable checks a full 2×2 with coupling against a
// closed-form reference computed from its eigen-decomposition.
func TestExpNonDiagonalizable(t *testing.T) {
	// A = [[1,2],[3,2]]; eigenvalues 4 and -1.
	// exp(A) via A = V diag(e^4,e^-1) V^-1.
	a := [][]float64{{1, 2}, {3, 2}}
	// Reference from eigendecomposition:
	// eigenvectors: for 4 -> (2,3), for -1 -> (1,-1).
	l1, l2 := 4.0, -1.0
	e1, e2 := math.Exp(l1), math.Exp(l2)
	// V = [[2,1],[3,-1]], det = -2-3 = -5, V^-1 = 1/-5 [[-1,-1],[-3,2]]
	// = [[0.2,0.2],[0.6,-0.4]]
	V := [][]float64{{2, 1}, {3, -1}}
	Vinv := [][]float64{{0.2, 0.2}, {0.6, -0.4}}
	D := [][]float64{{e1, 0}, {0, e2}}
	// want = V D Vinv
	want := matMul2(matMul2(V, D), Vinv)

	got, err := FromFloats(a).Exp()
	if err != nil {
		t.Fatalf("Exp: %v", err)
	}
	matExpEqual(t, got, want)
}

func matMul2(x, y [][]float64) [][]float64 {
	n := len(x)
	out := make([][]float64, n)
	for i := range out {
		out[i] = make([]float64, n)
		for j := 0; j < n; j++ {
			var s float64
			for k := 0; k < n; k++ {
				s += x[i][k] * y[k][j]
			}
			out[i][j] = s
		}
	}
	return out
}

// TestExpSemigroup verifies the identity exp(A)·exp(-A) = I.
func TestExpSemigroup(t *testing.T) {
	a := FromFloats([][]float64{{0.3, 1.2, -0.5}, {0.1, -0.4, 0.8}, {-0.7, 0.2, 0.6}})
	expA, err := a.Exp()
	if err != nil {
		t.Fatalf("Exp: %v", err)
	}
	expNegA, err := a.Neg().Exp()
	if err != nil {
		t.Fatalf("Exp neg: %v", err)
	}
	prod, err := expA.Mul(expNegA)
	if err != nil {
		t.Fatalf("Mul: %v", err)
	}
	matExpEqual(t, prod, [][]float64{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}})
}

func TestExpScaledKnownAnswers(t *testing.T) {
	a := FromFloats([][]float64{{1, 0}, {0, 2}})
	t.Run("t=0-is-identity", func(t *testing.T) {
		got, err := a.ExpScaled(0)
		if err != nil {
			t.Fatalf("ExpScaled: %v", err)
		}
		matExpEqual(t, got, [][]float64{{1, 0}, {0, 1}})
	})
	t.Run("t=0.5", func(t *testing.T) {
		got, err := a.ExpScaled(0.5)
		if err != nil {
			t.Fatalf("ExpScaled: %v", err)
		}
		matExpEqual(t, got, [][]float64{{math.Exp(0.5), 0}, {0, math.Exp(1.0)}})
	})
	t.Run("t=-2", func(t *testing.T) {
		got, err := a.ExpScaled(-2)
		if err != nil {
			t.Fatalf("ExpScaled: %v", err)
		}
		matExpEqual(t, got, [][]float64{{math.Exp(-2), 0}, {0, math.Exp(-4)}})
	})
}

// TestExpScaledConsistentWithExp checks ExpScaled(1) == Exp.
func TestExpScaledConsistentWithExp(t *testing.T) {
	a := FromFloats([][]float64{{0.5, -1.5}, {2.0, 0.25}})
	exp, err := a.Exp()
	if err != nil {
		t.Fatalf("Exp: %v", err)
	}
	scaled, err := a.ExpScaled(1)
	if err != nil {
		t.Fatalf("ExpScaled: %v", err)
	}
	matExpEqual(t, scaled, matExpFloats(t, exp))
}

func TestExpNotSquare(t *testing.T) {
	a := FromFloats([][]float64{{1, 2, 3}, {4, 5, 6}})
	if _, err := a.Exp(); !errors.Is(err, ErrNotSquare) {
		t.Fatalf("Exp: err = %v, want ErrNotSquare", err)
	}
	if _, err := a.ExpScaled(2); !errors.Is(err, ErrNotSquare) {
		t.Fatalf("ExpScaled: err = %v, want ErrNotSquare", err)
	}
}

func TestExpSymbolicUnsupported(t *testing.T) {
	sym := FromExpr([][]algebra.Expr{
		{algebra.Sym("a"), algebra.Int(1)},
		{algebra.Int(0), algebra.Int(2)},
	})
	if _, err := sym.Exp(); !errors.Is(err, ErrUnsupported) {
		t.Fatalf("Exp: err = %v, want ErrUnsupported", err)
	}
	if _, err := sym.ExpScaled(1); !errors.Is(err, ErrUnsupported) {
		t.Fatalf("ExpScaled: err = %v, want ErrUnsupported", err)
	}
}

// TestExpDeterministic checks that repeated calls produce bit-identical output.
func TestExpDeterministic(t *testing.T) {
	a := FromFloats([][]float64{{0.4, 1.1, -2.3}, {0.9, -0.2, 0.7}, {1.5, -0.8, 0.3}})
	first, err := a.Exp()
	if err != nil {
		t.Fatalf("Exp: %v", err)
	}
	f0 := matExpFloats(t, first)
	for iter := 0; iter < 3; iter++ {
		next, err := a.Exp()
		if err != nil {
			t.Fatalf("Exp: %v", err)
		}
		fn := matExpFloats(t, next)
		for i := range f0 {
			for j := range f0[i] {
				if f0[i][j] != fn[i][j] {
					t.Fatalf("non-deterministic at (%d,%d): %v vs %v", i, j, f0[i][j], fn[i][j])
				}
			}
		}
	}
}

func benchExpInput() *Matrix {
	return FromFloats([][]float64{
		{0.4, 1.1, -2.3, 0.6},
		{0.9, -0.2, 0.7, -1.4},
		{1.5, -0.8, 0.3, 0.2},
		{-0.5, 0.9, 1.2, -0.7},
	})
}

func BenchmarkExp(b *testing.B) {
	a := benchExpInput()
	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if _, err := a.Exp(); err != nil {
			b.Fatalf("Exp: %v", err)
		}
	}
}

func BenchmarkExpScaled(b *testing.B) {
	a := benchExpInput()
	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if _, err := a.ExpScaled(2.5); err != nil {
			b.Fatalf("ExpScaled: %v", err)
		}
	}
}
