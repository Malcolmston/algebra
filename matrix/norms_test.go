package matrix

import (
	"errors"
	"math"
	"testing"

	"github.com/malcolmston/algebra"
)

func TestNormFro(t *testing.T) {
	tests := []struct {
		name string
		m    *Matrix
		want float64 // numeric value of the exact norm
		exp  string  // expected exact string, "" to skip the structural check
	}{
		{"ones2x2", FromInts([][]int64{{1, 1}, {1, 1}}), 2, "2"},
		{"pythagorean", FromInts([][]int64{{3, 4}}), 5, "5"},
		{"radical8", FromInts([][]int64{{2, 2}, {0, 0}}), math.Sqrt(8), "2*sqrt(2)"},
		{"sqrt30", FromInts([][]int64{{1, 2}, {3, 4}}), math.Sqrt(30), "sqrt(30)"},
		{"empty", New(0, 0), 0, "0"},
		{"negatives", FromInts([][]int64{{-3, 0}, {0, -4}}), 5, "5"},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := tt.m.NormFro()
			numEqual(t, got, tt.want, 1e-12)
			if tt.exp != "" && got.String() != tt.exp {
				t.Fatalf("NormFro() = %q, want %q", got.String(), tt.exp)
			}
		})
	}
}

func TestNormFroSymbolic(t *testing.T) {
	a := algebra.Sym("a")
	b := algebra.Sym("b")
	m := FromExpr([][]algebra.Expr{{a, b}})
	got := m.NormFro()
	// sqrt(a^2 + b^2); confirm it is a sqrt of the sum of squares by evaluating
	// at a=3, b=4.
	val, err := algebra.Eval(got, map[string]float64{"a": 3, "b": 4})
	if err != nil {
		t.Fatalf("Eval: %v", err)
	}
	if math.Abs(val-5) > 1e-12 {
		t.Fatalf("NormFro at (3,4) = %v, want 5", val)
	}
}

func TestNumericNorms(t *testing.T) {
	tests := []struct {
		name                    string
		m                       *Matrix
		norm1, normInf, normMax float64
	}{
		{
			"mixed-signs",
			FromInts([][]int64{{1, -2}, {-3, 4}}),
			6, 7, 4,
		},
		{
			"2x3",
			FromInts([][]int64{{1, 2, 3}, {4, 5, 6}}),
			9, 15, 6,
		},
		{
			"identity3",
			Identity(3),
			1, 1, 1,
		},
		{
			"zeros",
			New(2, 2),
			0, 0, 0,
		},
		{
			"floats",
			FromFloats([][]float64{{1.5, -2.5}, {0.5, 0.5}}),
			3, 4, 2.5,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			n1, err := tt.m.Norm1()
			if err != nil {
				t.Fatalf("Norm1: %v", err)
			}
			if math.Abs(n1-tt.norm1) > 1e-12 {
				t.Fatalf("Norm1 = %v, want %v", n1, tt.norm1)
			}
			nInf, err := tt.m.NormInf()
			if err != nil {
				t.Fatalf("NormInf: %v", err)
			}
			if math.Abs(nInf-tt.normInf) > 1e-12 {
				t.Fatalf("NormInf = %v, want %v", nInf, tt.normInf)
			}
			nMax, err := tt.m.NormMax()
			if err != nil {
				t.Fatalf("NormMax: %v", err)
			}
			if math.Abs(nMax-tt.normMax) > 1e-12 {
				t.Fatalf("NormMax = %v, want %v", nMax, tt.normMax)
			}
		})
	}
}

func TestNumericNormsSymbolicUnsupported(t *testing.T) {
	m := FromExpr([][]algebra.Expr{{algebra.Sym("a"), algebra.Int(1)}})
	if _, err := m.Norm1(); !errors.Is(err, ErrUnsupported) {
		t.Fatalf("Norm1 err = %v, want ErrUnsupported", err)
	}
	if _, err := m.NormInf(); !errors.Is(err, ErrUnsupported) {
		t.Fatalf("NormInf err = %v, want ErrUnsupported", err)
	}
	if _, err := m.NormMax(); !errors.Is(err, ErrUnsupported) {
		t.Fatalf("NormMax err = %v, want ErrUnsupported", err)
	}
}

func TestHadamard(t *testing.T) {
	a := FromInts([][]int64{{1, 2}, {3, 4}})
	b := FromInts([][]int64{{5, 6}, {7, 8}})
	got, err := a.Hadamard(b)
	if err != nil {
		t.Fatal(err)
	}
	want := FromInts([][]int64{{5, 12}, {21, 32}})
	if !got.Equal(want) {
		t.Fatalf("Hadamard = %s, want %s", got, want)
	}
}

func TestHadamardSymbolic(t *testing.T) {
	x := algebra.Sym("x")
	y := algebra.Sym("y")
	a := FromExpr([][]algebra.Expr{{x, y}})
	got, err := a.Hadamard(a)
	if err != nil {
		t.Fatal(err)
	}
	want := FromExpr([][]algebra.Expr{{algebra.Pow(x, algebra.Int(2)), algebra.Pow(y, algebra.Int(2))}})
	if !got.Equal(want) {
		t.Fatalf("Hadamard = %s, want %s", got, want)
	}
}

func TestHadamardDimension(t *testing.T) {
	a := FromInts([][]int64{{1, 2}})
	b := FromInts([][]int64{{1, 2}, {3, 4}})
	if _, err := a.Hadamard(b); !errors.Is(err, ErrDimension) {
		t.Fatalf("err = %v, want ErrDimension", err)
	}
}

func TestCondP(t *testing.T) {
	tests := []struct {
		name string
		m    *Matrix
		p    int
		want float64
	}{
		{"identity-1", Identity(3), 1, 1},
		{"identity-inf", Identity(3), NormInfinity, 1},
		{"diag-1", FromInts([][]int64{{2, 0}, {0, 4}}), 1, 2},
		{"diag-inf", FromInts([][]int64{{2, 0}, {0, 4}}), NormInfinity, 2},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := tt.m.CondP(tt.p)
			if err != nil {
				t.Fatal(err)
			}
			if math.Abs(got-tt.want) > 1e-9 {
				t.Fatalf("CondP = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestCondPSingular(t *testing.T) {
	m := FromInts([][]int64{{1, 2}, {2, 4}})
	got, err := m.CondP(1)
	if !errors.Is(err, ErrSingular) {
		t.Fatalf("err = %v, want ErrSingular", err)
	}
	if !math.IsInf(got, 1) {
		t.Fatalf("CondP = %v, want +Inf", got)
	}
}

func TestCondPNotSquare(t *testing.T) {
	m := FromInts([][]int64{{1, 2, 3}})
	if _, err := m.CondP(1); !errors.Is(err, ErrNotSquare) {
		t.Fatalf("err = %v, want ErrNotSquare", err)
	}
}

func TestCondPUnsupportedSelector(t *testing.T) {
	m := Identity(2)
	if _, err := m.CondP(2); !errors.Is(err, ErrUnsupported) {
		t.Fatalf("err = %v, want ErrUnsupported", err)
	}
}

func benchMatrix(n int) *Matrix {
	vals := make([][]float64, n)
	for i := 0; i < n; i++ {
		vals[i] = make([]float64, n)
		for j := 0; j < n; j++ {
			vals[i][j] = float64((i*7+j*3)%11) - 5
		}
	}
	return FromFloats(vals)
}

func BenchmarkNorm1(b *testing.B) {
	m := benchMatrix(64)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if _, err := m.Norm1(); err != nil {
			b.Fatal(err)
		}
	}
}

func BenchmarkNormInf(b *testing.B) {
	m := benchMatrix(64)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if _, err := m.NormInf(); err != nil {
			b.Fatal(err)
		}
	}
}

func BenchmarkNormMax(b *testing.B) {
	m := benchMatrix(64)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if _, err := m.NormMax(); err != nil {
			b.Fatal(err)
		}
	}
}

func BenchmarkHadamard(b *testing.B) {
	m := benchMatrix(64)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if _, err := m.Hadamard(m); err != nil {
			b.Fatal(err)
		}
	}
}
