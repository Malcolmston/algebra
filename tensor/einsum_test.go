package tensor

import (
	"errors"
	"testing"
)

func TestEinsumMatMul(t *testing.T) {
	a, _ := NewWithData([]int{2, 2}, []float64{1, 2, 3, 4})
	b, _ := NewWithData([]int{2, 2}, []float64{5, 6, 7, 8})
	got, err := Einsum("ij,jk->ik", a, b)
	if err != nil {
		t.Fatal(err)
	}
	want, _ := NewWithData([]int{2, 2}, []float64{19, 22, 43, 50})
	if !got.AlmostEqual(want, tol) {
		t.Fatalf("einsum matmul = %v, want %v", got, want)
	}
}

func TestEinsumTraceAndDiagonal(t *testing.T) {
	a, _ := NewWithData([]int{3, 3}, []float64{1, 2, 3, 4, 5, 6, 7, 8, 9})
	tr, err := Einsum("ii->", a)
	if err != nil {
		t.Fatal(err)
	}
	if v, _ := tr.ScalarValue(); v != 15 {
		t.Fatalf("einsum trace = %v, want 15", v)
	}
	diag, err := Einsum("ii->i", a)
	if err != nil {
		t.Fatal(err)
	}
	if !diag.Equal(mustVec(1, 5, 9)) {
		t.Fatalf("einsum diagonal = %v, want [1 5 9]", diag)
	}
}

func TestEinsumTransposeAndSum(t *testing.T) {
	a, _ := NewWithData([]int{2, 3}, []float64{1, 2, 3, 4, 5, 6})
	tp, err := Einsum("ij->ji", a)
	if err != nil {
		t.Fatal(err)
	}
	if !tp.Equal(a.Transpose()) {
		t.Fatalf("einsum transpose = %v", tp)
	}
	total, err := Einsum("ij->", a)
	if err != nil {
		t.Fatal(err)
	}
	if v, _ := total.ScalarValue(); v != 21 {
		t.Fatalf("einsum full sum = %v, want 21", v)
	}
	// Implicit mode: "ba" reorders to output "ab", i.e. the transpose.
	imp, err := Einsum("ba", a)
	if err != nil {
		t.Fatal(err)
	}
	if !imp.Equal(a.Transpose()) {
		t.Fatalf("implicit transpose = %v", imp)
	}
}

func TestEinsumOuterAndBilinear(t *testing.T) {
	x := FromVector([]float64{1, 2, 3})
	y := FromVector([]float64{4, 5})
	o, err := Einsum("i,j->ij", x, y)
	if err != nil {
		t.Fatal(err)
	}
	if !o.Equal(Outer(x, y)) {
		t.Fatalf("einsum outer mismatch: %v", o)
	}
	// Bilinear form x^T M y with x=[1,2], M=[[1,2],[3,4]], y=[3,4] = 61.
	xx := FromVector([]float64{1, 2})
	m, _ := NewWithData([]int{2, 2}, []float64{1, 2, 3, 4})
	yy := FromVector([]float64{3, 4})
	bl, err := Einsum("i,ij,j->", xx, m, yy)
	if err != nil {
		t.Fatal(err)
	}
	if v, _ := bl.ScalarValue(); v != 61 {
		t.Fatalf("bilinear = %v, want 61", v)
	}
}

func TestEinsumErrors(t *testing.T) {
	a, _ := NewWithData([]int{2, 2}, []float64{1, 2, 3, 4})
	if _, err := Einsum("ij,jk->ik", a); !errors.Is(err, ErrSpec) {
		t.Fatalf("wrong term count err = %v, want ErrSpec", err)
	}
	if _, err := Einsum("ijk->i", a); !errors.Is(err, ErrSpec) {
		t.Fatalf("rank mismatch err = %v, want ErrSpec", err)
	}
	// Shared label bound to differing sizes.
	b, _ := NewWithData([]int{3, 3}, []float64{1, 2, 3, 4, 5, 6, 7, 8, 9})
	if _, err := Einsum("ij,jk->ik", a, b); !errors.Is(err, ErrShape) {
		t.Fatalf("size clash err = %v, want ErrShape", err)
	}
	if _, err := Einsum("ij->iz", a); !errors.Is(err, ErrSpec) {
		t.Fatalf("unknown output err = %v, want ErrSpec", err)
	}
}

func BenchmarkEinsumMatMul(b *testing.B) {
	const n = 24
	data := make([]float64, n*n)
	for i := range data {
		data[i] = float64(i%7) - 3
	}
	x, _ := NewWithData([]int{n, n}, data)
	y, _ := NewWithData([]int{n, n}, data)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if _, err := Einsum("ij,jk->ik", x, y); err != nil {
			b.Fatal(err)
		}
	}
}
