package tensor

import (
	"errors"
	"math"
	"testing"
)

func TestElementwise(t *testing.T) {
	a, _ := NewWithData([]int{2, 2}, []float64{1, 2, 3, 4})
	b, _ := NewWithData([]int{2, 2}, []float64{5, 6, 7, 8})
	add, _ := a.Add(b)
	sub, _ := a.Sub(b)
	mul, _ := a.Mul(b)
	div, _ := b.Div(a)
	wantAdd, _ := NewWithData([]int{2, 2}, []float64{6, 8, 10, 12})
	wantSub, _ := NewWithData([]int{2, 2}, []float64{-4, -4, -4, -4})
	wantMul, _ := NewWithData([]int{2, 2}, []float64{5, 12, 21, 32})
	wantDiv, _ := NewWithData([]int{2, 2}, []float64{5, 3, 7.0 / 3.0, 2})
	if !add.Equal(wantAdd) || !sub.Equal(wantSub) || !mul.Equal(wantMul) {
		t.Fatal("elementwise add/sub/mul mismatch")
	}
	if !div.AlmostEqual(wantDiv, tol) {
		t.Fatalf("div = %v, want %v", div, wantDiv)
	}
	if _, err := a.Add(FromVector([]float64{1, 2})); !errors.Is(err, ErrShape) {
		t.Fatalf("err = %v, want ErrShape", err)
	}
}

func TestScaleNegApply(t *testing.T) {
	a, _ := NewWithData([]int{3}, []float64{1, -2, 3})
	if !a.Scale(2).Equal(mustVec(2, -4, 6)) {
		t.Fatal("Scale wrong")
	}
	if !a.Neg().Equal(mustVec(-1, 2, -3)) {
		t.Fatal("Neg wrong")
	}
	if !a.AddScalar(1).Equal(mustVec(2, -1, 4)) {
		t.Fatal("AddScalar wrong")
	}
	sq := a.Apply(func(x float64) float64 { return x * x })
	if !sq.Equal(mustVec(1, 4, 9)) {
		t.Fatal("Apply wrong")
	}
	if !a.Abs().Equal(mustVec(1, 2, 3)) {
		t.Fatal("Abs wrong")
	}
}

func TestReductions(t *testing.T) {
	a, _ := NewWithData([]int{2, 3}, []float64{1, 2, 3, 4, 5, 6})
	if a.Sum() != 21 {
		t.Fatalf("Sum = %v", a.Sum())
	}
	if a.Mean() != 3.5 {
		t.Fatalf("Mean = %v", a.Mean())
	}
	if a.Product() != 720 {
		t.Fatalf("Product = %v", a.Product())
	}
	if a.Max() != 6 || a.Min() != 1 {
		t.Fatalf("Max/Min = %v/%v", a.Max(), a.Min())
	}
	// Frobenius norm sqrt(1+4+9+16+25+36) = sqrt(91).
	if !almostEqScalar(a.Norm(), math.Sqrt(91)) {
		t.Fatalf("Norm = %v, want sqrt(91)", a.Norm())
	}
}

func TestSumAxis(t *testing.T) {
	a, _ := NewWithData([]int{2, 3}, []float64{1, 2, 3, 4, 5, 6})
	s0, err := a.SumAxis(0)
	if err != nil {
		t.Fatal(err)
	}
	// Column sums: [1+4, 2+5, 3+6] = [5,7,9].
	if !s0.Equal(mustVec(5, 7, 9)) {
		t.Fatalf("SumAxis(0) = %v, want [5 7 9]", s0)
	}
	s1, err := a.SumAxis(1)
	if err != nil {
		t.Fatal(err)
	}
	// Row sums: [1+2+3, 4+5+6] = [6,15].
	if !s1.Equal(mustVec(6, 15)) {
		t.Fatalf("SumAxis(1) = %v, want [6 15]", s1)
	}
	if _, err := a.SumAxis(2); !errors.Is(err, ErrAxis) {
		t.Fatalf("err = %v, want ErrAxis", err)
	}
}

func mustVec(vs ...float64) *Tensor { return FromVector(vs) }
