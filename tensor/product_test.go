package tensor

import (
	"errors"
	"testing"
)

func TestOuter(t *testing.T) {
	a := FromVector([]float64{1, 2, 3})
	b := FromVector([]float64{4, 5})
	o := Outer(a, b)
	want, _ := NewWithData([]int{3, 2}, []float64{4, 5, 8, 10, 12, 15})
	if !o.Equal(want) {
		t.Fatalf("Outer = %v, want %v", o, want)
	}
	if !TensorProduct(a, b).Equal(want) {
		t.Fatal("TensorProduct should equal Outer")
	}
	// Outer of two scalars is their scalar product.
	s := Outer(FromScalar(3), FromScalar(4))
	if v, _ := s.ScalarValue(); v != 12 {
		t.Fatalf("scalar outer = %v, want 12", v)
	}
}

func TestKronecker(t *testing.T) {
	a, _ := NewWithData([]int{2, 2}, []float64{1, 2, 3, 4})
	b, _ := NewWithData([]int{2, 2}, []float64{0, 5, 6, 7})
	k, err := Kronecker(a, b)
	if err != nil {
		t.Fatal(err)
	}
	want, _ := NewWithData([]int{4, 4}, []float64{
		0, 5, 0, 10,
		6, 7, 12, 14,
		0, 15, 0, 20,
		18, 21, 24, 28,
	})
	if !k.Equal(want) {
		t.Fatalf("Kronecker = %v, want %v", k, want)
	}
	if _, err := Kronecker(a, FromVector([]float64{1})); !errors.Is(err, ErrRank) {
		t.Fatalf("err = %v, want ErrRank", err)
	}
}

func TestMatMulAndTensorDot(t *testing.T) {
	a, _ := NewWithData([]int{2, 2}, []float64{1, 2, 3, 4})
	b, _ := NewWithData([]int{2, 2}, []float64{5, 6, 7, 8})
	m, err := MatMul(a, b)
	if err != nil {
		t.Fatal(err)
	}
	want, _ := NewWithData([]int{2, 2}, []float64{19, 22, 43, 50})
	if !m.AlmostEqual(want, tol) {
		t.Fatalf("MatMul = %v, want %v", m, want)
	}
	// TensorDot over inner axes reproduces the matrix product.
	td, err := TensorDot(a, b, []int{1}, []int{0})
	if err != nil {
		t.Fatal(err)
	}
	if !td.AlmostEqual(want, tol) {
		t.Fatalf("TensorDot = %v, want %v", td, want)
	}
	// Non-square: (2x3)(3x2) = 2x2.
	p, _ := NewWithData([]int{2, 3}, []float64{1, 2, 3, 4, 5, 6})
	q, _ := NewWithData([]int{3, 2}, []float64{7, 8, 9, 10, 11, 12})
	pq, err := MatMul(p, q)
	if err != nil {
		t.Fatal(err)
	}
	wantPQ, _ := NewWithData([]int{2, 2}, []float64{58, 64, 139, 154})
	if !pq.AlmostEqual(wantPQ, tol) {
		t.Fatalf("MatMul nonsquare = %v, want %v", pq, wantPQ)
	}
	if _, err := MatMul(p, p); !errors.Is(err, ErrShape) {
		t.Fatalf("err = %v, want ErrShape", err)
	}
}

func TestTensorDotFullContraction(t *testing.T) {
	// Contracting all axes of two equal matrices gives the Frobenius inner
	// product sum_ij a_ij b_ij.
	a, _ := NewWithData([]int{2, 2}, []float64{1, 2, 3, 4})
	b, _ := NewWithData([]int{2, 2}, []float64{5, 6, 7, 8})
	full, err := TensorDot(a, b, []int{0, 1}, []int{0, 1})
	if err != nil {
		t.Fatal(err)
	}
	v, _ := full.ScalarValue()
	if v != 1*5+2*6+3*7+4*8 {
		t.Fatalf("full contraction = %v, want 70", v)
	}
}

func TestDot(t *testing.T) {
	a := FromVector([]float64{1, 2, 3})
	b := FromVector([]float64{4, 5, 6})
	d, err := Dot(a, b)
	if err != nil {
		t.Fatal(err)
	}
	if d != 32 {
		t.Fatalf("Dot = %v, want 32", d)
	}
	if _, err := Dot(a, FromVector([]float64{1, 2})); !errors.Is(err, ErrShape) {
		t.Fatalf("err = %v, want ErrShape", err)
	}
}

func TestTraceAndContract(t *testing.T) {
	m, _ := NewWithData([]int{2, 2}, []float64{1, 2, 3, 4})
	tr, err := m.Trace()
	if err != nil {
		t.Fatal(err)
	}
	if tr != 5 {
		t.Fatalf("Trace = %v, want 5", tr)
	}
	// Contract a rank-3 tensor over axes 0 and 1.
	x, _ := NewWithData([]int{2, 2, 2}, []float64{0, 1, 2, 3, 4, 5, 6, 7})
	c, err := x.Contract(0, 1)
	if err != nil {
		t.Fatal(err)
	}
	// out[k] = x[0,0,k] + x[1,1,k] = [0+6, 1+7] = [6,8].
	if !c.Equal(mustVec(6, 8)) {
		t.Fatalf("Contract(0,1) = %v, want [6 8]", c)
	}
	if _, err := m.Contract(0, 0); !errors.Is(err, ErrAxis) {
		t.Fatalf("err = %v, want ErrAxis", err)
	}
}
