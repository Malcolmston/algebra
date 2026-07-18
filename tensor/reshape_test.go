package tensor

import (
	"errors"
	"testing"
)

func TestTranspose(t *testing.T) {
	x, _ := NewWithData([]int{2, 3}, []float64{1, 2, 3, 4, 5, 6})
	xt := x.Transpose()
	want, _ := NewWithData([]int{3, 2}, []float64{1, 4, 2, 5, 3, 6})
	if !xt.Equal(want) {
		t.Fatalf("transpose = %v, want %v", xt, want)
	}
}

func TestPermute(t *testing.T) {
	// 2x2x2 tensor with values 0..7 in row-major order.
	x, _ := NewWithData([]int{2, 2, 2}, []float64{0, 1, 2, 3, 4, 5, 6, 7})
	p, err := x.Permute(2, 0, 1)
	if err != nil {
		t.Fatal(err)
	}
	// Result axis order (k,i,j): p[k,i,j] = x[i,j,k].
	for i := 0; i < 2; i++ {
		for j := 0; j < 2; j++ {
			for k := 0; k < 2; k++ {
				if p.At(k, i, j) != x.At(i, j, k) {
					t.Fatalf("permute mismatch at %d,%d,%d", i, j, k)
				}
			}
		}
	}
	if _, err := x.Permute(0, 0, 1); !errors.Is(err, ErrAxis) {
		t.Fatalf("err = %v, want ErrAxis", err)
	}
}

func TestSwapAndMoveAxis(t *testing.T) {
	x, _ := NewWithData([]int{2, 3}, []float64{1, 2, 3, 4, 5, 6})
	s, err := x.SwapAxes(0, 1)
	if err != nil {
		t.Fatal(err)
	}
	if !s.Equal(x.Transpose()) {
		t.Fatal("SwapAxes(0,1) should equal Transpose for a matrix")
	}
	y, _ := NewWithData([]int{2, 3, 4}, make([]float64, 24))
	m, err := y.MoveAxis(0, 2)
	if err != nil {
		t.Fatal(err)
	}
	if got := m.Shape(); got[0] != 3 || got[1] != 4 || got[2] != 2 {
		t.Fatalf("MoveAxis shape = %v, want [3 4 2]", got)
	}
}

func TestExpandSqueeze(t *testing.T) {
	x := FromVector([]float64{1, 2, 3})
	e, err := x.ExpandDims(0)
	if err != nil {
		t.Fatal(err)
	}
	if e.Rank() != 2 || e.Dim(0) != 1 || e.Dim(1) != 3 {
		t.Fatalf("ExpandDims shape = %v", e.Shape())
	}
	sq := e.Squeeze()
	if !sq.Equal(x) {
		t.Fatalf("Squeeze did not recover vector: %v", sq)
	}
	// A tensor of all length-1 axes squeezes to a scalar.
	one := New(1, 1, 1)
	one.Set(5, 0, 0, 0)
	if s := one.Squeeze(); !s.IsScalar() || s.At() != 5 {
		t.Fatalf("Squeeze to scalar failed: %v", s)
	}
}

func TestConcatenate(t *testing.T) {
	a, _ := NewWithData([]int{2, 2}, []float64{1, 2, 3, 4})
	b, _ := NewWithData([]int{1, 2}, []float64{5, 6})
	c, err := Concatenate(0, a, b)
	if err != nil {
		t.Fatal(err)
	}
	want, _ := NewWithData([]int{3, 2}, []float64{1, 2, 3, 4, 5, 6})
	if !c.Equal(want) {
		t.Fatalf("Concatenate axis 0 = %v, want %v", c, want)
	}
	d, _ := NewWithData([]int{2, 1}, []float64{7, 8})
	c2, err := Concatenate(1, a, d)
	if err != nil {
		t.Fatal(err)
	}
	want2, _ := NewWithData([]int{2, 3}, []float64{1, 2, 7, 3, 4, 8})
	if !c2.Equal(want2) {
		t.Fatalf("Concatenate axis 1 = %v, want %v", c2, want2)
	}
	if _, err := Concatenate(0, a, d); !errors.Is(err, ErrShape) {
		t.Fatalf("err = %v, want ErrShape", err)
	}
}

func TestStack(t *testing.T) {
	a := FromVector([]float64{1, 2, 3})
	b := FromVector([]float64{4, 5, 6})
	s, err := Stack(0, a, b)
	if err != nil {
		t.Fatal(err)
	}
	want, _ := NewWithData([]int{2, 3}, []float64{1, 2, 3, 4, 5, 6})
	if !s.Equal(want) {
		t.Fatalf("Stack axis 0 = %v, want %v", s, want)
	}
	s1, err := Stack(1, a, b)
	if err != nil {
		t.Fatal(err)
	}
	want1, _ := NewWithData([]int{3, 2}, []float64{1, 4, 2, 5, 3, 6})
	if !s1.Equal(want1) {
		t.Fatalf("Stack axis 1 = %v, want %v", s1, want1)
	}
}
