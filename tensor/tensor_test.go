package tensor

import (
	"errors"
	"math"
	"testing"
)

const tol = 1e-12

func TestNewAndFill(t *testing.T) {
	z := New(2, 3)
	if z.Rank() != 2 || z.Size() != 6 {
		t.Fatalf("rank/size = %d/%d, want 2/6", z.Rank(), z.Size())
	}
	for _, v := range z.Data() {
		if v != 0 {
			t.Fatalf("Zeros produced nonzero %v", v)
		}
	}
	o := Ones(2, 2)
	if o.Sum() != 4 {
		t.Fatalf("Ones sum = %v, want 4", o.Sum())
	}
	f := Full(2.5, 4)
	if f.Sum() != 10 {
		t.Fatalf("Full sum = %v, want 10", f.Sum())
	}
}

func TestNewWithDataErrors(t *testing.T) {
	if _, err := NewWithData([]int{2, 2}, []float64{1, 2, 3}); !errors.Is(err, ErrDataLength) {
		t.Fatalf("err = %v, want ErrDataLength", err)
	}
	if _, err := NewWithData([]int{2, 0}, []float64{}); !errors.Is(err, ErrShape) {
		t.Fatalf("err = %v, want ErrShape", err)
	}
	got, err := NewWithData([]int{2, 2}, []float64{1, 2, 3, 4})
	if err != nil {
		t.Fatal(err)
	}
	if got.At(1, 0) != 3 {
		t.Fatalf("At(1,0) = %v, want 3", got.At(1, 0))
	}
}

func TestAtSet(t *testing.T) {
	x := New(2, 3)
	x.Set(7, 1, 2)
	if x.At(1, 2) != 7 {
		t.Fatalf("At = %v, want 7", x.At(1, 2))
	}
	// Flat offset for (1,2) in a 2x3 row-major tensor is 1*3+2 = 5.
	if x.Data()[5] != 7 {
		t.Fatalf("flat data[5] = %v, want 7", x.Data()[5])
	}
}

func TestAtPanicsOutOfRange(t *testing.T) {
	defer func() {
		if r := recover(); r == nil {
			t.Fatal("expected panic on out-of-range index")
		}
	}()
	x := New(2, 2)
	_ = x.At(2, 0)
}

func TestScalar(t *testing.T) {
	s := FromScalar(3.14)
	if !s.IsScalar() || s.Rank() != 0 {
		t.Fatalf("scalar rank = %d", s.Rank())
	}
	v, err := s.ScalarValue()
	if err != nil || v != 3.14 {
		t.Fatalf("ScalarValue = %v, %v", v, err)
	}
	if s.At() != 3.14 {
		t.Fatalf("At() = %v", s.At())
	}
}

func TestCloneIndependent(t *testing.T) {
	a := FromVector([]float64{1, 2, 3})
	b := a.Clone()
	b.Set(99, 0)
	if a.At(0) != 1 {
		t.Fatalf("clone mutated original: %v", a.At(0))
	}
}

func TestEqualAndAlmost(t *testing.T) {
	a := FromVector([]float64{1, 2, 3})
	b := FromVector([]float64{1, 2, 3})
	c := FromVector([]float64{1, 2, 3.0000001})
	if !a.Equal(b) {
		t.Fatal("Equal should hold")
	}
	if a.Equal(c) {
		t.Fatal("Equal should fail on tiny diff")
	}
	if !a.AlmostEqual(c, 1e-6) {
		t.Fatal("AlmostEqual should hold within 1e-6")
	}
	if a.AlmostEqual(c, 1e-9) {
		t.Fatal("AlmostEqual should fail within 1e-9")
	}
}

func TestMatrixRoundTrip(t *testing.T) {
	m := [][]float64{{1, 2, 3}, {4, 5, 6}}
	tn, err := FromMatrix(m)
	if err != nil {
		t.Fatal(err)
	}
	back, err := tn.ToMatrix()
	if err != nil {
		t.Fatal(err)
	}
	for i := range m {
		for j := range m[i] {
			if back[i][j] != m[i][j] {
				t.Fatalf("round-trip mismatch at %d,%d", i, j)
			}
		}
	}
	if _, err := FromMatrix([][]float64{{1, 2}, {3}}); !errors.Is(err, ErrShape) {
		t.Fatalf("ragged err = %v, want ErrShape", err)
	}
}

func TestRavelUnravelRoundTrip(t *testing.T) {
	shape := []int{2, 3, 4}
	for flat := 0; flat < 24; flat++ {
		idx, err := UnravelIndex(shape, flat)
		if err != nil {
			t.Fatal(err)
		}
		back, err := RavelIndex(shape, idx)
		if err != nil {
			t.Fatal(err)
		}
		if back != flat {
			t.Fatalf("round trip %d -> %v -> %d", flat, idx, back)
		}
	}
	if _, err := UnravelIndex(shape, 24); !errors.Is(err, ErrIndex) {
		t.Fatalf("err = %v, want ErrIndex", err)
	}
}

func TestReshapeRavel(t *testing.T) {
	x, _ := NewWithData([]int{2, 3}, []float64{1, 2, 3, 4, 5, 6})
	r, err := x.Reshape(3, 2)
	if err != nil {
		t.Fatal(err)
	}
	if r.At(2, 1) != 6 {
		t.Fatalf("reshaped At(2,1) = %v, want 6", r.At(2, 1))
	}
	if _, err := x.Reshape(4, 2); !errors.Is(err, ErrShape) {
		t.Fatalf("err = %v, want ErrShape", err)
	}
	flat := x.Ravel()
	if flat.Rank() != 1 || flat.Size() != 6 || flat.At(5) != 6 {
		t.Fatalf("ravel wrong: %v", flat)
	}
}

func TestStringDeterministic(t *testing.T) {
	x := FromVector([]float64{1, 2})
	if x.String() != "Tensor(shape=[2], data=[1 2])" {
		t.Fatalf("String = %q", x.String())
	}
}

func almostEqScalar(a, b float64) bool { return math.Abs(a-b) <= tol }
