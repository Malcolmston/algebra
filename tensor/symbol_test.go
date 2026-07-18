package tensor

import (
	"testing"
)

func TestKroneckerDelta(t *testing.T) {
	if KroneckerDelta(2, 2) != 1 || KroneckerDelta(1, 3) != 0 {
		t.Fatal("KroneckerDelta wrong")
	}
	id := Identity(3)
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			if id.At(i, j) != KroneckerDelta(i, j) {
				t.Fatalf("Identity mismatch at %d,%d", i, j)
			}
		}
	}
}

func TestLeviCivita(t *testing.T) {
	cases := []struct {
		idx  []int
		want int
	}{
		{[]int{0, 1, 2}, 1},
		{[]int{1, 2, 0}, 1},
		{[]int{2, 0, 1}, 1},
		{[]int{2, 1, 0}, -1},
		{[]int{0, 2, 1}, -1},
		{[]int{1, 0, 2}, -1},
		{[]int{0, 0, 1}, 0},
		{[]int{0, 1}, 1}, // identity permutation -> +1
		{[]int{1, 0}, -1},
		{[]int{}, 1}, // empty product convention
	}
	for _, c := range cases {
		if got := LeviCivita(c.idx...); got != c.want {
			t.Fatalf("LeviCivita(%v) = %d, want %d", c.idx, got, c.want)
		}
	}
}

func TestLeviCivitaTensor(t *testing.T) {
	e := LeviCivitaTensor(3)
	if e.Rank() != 3 {
		t.Fatalf("rank = %d, want 3", e.Rank())
	}
	if e.At(0, 1, 2) != 1 || e.At(2, 1, 0) != -1 || e.At(1, 1, 2) != 0 {
		t.Fatalf("LeviCivitaTensor components wrong")
	}
	// Sum of all components must be zero and sum of absolute values must be 6
	// (there are 3! = 6 nonzero entries, each ±1).
	if e.Sum() != 0 {
		t.Fatalf("sum = %v, want 0", e.Sum())
	}
	if e.Abs().Sum() != 6 {
		t.Fatalf("abs sum = %v, want 6", e.Abs().Sum())
	}
}

func TestCrossProductViaLeviCivita(t *testing.T) {
	// c_i = eps_ijk a_j b_k should reproduce the 3-D cross product.
	e := LeviCivitaTensor(3)
	a := FromVector([]float64{1, 0, 0})
	b := FromVector([]float64{0, 1, 0})
	c, err := Einsum("ijk,j,k->i", e, a, b)
	if err != nil {
		t.Fatal(err)
	}
	if !c.Equal(mustVec(0, 0, 1)) {
		t.Fatalf("cross product = %v, want [0 0 1]", c)
	}
}

func TestMetricRaiseLower(t *testing.T) {
	// Euclidean metric leaves components unchanged.
	g := EuclideanMetric(3)
	v := FromVector([]float64{1, 2, 3})
	low, err := LowerIndex(v, g, 0)
	if err != nil {
		t.Fatal(err)
	}
	if !low.Equal(v) {
		t.Fatalf("Euclidean lower changed vector: %v", low)
	}
	// Minkowski (-,+,+,+) flips the sign of the time component.
	eta := MinkowskiMetric(4)
	u := FromVector([]float64{1, 2, 3, 4})
	lu, err := LowerIndex(u, eta, 0)
	if err != nil {
		t.Fatal(err)
	}
	if !lu.Equal(mustVec(-1, 2, 3, 4)) {
		t.Fatalf("Minkowski lower = %v, want [-1 2 3 4]", lu)
	}
	// Raising with the same (self-inverse) metric returns the original.
	back, err := RaiseIndex(lu, eta, 0)
	if err != nil {
		t.Fatal(err)
	}
	if !back.Equal(u) {
		t.Fatalf("raise(lower(u)) = %v, want %v", back, u)
	}
}

func TestLowerIndexMatrixAxis(t *testing.T) {
	// Lowering the first index of a rank-2 tensor with Minkowski metric negates
	// its first row.
	eta := MinkowskiMetric(2)
	m, _ := NewWithData([]int{2, 2}, []float64{1, 2, 3, 4})
	low, err := LowerIndex(m, eta, 0)
	if err != nil {
		t.Fatal(err)
	}
	want, _ := NewWithData([]int{2, 2}, []float64{-1, -2, 3, 4})
	if !low.Equal(want) {
		t.Fatalf("LowerIndex axis 0 = %v, want %v", low, want)
	}
}
