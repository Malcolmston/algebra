package settheory

import "testing"

func TestBellNumbers(t *testing.T) {
	// OEIS A000110.
	want := []int{1, 1, 2, 5, 15, 52, 203, 877, 4140, 21147, 115975}
	for n, w := range want {
		if got := BellNumber(n); got != w {
			t.Errorf("BellNumber(%d) = %d, want %d", n, got, w)
		}
	}
}

func TestStirlingSecondKnown(t *testing.T) {
	cases := []struct {
		n, k, want int
	}{
		{0, 0, 1},
		{1, 1, 1},
		{4, 2, 7},
		{5, 3, 25},
		{6, 3, 90},
		{5, 5, 1},
		{5, 0, 0},
		{3, 5, 0},
		{10, 4, 34105},
	}
	for _, c := range cases {
		if got := StirlingSecond(c.n, c.k); got != c.want {
			t.Errorf("StirlingSecond(%d,%d) = %d, want %d", c.n, c.k, got, c.want)
		}
	}
}

func TestBellEqualsSumOfStirling(t *testing.T) {
	// Bell(n) = sum_{k=0}^{n} S(n,k).
	for n := 0; n <= 9; n++ {
		sum := 0
		for k := 0; k <= n; k++ {
			sum += StirlingSecond(n, k)
		}
		if sum != BellNumber(n) {
			t.Errorf("sum S(%d,k) = %d, Bell = %d", n, sum, BellNumber(n))
		}
	}
}

func TestPartitionRoundTrip(t *testing.T) {
	set := NewIntSet(1, 2, 3, 4, 5)
	// Classes {1,3,5} (odd) and {2,4} (even) via mod-2 equivalence.
	r := make(Relation)
	for a := range set {
		for b := range set {
			if a%2 == b%2 {
				r.Add(a, b)
			}
		}
	}
	p := EquivalenceClasses(r, set)
	if !p.IsValidPartitionOf(set) {
		t.Fatalf("not a valid partition")
	}
	if p.Len() != 2 {
		t.Errorf("expected 2 blocks, got %d", p.Len())
	}
	// Reconstruct the relation from the partition; must match the original.
	if !RelationFromPartition(p).Equal(r) {
		t.Errorf("relation reconstructed from partition differs")
	}
}

func TestPartitionRefines(t *testing.T) {
	fine := Partition{NewIntSet(1), NewIntSet(2), NewIntSet(3)}
	coarse := Partition{NewIntSet(1, 2), NewIntSet(3)}
	if !fine.Refines(coarse) {
		t.Errorf("singleton partition should refine any partition")
	}
	if coarse.Refines(fine) {
		t.Errorf("coarse partition should not refine the finer one")
	}
}
