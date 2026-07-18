package groups

import "testing"

func TestComposeAndInverse(t *testing.T) {
	p := Perm{1, 2, 0}  // cycle (0 1 2)
	q := Perm{0, 2, 1}  // transposition (1 2)
	pq := Compose(p, q) // apply q then p
	// q: 0->0,1->2,2->1 ; then p: 0->1,2->0,1->2 => 0->1,1->0,2->2
	want := Perm{1, 0, 2}
	if !pq.Equal(want) {
		t.Errorf("Compose=%v want %v", pq, want)
	}
	if !Compose(p, p.Inverse()).IsIdentity() {
		t.Errorf("p * p^-1 != identity")
	}
}

func TestPermOrderSign(t *testing.T) {
	cases := []struct {
		p          Perm
		order, sig int
	}{
		{Perm{0, 1, 2}, 1, 1},       // identity
		{Perm{1, 0, 2}, 2, -1},      // transposition
		{Perm{1, 2, 0}, 3, 1},       // 3-cycle even
		{Perm{1, 2, 3, 0}, 4, -1},   // 4-cycle odd
		{Perm{1, 0, 3, 2}, 2, 1},    // (01)(23) even
		{Perm{2, 3, 4, 0, 1}, 5, 1}, // 5-cycle even
	}
	for _, c := range cases {
		if got := c.p.Order(); got != c.order {
			t.Errorf("Order(%v)=%d want %d", c.p, got, c.order)
		}
		if got := c.p.Sign(); got != c.sig {
			t.Errorf("Sign(%v)=%d want %d", c.p, got, c.sig)
		}
	}
}

func TestSignMatchesInversionParity(t *testing.T) {
	for _, p := range SymmetricGroup(5) {
		inv := p.NumInversions()
		sig := 1
		if inv%2 == 1 {
			sig = -1
		}
		if sig != p.Sign() {
			t.Errorf("sign/inversion mismatch for %v", p)
		}
	}
}

func TestCyclesAndType(t *testing.T) {
	p := Perm{1, 2, 0, 4, 3} // (0 1 2)(3 4)
	cyc := p.Cycles()
	if len(cyc) != 2 {
		t.Fatalf("expected 2 cycles got %v", cyc)
	}
	ct := p.CycleType()
	if len(ct) != 2 || ct[0] != 3 || ct[1] != 2 {
		t.Errorf("CycleType=%v want [3 2]", ct)
	}
	// round trip through PermFromCycles
	rebuilt := PermFromCycles(5, [][]int{{0, 1, 2}, {3, 4}})
	if !rebuilt.Equal(p) {
		t.Errorf("PermFromCycles=%v want %v", rebuilt, p)
	}
}

func TestPermPow(t *testing.T) {
	p := Perm{1, 2, 3, 0} // 4-cycle
	if !p.Pow(4).IsIdentity() {
		t.Errorf("p^4 should be identity")
	}
	if !p.Pow(-1).Equal(p.Inverse()) {
		t.Errorf("p^-1 mismatch")
	}
	if !p.Pow(0).IsIdentity() {
		t.Errorf("p^0 should be identity")
	}
}

func TestTranspositionAndFactorial(t *testing.T) {
	tr := Transposition(4, 1, 3)
	if tr[1] != 3 || tr[3] != 1 || tr[0] != 0 || tr[2] != 2 {
		t.Errorf("Transposition wrong: %v", tr)
	}
	if tr.Sign() != -1 {
		t.Errorf("transposition sign want -1")
	}
	facts := []int{1, 1, 2, 6, 24, 120, 720}
	for n, w := range facts {
		if Factorial(n) != w {
			t.Errorf("Factorial(%d)=%d want %d", n, Factorial(n), w)
		}
	}
}

func TestSymmetricAlternatingGroup(t *testing.T) {
	for n := 0; n <= 5; n++ {
		sg := SymmetricGroup(n)
		if len(sg) != SymmetricGroupOrder(n) {
			t.Errorf("|S_%d|=%d want %d", n, len(sg), SymmetricGroupOrder(n))
		}
		// validity + distinctness
		seen := map[string]bool{}
		for _, p := range sg {
			if !IsPermutation(p) {
				t.Errorf("invalid perm in S_%d: %v", n, p)
			}
			k := groupsPermKey(p)
			if seen[k] {
				t.Errorf("duplicate perm in S_%d", n)
			}
			seen[k] = true
		}
		ag := AlternatingGroup(n)
		if len(ag) != AlternatingGroupOrder(n) {
			t.Errorf("|A_%d|=%d want %d", n, len(ag), AlternatingGroupOrder(n))
		}
		for _, p := range ag {
			if p.Sign() != 1 {
				t.Errorf("odd perm in A_%d", n)
			}
		}
	}
}
