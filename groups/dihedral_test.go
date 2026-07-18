package groups

import "testing"

func TestDihedralBasics(t *testing.T) {
	n := 4
	if DihedralOrder(n) != 8 {
		t.Errorf("DihedralOrder(4) want 8")
	}
	if len(DihedralGroup(n)) != 8 {
		t.Errorf("DihedralGroup(4) size")
	}
	e := DihedralIdentity()
	r := DihedralRotation(n, 1)
	s := DihedralReflection(n, 0)
	// r has order n
	if DihedralElementOrder(n, r) != 4 {
		t.Errorf("rotation order want 4")
	}
	// s has order 2
	if DihedralElementOrder(n, s) != 2 {
		t.Errorf("reflection order want 2")
	}
	if DihedralElementOrder(n, e) != 1 {
		t.Errorf("identity order want 1")
	}
}

func TestDihedralRelations(t *testing.T) {
	n := 5
	r := DihedralRotation(n, 1)
	s := DihedralReflection(n, 0)
	// s*r = r^-1*s  (i.e. srs = r^-1)
	lhs := DihedralCompose(n, DihedralCompose(n, s, r), s)
	rhs := DihedralInverse(n, r)
	if lhs != rhs {
		t.Errorf("srs^-1 = r^-1 relation failed: %v vs %v", lhs, rhs)
	}
	// s^2 = e
	if DihedralCompose(n, s, s) != DihedralIdentity() {
		t.Errorf("s^2 != e")
	}
	// r^n = e
	rn := DihedralIdentity()
	for i := 0; i < n; i++ {
		rn = DihedralCompose(n, rn, r)
	}
	if rn != DihedralIdentity() {
		t.Errorf("r^n != e")
	}
}

func TestDihedralInverseAllElements(t *testing.T) {
	n := 6
	for _, a := range DihedralGroup(n) {
		inv := DihedralInverse(n, a)
		if DihedralCompose(n, a, inv) != DihedralIdentity() ||
			DihedralCompose(n, inv, a) != DihedralIdentity() {
			t.Errorf("inverse failed for %v", a)
		}
	}
}

func TestDihedralIsGroup(t *testing.T) {
	// Build a FiniteGroup from D_4 and check axioms.
	n := 4
	elems := DihedralGroup(n)
	index := map[DihedralElement]int{}
	for i, e := range elems {
		index[e] = i
	}
	g := NewFiniteGroup(len(elems), index[DihedralIdentity()], func(a, b int) int {
		return index[DihedralCompose(n, elems[a], elems[b])]
	})
	if !g.IsValid() {
		t.Errorf("D_4 fails group axioms")
	}
	if g.IsAbelian() {
		t.Errorf("D_4 should be non-abelian")
	}
}
