package groups

import "testing"

func TestModBasics(t *testing.T) {
	if Mod(-1, 7) != 6 {
		t.Errorf("Mod(-1,7) want 6")
	}
	if ModAdd(5, 4, 7) != 2 {
		t.Errorf("ModAdd want 2")
	}
	if ModSub(2, 5, 7) != 4 {
		t.Errorf("ModSub want 4")
	}
	if ModMul(6, 6, 7) != 1 {
		t.Errorf("ModMul 6*6 mod 7 want 1")
	}
	if ModNeg(3, 7) != 4 {
		t.Errorf("ModNeg want 4")
	}
}

func TestModPow(t *testing.T) {
	cases := []struct{ b, e, n, want int }{
		{2, 10, 1000, 24},
		{3, 0, 7, 1},
		{7, 4, 13, 9}, // 7^4 = 2401 = 184*13 + 9
		{5, 3, 13, 8}, // 125 mod 13 = 8
		{2, 5, 1, 0},
	}
	for _, c := range cases {
		if got := ModPow(c.b, c.e, c.n); got != c.want {
			t.Errorf("ModPow(%d,%d,%d)=%d want %d", c.b, c.e, c.n, got, c.want)
		}
	}
}

func TestModInverse(t *testing.T) {
	if inv, ok := ModInverse(3, 11); !ok || inv != 4 {
		t.Errorf("ModInverse(3,11)=%d,%v want 4,true", inv, ok)
	}
	if _, ok := ModInverse(6, 9); ok {
		t.Errorf("ModInverse(6,9) should fail (not coprime)")
	}
	// Verify a*inv == 1 for all units mod 17.
	for a := 1; a < 17; a++ {
		inv, ok := ModInverse(a, 17)
		if !ok || ModMul(a, inv, 17) != 1 {
			t.Errorf("inverse check failed for %d mod 17", a)
		}
	}
}

func TestElementOrderZn(t *testing.T) {
	// Additive order of a in Z/nZ is n/gcd(a,n).
	cases := []struct{ a, n, want int }{
		{1, 12, 12},
		{3, 12, 4},
		{4, 12, 3},
		{6, 12, 2},
		{0, 12, 1},
	}
	for _, c := range cases {
		if got := ElementOrderZn(c.a, c.n); got != c.want {
			t.Errorf("ElementOrderZn(%d,%d)=%d want %d", c.a, c.n, got, c.want)
		}
	}
}

func TestMultiplicativeOrder(t *testing.T) {
	cases := []struct {
		a, n, want int
		ok         bool
	}{
		{2, 7, 3, true},  // 2,4,1
		{3, 7, 6, true},  // primitive root
		{2, 15, 4, true}, // 2,4,8,1
		{6, 9, 0, false}, // not coprime
		{1, 5, 1, true},
	}
	for _, c := range cases {
		got, ok := MultiplicativeOrder(c.a, c.n)
		if got != c.want || ok != c.ok {
			t.Errorf("MultiplicativeOrder(%d,%d)=%d,%v want %d,%v", c.a, c.n, got, ok, c.want, c.ok)
		}
	}
}

func TestEulerTotient(t *testing.T) {
	cases := []struct{ n, want int }{
		{1, 1}, {2, 1}, {6, 2}, {9, 6}, {10, 4}, {12, 4}, {36, 12}, {97, 96},
	}
	for _, c := range cases {
		if got := EulerTotient(c.n); got != c.want {
			t.Errorf("EulerTotient(%d)=%d want %d", c.n, got, c.want)
		}
	}
}

func TestUnitsModN(t *testing.T) {
	u := UnitsModN(12)
	want := []int{1, 5, 7, 11}
	if len(u) != len(want) {
		t.Fatalf("UnitsModN(12)=%v want %v", u, want)
	}
	for i := range want {
		if u[i] != want[i] {
			t.Fatalf("UnitsModN(12)=%v want %v", u, want)
		}
	}
	if len(UnitsModN(12)) != EulerTotient(12) {
		t.Errorf("units count mismatch totient")
	}
}

func TestPrimitiveRoots(t *testing.T) {
	// Primitive roots mod 7 are 3 and 5.
	pr := PrimitiveRoots(7)
	want := map[int]bool{3: true, 5: true}
	if len(pr) != 2 {
		t.Fatalf("PrimitiveRoots(7)=%v", pr)
	}
	for _, g := range pr {
		if !want[g] {
			t.Errorf("unexpected primitive root %d", g)
		}
	}
	// 8 has no primitive root.
	if len(PrimitiveRoots(8)) != 0 {
		t.Errorf("PrimitiveRoots(8) should be empty")
	}
	if !IsPrimitiveRoot(3, 7) {
		t.Errorf("3 is a primitive root mod 7")
	}
	if IsPrimitiveRoot(2, 7) {
		t.Errorf("2 is not a primitive root mod 7")
	}
}
