package groups

import "testing"

func groupsIntSliceEq(a, b []int) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

func TestPolyModArith(t *testing.T) {
	p := 5
	a := []int{1, 2, 3} // 3x^2+2x+1
	b := []int{4, 4}    // 4x+4
	if !groupsIntSliceEq(PolyModAdd(a, b, p), []int{0, 1, 3}) {
		t.Errorf("PolyModAdd=%v", PolyModAdd(a, b, p))
	}
	// (x+1)(x+1) over GF(5) = x^2+2x+1
	sq := PolyModMul([]int{1, 1}, []int{1, 1}, p)
	if !groupsIntSliceEq(sq, []int{1, 2, 1}) {
		t.Errorf("PolyModMul=%v want [1 2 1]", sq)
	}
	if PolyModDegree(a, p) != 2 {
		t.Errorf("PolyModDegree want 2")
	}
}

func TestPolyModDivMod(t *testing.T) {
	p := 7
	// x^2 - 1 over GF(7): [6,0,1]; divide by x-1 = [6,1]
	a := []int{6, 0, 1}
	b := []int{6, 1}
	q, r := PolyModDivMod(a, b, p)
	// reconstruct
	recon := PolyModAdd(PolyModMul(q, b, p), r, p)
	if !groupsIntSliceEq(recon, PolyModAdd(a, []int{}, p)) {
		t.Errorf("a != q*b+r over GF(%d): %v", p, recon)
	}
	if len(r) != 0 {
		t.Errorf("x^2-1 divisible by x-1, remainder=%v", r)
	}
}

func TestPolyModGCD(t *testing.T) {
	p := 5
	// gcd((x-1)^2, x^2-1) = x-1, monic over GF(5): [4,1]
	a := PolyModMul([]int{4, 1}, []int{4, 1}, p) // (x-1)^2
	b := []int{4, 0, 1}                          // x^2-1
	g := PolyModGCD(a, b, p)
	if !groupsIntSliceEq(g, []int{4, 1}) {
		t.Errorf("PolyModGCD=%v want [4 1] (x-1)", g)
	}
	// gcd is monic
	if g[len(g)-1] != 1 {
		t.Errorf("gcd not monic: %v", g)
	}
}

func TestPolyModEval(t *testing.T) {
	p := 7
	a := []int{1, 2, 3} // 3x^2+2x+1
	// at x=2: 12+4+1=17 mod 7 = 3
	if PolyModEval(a, 2, p) != 3 {
		t.Errorf("PolyModEval want 3 got %d", PolyModEval(a, 2, p))
	}
	// root check: (x-3) evaluated at 3 over GF(7) is 0
	if PolyModEval([]int{4, 1}, 3, p) != 0 {
		t.Errorf("x-3 at 3 should be 0")
	}
}

func TestPolyModMonic(t *testing.T) {
	p := 7
	m := PolyModMonic([]int{4, 2}, p) // 2x+4 -> monic
	if m[len(m)-1] != 1 {
		t.Errorf("not monic: %v", m)
	}
}
