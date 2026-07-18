package crypto

import (
	"math/big"
	"testing"
)

func bi(n int64) *big.Int { return big.NewInt(n) }

func TestModExp(t *testing.T) {
	cases := []struct {
		base, exp, mod, want int64
	}{
		{2, 10, 1000, 24}, // 1024 mod 1000
		{3, 4, 5, 1},      // 81 mod 5
		{7, 0, 13, 1},     // anything^0
		{5, 3, 13, 8},     // 125 mod 13
		{2, 100, 7, 2},    // 2^100 mod 7
		{10, 5, 6, 4},     // 100000 mod 6
		{123, 456, 789, 699},
	}
	for _, c := range cases {
		got := ModExp(bi(c.base), bi(c.exp), bi(c.mod))
		if got.Int64() != c.want {
			t.Errorf("ModExp(%d,%d,%d)=%d want %d", c.base, c.exp, c.mod, got.Int64(), c.want)
		}
		// Cross-check against big.Int.Exp.
		ref := new(big.Int).Exp(bi(c.base), bi(c.exp), bi(c.mod))
		if got.Cmp(ref) != 0 {
			t.Errorf("ModExp disagrees with big.Int.Exp for %v", c)
		}
	}
}

func TestModExpNegativeExponent(t *testing.T) {
	// 3^-1 mod 7 = 5, so 3^-2 = 25 mod 7 = 4.
	got := ModExp(bi(3), bi(-2), bi(7))
	if got.Int64() != 4 {
		t.Errorf("ModExp(3,-2,7)=%d want 4", got.Int64())
	}
}

func TestModInverse(t *testing.T) {
	cases := []struct{ a, m, want int64 }{
		{3, 7, 5},
		{3, 11, 4},
		{10, 17, 12},
		{1, 5, 1},
	}
	for _, c := range cases {
		got, err := ModInverse(bi(c.a), bi(c.m))
		if err != nil {
			t.Fatalf("ModInverse(%d,%d) err %v", c.a, c.m, err)
		}
		if got.Int64() != c.want {
			t.Errorf("ModInverse(%d,%d)=%d want %d", c.a, c.m, got.Int64(), c.want)
		}
	}
	if _, err := ModInverse(bi(2), bi(4)); err != ErrNotInvertible {
		t.Errorf("ModInverse(2,4) expected ErrNotInvertible, got %v", err)
	}
}

func TestGCDLCM(t *testing.T) {
	cases := []struct{ a, b, g, l int64 }{
		{12, 18, 6, 36},
		{17, 5, 1, 85},
		{0, 9, 9, 0},
		{-12, 18, 6, 36},
	}
	for _, c := range cases {
		if g := GCD(bi(c.a), bi(c.b)); g.Int64() != c.g {
			t.Errorf("GCD(%d,%d)=%d want %d", c.a, c.b, g.Int64(), c.g)
		}
		if l := LCM(bi(c.a), bi(c.b)); l.Int64() != c.l {
			t.Errorf("LCM(%d,%d)=%d want %d", c.a, c.b, l.Int64(), c.l)
		}
	}
}

func TestExtendedGCD(t *testing.T) {
	cases := []struct{ a, b int64 }{{240, 46}, {17, 5}, {99, 78}}
	for _, c := range cases {
		g, x, y := ExtendedGCD(bi(c.a), bi(c.b))
		// Verify a*x + b*y == g.
		lhs := new(big.Int).Add(new(big.Int).Mul(bi(c.a), x), new(big.Int).Mul(bi(c.b), y))
		if lhs.Cmp(g) != 0 {
			t.Errorf("ExtendedGCD(%d,%d): a*x+b*y=%v want g=%v", c.a, c.b, lhs, g)
		}
		if g.Cmp(GCD(bi(c.a), bi(c.b))) != 0 {
			t.Errorf("ExtendedGCD(%d,%d) g mismatch", c.a, c.b)
		}
	}
}

func TestModMulAddSub(t *testing.T) {
	if ModMul(bi(7), bi(8), bi(10)).Int64() != 6 {
		t.Error("ModMul(7,8,10) want 6")
	}
	if ModAdd(bi(7), bi(8), bi(10)).Int64() != 5 {
		t.Error("ModAdd(7,8,10) want 5")
	}
	if ModSub(bi(3), bi(8), bi(10)).Int64() != 5 {
		t.Error("ModSub(3,8,10) want 5")
	}
}

func TestCRT(t *testing.T) {
	// Classic Sunzi problem: x≡2(3), x≡3(5), x≡2(7) -> 23.
	got, err := CRT([]*big.Int{bi(2), bi(3), bi(2)}, []*big.Int{bi(3), bi(5), bi(7)})
	if err != nil {
		t.Fatal(err)
	}
	if got.Int64() != 23 {
		t.Errorf("CRT Sunzi=%d want 23", got.Int64())
	}
	// Pair: x≡1(4), x≡2(9) -> 29 (mod 36).
	p, err := CRTPair(bi(1), bi(4), bi(2), bi(9))
	if err != nil {
		t.Fatal(err)
	}
	if p.Int64() != 29 {
		t.Errorf("CRTPair=%d want 29", p.Int64())
	}
	// Non-coprime moduli should error.
	if _, err := CRTPair(bi(1), bi(4), bi(2), bi(6)); err == nil {
		t.Error("CRTPair with non-coprime moduli should error")
	}
}
