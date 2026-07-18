package ecc

import (
	"math/big"
	"testing"
)

func bi(n int64) *big.Int { return big.NewInt(n) }

func TestModArithmetic(t *testing.T) {
	p := bi(17)
	if got := Mod(bi(-1), p); got.Cmp(bi(16)) != 0 {
		t.Errorf("Mod(-1,17) = %s, want 16", got)
	}
	if got := ModAdd(bi(10), bi(12), p); got.Cmp(bi(5)) != 0 {
		t.Errorf("ModAdd(10,12,17) = %s, want 5", got)
	}
	if got := ModSub(bi(3), bi(10), p); got.Cmp(bi(10)) != 0 {
		t.Errorf("ModSub(3,10,17) = %s, want 10", got)
	}
	if got := ModMul(bi(6), bi(6), p); got.Cmp(bi(2)) != 0 {
		t.Errorf("ModMul(6,6,17) = %s, want 2", got)
	}
	if got := ModNeg(bi(5), p); got.Cmp(bi(12)) != 0 {
		t.Errorf("ModNeg(5,17) = %s, want 12", got)
	}
	if got := ModExp(bi(2), bi(4), p); got.Cmp(bi(16)) != 0 {
		t.Errorf("ModExp(2,4,17) = %s, want 16", got)
	}
	// Negative exponent: 2^-1 mod 17 = 9.
	if got := ModExp(bi(2), bi(-1), p); got.Cmp(bi(9)) != 0 {
		t.Errorf("ModExp(2,-1,17) = %s, want 9", got)
	}
}

func TestModInverseAndDiv(t *testing.T) {
	p := bi(17)
	inv, ok := ModInverse(bi(4), p)
	if !ok || inv.Cmp(bi(13)) != 0 {
		t.Errorf("ModInverse(4,17) = %s ok=%v, want 13 true", inv, ok)
	}
	if _, ok := ModInverse(bi(0), p); ok {
		t.Errorf("ModInverse(0,17) should not exist")
	}
	q, ok := ModDiv(bi(1), bi(2), p)
	if !ok || q.Cmp(bi(9)) != 0 {
		t.Errorf("ModDiv(1,2,17) = %s ok=%v, want 9 true", q, ok)
	}
}

func TestModSqrtAndLegendre(t *testing.T) {
	p := bi(17)
	// 2 is a QR mod 17 (6^2 = 36 = 2), canonical root is min(6,11) = 6.
	if Legendre(bi(2), p) != 1 {
		t.Errorf("Legendre(2,17) want 1")
	}
	r, ok := ModSqrt(bi(2), p)
	if !ok || r.Cmp(bi(6)) != 0 {
		t.Errorf("ModSqrt(2,17) = %s ok=%v, want 6 true", r, ok)
	}
	// verify r^2 = 2
	if ModMul(r, r, p).Cmp(bi(2)) != 0 {
		t.Errorf("ModSqrt root does not square back to 2")
	}
	// 3 is a non-residue mod 17.
	if Legendre(bi(3), p) != -1 {
		t.Errorf("Legendre(3,17) want -1")
	}
	if _, ok := ModSqrt(bi(3), p); ok {
		t.Errorf("ModSqrt(3,17) should have no root")
	}
	if IsQuadraticResidue(bi(3), p) {
		t.Errorf("3 should not be a QR mod 17")
	}
	if !IsQuadraticResidue(bi(2), p) {
		t.Errorf("2 should be a QR mod 17")
	}
	if Legendre(bi(17), p) != 0 {
		t.Errorf("Legendre(17,17) want 0")
	}
}
