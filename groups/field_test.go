package groups

import "testing"

func TestIsPrime(t *testing.T) {
	// Reference primes up to 101 (closed-form known list).
	primeList := []int{2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
		53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101}
	primes := map[int]bool{}
	for _, p := range primeList {
		primes[p] = true
	}
	for n := 0; n <= 101; n++ {
		got := IsPrime(n)
		if primes[n] != got {
			t.Errorf("IsPrime(%d)=%v want %v", n, got, primes[n])
		}
	}
	if IsPrime(1) || IsPrime(0) || IsPrime(-7) {
		t.Errorf("non-primes flagged prime")
	}
}

func TestGFArith(t *testing.T) {
	p := 7
	if GFAdd(5, 4, p) != 2 {
		t.Errorf("GFAdd")
	}
	if GFSub(2, 5, p) != 4 {
		t.Errorf("GFSub")
	}
	if GFMul(3, 5, p) != 1 {
		t.Errorf("GFMul 3*5 mod 7 want 1")
	}
	if GFNeg(3, p) != 4 {
		t.Errorf("GFNeg")
	}
	if GFPow(3, 6, p) != 1 { // Fermat
		t.Errorf("GFPow Fermat")
	}
}

func TestGFInvDiv(t *testing.T) {
	p := 13
	for a := 1; a < p; a++ {
		inv := GFInv(a, p)
		if GFMul(a, inv, p) != 1 {
			t.Errorf("GFInv(%d) mod %d wrong", a, p)
		}
	}
	if GFDiv(6, 3, p) != 2 {
		t.Errorf("GFDiv 6/3 want 2")
	}
	// division reverses multiplication
	if GFDiv(GFMul(9, 5, p), 5, p) != 9 {
		t.Errorf("GFDiv inverse of GFMul")
	}
}

func TestGFInvPanicsOnZero(t *testing.T) {
	defer func() {
		if recover() == nil {
			t.Errorf("GFInv(0) should panic")
		}
	}()
	GFInv(0, 7)
}

func TestGFPanicsNonPrime(t *testing.T) {
	defer func() {
		if recover() == nil {
			t.Errorf("GFMul with composite modulus should panic")
		}
	}()
	GFMul(2, 3, 8)
}
