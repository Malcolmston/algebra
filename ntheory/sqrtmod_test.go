package ntheory

import (
	"math/big"
	"reflect"
	"testing"
)

func TestSqrtMod(t *testing.T) {
	tests := []struct {
		a, p   int64
		wantX  int64
		wantOk bool
	}{
		// p ≡ 3 (mod 4) closed-form branch.
		{a: 0, p: 7, wantX: 0, wantOk: true},
		{a: 1, p: 7, wantX: 1, wantOk: true},
		{a: 2, p: 7, wantX: 3, wantOk: true},
		{a: 4, p: 7, wantX: 2, wantOk: true},
		{a: 9, p: 7, wantX: 3, wantOk: true}, // 9 ≡ 2 (mod 7)
		{a: 3, p: 7, wantOk: false},          // non-residue
		{a: 5, p: 7, wantOk: false},
		// p == 2 special case.
		{a: 0, p: 2, wantX: 0, wantOk: true},
		{a: 1, p: 2, wantX: 1, wantOk: true},
		{a: 3, p: 2, wantX: 1, wantOk: true},
		// p ≡ 1 (mod 4) Tonelli-Shanks branch.
		{a: 4, p: 13, wantX: 2, wantOk: true},
		{a: 10, p: 13, wantX: 6, wantOk: true},
		{a: 3, p: 13, wantX: 4, wantOk: true},
		{a: 2, p: 13, wantOk: false},
		{a: 2, p: 17, wantX: 6, wantOk: true},
		{a: 8, p: 17, wantX: 5, wantOk: true},
		{a: 3, p: 17, wantOk: false},
	}
	for _, tc := range tests {
		x, ok := SqrtMod(tc.a, tc.p)
		if ok != tc.wantOk {
			t.Errorf("SqrtMod(%d,%d) ok=%v, want %v", tc.a, tc.p, ok, tc.wantOk)
			continue
		}
		if !ok {
			continue
		}
		if x != tc.wantX {
			t.Errorf("SqrtMod(%d,%d)=%d, want %d", tc.a, tc.p, x, tc.wantX)
		}
		// Structural check: x is really a root and lies in [0,p).
		if x < 0 || x >= tc.p {
			t.Errorf("SqrtMod(%d,%d)=%d out of range", tc.a, tc.p, x)
		}
		if mulMod(x, x, tc.p) != ((tc.a%tc.p)+tc.p)%tc.p {
			t.Errorf("SqrtMod(%d,%d)=%d not a valid root", tc.a, tc.p, x)
		}
	}
}

func TestSqrtModExhaustive(t *testing.T) {
	// For several primes, every residue must round-trip and every non-residue
	// must be rejected.
	for _, p := range []int64{3, 5, 7, 11, 13, 17, 19, 23, 29, 31} {
		residues := map[int64]bool{}
		for x := int64(0); x < p; x++ {
			residues[mulMod(x, x, p)] = true
		}
		for a := int64(0); a < p; a++ {
			x, ok := SqrtMod(a, p)
			if ok != residues[a] {
				t.Fatalf("SqrtMod(%d,%d) ok=%v, want %v", a, p, ok, residues[a])
			}
			if ok && mulMod(x, x, p) != a {
				t.Fatalf("SqrtMod(%d,%d)=%d not a root", a, p, x)
			}
		}
	}
}

func TestSqrtModBig(t *testing.T) {
	// p ≡ 3 (mod 4): 2^61 - 1 is prime.
	p3 := new(big.Int).Sub(new(big.Int).Lsh(big.NewInt(1), 61), big.NewInt(1))
	// p ≡ 1 (mod 4): a known prime.
	p1 := big.NewInt(1000000009)

	for _, p := range []*big.Int{big.NewInt(7), big.NewInt(13), big.NewInt(17), p3, p1} {
		for _, av := range []int64{0, 1, 2, 3, 4, 5, 6, 10} {
			a := big.NewInt(av)
			r, ok := SqrtModBig(a, p)
			if !ok {
				continue
			}
			sq := new(big.Int).Mod(new(big.Int).Mul(r, r), p)
			want := new(big.Int).Mod(a, p)
			if sq.Cmp(want) != 0 {
				t.Errorf("SqrtModBig(%v,%v)=%v not a root (got %v)", a, p, r, sq)
			}
			if r.Sign() < 0 || r.Cmp(p) >= 0 {
				t.Errorf("SqrtModBig(%v,%v)=%v out of range", a, p, r)
			}
		}
	}
}

func TestSqrtModBigNonResidue(t *testing.T) {
	// 2 is a non-residue modulo 5; 3 is a non-residue modulo 7.
	if r, ok := SqrtModBig(big.NewInt(2), big.NewInt(5)); ok {
		t.Errorf("SqrtModBig(2,5)=%v, want non-residue", r)
	}
	if r, ok := SqrtModBig(big.NewInt(3), big.NewInt(7)); ok {
		t.Errorf("SqrtModBig(3,7)=%v, want non-residue", r)
	}
}

func TestSqrtModBigAgreesWithInt(t *testing.T) {
	for _, p := range []int64{7, 13, 17, 23, 29} {
		for a := int64(0); a < p; a++ {
			xi, oki := SqrtMod(a, p)
			rb, okb := SqrtModBig(big.NewInt(a), big.NewInt(p))
			if oki != okb {
				t.Fatalf("residue disagreement for a=%d p=%d: int=%v big=%v", a, p, oki, okb)
			}
			if oki && rb.Int64() != xi {
				t.Fatalf("root disagreement a=%d p=%d: int=%d big=%d", a, p, xi, rb.Int64())
			}
		}
	}
}

func TestSqrtModPrimePower(t *testing.T) {
	tests := []struct {
		a, p   int64
		k      int
		wantX  int64
		wantOk bool
	}{
		{a: 2, p: 7, k: 2, wantX: 10, wantOk: true}, // 10^2 = 100 ≡ 2 (mod 49)
		{a: 1, p: 3, k: 3, wantX: 1, wantOk: true},  // mod 27
		{a: 0, p: 5, k: 2, wantX: 0, wantOk: true},  // mod 25
		{a: 3, p: 7, k: 2, wantOk: false},           // non-residue mod 7 -> mod 49
		{a: 1, p: 2, k: 3, wantX: 1, wantOk: true},  // mod 8, roots {1,3,5,7}
		{a: 4, p: 2, k: 3, wantX: 2, wantOk: true},  // mod 8, roots {2,6}
		{a: 3, p: 2, k: 3, wantOk: false},           // 3 is not a square mod 8
		{a: 9, p: 3, k: 2, wantX: 0, wantOk: true},  // 9 ≡ 0 (mod 9)
	}
	for _, tc := range tests {
		x, ok := SqrtModPrimePower(tc.a, tc.p, tc.k)
		if ok != tc.wantOk {
			t.Errorf("SqrtModPrimePower(%d,%d,%d) ok=%v, want %v", tc.a, tc.p, tc.k, ok, tc.wantOk)
			continue
		}
		if !ok {
			continue
		}
		if x != tc.wantX {
			t.Errorf("SqrtModPrimePower(%d,%d,%d)=%d, want %d", tc.a, tc.p, tc.k, x, tc.wantX)
		}
		q := ntheoryIntPow(tc.p, tc.k)
		if mulMod(x, x, q) != ((tc.a%q)+q)%q {
			t.Errorf("SqrtModPrimePower(%d,%d,%d)=%d not a valid root mod %d", tc.a, tc.p, tc.k, x, q)
		}
	}
}

func TestSqrtModPrimePowerExhaustive(t *testing.T) {
	type pk struct {
		p int64
		k int
	}
	for _, c := range []pk{{3, 3}, {5, 2}, {7, 2}, {2, 4}, {2, 5}, {11, 2}} {
		q := ntheoryIntPow(c.p, c.k)
		residues := map[int64]bool{}
		for x := int64(0); x < q; x++ {
			residues[mulMod(x, x, q)] = true
		}
		for a := int64(0); a < q; a++ {
			x, ok := SqrtModPrimePower(a, c.p, c.k)
			if ok != residues[a] {
				t.Fatalf("SqrtModPrimePower(%d,%d,%d) ok=%v, want %v", a, c.p, c.k, ok, residues[a])
			}
			if ok && mulMod(x, x, q) != a {
				t.Fatalf("SqrtModPrimePower(%d,%d,%d)=%d not a root", a, c.p, c.k, x)
			}
		}
	}
}

func TestAllSqrtModComposite(t *testing.T) {
	tests := []struct {
		a, m uint64
		want []uint64
	}{
		{a: 4, m: 15, want: []uint64{2, 7, 8, 13}},
		{a: 1, m: 12, want: []uint64{1, 5, 7, 11}},
		{a: 0, m: 12, want: []uint64{0, 6}},
		{a: 2, m: 7, want: []uint64{3, 4}},
		{a: 1, m: 8, want: []uint64{1, 3, 5, 7}},
		{a: 4, m: 15, want: []uint64{2, 7, 8, 13}},
		{a: 6, m: 10, want: []uint64{4, 6}},
		{a: 1, m: 1, want: []uint64{0}},
	}
	for _, tc := range tests {
		got := AllSqrtModComposite(tc.a, tc.m)
		if !reflect.DeepEqual(got, tc.want) {
			t.Errorf("AllSqrtModComposite(%d,%d)=%v, want %v", tc.a, tc.m, got, tc.want)
		}
	}
}

func TestAllSqrtModCompositeNonResidue(t *testing.T) {
	// 3 is not a square modulo 5, hence not modulo 15.
	if got := AllSqrtModComposite(3, 15); len(got) != 0 {
		t.Errorf("AllSqrtModComposite(3,15)=%v, want empty", got)
	}
	// 7 is not a square modulo 12.
	if got := AllSqrtModComposite(7, 12); len(got) != 0 {
		t.Errorf("AllSqrtModComposite(7,12)=%v, want empty", got)
	}
}

func TestAllSqrtModCompositeExhaustive(t *testing.T) {
	for _, m := range []uint64{6, 8, 12, 15, 16, 24, 30, 36, 45} {
		want := map[uint64][]uint64{}
		for x := uint64(0); x < m; x++ {
			r := (x * x) % m
			want[r] = append(want[r], x)
		}
		for a := uint64(0); a < m; a++ {
			got := AllSqrtModComposite(a, m)
			exp := ntheoryDedupSortU64(append([]uint64(nil), want[a]...))
			if len(got) == 0 {
				got = nil
			}
			if len(exp) == 0 {
				exp = nil
			}
			if !reflect.DeepEqual(got, exp) {
				t.Fatalf("AllSqrtModComposite(%d,%d)=%v, want %v", a, m, got, exp)
			}
		}
	}
}

func BenchmarkSqrtModCongruent3Mod4(b *testing.B) {
	// Exercises the closed-form fast path (p ≡ 3 mod 4).
	const p = int64(1000000007) // ≡ 3 (mod 4)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		SqrtMod(int64(i)%p*int64(i)%p%p, p)
	}
}

func BenchmarkSqrtModCongruent1Mod4(b *testing.B) {
	// Exercises the full Tonelli-Shanks path (p ≡ 1 mod 4).
	const p = int64(1000000009) // ≡ 1 (mod 4)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		SqrtMod(int64(i)%p*int64(i)%p%p, p)
	}
}

func BenchmarkSqrtModBig(b *testing.B) {
	p := new(big.Int).Sub(new(big.Int).Lsh(big.NewInt(1), 127), big.NewInt(1)) // 2^127-1 (prime)
	a := big.NewInt(2)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		SqrtModBig(a, p)
	}
}

func BenchmarkAllSqrtModComposite(b *testing.B) {
	const m = uint64(360360) // 2^3 * 3^2 * 5 * 7 * 11 * 13
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		AllSqrtModComposite(uint64(i)%m, m)
	}
}
