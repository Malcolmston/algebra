package ntheory

import (
	"math/big"
	"testing"
)

// ntheoryCheckFactorU64 multiplies the prime powers back out, confirms each base
// is prime, and checks the product reconstructs n.
func ntheoryCheckFactorU64(t *testing.T, n uint64, factors map[uint64]int) {
	t.Helper()
	product := uint64(1)
	for p, e := range factors {
		if !IsPrimeU64(p) {
			t.Errorf("factor %d of %d is not prime", p, n)
		}
		if e < 1 {
			t.Errorf("factor %d of %d has non-positive exponent %d", p, n, e)
		}
		for i := 0; i < e; i++ {
			product *= p
		}
	}
	if n >= 2 && product != n {
		t.Errorf("FactorizeU64(%d) product = %d, want %d", n, product, n)
	}
}

func TestFactorizeU64(t *testing.T) {
	tests := []struct {
		n    uint64
		want map[uint64]int
	}{
		{0, map[uint64]int{}},
		{1, map[uint64]int{}},
		{2, map[uint64]int{2: 1}},
		{3, map[uint64]int{3: 1}},
		{12, map[uint64]int{2: 2, 3: 1}},
		{360, map[uint64]int{2: 3, 3: 2, 5: 1}},
		{1000000, map[uint64]int{2: 6, 5: 6}},
		{9999991, map[uint64]int{9999991: 1}}, // prime
		{600851475143, map[uint64]int{71: 1, 839: 1, 1471: 1, 6857: 1}},
		{1 << 20, map[uint64]int{2: 20}},
	}
	for _, tt := range tests {
		got := FactorizeU64(tt.n)
		if got == nil {
			t.Errorf("FactorizeU64(%d) returned nil map", tt.n)
			continue
		}
		if len(got) != len(tt.want) {
			t.Errorf("FactorizeU64(%d) = %v, want %v", tt.n, got, tt.want)
			continue
		}
		for p, e := range tt.want {
			if got[p] != e {
				t.Errorf("FactorizeU64(%d)[%d] = %d, want %d", tt.n, p, got[p], e)
			}
		}
	}
}

func TestFactorizeU64HardSemiprime(t *testing.T) {
	// Product of two large primes near 2^32; naive trial division is very slow.
	const p, q = uint64(4294967291), uint64(4294967279)
	n := p * q
	factors := FactorizeU64(n)
	ntheoryCheckFactorU64(t, n, factors)
	if factors[p] != 1 || factors[q] != 1 {
		t.Errorf("FactorizeU64(%d) = %v, want %d and %d each once", n, factors, p, q)
	}
}

func TestFactorizeU64Deterministic(t *testing.T) {
	// The result must be reproducible across calls (no randomness).
	const n = uint64(4294967291) * uint64(4294967279)
	a := FactorizeU64(n)
	b := FactorizeU64(n)
	if len(a) != len(b) {
		t.Fatalf("FactorizeU64 not deterministic: %v vs %v", a, b)
	}
	for p, e := range a {
		if b[p] != e {
			t.Fatalf("FactorizeU64 not deterministic at %d: %d vs %d", p, e, b[p])
		}
	}
}

func TestFactorizeU64Sweep(t *testing.T) {
	// Deterministic pseudo-random sweep validated by reconstruction.
	x := uint64(0xdeadbeef)
	for i := 0; i < 400; i++ {
		x = x*6364136223846793005 + 1442695040888963407
		n := x % 1000000000000
		if n < 2 {
			continue
		}
		ntheoryCheckFactorU64(t, n, FactorizeU64(n))
	}
}

func TestPollardBrentU64(t *testing.T) {
	// Composites: expect a proper divisor.
	for _, n := range []uint64{15, 21, 8051, 600851475143, uint64(4294967291) * uint64(4294967279)} {
		d := PollardBrentU64(n)
		if d <= 1 || d >= n || n%d != 0 {
			t.Errorf("PollardBrentU64(%d) = %d is not a proper divisor", n, d)
		}
	}
	// Primes: expect n itself.
	for _, n := range []uint64{2, 3, 17, 9999991, 4294967291} {
		if d := PollardBrentU64(n); d != n {
			t.Errorf("PollardBrentU64(%d) = %d, want %d (prime)", n, d, n)
		}
	}
	// Even composite: expect 2.
	if d := PollardBrentU64(1000000); d != 2 {
		t.Errorf("PollardBrentU64(1000000) = %d, want 2", d)
	}
}

func TestPollardBrentU64Deterministic(t *testing.T) {
	const n = uint64(4294967291) * uint64(4294967279)
	if PollardBrentU64(n) != PollardBrentU64(n) {
		t.Errorf("PollardBrentU64(%d) is not deterministic", n)
	}
}

func ntheoryBigFromString(t *testing.T, s string) *big.Int {
	t.Helper()
	v, ok := new(big.Int).SetString(s, 10)
	if !ok {
		t.Fatalf("bad big.Int literal %q", s)
	}
	return v
}

// ntheoryCheckFactorBig reconstructs |n| from its factorization and confirms
// each base is (probably) prime.
func ntheoryCheckFactorBig(t *testing.T, n *big.Int, factors []PrimePowerBig) {
	t.Helper()
	product := big.NewInt(1)
	for _, pp := range factors {
		if !IsProbablePrimeBig(pp.Prime) {
			t.Errorf("factor %s of %s is not prime", pp.Prime, n)
		}
		if pp.Exponent < 1 {
			t.Errorf("factor %s of %s has non-positive exponent %d", pp.Prime, n, pp.Exponent)
		}
		term := new(big.Int).Exp(pp.Prime, big.NewInt(int64(pp.Exponent)), nil)
		product.Mul(product, term)
	}
	want := new(big.Int).Abs(n)
	if want.Cmp(big.NewInt(2)) >= 0 && product.Cmp(want) != 0 {
		t.Errorf("FactorizeBig(%s) product = %s, want %s", n, product, want)
	}
}

func TestFactorizeBig(t *testing.T) {
	tests := []struct {
		n    string
		want []PrimePowerBig
	}{
		{"0", nil},
		{"1", nil},
		{"-1", nil},
		{"2", []PrimePowerBig{{big.NewInt(2), 1}}},
		{"360", []PrimePowerBig{{big.NewInt(2), 3}, {big.NewInt(3), 2}, {big.NewInt(5), 1}}},
		{"-360", []PrimePowerBig{{big.NewInt(2), 3}, {big.NewInt(3), 2}, {big.NewInt(5), 1}}},
		{"1000000", []PrimePowerBig{{big.NewInt(2), 6}, {big.NewInt(5), 6}}},
	}
	for _, tt := range tests {
		n := ntheoryBigFromString(t, tt.n)
		got := FactorizeBig(n)
		if len(got) != len(tt.want) {
			t.Errorf("FactorizeBig(%s) = %v, want %v", tt.n, got, tt.want)
			continue
		}
		for i := range tt.want {
			if got[i].Prime.Cmp(tt.want[i].Prime) != 0 || got[i].Exponent != tt.want[i].Exponent {
				t.Errorf("FactorizeBig(%s)[%d] = {%s,%d}, want {%s,%d}", tt.n, i,
					got[i].Prime, got[i].Exponent, tt.want[i].Prime, tt.want[i].Exponent)
			}
		}
	}
}

func TestFactorizeBigLarge(t *testing.T) {
	// A semiprime beyond the uint64 range. Pollard rho costs about
	// sqrt(smallest factor), so we keep one factor modest (2^31-1, a Mersenne
	// prime) and the other a 20-digit prime; their product still exceeds 2^64,
	// exercising the big.Int path while remaining tractable for rho.
	p := ntheoryBigFromString(t, "2147483647")
	q := ntheoryBigFromString(t, "10000000000000000051")
	n := new(big.Int).Mul(p, q)
	factors := FactorizeBig(n)
	ntheoryCheckFactorBig(t, n, factors)
	if len(factors) != 2 ||
		factors[0].Prime.Cmp(p) != 0 || factors[0].Exponent != 1 ||
		factors[1].Prime.Cmp(q) != 0 || factors[1].Exponent != 1 {
		t.Errorf("FactorizeBig(%s) = %v, want %s and %s each once", n, factors, p, q)
	}
}

func TestFactorizeBigPrimePower(t *testing.T) {
	// 7^13, all factors above the small-prime trial cutoff after one strip.
	n := new(big.Int).Exp(big.NewInt(7), big.NewInt(13), nil)
	factors := FactorizeBig(n)
	ntheoryCheckFactorBig(t, n, factors)
	if len(factors) != 1 || factors[0].Prime.Cmp(big.NewInt(7)) != 0 || factors[0].Exponent != 13 {
		t.Errorf("FactorizeBig(7^13) = %v, want {7,13}", factors)
	}
}

func TestPollardRhoBig(t *testing.T) {
	// Composites: expect a proper divisor.
	for _, s := range []string{"15", "8051", "600851475143", "100000000000000000039000000000000000001"} {
		n := ntheoryBigFromString(t, s)
		d := PollardRhoBig(n)
		if d.Cmp(big.NewInt(1)) <= 0 || d.Cmp(n) >= 0 || new(big.Int).Mod(n, d).Sign() != 0 {
			t.Errorf("PollardRhoBig(%s) = %s is not a proper divisor", s, d)
		}
	}
	// Primes: expect a copy of n.
	for _, s := range []string{"2", "17", "9999991", "10000000000000000051"} {
		n := ntheoryBigFromString(t, s)
		if d := PollardRhoBig(n); d.Cmp(n) != 0 {
			t.Errorf("PollardRhoBig(%s) = %s, want %s (prime)", s, d, s)
		}
	}
}

func TestEulerPhiBig(t *testing.T) {
	tests := []struct {
		n, want string
	}{
		{"0", "0"},
		{"1", "1"},
		{"10", "4"},
		{"36", "12"},
		{"-36", "12"},
		{"1000000", "400000"},
		{"1000003", "1000002"}, // prime p -> p-1
	}
	for _, tt := range tests {
		n := ntheoryBigFromString(t, tt.n)
		want := ntheoryBigFromString(t, tt.want)
		if got := EulerPhiBig(n); got.Cmp(want) != 0 {
			t.Errorf("EulerPhiBig(%s) = %s, want %s", tt.n, got, tt.want)
		}
	}
	// Cross-check against the int64 EulerPhi over a deterministic range.
	for n := int64(1); n <= 200; n++ {
		got := EulerPhiBig(big.NewInt(n))
		if got.Int64() != EulerPhi(n) {
			t.Errorf("EulerPhiBig(%d) = %s, want %d", n, got, EulerPhi(n))
		}
	}
}

func TestCountDivisorsBig(t *testing.T) {
	tests := []struct {
		n, want string
	}{
		{"0", "0"},
		{"1", "1"},
		{"12", "6"},
		{"360", "24"},
		{"-360", "24"},
		{"1000003", "2"}, // prime
	}
	for _, tt := range tests {
		n := ntheoryBigFromString(t, tt.n)
		want := ntheoryBigFromString(t, tt.want)
		if got := CountDivisorsBig(n); got.Cmp(want) != 0 {
			t.Errorf("CountDivisorsBig(%s) = %s, want %s", tt.n, got, tt.want)
		}
	}
	// Cross-check against the int64 CountDivisors over a deterministic range.
	for n := int64(1); n <= 200; n++ {
		got := CountDivisorsBig(big.NewInt(n))
		if got.Int64() != CountDivisors(n) {
			t.Errorf("CountDivisorsBig(%d) = %s, want %d", n, got, CountDivisors(n))
		}
	}
}

func BenchmarkFactorizeU64Semiprime(b *testing.B) {
	const n = uint64(4294967291) * uint64(4294967279)
	for i := 0; i < b.N; i++ {
		FactorizeU64(n)
	}
}

func BenchmarkPollardBrentU64(b *testing.B) {
	const n = uint64(4294967291) * uint64(4294967279)
	for i := 0; i < b.N; i++ {
		PollardBrentU64(n)
	}
}

func BenchmarkFactorizeBigSemiprime(b *testing.B) {
	p, _ := new(big.Int).SetString("10000000000000000051", 10)
	q, _ := new(big.Int).SetString("10000000000000000147", 10)
	n := new(big.Int).Mul(p, q)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		FactorizeBig(n)
	}
}

func BenchmarkPollardRhoBig(b *testing.B) {
	p, _ := new(big.Int).SetString("10000000000000000051", 10)
	q, _ := new(big.Int).SetString("10000000000000000147", 10)
	n := new(big.Int).Mul(p, q)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		PollardRhoBig(n)
	}
}
