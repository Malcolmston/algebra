package ntheory

import (
	"math/big"
	"testing"
)

// Well-known primes and composites spanning the full uint64 range.
var (
	// Primes, including the two Mersenne primes M31 and M61 and the largest
	// prime below 2^64 (2^64 - 59).
	knownU64Primes = []uint64{
		2, 3, 5, 7, 11, 13, 61, 97, 7919, 104729, 1000003,
		2147483647,           // M31 = 2^31 - 1
		2305843009213693951,  // M61 = 2^61 - 1
		18446744073709551557, // largest prime < 2^64 = 2^64 - 59
	}
	// Composites, including Carmichael numbers and strong pseudoprimes to
	// several small bases, plus 2^64 - 1.
	knownU64Composites = []uint64{
		0, 1, 4, 9, 100, 561, 1105, 1729, 2465, 41041, // Carmichael numbers among these
		25326001,             // strong pseudoprime to bases 2, 3, 5
		3215031751,           // strong pseudoprime to bases 2, 3, 5, 7
		2152302898747,        // strong pseudoprime to bases 2, 3, 5, 7, 11
		2147483647 * 3,       // 3 * M31
		18446744073709551615, // 2^64 - 1
	}
)

func TestIsPrimeU64(t *testing.T) {
	for _, n := range knownU64Primes {
		if !IsPrimeU64(n) {
			t.Errorf("IsPrimeU64(%d) = false, want true", n)
		}
	}
	for _, n := range knownU64Composites {
		if IsPrimeU64(n) {
			t.Errorf("IsPrimeU64(%d) = true, want false", n)
		}
	}
}

// TestIsPrimeU64AgreesWithIsPrime cross-checks the uint64 test against the
// existing int64 IsPrime over a small contiguous range.
func TestIsPrimeU64AgreesWithIsPrime(t *testing.T) {
	for n := int64(0); n < 2000; n++ {
		if got, want := IsPrimeU64(uint64(n)), IsPrime(n); got != want {
			t.Errorf("IsPrimeU64(%d)=%v but IsPrime(%d)=%v", n, got, n, want)
		}
	}
}

func TestMillerRabinU64(t *testing.T) {
	// 3215031751 is a strong probable prime to bases 2, 3, 5 and 7 but is
	// composite; base 11 exposes it.
	cases := []struct {
		n, a uint64
		want bool
	}{
		{2, 2, true},
		{3, 2, true},
		{97, 2, true},
		{97, 96, true},
		{561, 2, false},  // Carmichael, composite
		{221, 174, true}, // strong liar: 221 = 13*17 passes base 174
		{221, 137, false},
		{3215031751, 2, true},
		{3215031751, 7, true},
		{3215031751, 11, false},
		{1, 2, false},
		{4, 3, false},
	}
	for _, c := range cases {
		if got := MillerRabinU64(c.n, c.a); got != c.want {
			t.Errorf("MillerRabinU64(%d,%d)=%v want %v", c.n, c.a, got, c.want)
		}
	}
}

func TestNextPrimeU64(t *testing.T) {
	cases := []struct{ n, want uint64 }{
		{0, 2}, {1, 2}, {2, 3}, {3, 5}, {4, 5}, {5, 7}, {6, 7}, {7, 11},
		{10, 11}, {11, 13}, {13, 17}, {100, 101}, {104728, 104729},
		{2147483647, 2147483659},
		{2305843009213693950, 2305843009213693951},
	}
	for _, c := range cases {
		if got := NextPrimeU64(c.n); got != c.want {
			t.Errorf("NextPrimeU64(%d)=%d want %d", c.n, got, c.want)
		}
	}
	// No prime greater than the largest uint64 prime fits: expect 0.
	if got := NextPrimeU64(18446744073709551557); got != 0 {
		t.Errorf("NextPrimeU64(largest prime)=%d want 0", got)
	}
}

func TestPrevPrimeU64(t *testing.T) {
	cases := []struct {
		n    uint64
		want uint64
		ok   bool
	}{
		{0, 0, false}, {1, 0, false}, {2, 0, false},
		{3, 2, true}, {4, 3, true}, {5, 3, true}, {6, 5, true},
		{7, 5, true}, {8, 7, true}, {11, 7, true}, {12, 11, true},
		{101, 97, true}, {104729, 104723, true},
		{2147483659, 2147483647, true},
		{2305843009213693951, 2305843009213693921, true},
	}
	for _, c := range cases {
		got, ok := PrevPrimeU64(c.n)
		if got != c.want || ok != c.ok {
			t.Errorf("PrevPrimeU64(%d)=(%d,%v) want (%d,%v)", c.n, got, ok, c.want, c.ok)
		}
	}
}

// TestNextPrevPrimeU64RoundTrip checks that stepping forward then back (and
// vice versa) returns to the original prime.
func TestNextPrevPrimeU64RoundTrip(t *testing.T) {
	for _, p := range []uint64{2, 3, 5, 97, 7919, 104729, 1000003} {
		if got := NextPrimeU64(p); got == 0 {
			t.Fatalf("NextPrimeU64(%d) unexpectedly 0", p)
		} else if back, ok := PrevPrimeU64(got); !ok || back != p {
			t.Errorf("PrevPrimeU64(NextPrimeU64(%d))=(%d,%v) want (%d,true)", p, back, ok, p)
		}
	}
}

func TestIsProbablePrimeBig(t *testing.T) {
	primes := []string{
		"2", "3", "97", "7919", "2305843009213693951",
		"18446744073709551557",
		// A 100-digit prime.
		"2074722246773485207821695222107608587480996474721117292752992589912196684750549658310084416732550077",
	}
	composites := []string{
		"0", "1", "4", "561", "41041", "3215031751",
		"18446744073709551615", // 2^64 - 1
		// Product of two large primes (nextprime(1e15) * nextprime(1e16)).
		"10000000000000431000000000002257",
	}
	for _, s := range primes {
		n, _ := new(big.Int).SetString(s, 10)
		if !IsProbablePrimeBig(n) {
			t.Errorf("IsProbablePrimeBig(%s) = false, want true", s)
		}
	}
	for _, s := range composites {
		n, _ := new(big.Int).SetString(s, 10)
		if IsProbablePrimeBig(n) {
			t.Errorf("IsProbablePrimeBig(%s) = true, want false", s)
		}
	}
	// Negative values are never prime.
	if IsProbablePrimeBig(big.NewInt(-7)) {
		t.Errorf("IsProbablePrimeBig(-7) = true, want false")
	}
}

// TestIsProbablePrimeBigAgreesWithStdlib cross-checks against the standard
// library over a small range.
func TestIsProbablePrimeBigAgreesWithStdlib(t *testing.T) {
	for i := int64(0); i < 3000; i++ {
		n := big.NewInt(i)
		if got, want := IsProbablePrimeBig(n), n.ProbablyPrime(20); got != want {
			t.Errorf("IsProbablePrimeBig(%d)=%v but ProbablyPrime=%v", i, got, want)
		}
	}
}

func TestNextPrimeBig(t *testing.T) {
	cases := []struct{ n, want string }{
		{"-5", "2"}, {"0", "2"}, {"1", "2"}, {"2", "3"}, {"3", "5"},
		{"6", "7"}, {"7", "11"}, {"100", "101"},
		{"2305843009213693950", "2305843009213693951"},
	}
	for _, c := range cases {
		n, _ := new(big.Int).SetString(c.n, 10)
		want, _ := new(big.Int).SetString(c.want, 10)
		if got := NextPrimeBig(n); got.Cmp(want) != 0 {
			t.Errorf("NextPrimeBig(%s)=%s want %s", c.n, got, c.want)
		}
	}
}

func TestPrevPrimeBig(t *testing.T) {
	cases := []struct {
		n, want string
		ok      bool
	}{
		{"2", "", false}, {"1", "", false}, {"0", "", false}, {"-3", "", false},
		{"3", "2", true}, {"5", "3", true}, {"8", "7", true}, {"101", "97", true},
		{"2305843009213693951", "2305843009213693921", true},
	}
	for _, c := range cases {
		n, _ := new(big.Int).SetString(c.n, 10)
		got, ok := PrevPrimeBig(n)
		if ok != c.ok {
			t.Errorf("PrevPrimeBig(%s) ok=%v want %v", c.n, ok, c.ok)
			continue
		}
		if ok {
			want, _ := new(big.Int).SetString(c.want, 10)
			if got.Cmp(want) != 0 {
				t.Errorf("PrevPrimeBig(%s)=%s want %s", c.n, got, c.want)
			}
		}
	}
}

func BenchmarkIsPrimeU64(b *testing.B) {
	const n = 18446744073709551557 // largest prime < 2^64
	for i := 0; i < b.N; i++ {
		if !IsPrimeU64(n) {
			b.Fatal("expected prime")
		}
	}
}

func BenchmarkNextPrimeU64(b *testing.B) {
	const start = 1 << 40
	for i := 0; i < b.N; i++ {
		NextPrimeU64(start)
	}
}

func BenchmarkPrevPrimeU64(b *testing.B) {
	const start = 1 << 40
	for i := 0; i < b.N; i++ {
		PrevPrimeU64(start)
	}
}

func BenchmarkIsProbablePrimeBig(b *testing.B) {
	n, _ := new(big.Int).SetString(
		"2074722246773485207821695222107608587480996474721117292752992589912196684750549658310084416732550077", 10)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if !IsProbablePrimeBig(n) {
			b.Fatal("expected prime")
		}
	}
}

func BenchmarkNextPrimeBig(b *testing.B) {
	n, _ := new(big.Int).SetString("2305843009213693950", 10)
	for i := 0; i < b.N; i++ {
		NextPrimeBig(n)
	}
}
