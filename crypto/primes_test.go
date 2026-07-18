package crypto

import (
	"math/big"
	"math/rand"
	"testing"
)

func TestIsPrimeKnown(t *testing.T) {
	primes := map[int64]bool{
		2: true, 3: true, 5: true, 7: true, 11: true, 13: true, 97: true,
		561:  false,                           // Carmichael number
		1105: false, 1729: false, 2465: false, // Carmichael numbers
		7919: true, 104729: true, 1000003: true,
		1: false, 0: false, 4: false, 100: false, 999: false,
	}
	for n, want := range primes {
		if got := IsPrime(bi(n)); got != want {
			t.Errorf("IsPrime(%d)=%v want %v", n, got, want)
		}
		// TrialDivision must agree for these small inputs.
		if got := TrialDivision(bi(n)); got != want {
			t.Errorf("TrialDivision(%d)=%v want %v", n, got, want)
		}
	}
	// A large known Mersenne prime 2^61-1.
	m61 := new(big.Int).Sub(new(big.Int).Lsh(big.NewInt(1), 61), big.NewInt(1))
	if !IsPrime(m61) {
		t.Error("IsPrime(2^61-1) want true")
	}
}

func TestMillerRabinCarmichael(t *testing.T) {
	rng := rand.New(rand.NewSource(42))
	// 561 is a Carmichael number: Fermat can be fooled, Miller-Rabin must not be.
	if MillerRabinDeterministic(bi(561)) {
		t.Error("MillerRabinDeterministic(561) want false")
	}
	if MillerRabin(bi(561), 20, rng) {
		t.Error("MillerRabin(561) want false")
	}
	if !MillerRabin(bi(7919), 20, rng) {
		t.Error("MillerRabin(7919) want true")
	}
}

func TestFermatTest(t *testing.T) {
	rng := rand.New(rand.NewSource(7))
	if !FermatTest(bi(7919), 10, rng) {
		t.Error("FermatTest(7919) want true")
	}
	if FermatTest(bi(15), 10, rng) {
		t.Error("FermatTest(15) want false")
	}
}

func TestSieve(t *testing.T) {
	got := SieveOfEratosthenes(30)
	want := []int64{2, 3, 5, 7, 11, 13, 17, 19, 23, 29}
	if len(got) != len(want) {
		t.Fatalf("Sieve(30) len=%d want %d", len(got), len(want))
	}
	for i := range want {
		if got[i] != want[i] {
			t.Errorf("Sieve(30)[%d]=%d want %d", i, got[i], want[i])
		}
	}
	// Count of primes below 100 is 25.
	if n := len(PrimesUpTo(100)); n != 25 {
		t.Errorf("pi(100)=%d want 25", n)
	}
}

func TestNextPrevPrime(t *testing.T) {
	cases := []struct{ n, next, prev int64 }{
		{14, 17, 13},
		{100, 101, 97},
		{7, 11, 5},
	}
	for _, c := range cases {
		if np := NextPrime(bi(c.n)); np.Int64() != c.next {
			t.Errorf("NextPrime(%d)=%d want %d", c.n, np.Int64(), c.next)
		}
		if pp := PrevPrime(bi(c.n)); pp.Int64() != c.prev {
			t.Errorf("PrevPrime(%d)=%d want %d", c.n, pp.Int64(), c.prev)
		}
	}
}

func TestRandomPrimeDeterministic(t *testing.T) {
	rng1 := rand.New(rand.NewSource(123))
	rng2 := rand.New(rand.NewSource(123))
	for i := 0; i < 3; i++ {
		p1 := RandomPrime(64, rng1)
		p2 := RandomPrime(64, rng2)
		if p1.Cmp(p2) != 0 {
			t.Errorf("RandomPrime not deterministic: %v vs %v", p1, p2)
		}
		if !IsPrime(p1) {
			t.Errorf("RandomPrime returned composite %v", p1)
		}
		if p1.BitLen() != 64 {
			t.Errorf("RandomPrime bitlen=%d want 64", p1.BitLen())
		}
	}
}

func TestSafePrime(t *testing.T) {
	// 11 is safe ((11-1)/2=5 prime); 13 is not ((13-1)/2=6).
	if !IsSafePrime(bi(11)) {
		t.Error("IsSafePrime(11) want true")
	}
	if IsSafePrime(bi(13)) {
		t.Error("IsSafePrime(13) want false")
	}
	rng := rand.New(rand.NewSource(5))
	sp := GenerateSafePrime(16, rng)
	if !IsSafePrime(sp) {
		t.Errorf("GenerateSafePrime produced non-safe prime %v", sp)
	}
	if sp.BitLen() != 16 {
		t.Errorf("GenerateSafePrime bitlen=%d want 16", sp.BitLen())
	}
}
