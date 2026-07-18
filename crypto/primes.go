package crypto

import (
	"math/big"
	"math/rand"
)

// cryptoSmallPrimes holds the primes below 100, used for fast trial division
// before the more expensive probabilistic tests are run.
var cryptoSmallPrimes = []int64{
	2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
	53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
}

// cryptoDeterministicWitnesses is the set of Miller-Rabin bases that together
// prove primality deterministically for every n < 3.3 * 10^24. See Sorenson &
// Webster, "Strong pseudoprimes to twelve prime bases".
var cryptoDeterministicWitnesses = []int64{2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}

// millerRabinPass reports whether n passes one round of the Miller-Rabin strong
// probable-prime test to the given base a, where n-1 = d * 2^s with d odd. It
// assumes 2 <= a < n and n odd.
func cryptoMillerRabinPass(n, d *big.Int, s int, a *big.Int) bool {
	x := ModExp(a, d, n)
	nMinus1 := new(big.Int).Sub(n, cryptoOne)
	if x.Cmp(cryptoOne) == 0 || x.Cmp(nMinus1) == 0 {
		return true
	}
	for i := 1; i < s; i++ {
		x.Mul(x, x)
		x.Mod(x, n)
		if x.Cmp(nMinus1) == 0 {
			return true
		}
	}
	return false
}

// cryptoDecompose writes n-1 as d * 2^s with d odd and returns d and s.
func cryptoDecompose(n *big.Int) (d *big.Int, s int) {
	d = new(big.Int).Sub(n, cryptoOne)
	for d.Bit(0) == 0 {
		d.Rsh(d, 1)
		s++
	}
	return d, s
}

// MillerRabin runs the probabilistic Miller-Rabin primality test on n using the
// given number of independently chosen random bases (rounds) drawn from rng. It
// returns true if n is a probable prime and false if n is definitely composite.
// The probability that a composite passes all rounds is at most 4^-rounds.
// Inputs n <= 3 and even n are handled directly. It panics if rounds < 1 or rng
// is nil for n large enough to require random bases.
func MillerRabin(n *big.Int, rounds int, rng *rand.Rand) bool {
	if rounds < 1 {
		panic("crypto: MillerRabin requires rounds >= 1")
	}
	if n.Cmp(cryptoTwo) < 0 {
		return false
	}
	if n.Cmp(cryptoThree) <= 0 {
		return true
	}
	if n.Bit(0) == 0 {
		return false
	}
	d, s := cryptoDecompose(n)
	// Choose bases uniformly from [2, n-2].
	lo := big.NewInt(2)
	hi := new(big.Int).Sub(n, cryptoTwo)
	for i := 0; i < rounds; i++ {
		a := cryptoRandRange(rng, lo, hi)
		if !cryptoMillerRabinPass(n, d, s, a) {
			return false
		}
	}
	return true
}

// MillerRabinDeterministic performs a deterministic Miller-Rabin test using a
// fixed set of small prime bases. Its verdict is provably correct for every
// n < 3.3 * 10^24; for larger n it is an extremely strong probable-prime test
// (a composite would have to be a strong pseudoprime to all twelve of the first
// prime bases). It requires no randomness.
func MillerRabinDeterministic(n *big.Int) bool {
	if n.Cmp(cryptoTwo) < 0 {
		return false
	}
	for _, p := range cryptoSmallPrimes {
		bp := big.NewInt(p)
		if n.Cmp(bp) == 0 {
			return true
		}
		if new(big.Int).Mod(n, bp).Sign() == 0 {
			return false
		}
	}
	d, s := cryptoDecompose(n)
	for _, w := range cryptoDeterministicWitnesses {
		a := big.NewInt(w)
		if a.Cmp(n) >= 0 {
			continue
		}
		if !cryptoMillerRabinPass(n, d, s, a) {
			return false
		}
	}
	return true
}

// FermatTest runs the Fermat probable-prime test on n using the given number of
// random bases a drawn from rng, checking a^(n-1) ≡ 1 (mod n). It returns true
// if n passes every round. Note that Carmichael numbers pass the Fermat test to
// all coprime bases, so a positive result is weaker than Miller-Rabin. It
// panics if rounds < 1.
func FermatTest(n *big.Int, rounds int, rng *rand.Rand) bool {
	if rounds < 1 {
		panic("crypto: FermatTest requires rounds >= 1")
	}
	if n.Cmp(cryptoTwo) < 0 {
		return false
	}
	if n.Cmp(cryptoThree) <= 0 {
		return true
	}
	if n.Bit(0) == 0 {
		return false
	}
	nMinus1 := new(big.Int).Sub(n, cryptoOne)
	lo := big.NewInt(2)
	hi := new(big.Int).Sub(n, cryptoTwo)
	for i := 0; i < rounds; i++ {
		a := cryptoRandRange(rng, lo, hi)
		if GCD(a, n).Cmp(cryptoOne) != 0 {
			return false
		}
		if ModExp(a, nMinus1, n).Cmp(cryptoOne) != 0 {
			return false
		}
	}
	return true
}

// IsPrime reports whether n is prime. Small n are decided by trial division
// against the primes below 100; the rest is settled by the deterministic
// Miller-Rabin witness set, which is exact for all n < 3.3 * 10^24 and an
// exceptionally strong test beyond. It needs no randomness and is the
// recommended general-purpose primality predicate in this package.
func IsPrime(n *big.Int) bool {
	return MillerRabinDeterministic(n)
}

// IsProbablePrime reports whether n is a probable prime, combining fast trial
// division by small primes with the deterministic Miller-Rabin witness set and
// then rounds additional random Miller-Rabin bases drawn from rng for extra
// assurance on very large inputs. rounds may be zero, in which case it reduces
// to IsPrime. It panics if rounds < 0.
func IsProbablePrime(n *big.Int, rounds int, rng *rand.Rand) bool {
	if rounds < 0 {
		panic("crypto: IsProbablePrime requires rounds >= 0")
	}
	if !MillerRabinDeterministic(n) {
		return false
	}
	if rounds == 0 {
		return true
	}
	return MillerRabin(n, rounds, rng)
}

// TrialDivision reports whether n is prime by trial division against every
// candidate up to sqrt(n). It is exact for all n but only practical for small
// inputs; use IsPrime for anything large. n < 2 is not prime.
func TrialDivision(n *big.Int) bool {
	if n.Cmp(cryptoTwo) < 0 {
		return false
	}
	if n.Cmp(cryptoThree) <= 0 {
		return true
	}
	if n.Bit(0) == 0 {
		return false
	}
	i := big.NewInt(3)
	iSq := new(big.Int)
	for iSq.Mul(i, i); iSq.Cmp(n) <= 0; iSq.Mul(i, i) {
		if new(big.Int).Mod(n, i).Sign() == 0 {
			return false
		}
		i.Add(i, cryptoTwo)
	}
	return true
}

// SieveOfEratosthenes returns all primes p with 2 <= p <= limit as int64
// values, computed with the classic sieve. It panics if limit < 0. For limit <
// 2 the result is empty.
func SieveOfEratosthenes(limit int64) []int64 {
	if limit < 0 {
		panic("crypto: SieveOfEratosthenes requires limit >= 0")
	}
	if limit < 2 {
		return []int64{}
	}
	composite := make([]bool, limit+1)
	for i := int64(2); i*i <= limit; i++ {
		if !composite[i] {
			for j := i * i; j <= limit; j += i {
				composite[j] = true
			}
		}
	}
	var primes []int64
	for i := int64(2); i <= limit; i++ {
		if !composite[i] {
			primes = append(primes, i)
		}
	}
	return primes
}

// PrimesUpTo is an alias for SieveOfEratosthenes returning every prime not
// exceeding limit.
func PrimesUpTo(limit int64) []int64 {
	return SieveOfEratosthenes(limit)
}

// NextPrime returns the smallest prime strictly greater than n. For n < 2 it
// returns 2. The search uses the deterministic primality test.
func NextPrime(n *big.Int) *big.Int {
	if n.Cmp(cryptoTwo) < 0 {
		return big.NewInt(2)
	}
	cand := new(big.Int).Add(n, cryptoOne)
	if cand.Bit(0) == 0 {
		cand.Add(cand, cryptoOne)
	}
	for !IsPrime(cand) {
		cand.Add(cand, cryptoTwo)
	}
	return cand
}

// PrevPrime returns the largest prime strictly less than n. It returns nil when
// no such prime exists (n <= 2).
func PrevPrime(n *big.Int) *big.Int {
	if n.Cmp(cryptoThree) <= 0 {
		if n.Cmp(cryptoThree) == 0 {
			return big.NewInt(2)
		}
		return nil
	}
	cand := new(big.Int).Sub(n, cryptoOne)
	if cand.Bit(0) == 0 {
		cand.Sub(cand, cryptoOne)
	}
	for cand.Cmp(cryptoTwo) > 0 && !IsPrime(cand) {
		cand.Sub(cand, cryptoTwo)
	}
	return cand
}

// RandomPrime returns a random prime with exactly the requested bit length
// (both the top and bottom bits are forced to 1, guaranteeing an odd number of
// the requested size), using rng as the randomness source. bits must be at
// least 2. The determinism of rng makes the result reproducible.
func RandomPrime(bits int, rng *rand.Rand) *big.Int {
	if bits < 2 {
		panic("crypto: RandomPrime requires bits >= 2")
	}
	topBit := uint(bits - 1)
	for {
		cand := cryptoRandBits(rng, bits)
		cand.SetBit(cand, 0, 1)           // force odd
		cand.SetBit(cand, int(topBit), 1) // force exact bit length
		if IsPrime(cand) {
			return cand
		}
	}
}

// IsSafePrime reports whether p is a safe prime, that is a prime for which
// (p-1)/2 is also prime. Safe primes are preferred for Diffie-Hellman because
// the corresponding group has a large prime-order subgroup.
func IsSafePrime(p *big.Int) bool {
	if !IsPrime(p) {
		return false
	}
	q := new(big.Int).Sub(p, cryptoOne)
	q.Rsh(q, 1)
	return IsPrime(q)
}

// GenerateSafePrime returns a random safe prime p (one for which (p-1)/2 is
// also prime) with exactly the requested bit length, drawing randomness from
// rng. bits must be at least 3. Because safe primes are comparatively rare this
// may test many candidates.
func GenerateSafePrime(bits int, rng *rand.Rand) *big.Int {
	if bits < 3 {
		panic("crypto: GenerateSafePrime requires bits >= 3")
	}
	for {
		q := RandomPrime(bits-1, rng)
		p := new(big.Int).Lsh(q, 1)
		p.Add(p, cryptoOne)
		if p.BitLen() == bits && IsPrime(p) {
			return p
		}
	}
}
