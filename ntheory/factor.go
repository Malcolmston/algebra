package ntheory

import (
	"math/big"
	"sort"
)

// This file adds subexponential integer factorization on top of the package's
// pure O(sqrt n) trial-division [Factorize]. Trial division is int64-only and
// hopeless on large semiprimes; the routines here strip small primes first and
// then split the composite cofactor with Pollard's rho (Brent's cycle-finding
// variant), which runs in roughly O(n^1/4) time per prime factor. The uint64
// path multiplies through Montgomery reduction (MulMont from fastmod.go) and
// batches the greatest-common-divisor steps, and there is a parallel math/big
// path for arbitrary-precision inputs that reuses the memoized base primes from
// sieve.go.
//
// Everything is deterministic: the polynomial constant and start point follow a
// fixed schedule (c = 1, 2, 3, ... with x0 = 2), so repeated calls on the same
// input always return the same answer. No randomness is used anywhere.

// ntheoryBigTrialLimit is the upper bound on the small primes trial-divided out
// before [FactorizeBig] falls back to the heavier [PollardRhoBig]. The primes up
// to this bound are taken from the memoized base-prime cache in sieve.go.
const ntheoryBigTrialLimit uint64 = 1 << 16

// ntheoryBrentStep advances the Pollard-Brent iteration one step in the
// Montgomery domain, returning f(y) = y^2 + c mod n. Both y and cm are
// Montgomery-form residues in [0, n); the returned value is again a
// Montgomery-form residue. The squaring uses [Montgomery.MulMont] (division-free
// REDC) and the constant is added with [AddModU64] so the sum never overflows a
// uint64 even when the modulus is close to 2^64.
func ntheoryBrentStep(mont *Montgomery, y, cm uint64) uint64 {
	y = mont.MulMont(y, y)
	return AddModU64(y, cm, mont.Modulus())
}

// PollardBrentU64 returns a non-trivial factor of the composite integer n using
// Brent's improved cycle-finding variant of Pollard's rho algorithm. The
// returned factor is strictly between 1 and n but is not necessarily prime.
//
// If n is prime it returns n itself, and for an even n it returns 2 (so
// PollardBrentU64(2) == 2 consistently reports the prime). Values below 2 are
// returned unchanged.
//
// The iteration f(x) = x^2 + c mod n runs entirely in the Montgomery domain via
// [Montgomery.MulMont] from fastmod.go, and the differences |x - y| are
// accumulated into a running product whose greatest common divisor with n is
// taken only once per batch, so a hardware division is avoided on the hot path.
// The search restarts with the next constant c on failure. The seed schedule is
// fixed (c = 1, 2, 3, ... with x0 = 2), making the output reproducible.
func PollardBrentU64(n uint64) uint64 {
	if n < 2 {
		return n
	}
	if n%2 == 0 {
		return 2
	}
	if IsPrimeU64(n) {
		return n
	}
	// n is now an odd composite >= 9, so a Montgomery context is valid.
	mont := NewMontgomery(n)
	one := mont.one // Montgomery form of 1, used as the product accumulator seed.
	for c := uint64(1); c < n; c++ {
		cm := mont.ToMont(c % n)
		y := mont.ToMont(2 % n) // x0 = 2
		m := uint64(128)
		var g, r uint64 = 1, 1
		q := one
		var x, ys uint64
		for g == 1 {
			x = y
			for i := uint64(0); i < r; i++ {
				y = ntheoryBrentStep(mont, y, cm)
			}
			for k := uint64(0); k < r && g == 1; k += m {
				ys = y
				lim := m
				if r-k < m {
					lim = r - k
				}
				for i := uint64(0); i < lim; i++ {
					y = ntheoryBrentStep(mont, y, cm)
					// The Montgomery-domain difference |x - y| equals the true
					// difference times R mod n; since gcd(R, n) == 1, folding it
					// into the product via MulMont preserves the gcd with n.
					q = mont.MulMont(q, ntheoryAbsDiffU64(x, y))
				}
				g = ntheoryGCDU64(q, n)
			}
			r *= 2
		}
		if g == n {
			// The batch overshot (a whole factor fell into the product). Retrace
			// one step at a time from the saved position to isolate it.
			g = 1
			for g == 1 {
				ys = ntheoryBrentStep(mont, ys, cm)
				g = ntheoryGCDU64(ntheoryAbsDiffU64(x, ys), n)
			}
		}
		if g != n && g != 1 {
			return g
		}
	}
	return n
}

// ntheoryFactorU64Rec recursively factors n into prime powers, accumulating
// exponents into factors. n must be > 1. It tests primality with [IsPrimeU64]
// and splits composites with [PollardBrentU64].
func ntheoryFactorU64Rec(n uint64, factors map[uint64]int) {
	if n == 1 {
		return
	}
	if IsPrimeU64(n) {
		factors[n]++
		return
	}
	d := PollardBrentU64(n)
	if d <= 1 || d == n {
		// Should not happen for a composite, but guard against infinite recursion.
		factors[n]++
		return
	}
	ntheoryFactorU64Rec(d, factors)
	ntheoryFactorU64Rec(n/d, factors)
}

// FactorizeU64 returns the complete prime factorization of n as a map from each
// prime factor to its (positive) multiplicity.
//
// It first strips small primes by trial division, tests the residual with
// [IsPrimeU64], and splits any remaining composite with [PollardBrentU64],
// recursing until only primes remain. This makes it subexponential and
// dramatically faster than the pure trial-division [Factorize] on large
// semiprimes, while covering the full uint64 range. FactorizeU64(0) and
// FactorizeU64(1) return an empty (non-nil) map.
func FactorizeU64(n uint64) map[uint64]int {
	factors := make(map[uint64]int)
	if n < 2 {
		return factors
	}
	for _, p := range ntheoryTrialPrimes {
		if p*p > n {
			break
		}
		for n%p == 0 {
			factors[p]++
			n /= p
		}
	}
	if n > 1 {
		ntheoryFactorU64Rec(n, factors)
	}
	return factors
}

// PrimePowerBig pairs an arbitrary-precision prime with its exponent in a
// factorization. It is the [math/big] analogue of [PrimePower].
type PrimePowerBig struct {
	// Prime is the prime factor.
	Prime *big.Int
	// Exponent is the multiplicity of Prime in the factorization (>= 1).
	Exponent int
}

// PollardRhoBig returns a non-trivial factor of the composite integer n using
// Brent's variant of Pollard's rho over [math/big], with the greatest common
// divisor taken periodically over a batched product of differences. The sign of
// n is ignored and the returned factor lies strictly between 1 and |n|.
//
// If |n| is (probably) prime it returns a copy of |n|; for an even |n| it
// returns 2. The polynomial constant follows the fixed schedule c = 1, 2, 3,
// ... with start point x0 = 2, so the result is deterministic. Primality is
// checked with [IsProbablePrimeBig].
func PollardRhoBig(n *big.Int) *big.Int {
	m := new(big.Int).Abs(n)
	two := big.NewInt(2)
	if m.Cmp(two) < 0 {
		return m
	}
	if m.Bit(0) == 0 {
		return big.NewInt(2)
	}
	if IsProbablePrimeBig(m) {
		return m
	}
	one := big.NewInt(1)
	x := new(big.Int)
	y := new(big.Int)
	ys := new(big.Int)
	diff := new(big.Int)
	q := new(big.Int)
	g := new(big.Int)
	c := new(big.Int)
	step := func(z *big.Int) {
		z.Mul(z, z)
		z.Add(z, c)
		z.Mod(z, m)
	}
	for cInt := int64(1); ; cInt++ {
		c.SetInt64(cInt)
		y.SetInt64(2) // x0 = 2
		q.SetInt64(1)
		g.SetInt64(1)
		var r int64 = 1
		const mBatch int64 = 128
		for g.Cmp(one) == 0 {
			x.Set(y)
			for i := int64(0); i < r; i++ {
				step(y)
			}
			for k := int64(0); k < r && g.Cmp(one) == 0; k += mBatch {
				ys.Set(y)
				lim := mBatch
				if r-k < mBatch {
					lim = r - k
				}
				for i := int64(0); i < lim; i++ {
					step(y)
					diff.Sub(x, y)
					diff.Abs(diff)
					q.Mul(q, diff)
					q.Mod(q, m)
				}
				g.GCD(nil, nil, q, m)
			}
			r *= 2
		}
		if g.Cmp(m) == 0 {
			// The batch overshot; retrace one step at a time to isolate a factor.
			g.SetInt64(1)
			for g.Cmp(one) == 0 {
				step(ys)
				diff.Sub(x, ys)
				diff.Abs(diff)
				g.GCD(nil, nil, diff, m)
			}
		}
		if g.Cmp(m) != 0 && g.Cmp(one) != 0 {
			return g
		}
	}
}

// ntheoryFactorBigRec recursively factors n, feeding each discovered prime and
// its multiplicity to add. n must be > 1. It tests primality with
// [IsProbablePrimeBig] and splits composites with [PollardRhoBig].
func ntheoryFactorBigRec(n *big.Int, add func(p *big.Int, e int)) {
	if n.Cmp(big.NewInt(1)) == 0 {
		return
	}
	if IsProbablePrimeBig(n) {
		add(n, 1)
		return
	}
	d := PollardRhoBig(n)
	if d.Cmp(big.NewInt(1)) <= 0 || d.Cmp(n) == 0 {
		// Should not happen for a composite, but guard against infinite recursion.
		add(n, 1)
		return
	}
	ntheoryFactorBigRec(d, add)
	ntheoryFactorBigRec(new(big.Int).Quo(n, d), add)
}

// FactorizeBig returns the prime factorization of |n| as a slice of
// [PrimePowerBig] values sorted by ascending prime. The sign of n is ignored.
//
// It trial-divides by the cached small primes from sieve.go (up to
// [ntheoryBigTrialLimit]), tests the residual with [IsProbablePrimeBig], and
// splits any remaining composite with [PollardRhoBig], recursing until only
// primes remain. It is the arbitrary-precision counterpart of [FactorList] and
// returns nil for |n| < 2.
func FactorizeBig(n *big.Int) []PrimePowerBig {
	m := new(big.Int).Abs(n)
	if m.Cmp(big.NewInt(2)) < 0 {
		return nil
	}
	counts := make(map[string]int)
	primesByKey := make(map[string]*big.Int)
	add := func(p *big.Int, e int) {
		k := p.String()
		if _, ok := primesByKey[k]; !ok {
			primesByKey[k] = new(big.Int).Set(p)
		}
		counts[k] += e
	}

	// Strip small primes using the memoized base-prime cache from sieve.go.
	smallPrimes := ntheorySieveBasePrimes(ntheoryBigTrialLimit)
	bp := new(big.Int)
	sq := new(big.Int)
	rem := new(big.Int)
	for _, p := range smallPrimes {
		bp.SetUint64(p)
		if sq.Mul(bp, bp).Cmp(m) > 0 {
			break
		}
		for {
			rem.Mod(m, bp)
			if rem.Sign() != 0 {
				break
			}
			add(bp, 1)
			m.Quo(m, bp)
		}
	}
	if m.Cmp(big.NewInt(1)) > 0 {
		ntheoryFactorBigRec(m, add)
	}

	list := make([]PrimePowerBig, 0, len(counts))
	for k, e := range counts {
		list = append(list, PrimePowerBig{Prime: primesByKey[k], Exponent: e})
	}
	sort.Slice(list, func(i, j int) bool {
		return list[i].Prime.Cmp(list[j].Prime) < 0
	})
	return list
}

// EulerPhiBig returns Euler's totient phi(|n|): the count of integers in [1, |n|]
// that are coprime to |n|. The sign of n is ignored. By convention
// EulerPhiBig(0) == 0 and EulerPhiBig(1) == 1. It is the arbitrary-precision
// counterpart of [EulerPhi] and derives from [FactorizeBig] via the product
// phi(n) = prod p^(e-1) * (p - 1) over the prime powers p^e dividing n.
func EulerPhiBig(n *big.Int) *big.Int {
	m := new(big.Int).Abs(n)
	if m.Sign() == 0 {
		return big.NewInt(0)
	}
	result := big.NewInt(1)
	for _, pp := range FactorizeBig(m) {
		// term = p^(e-1) * (p - 1)
		term := new(big.Int).Exp(pp.Prime, big.NewInt(int64(pp.Exponent-1)), nil)
		term.Mul(term, new(big.Int).Sub(pp.Prime, big.NewInt(1)))
		result.Mul(result, term)
	}
	return result
}

// CountDivisorsBig returns tau(|n|), the number of positive divisors of |n|. The
// sign of n is ignored. CountDivisorsBig(0) == 0 and CountDivisorsBig(1) == 1.
// It is the arbitrary-precision counterpart of [CountDivisors] and derives from
// [FactorizeBig] via the product tau(n) = prod (e + 1) over the exponents e in
// the factorization.
func CountDivisorsBig(n *big.Int) *big.Int {
	m := new(big.Int).Abs(n)
	if m.Sign() == 0 {
		return big.NewInt(0)
	}
	result := big.NewInt(1)
	for _, pp := range FactorizeBig(m) {
		result.Mul(result, big.NewInt(int64(pp.Exponent+1)))
	}
	return result
}
