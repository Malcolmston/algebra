package ntheory

import (
	"math/big"
	"math/bits"
)

// This file complements the int64-only [IsPrime], [IsPrimeBig] and [NextPrime]
// with primality routines that cover the full uint64 range and add the
// previously missing "previous prime" queries (the SymPy prevprime gap).
//
// The performance-sensitive uint64 powering is done with Montgomery
// multiplication supplied by fastmod.go (NewMontgomery/PowMont), falling back
// to its scalar ModPowU64 for the even-modulus path. Candidate stepping for
// the Next*/Prev* helpers uses a mod-30 wheel that skips every multiple of 2,
// 3 and 5, testing only the eight residues coprime to 30.
//
// Everything here is deterministic: no randomness is used anywhere, so repeated
// calls on the same input always return the same answer.

// ntheorySmallPrimes64 lists every prime below 64. IsPrimeU64 trial-divides by
// these before running the more expensive strong-probable-prime tests, which
// resolves the overwhelming majority of composite inputs cheaply.
var ntheorySmallPrimes64 = []uint64{
	2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61,
}

// ntheoryU64Witnesses is the seven-base set found by Jim Sinclair (a
// Jaeschke-style deterministic set). A strong-probable-prime test to every one
// of these bases is a proven, deterministic primality test for all
// n < 2^64, so it is exact across the entire uint64 range. (The classic
// minimal set {2,3,5,7,11,13,17,19,23,29,31,37} — proven for
// n < 3,317,044,064,679,887,385,961,981 — would also work but needs more
// rounds.)
var ntheoryU64Witnesses = []uint64{
	2, 325, 9375, 28178, 450775, 9780504, 1795265022,
}

// ntheoryWheelInc30 and ntheoryWheelDec30 are the forward and backward step
// tables of a mod-30 wheel. For a residue r in [0,30), ntheoryWheelInc30[r] is
// the distance to the next larger integer whose residue mod 30 is coprime to
// 30, and ntheoryWheelDec30[r] the distance to the next smaller such integer.
// Advancing by these amounts visits exactly the eight residues {1,7,11,13,17,
// 19,23,29}, skipping every multiple of 2, 3 and 5.
var ntheoryWheelInc30, ntheoryWheelDec30 = func() ([30]uint64, [30]uint64) {
	coprime := func(x int) bool {
		x = ((x % 30) + 30) % 30
		return x%2 != 0 && x%3 != 0 && x%5 != 0
	}
	var inc, dec [30]uint64
	for r := 0; r < 30; r++ {
		d := 1
		for !coprime(r + d) {
			d++
		}
		inc[r] = uint64(d)
		d = 1
		for !coprime(r - d) {
			d++
		}
		dec[r] = uint64(d)
	}
	return inc, dec
}()

// ntheoryMulModU64 returns (a*b) mod n for a, b < n and n > 0, using a 128-bit
// intermediate product so nothing overflows. The precondition a, b < n
// guarantees the high word of the product is below n, so bits.Div64 never
// panics.
func ntheoryMulModU64(a, b, n uint64) uint64 {
	hi, lo := bits.Mul64(a, b)
	_, rem := bits.Div64(hi, lo, n)
	return rem
}

// ntheoryPowMod returns base^exp mod n. For an odd modulus it uses Montgomery
// multiplication from fastmod.go (NewMontgomery/PowMont); for an even modulus
// it falls back to the scalar ModPowU64. base is assumed already reduced modulo
// n.
func ntheoryPowMod(n, base, exp uint64) uint64 {
	if n&1 == 1 {
		m := NewMontgomery(n)
		return m.PowMont(base, exp)
	}
	return ModPowU64(base, exp, n)
}

// MillerRabinU64 runs a single strong-probable-prime (Miller-Rabin) test of n
// to base a and reports whether n passes it. A return of false proves n is
// composite; a return of true means n is a strong probable prime to base a
// (definitely prime only after enough well-chosen bases pass — see
// [IsPrimeU64]).
//
// The base is reduced modulo n; a base that is a multiple of n carries no
// information and is reported as passing. n < 2 fails and n = 2 passes; every
// other even n fails. The inner powering is Montgomery-accelerated, so the test
// is cheap enough to reuse from factoring code such as factor.go.
func MillerRabinU64(n, a uint64) bool {
	if n < 2 {
		return false
	}
	if n%2 == 0 {
		return n == 2
	}
	if a >= n {
		a %= n
	}
	if a == 0 {
		return true // base ≡ 0 (mod n): not a usable witness, treat as passing.
	}
	// Write n-1 = d * 2^s with d odd.
	d := n - 1
	s := bits.TrailingZeros64(d)
	d >>= uint(s)

	x := ntheoryPowMod(n, a, d)
	if x == 1 || x == n-1 {
		return true
	}
	for i := 1; i < s; i++ {
		x = ntheoryMulModU64(x, x, n)
		if x == n-1 {
			return true
		}
	}
	return false
}

// IsPrimeU64 reports whether n is prime. It is a deterministic test that is
// exact for every value of a uint64 (never merely probable).
//
// It first trial-divides by the primes below 64, which rejects most composites
// immediately, then runs strong-probable-prime tests to the fixed base set
// [ntheoryU64Witnesses], which is proven to give the correct answer for all
// n < 2^64. The inner powering uses Montgomery multiplication from fastmod.go
// whenever the modulus is odd. Values below 2 are not prime.
//
// This is the uint64 counterpart of the int64-only [IsPrime]; use it when
// operands can exceed 2^63.
func IsPrimeU64(n uint64) bool {
	if n < 2 {
		return false
	}
	for _, p := range ntheorySmallPrimes64 {
		if n == p {
			return true
		}
		if n%p == 0 {
			return false
		}
	}
	// No factor below 64; the seven-base set is a proven test for n < 2^64.
	for _, a := range ntheoryU64Witnesses {
		if !MillerRabinU64(n, a) {
			return false
		}
	}
	return true
}

// NextPrimeU64 returns the least prime strictly greater than n.
//
// Unlike the int64 [NextPrime] it works correctly for values near and above
// 2^63 without overflow. Candidate stepping uses a mod-30 wheel so only the
// eight residues coprime to 30 are tested. If no prime greater than n is
// representable in a uint64 (that is, n is at least the largest uint64 prime
// 18446744073709551557), it returns 0.
func NextPrimeU64(n uint64) uint64 {
	if n < 7 {
		// The primes 2, 3, 5 lie on residues that the mod-30 wheel skips, so
		// resolve the small range directly.
		for _, p := range []uint64{2, 3, 5, 7} {
			if p > n {
				return p
			}
		}
	}
	// Every prime greater than 7 is coprime to 30. Walk the wheel; the loop
	// condition also detects the uint64 overflow that ends the search when no
	// larger prime fits.
	for c := n + 1; c > n; c += ntheoryWheelInc30[c%30] {
		if c%2 != 0 && c%3 != 0 && c%5 != 0 && IsPrimeU64(c) {
			return c
		}
	}
	return 0
}

// PrevPrimeU64 returns the greatest prime strictly less than n together with
// ok == true. When no prime is smaller than n (that is, n <= 2) it returns
// (0, false).
//
// This fills the "previous prime" gap (SymPy's prevprime) that the package
// otherwise lacked entirely. Candidate stepping uses the descending mod-30
// wheel for speed.
func PrevPrimeU64(n uint64) (uint64, bool) {
	if n <= 2 {
		return 0, false
	}
	if n <= 7 {
		// 2, 3, 5 are skipped by the wheel, so handle the small range directly.
		for _, p := range []uint64{5, 3, 2} {
			if p < n {
				return p, true
			}
		}
		return 0, false
	}
	// c stays >= 7 throughout, and the maximum backward step is 6, so the
	// subtraction never underflows.
	for c := n - 1; c >= 7; c -= ntheoryWheelDec30[c%30] {
		if c%2 != 0 && c%3 != 0 && c%5 != 0 && IsPrimeU64(c) {
			return c, true
		}
	}
	for _, p := range []uint64{5, 3, 2} {
		if p < n {
			return p, true
		}
	}
	return 0, false
}

// ntheorySprpBase2Big reports whether n is a strong probable prime to base 2.
// n must be odd and greater than 2.
func ntheorySprpBase2Big(n *big.Int) bool {
	one := big.NewInt(1)
	nm1 := new(big.Int).Sub(n, one)
	s := nm1.TrailingZeroBits()
	d := new(big.Int).Rsh(nm1, s)

	x := new(big.Int).Exp(big.NewInt(2), d, n)
	if x.Cmp(one) == 0 || x.Cmp(nm1) == 0 {
		return true
	}
	for i := uint(1); i < s; i++ {
		x.Mul(x, x)
		x.Mod(x, n)
		if x.Cmp(nm1) == 0 {
			return true
		}
	}
	return false
}

// ntheoryHalveMod sets x to (x / 2) mod n for odd n, interpreting the division
// modularly: if x is odd it is first shifted up by n (which is even after the
// addition, since n is odd) before halving. x must already be in [0, n).
func ntheoryHalveMod(x, n *big.Int) {
	if x.Bit(0) == 1 {
		x.Add(x, n)
	}
	x.Rsh(x, 1)
	if x.Cmp(n) >= 0 {
		x.Mod(x, n)
	}
}

// ntheoryStrongLucasBig reports whether n passes the strong Lucas-Selfridge
// probable-prime test. n must be odd and greater than 2. A perfect-square n is
// rejected first so the Selfridge parameter search always terminates.
func ntheoryStrongLucasBig(n *big.Int) bool {
	// A perfect square never yields a Jacobi symbol of -1, which would make the
	// parameter search below loop forever; such n is composite.
	sqrt := new(big.Int).Sqrt(n)
	if new(big.Int).Mul(sqrt, sqrt).Cmp(n) == 0 {
		return false
	}

	// Selfridge's method: try D = 5, -7, 9, -11, 13, ... until (D/n) = -1.
	D := big.NewInt(5)
	two := big.NewInt(2)
	for {
		j := big.Jacobi(D, n)
		if j == -1 {
			break
		}
		if j == 0 {
			// gcd(D, n) > 1 with |D| < n means n has a small factor.
			absD := new(big.Int).Abs(D)
			return absD.Cmp(n) == 0
		}
		if D.Sign() > 0 {
			D.Add(D, two)
			D.Neg(D)
		} else {
			D.Neg(D)
			D.Add(D, two)
		}
	}

	// P = 1, Q = (1 - D) / 4 (exact, since D ≡ 1 (mod 4)).
	P := big.NewInt(1)
	Q := new(big.Int).Sub(big.NewInt(1), D)
	Q.Quo(Q, big.NewInt(4))

	// Reduce the parameters modulo n into [0, n).
	Pm := new(big.Int).Mod(P, n)
	Dm := new(big.Int).Mod(D, n)
	Qk := new(big.Int).Mod(Q, n)

	// n + 1 = d * 2^s with d odd.
	np1 := new(big.Int).Add(n, big.NewInt(1))
	s := np1.TrailingZeroBits()
	d := new(big.Int).Rsh(np1, s)

	// Lucas sequences U_k, V_k and Q^k, starting at k = 1.
	U := big.NewInt(1)
	V := new(big.Int).Set(Pm)

	tmp := new(big.Int)
	for i := d.BitLen() - 2; i >= 0; i-- {
		// Doubling: U_2k = U_k V_k; V_2k = V_k^2 - 2 Q^k; Q^2k = (Q^k)^2.
		U.Mul(U, V)
		U.Mod(U, n)

		V.Mul(V, V)
		tmp.Lsh(Qk, 1)
		V.Sub(V, tmp)
		V.Mod(V, n)

		Qk.Mul(Qk, Qk)
		Qk.Mod(Qk, n)

		if d.Bit(i) == 1 {
			// Indexing up by one:
			//   U_{2k+1} = (P U_2k + V_2k) / 2
			//   V_{2k+1} = (D U_2k + P V_2k) / 2
			newU := new(big.Int).Mul(Pm, U)
			newU.Add(newU, V)
			newU.Mod(newU, n)
			ntheoryHalveMod(newU, n)

			newV := new(big.Int).Mul(Dm, U)
			tmp.Mul(Pm, V)
			newV.Add(newV, tmp)
			newV.Mod(newV, n)
			ntheoryHalveMod(newV, n)

			U, V = newU, newV

			Qk.Mul(Qk, Q)
			Qk.Mod(Qk, n)
		}
	}

	// Strong test: probable prime iff U_d ≡ 0, or V_{d·2^r} ≡ 0 for some
	// 0 <= r < s.
	if U.Sign() == 0 {
		return true
	}
	for r := uint(0); r < s; r++ {
		if V.Sign() == 0 {
			return true
		}
		V.Mul(V, V)
		tmp.Lsh(Qk, 1)
		V.Sub(V, tmp)
		V.Mod(V, n)

		Qk.Mul(Qk, Qk)
		Qk.Mod(Qk, n)
	}
	return false
}

// IsProbablePrimeBig reports whether the arbitrary-precision n is (very
// probably) prime using a Baillie-PSW test: a strong probable-prime test to
// base 2 together with a strong Lucas-Selfridge test, then additional
// Miller-Rabin rounds delegated to [math/big.Int.ProbablyPrime]. No composite
// is known to pass, and because none of the steps use randomness the result is
// deterministic for a given n. Values below 2 are not prime.
func IsProbablePrimeBig(n *big.Int) bool {
	if n.Sign() <= 0 || n.BitLen() < 2 {
		return false // n <= 1
	}
	if n.Bit(0) == 0 {
		return n.Cmp(big.NewInt(2)) == 0 // even: prime only if 2.
	}
	for _, p := range []int64{3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37} {
		bp := big.NewInt(p)
		switch n.Cmp(bp) {
		case 0:
			return true
		case 1:
			if new(big.Int).Mod(n, bp).Sign() == 0 {
				return false
			}
		}
	}
	// n is odd, greater than 37 and free of small factors: run Baillie-PSW.
	if !ntheorySprpBase2Big(n) {
		return false
	}
	if !ntheoryStrongLucasBig(n) {
		return false
	}
	// Extra Miller-Rabin rounds from the standard library (itself deterministic
	// for a fixed count and input).
	return n.ProbablyPrime(20)
}

// NextPrimeBig returns the least prime strictly greater than n, using
// [IsProbablePrimeBig] to test candidates and a mod-30 wheel to skip multiples
// of 2, 3 and 5. For n < 2 it returns 2. The returned value is a fresh
// *big.Int; n is not modified.
func NextPrimeBig(n *big.Int) *big.Int {
	if n.Cmp(big.NewInt(2)) < 0 {
		return big.NewInt(2)
	}
	c := new(big.Int).Add(n, big.NewInt(1))
	// The primes 2, 3, 5, 7 sit on residues the wheel skips (or below its
	// start), so resolve them explicitly.
	for _, p := range []int64{2, 3, 5, 7} {
		bp := big.NewInt(p)
		if c.Cmp(bp) <= 0 && bp.Cmp(n) > 0 {
			return bp
		}
	}
	thirty := big.NewInt(30)
	mod := new(big.Int)
	for {
		r := mod.Mod(c, thirty).Uint64()
		if r%2 != 0 && r%3 != 0 && r%5 != 0 && IsProbablePrimeBig(c) {
			return new(big.Int).Set(c)
		}
		c.Add(c, new(big.Int).SetUint64(ntheoryWheelInc30[r]))
	}
}

// PrevPrimeBig returns the greatest prime strictly less than n together with
// ok == true, using [IsProbablePrimeBig] and a descending mod-30 wheel. When no
// prime is smaller than n (that is, n <= 2) it returns (nil, false). The
// returned value is a fresh *big.Int; n is not modified.
func PrevPrimeBig(n *big.Int) (*big.Int, bool) {
	if n.Cmp(big.NewInt(2)) <= 0 {
		return nil, false
	}
	seven := big.NewInt(7)
	thirty := big.NewInt(30)
	mod := new(big.Int)
	c := new(big.Int).Sub(n, big.NewInt(1))
	for c.Cmp(seven) >= 0 {
		r := mod.Mod(c, thirty).Uint64()
		if r%2 != 0 && r%3 != 0 && r%5 != 0 && IsProbablePrimeBig(c) {
			return new(big.Int).Set(c), true
		}
		c.Sub(c, new(big.Int).SetUint64(ntheoryWheelDec30[r]))
	}
	// Below 7 the only primes are 5, 3, 2, which the wheel skips.
	for _, p := range []int64{5, 3, 2} {
		bp := big.NewInt(p)
		if bp.Cmp(n) < 0 {
			return bp, true
		}
	}
	return nil, false
}
