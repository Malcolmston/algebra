package ntheory

import (
	"math/big"
	"sort"
)

// SqrtMod returns a square root x in [0, p) of a modulo the prime p, i.e. an x
// with x*x ≡ a (mod p), together with ok == true. When a is a quadratic
// non-residue modulo p (no square root exists) it returns (0, false).
//
// The residue is validated with [LegendreSymbol]. For p ≡ 3 (mod 4) it uses the
// direct closed form x = a^((p+1)/4), avoiding the Tonelli-Shanks inner loop;
// for the remaining odd primes it runs the full Tonelli-Shanks algorithm. All
// powering is performed with [ModPow] so no intermediate value overflows. The
// special cases a ≡ 0 (returning (0, true)) and p == 2 are handled directly.
//
// p must be an odd prime or 2; a is reduced modulo p before use. When the two
// roots x and p-x both exist the smaller one is returned, so the result is
// deterministic. The behavior is undefined when p is not prime.
func SqrtMod(a, p int64) (x int64, ok bool) {
	if p <= 1 {
		panic("ntheory: SqrtMod requires an odd prime p (or 2)")
	}
	a = ((a % p) + p) % p
	if a == 0 {
		return 0, true
	}
	if p == 2 {
		return a % 2, true
	}
	// A root exists only for quadratic residues.
	if LegendreSymbol(a, p) != 1 {
		return 0, false
	}
	// Closed form when p ≡ 3 (mod 4): x = a^((p+1)/4).
	if p%4 == 3 {
		root := ModPow(a, (p+1)/4, p)
		return ntheorySmallerRoot(root, p), true
	}
	// Tonelli-Shanks: write p-1 = q * 2^s with q odd.
	q := p - 1
	s := int64(0)
	for q%2 == 0 {
		q /= 2
		s++
	}
	// Find a quadratic non-residue z to seed the search.
	z := int64(2)
	for LegendreSymbol(z, p) != -1 {
		z++
	}
	m := s
	c := ModPow(z, q, p)
	t := ModPow(a, q, p)
	root := ModPow(a, (q+1)/2, p)
	for t != 1 {
		// Least i in (0, m) with t^(2^i) == 1.
		i := int64(0)
		t2 := t
		for t2 != 1 {
			t2 = mulMod(t2, t2, p)
			i++
			if i == m {
				// Unreachable for a genuine residue modulo a prime.
				return 0, false
			}
		}
		// b = c^(2^(m-i-1)).
		b := c
		for j := int64(0); j < m-i-1; j++ {
			b = mulMod(b, b, p)
		}
		m = i
		c = mulMod(b, b, p)
		t = mulMod(t, c, p)
		root = mulMod(root, b, p)
	}
	return ntheorySmallerRoot(root, p), true
}

// SqrtModBig returns a square root of a modulo the prime p using the
// Tonelli-Shanks algorithm over math/big, suitable for primes far larger than
// int64. It returns (nil, false) when a is a quadratic non-residue modulo p.
//
// Residue membership is decided with the Euler criterion a^((p-1)/2) ≡ 1
// (mod p), evaluated with [big.Int.Exp] (the same square-and-multiply used by
// [JacobiSymbol]-style checks). For p ≡ 3 (mod 4) the direct form
// a^((p+1)/4) is used. The cases a ≡ 0 (returning 0) and p == 2 are handled
// directly. a is reduced modulo p, and the smaller of the two roots is
// returned, making the result deterministic. p must be an odd prime or 2.
func SqrtModBig(a, p *big.Int) (*big.Int, bool) {
	one := big.NewInt(1)
	two := big.NewInt(2)
	three := big.NewInt(3)
	if p.Cmp(two) == 0 {
		return new(big.Int).Mod(a, two), true
	}
	aa := new(big.Int).Mod(a, p)
	if aa.Sign() == 0 {
		return big.NewInt(0), true
	}
	pm1 := new(big.Int).Sub(p, one)
	half := new(big.Int).Rsh(pm1, 1) // (p-1)/2
	// Euler criterion via square-and-multiply exponentiation.
	if new(big.Int).Exp(aa, half, p).Cmp(one) != 0 {
		return nil, false
	}
	// Closed form when p ≡ 3 (mod 4).
	if new(big.Int).And(p, three).Cmp(three) == 0 {
		exp := new(big.Int).Rsh(new(big.Int).Add(p, one), 2) // (p+1)/4
		root := new(big.Int).Exp(aa, exp, p)
		return ntheoryBigSmallerRoot(root, p), true
	}
	// Tonelli-Shanks: p-1 = q * 2^s with q odd.
	q := new(big.Int).Set(pm1)
	s := 0
	for q.Bit(0) == 0 {
		q.Rsh(q, 1)
		s++
	}
	// Find a non-residue z (z^((p-1)/2) ≡ -1).
	z := big.NewInt(2)
	for new(big.Int).Exp(z, half, p).Cmp(pm1) != 0 {
		z.Add(z, one)
	}
	m := s
	c := new(big.Int).Exp(z, q, p)
	qp1o2 := new(big.Int).Rsh(new(big.Int).Add(q, one), 1) // (q+1)/2
	t := new(big.Int).Exp(aa, q, p)
	root := new(big.Int).Exp(aa, qp1o2, p)
	for t.Cmp(one) != 0 {
		// Least i in (0, m) with t^(2^i) == 1.
		i := 0
		t2 := new(big.Int).Set(t)
		for t2.Cmp(one) != 0 {
			t2.Mul(t2, t2)
			t2.Mod(t2, p)
			i++
			if i == m {
				return nil, false
			}
		}
		// b = c^(2^(m-i-1)).
		b := new(big.Int).Set(c)
		for j := 0; j < m-i-1; j++ {
			b.Mul(b, b)
			b.Mod(b, p)
		}
		m = i
		c.Mul(b, b)
		c.Mod(c, p)
		t.Mul(t, c)
		t.Mod(t, p)
		root.Mul(root, b)
		root.Mod(root, p)
	}
	return ntheoryBigSmallerRoot(root, p), true
}

// SqrtModPrimePower returns a square root x in [0, p^k) of a modulo the prime
// power p^k, together with ok == true, or (0, false) when a is a non-residue.
// A root modulo p is first obtained with [SqrtMod] and then lifted to modulus
// p^k by Hensel lifting (Newton iteration on f(x) = x^2 - a). The p == 2 branch
// and the non-unit cases (where p divides a) are handled by the same complete
// prime-power solver used by [AllSqrtModComposite], and the smallest root is
// returned for determinism.
//
// p must be prime and k must be >= 1. The behavior is undefined when p is not
// prime, and it panics when k < 1.
func SqrtModPrimePower(a, p int64, k int) (x int64, ok bool) {
	if k < 1 {
		panic("ntheory: SqrtModPrimePower requires k >= 1")
	}
	if p == 2 {
		roots := ntheorySqrtMod2kAll(a, k)
		if len(roots) == 0 {
			return 0, false
		}
		return roots[0], true
	}
	if k == 1 {
		return SqrtMod(a, p)
	}
	q := ntheoryIntPow(p, k)
	a = ((a % q) + q) % q
	if a == 0 {
		return 0, true
	}
	roots := ntheorySqrtModOddPrimePowerAll(a, p, uint(k))
	if len(roots) == 0 {
		return 0, false
	}
	return roots[0], true
}

// AllSqrtModComposite returns every square root of a modulo the composite
// modulus m, sorted ascending. The modulus is factored with [FactorListUint64],
// each prime-power component is solved for its complete set of roots, and the
// component roots are recombined across all combinations with [CRT]. The result
// is empty when a is a non-residue modulo some component (and therefore has no
// square root modulo m).
//
// m must be greater than 0; it panics when m == 0. a is reduced modulo m before
// use.
func AllSqrtModComposite(a, m uint64) []uint64 {
	if m == 0 {
		panic("ntheory: AllSqrtModComposite requires m > 0")
	}
	if m == 1 {
		return []uint64{0}
	}
	a %= m
	comps := FactorListUint64(m)
	moduli := make([]int64, 0, len(comps))
	rootSets := make([][]int64, 0, len(comps))
	for _, pp := range comps {
		q := ntheoryIntPow(int64(pp.Prime), int(pp.Exponent))
		ra := int64(a % uint64(q))
		var roots []int64
		if pp.Prime == 2 {
			roots = ntheorySqrtMod2kAll(ra, int(pp.Exponent))
		} else {
			roots = ntheorySqrtModOddPrimePowerAll(ra, int64(pp.Prime), pp.Exponent)
		}
		if len(roots) == 0 {
			return nil
		}
		moduli = append(moduli, q)
		rootSets = append(rootSets, roots)
	}
	// Recombine every combination of component roots with CRT.
	results := make([]uint64, 0)
	idx := make([]int, len(rootSets))
	for {
		rs := make([]int64, len(rootSets))
		for i := range rootSets {
			rs[i] = rootSets[i][idx[i]]
		}
		if x, _, ok := CRT(rs, moduli); ok {
			results = append(results, uint64(x))
		}
		// Advance the odometer over the Cartesian product of root sets.
		i := 0
		for ; i < len(idx); i++ {
			idx[i]++
			if idx[i] < len(rootSets[i]) {
				break
			}
			idx[i] = 0
		}
		if i == len(idx) {
			break
		}
	}
	return ntheoryDedupSortU64(results)
}

// ntheoryIntPow returns base**exp for a non-negative exp using repeated
// multiplication. The caller is responsible for keeping the result within int64
// range.
func ntheoryIntPow(base int64, exp int) int64 {
	r := int64(1)
	for i := 0; i < exp; i++ {
		r *= base
	}
	return r
}

// ntheoryBigSmallerRoot returns min(root, p-root), the deterministic
// representative of a square-root pair modulo p.
func ntheoryBigSmallerRoot(root, p *big.Int) *big.Int {
	other := new(big.Int).Sub(p, root)
	if other.Cmp(root) < 0 {
		return other
	}
	return root
}

// ntheorySqrtModPrimePowerUnit lifts a square root of the unit a (gcd(a, p) == 1)
// from modulus p to modulus p^k using Newton/Hensel iteration on f(x)=x^2-a. It
// returns (root, true) or (0, false) when a is a non-residue modulo p. p must be
// an odd prime.
func ntheorySqrtModPrimePowerUnit(a, p int64, k int) (int64, bool) {
	r, ok := SqrtMod(a, p)
	if !ok {
		return 0, false
	}
	mod := p
	for j := 1; j < k; j++ {
		next := mod * p
		ra := ((a % next) + next) % next
		f := ((mulMod(r, r, next)-ra)%next + next) % next
		// 2r is a unit modulo the odd prime power, so the inverse exists.
		inv, iok := ModInverse(mulMod(2, r, next), next)
		if !iok {
			return 0, false
		}
		r = ((r-mulMod(f, inv, next))%next + next) % next
		mod = next
	}
	return r, true
}

// ntheorySqrtModOddPrimePowerAll returns every square root of a modulo q = p^e
// for an odd prime p, sorted ascending. It handles the unit case, the
// a ≡ 0 case, and the intermediate non-unit case via p-adic valuation
// reduction. The slice is empty when a is a non-residue modulo q.
func ntheorySqrtModOddPrimePowerAll(a, p int64, e uint) []int64 {
	q := ntheoryIntPow(p, int(e))
	a = ((a % q) + q) % q
	if a == 0 {
		// Roots of x^2 ≡ 0: the multiples of p^ceil(e/2).
		step := ntheoryIntPow(p, int((e+1)/2))
		res := make([]int64, 0, q/step)
		for x := int64(0); x < q; x += step {
			res = append(res, x)
		}
		return res
	}
	// p-adic valuation t of a, with a = p^t * u and u a unit.
	t := uint(0)
	u := a
	for u%p == 0 {
		u /= p
		t++
	}
	if t%2 == 1 {
		// Odd valuation: no square root (squares have even valuation).
		return nil
	}
	s := t / 2
	reducedExp := int(e - t)
	r0, ok := ntheorySqrtModPrimePowerUnit(u, p, reducedExp)
	if !ok {
		return nil
	}
	reducedMod := ntheoryIntPow(p, reducedExp)
	unitRoots := []int64{r0}
	if other := reducedMod - r0; other != r0 {
		unitRoots = append(unitRoots, other)
	}
	// Each unit root y contributes x = p^s*y + j*p^(e-s) for j in [0, p^s).
	ps := ntheoryIntPow(p, int(s))
	topStep := ntheoryIntPow(p, int(e)-int(s))
	res := make([]int64, 0, len(unitRoots)*int(ps))
	for _, y := range unitRoots {
		base := mulMod(ps, y%q, q)
		for j := int64(0); j < ps; j++ {
			res = append(res, (base+j*topStep)%q)
		}
	}
	return ntheoryDedupSort(res)
}

// ntheorySqrtMod2kAll returns every square root of a modulo 2^k, sorted
// ascending, computed by lifting the root set one bit at a time. The slice is
// empty when a is a non-residue modulo 2^k. k must satisfy 1 <= k < 63.
func ntheorySqrtMod2kAll(a int64, k int) []int64 {
	mod := int64(1) << uint(k)
	a = ((a % mod) + mod) % mod
	// Roots modulo 2.
	roots := make([]int64, 0, 4)
	for x := int64(0); x < 2; x++ {
		if mulMod(x, x, 2) == a%2 {
			roots = append(roots, x)
		}
	}
	curMod := int64(2)
	for j := 2; j <= k; j++ {
		nextMod := curMod << 1
		ra := ((a % nextMod) + nextMod) % nextMod
		next := make([]int64, 0, len(roots)*2)
		for _, r := range roots {
			for _, cand := range [2]int64{r, r + curMod} {
				if mulMod(cand, cand, nextMod) == ra {
					next = append(next, cand)
				}
			}
		}
		roots = next
		curMod = nextMod
	}
	return ntheoryDedupSort(roots)
}

// ntheoryDedupSort sorts a slice of int64 ascending and removes duplicates.
func ntheoryDedupSort(v []int64) []int64 {
	if len(v) == 0 {
		return v
	}
	sort.Slice(v, func(i, j int) bool { return v[i] < v[j] })
	out := v[:1]
	for _, x := range v[1:] {
		if x != out[len(out)-1] {
			out = append(out, x)
		}
	}
	return out
}

// ntheoryDedupSortU64 sorts a slice of uint64 ascending and removes duplicates.
func ntheoryDedupSortU64(v []uint64) []uint64 {
	if len(v) == 0 {
		return v
	}
	sort.Slice(v, func(i, j int) bool { return v[i] < v[j] })
	out := v[:1]
	for _, x := range v[1:] {
		if x != out[len(out)-1] {
			out = append(out, x)
		}
	}
	return out
}
