package galois

import (
	"errors"
	"math/big"
	"math/rand"
	"sort"
)

// distinctPrimeDivisorsInt returns the distinct prime divisors of the small
// positive integer n.
func distinctPrimeDivisorsInt(n int) []int {
	var ps []int
	m := n
	for d := 2; d*d <= m; d++ {
		if m%d == 0 {
			ps = append(ps, d)
			for m%d == 0 {
				m /= d
			}
		}
	}
	if m > 1 {
		ps = append(ps, m)
	}
	return ps
}

// frobeniusIterateX applies the Frobenius map t ↦ t^p modulo f to the residue
// of x exactly k times, returning x^(p^k) mod f.
func frobeniusIterateX(p *big.Int, k int, f *Poly) *Poly {
	h := PolyX(p)
	_, h, _ = h.DivMod(f)
	for i := 0; i < k; i++ {
		h, _ = h.PowMod(p, f)
	}
	return h
}

// IsIrreducible reports whether f is irreducible over GF(p) using Rabin's test.
// Constants and the zero polynomial are not irreducible; every degree-1
// polynomial is irreducible.
func IsIrreducible(f *Poly) bool {
	n := f.Degree()
	if n <= 0 {
		return false
	}
	if n == 1 {
		return true
	}
	fm := f.Monic()
	x := PolyX(f.P)
	for _, d := range distinctPrimeDivisorsInt(n) {
		h := frobeniusIterateX(f.P, n/d, fm)
		if fm.Gcd(h.Sub(x)).Degree() != 0 {
			return false
		}
	}
	h := frobeniusIterateX(f.P, n, fm)
	_, r, _ := h.Sub(x).DivMod(fm)
	return r.IsZero()
}

// IsPrimitivePoly reports whether f is a primitive polynomial over GF(p): a
// monic irreducible polynomial whose root x generates the whole multiplicative
// group of GF(p^deg f), i.e. x has order p^(deg f) - 1.
func IsPrimitivePoly(f *Poly) bool {
	if !IsIrreducible(f) {
		return false
	}
	n := f.Degree()
	fm := f.Monic()
	e := new(big.Int).Sub(IntPow(f.P, n), big1)
	x := PolyX(f.P)
	for _, r := range PrimeFactors(e) {
		exp := new(big.Int).Div(e, r)
		h, _ := x.PowMod(exp, fm)
		if h.IsOne() {
			return false
		}
	}
	return true
}

// monicFromIndex builds the monic degree-n polynomial whose lower coefficients
// are the base-p digits of idx (little-endian).
func monicFromIndex(p *big.Int, n int, idx *big.Int) *Poly {
	coeffs := make([]*big.Int, n+1)
	k := clone(idx)
	for i := 0; i < n; i++ {
		q, r := new(big.Int).DivMod(k, p, new(big.Int))
		coeffs[i] = r
		k = q
	}
	coeffs[n] = big.NewInt(1)
	return NewPolyBig(p, coeffs)
}

// FindIrreducible returns the first monic irreducible polynomial of degree n
// over GF(p), scanning candidates by ascending coefficient index.
func FindIrreducible(p *big.Int, n int) (*Poly, error) {
	if n < 1 {
		return nil, errors.New("galois: degree must be at least 1")
	}
	if !IsPrimeInt(p) {
		return nil, errors.New("galois: characteristic must be prime")
	}
	total := IntPow(p, n)
	for idx := big.NewInt(0); idx.Cmp(total) < 0; idx.Add(idx, big1) {
		f := monicFromIndex(p, n, idx)
		if IsIrreducible(f) {
			return f, nil
		}
	}
	return nil, errors.New("galois: no irreducible polynomial found")
}

// FindPrimitivePoly returns the first primitive polynomial of degree n over
// GF(p), scanning candidates by ascending coefficient index.
func FindPrimitivePoly(p *big.Int, n int) (*Poly, error) {
	if n < 1 {
		return nil, errors.New("galois: degree must be at least 1")
	}
	if !IsPrimeInt(p) {
		return nil, errors.New("galois: characteristic must be prime")
	}
	total := IntPow(p, n)
	for idx := big.NewInt(0); idx.Cmp(total) < 0; idx.Add(idx, big1) {
		f := monicFromIndex(p, n, idx)
		if IsPrimitivePoly(f) {
			return f, nil
		}
	}
	return nil, errors.New("galois: no primitive polynomial found")
}

// CanonicalIrreduciblePoly returns the deterministic canonical (lexicographically
// least by coefficient index) monic irreducible polynomial of degree n over
// GF(p). It is identical to FindIrreducible.
func CanonicalIrreduciblePoly(p *big.Int, n int) (*Poly, error) {
	return FindIrreducible(p, n)
}

// ConwayStylePoly returns a deterministic Conway-style defining polynomial for
// GF(p^n): the lexicographically least primitive polynomial of degree n over
// GF(p). Unlike true Conway polynomials it does not enforce norm-compatibility
// across subfields, but it is a well-defined, reproducible canonical choice.
func ConwayStylePoly(p *big.Int, n int) (*Poly, error) {
	return FindPrimitivePoly(p, n)
}

// RandomIrreducible returns a random monic irreducible polynomial of degree n
// over GF(p), drawing candidate coefficient vectors from r until one is
// irreducible.
func RandomIrreducible(p *big.Int, n int, r *rand.Rand) (*Poly, error) {
	if n < 1 {
		return nil, errors.New("galois: degree must be at least 1")
	}
	total := IntPow(p, n)
	for tries := 0; tries < 100000; tries++ {
		idx := new(big.Int).Rand(r, total)
		f := monicFromIndex(p, n, idx)
		if IsIrreducible(f) {
			return f, nil
		}
	}
	return nil, errors.New("galois: failed to find a random irreducible polynomial")
}

// AllIrreduciblePolys returns every monic irreducible polynomial of degree n
// over GF(p). It enumerates all p^n candidates and is intended for small fields.
func AllIrreduciblePolys(p *big.Int, n int) []*Poly {
	total := IntPow(p, n)
	if !total.IsInt64() {
		return nil
	}
	var out []*Poly
	for idx := big.NewInt(0); idx.Cmp(total) < 0; idx.Add(idx, big1) {
		f := monicFromIndex(p, n, idx)
		if IsIrreducible(f) {
			out = append(out, f)
		}
	}
	return out
}

// NumberOfIrreducibles returns the number of monic irreducible polynomials of
// degree n over GF(p), via the necklace formula
// (1/n)·Σ_{d|n} μ(d)·p^(n/d).
func NumberOfIrreducibles(p *big.Int, n int) *big.Int {
	sum := big.NewInt(0)
	for _, d := range divisorsInt64(n) {
		term := IntPow(p, n/d)
		mu := MobiusMu(big.NewInt(int64(d)))
		switch mu {
		case 1:
			sum.Add(sum, term)
		case -1:
			sum.Sub(sum, term)
		}
	}
	return sum.Div(sum, big.NewInt(int64(n)))
}

// CyclotomicCoset returns the q-cyclotomic coset of s modulo q^n - 1, that is
// the sorted set {s, s·q, s·q², …} reduced modulo q^n - 1.
func CyclotomicCoset(q *big.Int, s, n int) []int64 {
	N := new(big.Int).Sub(IntPow(q, n), big1).Int64()
	qq := q.Int64()
	seen := map[int64]bool{}
	var coset []int64
	v := int64(((s % int(N)) + int(N))) % N
	for !seen[v] {
		seen[v] = true
		coset = append(coset, v)
		v = (v * qq) % N
	}
	sort.Slice(coset, func(i, j int) bool { return coset[i] < coset[j] })
	return coset
}

// CyclotomicCosets returns the partition of {0, 1, …, q^n-2} into q-cyclotomic
// cosets modulo q^n - 1, each coset sorted and the list ordered by smallest
// representative.
func CyclotomicCosets(q *big.Int, n int) [][]int64 {
	N := new(big.Int).Sub(IntPow(q, n), big1).Int64()
	used := make([]bool, N)
	var cosets [][]int64
	for s := int64(0); s < N; s++ {
		if used[s] {
			continue
		}
		c := CyclotomicCoset(q, int(s), n)
		for _, v := range c {
			used[v] = true
		}
		cosets = append(cosets, c)
	}
	return cosets
}
