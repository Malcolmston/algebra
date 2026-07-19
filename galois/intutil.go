package galois

import (
	"math/big"
	"sort"
)

// small reusable constants; never mutated in place.
var (
	big0 = big.NewInt(0)
	big1 = big.NewInt(1)
	big2 = big.NewInt(2)
)

func bi(n int64) *big.Int       { return big.NewInt(n) }
func clone(a *big.Int) *big.Int { return new(big.Int).Set(a) }

// PrimePower records a prime together with its exponent in a factorisation,
// so that n = Prime^Exp times the other factors.
type PrimePower struct {
	Prime *big.Int
	Exp   int
}

// IsPrimeInt reports whether n is a prime integer. It uses a deterministic
// small-factor screen followed by a strong probabilistic (Baillie–PSW style)
// test via math/big, which has no known composite false positives.
func IsPrimeInt(n *big.Int) bool {
	if n.Sign() <= 0 {
		return false
	}
	if n.Cmp(big2) < 0 {
		return false
	}
	return n.ProbablyPrime(40)
}

// NextPrime returns the smallest prime strictly greater than n.
func NextPrime(n *big.Int) *big.Int {
	c := new(big.Int).Add(n, big1)
	if c.Cmp(big2) < 0 {
		return big.NewInt(2)
	}
	if c.Bit(0) == 0 && c.Cmp(big2) != 0 {
		c.Add(c, big1)
	}
	for !c.ProbablyPrime(40) {
		if c.Cmp(big2) == 0 {
			break
		}
		c.Add(c, big2)
	}
	return c
}

// Gcd returns the non-negative greatest common divisor of a and b.
func Gcd(a, b *big.Int) *big.Int {
	return new(big.Int).GCD(nil, nil, new(big.Int).Abs(a), new(big.Int).Abs(b))
}

// ExtendedGcd returns g = gcd(a, b) together with Bézout coefficients x and y
// satisfying a*x + b*y = g.
func ExtendedGcd(a, b *big.Int) (g, x, y *big.Int) {
	x = new(big.Int)
	y = new(big.Int)
	g = new(big.Int).GCD(x, y, a, b)
	return g, x, y
}

// Lcm returns the non-negative least common multiple of a and b. Lcm(0, x) is 0.
func Lcm(a, b *big.Int) *big.Int {
	if a.Sign() == 0 || b.Sign() == 0 {
		return big.NewInt(0)
	}
	g := Gcd(a, b)
	l := new(big.Int).Div(new(big.Int).Abs(a), g)
	l.Mul(l, new(big.Int).Abs(b))
	return l
}

// IntPow returns base raised to the power exp (exp >= 0) as an exact integer.
func IntPow(base *big.Int, exp int) *big.Int {
	if exp < 0 {
		return big.NewInt(0)
	}
	return new(big.Int).Exp(base, big.NewInt(int64(exp)), nil)
}

// pollardRho returns a non-trivial factor of the composite odd integer n.
func pollardRho(n *big.Int) *big.Int {
	if n.Bit(0) == 0 {
		return big.NewInt(2)
	}
	x := big.NewInt(2)
	y := big.NewInt(2)
	c := big.NewInt(1)
	d := big.NewInt(1)
	f := func(v *big.Int) *big.Int {
		r := new(big.Int).Mul(v, v)
		r.Add(r, c)
		r.Mod(r, n)
		return r
	}
	for d.Cmp(big1) == 0 {
		x = f(x)
		y = f(f(y))
		diff := new(big.Int).Sub(x, y)
		diff.Abs(diff)
		if diff.Sign() == 0 {
			// cycle without a factor; bump the increment and restart.
			c.Add(c, big1)
			x = big.NewInt(2)
			y = big.NewInt(2)
			d = big.NewInt(1)
			continue
		}
		d = new(big.Int).GCD(nil, nil, diff, n)
	}
	return d
}

func factorRec(n *big.Int, out map[string]int) {
	if n.Cmp(big1) == 0 {
		return
	}
	if IsPrimeInt(n) {
		out[n.String()]++
		return
	}
	// strip small primes quickly.
	for _, sp := range []int64{2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37} {
		p := big.NewInt(sp)
		for new(big.Int).Mod(n, p).Sign() == 0 {
			out[p.String()]++
			n = new(big.Int).Div(n, p)
		}
	}
	if n.Cmp(big1) == 0 {
		return
	}
	if IsPrimeInt(n) {
		out[n.String()]++
		return
	}
	d := pollardRho(n)
	factorRec(d, out)
	factorRec(new(big.Int).Div(n, d), out)
}

// FactorInt returns the prime factorisation of |n| as a slice of PrimePower
// entries sorted by increasing prime. FactorInt(0) and FactorInt(±1) return an
// empty slice.
func FactorInt(n *big.Int) []PrimePower {
	m := new(big.Int).Abs(n)
	out := map[string]int{}
	if m.Cmp(big1) > 0 {
		factorRec(m, out)
	}
	res := make([]PrimePower, 0, len(out))
	for s, e := range out {
		p, _ := new(big.Int).SetString(s, 10)
		res = append(res, PrimePower{Prime: p, Exp: e})
	}
	sort.Slice(res, func(i, j int) bool { return res[i].Prime.Cmp(res[j].Prime) < 0 })
	return res
}

// PrimeFactors returns the distinct prime factors of |n| in increasing order.
func PrimeFactors(n *big.Int) []*big.Int {
	f := FactorInt(n)
	res := make([]*big.Int, len(f))
	for i, pp := range f {
		res[i] = pp.Prime
	}
	return res
}

// DivisorsInt returns all positive divisors of |n| in increasing order.
// DivisorsInt(0) returns an empty slice.
func DivisorsInt(n *big.Int) []*big.Int {
	m := new(big.Int).Abs(n)
	if m.Sign() == 0 {
		return nil
	}
	divs := []*big.Int{big.NewInt(1)}
	for _, pp := range FactorInt(m) {
		var next []*big.Int
		pk := big.NewInt(1)
		for e := 0; e <= pp.Exp; e++ {
			for _, d := range divs {
				next = append(next, new(big.Int).Mul(d, pk))
			}
			pk = new(big.Int).Mul(pk, pp.Prime)
		}
		divs = next
	}
	sort.Slice(divs, func(i, j int) bool { return divs[i].Cmp(divs[j]) < 0 })
	return divs
}

// NumDivisors returns the number of positive divisors of |n|.
func NumDivisors(n *big.Int) int {
	count := 1
	for _, pp := range FactorInt(n) {
		count *= pp.Exp + 1
	}
	return count
}

// EulerPhi returns Euler's totient φ(n): the count of integers in [1, n]
// coprime to n. EulerPhi(1) == 1.
func EulerPhi(n *big.Int) *big.Int {
	m := new(big.Int).Abs(n)
	if m.Cmp(big1) == 0 {
		return big.NewInt(1)
	}
	result := clone(m)
	for _, pp := range FactorInt(m) {
		// result *= (1 - 1/p)  ==  result / p * (p-1)
		result.Div(result, pp.Prime)
		result.Mul(result, new(big.Int).Sub(pp.Prime, big1))
	}
	return result
}

// MobiusMu returns the Möbius function μ(n): 0 if n is divisible by a square
// greater than 1, otherwise (-1)^k where k is the number of distinct primes.
func MobiusMu(n *big.Int) int {
	m := new(big.Int).Abs(n)
	if m.Cmp(big1) == 0 {
		return 1
	}
	if m.Sign() == 0 {
		return 0
	}
	k := 0
	for _, pp := range FactorInt(m) {
		if pp.Exp > 1 {
			return 0
		}
		k++
	}
	if k%2 == 0 {
		return 1
	}
	return -1
}

// divisorsInt64 is an internal helper returning divisors of a small positive n.
func divisorsInt64(n int) []int {
	var ds []int
	for d := 1; d*d <= n; d++ {
		if n%d == 0 {
			ds = append(ds, d)
			if d != n/d {
				ds = append(ds, n/d)
			}
		}
	}
	sort.Ints(ds)
	return ds
}
