package designs

import (
	"errors"
	"sort"
)

// Gcd returns the non-negative greatest common divisor of a and b. Gcd(0,0) is
// 0. Signs of the inputs are ignored.
func Gcd(a, b int) int {
	if a < 0 {
		a = -a
	}
	if b < 0 {
		b = -b
	}
	for b != 0 {
		a, b = b, a%b
	}
	return a
}

// ExtendedGcd returns g = gcd(a,b) together with Bezout coefficients x, y such
// that a*x + b*y = g.
func ExtendedGcd(a, b int) (g, x, y int) {
	oldR, r := a, b
	oldS, s := 1, 0
	oldT, t := 0, 1
	for r != 0 {
		q := oldR / r
		oldR, r = r, oldR-q*r
		oldS, s = s, oldS-q*s
		oldT, t = t, oldT-q*t
	}
	if oldR < 0 {
		return -oldR, -oldS, -oldT
	}
	return oldR, oldS, oldT
}

// Lcm returns the least common multiple of a and b. Lcm with a zero argument is
// 0. The result is non-negative.
func Lcm(a, b int) int {
	if a == 0 || b == 0 {
		return 0
	}
	g := Gcd(a, b)
	r := (a / g) * b
	if r < 0 {
		return -r
	}
	return r
}

// ModInverse returns the multiplicative inverse of a modulo m as a value in
// [0,m). It reports an error when m<=0 or gcd(a,m) != 1.
func ModInverse(a, m int) (int, error) {
	if m <= 0 {
		return 0, errors.New("designs: modulus must be positive")
	}
	g, x, _ := ExtendedGcd(((a%m)+m)%m, m)
	if g != 1 {
		return 0, errors.New("designs: value not invertible modulo m")
	}
	return ((x % m) + m) % m, nil
}

// ModExp returns base**exp modulo m, computed by binary exponentiation, as a
// value in [0,m). It requires m>0 and exp>=0.
func ModExp(base, exp, m int) int {
	if m == 1 {
		return 0
	}
	result := 1
	base = ((base % m) + m) % m
	for exp > 0 {
		if exp&1 == 1 {
			result = result * base % m
		}
		exp >>= 1
		base = base * base % m
	}
	return result
}

// IsPrime reports whether n is a prime number using trial division.
func IsPrime(n int) bool {
	if n < 2 {
		return false
	}
	if n%2 == 0 {
		return n == 2
	}
	if n%3 == 0 {
		return n == 3
	}
	for i := 5; i*i <= n; i += 6 {
		if n%i == 0 || n%(i+2) == 0 {
			return false
		}
	}
	return true
}

// NextPrime returns the smallest prime strictly greater than n.
func NextPrime(n int) int {
	c := n + 1
	if c < 2 {
		return 2
	}
	for !IsPrime(c) {
		c++
	}
	return c
}

// PrevPrime returns the largest prime strictly less than n, or an error when no
// such prime exists (n<=2).
func PrevPrime(n int) (int, error) {
	c := n - 1
	for c >= 2 {
		if IsPrime(c) {
			return c, nil
		}
		c--
	}
	return 0, errors.New("designs: no prime below n")
}

// PrimesUpTo returns all primes p with p<=n in increasing order, using the
// sieve of Eratosthenes.
func PrimesUpTo(n int) []int {
	if n < 2 {
		return nil
	}
	sieve := make([]bool, n+1)
	var primes []int
	for i := 2; i <= n; i++ {
		if !sieve[i] {
			primes = append(primes, i)
			for j := i * i; j <= n; j += i {
				sieve[j] = true
			}
		}
	}
	return primes
}

// Factorize returns the prime factorization of n>=1 as sorted prime bases with
// their exponents. Factorize(1) returns empty slices.
func Factorize(n int) (primes, exponents []int) {
	if n < 1 {
		return nil, nil
	}
	for p := 2; p*p <= n; p++ {
		if n%p == 0 {
			e := 0
			for n%p == 0 {
				n /= p
				e++
			}
			primes = append(primes, p)
			exponents = append(exponents, e)
		}
	}
	if n > 1 {
		primes = append(primes, n)
		exponents = append(exponents, 1)
	}
	return primes, exponents
}

// IsPrimePower reports whether q equals p**k for a prime p and integer k>=1.
func IsPrimePower(q int) bool {
	if q < 2 {
		return false
	}
	primes, _ := Factorize(q)
	return len(primes) == 1
}

// PrimePowerBase returns the prime p when q = p**k is a prime power, together
// with an error otherwise.
func PrimePowerBase(q int) (int, error) {
	p, _, err := FactorPrimePower(q)
	return p, err
}

// PrimePowerExponent returns the exponent k when q = p**k is a prime power,
// together with an error otherwise.
func PrimePowerExponent(q int) (int, error) {
	_, k, err := FactorPrimePower(q)
	return k, err
}

// FactorPrimePower writes q = p**k with p prime and k>=1, returning p and k. It
// reports an error when q is not a prime power.
func FactorPrimePower(q int) (p, k int, err error) {
	if q < 2 {
		return 0, 0, errors.New("designs: not a prime power")
	}
	primes, exps := Factorize(q)
	if len(primes) != 1 {
		return 0, 0, errors.New("designs: not a prime power")
	}
	return primes[0], exps[0], nil
}

// Divisors returns the positive divisors of n>=1 in increasing order.
func Divisors(n int) []int {
	if n < 1 {
		return nil
	}
	var small, large []int
	for i := 1; i*i <= n; i++ {
		if n%i == 0 {
			small = append(small, i)
			if i != n/i {
				large = append(large, n/i)
			}
		}
	}
	for i := len(large) - 1; i >= 0; i-- {
		small = append(small, large[i])
	}
	return small
}

// NumDivisors returns the number of positive divisors of n>=1.
func NumDivisors(n int) int {
	if n < 1 {
		return 0
	}
	_, exps := Factorize(n)
	c := 1
	for _, e := range exps {
		c *= e + 1
	}
	return c
}

// EulerPhi returns Euler's totient of n>=1, the count of integers in [1,n]
// coprime to n.
func EulerPhi(n int) int {
	if n < 1 {
		return 0
	}
	result := n
	primes, _ := Factorize(n)
	for _, p := range primes {
		result -= result / p
	}
	return result
}

// LegendreSymbol returns the Legendre symbol (a|p) for an odd prime p: it is 0
// when p divides a, +1 when a is a non-zero quadratic residue modulo p, and -1
// when a is a non-residue. It reports an error when p is not an odd prime.
func LegendreSymbol(a, p int) (int, error) {
	if p < 3 || !IsPrime(p) {
		return 0, errors.New("designs: p must be an odd prime")
	}
	a = ((a % p) + p) % p
	if a == 0 {
		return 0, nil
	}
	r := ModExp(a, (p-1)/2, p)
	if r == 1 {
		return 1, nil
	}
	return -1, nil
}

// IsQuadraticResidue reports whether a is a non-zero quadratic residue modulo
// the odd prime p.
func IsQuadraticResidue(a, p int) bool {
	s, err := LegendreSymbol(a, p)
	return err == nil && s == 1
}

// QuadraticResidues returns the non-zero quadratic residues modulo the odd
// prime p in increasing order.
func QuadraticResidues(p int) []int {
	if p < 3 || !IsPrime(p) {
		return nil
	}
	seen := make(map[int]bool)
	for x := 1; x < p; x++ {
		seen[x*x%p] = true
	}
	var out []int
	for r := range seen {
		out = append(out, r)
	}
	sort.Ints(out)
	return out
}

// QuadraticNonResidues returns the quadratic non-residues modulo the odd prime
// p in increasing order.
func QuadraticNonResidues(p int) []int {
	if p < 3 || !IsPrime(p) {
		return nil
	}
	res := make(map[int]bool)
	for _, r := range QuadraticResidues(p) {
		res[r] = true
	}
	var out []int
	for x := 1; x < p; x++ {
		if !res[x] {
			out = append(out, x)
		}
	}
	return out
}

// MultiplicativeOrder returns the multiplicative order of a modulo m, the least
// e>=1 with a**e == 1 (mod m). It reports an error when gcd(a,m) != 1 or m<=1.
func MultiplicativeOrder(a, m int) (int, error) {
	if m <= 1 {
		return 0, errors.New("designs: modulus must exceed 1")
	}
	a = ((a % m) + m) % m
	if Gcd(a, m) != 1 {
		return 0, errors.New("designs: a and m are not coprime")
	}
	e := 1
	cur := a % m
	for cur != 1 {
		cur = cur * a % m
		e++
		if e > m {
			return 0, errors.New("designs: order not found")
		}
	}
	return e, nil
}

// PrimitiveRoot returns the smallest primitive root modulo the prime p (a
// generator of the multiplicative group). It reports an error when p is not
// prime.
func PrimitiveRoot(p int) (int, error) {
	if !IsPrime(p) {
		return 0, errors.New("designs: p must be prime")
	}
	if p == 2 {
		return 1, nil
	}
	phi := p - 1
	primes, _ := Factorize(phi)
	for g := 2; g < p; g++ {
		ok := true
		for _, q := range primes {
			if ModExp(g, phi/q, p) == 1 {
				ok = false
				break
			}
		}
		if ok {
			return g, nil
		}
	}
	return 0, errors.New("designs: no primitive root found")
}

// Binomial returns the binomial coefficient C(n,k) computed without overflow
// for moderate arguments. It returns 0 when k<0 or k>n.
func Binomial(n, k int) int {
	if k < 0 || k > n {
		return 0
	}
	if k > n-k {
		k = n - k
	}
	result := 1
	for i := 0; i < k; i++ {
		result = result * (n - i) / (i + 1)
	}
	return result
}

// Factorial returns n! for 0<=n. It returns an error when n is negative.
func Factorial(n int) (int, error) {
	if n < 0 {
		return 0, errors.New("designs: factorial of negative number")
	}
	r := 1
	for i := 2; i <= n; i++ {
		r *= i
	}
	return r, nil
}

// IntSqrt returns the integer square root floor(sqrt(n)) for n>=0.
func IntSqrt(n int) int {
	if n < 0 {
		return 0
	}
	if n == 0 {
		return 0
	}
	x := n
	y := (x + 1) / 2
	for y < x {
		x = y
		y = (x + n/x) / 2
	}
	return x
}

// IsPerfectSquare reports whether n is a perfect square (n>=0).
func IsPerfectSquare(n int) bool {
	if n < 0 {
		return false
	}
	r := IntSqrt(n)
	return r*r == n
}
