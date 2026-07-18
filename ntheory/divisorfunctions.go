package ntheory

// DivisorSigma returns σ_k(n), the sum of the k-th powers of the positive
// divisors of n: Σ_{d | n} d^k. The sign of n is ignored.
//
// Special cases: σ_0(n) is the number of divisors τ(n), and σ_1(n) is the sum of
// divisors σ(n). k must be non-negative. DivisorSigma(k, 0) == 0, and for n >= 1
// the value is computed from the prime factorization via the multiplicative
// formula σ_k(p^a) = (p^(k(a+1)) - 1)/(p^k - 1) for k > 0, and (a+1) for k == 0.
func DivisorSigma(k, n int64) int64 {
	if k < 0 {
		panic("ntheory: DivisorSigma requires k >= 0")
	}
	n = abs64(n)
	if n == 0 {
		return 0
	}
	if k == 0 {
		return CountDivisors(n)
	}
	result := int64(1)
	for p, a := range Factorize(n) {
		pk := int64(1)
		for i := int64(0); i < k; i++ {
			pk *= p // p^k
		}
		// term = 1 + p^k + p^2k + ... + p^ak
		term := int64(1)
		cur := int64(1)
		for i := 0; i < a; i++ {
			cur *= pk
			term += cur
		}
		result *= term
	}
	return result
}

// ntheorySpfSieve returns a slice spf of length n+1 where spf[i] is the smallest
// prime factor of i for i >= 2 (spf[0] and spf[1] are 0). It underpins the fast
// linear-time multiplicative-function sieves below. n must be >= 0.
func ntheorySpfSieve(n int64) []int64 {
	spf := make([]int64, n+1)
	for i := int64(2); i <= n; i++ {
		if spf[i] == 0 {
			for j := i; j <= n; j += i {
				if spf[j] == 0 {
					spf[j] = i
				}
			}
		}
	}
	return spf
}

// TotientSieve returns a slice phi of length n+1 where phi[i] == φ(i), Euler's
// totient of i, for every 0 <= i <= n. By convention phi[0] == 0 and phi[1] == 1.
//
// It computes all values in a single sieve pass, which is far faster than
// calling [EulerPhi] on each i independently when the whole table is needed.
// n must be non-negative.
func TotientSieve(n int64) []int64 {
	if n < 0 {
		panic("ntheory: TotientSieve requires n >= 0")
	}
	phi := make([]int64, n+1)
	for i := int64(0); i <= n; i++ {
		phi[i] = i
	}
	for i := int64(2); i <= n; i++ {
		if phi[i] == i { // i is prime
			for j := i; j <= n; j += i {
				phi[j] -= phi[j] / i
			}
		}
	}
	if n >= 0 {
		phi[0] = 0
	}
	return phi
}

// MobiusSieve returns a slice mu of length n+1 where mu[i] == μ(i), the Möbius
// function of i, for every 0 <= i <= n. By convention mu[0] == 0 and mu[1] == 1.
//
// All values are produced in one sieve pass using smallest-prime factorization,
// which is much faster than evaluating [MobiusMu] at each i separately. n must be
// non-negative.
func MobiusSieve(n int64) []int8 {
	if n < 0 {
		panic("ntheory: MobiusSieve requires n >= 0")
	}
	mu := make([]int8, n+1)
	if n >= 1 {
		mu[1] = 1
	}
	spf := ntheorySpfSieve(n)
	for i := int64(2); i <= n; i++ {
		p := spf[i]
		m := i / p
		if m%p == 0 {
			// p^2 divides i, so μ(i) = 0.
			mu[i] = 0
		} else {
			mu[i] = -mu[m]
		}
	}
	return mu
}

// MertensFunction returns M(n) = Σ_{k=1}^{n} μ(k), the Mertens function, the
// running sum of the Möbius function up to n. MertensFunction(0) == 0. It sieves
// the Möbius values once and accumulates them. n must be non-negative.
func MertensFunction(n int64) int64 {
	if n < 0 {
		panic("ntheory: MertensFunction requires n >= 0")
	}
	if n == 0 {
		return 0
	}
	mu := MobiusSieve(n)
	var sum int64
	for i := int64(1); i <= n; i++ {
		sum += int64(mu[i])
	}
	return sum
}
