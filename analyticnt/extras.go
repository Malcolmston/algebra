package analyticnt

import (
	"math"
	"math/cmplx"
)

// GCD returns the greatest common divisor of a and b (non-negative).
func GCD(a, b int64) int64 { return gcdInt64(a, b) }

// LCM returns the least common multiple of a and b. LCM(0, x) = 0.
func LCM(a, b int64) int64 {
	if a == 0 || b == 0 {
		return 0
	}
	g := gcdInt64(a, b)
	r := a / g * b
	if r < 0 {
		return -r
	}
	return r
}

// IsCoprime reports whether a and b are coprime (gcd 1).
func IsCoprime(a, b int64) bool { return gcdInt64(a, b) == 1 }

// FactorInt returns the prime factorization of n >= 1 as a map from prime to
// exponent. FactorInt(1) is the empty map.
func FactorInt(n int64) map[int64]int {
	if n < 1 {
		panic("analyticnt: FactorInt requires n >= 1")
	}
	out := make(map[int64]int)
	for p := int64(2); p*p <= n; p++ {
		for n%p == 0 {
			out[p]++
			n /= p
		}
	}
	if n > 1 {
		out[n]++
	}
	return out
}

// PrimeFactors returns the distinct prime factors of n in increasing order.
func PrimeFactors(n int64) []int64 {
	if n < 1 {
		panic("analyticnt: PrimeFactors requires n >= 1")
	}
	return distinctPrimeFactors(n)
}

// PrimeFactorsWithMultiplicity returns the prime factors of n listed with
// repetition, in increasing order, so that their product is n.
func PrimeFactorsWithMultiplicity(n int64) []int64 {
	if n < 1 {
		panic("analyticnt: PrimeFactorsWithMultiplicity requires n >= 1")
	}
	var out []int64
	for p := int64(2); p*p <= n; p++ {
		for n%p == 0 {
			out = append(out, p)
			n /= p
		}
	}
	if n > 1 {
		out = append(out, n)
	}
	return out
}

// IsSophieGermain reports whether p and 2p+1 are both prime (p is a Sophie
// Germain prime).
func IsSophieGermain(p int64) bool {
	return IsPrime(p) && IsPrime(2*p+1)
}

// IsSafePrime reports whether p is a safe prime, i.e. (p−1)/2 is also prime.
func IsSafePrime(p int64) bool {
	return IsPrime(p) && p%2 == 1 && IsPrime((p-1)/2)
}

// SophieGermainPrimesUpTo returns all Sophie Germain primes p (p and 2p+1 both
// prime) with p ≤ x.
func SophieGermainPrimesUpTo(x int) []int64 {
	var out []int64
	for _, p := range PrimesUpTo(x) {
		if IsPrime(2*p + 1) {
			out = append(out, p)
		}
	}
	return out
}

// SafePrimesUpTo returns all safe primes q ≤ x, those for which (q−1)/2 is
// prime.
func SafePrimesUpTo(x int) []int64 {
	var out []int64
	for _, q := range PrimesUpTo(x) {
		if q%2 == 1 && IsPrime((q-1)/2) {
			out = append(out, q)
		}
	}
	return out
}

// BertrandPrime returns a prime p with n < p < 2n for n >= 1 (guaranteed to
// exist by Bertrand's postulate). It returns the smallest such prime.
func BertrandPrime(n int64) int64 {
	if n < 1 {
		panic("analyticnt: BertrandPrime requires n >= 1")
	}
	return NextPrime(n)
}

// FactorialMod returns n! mod m using iterative multiplication.
func FactorialMod(n, m int64) int64 {
	if m == 1 {
		return 0
	}
	r := int64(1 % m)
	for i := int64(2); i <= n; i++ {
		r = mulMod(r, i%m, m)
	}
	return r
}

// IsWilsonPrime reports whether p is a Wilson prime, i.e. p² divides (p−1)! + 1.
// Only three are known (5, 13, 563); this checks by direct computation and is
// suitable for small p.
func IsWilsonPrime(p int64) bool {
	if !IsPrime(p) {
		return false
	}
	m := p * p
	fact := FactorialMod(p-1, m)
	return (fact+1)%m == 0
}

// Bernoulli returns the n-th Bernoulli number B_n as a float64, computed by the
// Akiyama–Tanigawa algorithm. This uses the convention B_1 = +1/2.
func Bernoulli(n int) float64 {
	if n < 0 {
		panic("analyticnt: Bernoulli requires n >= 0")
	}
	a := make([]float64, n+1)
	for m := 0; m <= n; m++ {
		a[m] = 1 / float64(m+1)
		for j := m; j >= 1; j-- {
			a[j-1] = float64(j) * (a[j-1] - a[j])
		}
	}
	return a[0]
}

// Gamma returns the gamma function Γ(x) for real x, delegating to the standard
// library implementation.
func Gamma(x float64) float64 { return math.Gamma(x) }

// LogGamma returns the natural logarithm of |Γ(x)| for real x.
func LogGamma(x float64) float64 {
	lg, _ := math.Lgamma(x)
	return lg
}

// Digamma returns the digamma function ψ(x) = Γ'(x)/Γ(x) for real x > 0, using
// the recurrence to shift the argument upward followed by the standard
// asymptotic series.
func Digamma(x float64) float64 {
	if x <= 0 {
		panic("analyticnt: Digamma requires x > 0")
	}
	result := 0.0
	for x < 6 {
		result -= 1 / x
		x++
	}
	f := 1 / (x * x)
	result += math.Log(x) - 0.5/x +
		f*(-1.0/12+f*(1.0/120+f*(-1.0/252+f*(1.0/240))))
	return result
}

// Trigamma returns the trigamma function ψ₁(x) = d²/dx² ln Γ(x) for real x > 0.
func Trigamma(x float64) float64 {
	if x <= 0 {
		panic("analyticnt: Trigamma requires x > 0")
	}
	result := 0.0
	for x < 6 {
		result += 1 / (x * x)
		x++
	}
	f := 1 / (x * x)
	result += 1/x + 0.5*f + f/x*(1.0/6+f*(-1.0/30+f*(1.0/42)))
	return result
}

// RiemannXiComplex returns the Riemann ξ-function
// ξ(s) = ½ s(s−1) π^{−s/2} Γ(s/2) ζ(s), an entire function whose zeros are
// exactly the nontrivial zeros of ζ and which satisfies ξ(s) = ξ(1−s).
func RiemannXiComplex(s complex128) complex128 {
	half := s / 2
	logG := LogGammaComplex(half)
	pref := 0.5 * s * (s - 1) * cmplx.Exp(logG-half*complex(math.Log(math.Pi), 0))
	return pref * ZetaComplex(s)
}

// RiemannXiReal returns ξ(s) for a real argument s as a complex128 (the value is
// real).
func RiemannXiReal(s float64) complex128 {
	return RiemannXiComplex(complex(s, 0))
}

// CompletedZeta returns the completed zeta function
// ξ̂(s) = π^{−s/2} Γ(s/2) ζ(s), which satisfies the symmetric functional
// equation ξ̂(s) = ξ̂(1−s).
func CompletedZeta(s complex128) complex128 {
	half := s / 2
	logG := LogGammaComplex(half)
	return cmplx.Exp(logG-half*complex(math.Log(math.Pi), 0)) * ZetaComplex(s)
}

// ZetaZeroRealPart returns the common real part 1/2 of every nontrivial zero on
// the critical line; it is a convenience constant-returning helper.
func ZetaZeroRealPart() float64 { return 0.5 }

// PrimePowersUpTo returns every prime power p^k (k >= 1) that is ≤ x, in
// increasing order.
func PrimePowersUpTo(x int64) []int64 {
	if x < 2 {
		return []int64{}
	}
	var out []int64
	for _, p := range PrimesUpTo(int(x)) {
		pk := p
		for pk <= x {
			out = append(out, pk)
			if pk > x/p {
				break
			}
			pk *= p
		}
	}
	sortInt64(out)
	return out
}

// PrimePowerCount returns the number of prime powers p^k (k >= 1) not exceeding
// x. This is the summatory function of the prime-power indicator.
func PrimePowerCount(x int64) int64 {
	return int64(len(PrimePowersUpTo(x)))
}

// Semiprimes returns all semiprimes (products of exactly two primes, counted
// with multiplicity, so squares of primes qualify) that are ≤ x.
func Semiprimes(x int64) []int64 {
	if x < 4 {
		return []int64{}
	}
	var out []int64
	for n := int64(4); n <= x; n++ {
		if BigOmega(n) == 2 {
			out = append(out, n)
		}
	}
	return out
}

// SemiprimeCount returns the number of semiprimes ≤ x.
func SemiprimeCount(x int64) int64 {
	return int64(len(Semiprimes(x)))
}

// SquareFreeCount returns the number of square-free integers in [1, x],
// computed via Σ_{d≤√x} μ(d)⌊x/d²⌋. Its density approaches 6/π².
func SquareFreeCount(x int64) int64 {
	if x < 1 {
		return 0
	}
	root := int64(math.Sqrt(float64(x)))
	for (root+1)*(root+1) <= x {
		root++
	}
	var sum int64
	for d := int64(1); d <= root; d++ {
		sum += int64(MobiusMu(d)) * (x / (d * d))
	}
	return sum
}

// CountPrimesCongruent returns the number of primes p ≤ x with p ≡ a (mod q).
// By Dirichlet's theorem, for gcd(a,q)=1 these are asymptotically equidistributed
// across the φ(q) admissible residue classes.
func CountPrimesCongruent(a, q int64, x int) int64 {
	if q < 1 {
		panic("analyticnt: CountPrimesCongruent requires q >= 1")
	}
	a %= q
	if a < 0 {
		a += q
	}
	var count int64
	for _, p := range PrimesUpTo(x) {
		if p%q == a {
			count++
		}
	}
	return count
}

// PrimesInResidueClass returns the primes p ≤ x with p ≡ a (mod q).
func PrimesInResidueClass(a, q int64, x int) []int64 {
	if q < 1 {
		panic("analyticnt: PrimesInResidueClass requires q >= 1")
	}
	a %= q
	if a < 0 {
		a += q
	}
	var out []int64
	for _, p := range PrimesUpTo(x) {
		if p%q == a {
			out = append(out, p)
		}
	}
	return out
}

// ChebyshevBiasCount compares the number of primes ≤ x in the non-residue class
// versus the residue class modulo q (for q = 4, classes 3 and 1), returning
// their difference N₃ − N₁, which is positive far more often than not (the
// Chebyshev bias).
func ChebyshevBiasCount(x int) int64 {
	return CountPrimesCongruent(3, 4, x) - CountPrimesCongruent(1, 4, x)
}

// InverseLi returns the value x such that LiOffset(x) = y, obtained by Newton
// iteration on Li. It inverts the standard prime-counting estimate and is the
// basis of accurate n-th-prime approximations. y must be > 0.
func InverseLi(y float64) float64 {
	if y <= 0 {
		panic("analyticnt: InverseLi requires y > 0")
	}
	// Initial guess x ≈ y ln y.
	x := y * math.Log(y+2)
	if x < 3 {
		x = 3
	}
	for i := 0; i < 100; i++ {
		f := LiOffset(x) - y
		fp := 1 / math.Log(x)
		nx := x - f/fp
		if nx < 2.001 {
			nx = 2.001
		}
		if math.Abs(nx-x) < 1e-6*x {
			x = nx
			break
		}
		x = nx
	}
	return x
}

// sortInt64 sorts a slice of int64 in place (ascending).
func sortInt64(a []int64) {
	// Simple insertion sort adequate for the short, nearly-sorted slices here.
	for i := 1; i < len(a); i++ {
		v := a[i]
		j := i - 1
		for j >= 0 && a[j] > v {
			a[j+1] = a[j]
			j--
		}
		a[j+1] = v
	}
}
