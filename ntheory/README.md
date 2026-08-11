# ntheory — self-contained number-theory and combinatorics toolkit built entirely on the Go

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/ntheory.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/ntheory)

Package ntheory is a self-contained number-theory and combinatorics toolkit
built entirely on the Go standard library (only math and math/big).

It is a sibling subpackage of github.com/malcolmston/algebra but does not
depend on it: every routine here works directly on machine integers or on
arbitrary-precision math/big values.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/ntheory"
```

## Usage

This is the package's own `ExampleExtendedGCD`, so it compiles and its output is
asserted on every `go test ./ntheory/`.

```go
g, x, y := ntheory.ExtendedGCD(240, 46)
	fmt.Printf("gcd=%d, 240*%d + 46*%d = %d\n", g, x, y, 240*x+46*y)
```

```
gcd=2, 240*-9 + 46*47 = 2
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func AddModU64(a, b, m uint64) uint64` | AddModU64 returns (a + b) mod m without overflowing, even when a + b does not fit in a uint64. |
| `func AllSqrtModComposite(a, m uint64) []uint64` | AllSqrtModComposite returns every square root of a modulo the composite modulus m, sorted ascending. |
| `func Bernoulli(n int64) *big.Rat` | Bernoulli returns the n-th Bernoulli number Bₙ as an exact rational (math/big.Rat), using the convention B₁ = -1/2. |
| `func Binomial(n, k int64) *big.Int` | Binomial returns the binomial coefficient C(n, k) = n! / (k!·(n-k)!) as an arbitrary-precision integer. |
| `func CRT(r, n []int64) (x, mod int64, ok bool)` | CRT solves a system of simultaneous congruences via the Chinese Remainder Theorem. |
| `func Carmichael(n uint64) uint64` | Carmichael returns λ(n), the Carmichael function: the exponent of the multiplicative group modulo n, i.e. |
| `func CatalanNumber(n int64) *big.Int` | CatalanNumber returns the n-th Catalan number Cₙ = C(2n, n) / (n+1) as an arbitrary-precision integer. |
| `func ContinuedFraction(p, q int64) []int64` | ContinuedFraction returns the finite continued-fraction expansion [a0; a1, a2, ...] of the rational p/q (q != 0) as a slice of int64 partial… |
| `func ContinuedFractionRat(r *big.Rat) []int64` | ContinuedFractionRat is the big.Rat form of ContinuedFraction for exact large rationals. |
| `func Convergents(cf []int64) []*big.Rat` | Convergents returns the successive convergents of the partial-quotient slice cf as reduced big.Rat values. |
| `func CountDivisors(n int64) int64` | CountDivisors returns τ(n), the number of positive divisors of n. |
| `func CountDivisorsBig(n *big.Int) *big.Int` | CountDivisorsBig returns tau(\|n\|), the number of positive divisors of \|n\|. |
| `func DiscreteLog(g, h, m int64) (x int64, ok bool)` | DiscreteLog returns the smallest non-negative x with g**x ≡ h (mod m), found via the baby-step giant-step algorithm, together with ok == true. |
| `func DivisorSigma(k, n int64) int64` | DivisorSigma returns σ_k(n), the sum of the k-th powers of the positive divisors of n: Σ_{d \| n} d^k. |
| `func Divisors(n int64) []int64` | Divisors returns all positive divisors of n in ascending order. |
| `func DoubleFactorial(n int64) *big.Int` | DoubleFactorial returns n!! = n·(n-2)·(n-4)·…, the product of every second integer down to 1 or 2. |
| `func EulerPhi(n int64) int64` | EulerPhi returns Euler's totient φ(n): the count of integers in [1, n] that are coprime to n. |
| `func EulerPhiBig(n *big.Int) *big.Int` | EulerPhiBig returns Euler's totient phi(\|n\|): the count of integers in [1, \|n\|] that are coprime to \|n\|. |
| `func ExtendedGCD(a, b int64) (g, x, y int64)` | ExtendedGCD returns g = GCD(a, b) together with Bézout coefficients x and y satisfying a*x + b*y == g. |
| `func Factorial(n int64) *big.Int` | Factorial returns n! as an arbitrary-precision integer. |
| `func Factorize(n int64) map[int64]int` | Factorize returns the prime factorization of \|n\| as a map from each prime factor to its exponent. |
| `func FactorizeU64(n uint64) map[uint64]int` | FactorizeU64 returns the complete prime factorization of n as a map from each prime factor to its (positive) multiplicity. |
| `func Fibonacci(n int64) *big.Int` | Fibonacci returns the n-th Fibonacci number Fₙ, with F₀ = 0, F₁ = 1 and Fₙ = Fₙ₋₁ + Fₙ₋₂, as an arbitrary-precision integer. |
| `func GCD(a, b int64) int64` | GCD returns the greatest common divisor of a and b. |
| `func IsPerfect(n int64) bool` | IsPerfect reports whether n is a perfect number, i.e. |
| `func IsPrime(n int64) bool` | IsPrime reports whether n is prime using a deterministic Miller-Rabin test. |
| `func IsPrimeBig(n *big.Int) bool` | IsPrimeBig reports whether n is (very probably) prime using math/big.Int.ProbablyPrime with 20 Miller-Rabin rounds plus a Baillie-PSW test. |
| `func IsPrimeU64(n uint64) bool` | IsPrimeU64 reports whether n is prime. |
| `func IsPrimitiveRoot(g, m uint64) bool` | IsPrimitiveRoot reports whether g is a primitive root modulo m, i.e. |
| `func IsProbablePrimeBig(n *big.Int) bool` | IsProbablePrimeBig reports whether the arbitrary-precision n is (very probably) prime using a Baillie-PSW test: a strong probable-prime test to base… |
| `func IsQuadraticResidue(a, p int64) bool` | IsQuadraticResidue reports whether a is a quadratic residue modulo the prime p, i.e. |
| `func IsSquare(n int64) bool` | IsSquare reports whether n is a perfect square (n == k*k for some integer k). |
| `func IsqrtBig(n *big.Int) *big.Int` | IsqrtBig returns the integer square root of n: the greatest integer r with r*r <= n, using math/big's exact Sqrt. |
| `func JacobiSymbol(a, n int64) int` | JacobiSymbol returns the Jacobi symbol (a/n) for an odd positive integer n. |
| `func LCM(a, b int64) int64` | LCM returns the least common multiple of a and b. |
| `func LegendreSymbol(a, p int64) int` | LegendreSymbol returns the Legendre symbol (a/p) for an odd prime p: |
| `func Lucas(n int64) *big.Int` | Lucas returns the n-th Lucas number Lₙ, with L₀ = 2, L₁ = 1 and Lₙ = Lₙ₋₁ + Lₙ₋₂, as an arbitrary-precision integer. |
| `func MertensFunction(n int64) int64` | MertensFunction returns M(n) = Σ_{k=1}^{n} μ(k), the Mertens function, the running sum of the Möbius function up to n. |
| `func MillerRabinU64(n, a uint64) bool` | MillerRabinU64 runs a single strong-probable-prime (Miller-Rabin) test of n to base a and reports whether n passes it. |
| `func MobiusMu(n int64) int` | MobiusMu returns the Möbius function μ(n): 0 if n is divisible by a square greater than 1, otherwise (-1)^k where k is the number of distinct prime… |
| `func MobiusSieve(n int64) []int8` | MobiusSieve returns a slice mu of length n+1 where mu`i` == μ(i), the Möbius function of i, for every 0 <= i <= n. |
| `func ModInverse(a, m int64) (inv int64, ok bool)` | ModInverse returns the modular multiplicative inverse of a modulo m: the unique x in [0, m) with a*x ≡ 1 (mod m), together with ok == true. |
| `func ModPow(base, exp, m int64) int64` | ModPow returns base**exp mod m computed with math/big internally so no intermediate value overflows. |
| `func ModPowU64(base, exp, m uint64) uint64` | ModPowU64 returns base**exp mod m by square-and-multiply. |
| `func MulModU64(a, b, m uint64) uint64` | MulModU64 returns a*b mod m computed from the full 128-bit product with math/bits (bits.Mul64 followed by bits.Div64), so no intermediate value… |
| `func Multinomial(ks ...int64) *big.Int` | Multinomial returns the multinomial coefficient (Σkᵢ)! / (k₀!·k₁!·…) as an arbitrary-precision integer. |
| `func MultiplicativeOrder(a, m uint64) (uint64, bool)` | MultiplicativeOrder returns the multiplicative order of a modulo m — the smallest k > 0 with a^k ≡ 1 (mod m) — together with ok == true. |
| `func NextPrime(n int64) int64` | NextPrime returns the smallest prime strictly greater than n. |
| `func NextPrimeBig(n *big.Int) *big.Int` | NextPrimeBig returns the least prime strictly greater than n, using `IsProbablePrimeBig` to test candidates and a mod-30 wheel to skip multiples of… |
| `func NextPrimeU64(n uint64) uint64` | NextPrimeU64 returns the least prime strictly greater than n. |
| `func NthPrime(n uint64) uint64` | NthPrime returns the n-th prime, 1-indexed, so NthPrime(1) == 2, NthPrime(2) == 3 and so on. |
| `func Order(a, m int64) (ord int64, ok bool)` | Order returns the multiplicative order of a modulo m: the smallest positive integer k with a**k ≡ 1 (mod m), together with ok == true. |
| `func Partition(n int64) *big.Int` | Partition returns p(n), the number of distinct ways to write the non-negative integer n as a sum of positive integers where order does not matter, as… |
| `func PellFundamental(n uint64) (x, y *big.Int)` | PellFundamental returns the fundamental solution (x, y) of the Pell equation x^2 - n*y^2 = 1 as arbitrary-precision integers, derived from the… |
| `func Permutations(n, k int64) *big.Int` | Permutations returns the number of ordered k-permutations of n objects, nPr = n! / (n-k)!, as an arbitrary-precision integer. |
| `func PollardBrentU64(n uint64) uint64` | PollardBrentU64 returns a non-trivial factor of the composite integer n using Brent's improved cycle-finding variant of Pollard's rho algorithm. |
| `func PollardRhoBig(n *big.Int) *big.Int` | PollardRhoBig returns a non-trivial factor of the composite integer n using Brent's variant of Pollard's rho over math/big, with the greatest common… |
| `func PrevPrime(n int64) (prime int64, ok bool)` | PrevPrime returns the greatest prime strictly less than n together with ok == true. |
| `func PrevPrimeBig(n *big.Int) (*big.Int, bool)` | PrevPrimeBig returns the greatest prime strictly less than n together with ok == true, using `IsProbablePrimeBig` and a descending mod-30 wheel. |
| `func PrevPrimeU64(n uint64) (uint64, bool)` | PrevPrimeU64 returns the greatest prime strictly less than n together with ok == true. |
| `func PrimePi(n int64) int64` | PrimePi returns π(n), the number of primes p with p <= n. |
| `func PrimePiRange(lo, hi uint64) uint64` | PrimePiRange counts the primes p with lo <= p <= hi without materializing them, sieving one segment at a time and tallying survivors. |
| `func PrimesInRange(lo, hi uint64) []uint64` | PrimesInRange returns all primes in the half-open interval [lo, hi) in ascending order. |
| `func PrimesUpTo(n int64) []int64` | PrimesUpTo returns all primes p with p <= n in ascending order, computed with the sieve of Eratosthenes. |
| `func PrimitiveRoot(m uint64) (uint64, bool)` | PrimitiveRoot returns the smallest primitive root modulo m, together with ok == true. |
| `func PrimitiveRoots(m uint64) []uint64` | PrimitiveRoots returns all primitive roots modulo m in ascending order, or nil when none exist. |
| `func Radical(n int64) int64` | Radical returns the radical of n: the product of its distinct prime factors. |
| `func RatFromContinuedFraction(cf []int64) *big.Rat` | RatFromContinuedFraction reconstructs the exact big.Rat value of the partial-quotient slice cf, i.e. |
| `func SegmentedSieve(lo, hi uint64) []uint64` | SegmentedSieve returns all primes p with lo <= p <= hi in ascending order. |
| `func SqrtContinuedFraction(n uint64) (a0 int64, period []int64)` | SqrtContinuedFraction returns the continued-fraction expansion of the square root of n for a positive n. |
| `func SqrtMod(a, p int64) (x int64, ok bool)` | SqrtMod returns a square root x in [0, p) of a modulo the prime p, i.e. |
| `func SqrtModBig(a, p *big.Int) (*big.Int, bool)` | SqrtModBig returns a square root of a modulo the prime p using the Tonelli-Shanks algorithm over math/big, suitable for primes far larger than int64. |
| `func SqrtModPrimePower(a, p int64, k int) (x int64, ok bool)` | SqrtModPrimePower returns a square root x in [0, p^k) of a modulo the prime power p^k, together with ok == true, or (0, false) when a is a… |
| `func StirlingSecond(n, k int64) *big.Int` | StirlingSecond returns the Stirling number of the second kind S(n, k): the number of ways to partition a set of n labeled elements into k non-empty… |
| `func SubModU64(a, b, m uint64) uint64` | SubModU64 returns (a - b) mod m without underflowing. |
| `func SumDivisors(n int64) int64` | SumDivisors returns σ(n), the sum of all positive divisors of n (including 1 and n itself). |
| `func TotientSieve(n int64) []int64` | TotientSieve returns a slice phi of length n+1 where phi`i` == φ(i), Euler's totient of i, for every 0 <= i <= n. |
| `func Tribonacci(n int64) *big.Int` | Tribonacci returns the n-th Tribonacci number Tₙ, with T₀ = 0, T₁ = 1, T₂ = 1 and Tₙ = Tₙ₋₁ + Tₙ₋₂ + Tₙ₋₃, as an… |

### Types

| Type | What it is |
| --- | --- |
| `Barrett` | Barrett holds a precomputed reciprocal for repeated reduction by a fixed modulus. |
| `Montgomery` | Montgomery holds the precomputed constants for Montgomery reduction with a fixed odd modulus. |
| `PrimePower` | PrimePower pairs a prime with its exponent in a factorization. |
| `PrimePowerBig` | PrimePowerBig pairs an arbitrary-precision prime with its exponent in a factorization. |
| `PrimePowerU64` | PrimePowerU64 pairs a prime with its exponent in a uint64 factorization. |
| `PrimeSieve` | PrimeSieve is a stateful generator that yields the primes in ascending order, one per call to `PrimeSieve.Next`, without ever holding them all in… |

<details>
<summary><code>Barrett</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewBarrett(m uint64) *Barrett` | NewBarrett returns a Barrett context for the modulus m. |
| `func (b *Barrett) Modulus() uint64` | Modulus returns the modulus this context reduces by. |
| `func (b *Barrett) MulMod(a, c uint64) uint64` | MulMod returns (a*c) mod m for the context's modulus m. |
| `func (b *Barrett) Reduce(x uint64) uint64` | Reduce returns x mod m for the context's modulus m. |

</details>

<details>
<summary><code>Montgomery</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewMontgomery(m uint64) *Montgomery` | NewMontgomery returns a Montgomery context for the odd modulus m. |
| `func (mo *Montgomery) FromMont(a uint64) uint64` | FromMont maps a Montgomery-domain value back to normal form, returning the Montgomery reduction a*R**-1 mod n. |
| `func (mo *Montgomery) Modulus() uint64` | Modulus returns the odd modulus this context reduces by. |
| `func (mo *Montgomery) MulMont(a, b uint64) uint64` | MulMont returns REDC(a*b), the Montgomery product of a and b. |
| `func (mo *Montgomery) PowMont(base, exp uint64) uint64` | PowMont returns base**exp mod n. |
| `func (mo *Montgomery) ToMont(a uint64) uint64` | ToMont maps a into the Montgomery domain, returning a*R mod n. |

</details>

<details>
<summary><code>PrimePower</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func FactorList(n int64) []PrimePower` | FactorList returns the prime factorization of \|n\| as a slice of `PrimePower` values sorted by ascending prime. |

</details>

<details>
<summary><code>PrimePowerBig</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func FactorizeBig(n *big.Int) []PrimePowerBig` | FactorizeBig returns the prime factorization of \|n\| as a slice of `PrimePowerBig` values sorted by ascending prime. |

</details>

<details>
<summary><code>PrimePowerU64</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func FactorListUint64(n uint64) []PrimePowerU64` | FactorListUint64 returns the prime factorization of n as a slice of `PrimePowerU64` sorted by ascending prime. |

</details>

<details>
<summary><code>PrimeSieve</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewPrimeSieve() *PrimeSieve` | NewPrimeSieve returns a PrimeSieve positioned before the first prime, so the first call to Next returns 2. |
| `func (s *PrimeSieve) Next() uint64` | Next returns the next prime in ascending order, advancing to a further segment only when the current one is exhausted. |

</details>

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/ntheory).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
