# analyticnt — self-contained analytic number theory toolkit built entirely on the Go standard library

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/analyticnt.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/analyticnt)

Package analyticnt is a self-contained analytic number theory toolkit built
entirely on the Go standard library.

It is a sibling subpackage of github.com/malcolmston/algebra but does not
depend on it. Everything here is implemented directly on machine floats,
integers, and math/big / math/cmplx values.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/analyticnt"
```

## Usage

This is the package's own `ExampleMobiusMu`, so it compiles and its output is
asserted on every `go test ./analyticnt/`.

```go
fmt.Println(ant.MobiusMu(30), ant.MobiusMu(12))
```

```
-1 0
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func AverageGap(x int) float64` | AverageGap returns the mean gap between consecutive primes up to x, which by the prime number theorem is asymptotic to ln x. |
| `func Bernoulli(n int) float64` | Bernoulli returns the n-th Bernoulli number B_n as a float64, computed by the Akiyama–Tanigawa algorithm. |
| `func BertrandPrime(n int64) int64` | BertrandPrime returns a prime p with n < p < 2n for n >= 1 (guaranteed to exist by Bertrand's postulate). |
| `func BigOmega(n int64) int` | BigOmega returns Ω(n), the number of prime factors of n counted with multiplicity. |
| `func CarmichaelLambda(n int64) int64` | CarmichaelLambda returns the Carmichael function λ(n), the exponent of the group (Z/nZ)*: the least k >= 1 with a^k ≡ 1 (mod n) for every a… |
| `func ChebyshevBiasCount(x int) int64` | ChebyshevBiasCount compares the number of primes ≤ x in the non-residue class versus the residue class modulo q (for q = 4, classes 3 and 1),… |
| `func ChebyshevPsi(x float64) float64` | ChebyshevPsi returns the second Chebyshev function ψ(x) = Σ_{n ≤ x} Λ(n) = Σ_{p^k ≤ x} ln p. |
| `func ChebyshevTheta(x float64) float64` | ChebyshevTheta returns the first Chebyshev function θ(x) = Σ_{p ≤ x} ln p, summing the natural logarithm over primes p at most x. |
| `func ChebyshevThetaExp(x float64) float64` | ChebyshevThetaExp returns e^{θ(x)} = Π_{p ≤ x} p, the primorial of x, as a float64. |
| `func ClassNumberDirichlet(p int64) int64` | ClassNumberDirichlet returns the class number h(−p) of the imaginary quadratic field Q(√−p) for a prime p ≡ 3 (mod 4), computed from the… |
| `func CompletedZeta(s complex128) complex128` | CompletedZeta returns the completed zeta function ξ̂(s) = π^{−s/2} Γ(s/2) ζ(s), which satisfies the symmetric functional equation ξ̂(s) =… |
| `func CountPrimesCongruent(a, q int64, x int) int64` | CountPrimesCongruent returns the number of primes p ≤ x with p ≡ a (mod q). |
| `func CousinPrimesUpTo(x int) [][2]int64` | CousinPrimesUpTo returns all cousin-prime pairs (p, p+4) with p+4 ≤ x. |
| `func CramerGap(p float64) float64` | CramerGap returns Cramér's heuristic estimate (ln p)² for the gap following a prime p, the conjectured order of the maximal prime gap near p. |
| `func Digamma(x float64) float64` | Digamma returns the digamma function ψ(x) = Γ'(x)/Γ(x) for real x > 0, using the recurrence to shift the argument upward followed by the standard… |
| `func DirichletBeta(s float64) float64` | DirichletBeta returns the Dirichlet beta function β(s) = Σ_{n≥0} (−1)^n (2n+1)^{-s}, the L-function of the non-principal character modulo 4. |
| `func DirichletEta(s float64) float64` | DirichletEta returns the Dirichlet eta function η(s) = Σ (-1)^{n-1} n^{-s} = (1 − 2^{1−s}) ζ(s) for real s. |
| `func DirichletEtaComplex(s complex128) complex128` | DirichletEtaComplex returns η(s) = (1 − 2^{1−s}) ζ(s) for complex s, with the removable value η(1) = ln 2. |
| `func DirichletL(s complex128, chi DirichletCharacter) complex128` | DirichletL returns the Dirichlet L-function L(s, χ) for a complex argument s (with s ≠ 1 when χ is principal), using the Hurwitz-zeta… |
| `func DirichletLSeries(s float64, chi DirichletCharacter, terms int) complex128` | DirichletLSeries returns the truncated Dirichlet series Σ_{n=1}^{terms} χ(n) n^{-s} for real s > 1. |
| `func DirichletLambda(s float64) float64` | DirichletLambda returns the Dirichlet lambda function λ(s) = Σ_{n≥0} (2n+1)^{-s} = (1 − 2^{-s}) ζ(s) for real s > 1. |
| `func DivisorCount(n int64) int64` | DivisorCount returns τ(n) = d(n), the number of positive divisors of n. |
| `func DivisorSigma(k int, n int64) int64` | DivisorSigma returns σ_k(n), the sum of the k-th powers of the divisors of n. |
| `func E1(x float64) float64` | E1 returns the exponential integral E1(x) = ∫_x^∞ e^{-t}/t dt for x > 0. |
| `func Ei(x float64) float64` | Ei returns the exponential integral Ei(x) = -PV ∫_{-x}^∞ e^{-t}/t dt, the principal value for x > 0 and the ordinary integral for x < 0. |
| `func EulerGammaExp() float64` | EulerGammaExp returns e^{−γ} ≈ 0.5614594836, the constant appearing in Mertens' third theorem. |
| `func EulerPhi(n int64) int64` | EulerPhi returns Euler's totient φ(n), the number of integers in [1, n] coprime to n. |
| `func ExplicitPsi(x float64, numZeros int) float64` | ExplicitPsi returns Riemann's explicit-formula approximation to ψ(x) using the first numZeros nontrivial zeros ρ = 1/2 + iγ: ψ(x) ≈ x − Σ_ρ… |
| `func FactorInt(n int64) map[int64]int` | FactorInt returns the prime factorization of n >= 1 as a map from prime to exponent. |
| `func FactorialMod(n, m int64) int64` | FactorialMod returns n! mod m using iterative multiplication. |
| `func GCD(a, b int64) int64` | GCD returns the greatest common divisor of a and b (non-negative). |
| `func Gamma(x float64) float64` | Gamma returns the gamma function Γ(x) for real x, delegating to the standard library implementation. |
| `func GaussSum(chi DirichletCharacter) complex128` | GaussSum returns the Gauss sum τ(χ) = Σ_{a=1}^{q} χ(a) e^{2πi a/q} associated with a Dirichlet character χ modulo q. |
| `func GoldbachPair(n int) (p, q int64, ok bool)` | GoldbachPair returns a pair of primes summing to the even integer n, choosing the pair with the smallest first prime. |
| `func GoldbachPartitions(n int) int64` | GoldbachPartitions returns the number of ways to write the even integer n as an unordered sum of two primes p + q with p ≤ q (the Goldbach… |
| `func GramPoint(n int) float64` | GramPoint returns the n-th Gram point g_n, the solution of θ(g_n) = nπ, found by Newton iteration. |
| `func HardyLittlewoodTwinEstimate(x float64) float64` | HardyLittlewoodTwinEstimate returns the Hardy–Littlewood asymptotic estimate 2·C₂·x/(ln x)² for the number of twin primes up to x. |
| `func HardyZ(t float64) float64` | HardyZ returns Hardy's Z-function Z(t) = e^{iθ(t)} ζ(1/2 + it), which is real for real t and shares its real zeros with ζ on the critical line. |
| `func HurwitzZeta(s, a float64) float64` | HurwitzZeta returns the Hurwitz zeta function ζ(s, a) = Σ_{n≥0} (n+a)^{-s} for real s > 1 and a > 0, via Euler–Maclaurin summation (which also… |
| `func HurwitzZetaComplex(s complex128, a float64) complex128` | HurwitzZetaComplex returns ζ(s, a) for complex s (s ≠ 1) and real a > 0. |
| `func InverseLi(y float64) float64` | InverseLi returns the value x such that LiOffset(x) = y, obtained by Newton iteration on Li. |
| `func IsCoprime(a, b int64) bool` | IsCoprime reports whether a and b are coprime (gcd 1). |
| `func IsPrime(n int64) bool` | IsPrime reports whether n is prime using deterministic trial division backed by a 6k±1 wheel. |
| `func IsPrimePower(n int64) (ok bool, p int64, k int)` | IsPrimePower reports whether n = p^k for some prime p and k >= 1, and if so returns the base prime p and exponent k. |
| `func IsSafePrime(p int64) bool` | IsSafePrime reports whether p is a safe prime, i.e. |
| `func IsSophieGermain(p int64) bool` | IsSophieGermain reports whether p and 2p+1 are both prime (p is a Sophie Germain prime). |
| `func IsSquareFree(n int64) bool` | IsSquareFree reports whether n has no repeated prime factor. |
| `func IsWilsonPrime(p int64) bool` | IsWilsonPrime reports whether p is a Wilson prime, i.e. |
| `func JacobiSymbol(a, n int64) int` | JacobiSymbol returns the Jacobi symbol (a\|n) for odd n >= 1, generalizing the Legendre symbol multiplicatively over the prime factorization of n. |
| `func JordanTotient(k int, n int64) int64` | JordanTotient returns the Jordan totient J_k(n) = n^k Π_{p\|n}(1 − p^{-k}), generalizing Euler's totient (k = 1). |
| `func KroneckerSymbol(a, n int64) int` | KroneckerSymbol returns the Kronecker symbol (a\|n), the completely multiplicative extension of the Jacobi symbol to every integer n, including even… |
| `func LCM(a, b int64) int64` | LCM returns the least common multiple of a and b. |
| `func LFunctionReal(s float64, chi DirichletCharacter) complex128` | LFunctionReal returns L(s, χ) for a real argument s, returning the value as a complex128 (which is real for real characters). |
| `func LegendrePhi(x int64, a int) int64` | LegendrePhi returns Legendre's φ(x, a): the count of integers in [1, x] not divisible by any of the first a primes. |
| `func LegendreSymbol(a, p int64) int` | LegendreSymbol returns the Legendre symbol (a\|p) for an odd prime p using Euler's criterion, evaluated as a^{(p−1)/2} mod p mapped to {−1, 0, 1}. |
| `func Li(x float64) float64` | Li returns the logarithmic integral li(x) = PV ∫_0^x dt/ln t for x > 0, x ≠ 1. |
| `func LiOffset(x float64) float64` | LiOffset returns the offset logarithmic integral Li(x) = li(x) − li(2) = ∫_2^x dt/ln t. |
| `func LiValueAt2() float64` | LiValueAt2 returns li(2) ≈ 1.045163780117493. |
| `func Liouville(n int64) int` | Liouville returns the Liouville function λ(n) = (−1)^{Ω(n)}, where Ω counts prime factors with multiplicity. |
| `func LogGamma(x float64) float64` | LogGamma returns the natural logarithm of \|Γ(x)\| for real x. |
| `func LogGammaComplex(z complex128) complex128` | LogGammaComplex returns log Γ(z) (principal branch) for complex z with positive real part sufficiently large, using the Lanczos approximation. |
| `func MangoldtSummatory(N int64) float64` | MangoldtSummatory returns Σ_{n≤N} Λ(n)/n as a partial sum related to the logarithmic derivative of ζ; it is a small utility for PNT experiments. |
| `func MaxPrimeGap(x int) (gap int64, startingPrime int64)` | MaxPrimeGap returns the largest gap between consecutive primes up to x, along with the prime that begins it. |
| `func MaximalPrimeGaps(x int) [][2]int64` | MaximalPrimeGaps returns the record ("maximal") prime gaps up to x: each entry is a gap strictly larger than every earlier gap, paired with the prime… |
| `func MeritOfGap(g int64, p int64) float64` | MeritOfGap returns the merit of a prime gap g following prime p, defined as g/ln p; merits above about 30 correspond to exceptionally large gaps. |
| `func MertensFunction(n int) int64` | MertensFunction returns the Mertens function M(n) = Σ_{k=1}^{n} μ(k). |
| `func MertensProduct(x float64) float64` | MertensProduct returns Π_{p ≤ x} (1 − 1/p). |
| `func MertensSequence(n int) []int64` | MertensSequence returns the running Mertens values M(0), M(1), …, M(n) as a slice of length n+1. |
| `func MobiusMu(n int64) int` | MobiusMu returns the Möbius function μ(n): 1 if n is a square-free positive integer with an even number of prime factors, −1 if square-free with… |
| `func MobiusSieve(n int) []int` | MobiusSieve returns μ(k) for every k from 0 to n as a slice of length n+1, with entry 0 left as 0. |
| `func MobiusSummatoryReciprocal(N int) float64` | MobiusSummatoryReciprocal returns Σ_{n≤N} μ(n)/n, which tends to 0 and whose value is bounded by 1 in absolute value. |
| `func MultiplicativeOrder(a, n int64) int64` | MultiplicativeOrder returns the multiplicative order of a modulo n: the least k >= 1 with a^k ≡ 1 (mod n). |
| `func NextPrime(n int64) int64` | NextPrime returns the smallest prime strictly greater than n. |
| `func NthPrime(n int) int64` | NthPrime returns the n-th prime (1-indexed): NthPrime(1) == 2, NthPrime(2) == 3, and so on. |
| `func NthPrimeApprox(n int) float64` | NthPrimeApprox returns an asymptotic estimate of the n-th prime using the expansion p_n ≈ n(ln n + ln ln n − 1 + (ln ln n − 2)/ln n). |
| `func NthPrimeInverseR(n int) float64` | NthPrimeInverseR returns the n-th prime estimated by inverting the Riemann R function, i.e. |
| `func Omega(n int64) int` | Omega returns ω(n), the number of distinct prime factors of n. |
| `func PiError(x int64) float64` | PiError returns the signed error li(x) − π(x) of the offset-logarithmic integral estimate of the prime-counting function at integer x. |
| `func PiErrorRiemann(x int64) float64` | PiErrorRiemann returns the signed error R(x) − π(x) of the Riemann estimate. |
| `func PrevPrime(n int64) int64` | PrevPrime returns the largest prime strictly less than n, or 0 if there is none (n <= 2). |
| `func PrimeCountBetween(lo, hi int64) int64` | PrimeCountBetween returns the number of primes p with lo < p <= hi. |
| `func PrimeCountingBias(x int64) (bias float64, liOverestimates bool)` | PrimeCountingBias returns li(x) − π(x) as a float and reports whether the value is positive (li over-estimates π), illustrating the Chebyshev… |
| `func PrimeDensity(x float64) float64` | PrimeDensity returns the prime density 1/ln x, the probability that a number near x is prime according to the prime number theorem. |
| `func PrimeFactors(n int64) []int64` | PrimeFactors returns the distinct prime factors of n in increasing order. |
| `func PrimeFactorsWithMultiplicity(n int64) []int64` | PrimeFactorsWithMultiplicity returns the prime factors of n listed with repetition, in increasing order, so that their product is n. |
| `func PrimeGapHistogram(x int) map[int64]int` | PrimeGapHistogram returns a map from gap size to its frequency among consecutive primes up to x. |
| `func PrimeGaps(x int) []int64` | PrimeGaps returns the sequence of gaps p_{n+1} − p_n between consecutive primes up to x. |
| `func PrimeHarmonic(x float64) float64` | PrimeHarmonic returns Σ_{p ≤ x} 1/p, the sum of reciprocals of primes up to x. |
| `func PrimeIndex(p int64) int64` | PrimeIndex returns the 1-based index of prime p in the sequence of primes (so PrimeIndex(2) == 1). |
| `func PrimeLogHarmonic(x float64) float64` | PrimeLogHarmonic returns Σ_{p ≤ x} (ln p)/p, which by Mertens' first theorem is asymptotic to ln x. |
| `func PrimePi(x int64) int64` | PrimePi returns the exact prime-counting function pi(x): the number of primes p with p <= x. |
| `func PrimePiApprox(x float64) float64` | PrimePiApprox returns the recommended real-valued estimate of π(x), namely the Riemann R function for x >= 2 and 0 below that. |
| `func PrimePiCheck(x int64) (value int64, agree bool)` | PrimePiCheck cross-checks all three exact counters (sieve, Legendre, Meissel) at x and reports whether they agree, returning the common value and the… |
| `func PrimePiGauss(x float64) float64` | PrimePiGauss returns Gauss's original estimate x/ln x of π(x). |
| `func PrimePiLegendre(x int64) int64` | PrimePiLegendre returns π(x) via Legendre's formula π(x) = φ(x, a) + a − 1, where a = π(⌊√x⌋). |
| `func PrimePiLegendreApprox(x float64) float64` | PrimePiLegendreApprox returns Legendre's classical estimate x/(ln x − B) with B = 1.08366 (Legendre's constant). |
| `func PrimePiLi(x float64) float64` | PrimePiLi returns the offset logarithmic-integral estimate Li(x) of π(x) for x >= 2. |
| `func PrimePiMeissel(x int64) int64` | PrimePiMeissel returns π(x) via the Meissel formula, which reduces the work relative to Legendre's method. |
| `func PrimePiSieve(n int) []int64` | PrimePiSieve returns pi(k) for every k from 0 to n as a cumulative slice of length n+1, where result`k` is the number of primes <= k. |
| `func PrimePowerCount(x int64) int64` | PrimePowerCount returns the number of prime powers p^k (k >= 1) not exceeding x. |
| `func PrimePowersUpTo(x int64) []int64` | PrimePowersUpTo returns every prime power p^k (k >= 1) that is ≤ x, in increasing order. |
| `func PrimeTripletsUpTo(x int) [][3]int64` | PrimeTripletsUpTo returns prime triplets (p, p+2, p+6) or (p, p+4, p+6) — the two admissible constellations of diameter 6 — with the largest… |
| `func PrimeZeta(s float64) float64` | PrimeZeta returns the prime zeta function P(s) = Σ_p p^{-s} for real s > 1, computed by Möbius inversion P(s) = Σ_{n≥1} μ(n)/n · ln ζ(ns). |
| `func PrimesInRange(lo, hi int) []int64` | PrimesInRange returns all primes p with lo <= p <= hi in increasing order. |
| `func PrimesInResidueClass(a, q int64, x int) []int64` | PrimesInResidueClass returns the primes p ≤ x with p ≡ a (mod q). |
| `func PrimesUpTo(n int) []int64` | PrimesUpTo returns all primes p with p <= n in increasing order. |
| `func PrimitiveRoot(n int64) (int64, bool)` | PrimitiveRoot returns a primitive root modulo n (a generator of the cyclic group (Z/nZ)*) together with true, or (0, false) when n has no primitive… |
| `func PrimorialPrimes(k int) []int64` | PrimorialPrimes returns the first k primes whose product is the k-th primorial. |
| `func PsiError(x float64) float64` | PsiError returns the normalized error (ψ(x) − x)/√x of the second Chebyshev function, the quantity whose boundedness is equivalent to the… |
| `func PsiSchoenfeldBound(x float64) float64` | PsiSchoenfeldBound returns Schoenfeld's RH-conditional bound \|ψ(x) − x\| < (1/(8π))·√x·(ln x)², valid for x >= 73.2. |
| `func Radical(n int64) int64` | Radical returns rad(n), the product of the distinct prime factors of n. |
| `func RandomPrimeBelow(n int, seed int64) int64` | RandomPrimeBelow returns a uniformly random prime strictly less than n using a deterministic math/rand source seeded by seed. |
| `func RiemannHypothesisBound(x float64) float64` | RiemannHypothesisBound returns the RH-conditional error bound for \|π(x) − li(x)\|, namely (1/(8π))·√x·ln x, which holds for x >= 2657 under… |
| `func RiemannR(x float64) float64` | RiemannR returns the Riemann prime-counting approximation R(x), evaluated via the Gram series R(x) = 1 + Σ_{k≥1} (ln x)^k / (k · k! · ζ(k+1)). |
| `func RiemannRPrime(x float64) float64` | RiemannRPrime returns the derivative R'(x) via the series R'(x) = Σ_{k≥1} (ln x)^{k-1} / (x · k! · ζ(k+1)), useful for inverting R. |
| `func RiemannRScaled(x float64) float64` | RiemannRScaled returns R(x) computed with the Möbius-inversion form R(x) = Σ_{n≥1} μ(n)/n · li(x^{1/n}). |
| `func RiemannSiegelTheta(t float64) float64` | RiemannSiegelTheta returns the Riemann–Siegel theta function θ(t) = arg Γ(1/4 + it/2) − (t/2) ln π, real for real t. |
| `func RiemannSiegelThetaAsymptotic(t float64) float64` | RiemannSiegelThetaAsymptotic returns the standard asymptotic expansion of θ(t) for large t; it is faster than RiemannSiegelTheta but only accurate… |
| `func RiemannSiegelZ(t float64) float64` | RiemannSiegelZ is a synonym for HardyZ, the real-valued function whose sign changes locate the nontrivial zeros of ζ on the critical line. |
| `func RiemannVonMangoldtN(T float64) float64` | RiemannVonMangoldtN returns the counting function N(T), the approximate number of nontrivial zeros of ζ with imaginary part in (0, T], given by the… |
| `func RiemannXiComplex(s complex128) complex128` | RiemannXiComplex returns the Riemann ξ-function ξ(s) = ½ s(s−1) π^{−s/2} Γ(s/2) ζ(s), an entire function whose zeros are exactly the… |
| `func RiemannXiReal(s float64) complex128` | RiemannXiReal returns ξ(s) for a real argument s as a complex128 (the value is real). |
| `func SafePrimesUpTo(x int) []int64` | SafePrimesUpTo returns all safe primes q ≤ x, those for which (q−1)/2 is prime. |
| `func SecondChebyshev(x float64) float64` | SecondChebyshev is a synonym for ChebyshevPsi. |
| `func SemiprimeCount(x int64) int64` | SemiprimeCount returns the number of semiprimes ≤ x. |
| `func Semiprimes(x int64) []int64` | Semiprimes returns all semiprimes (products of exactly two primes, counted with multiplicity, so squares of primes qualify) that are ≤ x. |
| `func SexyPrimesUpTo(x int) [][2]int64` | SexyPrimesUpTo returns all sexy-prime pairs (p, p+6) with p+6 ≤ x. |
| `func Sieve(n int) []bool` | Sieve returns a boolean slice s of length n+1 where s`i` reports whether i is prime, for 0 <= i <= n. |
| `func SkewesNumberLog() float64` | SkewesNumberLog returns the natural logarithm of the smallest x at which li(x) − π(x) is believed to change sign. |
| `func SophieGermainPrimesUpTo(x int) []int64` | SophieGermainPrimesUpTo returns all Sophie Germain primes p (p and 2p+1 both prime) with p ≤ x. |
| `func SquareFreeCount(x int64) int64` | SquareFreeCount returns the number of square-free integers in [1, x], computed via Σ_{d≤√x} μ(d)⌊x/d²⌋. |
| `func SquareFreeDensityEstimate(n float64) float64` | SquareFreeDensityEstimate returns the asymptotic count 6n/π² = n/ζ(2) of square-free integers in [1, n]. |
| `func ThetaError(x float64) float64` | ThetaError returns the normalized error (θ(x) − x)/√x of the first Chebyshev function. |
| `func TotientSummatory(n int) int64` | TotientSummatory returns Φ(n) = Σ_{k=1}^{n} φ(k), whose asymptotic is 3n²/π². |
| `func Trigamma(x float64) float64` | Trigamma returns the trigamma function ψ₁(x) = d²/dx² ln Γ(x) for real x > 0. |
| `func TwinPrimeCount(x int) int64` | TwinPrimeCount returns the number of twin-prime pairs (p, p+2) with p+2 ≤ x. |
| `func TwinPrimeDensity(x float64) float64` | TwinPrimeDensity returns the Hardy–Littlewood estimate of the number of twin prime pairs (p, p+2) with p ≤ x, namely 2·C₂·∫_2^x dt/(ln t)². |
| `func TwinPrimesUpTo(x int) [][2]int64` | TwinPrimesUpTo returns all twin-prime pairs (p, p+2) with p+2 ≤ x, as a slice of [p, p+2] pairs. |
| `func VonMangoldt(n int64) float64` | VonMangoldt returns the von Mangoldt function Λ(n): ln p if n = p^k for a prime p and integer k >= 1, and 0 otherwise. |
| `func VonMangoldtLambda(n int64) float64` | VonMangoldtLambda is a synonym for VonMangoldt. |
| `func VonMangoldtSum(N int64) float64` | VonMangoldtSum returns Σ_{n ≤ N} Λ(n), which equals ψ(N). |
| `func Zeta(s float64) float64` | Zeta returns the Riemann zeta function ζ(s) for a real argument s ≠ 1, using the Euler–Maclaurin summation formula, which analytically continues… |
| `func ZetaAtEvenInteger(k int) float64` | ZetaAtEvenInteger returns the exact value ζ(2k) = (−1)^{k+1} B_{2k} (2π)^{2k} / (2·(2k)!) for k >= 1, using tabulated Bernoulli numbers. |
| `func ZetaComplex(s complex128) complex128` | ZetaComplex returns ζ(s) for a complex argument s ≠ 1 via Euler–Maclaurin summation. |
| `func ZetaPrime(s float64) float64` | ZetaPrime returns the derivative ζ'(s) for real s ≠ 1 via a central finite difference of the Euler–Maclaurin evaluation. |
| `func ZetaZero(n int) float64` | ZetaZero returns the imaginary part of the n-th nontrivial zero of ζ on the critical line (n >= 1), located by bracketing sign changes of HardyZ and… |
| `func ZetaZeroRealPart() float64` | ZetaZeroRealPart returns the common real part 1/2 of every nontrivial zero on the critical line; it is a convenience constant-returning helper. |
| `func ZetaZeros(n int) []float64` | ZetaZeros returns the imaginary parts of the first n nontrivial zeros of ζ on the critical line, in increasing order. |

### Types

| Type | What it is |
| --- | --- |
| `DirichletCharacter` | DirichletCharacter represents a Dirichlet character χ modulo q: a completely multiplicative, q-periodic function on the integers, defined by its… |

<details>
<summary><code>DirichletCharacter</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func AllCharactersModulo(q int64) []DirichletCharacter` | AllCharactersModulo returns the full group of φ(q) Dirichlet characters modulo q when (Z/qZ)* is cyclic. |
| `func CharacterFromPrimitiveRoot(q int64, j int) DirichletCharacter` | CharacterFromPrimitiveRoot builds a Dirichlet character modulo q with cyclic group (Z/qZ)* by sending a fixed primitive root g to e^{2πi·j/φ(q)}. |
| `func LegendreCharacter(p int64) DirichletCharacter` | LegendreCharacter returns the quadratic character modulo an odd prime p given by the Legendre symbol (a\|p): +1 for quadratic residues, −1 for… |
| `func PrincipalCharacter(q int64) DirichletCharacter` | PrincipalCharacter returns the principal Dirichlet character χ₀ modulo q. |
| `func (c DirichletCharacter) Conjugate() DirichletCharacter` | Conjugate returns the complex-conjugate character χ̄, which is also the inverse character χ^{-1}. |
| `func (c DirichletCharacter) Eval(n int64) complex128` | Eval returns χ(n), reducing n modulo q. |
| `func (c DirichletCharacter) IsPrincipal() bool` | IsPrincipal reports whether χ is the principal character χ₀ modulo q, which is 1 on residues coprime to q and 0 otherwise. |
| `func (c DirichletCharacter) IsReal() bool` | IsReal reports whether χ takes only real values (±1, 0), i.e. |
| `func (c DirichletCharacter) Order() int` | Order returns the multiplicative order of χ in the character group, i.e. |

</details>

### Constants

`MertensConstant`, `TwinPrimeConstant`, `BrunConstant`, `ArtinConstant`, `LandauRamanujanConstant`, `FellerTornierConstant`, `SarnakConstant`, `GlaisherKinkelin`, `CatalanConstant`, `StieltjesGamma0`, `StieltjesGamma1`, `StieltjesGamma2`, `LegendreConstant`, `AperyConstant`, `EulerGamma`, `SoldnerConstant`

### Variables

`ErrNonPositive`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/analyticnt).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
