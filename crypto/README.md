# crypto — implements the mathematical primitives that underpin classical public-key cryptography

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/crypto.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/crypto)

Package crypto implements the mathematical primitives that underpin classical
public-key cryptography, built exclusively on top of math/big.Int and the Go
standard library.

The package is deliberately educational and self-contained: it exposes the
number-theoretic machinery directly rather than hiding it behind opaque,
hardened implementations. Nothing here draws on crypto/rand, crypto/rsa or any
other package under crypto/*; every routine is expressed in terms of
arbitrary-precision integer arithmetic so that the algorithms can be read,
tested against closed-form values and reused for teaching.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/crypto"
```

## Usage

This is the package's own `ExampleBabyStepGiantStep`, so it compiles and its output is
asserted on every `go test ./crypto/`.

```go
x, err := crypto.DiscreteLog(big.NewInt(3), big.NewInt(13), big.NewInt(17))
	fmt.Println(x, err)
	fmt.Println(crypto.ModExp(big.NewInt(3), x, big.NewInt(17)))
```

```
4 <nil>
13
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func BabyStepGiantStep(g, h, p, order *big.Int) (*big.Int, error)` | BabyStepGiantStep solves the discrete logarithm g^x ≡ h (mod p) for x in [0, order), using Shanks' baby-step giant-step algorithm in O(sqrt(order))… |
| `func CRT(residues, moduli []*big.Int) (*big.Int, error)` | CRT solves a system of simultaneous congruences x ≡ residues`i` (mod moduli`i`) via the Chinese Remainder Theorem and returns the unique solution x… |
| `func CRTPair(a1, n1, a2, n2 *big.Int) (*big.Int, error)` | CRTPair solves the pair of congruences x ≡ a1 (mod n1) and x ≡ a2 (mod n2) and returns the unique solution in [0, n1*n2). |
| `func CarmichaelLambda(n *big.Int) *big.Int` | CarmichaelLambda returns the Carmichael function λ(n): the smallest positive integer m such that a^m ≡ 1 (mod n) for every a coprime to n. |
| `func CombineShares(shares []Share, prime *big.Int) (*big.Int, error)` | CombineShares reconstructs the secret from a set of shares using Lagrange interpolation at x = 0 over the prime field Z/prime. |
| `func DiscreteLog(g, h, p *big.Int) (*big.Int, error)` | DiscreteLog solves g^x ≡ h (mod p) for a prime modulus p, searching the full range [0, p-1). |
| `func Divisors(n *big.Int) []*big.Int` | Divisors returns every positive divisor of n in ascending order. |
| `func ElGamalDecrypt(priv *ElGamalPrivateKey, c1, c2 *big.Int) (*big.Int, error)` | ElGamalDecrypt recovers the message m from a ciphertext pair (c1, c2) using the private key, computing m = c2 · (c1^x)^-1 mod P. It returns an error… |
| `func ElGamalEncrypt(pub ElGamalPublicKey, m *big.Int, rng *rand.Rand) (c1, c2 *big.Int, err error)` | ElGamalEncrypt encrypts the integer message m (0 <= m < P) under the public key, using the ephemeral randomness k drawn from rng. |
| `func ElGamalVerify(pub ElGamalPublicKey, m *big.Int, sig ElGamalSignature) bool` | ElGamalVerify reports whether sig is a valid ElGamal signature on the integer message m under the public key. |
| `func EulerTotient(n *big.Int) *big.Int` | EulerTotient returns Euler's totient φ(n): the count of integers in [1, n] that are coprime to n. |
| `func EvaluatePolynomial(coeffs []*big.Int, x, m *big.Int) *big.Int` | EvaluatePolynomial evaluates the polynomial whose coefficients are given in ascending order (coeffs[0] is the constant term) at the point x modulo m,… |
| `func ExtendedGCD(a, b *big.Int) (g, x, y *big.Int)` | ExtendedGCD returns the greatest common divisor g of a and b together with Bézout coefficients x and y satisfying a*x + b*y = g. |
| `func FermatFactor(n *big.Int, maxSteps int) *big.Int` | FermatFactor returns a non-trivial factor of the odd composite n using Fermat's method, which searches for a representation n = a² − b² and is… |
| `func FermatTest(n *big.Int, rounds int, rng *rand.Rand) bool` | FermatTest runs the Fermat probable-prime test on n using the given number of random bases a drawn from rng, checking a^(n-1) ≡ 1 (mod n). |
| `func GCD(a, b *big.Int) *big.Int` | GCD returns the non-negative greatest common divisor of a and b. |
| `func GenerateSafePrime(bits int, rng *rand.Rand) *big.Int` | GenerateSafePrime returns a random safe prime p (one for which (p-1)/2 is also prime) with exactly the requested bit length, drawing randomness from… |
| `func IsCarmichaelNumber(n *big.Int) bool` | IsCarmichaelNumber reports whether n is a Carmichael number: a composite n that nevertheless satisfies a^(n−1) ≡ 1 (mod n) for every a coprime to… |
| `func IsPrime(n *big.Int) bool` | IsPrime reports whether n is prime. |
| `func IsPrimitiveRoot(g, p *big.Int) bool` | IsPrimitiveRoot reports whether g is a primitive root modulo the odd prime p, i.e. |
| `func IsProbablePrime(n *big.Int, rounds int, rng *rand.Rand) bool` | IsProbablePrime reports whether n is a probable prime, combining fast trial division by small primes with the deterministic Miller-Rabin witness set… |
| `func IsQuadraticResidue(a, p *big.Int) bool` | IsQuadraticResidue reports whether a is a quadratic residue modulo the odd prime p, i.e. |
| `func IsSafePrime(p *big.Int) bool` | IsSafePrime reports whether p is a safe prime, that is a prime for which (p-1)/2 is also prime. |
| `func IsSquareFree(n *big.Int) bool` | IsSquareFree reports whether n is divisible by no square of a prime. |
| `func Jacobi(a, n *big.Int) int` | Jacobi returns the Jacobi symbol (a/n) for an odd positive integer n. |
| `func LCM(a, b *big.Int) *big.Int` | LCM returns the non-negative least common multiple of a and b. |
| `func Legendre(a, p *big.Int) int` | Legendre returns the Legendre symbol (a/p) for an odd prime p: +1 if a is a non-zero quadratic residue mod p, -1 if it is a non-residue, and 0 if p… |
| `func MillerRabin(n *big.Int, rounds int, rng *rand.Rand) bool` | MillerRabin runs the probabilistic Miller-Rabin primality test on n using the given number of independently chosen random bases (rounds) drawn from… |
| `func MillerRabinDeterministic(n *big.Int) bool` | MillerRabinDeterministic performs a deterministic Miller-Rabin test using a fixed set of small prime bases. |
| `func Mobius(n *big.Int) int` | Mobius returns the Möbius function μ(n): 0 if n is divisible by the square of a prime, otherwise (−1)^k where k is the number of distinct prime… |
| `func ModAdd(a, b, m *big.Int) *big.Int` | ModAdd returns a+b mod m. |
| `func ModExp(base, exp, m *big.Int) *big.Int` | ModExp returns base**exp mod m using right-to-left binary (square-and-multiply) exponentiation. |
| `func ModInverse(a, m *big.Int) (*big.Int, error)` | ModInverse returns the multiplicative inverse of a modulo m, i.e. |
| `func ModMul(a, b, m *big.Int) *big.Int` | ModMul returns a*b mod m. |
| `func ModSqrt(a, p *big.Int) (*big.Int, error)` | ModSqrt returns a square root r of a modulo the odd prime p (so that r^2 ≡ a (mod p)), using the Tonelli-Shanks algorithm. |
| `func ModSub(a, b, m *big.Int) *big.Int` | ModSub returns a-b mod m, normalised into [0, m). |
| `func MultiplicativeOrder(a, n *big.Int) (*big.Int, error)` | MultiplicativeOrder returns the multiplicative order of a modulo n, that is the smallest positive integer k with a^k ≡ 1 (mod n). |
| `func NextPrime(n *big.Int) *big.Int` | NextPrime returns the smallest prime strictly greater than n. |
| `func NumberOfDivisors(n *big.Int) *big.Int` | NumberOfDivisors returns d(n), the number of positive divisors of n, computed as the product of (e+1) over the prime-power exponents e. |
| `func PohligHellman(g, h, p, order *big.Int) (*big.Int, error)` | PohligHellman solves the discrete logarithm g^x ≡ h (mod p) when the order of g is smooth, by solving the problem separately in each prime-power… |
| `func PollardPMinus1(n *big.Int, b int64) *big.Int` | PollardPMinus1 returns a non-trivial factor of n found by Pollard's p−1 method, which succeeds when some prime factor p of n has p−1 built… |
| `func PollardRho(n *big.Int) *big.Int` | PollardRho returns a non-trivial factor of the composite integer n using Pollard's rho algorithm with the iteration x -> x^2 + c. |
| `func PollardRhoBrent(n *big.Int) *big.Int` | PollardRhoBrent returns a non-trivial factor of the composite integer n using Brent's improvement to Pollard's rho, which reduces the number of gcd… |
| `func PrevPrime(n *big.Int) *big.Int` | PrevPrime returns the largest prime strictly less than n. |
| `func PrimeFactors(n *big.Int) []*big.Int` | PrimeFactors returns the distinct-with-multiplicity prime factors of n in ascending order (each prime appears as many times as it divides n). |
| `func PrimesUpTo(limit int64) []int64` | PrimesUpTo is an alias for SieveOfEratosthenes returning every prime not exceeding limit. |
| `func PrimitiveRoot(p *big.Int) *big.Int` | PrimitiveRoot returns the smallest primitive root modulo the odd prime p, i.e. |
| `func RSADecrypt(priv *RSAPrivateKey, c *big.Int) (*big.Int, error)` | RSADecrypt performs the raw RSA decryption primitive m = c^D mod N using the private exponent directly. |
| `func RSADecryptCRT(priv *RSAPrivateKey, c *big.Int) (*big.Int, error)` | RSADecryptCRT performs RSA decryption using the Chinese Remainder Theorem, exponentiating separately modulo P and Q and recombining. |
| `func RSAEncrypt(pub RSAPublicKey, m *big.Int) (*big.Int, error)` | RSAEncrypt performs the raw (textbook) RSA encryption primitive c = m^E mod N. The message m, interpreted as an integer, must satisfy 0 <= m < N;… |
| `func RSASign(priv *RSAPrivateKey, m *big.Int) (*big.Int, error)` | RSASign produces a raw RSA signature s = m^D mod N over the integer message m (0 <= m < N). |
| `func RSAVerify(pub RSAPublicKey, m, s *big.Int) bool` | RSAVerify checks a raw RSA signature by computing s^E mod N and reporting whether it equals the original message m. |
| `func RandomPrime(bits int, rng *rand.Rand) *big.Int` | RandomPrime returns a random prime with exactly the requested bit length (both the top and bottom bits are forced to 1, guaranteeing an odd number of… |
| `func SieveOfEratosthenes(limit int64) []int64` | SieveOfEratosthenes returns all primes p with 2 <= p <= limit as int64 values, computed with the classic sieve. |
| `func SmallestPrimeFactor(n *big.Int) *big.Int` | SmallestPrimeFactor returns the least prime factor of n. |
| `func SumOfDivisors(n *big.Int) *big.Int` | SumOfDivisors returns σ(n), the sum of the positive divisors of n, computed from the factorization as ∏ (p^(e+1) − 1)/(p − 1). |
| `func TrialDivision(n *big.Int) bool` | TrialDivision reports whether n is prime by trial division against every candidate up to sqrt(n). |

### Types

| Type | What it is |
| --- | --- |
| `DHParams` | DHParams holds the public Diffie-Hellman domain parameters shared by both parties: a prime modulus P and a generator G of a (large) subgroup of the… |
| `ElGamalPrivateKey` | ElGamalPrivateKey is an ElGamal private key: the secret exponent X together with the matching public key. |
| `ElGamalPublicKey` | ElGamalPublicKey is an ElGamal public key over the multiplicative group modulo the prime P with generator G, where Y = G^x mod P for the secret x. |
| `ElGamalSignature` | ElGamalSignature is a signature produced by ElGamalSign: the pair (R, S) with R = G^k mod P and S = k^-1 (m − X·R) mod (P−1) for a per-signature… |
| `Factor` | Factor represents a prime power in an integer's factorization: Prime raised to the given positive Exponent. |
| `RSAPrivateKey` | RSAPrivateKey is an RSA private key. |
| `RSAPublicKey` | RSAPublicKey is an RSA public key: the modulus N = p*q and the public exponent E. Encryption and signature verification use only these values. |
| `Share` | Share is a single point (X, Y) on the secret-sharing polynomial produced by SplitSecret. |

<details>
<summary><code>DHParams</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func GenerateDHParams(bits int, rng *rand.Rand) DHParams` | GenerateDHParams produces Diffie-Hellman domain parameters using a freshly generated safe prime P of the requested bit length and the generator G =… |
| `func (params DHParams) DHGeneratePrivate(rng *rand.Rand) *big.Int` | DHGeneratePrivate returns a random Diffie-Hellman private exponent in the range [2, P-2] drawn from rng. |
| `func (params DHParams) DHPublicKey(priv *big.Int) *big.Int` | DHPublicKey computes the public value G^priv mod P that a party transmits to its peer. |
| `func (params DHParams) DHSharedSecret(priv, peerPub *big.Int) (*big.Int, error)` | DHSharedSecret computes the shared secret peerPub^priv mod P from this party's private exponent and the peer's public value. |

</details>

<details>
<summary><code>ElGamalPrivateKey</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func GenerateElGamalKey(bits int, rng *rand.Rand) *ElGamalPrivateKey` | GenerateElGamalKey builds an ElGamal key pair over a fresh safe-prime group of the requested bit length, using generator G = 2 and a random secret… |

</details>

<details>
<summary><code>ElGamalSignature</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func ElGamalSign(priv *ElGamalPrivateKey, m *big.Int, rng *rand.Rand) (ElGamalSignature, error)` | ElGamalSign signs the integer message m under the private key using the ElGamal signature scheme, drawing the per-signature ephemeral exponent k from… |

</details>

<details>
<summary><code>Factor</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Factorization(n *big.Int) []Factor` | Factorization returns the prime factorization of n as a slice of Factor (prime, exponent) pairs sorted by ascending prime. |

</details>

<details>
<summary><code>RSAPrivateKey</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func GenerateRSAKey(bits int, e *big.Int, rng *rand.Rand) (*RSAPrivateKey, error)` | GenerateRSAKey generates an RSA private key whose modulus has approximately the requested bit length, using the public exponent e (65537 is standard… |
| `func (k *RSAPrivateKey) PublicKey() RSAPublicKey` | PublicKey returns the public half of the private key. |

</details>

<details>
<summary><code>Share</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func SplitSecret(secret *big.Int, threshold, participants int, prime *big.Int, rng *rand.Rand) ([]Share, error)` | SplitSecret splits secret into shares participants using Shamir's threshold scheme over the prime field Z/prime. |

</details>

### Variables

`ErrNoDiscreteLog`, `ErrNotInvertible`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/crypto).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
