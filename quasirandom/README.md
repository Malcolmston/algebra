# quasirandom — implements low-discrepancy (quasi-random) sequences and the numerical machinery

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/quasirandom.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/quasirandom)

Package quasirandom implements low-discrepancy (quasi-random) sequences and
the numerical machinery that surrounds them.

Quasi-random, or quasi-Monte-Carlo (QMC), point sets fill the unit cube
[0,1)^d far more evenly than independent pseudo-random samples. Their defining
quantity is the discrepancy, a measure of how far the empirical distribution
of the points departs from the uniform distribution. The Koksma–Hlawka
inequality bounds the error of an equal-weight cubature rule by the product of
the star discrepancy of the nodes and the variation of the integrand, so
sequences with asymptotically small discrepancy — of order (log N)^d / N
rather than the 1/sqrt(N) of plain Monte Carlo — give faster, deterministic
convergence.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/quasirandom"
```

## Usage

This is the package's own `ExampleDigits`, so it compiles and its output is
asserted on every `go test ./quasirandom/`.

```go
for _, tc := range []struct {
		n    uint64
		base int
	}{{1000, 10}, {255, 16}, {13, 2}} {
		d, err := quasirandom.Digits(tc.n, tc.base)
		if err != nil {
			panic(err)
		}
		back, err := quasirandom.FromDigits(d, tc.base)
		if err != nil {
			panic(err)
		}
		fmt.Printf("%d base %d -> %v -> %d\n", tc.n, tc.base, d, back)
	}
```

```
1000 base 10 -> [0 0 0 1] -> 1000
255 base 16 -> [15 15] -> 255
13 base 2 -> [1 0 1 1] -> 13
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func AdditiveRecurrenceSequence(seed, alpha float64, count int) ([]float64, error)` | AdditiveRecurrenceSequence returns count terms of the one-dimensional additive recurrence x_{n} = frac(seed + n*alpha), a shifted Weyl sequence. |
| `func AffinePermutation(n, a, b int) []int` | AffinePermutation returns the permutation i -> (a*i + b) mod n over the digit set {0,...,n-1}. |
| `func ApplyPermutation(xs []float64, p []int) []float64` | ApplyPermutation returns a new slice with entry i taken from xs[p`i`], reordering xs by the permutation p. |
| `func AreCoprimeBases(bases []int) bool` | AreCoprimeBases reports whether every pair of the supplied bases is coprime, the condition under which a multi-dimensional Halton construction is… |
| `func Bernoulli2(x float64) float64` | Bernoulli2 returns the second Bernoulli polynomial B2(x) = x^2 - x + 1/6, the building block of the diaphony's reproducing kernel. |
| `func BigBinomial(n, k int) (*big.Int, error)` | BigBinomial returns the binomial coefficient C(n,k) as an exact big.Int. |
| `func BigFactorial(n int) (*big.Int, error)` | BigFactorial returns n! as a big.Int, a small helper occasionally useful when working with the exact combinatorics behind Faure generator matrices. |
| `func BitLength(n uint64) int` | BitLength returns the number of bits required to represent n, i.e. |
| `func BoundingBox(points [][]float64) (lo, hi []float64, err error)` | BoundingBox returns the per-coordinate minimum and maximum of a point set. |
| `func BoxToUnit(points [][]float64, lower, upper []float64) ([][]float64, error)` | BoxToUnit maps a point set from the axis-aligned box [lower,upper] back into the unit cube, the inverse of ScaleToBox. |
| `func CenteredL2Discrepancy(points [][]float64) (float64, error)` | CenteredL2Discrepancy returns Hickernell's centered L2-discrepancy, which is invariant under reflection of the point set about the center of the cube… |
| `func CenteredL2DiscrepancySquared(points [][]float64) (float64, error)` | CenteredL2DiscrepancySquared returns the square of the centered L2-discrepancy. |
| `func Centroid(points [][]float64) ([]float64, error)` | Centroid returns the coordinate-wise mean of a point set. |
| `func ClampUnit(x float64) float64` | ClampUnit clamps x into the closed unit interval [0,1]. |
| `func ComposePermutations(p, q []int) []int` | ComposePermutations returns the permutation r with r`i` = p[q`i`], the effect of applying q and then p. |
| `func CoordinateMean(points [][]float64, k int) (float64, error)` | CoordinateMean returns the mean of coordinate k over the point set. |
| `func CoordinateMeans(points [][]float64) ([]float64, error)` | CoordinateMeans returns the vector of per-coordinate means of a point set. |
| `func CoordinateVariance(points [][]float64, k int) (float64, error)` | CoordinateVariance returns the population variance of coordinate k over the point set. |
| `func CoordinateVariances(points [][]float64) ([]float64, error)` | CoordinateVariances returns the vector of per-coordinate population variances of a point set. |
| `func CountInBox(points [][]float64, t []float64) (int, error)` | CountInBox returns the number of points lying in the half-open anchored box [0,t) whose upper corner is t. |
| `func CountPrimesBelow(n int) int` | CountPrimesBelow returns the number of primes strictly less than n (the prime-counting function pi evaluated just below n). |
| `func Diaphony(points [][]float64) (float64, error)` | Diaphony returns the classical (Zinterhof) diaphony of a point set, a Fourier-based uniformity measure with the closed form F_N^2 = (1/N^2) sum_{i,j}… |
| `func DigitAt(n uint64, base, pos int) (int, error)` | DigitAt returns the base-b digit of n at position pos (position zero is the least significant digit). |
| `func DigitCount(n uint64, base int) (int, error)` | DigitCount returns the number of base-b digits required to represent n, with zero requiring zero digits. |
| `func DigitScramble(base int, n uint64, perm []int) (float64, error)` | DigitScramble returns the value obtained by scrambling the base-b digits of the radical inverse of n with the single permutation perm applied to… |
| `func DigitSum(n uint64, base int) (int, error)` | DigitSum returns the sum of the base-b digits of n. |
| `func Digits(n uint64, base int) ([]int, error)` | Digits returns the base-b representation of the non-negative integer n as a slice of digits ordered from least significant to most significant. |
| `func DigitsPadded(n uint64, base, width int) ([]int, error)` | DigitsBig behaves like Digits but always returns a slice of length at least width, left-padding the most-significant end with zeros. |
| `func DirectionNumbers(dim int) ([]float64, error)` | DirectionNumbers returns the fractional direction numbers m_i / 2^i for the given Sobol dimension, one value per bit of resolution, as computed from… |
| `func EuclideanDistance(a, b []float64) (float64, error)` | EuclideanDistance returns the Euclidean distance between two equal-length points. |
| `func EvaluateSamples(f Integrand, points [][]float64) []float64` | EvaluateSamples returns the slice of function values f(p) for each point p, the raw samples underlying the QMC average and its error estimates. |
| `func Faure(dim int, n uint64) ([]float64, error)` | Faure returns the n-th point (zero-based) of the Faure sequence in the given dimension. |
| `func FaureBase(dim int) (int, error)` | FaureBase returns the base used by the Faure sequence in the given dimension: the smallest prime greater than or equal to the dimension (and at least… |
| `func FaureGeneratorMatrix(base, power, size int) ([][]int, error)` | FaureGeneratorMatrix returns the power-th power of the upper-triangular Pascal generator matrix modulo base, truncated to size rows and columns. |
| `func FaurePermutation(base int) []int` | FaurePermutation returns Faure's deterministic scrambling permutation of the digit set {0,...,base-1}. |
| `func FaureSequence(dim, count int) ([][]float64, error)` | FaureSequence returns the first count points of the dim-dimensional Faure sequence beginning at index zero. |
| `func FaureSequenceOffset(dim, start, count int) ([][]float64, error)` | FaureSequenceOffset returns count Faure points starting at index start. |
| `func FaureWithBase(dim, base int, n uint64) ([]float64, error)` | FaureWithBase returns the n-th Faure point using an explicit prime base, which must be at least as large as dim for the construction to be a genuine… |
| `func Fibonacci(n int) int` | Fibonacci returns the n-th Fibonacci number with Fibonacci(0)==0, Fibonacci(1)==1. |
| `func FibonacciLattice(m int) ([][]float64, error)` | FibonacciLattice returns the two-dimensional Fibonacci lattice with Fibonacci(m) points, the rank-1 lattice with generating vector (1,… |
| `func FillDistance(points [][]float64, res int) (float64, error)` | FillDistance returns the fill distance (covering radius) of a point set with respect to a regular test grid of res divisions per dimension: the… |
| `func FoldedRadicalInverse(base int, n uint64) (float64, error)` | FoldedRadicalInverse returns the folded (Hammersley–Kronecker) radical inverse, in which digit position k contributes (d_k + k) mod base rather… |
| `func Frac(x float64) float64` | Frac returns the fractional part of x in the half-open interval [0,1), computed as x - floor(x) so that negative inputs also map into [0,1). |
| `func FromDigits(digits []int, base int) (uint64, error)` | FromDigits reconstructs the integer whose base-b digits (least significant first) are given. |
| `func GCD(a, b int) int` | GCD returns the greatest common divisor of a and b, with GCD(0,0)==0. |
| `func GoldenRatioSequence(count int) ([]float64, error)` | GoldenRatioSequence returns count terms of the Weyl sequence generated by the golden ratio, the best-distributed one-dimensional additive sequence. |
| `func GrayCode(n uint64) uint64` | GrayCode returns the reflected binary Gray code of n, namely n XOR (n>>1). |
| `func Halton(dim int, n uint64) ([]float64, error)` | Halton returns the n-th point (zero-based) of the Halton sequence in the given dimension, using the first dim prime numbers as the coordinate bases. |
| `func HaltonBig(dim int, n *big.Int) ([]float64, error)` | HaltonBig returns the n-th Halton point for an index given as a big.Int, using the first dim primes as bases and returning each coordinate as a… |
| `func HaltonCoordinate(dim, k int, n uint64) (float64, error)` | HaltonCoordinate returns coordinate k (zero-based) of the n-th Halton point in the given dimension, the radical inverse of n in the (k+1)-th prime. |
| `func HaltonSequence(dim, count int) ([][]float64, error)` | HaltonSequence returns the first count points of the dim-dimensional Halton sequence starting from index zero as a slice of points. |
| `func HaltonSequenceOffset(dim, start, count int) ([][]float64, error)` | HaltonSequenceOffset returns count Halton points beginning at index start. |
| `func HaltonWithBases(bases []int, n uint64) ([]float64, error)` | HaltonWithBases returns the n-th Halton point using the supplied coordinate bases, one per dimension. |
| `func Hammersley(dim, i, total int) ([]float64, error)` | Hammersley returns the i-th point (zero-based, 0 <= i < total) of the dim-dimensional Hammersley point set of total points. |
| `func HammersleySet(dim, total int) ([][]float64, error)` | HammersleySet returns the complete dim-dimensional Hammersley point set of total points as a slice of points. |
| `func HammersleyWithBases(bases []int, i, total int) ([]float64, error)` | HammersleyWithBases returns the i-th Hammersley point whose trailing coordinates use the supplied bases (one per trailing dimension, so the point has… |
| `func IdentityPermutation(n int) []int` | IdentityPermutation returns the identity permutation (0,1,...,n-1). |
| `func InitialDirectionNumbers(dim int) ([]int, error)` | InitialDirectionNumbers returns the initial direction integers m_1..m_s for the given Sobol dimension. |
| `func IntegrateBox(f Integrand, lower, upper []float64, points [][]float64) (float64, error)` | IntegrateBox estimates the integral of f over the axis-aligned box [lower,upper] by mapping the unit-cube points into the box and multiplying the… |
| `func IntegrateUnitCube(f Integrand, points [][]float64) (float64, error)` | IntegrateUnitCube estimates the integral of f over the unit cube [0,1]^d as the equal-weight average of f over the supplied quasi-random points. |
| `func InverseGrayCode(g uint64) uint64` | InverseGrayCode returns the integer whose Gray code is g, inverting GrayCode. |
| `func InvertPermutation(p []int) []int` | InvertPermutation returns the inverse of the permutation p, so that InvertPermutation(p)[p`i`] == i. |
| `func IsPermutation(p []int) bool` | IsPermutation reports whether p is a permutation of {0,...,len(p)-1}: each value appears exactly once and lies in range. |
| `func IsPowerOfTwo(n uint64) bool` | IsPowerOfTwo reports whether n is a positive power of two. |
| `func IsPrime(n int) bool` | IsPrime reports whether n is a prime number. |
| `func KoksmaHlawkaBound(variation, starDiscrepancy float64) float64` | KoksmaHlawkaBound returns the Koksma–Hlawka upper bound on the absolute error of an equal-weight cubature rule, namely variation * discrepancy,… |
| `func KorobovGenerator(dim, a, total int) ([]int, error)` | KorobovGenerator returns the Korobov generating vector (1, a, a^2, ..., a^{dim-1}) reduced modulo total. |
| `func KorobovLattice(dim, a, total int) ([][]float64, error)` | KorobovLattice returns all total points of the Korobov lattice rule of the given dimension with parameter a. |
| `func KroneckerPoint(alpha []float64, n uint64) ([]float64, error)` | KroneckerPoint returns the n-th point (zero-based) of the multi-dimensional Kronecker sequence with generating vector alpha, the fractional part of… |
| `func KroneckerSequence(alpha []float64, count int) ([][]float64, error)` | KroneckerSequence returns the first count points of the Kronecker sequence with generating vector alpha. |
| `func L2StarDiscrepancy(points [][]float64) (float64, error)` | L2StarDiscrepancy returns the L2 star discrepancy of the point set using Warnock's closed-form formula, an O(d N^2) exact evaluation of the L2 norm… |
| `func L2StarDiscrepancySquared(points [][]float64) (float64, error)` | L2StarDiscrepancySquared returns the square of the L2 star discrepancy, avoiding the final square root when only the squared quantity is needed. |
| `func LeapedHalton(dim int, leap, n uint64) ([]float64, error)` | LeapedHalton returns the n-th point of a leaped Halton sequence, which takes every leap-th point of the ordinary sequence (index n*leap). |
| `func LocalDiscrepancy(points [][]float64, t []float64) (float64, error)` | LocalDiscrepancy returns the signed local discrepancy of a point set at the anchored box [0,t): the fraction of points inside minus the box volume. |
| `func LowestZeroBit(n uint64) int` | LowestZeroBit returns the position (from zero) of the lowest zero bit of n, which drives the Antonov–Saleev Gray-code Sobol recurrence. |
| `func MaxSobolDimension() int` | MaxSobolDimension returns the largest dimension for which built-in direction numbers are available (dimensions 1 through this value inclusive). |
| `func MinimumDistance(points [][]float64) (float64, error)` | MinimumDistance returns the smallest pairwise Euclidean distance in a point set (twice the separation distance). |
| `func ModifiedL2StarDiscrepancy(points [][]float64) (float64, error)` | ModifiedL2StarDiscrepancy returns Hickernell's modified L2-star discrepancy, which weights the anchored boxes so that low-dimensional projections… |
| `func MonteCarloStandardError(values []float64) float64` | MonteCarloStandardError returns the estimated standard error of the sample mean, sqrt(SampleVariance/n). |
| `func NextPowerOfTwo(n uint64) uint64` | NextPowerOfTwo returns the smallest power of two greater than or equal to n, with NextPowerOfTwo(0)==1. |
| `func NextPrime(n int) int` | NextPrime returns the smallest prime strictly greater than n. |
| `func NextPrimeGE(n int) int` | NextPrimeGE returns the smallest prime greater than or equal to n. |
| `func PlasticConstant(dim int) (float64, error)` | PlasticConstant returns the generalized golden ratio for the given dimension: the unique real root greater than one of x^(dim+1) = x + 1. |
| `func PopCount(n uint64) int` | PopCount returns the number of one bits in n. |
| `func PrevPrime(n int) int` | PrevPrime returns the largest prime strictly smaller than n, or zero when no such prime exists (n <= 2). |
| `func Prime(i int) (int, error)` | Prime returns the i-th prime number using one-based indexing, so Prime(1)==2, Prime(2)==3 and so on. |
| `func PrimeBases(dim int) ([]int, error)` | PrimeBases returns the first dim prime numbers, the canonical choice of pairwise-coprime bases for a Halton sequence of the given dimension. |
| `func PrimeIndex(p int) int` | PrimeIndex returns the one-based index of the prime p, or zero when p is not prime. |
| `func PrimeSieve(limit int) []int` | PrimeSieve returns all primes less than or equal to limit in increasing order using the sieve of Eratosthenes. |
| `func Primes(n int) ([]int, error)` | Primes returns the first n prime numbers (2, 3, 5, ...) in increasing order. |
| `func PrimitivePolynomialCoefficients(dim int) ([]int, error)` | PrimitivePolynomialCoefficients returns the full coefficient vector of the primitive polynomial for the given Sobol dimension, ordered from the… |
| `func PrimitivePolynomialDegree(dim int) (int, error)` | PrimitivePolynomialDegree returns the degree of the primitive polynomial used for the given Sobol dimension. |
| `func PrimitivePolynomialEncoding(dim int) (int, error)` | PrimitivePolynomialEncoding returns the integer that encodes the interior coefficients a_1..a_{s-1} of the primitive polynomial for the given Sobol… |
| `func Project(points [][]float64, dims []int) ([][]float64, error)` | Project returns the projection of a point set onto the coordinate axes listed in dims, producing a lower-dimensional point set. |
| `func QMCIntegrateFaure(f Integrand, dim, n int) (float64, error)` | QMCIntegrateFaure estimates the integral of f over [0,1]^dim using the first n points of the Faure sequence. |
| `func QMCIntegrateHalton(f Integrand, dim, n int) (float64, error)` | QMCIntegrateHalton estimates the integral of f over [0,1]^dim using the first n points of the Halton sequence. |
| `func QMCIntegrateHammersley(f Integrand, dim, n int) (float64, error)` | QMCIntegrateHammersley estimates the integral of f over [0,1]^dim using the n-point Hammersley set. |
| `func QMCIntegrateSobol(f Integrand, dim, n int) (float64, error)` | QMCIntegrateSobol estimates the integral of f over [0,1]^dim using the first n Sobol points after skipping the origin. |
| `func RadicalInverse(base int, n uint64) (float64, error)` | RadicalInverse returns the radical inverse phi_base(n): the value obtained by reflecting the base-b digits of n about the radix point, mapping the… |
| `func RadicalInverseBase2(n uint64) float64` | RadicalInverseBase2 returns the radical inverse in base two using pure bit reversal, faster than the general routine and exact for indices below 2^52. |
| `func RadicalInverseBig(base int, n *big.Int) (*big.Rat, error)` | RadicalInverseBig returns the radical inverse in the given base of a non-negative index supplied as a big.Int, as an exact big.Rat in [0,1). |
| `func RadicalInverseBigFloat(base int, n *big.Int) (float64, error)` | RadicalInverseBigFloat returns RadicalInverseBig evaluated to a float64, convenient when only ordinary precision is needed but the index exceeds… |
| `func RadicalInverseDigits(base int, n uint64) ([]int, error)` | RadicalInverseDigits returns the reflected base-b digit vector of n, ordered from the most significant original digit to the least; interpreting the… |
| `func RankOneLattice(gen []int, total int) ([][]float64, error)` | RankOneLattice returns all total points of the rank-1 lattice rule with the given generating vector. |
| `func RankOneLatticePoint(gen []int, total, i int) ([]float64, error)` | RankOneLatticePoint returns the i-th point (zero-based) of the rank-1 lattice rule with integer generating vector gen and total nodes, namely the… |
| `func ReflectPoint(p []float64) []float64` | ReflectPoint returns the reflection of a single point about the center of the unit cube, mapping each coordinate x to 1-x. |
| `func ReflectPoints(points [][]float64) [][]float64` | ReflectPoints returns the reflection of every point of a set about the center of the unit cube. |
| `func ReverseBits(n uint64, width int) uint64` | ReverseBits returns n with its low width bits reversed, the base-two analogue of ReverseDigits used implicitly by the base-two radical inverse. |
| `func ReverseDigits(n uint64, base, width int) (uint64, error)` | ReverseDigits returns the integer obtained by reversing the base-b digit string of n over exactly width digits. |
| `func ReverseHalton(dim int, n uint64) ([]float64, error)` | ReverseHalton returns the n-th point of the reverse (Faure–Tezuka) Halton sequence, using the deterministic per-base scrambling permutations… |
| `func ReversePermutation(n int) []int` | ReversePermutation returns the permutation (n-1,...,1,0) that reverses digit order, a simple scrambling of the digits {0,...,n-1}. |
| `func RobertsAlpha(dim int) ([]float64, error)` | RobertsAlpha returns the generating vector of Roberts' R_d low-discrepancy additive sequence in the given dimension: alpha_j = frac(g^-j) where g is… |
| `func RobertsPoint(dim int, n uint64) ([]float64, error)` | RobertsPoint returns the n-th point (zero-based) of Roberts' R_d additive low-discrepancy sequence in the given dimension, using the seed 1/2 so that… |
| `func RobertsSequence(dim, count int) ([][]float64, error)` | RobertsSequence returns the first count points of Roberts' R_d sequence in the given dimension. |
| `func SampleMean(values []float64) float64` | SampleMean returns the arithmetic mean of the values, or zero for an empty slice. |
| `func SampleVariance(values []float64) float64` | SampleVariance returns the unbiased (Bessel-corrected) sample variance of the values, or zero when fewer than two values are given. |
| `func ScaleToBox(points [][]float64, lower, upper []float64) ([][]float64, error)` | ScaleToBox maps a point set from the unit cube into the axis-aligned box [lower,upper]. |
| `func ScrambledHalton(perms [][]int, n uint64) ([]float64, error)` | ScrambledHalton returns the n-th Halton point in which each coordinate is a permutation-scrambled radical inverse. |
| `func ScrambledRadicalInverse(base int, n uint64, perm []int) (float64, error)` | ScrambledRadicalInverse returns the radical inverse of n in the given base with each digit d replaced by perm`d` before reflection. |
| `func SeparationDistance(points [][]float64) (float64, error)` | SeparationDistance returns half the minimum pairwise distance, the largest radius of non-overlapping balls centered at the points. |
| `func ShiftedRankOneLatticePoint(gen []int, shift []float64, total, i int) ([]float64, error)` | ShiftedRankOneLatticePoint returns the i-th point of the rank-1 lattice rule shifted by the vector shift (each coordinate taken modulo one), a cheap… |
| `func SobolPoint(dim int, n uint64) ([]float64, error)` | SobolPoint returns the n-th point (zero-based) of the dim-dimensional Sobol sequence directly. |
| `func SobolSequence(dim, count int) ([][]float64, error)` | SobolSequence returns the first count points of the dim-dimensional Sobol sequence beginning at the origin (index zero). |
| `func SobolSequenceSkip(dim, skip, count int) ([][]float64, error)` | SobolSequenceSkip returns count Sobol points after first skipping the leading skip points (a common way to discard the origin and low-index samples). |
| `func StarDiscrepancy(points [][]float64) (float64, error)` | StarDiscrepancy returns the exact star discrepancy D*_N of the point set, computed by the coordinate-grid method: the supremum of the local… |
| `func StarDiscrepancy1D(points [][]float64) (float64, error)` | StarDiscrepancy1D returns the exact star discrepancy of a one-dimensional point set (each point a length-one slice) using the closed form 1/(2N) +… |
| `func StarDiscrepancyGrid(points [][]float64, res int) (float64, error)` | StarDiscrepancyGrid returns an approximation to the star discrepancy obtained by sampling the local discrepancy on a regular grid of res divisions… |
| `func SymmetricL2Discrepancy(points [][]float64) (float64, error)` | SymmetricL2Discrepancy returns Hickernell's symmetric L2-discrepancy, which is invariant under reflection of any single coordinate. |
| `func ToroidalDistance(a, b []float64) (float64, error)` | ToroidalDistance returns the Euclidean distance between two points on the unit torus, where each coordinate difference is taken to be at most one… |
| `func TrailingZeroBits(n uint64) int` | TrailingZeroBits returns the number of trailing zero bits of n, i.e. |
| `func VanDerCorput(base int, n uint64) (float64, error)` | VanDerCorput returns the n-th term (zero-based) of the van der Corput sequence in the given base, which is exactly RadicalInverse(base, n). |
| `func VanDerCorputSequence(base, count int) ([]float64, error)` | VanDerCorputSequence returns the first count terms of the van der Corput sequence in the given base, starting from index zero (whose value is 0). |
| `func VanDerCorputSequenceOffset(base, start, count int) ([]float64, error)` | VanDerCorputSequenceOffset returns count terms of the van der Corput sequence in the given base beginning at index start rather than zero. |
| `func WeylSequence(alpha float64, n uint64) float64` | WeylSequence returns the n-th term (zero-based) of the one-dimensional Weyl (Kronecker) sequence with irrational parameter alpha, namely the… |
| `func WrapAroundL2Discrepancy(points [][]float64) (float64, error)` | WrapAroundL2Discrepancy returns Hickernell's wrap-around L2-discrepancy, which treats the unit cube as a torus and so is invariant under coordinate… |

### Types

| Type | What it is |
| --- | --- |
| `Integrand` | Integrand is a scalar function of a point in R^d used by the integration helpers. |
| `Sobol` | Sobol is a stateful generator of the multi-dimensional Sobol sequence. |

<details>
<summary><code>Sobol</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewSobol(dim int) (*Sobol, error)` | NewSobol returns a Sobol generator for the given dimension using the built-in direction numbers. |
| `func (s *Sobol) Dimension() int` | Dimension returns the number of coordinates the generator produces. |
| `func (s *Sobol) Index() uint64` | Index returns the zero-based index of the point that the next call to Next will return. |
| `func (s *Sobol) Next() []float64` | Next returns the current Sobol point and advances the generator by one index using the Antonov–Saleev Gray-code recurrence. |
| `func (s *Sobol) Point(n uint64) []float64` | Point returns the n-th Sobol point (zero-based) directly, without disturbing the generator's streaming state, by expanding the Gray code of n. |
| `func (s *Sobol) Reset()` | Reset returns the generator to its initial state, so the next call to Next again yields the origin. |
| `func (s *Sobol) Skip(k uint64)` | Skip advances the generator by k indices without returning the intervening points, leaving it positioned so that the next Next returns point… |

</details>

### Constants

`GoldenRatio`, `GoldenRatioConjugate`

### Variables

`ErrBadBase`, `ErrDimension`, `ErrEmptyPointSet`, `ErrNonPositive`, `ErrOutOfUnitCube`, `ErrRaggedPointSet`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/quasirandom).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
