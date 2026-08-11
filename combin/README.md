# combin — implements combinatorics: counting, arrangements and combinatorial numbers

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/combin.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/combin)

Package combin implements combinatorics: counting, arrangements and
combinatorial numbers.

The package collects the classical counting functions used throughout discrete
mathematics. Where results grow without bound the routines return
arbitrary-precision *big.Int values so that answers stay exact; parallel
float64 and natural-logarithm variants are provided for the cases where an
approximate but fast answer (or one that survives magnitudes far beyond the
range of a machine integer) is preferable.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/combin"
```

## Usage

This is the package's own `ExampleBernoulliNumber`, so it compiles and its output is
asserted on every `go test ./combin/`.

```go
for n := 0; n <= 6; n++ {
		fmt.Print(combin.BernoulliNumber(n).RatString(), " ")
	}
	fmt.Println()
```

```
1 -1/2 1/6 0 -1/30 0 1/42
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func AllPermutations(items []int) [][]int` | AllPermutations returns all len(items)! permutations of items (treating equal values as distinct by position), each as a fresh slice. |
| `func Bell(n int) *big.Int` | Bell returns the nth Bell number, the number of partitions of a set of n elements. |
| `func BellTriangleFull(n int) [][]*big.Int` | BellTriangleFull returns rows 0 through n of the Bell (Aitken/Peirce) triangle. |
| `func BellTriangleRow(n int) []*big.Int` | BellTriangleRow returns row n of the Bell (Aitken) triangle. |
| `func BernoulliNumber(n int) *big.Rat` | BernoulliNumber returns the n-th Bernoulli number B(n) as an exact rational, using the convention B(1) = -1/2. |
| `func BernoulliNumberFloat(n int) float64` | BernoulliNumberFloat returns the n-th Bernoulli number as a float64. |
| `func BernoulliSequence(n int) []*big.Rat` | BernoulliSequence returns the Bernoulli numbers B(0), B(1), ..., B(n) as exact rationals, using the convention B(1) = -1/2. |
| `func Binomial(n, k int) *big.Int` | Binomial returns the binomial coefficient C(n, k) = n!/(k!·(n-k)!) using a multiplicative scheme that keeps every intermediate value integral. |
| `func BinomialFloat(n, k int) float64` | BinomialFloat returns the binomial coefficient C(n, k) as a float64. |
| `func BinomialRow(n int) []*big.Int` | BinomialRow is a synonym for PascalRow: it returns the binomial coefficients C(n, 0), C(n, 1), …, C(n, n). |
| `func Catalan(n int) *big.Int` | Catalan returns the nth Catalan number, C(2n, n)/(n+1). |
| `func CatalanFloat(n int) float64` | CatalanFloat returns the nth Catalan number as a float64. |
| `func CatalanSequence(n int) []*big.Int` | CatalanSequence returns the Catalan numbers C(0), C(1), ..., C(n), using C(i+1) = C(i) (2)(2i+1) / (i+2). |
| `func CatalanTriangle(n, k int) *big.Int` | CatalanTriangle returns the entry C(n, k) of Catalan's triangle, the number of lattice paths that never rise above the diagonal, equal to binom(n+k,… |
| `func CentralBinomial(n int) *big.Int` | CentralBinomial returns the central binomial coefficient C(2n, n). |
| `func CombinationList(n, k int) [][]int` | CombinationList returns all k-element index subsets of {0, 1, ..., n-1}, in lexicographic order, each as a strictly increasing slice. |
| `func CombinationListOf(items []int, k int) [][]int` | CombinationListOf returns all k-element combinations of items, in the order induced by index position, each as a fresh slice. |
| `func CombinationRank(c []int) *big.Int` | CombinationRank returns the colexicographic rank of the strictly increasing combination c in the combinatorial number system: sum of C(c`i`, i+1). |
| `func CombinationUnrank(k int, rank *big.Int) []int` | CombinationUnrank returns the k-combination (strictly increasing) whose colexicographic rank is rank, inverting CombinationRank. |
| `func Combinations(n, r int) *big.Int` | Combinations returns the number of unordered selections of k items from n, nCr = C(n, k). |
| `func CombinationsFloat(n, r int) float64` | CombinationsFloat returns C(n, k) as a float64, evaluated via logarithms so that it stays finite for large arguments. |
| `func CombinationsWithRepetition(n, r int) *big.Int` | CombinationsWithRepetition returns the number of multisets of size r drawn from n distinct items, C(n+r-1, r). |
| `func CompositionList(n int) [][]int` | CompositionList returns all compositions of n into positive parts, in lexicographic order. |
| `func CompositionListInto(n, k int) [][]int` | CompositionListInto returns all compositions of n into exactly k positive parts. |
| `func CompositionNumber(n int) int` | CompositionNumber returns the number of compositions (ordered partitions) of n into positive parts, which is 2^(n-1) for n >= 1 and 1 for n == 0. |
| `func CompositionNumberInto(n, k int) int` | CompositionNumberInto returns the number of compositions of n into exactly k positive parts, which is C(n-1, k-1). |
| `func Compositions(n int) *big.Int` | Compositions returns the number of ordered compositions of n into positive parts, 2^(n-1) for n >= 1 and 1 for n = 0. |
| `func DelannoyNumber(m, n int) *big.Int` | DelannoyNumber returns the Delannoy number D(m, n), the number of lattice paths from (0,0) to (m,n) using steps east, north and north-east. |
| `func Derangement(n int) *big.Int` | Derangement is a synonym for Subfactorial: it returns the number of permutations of n elements that leave no element fixed. |
| `func DerangementList(n int) [][]int` | DerangementList returns every derangement of {0, 1, ..., n-1} as a slice, in lexicographic order. |
| `func DerangementNumber(n int) int` | DerangementNumber returns !n, the number of permutations of n elements with no fixed point, using the recurrence !n = (n-1)(!(n-1) + !(n-2)). |
| `func DerangementProbability(n int) float64` | DerangementProbability returns !n / n!, the probability that a uniformly random permutation of n elements is a derangement. |
| `func DoubleFactorial(n int) *big.Int` | DoubleFactorial returns the double factorial n!! = n·(n-2)·(n-4)·…, stopping at 2 or 1. |
| `func DoubleFactorialFloat(n int) float64` | DoubleFactorialFloat returns the double factorial n!! as a float64. |
| `func EulerNumber(n int) *big.Int` | EulerNumber returns the n-th Euler (secant) number E(n). |
| `func EulerNumberSequence(n int) []*big.Int` | EulerNumberSequence returns the Euler numbers E(0), E(1), ..., E(n). |
| `func EulerianNumber(n, k int) *big.Int` | EulerianNumber returns the Eulerian number A(n, k), the number of permutations of n elements with exactly k ascents. |
| `func EulerianRow(n int) []*big.Int` | EulerianRow returns row n of the Eulerian triangle: the counts of permutations of n elements with 0, 1, ..., n-1 ascents. |
| `func EulerianSecondOrder(n, k int) *big.Int` | EulerianSecondOrder returns the second-order Eulerian number <<n, k>>, satisfying <<n, k>> = (k+1)·<<n-1, k>> + (2n-1-k)·<<n-1, k-1>> with <<0, 0>>… |
| `func EulerianTriangle(n int) [][]*big.Int` | EulerianTriangle returns rows 0 through n of the Eulerian triangle. |
| `func Factorial(n int) *big.Int` | Factorial returns n! = 1·2·…·n as an exact integer. |
| `func FactorialFloat(n int) float64` | FactorialFloat returns n! as a float64, using the Gamma function so that the result stays finite up to n == 170. |
| `func FactorialNumberSystem(rank *big.Int, n int) []int` | FactorialNumberSystem returns the factoradic (Lehmer) digits of rank as a length-n slice; digit i lies in [0, n-1-i] and the final digit is always 0. |
| `func FallingFactorial(x, n int) *big.Int` | FallingFactorial returns the falling factorial (x)_n = x·(x-1)·…·(x-n+1), with the empty product FallingFactorial(x, 0) = 1. |
| `func FallingFactorialFloat(x float64, n int) float64` | FallingFactorialFloat returns the falling factorial (x)_n for a real base x. |
| `func Fibonacci(n int) *big.Int` | Fibonacci returns the nth Fibonacci number, F(0) = 0, F(1) = 1 and F(n) = F(n-1) + F(n-2). |
| `func GammaBinomial(x, y float64) float64` | GammaBinomial returns the generalized binomial coefficient for real arguments, C(x, y) = Γ(x+1) / (Γ(y+1)·Γ(x-y+1)). |
| `func GeneralizedPentagonal(i int) int` | GeneralizedPentagonal returns the i-th generalized pentagonal number in the order 0, 1, 2, 5, 7, 12, 15, ... |
| `func GrayCodeAt(rank uint64) uint64` | GrayCodeAt returns the Gray code at position rank in the standard reflected Gray sequence; it is an alias for GrayEncode. |
| `func GrayCodeRank(g uint64) uint64` | GrayCodeRank returns the position (rank) of the Gray code g within the standard reflected Gray sequence; it is an alias for GrayDecode. |
| `func GrayCodeSequence(bits int) []uint64` | GrayCodeSequence returns the 2^bits reflected Gray codes in order. |
| `func GrayDecode(g uint64) uint64` | GrayDecode returns the integer whose Gray code is g, inverting GrayEncode. |
| `func GrayEncode(n uint64) uint64` | GrayEncode returns the reflected binary Gray code of n (n XOR (n >> 1)). |
| `func HarmonicNumber(n int) float64` | HarmonicNumber returns the nth harmonic number H(n) = 1 + 1/2 + … + 1/n as a float64, with H(0) = 0. |
| `func HarmonicNumberGeneralized(n, m int) float64` | HarmonicNumberGeneralized returns the generalized harmonic number H(n, m) = 1 + 1/2^m + … + 1/n^m as a float64, with H(0, m) = 0. |
| `func HyperFactorial(n int) *big.Int` | HyperFactorial returns the hyperfactorial H(n) = 1^1·2^2·…·n^n, with HyperFactorial(0) = 1. |
| `func InvolutionList(n int) [][]int` | InvolutionList returns every involution of {0, 1, ..., n-1} as a slice, in lexicographic order. |
| `func InvolutionNumber(n int) int` | InvolutionNumber returns the number of involutions (self-inverse permutations) of n elements, using I(n) = I(n-1) + (n-1) I(n-2). |
| `func LahNumber(n, k int) *big.Int` | LahNumber returns the unsigned Lah number L(n, k) = C(n-1, k-1)·n!/k!, the number of ways to partition n elements into k non-empty ordered lists. |
| `func LogBinomial(n, k int) float64` | LogBinomial returns the natural logarithm of the binomial coefficient C(n, k). |
| `func LogCombinations(n, r int) float64` | LogCombinations returns the natural logarithm of C(n, k). |
| `func LogFactorial(n int) float64` | LogFactorial returns the natural logarithm of n!, computed from the log-Gamma function so that it is accurate for very large n. |
| `func LogMultinomial(ks ...int) float64` | LogMultinomial returns the natural logarithm of the multinomial coefficient for the given part sizes. |
| `func LogPermutations(n, r int) float64` | LogPermutations returns the natural logarithm of nPr. |
| `func Lucas(n int) *big.Int` | Lucas returns the nth Lucas number, L(0) = 2, L(1) = 1 and L(n) = L(n-1) + L(n-2). |
| `func LyndonWordNumber(n, k int) *big.Int` | LyndonWordNumber returns the number of Lyndon words of length n over a k-letter alphabet: (1/n) sum_{d\|n} mu(d) k^(n/d). |
| `func MotzkinNumber(n int) *big.Int` | MotzkinNumber returns the nth Motzkin number, counting the ways of drawing non-intersecting chords between n points on a circle. |
| `func MotzkinSequence(n int) []*big.Int` | MotzkinSequence returns the Motzkin numbers M(0), M(1), ..., M(n), using M(i) = M(i-1) + sum_{j=0}^{i-2} M(j) M(i-2-j). |
| `func Multifactorial(n, k int) *big.Int` | Multifactorial returns the k-fold factorial n·(n-k)·(n-2k)·… down to the smallest positive term. |
| `func Multinomial(ks ...int) *big.Int` | Multinomial returns the multinomial coefficient (k1+k2+…)! / (k1!·k2!·…), the number of ways to partition a set of that many items into… |
| `func MultinomialFloat(ks ...int) float64` | MultinomialFloat returns the multinomial coefficient as a float64, evaluated through logarithms. |
| `func MultisetCoefficient(n, k int) *big.Int` | MultisetCoefficient returns the number of multisets of size k drawn from n distinct kinds, i.e. |
| `func MultisetPermutationList(counts []int) [][]int` | MultisetPermutationList expands the count vector counts into a multiset (value i repeated counts`i` times) and returns all its distinct permutations. |
| `func MultisetPermutationNumber(counts []int) *big.Int` | MultisetPermutationNumber returns the number of distinct permutations of a multiset whose element multiplicities are counts, i.e. |
| `func MultisetPermutations(items []int) [][]int` | MultisetPermutations returns all distinct permutations of items (a multiset), in lexicographic order, each as a fresh slice. |
| `func NarayanaNumber(n, k int) *big.Int` | NarayanaNumber returns the Narayana number N(n, k) = C(n, k)·C(n, k-1)/n, counting Dyck paths of length 2n with exactly k peaks. |
| `func NarayanaTriangle(n int) [][]*big.Int` | NarayanaTriangle returns rows 0 through n of the Narayana triangle. |
| `func NarayanaTriangleRow(n int) []*big.Int` | NarayanaTriangleRow returns row n of the Narayana triangle, the values N(n,1), ..., N(n,n) where N(n,k) = C(n,k) C(n,k-1) / n. |
| `func NecklaceNumber(n, k int) *big.Int` | NecklaceNumber returns the number of distinct necklaces of length n using up to k colors, counted up to rotation: (1/n) sum_{d\|n} phi(d) k^(n/d). |
| `func NextCombination(c []int, n int) bool` | NextCombination advances c, a strictly increasing k-combination of indices drawn from {0, ..., n-1}, to the next combination in lexicographic order. |
| `func NextPartition(p []int) ([]int, bool)` | NextPartition advances p (a non-increasing partition) in place to the next partition of the same integer in reverse lexicographic order. |
| `func NextPermutation(a []int) bool` | NextPermutation rearranges a into the lexicographically next greater permutation of its elements and returns true. |
| `func PartitionConjugate(p []int) []int` | PartitionConjugate returns the conjugate (transpose) of the partition p, which must be given in non-increasing order. |
| `func PartitionCount(n int) *big.Int` | PartitionCount returns p(n), the number of unordered ways to write n as a sum of positive integers, computed with Euler's pentagonal-number… |
| `func PartitionDurfeeSquare(p []int) int` | PartitionDurfeeSquare returns the side length of the Durfee square of the partition p (the largest d with at least d parts each of size >= d). |
| `func PartitionNumber(n int) int` | PartitionNumber returns p(n), the number of integer partitions of n, computed by an O(n^2) dynamic program. |
| `func PartitionNumberAtMost(n, k int) int` | PartitionNumberAtMost returns the number of partitions of n into at most k parts, equivalently the number of partitions of n whose parts are all <= k. |
| `func PartitionNumberDistinct(n int) int` | PartitionNumberDistinct returns the number of partitions of n into distinct (unequal) positive parts. |
| `func PartitionNumberInto(n, k int) int` | PartitionNumberInto returns the number of partitions of n into exactly k positive parts, using the recurrence P(n,k) = P(n-1,k-1) + P(n-k,k). |
| `func PartitionNumberOdd(n int) int` | PartitionNumberOdd returns the number of partitions of n into odd parts. |
| `func PartitionNumberTable(n int) []int` | PartitionNumberTable returns the slice [p(0), p(1), ..., p(n)]. |
| `func PartitionRank(p []int) int` | PartitionRank returns the rank of the partition p (largest part minus the number of parts). |
| `func Partitions(n int) [][]int` | Partitions returns all integer partitions of n as non-increasing slices, in reverse lexicographic order (starting with {n} and ending with all ones). |
| `func PartitionsDistinct(n int) [][]int` | PartitionsDistinct returns all partitions of n into distinct positive parts, each as a strictly decreasing slice. |
| `func PartitionsInto(n, k int) [][]int` | PartitionsInto returns all partitions of n into exactly k positive parts, each as a non-increasing slice. |
| `func PartitionsIntoKParts(n, k int) *big.Int` | PartitionsIntoKParts returns the number of partitions of n into exactly k positive parts, using the recurrence p(n, k) = p(n-1, k-1) + p(n-k, k) with… |
| `func PascalRow(n int) []*big.Int` | PascalRow returns row n of Pascal's triangle, the binomial coefficients C(n, 0), C(n, 1), …, C(n, n). |
| `func PascalTriangle(n int) [][]*big.Int` | PascalTriangle returns rows 0 through n of Pascal's triangle. |
| `func PentagonalNumber(k int) int` | PentagonalNumber returns the k-th (signed) pentagonal number k(3k-1)/2. |
| `func PermutationParity(a []int) int` | PermutationParity returns the sign of the permutation a of {0,...,n-1}: +1 if it has an even number of inversions, -1 if odd. |
| `func PermutationRank(a []int) *big.Int` | PermutationRank returns the lexicographic rank (0-based) of the permutation a, which must be a permutation of {0, 1, ..., len(a)-1}. |
| `func PermutationUnrank(n int, rank *big.Int) []int` | PermutationUnrank returns the permutation of {0, 1, ..., n-1} whose lexicographic rank (0-based) is rank. |
| `func Permutations(n, r int) *big.Int` | Permutations returns the number of ordered arrangements of r items chosen from n, nPr = n!/(n-r)!. |
| `func PermutationsFloat(n, r int) float64` | PermutationsFloat returns nPr as a float64, evaluated through logarithms so that it remains finite for large arguments. |
| `func PermutationsWithRepetition(n, r int) *big.Int` | PermutationsWithRepetition returns n^r, the number of ordered sequences of length r drawn from n symbols when repetition is allowed. |
| `func PowerSet(items []int) [][]int` | PowerSet returns every subset of items as a fresh slice, ordered by the binary value of the inclusion mask (the empty set first). |
| `func PowerSetSize(n int) *big.Int` | PowerSetSize returns 2^n, the number of subsets of an n-element set. |
| `func PrevPermutation(a []int) bool` | PrevPermutation rearranges a into the lexicographically previous permutation of its elements and returns true. |
| `func Primorial(n int) *big.Int` | Primorial returns n# — the product of all prime numbers not exceeding n. |
| `func RencontresNumber(n, k int) *big.Int` | RencontresNumber returns the number of permutations of n elements with exactly k fixed points, C(n, k)·!(n-k). |
| `func RisingFactorial(x, n int) *big.Int` | RisingFactorial returns the rising factorial (Pochhammer symbol) x^(n) = x·(x+1)·…·(x+n-1), with the empty product RisingFactorial(x, 0) = 1. |
| `func RisingFactorialFloat(x float64, n int) float64` | RisingFactorialFloat returns the rising factorial x^(n) for a real base x. |
| `func SchroederNumber(n int) *big.Int` | SchroederNumber returns the nth large Schröder number, counting the lattice paths from (0,0) to (n,n) using steps east, north and north-east that… |
| `func SecantNumber(n int) *big.Int` | SecantNumber returns the n-th secant number \|E(2n)\| = A(2n), the number of alternating permutations of even length 2n. |
| `func SignedLahNumber(n, k int) *big.Int` | SignedLahNumber returns the signed Lah number (-1)^n·L(n, k), the coefficients that express rising factorials in terms of falling factorials. |
| `func StirlingApproximation(n int) float64` | StirlingApproximation returns Stirling's asymptotic estimate of n!, sqrt(2·π·n)·(n/e)^n, as a float64. |
| `func StirlingFirst(n, k int) *big.Int` | StirlingFirst returns the signed Stirling number of the first kind s(n, k) = (-1)^(n-k)·c(n, k), the coefficients relating falling factorials to… |
| `func StirlingFirstRow(n int) []*big.Int` | StirlingFirstRow returns the unsigned Stirling numbers of the first kind c(n, 0), c(n, 1), …, c(n, n). |
| `func StirlingFirstUnsigned(n, k int) *big.Int` | StirlingFirstUnsigned returns the unsigned Stirling number of the first kind c(n, k), the number of permutations of n elements with exactly k… |
| `func StirlingSecond(n, k int) *big.Int` | StirlingSecond returns the Stirling number of the second kind S(n, k), the number of ways to partition n elements into k non-empty unlabelled subsets. |
| `func StirlingSecondRow(n int) []*big.Int` | StirlingSecondRow returns the Stirling numbers of the second kind S(n, 0), S(n, 1), …, S(n, n). |
| `func Subfactorial(n int) *big.Int` | Subfactorial returns the number of derangements !n, i.e. |
| `func SubfactorialFloat(n int) float64` | SubfactorialFloat returns the number of derangements !n as a float64, using the closed form round(n!/e). |
| `func SubsetList(n int) [][]int` | SubsetList returns every subset of {0, 1, ..., n-1} as an index slice, ordered by the binary value of the inclusion mask. |
| `func Superfactorial(n int) *big.Int` | Superfactorial returns the product of the first n factorials, 1!·2!·…·n!, with Superfactorial(0) = 1. |
| `func TangentNumber(n int) *big.Int` | TangentNumber returns the n-th tangent number A(2n+1), the number of alternating permutations of odd length 2n+1. |
| `func TelephoneNumber(n int) *big.Int` | TelephoneNumber returns the nth telephone (involution) number, the number of self-inverse permutations of n elements. |
| `func Tribonacci(n int) *big.Int` | Tribonacci returns the nth tribonacci number with T(0) = T(1) = 0, T(2) = 1 and T(n) = T(n-1) + T(n-2) + T(n-3). |
| `func WeakCompositionList(n, k int) [][]int` | WeakCompositionList returns all weak compositions of n into exactly k non-negative parts, in lexicographic order. |
| `func WeakCompositionNumber(n, k int) int` | WeakCompositionNumber returns the number of weak compositions of n into exactly k non-negative parts, which is C(n+k-1, k-1). |
| `func WeakCompositions(n, k int) *big.Int` | WeakCompositions returns the number of ways to write n as an ordered sum of k non-negative integers, C(n+k-1, k-1). |
| `func ZigzagNumber(n int) *big.Int` | ZigzagNumber returns the n-th Euler up/down (zigzag) number A(n): the number of alternating permutations of n elements. |

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/combin).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
