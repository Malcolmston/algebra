# discretemath — broad collection of discrete-mathematics primitives built entirely on the Go standard

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/discretemath.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/discretemath)

Package discretemath implements a broad collection of discrete-mathematics
primitives built entirely on the Go standard library.

The package is deliberately self-contained and deterministic: identical inputs
always yield identical outputs, no global state is mutated after
initialization, and no third-party dependency is required.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/discretemath"
```

## Usage

This is the package's own `ExampleCartesianProductN`, so it compiles and its output is
asserted on every `go test ./discretemath/`.

```go
for _, t := range discretemath.CartesianProductN([]int{0, 1}, []int{7, 8}) {
		fmt.Printf("%v ", t)
	}
	fmt.Println()
```

```
[0 7] [0 8] [1 7] [1 8]
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func BitString(x uint64, width int) string` | BitString returns the base-2 representation of x, left-padded with zeros to at least width characters. |
| `func CartesianProductN[T any](sets ...[]T) [][]T` | CartesianProductN returns the Cartesian product of any number of homogeneous slices as a slice of tuples, each tuple having one element from each… |
| `func ClearBit(x uint64, i uint) uint64` | ClearBit returns x with bit i set to zero. |
| `func ClearLowestSetBit(x uint64) uint64` | ClearLowestSetBit returns x with its lowest set bit cleared. |
| `func DeBruijnSequence(k, n int) ([]int, error)` | DeBruijnSequence returns a De Bruijn sequence B(k, n): a sequence of symbols drawn from the alphabet {0, 1, ..., k-1} of length k**n in which every… |
| `func DeBruijnString(alphabet string, n int) (string, error)` | DeBruijnString returns a De Bruijn sequence of order n over the given alphabet as a string. |
| `func FromBase(s string, base int) (int64, error)` | FromBase parses the string s as a signed integer in the given base, which must be in the range 2..36. |
| `func FromBaseUint(s string, base int) (uint64, error)` | FromBaseUint parses the string s as an unsigned integer in the given base, which must be in the range 2..36. |
| `func GrayDecode(g uint64) uint64` | GrayDecode converts a reflected binary Gray code g back into the ordinary binary value it represents. |
| `func GrayEncode(x uint64) uint64` | GrayEncode converts the binary value x into its reflected binary Gray code, in which successive integers differ in exactly one bit. |
| `func HammingDistance(a, b uint64) int` | HammingDistance returns the number of bit positions at which a and b differ, which is the population count of their exclusive-or. |
| `func HammingDistanceBytes(a, b []byte) (int, error)` | HammingDistanceBytes returns the number of differing bits between two equal length byte slices. |
| `func HammingDistanceString(a, b string) (int, error)` | HammingDistanceString returns the number of positions at which two equal length strings differ, comparing by rune. |
| `func IntToRoman(n int) (string, error)` | IntToRoman converts an integer in the range 1..3999 into its canonical Roman numeral using subtractive notation (for example 1994 becomes "MCMXCIV"). |
| `func IntToWords(n int64) string` | IntToWords converts a signed 64-bit integer into its English cardinal words, covering the full int64 range up to quintillions. |
| `func IntToWordsOrdinal(n int64) string` | IntToWordsOrdinal converts a signed 64-bit integer into its English ordinal words (for example 21 becomes "twenty-first" and 100 becomes "one… |
| `func IsPowerOfTwo(x uint64) bool` | IsPowerOfTwo reports whether x is a positive power of two (1, 2, 4, 8, ...). |
| `func IsValidRoman(s string) bool` | IsValidRoman reports whether s is a well-formed canonical Roman numeral in the range 1..3999. |
| `func LeadingZeros(x uint64) int` | LeadingZeros returns the number of consecutive zero bits at the high end of x. |
| `func LevenshteinDistance(a, b string) int` | LevenshteinDistance returns the Levenshtein edit distance between two strings: the minimum number of single-character insertions, deletions and… |
| `func Log2Ceil(x uint64) int` | Log2Ceil returns the ceiling of the base-2 logarithm of x, equivalently the smallest e such that 2**e >= x. |
| `func Log2Floor(x uint64) int` | Log2Floor returns the floor of the base-2 logarithm of x, equivalently the index of its highest set bit. |
| `func LongestCommonSubsequence(a, b string) string` | LongestCommonSubsequence returns one longest common subsequence of two strings as a string, comparing by rune. |
| `func LongestCommonSubsequenceLength(a, b string) int` | LongestCommonSubsequenceLength returns the length of the longest common subsequence of two strings, comparing by rune. |
| `func LowestSetBit(x uint64) uint64` | LowestSetBit returns x with all bits cleared except its lowest set bit. |
| `func MortonDecode2D(m uint64) (x, y uint32)` | MortonDecode2D recovers the two coordinates x and y from a 2-D Morton code m. |
| `func MortonDecode3D(m uint64) (x, y, z uint32)` | MortonDecode3D recovers the three coordinates x, y and z from a 3-D Morton code m. |
| `func MortonEncode2D(x, y uint32) uint64` | MortonEncode2D interleaves the bits of x and y into a single Morton (Z-order) code, with the bits of x occupying the even positions and the bits of y… |
| `func MortonEncode3D(x, y, z uint32) uint64` | MortonEncode3D interleaves the low 21 bits of x, y and z into a single 3-D Morton (Z-order) code. |
| `func NextPowerOfTwo(x uint64) uint64` | NextPowerOfTwo returns the smallest power of two that is greater than or equal to x. |
| `func OptimalStringAlignmentDistance(a, b string) int` | OptimalStringAlignmentDistance returns the optimal string alignment distance between two strings, the restricted form of the Damerau-Levenshtein… |
| `func Parity(x uint64) int` | Parity returns the parity of x: 1 if the number of set bits is odd and 0 if it is even. |
| `func PopCount(x uint64) int` | PopCount returns the number of one bits (the population count) in x. |
| `func PowerSet[T any](elements []T) [][]T` | PowerSet returns every subset of the given elements as a slice of slices. |
| `func PrevPowerOfTwo(x uint64) uint64` | PrevPowerOfTwo returns the largest power of two that is less than or equal to x. |
| `func ReverseBits(x uint64) uint64` | ReverseBits returns x with the order of its 64 bits reversed, so that bit i of the result equals bit 63-i of x. |
| `func ReverseBits32(x uint32) uint32` | ReverseBits32 returns x with the order of its 32 bits reversed. |
| `func RomanToInt(s string) (int, error)` | RomanToInt converts a Roman numeral into its integer value. |
| `func RotateLeft(x uint64, k int) uint64` | RotateLeft returns x rotated left by k bits within a 64-bit word. |
| `func RotateRight(x uint64, k int) uint64` | RotateRight returns x rotated right by k bits within a 64-bit word. |
| `func RunLengthDecode[T comparable](runs []Run[T]) []T` | RunLengthDecode expands a sequence of runs back into the original slice. |
| `func RunLengthDecodeString(s string) (string, error)` | RunLengthDecodeString expands a string produced by RunLengthEncodeString back into its original form. |
| `func RunLengthEncodeString(s string) string` | RunLengthEncodeString compresses a string into the form "<count><rune>" for each maximal run (for example "aaabbc" becomes "3a2b1c"). |
| `func SetBit(x uint64, i uint) uint64` | SetBit returns x with bit i set to one. |
| `func TestBit(x uint64, i uint) bool` | TestBit reports whether bit i of x is set. |
| `func ToBase(n int64, base int) (string, error)` | ToBase returns the representation of the signed integer n in the given base, which must be in the range 2..36. |
| `func ToBaseUint(n uint64, base int) (string, error)` | ToBaseUint returns the representation of n in the given base, which must be in the range 2..36. |
| `func ToggleBit(x uint64, i uint) uint64` | ToggleBit returns x with bit i flipped. |
| `func TrailingZeros(x uint64) int` | TrailingZeros returns the number of consecutive zero bits at the low end of x, counted using a De Bruijn perfect-hash table. |

### Types

| Type | What it is |
| --- | --- |
| `Pair` | Pair is an ordered two-element tuple, used as the element type of a Cartesian product of two (possibly heterogeneous) slices. |
| `Run` | Run is a single maximal run of an equal value produced by run-length encoding: Value repeated Count times. |
| `Set` | Set is a mathematical set of distinct comparable elements, represented as a map whose keys are the members. |

<details>
<summary><code>Pair</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func CartesianProduct[T, U any](a []T, b []U) []Pair[T, U]` | CartesianProduct returns the Cartesian product of a and b: every ordered pair (x, y) with x drawn from a and y from b. |

</details>

<details>
<summary><code>Run</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func RunLengthEncode[T comparable](s []T) []Run[T]` | RunLengthEncode compresses a slice into a sequence of runs, each recording a value and the number of times it repeats consecutively. |

</details>

<details>
<summary><code>Set</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Difference[T comparable](a, b Set[T]) Set[T]` | Difference returns a new set containing the elements of a that are not in b (the relative complement a \ b). |
| `func Intersection[T comparable](sets ...Set[T]) Set[T]` | Intersection returns a new set containing exactly the elements common to all of the given sets. |
| `func NewSet[T comparable](items ...T) Set[T]` | NewSet returns a new Set containing the given items, with duplicates removed. |
| `func SymmetricDifference[T comparable](a, b Set[T]) Set[T]` | SymmetricDifference returns a new set containing the elements that are in exactly one of a and b. |
| `func Union[T comparable](sets ...Set[T]) Set[T]` | Union returns a new set containing every element that appears in any of the given sets. |
| `func (s Set[T]) Add(items ...T)` | Add inserts the given items into the set. |
| `func (s Set[T]) Clone() Set[T]` | Clone returns an independent copy of the set. |
| `func (s Set[T]) Contains(item T) bool` | Contains reports whether item is a member of the set. |
| `func (s Set[T]) Elements() []T` | Elements returns the members of the set as a slice in unspecified order. |
| `func (s Set[T]) Equal(other Set[T]) bool` | Equal reports whether the set has exactly the same members as other. |
| `func (s Set[T]) IsDisjoint(other Set[T]) bool` | IsDisjoint reports whether the receiver and other share no common element. |
| `func (s Set[T]) IsSubset(other Set[T]) bool` | IsSubset reports whether every element of the receiver is also in other. |
| `func (s Set[T]) IsSuperset(other Set[T]) bool` | IsSuperset reports whether the receiver contains every element of other. |
| `func (s Set[T]) Len() int` | Len returns the number of elements in the set (its cardinality). |
| `func (s Set[T]) Remove(items ...T)` | Remove deletes the given items from the set. |

</details>

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/discretemath).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
