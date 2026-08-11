# seq — implements classical integer sequences and figurate numbers

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/seq.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/seq)

Package seq implements classical integer sequences and figurate numbers.

- Linear recurrences and their relatives: the Fibonacci numbers computed by
the fast-doubling method, together with the Lucas, Pell, Pell-Lucas,
Jacobsthal, Jacobsthal-Lucas, Tribonacci, Tetranacci, Padovan and Perrin
sequences, a fully general constant-coefficient LinearRecurrence engine and
the Zeckendorf representation.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/seq"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func CenteredHexagonal(n int) uint64` | CenteredHexagonal returns the n-th centered hexagonal number 3n(n+1) + 1, also known as a hex number: 1, 7, 19, 37, 61, … n must be non-negative. |
| `func CenteredPolygonal(s, n int) uint64` | CenteredPolygonal returns the n-th centered s-gonal number s·n(n+1)/2 + 1, the count of dots in a centered polygonal arrangement with n rings around… |
| `func CenteredSquare(n int) uint64` | CenteredSquare returns the n-th centered square number 2n(n+1) + 1: 1, 5, 13, 25, 41, … n must be non-negative. |
| `func CenteredTriangular(n int) uint64` | CenteredTriangular returns the n-th centered triangular number 3n(n+1)/2 + 1: 1, 4, 10, 19, 31, … n must be non-negative. |
| `func CollatzMax(n uint64) uint64` | CollatzMax returns the largest value attained along the Collatz trajectory of n, including n itself. |
| `func CollatzNext(n uint64) uint64` | CollatzNext returns the successor of n under the Collatz (3n+1) map: n/2 when n is even and 3n+1 when n is odd. |
| `func CollatzSequence(n uint64) []uint64` | CollatzSequence returns the Collatz trajectory of n: the sequence starting at n and repeatedly applying CollatzNext until it reaches 1, inclusive of… |
| `func CollatzSteps(n uint64) int` | CollatzSteps returns the total stopping time of n: the number of applications of the Collatz map required to reach 1. |
| `func Cube(n int) uint64` | Cube returns the n-th cube number n³. |
| `func Fibonacci(n int) uint64` | Fibonacci returns the n-th Fibonacci number Fₙ, defined by F₀ = 0, F₁ = 1 and Fₙ = Fₙ₋₁ + Fₙ₋₂. |
| `func FibonacciBig(n int) *big.Int` | FibonacciBig returns the n-th Fibonacci number Fₙ as an arbitrary-precision integer, computed by the fast-doubling method in O(log n) big-integer… |
| `func FibonacciSequence(n int) []uint64` | FibonacciSequence returns the first n Fibonacci numbers F₀ … Fₙ₋₁. |
| `func GeneralizedPentagonal(n int) int64` | GeneralizedPentagonal returns n(3n−1)/2 evaluated at the signed index n, so that ranging n over 0, 1, −1, 2, −2, … yields the generalized… |
| `func Gnomonic(n int) uint64` | Gnomonic returns the n-th gnomonic number 2n−1, the n-th odd number and the difference between consecutive square numbers: 1, 3, 5, 7, 9, … n… |
| `func HappyNumbers(count int) []uint64` | HappyNumbers returns the first count happy numbers in increasing order, beginning with 1, 7, 10, 13, 19, 23, 28, … count must be non-negative. |
| `func Heptagonal(n int) uint64` | Heptagonal returns the n-th heptagonal number n(5n−3)/2: 0, 1, 7, 18, 34, 55, 81, … n must be non-negative. |
| `func Hexagonal(n int) uint64` | Hexagonal returns the n-th hexagonal number Hₙ = n(2n−1): 0, 1, 6, 15, 28, 45, 66, … n must be non-negative. |
| `func IntSqrt(v uint64) uint64` | IntSqrt returns the integer square root of v, that is the largest integer r with r*r <= v. |
| `func IsFibonacci(n uint64) bool` | IsFibonacci reports whether n is a Fibonacci number. |
| `func IsHappy(n uint64) bool` | IsHappy reports whether n is a happy number: starting from n and repeatedly replacing the value by the sum of the squares of its decimal digits, the… |
| `func IsHexagonal(x uint64) bool` | IsHexagonal reports whether x is a hexagonal number, that is whether (1+√(8x+1))/4 is a positive integer. |
| `func IsKaprekar(n uint64) bool` | IsKaprekar reports whether n is a Kaprekar number: the decimal representation of n² can be split into a left part and a non-zero right part that add… |
| `func IsPentagonal(x uint64) bool` | IsPentagonal reports whether x is a (non-negative-index) pentagonal number, that is whether (1+√(24x+1))/6 is a positive integer. |
| `func IsPerfectSquare(v uint64) bool` | IsPerfectSquare reports whether v is a perfect square. |
| `func IsPolygonal(s int, x uint64) bool` | IsPolygonal reports whether x is an s-gonal number for some non-negative index, i.e. |
| `func IsPronic(x uint64) bool` | IsPronic reports whether x is a pronic number, that is whether x = n(n+1) for some non-negative integer n. |
| `func IsSquare(x uint64) bool` | IsSquare reports whether x is a perfect square. |
| `func IsTriangular(x uint64) bool` | IsTriangular reports whether x is a triangular number, that is whether 8x+1 is a perfect square. |
| `func Jacobsthal(n int) uint64` | Jacobsthal returns the n-th Jacobsthal number Jₙ, defined by J₀ = 0, J₁ = 1 and Jₙ = Jₙ₋₁ + 2·Jₙ₋₂: 0, 1, 1, 3, 5, 11, 21, …… |
| `func JacobsthalLucas(n int) uint64` | JacobsthalLucas returns the n-th Jacobsthal-Lucas number, defined by the same recurrence with initial values 2 and 1: 2, 1, 5, 7, 17, 31, …… |
| `func JacobsthalSequence(n int) []uint64` | JacobsthalSequence returns the first n Jacobsthal numbers J₀ … Jₙ₋₁. |
| `func KaprekarNumbers(limit uint64) []uint64` | KaprekarNumbers returns every Kaprekar number in the closed interval [1, limit], in increasing order: 1, 9, 45, 55, 99, 297, 703, 999, … |
| `func KaprekarSequence(n int) []int` | KaprekarSequence returns the trajectory of the four-digit Kaprekar routine starting from n, inclusive of n, continuing until it reaches the fixed… |
| `func KaprekarStep(n int) int` | KaprekarStep applies one step of Kaprekar's routine to a four-digit number: the digits of n (zero-padded to four places) are arranged into the… |
| `func LinearRecurrence(coeffs, init []int64, n int) int64` | LinearRecurrence evaluates the term aₙ of the constant-coefficient linear recurrence |
| `func LinearRecurrenceBig(coeffs, init []*big.Int, n int) *big.Int` | LinearRecurrenceBig evaluates term aₙ of the constant-coefficient linear recurrence described by LinearRecurrence using arbitrary-precision integer… |
| `func LinearRecurrenceSequence(coeffs, init []int64, m int) []int64` | LinearRecurrenceSequence returns the first m terms a₀ … a₍ₘ₋₁₎ of the linear recurrence described by LinearRecurrence. |
| `func LookAndSaySequence(seed string, n int) []string` | LookAndSaySequence returns the first n terms of the look-and-say sequence generated from the seed string, inclusive of the seed. |
| `func LookAndSayStep(s string) string` | LookAndSayStep returns the look-and-say successor of the decimal digit string s: each maximal run of identical digits is replaced by its length… |
| `func Lucas(n int) uint64` | Lucas returns the n-th Lucas number Lₙ, defined by L₀ = 2, L₁ = 1 and Lₙ = Lₙ₋₁ + Lₙ₋₂. |
| `func LucasBig(n int) *big.Int` | LucasBig returns the n-th Lucas number Lₙ as an arbitrary-precision integer. |
| `func LucasSequence(n int) []uint64` | LucasSequence returns the first n Lucas numbers L₀ … Lₙ₋₁. |
| `func Octagonal(n int) uint64` | Octagonal returns the n-th octagonal number n(3n−2): 0, 1, 8, 21, 40, 65, 96, … n must be non-negative. |
| `func Padovan(n int) uint64` | Padovan returns the n-th Padovan number, defined by P₀ = P₁ = P₂ = 1 and Pₙ = Pₙ₋₂ + Pₙ₋₃: 1, 1, 1, 2, 2, 3, 4, 5, 7, 9, 12, 16,… |
| `func PadovanSequence(n int) []uint64` | PadovanSequence returns the first n Padovan numbers P₀ … Pₙ₋₁. |
| `func Pell(n int) uint64` | Pell returns the n-th Pell number Pₙ, defined by P₀ = 0, P₁ = 1 and Pₙ = 2·Pₙ₋₁ + Pₙ₋₂. |
| `func PellBig(n int) *big.Int` | PellBig returns the n-th Pell number Pₙ as an arbitrary-precision integer. |
| `func PellLucas(n int) uint64` | PellLucas returns the n-th Pell-Lucas (companion Pell) number Qₙ, defined by Q₀ = 2, Q₁ = 2 and Qₙ = 2·Qₙ₋₁ + Qₙ₋₂: 2, 2, 6, 14,… |
| `func PellSequence(n int) []uint64` | PellSequence returns the first n Pell numbers P₀ … Pₙ₋₁. |
| `func Pentagonal(n int) uint64` | Pentagonal returns the n-th pentagonal number Pₙ = n(3n−1)/2: 0, 1, 5, 12, 22, 35, 51, … n must be non-negative. |
| `func Pentatope(n int) uint64` | Pentatope returns the n-th pentatope (4-simplex) number n(n+1)(n+2)(n+3)/24, the four-dimensional analogue of the triangular and tetrahedral numbers:… |
| `func Perrin(n int) uint64` | Perrin returns the n-th Perrin number, defined by P₀ = 3, P₁ = 0, P₂ = 2 and Pₙ = Pₙ₋₂ + Pₙ₋₃: 3, 0, 2, 3, 2, 5, 5, 7, 10, 12,… |
| `func PerrinSequence(n int) []uint64` | PerrinSequence returns the first n Perrin numbers P₀ … Pₙ₋₁. |
| `func Polygonal(s, n int) uint64` | Polygonal returns the n-th s-gonal number |
| `func Pronic(n int) uint64` | Pronic returns the n-th pronic (oblong) number n(n+1), the product of two consecutive integers: 0, 2, 6, 12, 20, 30, … n must be non-negative. |
| `func RecamanSequence(n int) []int64` | RecamanSequence returns the first n terms of Recaman's sequence, defined by a₀ = 0 and, for k ≥ 1, |
| `func Square(n int) uint64` | Square returns the n-th square number n². |
| `func SquarePyramidal(n int) uint64` | SquarePyramidal returns the n-th square-pyramidal number n(n+1)(2n+1)/6, the number of spheres in a pyramid with a square base of side n: 0, 1, 5,… |
| `func StarNumber(n int) uint64` | StarNumber returns the n-th star number (centered 12-gonal number) 6n(n−1) + 1: 1, 13, 37, 73, 121, … n must be at least 1. |
| `func Tetrahedral(n int) uint64` | Tetrahedral returns the n-th tetrahedral (triangular-pyramidal) number n(n+1)(n+2)/6, the number of spheres in a tetrahedral stack of n layers: 0, 1,… |
| `func Tetranacci(n int) uint64` | Tetranacci returns the n-th Tetranacci number, defined by initial values 0, 0, 0, 1 and the four-term sum recurrence: 0, 0, 0, 1, 1, 2, 4, 8, 15, 29,… |
| `func Triangular(n int) uint64` | Triangular returns the n-th triangular number Tₙ = n(n+1)/2, the number of dots in a triangular arrangement with n dots on a side: 0, 1, 3, 6, 10,… |
| `func Tribonacci(n int) uint64` | Tribonacci returns the n-th Tribonacci number Tₙ, defined by T₀ = 0, T₁ = 0, T₂ = 1 and Tₙ = Tₙ₋₁ + Tₙ₋₂ + Tₙ₋₃: 0, 0, 1,… |
| `func TribonacciSequence(n int) []uint64` | TribonacciSequence returns the first n Tribonacci numbers T₀ … Tₙ₋₁. |
| `func Zeckendorf(n uint64) []uint64` | Zeckendorf returns the Zeckendorf representation of n: the unique set of non-consecutive Fibonacci numbers (each ≥ 1, drawn from 1, 2, 3, 5, 8,… |

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/seq).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
