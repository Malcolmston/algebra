# contfrac — self-contained toolkit for continued fractions and Diophantine approximation

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/contfrac.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/contfrac)

Package contfrac is a self-contained toolkit for continued fractions and
Diophantine approximation, built entirely on the Go standard library (math,
math/big, sort, errors, fmt, strings).

- `Frac` is an exact rational with int64 numerator and denominator, with the
usual field of arithmetic and ordering operations plus conversions to and from
continued fractions. - `CF` is a finite continued fraction [a0; a1, a2, ...]
stored as a slice of partial quotients. The Euclidean/floor convention is
used, so a0 may be zero or negative but every later term is a positive
integer. - `QuadraticSurd` represents (P + sqrt(D))/Q and `PeriodicCF`
represents an eventually periodic continued fraction; the two convert to one
another, realising the classical theorem of Lagrange.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/contfrac"
```

## Usage

This is the package's own `ExampleBestApproximation`, so it compiles and its output is
asserted on every `go test ./contfrac/`.

```go
p, q := contfrac.BestApproximation(3.141592653589793, 113)
	fmt.Printf("%d/%d\n", p, q)
```

```
355/113
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func AreFareyNeighbors(a, b, c, d int64) bool` | AreFareyNeighbors reports whether a/b and c/d are adjacent in some Farey sequence, i.e. |
| `func BestApproximation(x float64, maxDen int64) (int64, int64)` | BestApproximation returns the rational p/q closest to x whose denominator q satisfies 1 <= q <= maxDen. |
| `func BestApproximationRat(target *big.Rat, maxDen int64) *big.Rat` | BestApproximationRat returns the rational closest to the target *big.Rat whose denominator does not exceed maxDen, computed exactly from the… |
| `func CFToPath(cf CF) string` | CFToPath converts a continued fraction of a positive value into the corresponding Stern-Brocot L/R path. |
| `func Continuant(a ...int64) int64` | Continuant returns the continuant K(a0, a1, ..., a_{n-1}), the polynomial that equals the numerator of the convergent of the continued fraction with… |
| `func ContinuantBig(a []int64) *big.Int` | ContinuantBig is the arbitrary-precision form of `Continuant`. |
| `func ConvergentError(cf CF, k int) float64` | ConvergentError returns the absolute difference between the value of the full continued fraction and its k-th convergent, evaluated in floating point. |
| `func ECFTerm(k int) int64` | ECFTerm returns the k-th partial quotient (0-indexed) of the continued fraction of Euler's number e = [2; 1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8, ...]. |
| `func EgyptianFraction(p, q int64) ([]int64, error)` | EgyptianFraction returns the denominators of a sum of distinct unit fractions equal to the proper fraction p/q (with 0 < p/q < 1) using Fibonacci's… |
| `func EgyptianFractionBig(r *big.Rat) ([]*big.Int, error)` | EgyptianFractionBig returns the denominators of a distinct-unit-fraction decomposition of the proper fraction r using the greedy algorithm, with… |
| `func EgyptianFractionFrac(f Frac) ([]int64, error)` | EgyptianFractionFrac is the `Frac` form of `EgyptianFraction`. |
| `func EulerPhi(n int64) int64` | EulerPhi returns Euler's totient of n for n >= 1: the count of integers in [1, n] coprime to n. |
| `func Evaluate(cf CF) float64` | Evaluate is a free-function form of `CF.Value`. |
| `func EvaluateRat(cf CF) *big.Rat` | EvaluateRat is a free-function form of `CF.Rat`. |
| `func FareyLength(n int64) int64` | FareyLength returns the number of terms in the Farey sequence F_n, which is 1 + sum_{k=1}^{n} phi(k). |
| `func GCD(a, b int64) int64` | GCD returns the non-negative greatest common divisor of a and b using the Euclidean algorithm. |
| `func GoldenRatio() float64` | GoldenRatio returns the golden ratio (1 + sqrt(5))/2 as a float64. |
| `func IsPellSolution(D, x, y, n int64) bool` | IsPellSolution reports whether (x, y) satisfies x^2 - D*y^2 == n. |
| `func IsPerfectSquare(n int64) bool` | IsPerfectSquare reports whether n is a perfect square. |
| `func IsPeriodPalindrome(period []int64) bool` | IsPeriodPalindrome reports whether the given period (excluding its final element) reads the same forwards and backwards, the classical structure of… |
| `func IsUnitFraction(p, q int64) bool` | IsUnitFraction reports whether p/q reduces to 1/n for some positive integer n. |
| `func Isqrt(n int64) int64` | Isqrt returns the integer square root floor(sqrt(n)) for n >= 0 using a Newton iteration that is exact for every non-negative int64. |
| `func LCM(a, b int64) int64` | LCM returns the least common multiple of a and b as a non-negative value. |
| `func Mediant(a, b, c, d int64) (int64, int64)` | Mediant returns the mediant (a+c)/(b+d) of the fractions a/b and c/d as an unreduced numerator and denominator. |
| `func MediantRat(f, g *big.Rat) *big.Rat` | MediantRat returns the mediant of two *big.Rat values. |
| `func PellFundamental(D int64) (x, y *big.Int, ok bool)` | PellFundamental returns the fundamental (smallest positive) solution (x, y) of the Pell equation x^2 - D*y^2 = 1 for a non-square D > 1, together… |
| `func PellFundamentalInt(D int64) (x, y int64, ok bool)` | PellFundamentalInt returns the fundamental solution of x^2 - D*y^2 = 1 as int64 values when they fit, together with ok. |
| `func PellNegative(D int64) (x, y *big.Int, ok bool)` | PellNegative returns the fundamental solution (x, y) of the negative Pell equation x^2 - D*y^2 = -1, together with ok. |
| `func PellNegativeSolvable(D int64) bool` | PellNegativeSolvable reports whether x^2 - D*y^2 = -1 has a solution, i.e. |
| `func PellNthSolution(D int64, n int) (x, y *big.Int, ok bool)` | PellNthSolution returns the n-th positive solution (n >= 1) of x^2 - D*y^2 = 1, with n == 1 the fundamental solution. |
| `func PellSolutions(D int64, count int) [][2]*big.Int` | PellSolutions returns the first count positive solutions of x^2 - D*y^2 = 1, generated from the fundamental solution by the recurrence x_{k+1} =… |
| `func PiFloat(digits int) *big.Float` | PiFloat returns pi computed to roughly the requested number of decimal digits as a *big.Float. |
| `func ReduceFraction(p, q int64) (int64, int64)` | ReduceFraction returns p/q in lowest terms with a positive denominator. |
| `func SqrtCFPeriod(n int64) []int64` | SqrtCFPeriod returns just the repeating block of the continued fraction of sqrt(n) (nil for a perfect square). |
| `func SqrtCFPeriodLength(n int64) int` | SqrtCFPeriodLength returns the length of the period of the continued fraction of sqrt(n) (0 for a perfect square). |
| `func SternBrocotDepth(p, q int64) int` | SternBrocotDepth returns the depth of p/q in the Stern-Brocot tree (the root 1/1 has depth 0). |
| `func SternBrocotPath(p, q int64) string` | SternBrocotPath returns the path from the root 1/1 of the Stern-Brocot tree to the positive reduced fraction p/q, encoded as a string of 'L' (left)… |
| `func SumUnitFractions(dens []int64) *big.Rat` | SumUnitFractions returns the exact sum of the unit fractions 1/d for each d in dens as a *big.Rat. |

### Types

| Type | What it is |
| --- | --- |
| `CF` | CF is a finite continued fraction [a0; a1, a2, ...] represented by its slice of partial quotients. |
| `Frac` | Frac is an exact rational number with an int64 numerator and denominator. |
| `PeriodicCF` | PeriodicCF is an eventually periodic continued fraction consisting of a non-repeating Head followed by an infinitely repeated Period. |
| `QuadraticSurd` | QuadraticSurd represents the real number (P + sqrt(D))/Q with integer P, Q, D and D >= 0, Q != 0. |

<details>
<summary><code>CF</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func ECF(n int) CF` | ECF returns the first n partial quotients of the continued fraction of e. |
| `func FromBigFloat(x *big.Float, n int) CF` | FromBigFloat returns up to n partial quotients of the continued fraction of a *big.Float, using the float's own precision for each reciprocal step. |
| `func FromFloat(x float64, maxTerms int) CF` | FromFloat returns up to maxTerms partial quotients of the continued fraction of x, stopping early once the fractional part is negligible. |
| `func FromFloatEps(x float64, eps float64, maxTerms int) CF` | FromFloatEps returns up to maxTerms partial quotients of the continued fraction of x. |
| `func FromFrac(f Frac) CF` | FromFrac returns the continued-fraction expansion of f. |
| `func FromRational(p, q int64) CF` | FromRational returns the canonical continued-fraction expansion of the rational p/q using the Euclidean algorithm. |
| `func FromRationalBig(r *big.Rat) CF` | FromRationalBig returns the continued-fraction expansion of a *big.Rat as a CF of int64 terms. |
| `func GoldenRatioCF(n int) CF` | GoldenRatioCF returns the first n partial quotients of the continued fraction of the golden ratio phi = (1 + sqrt(5))/2 = [1; 1, 1, 1, ...], which is… |
| `func ParseCF(s string) (CF, error)` | ParseCF parses the bracket notation produced by `CF.String`, for example "[3; 7, 15, 1]" or "[3]". |
| `func PathToCF(path string) CF` | PathToCF converts a Stern-Brocot L/R path into the continued fraction of the fraction it names. |
| `func PiCF(n int) CF` | PiCF returns the first n partial quotients of the continued fraction of pi, which begins [3; 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, ...]. |
| `func SqrtCFExpand(n int64, count int) CF` | SqrtCFExpand returns the first count partial quotients of the continued fraction of sqrt(n). |
| `func SqrtTwoCF(n int) CF` | SqrtTwoCF returns the first n partial quotients of the continued fraction of sqrt(2) = [1; 2, 2, 2, ...]. |
| `func (c CF) Canonical() CF` | Canonical returns the canonical form of the continued fraction: a trailing term equal to one is absorbed into the preceding term ([..., x, 1] becomes… |
| `func (c CF) Clone() CF` | Clone returns an independent copy of c. |
| `func (c CF) Compare(d CF) int` | Compare compares the values of c and d, returning -1, 0 or +1. |
| `func (c CF) Convergent(k int) Frac` | Convergent returns the k-th convergent (0-indexed) as a reduced `Frac`. |
| `func (c CF) Convergents() []Frac` | Convergents returns all convergents p0/q0, p1/q1, ... |
| `func (c CF) ConvergentsBig() []*big.Rat` | ConvergentsBig returns all convergents as *big.Rat, correct for expansions of any length. |
| `func (c CF) Denominators() []int64` | Denominators returns the denominators of the successive convergents. |
| `func (c CF) Equal(d CF) bool` | Equal reports whether c and d have identical partial quotients. |
| `func (c CF) EqualValue(d CF) bool` | EqualValue reports whether c and d denote the same rational number, comparing their canonical forms. |
| `func (c CF) Frac() Frac` | Frac evaluates the continued fraction exactly to a reduced `Frac` via the convergent recurrence. |
| `func (c CF) IsValid() bool` | IsValid reports whether c is a well-formed expansion: it is non-empty and every partial quotient after the first is a positive integer. |
| `func (c CF) Len() int` | Len returns the number of partial quotients. |
| `func (c CF) Numerators() []int64` | Numerators returns the numerators of the successive convergents. |
| `func (c CF) Rat() *big.Rat` | Rat evaluates the continued fraction exactly to a *big.Rat, correct for expansions of any length. |
| `func (c CF) Reverse() CF` | Reverse returns the continued fraction with its partial quotients in reverse order. |
| `func (c CF) Semiconvergents() []Frac` | Semiconvergents returns the sequence of intermediate fractions generated along the way to the value: for each partial quotient a_i and each m = 1,… |
| `func (c CF) String() string` | String renders the continued fraction in the usual bracket notation "[a0; a1, a2, ...]" (or "`a0`" for a single term). |
| `func (c CF) Terms() []int64` | Terms returns a copy of the underlying slice of partial quotients. |
| `func (c CF) Truncate(n int) CF` | Truncate returns the continued fraction consisting of the first n partial quotients. |
| `func (c CF) Value() float64` | Value evaluates the continued fraction as a float64 using stable backward recurrence, which avoids the overflow that forward convergent recurrence… |

</details>

<details>
<summary><code>Frac</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func BestApproximationFrac(x float64, maxDen int64) Frac` | BestApproximationFrac returns `BestApproximation` as a reduced `Frac`. |
| `func BestApproximations(x float64, maxDen int64) []Frac` | BestApproximations returns the sequence of convergents of x whose denominators do not exceed maxDen. |
| `func Convergent(cf CF, k int) Frac` | Convergent is a free-function form of `CF.Convergent`. |
| `func Convergents(cf CF) []Frac` | Convergents is a free-function form of `CF.Convergents`. |
| `func EApprox(n int) Frac` | EApprox returns the n-th convergent of e (using n partial quotients) as a single best rational approximation. |
| `func EConvergents(n int) []Frac` | EConvergents returns the first n convergents of e as reduced fractions. |
| `func FareyNeighbors(a, b, n int64) (prev, next Frac)` | FareyNeighbors returns the predecessor and successor of a/b in F_n. |
| `func FareyNext(a, b, c, d, n int64) Frac` | FareyNext returns the term after c/d given the two consecutive Farey terms a/b and c/d of F_n, using the standard mediant recurrence k =… |
| `func FareyPredecessor(a, b, n int64) Frac` | FareyPredecessor returns the term immediately before a/b in the Farey sequence F_n. |
| `func FareySequence(n int64) []Frac` | FareySequence returns the Farey sequence F_n: all reduced fractions in the closed interval [0, 1] with denominator at most n, listed in ascending… |
| `func FareySuccessor(a, b, n int64) Frac` | FareySuccessor returns the term immediately after a/b in the Farey sequence F_n. |
| `func FracFromFloat(x float64) (Frac, bool)` | FracFromFloat returns the exact Frac whose value equals x, provided x is a dyadic rational representable in an int64 numerator and denominator,… |
| `func MediantFrac(f, g Frac) Frac` | MediantFrac returns the mediant of two fractions as a `Frac` (not reduced, mirroring the Stern-Brocot construction). |
| `func NewFrac(n, d int64) Frac` | NewFrac returns the fraction n/d reduced to lowest terms with a positive denominator. |
| `func PiApprox(n int) Frac` | PiApprox returns the continued-fraction convergent of pi obtained from n partial quotients as a single best rational approximation. |
| `func PiConvergents(n int) []Frac` | PiConvergents returns the first n convergents of pi as reduced fractions, including the classical 22/7 and 355/113. |
| `func Rationalize(x float64, maxDen int64) Frac` | Rationalize is an alias for `BestApproximationFrac`: it returns the simplest good rational approximation of x with denominator at most maxDen. |
| `func SmallestDenominatorWithin(x, eps float64) Frac` | SmallestDenominatorWithin returns the fraction with the smallest possible denominator whose value lies within eps of x. |
| `func SqrtConvergent(n int64, k int) Frac` | SqrtConvergent returns the k-th convergent (0-indexed) of sqrt(n). |
| `func SqrtConvergents(n int64, count int) []Frac` | SqrtConvergents returns the first count convergents of sqrt(n) as reduced fractions. |
| `func SternBrocotAncestors(p, q int64) []Frac` | SternBrocotAncestors returns the fractions on the path from the root down to (but not including) p/q, in order from the root. |
| `func SternBrocotChildren(p, q int64) (left, right Frac)` | SternBrocotChildren returns both children (left, right) of p/q. |
| `func SternBrocotFromPath(path string) Frac` | SternBrocotFromPath returns the fraction reached by following the given L/R path from the root of the Stern-Brocot tree. |
| `func SternBrocotLeftChild(p, q int64) Frac` | SternBrocotLeftChild returns the left child of p/q in the Stern-Brocot tree. |
| `func SternBrocotParent(p, q int64) Frac` | SternBrocotParent returns the parent of p/q in the Stern-Brocot tree. |
| `func SternBrocotRightChild(p, q int64) Frac` | SternBrocotRightChild returns the right child of p/q in the Stern-Brocot tree. |
| `func UnitFraction(d int64) Frac` | UnitFraction returns the unit fraction 1/d as a `Frac`. |
| `func (f Frac) Abs() Frac` | Abs returns \|f\|. |
| `func (f Frac) Add(g Frac) Frac` | Add returns the reduced sum f + g. |
| `func (f Frac) Ceil() int64` | Ceil returns the least integer not less than f. |
| `func (f Frac) Cmp(g Frac) int` | Cmp compares f and g, returning -1 if f < g, 0 if they are equal and +1 if f > g. |
| `func (f Frac) Div(g Frac) Frac` | Div returns the reduced quotient f / g. |
| `func (f Frac) Equal(g Frac) bool` | Equal reports whether f and g represent the same rational number. |
| `func (f Frac) Float() float64` | Float returns the value of f as a float64. |
| `func (f Frac) Floor() int64` | Floor returns the greatest integer not exceeding f. |
| `func (f Frac) Inv() Frac` | Inv returns the reciprocal 1/f. |
| `func (f Frac) IsInteger() bool` | IsInteger reports whether f is an integer. |
| `func (f Frac) IsZero() bool` | IsZero reports whether f == 0. |
| `func (f Frac) Less(g Frac) bool` | Less reports whether f < g. |
| `func (f Frac) Mediant(g Frac) Frac` | Mediant returns the mediant (f.Num+g.Num)/(f.Den+g.Den) of f and g. |
| `func (f Frac) Mul(g Frac) Frac` | Mul returns the reduced product f * g. |
| `func (f Frac) Neg() Frac` | Neg returns -f. |
| `func (f Frac) Pow(n int) Frac` | Pow returns f raised to the integer power n (which may be negative). |
| `func (f Frac) Rat() *big.Rat` | Rat returns f as a *big.Rat. |
| `func (f Frac) Reduce() Frac` | Reduce returns an equivalent fraction in lowest terms with a positive denominator. |
| `func (f Frac) Round() int64` | Round returns f rounded to the nearest integer, rounding halves away from zero. |
| `func (f Frac) Sign() int` | Sign returns -1, 0 or +1 according to the sign of f (assuming a positive denominator, as produced by the constructors). |
| `func (f Frac) String() string` | String renders f as "p/q", or just "p" when the denominator is 1. |
| `func (f Frac) Sub(g Frac) Frac` | Sub returns the reduced difference f - g. |

</details>

<details>
<summary><code>PeriodicCF</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func SqrtCF(n int64) PeriodicCF` | SqrtCF returns the continued fraction of sqrt(n) for n >= 0 as a PeriodicCF. |
| `func (pc PeriodicCF) Convergents(n int) []Frac` | Convergents returns the first n convergents of the expansion. |
| `func (pc PeriodicCF) Expand(n int) CF` | Expand returns the first n partial quotients of the (infinite) continued fraction: the head followed by as many copies of the period as needed. |
| `func (pc PeriodicCF) HeadLength() int` | HeadLength returns the length of the non-repeating prefix. |
| `func (pc PeriodicCF) IsPurelyPeriodic() bool` | IsPurelyPeriodic reports whether the expansion has an empty head and a non-empty period, so the whole continued fraction repeats from the start. |
| `func (pc PeriodicCF) IsRational() bool` | IsRational reports whether the expansion terminates (empty period). |
| `func (pc PeriodicCF) PeriodLength() int` | PeriodLength returns the length of the repeating block (0 for a rational). |
| `func (pc PeriodicCF) String() string` | String renders the expansion as "[h0; h1, ..., (p0, p1, ...)]" with the period parenthesised. |
| `func (pc PeriodicCF) Surd() QuadraticSurd` | Surd returns the quadratic surd represented by the periodic continued fraction, obtained by solving the quadratic fixed-point equation of the period… |
| `func (pc PeriodicCF) Terms(n int) []int64` | Terms is an alias for Expand returning a plain []int64. |
| `func (pc PeriodicCF) Value() float64` | Value returns the numeric value of the periodic continued fraction as a float64. |

</details>

<details>
<summary><code>QuadraticSurd</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewSurd(P, Q, D int64) QuadraticSurd` | NewSurd returns the quadratic surd (P + sqrt(D))/Q. It panics if Q == 0 or D < 0. |
| `func (s QuadraticSurd) CF() PeriodicCF` | CF returns the eventually periodic continued fraction of the surd. |
| `func (s QuadraticSurd) Conjugate() QuadraticSurd` | Conjugate returns the algebraic conjugate (P - sqrt(D))/Q, written in the canonical form with a non-negative radical by negating both P and Q. |
| `func (s QuadraticSurd) Floor() int64` | Floor returns the greatest integer not exceeding the surd, computed exactly with integer arithmetic (no floating-point rounding error). |
| `func (s QuadraticSurd) IsRational() bool` | IsRational reports whether the surd is actually a rational number, i.e. |
| `func (s QuadraticSurd) String() string` | String renders the surd as "(P + sqrt(D))/Q". |
| `func (s QuadraticSurd) Value() float64` | Value returns the surd as a float64. |

</details>

### Variables

`ErrZeroDenominator`, `ErrPerfectSquare`, `ErrNonPositive`, `ErrNotProper`, `ErrOverflow`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/contfrac).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
