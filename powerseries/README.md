# powerseries — implements truncated formal power series over float64 coefficients together

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/powerseries.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/powerseries)

Package powerseries implements truncated formal power series over float64
coefficients together with the classical generating-function toolkit.

A Series holds a finite slice of coefficients c[0], c[1], … , c`n-1` and is
interpreted as the polynomial

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/powerseries"
```

## Usage

This is the package's own `ExampleLagrangeInversion`, so it compiles and its output is
asserted on every `go test ./powerseries/`.

```go
phi := powerseries.New(1, 2, 1, 0, 0, 0, 0)
	w := powerseries.LagrangeInversion(phi, 7)
	for i := 1; i < w.Len(); i++ {
		fmt.Printf("%.0f ", w.Coeff(i))
	}
	fmt.Println()
```

```
1 2 5 14 42 132
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func Bell(n int) float64` | Bell returns the n-th Bell number as a float64. |
| `func BellNumbers(n int) []float64` | BellNumbers returns the first n Bell numbers B_0, B_1, … as float64 values, computed with the Bell triangle. |
| `func Bernoulli(n int) float64` | Bernoulli returns the n-th Bernoulli number as a float64, using B_1 = −1/2. |
| `func BernoulliNumbers(n int) []float64` | BernoulliNumbers returns the first n Bernoulli numbers B_0, B_1, … as float64 values, using the convention B_1 = −1/2. |
| `func Catalan(n int) float64` | Catalan returns the n-th Catalan number as a float64. |
| `func CatalanNumbers(n int) []float64` | CatalanNumbers returns the first n Catalan numbers C_0, C_1, … as float64 values, where C_k = C(2k, k)/(k+1). |
| `func Derangement(n int) float64` | Derangement returns the n-th derangement number as a float64. |
| `func DerangementNumbers(n int) []float64` | DerangementNumbers returns the first n derangement counts D_0, D_1, … as float64 values, where D_k is the number of fixed-point-free permutations… |
| `func Fibonacci(n int) float64` | Fibonacci returns the n-th Fibonacci number as a float64 with F_0 = 0. |
| `func FibonacciNumbers(n int) []float64` | FibonacciNumbers returns the first n Fibonacci numbers F_0, F_1, … as float64 values with F_0 = 0 and F_1 = 1. |
| `func MotzkinNumbers(n int) []float64` | MotzkinNumbers returns the first n Motzkin numbers M_0, M_1, … as float64 values, using the recurrence (n+2)·M_n = (2n+1)·M_{n−1} +… |
| `func PartitionCount(n int) float64` | PartitionCount returns p(n), the number of integer partitions of n, as a float64. |
| `func PartitionNumbers(n int) []float64` | PartitionNumbers returns the first n integer-partition counts p(0), p(1), … as float64 values, computed with Euler's pentagonal-number recurrence. |

### Types

| Type | What it is |
| --- | --- |
| `Series` | Series is a truncated formal power series with float64 coefficients. |

<details>
<summary><code>Series</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func BellGF(n int) Series` | BellGF returns the length-n exponential generating function of the Bell numbers, exp(e^x − 1), whose degree-k coefficient is B_k/k! where B_k… |
| `func BernoulliGF(n int) Series` | BernoulliGF returns the length-n exponential generating function of the Bernoulli numbers, x/(e^x − 1) = Σ B_k·x^k/k!, using the convention B_1 =… |
| `func BinomialGF(alpha float64, n int) Series` | BinomialGF returns the length-n binomial series (1+x)^alpha, whose degree-k coefficient is the generalized binomial coefficient C(alpha, k). |
| `func CatalanGF(n int) Series` | CatalanGF returns the length-n ordinary generating function of the Catalan numbers, C(x) = (1 − √(1−4x))/(2x), computed with the series square… |
| `func Constant(c float64, n int) Series` | Constant returns the constant series c carrying n coefficients. |
| `func DerangementGF(n int) Series` | DerangementGF returns the length-n exponential generating function of the derangement numbers, e^{−x}/(1−x), whose degree-k coefficient is D_k/k!. |
| `func EGFFromSequence(a []float64) Series` | EGFFromSequence returns the exponential generating function whose degree-n coefficient is a`n`/n!. |
| `func EGFtoOGF(s Series) Series` | EGFtoOGF reinterprets the exponential generating function s, whose degree-n coefficient is a_n/n!, as the ordinary generating function of the same… |
| `func ExpGF(n int) Series` | ExpGF returns the length-n exponential series e^x = Σ x^k/k!, the exponential generating function of the all-ones sequence. |
| `func FibonacciGF(n int) Series` | FibonacciGF returns the length-n ordinary generating function of the Fibonacci numbers, x/(1 − x − x²), whose degree-k coefficient is F_k with… |
| `func FromFunc(f func(i int) float64, n int) Series` | FromFunc returns the length-n series whose degree-i coefficient is f(i). |
| `func FromSlice(coeffs []float64) Series` | FromSlice returns a series that copies the coefficients of the supplied slice. |
| `func GeometricGF(n int) Series` | GeometricGF returns the length-n ordinary generating function 1/(1−x) = Σ x^k, whose coefficients are all one. |
| `func GeometricRatioGF(r float64, n int) Series` | GeometricRatioGF returns the length-n ordinary generating function 1/(1−r·x) = Σ r^k·x^k. |
| `func HarmonicGF(n int) Series` | HarmonicGF returns the length-n ordinary generating function of the harmonic numbers, −log(1−x)/(1−x), whose degree-k coefficient is H_k =… |
| `func Ident(n int) Series` | Ident returns the identity series x carrying n coefficients. |
| `func LagrangeInversion(phi Series, n int) Series` | LagrangeInversion returns the length-n power series w(t) with zero constant term that solves the functional equation w = t·φ(w), where phi is the… |
| `func LagrangeInversionApply(phi, h Series, n int) Series` | LagrangeInversionApply returns the length-n series H(w(t)), where w is the solution of w = t·φ(w) with zero constant term. |
| `func Monomial(c float64, degree, n int) Series` | Monomial returns the series c·x^degree carrying n coefficients. |
| `func MotzkinGF(n int) Series` | MotzkinGF returns the length-n ordinary generating function of the Motzkin numbers, (1 − x − √(1−2x−3x²))/(2x²), whose degree-k… |
| `func New(coeffs ...float64) Series` | New returns the series whose coefficients are the supplied values in order of increasing degree, so New(1, 2, 3) is 1 + 2x + 3x². |
| `func OGFFromSequence(a []float64) Series` | OGFFromSequence returns the ordinary generating function whose degree-n coefficient is a`n`. |
| `func OGFtoEGF(s Series) Series` | OGFtoEGF reinterprets the ordinary generating function s, whose degree-n coefficient is the sequence value a_n, as the exponential generating… |
| `func One(n int) Series` | One returns the constant series 1 carrying n coefficients. |
| `func PartitionGF(n int) Series` | PartitionGF returns the length-n ordinary generating function of the integer partition numbers, Π_{k≥1} 1/(1−x^k), whose degree-m coefficient is… |
| `func Zero(n int) Series` | Zero returns the zero series carrying n coefficients (precision n). |
| `func (s Series) Add(other Series) Series` | Add returns the sum s + other. |
| `func (s Series) Atan() Series` | Atan returns arctan(s) as a truncated power series, valid for an arbitrary constant term. |
| `func (s Series) Clone() Series` | Clone returns an independent copy of the series. |
| `func (s Series) Coeff(i int) float64` | Coeff returns the coefficient of x^i. |
| `func (s Series) Coeffs() []float64` | Coeffs returns a fresh copy of the coefficient slice in order of increasing degree. |
| `func (s Series) Compose(other Series) Series` | Compose returns the composition s(other(x)), truncated to the precision of other. |
| `func (s Series) Cos() Series` | Cos returns cos(s) as a truncated power series, valid for an arbitrary constant term. |
| `func (s Series) Cosh() Series` | Cosh returns cosh(s) as a truncated power series, valid for an arbitrary constant term. |
| `func (s Series) Derivative() Series` | Derivative returns the formal derivative d/dx of the series. |
| `func (s Series) Div(other Series) Series` | Div returns the quotient s/other as a truncated power series. |
| `func (s Series) Equal(other Series, tol float64) bool` | Equal reports whether s and other agree, coefficient by coefficient, to within an absolute tolerance of tol. |
| `func (s Series) Evaluate(x float64) float64` | Evaluate returns the value of the truncated series at x using Horner's rule. |
| `func (s Series) Exp() Series` | Exp returns exp(s) as a truncated power series. |
| `func (s Series) Extend(n int) Series` | Extend returns the series padded with zero coefficients up to length n. |
| `func (s Series) Hadamard(other Series) Series` | Hadamard returns the coefficient-wise (Hadamard) product of s and other, whose degree-i coefficient is s.Coeff(i)·other.Coeff(i). |
| `func (s Series) Integral() Series` | Integral returns the antiderivative of the series whose constant term is zero. |
| `func (s Series) IntegralConst(c float64) Series` | IntegralConst returns the antiderivative of the series whose constant term is c. |
| `func (s Series) Inverse() Series` | Inverse returns the multiplicative inverse 1/s as a truncated power series at the receiver's precision. |
| `func (s Series) IsZero(tol float64) bool` | IsZero reports whether every coefficient has magnitude at most tol. |
| `func (s Series) Len() int` | Len reports the number of stored coefficients, that is the working precision of the series. |
| `func (s Series) Log() Series` | Log returns the natural logarithm log(s) as a truncated power series. |
| `func (s Series) Mul(other Series) Series` | Mul returns the Cauchy product s · other, truncated to the precision of the longer operand. |
| `func (s Series) Neg() Series` | Neg returns −s. |
| `func (s Series) Order(tol float64) int` | Order returns the valuation of the series, that is the index of its lowest non-zero coefficient. |
| `func (s Series) Pow(k int) Series` | Pow returns s raised to the non-negative integer power k, using exponentiation by squaring. |
| `func (s Series) PowReal(p float64) Series` | PowReal returns s raised to the real power p as a truncated power series. |
| `func (s Series) Reversion() Series` | Reversion returns the compositional inverse g of the series, the unique series with zero constant term satisfying s(g(x)) = g(s(x)) = x. |
| `func (s Series) Scale(c float64) Series` | Scale returns c·s, multiplying every coefficient by the scalar c. |
| `func (s Series) SequenceFromEGF() []float64` | SequenceFromEGF returns the underlying sequence a_0, a_1, … of the series read as an exponential generating function, multiplying the degree-n… |
| `func (s Series) SequenceFromOGF() []float64` | SequenceFromOGF returns the underlying sequence a_0, a_1, … of the series read as an ordinary generating function, that is a plain copy of the… |
| `func (s Series) Shift(k int) Series` | Shift returns x^k · s, moving every coefficient up by k degrees and dropping terms pushed beyond the current precision. |
| `func (s Series) Sin() Series` | Sin returns sin(s) as a truncated power series, valid for an arbitrary constant term. |
| `func (s Series) Sinh() Series` | Sinh returns sinh(s) as a truncated power series, valid for an arbitrary constant term. |
| `func (s Series) Sqrt() Series` | Sqrt returns the principal square root of the series as a truncated power series. |
| `func (s Series) String() string` | String renders the series in the conventional a + b·x + c·x² notation, ending with an O(x^n) remainder term recording the precision. |
| `func (s Series) Sub(other Series) Series` | Sub returns the difference s − other. |
| `func (s Series) Tan() Series` | Tan returns tan(s) as a truncated power series, computed as Sin/Cos. |
| `func (s Series) Truncate(n int) Series` | Truncate returns the series keeping only the coefficients of degree below n. |

</details>

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/powerseries).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
