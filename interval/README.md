# interval — implements rigorous interval arithmetic over IEEE-754 double precision using outward

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/interval.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/interval)

Package interval implements rigorous interval arithmetic over IEEE-754 double
precision using outward directed rounding.

An `Interval` represents the set of real numbers {x : Lo <= x <= Hi}. Every
operation returns an enclosure: an interval that is guaranteed to contain the
true result of the corresponding real operation applied to every point of the
inputs. This makes the package suitable for verified numerics, where a
computed answer must come with a mathematically sound error bound rather than
a heuristic estimate.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/interval"
```

## Usage

This is the package's own `ExampleInterval_Intersect`, so it compiles and its output is
asserted on every `go test ./interval/`.

```go
fmt.Println(interval.New(1, 5).Intersect(interval.New(3, 9)))
	fmt.Println(interval.New(1, 2).Intersect(interval.New(5, 6)))
```

```
[3, 5]
[empty]
```

## Exported surface

### Types

| Type | What it is |
| --- | --- |
| `Interval` | Interval is a closed real interval [Lo, Hi] interpreted as the set of real numbers x with Lo <= x <= Hi. |
| `Matrix` | Matrix is a dense matrix whose entries are intervals. |

<details>
<summary><code>Interval</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Add(a, b Interval) Interval` | Add returns an enclosure of a + b. |
| `func Bisect(x Interval) (left, right Interval)` | Bisect splits x at its midpoint and returns the closed left and right halves, which share the midpoint as a common endpoint. |
| `func Cos(a Interval) Interval` | Cos returns an enclosure of cos(a). |
| `func Div(a, b Interval) Interval` | Div returns an enclosure of a / b. |
| `func DotVec(x, y []Interval) (Interval, error)` | DotVec returns an outward rounded enclosure of the dot product of interval vectors x and y. |
| `func Empty() Interval` | Empty returns the empty interval, which contains no real numbers. |
| `func Entire() Interval` | Entire returns the interval [-Inf, +Inf] representing the whole real line. |
| `func Exp(a Interval) Interval` | Exp returns an enclosure of e**a. |
| `func Hull(a, b Interval) Interval` | Hull returns the convex hull of a and b. |
| `func Intersect(a, b Interval) Interval` | Intersect returns the intersection of a and b. |
| `func Log(a Interval) Interval` | Log returns an enclosure of the natural logarithm of a. |
| `func Mul(a, b Interval) Interval` | Mul returns an enclosure of a * b. |
| `func New(lo, hi float64) Interval` | New returns the interval [lo, hi]. |
| `func Newton(f, df func(Interval) Interval, x Interval, tol float64, maxIter int) (Interval, bool)` | Newton runs the interval Newton iteration to enclose a root of f within the starting interval x. |
| `func Point(x float64) Interval` | Point returns the degenerate interval [x, x] containing only the single value x. |
| `func Pow(a Interval, y float64) Interval` | Pow returns an enclosure of a**y for a real exponent y. |
| `func Sin(a Interval) Interval` | Sin returns an enclosure of sin(a). |
| `func Sqrt(a Interval) Interval` | Sqrt returns an enclosure of the square root of a. |
| `func Sub(a, b Interval) Interval` | Sub returns an enclosure of a - b. |
| `func (a Interval) Abs() Interval` | Abs returns the interval of absolute values {\|x\| : x in a}. |
| `func (a Interval) Add(b Interval) Interval` | Add returns an enclosure of a + b, with the lower bound rounded down and the upper bound rounded up so the result rigorously contains every sum x + y… |
| `func (a Interval) AddFloat(f float64) Interval` | AddFloat returns an enclosure of a + f, outward rounded. |
| `func (a Interval) Atan() Interval` | Atan returns an enclosure of the arctangent of a. |
| `func (a Interval) Contains(x float64) bool` | Contains reports whether the real number x lies within the interval. |
| `func (a Interval) ContainsInterval(b Interval) bool` | ContainsInterval reports whether b is a subset of a (every point of b lies in a). |
| `func (a Interval) Cos() Interval` | Cos returns an enclosure of the cosine of a. |
| `func (a Interval) Cosh() Interval` | Cosh returns an enclosure of the hyperbolic cosine of a. |
| `func (a Interval) Div(b Interval) Interval` | Div returns an enclosure of a / b, outward rounded. |
| `func (a Interval) DivFloat(f float64) Interval` | DivFloat returns an enclosure of a / f, outward rounded. |
| `func (a Interval) Equal(b Interval) bool` | Equal reports whether a and b denote exactly the same set. |
| `func (a Interval) Exp() Interval` | Exp returns an enclosure of e raised to the power of a. |
| `func (a Interval) Hull(b Interval) Interval` | Hull returns the convex hull (union enclosure) of a and b, the smallest interval containing both. |
| `func (a Interval) IntPow(n int) Interval` | IntPow returns an enclosure of a raised to the integer power n. |
| `func (a Interval) Intersect(b Interval) Interval` | Intersect returns the intersection of a and b, the largest interval contained in both. |
| `func (a Interval) IsBounded() bool` | IsBounded reports whether both bounds of the interval are finite. |
| `func (a Interval) IsEmpty() bool` | IsEmpty reports whether the interval is empty, i.e. |
| `func (a Interval) IsEntire() bool` | IsEntire reports whether the interval is the whole real line [-Inf, +Inf]. |
| `func (a Interval) IsPoint() bool` | IsPoint reports whether the interval is degenerate, containing exactly one finite value (Lo == Hi and finite). |
| `func (a Interval) Log() Interval` | Log returns an enclosure of the natural logarithm of a. |
| `func (a Interval) Mag() float64` | Mag returns the magnitude of the interval, the supremum of \|x\| over all x in the interval (also called the norm). |
| `func (a Interval) Max(b Interval) Interval` | Max returns the elementwise maximum enclosure {max(x, y) : x in a, y in b}. |
| `func (a Interval) Midpoint() float64` | Midpoint returns a representable point near the centre (Lo + Hi) / 2. |
| `func (a Interval) Mig() float64` | Mig returns the mignitude of the interval, the infimum of \|x\| over all x in the interval. |
| `func (a Interval) Min(b Interval) Interval` | Min returns the elementwise minimum enclosure {min(x, y) : x in a, y in b}. |
| `func (a Interval) Mul(b Interval) Interval` | Mul returns an enclosure of a * b, outward rounded. |
| `func (a Interval) MulFloat(f float64) Interval` | MulFloat returns an enclosure of a * f, outward rounded. |
| `func (a Interval) Neg() Interval` | Neg returns -a, the reflection of the interval through the origin. |
| `func (a Interval) Overlaps(b Interval) bool` | Overlaps reports whether a and b share at least one common point. |
| `func (a Interval) Pow(y float64) Interval` | Pow returns an enclosure of a raised to the real power y, defined for a with positive lower bound via exp(y * log(a)). |
| `func (a Interval) Radius() float64` | Radius returns half the width of the interval, rounded up. |
| `func (a Interval) Recip() Interval` | Recip returns an enclosure of 1 / a. |
| `func (a Interval) Sin() Interval` | Sin returns an enclosure of the sine of a. |
| `func (a Interval) Sinh() Interval` | Sinh returns an enclosure of the hyperbolic sine of a. |
| `func (a Interval) Sqrt() Interval` | Sqrt returns an enclosure of the square root of a. |
| `func (a Interval) Square() Interval` | Square returns an enclosure of a*a. |
| `func (a Interval) String() string` | String formats the interval as "[Lo, Hi]" using the shortest decimal that round-trips each bound. |
| `func (a Interval) Sub(b Interval) Interval` | Sub returns an enclosure of a - b, outward rounded. |
| `func (a Interval) SubFloat(f float64) Interval` | SubFloat returns an enclosure of a - f, outward rounded. |
| `func (a Interval) Tan() Interval` | Tan returns an enclosure of the tangent of a, computed as Sin(a) / Cos(a). |
| `func (a Interval) Tanh() Interval` | Tanh returns an enclosure of the hyperbolic tangent of a. |
| `func (a Interval) Width() float64` | Width returns Hi - Lo, the length of the interval, rounded up so the result never underestimates the true width. |

</details>

<details>
<summary><code>Matrix</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Identity(n int) *Matrix` | Identity returns the n-by-n interval identity matrix, with degenerate 1s on the diagonal and 0s elsewhere. |
| `func MatrixFromFloat(rows, cols int, values []float64) (*Matrix, error)` | MatrixFromFloat builds a rows-by-cols interval matrix whose entries are the degenerate (point) intervals of the given row-major values. |
| `func NewMatrix(rows, cols int) *Matrix` | NewMatrix returns a rows-by-cols interval matrix with every entry set to the degenerate interval [0, 0]. |
| `func (m *Matrix) Add(n *Matrix) (*Matrix, error)` | Add returns the entrywise sum m + n as a new matrix, each entry an outward rounded interval sum. |
| `func (m *Matrix) At(i, j int) Interval` | At returns the interval stored at row i, column j (both 0-indexed). |
| `func (m *Matrix) Clone() *Matrix` | Clone returns a deep copy of the matrix that shares no storage with the receiver. |
| `func (m *Matrix) Cols() int` | Cols returns the number of columns in the matrix. |
| `func (m *Matrix) Equal(n *Matrix) bool` | Equal reports whether m and n have the same shape and every corresponding entry denotes the same interval. |
| `func (m *Matrix) MaxWidth() float64` | MaxWidth returns the largest entry width in the matrix, a scalar measure of how tight the enclosure is. |
| `func (m *Matrix) Midpoint() []float64` | Midpoint returns the real matrix of entrywise midpoints as a row-major slice of length Rows*Cols. |
| `func (m *Matrix) Mul(n *Matrix) (*Matrix, error)` | Mul returns the matrix product m * n as a new matrix. |
| `func (m *Matrix) MulVec(x []Interval) ([]Interval, error)` | MulVec returns the matrix-vector product m * x as a new interval vector. |
| `func (m *Matrix) Rows() int` | Rows returns the number of rows in the matrix. |
| `func (m *Matrix) Scale(s Interval) *Matrix` | Scale returns a new matrix with every entry multiplied by the interval s. |
| `func (m *Matrix) Set(i, j int, v Interval)` | Set stores the interval v at row i, column j (both 0-indexed). |
| `func (m *Matrix) Sub(n *Matrix) (*Matrix, error)` | Sub returns the entrywise difference m - n as a new matrix. |
| `func (m *Matrix) Transpose() *Matrix` | Transpose returns a new matrix that is the transpose of m. |

</details>

### Variables

`ErrDimension`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/interval).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
