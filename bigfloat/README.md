# bigfloat — implements arbitrary-precision real analysis on top of math/big.Float

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/bigfloat.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/bigfloat)

Package bigfloat implements arbitrary-precision real analysis on top of
math/big.Float.

Everything in the package is computed to a caller-requested precision,
expressed in bits of mantissa. Each routine accepts a prec argument and
returns a freshly allocated *big.Float carrying exactly that precision,
rounded to nearest. Internally the routines work at a slightly higher "guard"
precision and, where the argument magnitude demands it, add extra bits so that
argument reduction does not lose significance; the final result is then
rounded down to the requested precision.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/bigfloat"
```

## Usage

This is the package's own `ExampleEllipticK`, so it compiles and its output is
asserted on every `go test ./bigfloat/`.

```go
k, err := bf.EllipticK(bf.Zero(128), 128)
	if err != nil {
		panic(err)
	}

	fmt.Println(k.Text('f', 20))
	fmt.Println(bf.HalfPi(128).Text('f', 20))
```

```
1.57079632679489661923
1.57079632679489661923
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func AGM(a, b *big.Float, prec uint) (*big.Float, error)` | AGM returns the arithmetic-geometric mean of a and b to prec bits. |
| `func Abs(x *big.Float, prec uint) *big.Float` | Abs returns \|x\| rounded to prec bits. |
| `func Acos(x *big.Float, prec uint) (*big.Float, error)` | Acos returns the arccosine of x to prec bits, in [0, pi]. |
| `func Acosh(x *big.Float, prec uint) (*big.Float, error)` | Acosh returns the inverse hyperbolic cosine log(x + sqrt(x^2-1)) to prec bits. |
| `func Acot(x *big.Float, prec uint) *big.Float` | Acot returns the arccotangent atan(1/x) to prec bits (principal value in (0, pi)); Acot(0) == pi/2. |
| `func Acoth(x *big.Float, prec uint) (*big.Float, error)` | Acoth returns the inverse hyperbolic cotangent (1/2)*log((x+1)/(x-1)) to prec bits. |
| `func Acsc(x *big.Float, prec uint) (*big.Float, error)` | Acsc returns the arccosecant asin(1/x) to prec bits. |
| `func Acsch(x *big.Float, prec uint) (*big.Float, error)` | Acsch returns the inverse hyperbolic cosecant asinh(1/x) to prec bits. |
| `func Add(x, y *big.Float, prec uint) *big.Float` | Add returns x + y rounded to prec bits. |
| `func AlmostEqual(x, y, tol *big.Float) bool` | AlmostEqual reports whether \|x-y\| <= tol. |
| `func Apery(prec uint) *big.Float` | Apery returns Apery's constant zeta(3) to prec bits using |
| `func Asec(x *big.Float, prec uint) (*big.Float, error)` | Asec returns the arcsecant acos(1/x) to prec bits. |
| `func Asech(x *big.Float, prec uint) (*big.Float, error)` | Asech returns the inverse hyperbolic secant acosh(1/x) to prec bits. |
| `func Asin(x *big.Float, prec uint) (*big.Float, error)` | Asin returns the arcsine of x to prec bits, in [-pi/2, pi/2]. |
| `func Asinh(x *big.Float, prec uint) *big.Float` | Asinh returns the inverse hyperbolic sine log(x + sqrt(x^2+1)) to prec bits, for any real x. |
| `func Atan(x *big.Float, prec uint) *big.Float` | Atan returns the arctangent of x to prec bits, in (-pi/2, pi/2). |
| `func Atan2(y, x *big.Float, prec uint) *big.Float` | Atan2 returns the angle of the point (x, y) in (-pi, pi], to prec bits. |
| `func Atanh(x *big.Float, prec uint) (*big.Float, error)` | Atanh returns the inverse hyperbolic tangent (1/2)*log((1+x)/(1-x)) to prec bits. |
| `func Bernoulli(n uint) *big.Rat` | Bernoulli returns the n-th Bernoulli number B_n as an exact *big.Rat, using the convention B_1 = -1/2. |
| `func BernoulliFloat(n uint, prec uint) *big.Float` | BernoulliFloat returns the n-th Bernoulli number rounded to prec bits. |
| `func Beta(a, b *big.Float, prec uint) (*big.Float, error)` | Beta returns the beta function B(a,b) = Gamma(a)Gamma(b)/Gamma(a+b) to prec bits, for a > 0 and b > 0. |
| `func Binomial(n, k uint, prec uint) *big.Float` | Binomial returns the binomial coefficient C(n,k) rounded to prec bits. |
| `func BinomialBig(n, k uint) *big.Int` | BinomialBig returns the binomial coefficient C(n,k) as an exact *big.Int. |
| `func BinomialReal(x, y *big.Float, prec uint) (*big.Float, error)` | BinomialReal returns the generalised binomial coefficient Gamma(x+1)/(Gamma(y+1)*Gamma(x-y+1)) to prec bits, for real x and y where all three gamma… |
| `func Catalan(prec uint) *big.Float` | Catalan returns Catalan's constant G to prec bits using the rapidly convergent identity |
| `func Cbrt(x *big.Float, prec uint) *big.Float` | Cbrt returns the real cube root of x to prec bits, for any real x (negative arguments yield negative roots). |
| `func Ceil(x *big.Float, prec uint) *big.Float` | Ceil returns the least integer >= x, to prec bits. |
| `func Clamp(x, lo, hi *big.Float, prec uint) *big.Float` | Clamp returns x confined to the interval [lo, hi], rounded to prec bits. |
| `func Cmp(x, y *big.Float) int` | Cmp compares x and y and returns -1, 0, or +1. |
| `func Copysign(x, y *big.Float, prec uint) *big.Float` | Copysign returns a value with the magnitude of x and the sign of y. |
| `func Cos(x *big.Float, prec uint) *big.Float` | Cos returns the cosine of x (radians) to prec bits. |
| `func Cosh(x *big.Float, prec uint) *big.Float` | Cosh returns the hyperbolic cosine of x to prec bits. |
| `func Cospi(x *big.Float, prec uint) *big.Float` | Cospi returns cos(pi*x) to prec bits, reducing x modulo 2 exactly. |
| `func Cot(x *big.Float, prec uint) *big.Float` | Cot returns the cotangent of x (radians) to prec bits. |
| `func Coth(x *big.Float, prec uint) (*big.Float, error)` | Coth returns the hyperbolic cotangent of x to prec bits. |
| `func Coversin(x *big.Float, prec uint) *big.Float` | Coversin returns the coversed sine 1 - sin(x) to prec bits. |
| `func Csc(x *big.Float, prec uint) *big.Float` | Csc returns the cosecant 1/sin(x) to prec bits. |
| `func Csch(x *big.Float, prec uint) (*big.Float, error)` | Csch returns the hyperbolic cosecant 1/sinh(x) to prec bits. |
| `func Cube(x *big.Float, prec uint) *big.Float` | Cube returns x*x*x rounded to prec bits. |
| `func Deg2Rad(x *big.Float, prec uint) *big.Float` | Deg2Rad converts x degrees to radians to prec bits. |
| `func Degree(prec uint) *big.Float` | Degree returns the size of one degree in radians, pi/180, to prec bits. |
| `func Digamma(x *big.Float, prec uint) (*big.Float, error)` | Digamma returns the digamma function psi(x) = Gamma'(x)/Gamma(x) to prec bits, for any real x that is not a non-positive integer, for which it… |
| `func Dim(x, y *big.Float, prec uint) *big.Float` | Dim returns max(x-y, 0) rounded to prec bits. |
| `func Div(x, y *big.Float, prec uint) *big.Float` | Div returns x / y rounded to prec bits. |
| `func DoubleFactorial(n uint, prec uint) *big.Float` | DoubleFactorial returns the double factorial n!! rounded to prec bits. |
| `func DoubleFactorialBig(n uint) *big.Int` | DoubleFactorialBig returns the double factorial n!! as an exact *big.Int. |
| `func E(prec uint) *big.Float` | E returns Euler's number e = exp(1) to prec bits. |
| `func EllipticE(k *big.Float, prec uint) (*big.Float, error)` | EllipticE returns the complete elliptic integral of the second kind E(k) to prec bits, via the descending AGM sequence |
| `func EllipticK(k *big.Float, prec uint) (*big.Float, error)` | EllipticK returns the complete elliptic integral of the first kind K(k) to prec bits, using K(k) = pi / (2 * AGM(1, sqrt(1-k^2))). |
| `func Equal(x, y *big.Float) bool` | Equal reports whether x and y are exactly equal. |
| `func EulerGamma(prec uint) *big.Float` | EulerGamma returns the Euler-Mascheroni constant gamma to prec bits, computed with the Brent-McMillan series |
| `func Exp(x *big.Float, prec uint) *big.Float` | Exp returns e**x to prec bits. |
| `func Exp10(x *big.Float, prec uint) *big.Float` | Exp10 returns 10**x to prec bits. |
| `func Exp2(x *big.Float, prec uint) *big.Float` | Exp2 returns 2**x to prec bits. |
| `func Expm1(x *big.Float, prec uint) *big.Float` | Expm1 returns e**x - 1 to prec bits, accurate even when x is small. |
| `func Exsec(x *big.Float, prec uint) *big.Float` | Exsec returns the exsecant sec(x) - 1 to prec bits. |
| `func FMA(x, y, z *big.Float, prec uint) *big.Float` | FMA returns x*y + z rounded to prec bits (computed at higher internal precision to reduce intermediate rounding). |
| `func Factorial(n uint, prec uint) *big.Float` | Factorial returns n! rounded to prec bits. |
| `func FactorialBig(n uint) *big.Int` | FactorialBig returns n! as an exact *big.Int. |
| `func FallingFactorial(x *big.Float, n uint, prec uint) *big.Float` | FallingFactorial returns x*(x-1)*...*(x-n+1) to prec bits. |
| `func Float64(x *big.Float) (float64, big.Accuracy)` | Float64 returns the value of x as a float64, together with the accuracy of the conversion. |
| `func Floor(x *big.Float, prec uint) *big.Float` | Floor returns the greatest integer <= x, to prec bits. |
| `func Frac(x *big.Float, prec uint) *big.Float` | Frac returns the fractional part x - Trunc(x), to prec bits. |
| `func Frexp(x *big.Float, prec uint) (frac *big.Float, exp int)` | Frexp breaks x into a normalised fraction f in [0.5,1) and exponent e such that x == f * 2**e. |
| `func FromBig(n *big.Int, prec uint) *big.Float` | FromBig returns the big.Int n as a *big.Float with the given precision. |
| `func FromFloat64(x float64, prec uint) *big.Float` | FromFloat64 returns the float64 x as a *big.Float with the given precision. |
| `func FromInt(n int64, prec uint) *big.Float` | FromInt returns the integer n as a *big.Float with the given precision. |
| `func FromRat(q *big.Rat, prec uint) *big.Float` | FromRat returns the rational q as a *big.Float with the given precision. |
| `func Gamma(x *big.Float, prec uint) (*big.Float, error)` | Gamma returns the gamma function Gamma(x) to prec bits, for any real x that is not a non-positive integer (a pole), for which it returns ErrPole. |
| `func GoldenRatio(prec uint) *big.Float` | GoldenRatio returns the golden ratio (1+sqrt(5))/2 to prec bits. |
| `func Half(prec uint) *big.Float` | Half returns 1/2 with the given precision. |
| `func HalfPi(prec uint) *big.Float` | HalfPi returns pi/2 to prec bits. |
| `func Haversin(x *big.Float, prec uint) *big.Float` | Haversin returns the haversine (1 - cos(x))/2 to prec bits. |
| `func Hypot(x, y *big.Float, prec uint) *big.Float` | Hypot returns sqrt(x*x + y*y), computed without undue overflow, to prec bits. |
| `func InvE(prec uint) *big.Float` | InvE returns 1/e to prec bits. |
| `func InvGoldenRatio(prec uint) *big.Float` | InvGoldenRatio returns 1/phi = (sqrt(5)-1)/2 to prec bits. |
| `func InvPi(prec uint) *big.Float` | InvPi returns 1/pi to prec bits. |
| `func InvSqrt2(prec uint) *big.Float` | InvSqrt2 returns 1/sqrt(2) to prec bits. |
| `func InvSqrtPi(prec uint) *big.Float` | InvSqrtPi returns 1/sqrt(pi) to prec bits. |
| `func IsInf(x *big.Float) bool` | IsInf reports whether x is an infinity. |
| `func IsInteger(x *big.Float) bool` | IsInteger reports whether x has no fractional part. |
| `func IsZero(x *big.Float) bool` | IsZero reports whether x is zero. |
| `func LambertW(x *big.Float, prec uint) (*big.Float, error)` | LambertW returns the principal branch W0(x) of the Lambert W function to prec bits, i.e. |
| `func Ldexp(x *big.Float, n int, prec uint) *big.Float` | Ldexp returns x * 2**n, to prec bits. |
| `func Lgamma(x *big.Float, prec uint) (val *big.Float, sign int, err error)` | Lgamma returns the natural logarithm of the absolute value of Gamma(x) and the sign of Gamma(x) (+1 or -1), to prec bits. |
| `func Ln10(prec uint) *big.Float` | Ln10 returns the natural logarithm of 10 to prec bits. |
| `func Ln2(prec uint) *big.Float` | Ln2 returns the natural logarithm of 2 to prec bits. |
| `func Ln2Pi(prec uint) *big.Float` | Ln2Pi returns the natural logarithm of 2*pi to prec bits. |
| `func LnPi(prec uint) *big.Float` | LnPi returns the natural logarithm of pi to prec bits. |
| `func Log(x *big.Float, prec uint) (*big.Float, error)` | Log returns the natural logarithm of x to prec bits. |
| `func Log10(x *big.Float, prec uint) (*big.Float, error)` | Log10 returns the base-10 logarithm of x to prec bits. |
| `func Log10E(prec uint) *big.Float` | Log10E returns log base 10 of e, i.e. |
| `func Log1p(x *big.Float, prec uint) (*big.Float, error)` | Log1p returns log(1+x) to prec bits, accurate for small x. |
| `func Log2(x *big.Float, prec uint) (*big.Float, error)` | Log2 returns the base-2 logarithm of x to prec bits. |
| `func Log2E(prec uint) *big.Float` | Log2E returns log base 2 of e, i.e. |
| `func LogAddExp(x, y *big.Float, prec uint) *big.Float` | LogAddExp returns log(e**x + e**y) to prec bits, computed stably. |
| `func LogBase(x, base *big.Float, prec uint) (*big.Float, error)` | LogBase returns the logarithm of x in the given base to prec bits. |
| `func LogBeta(a, b *big.Float, prec uint) (*big.Float, error)` | LogBeta returns log(B(a,b)) to prec bits, for a > 0 and b > 0. |
| `func Logit(x *big.Float, prec uint) (*big.Float, error)` | Logit returns log(x/(1-x)) to prec bits, the inverse of Sigmoid. |
| `func Max(x, y *big.Float, prec uint) *big.Float` | Max returns the larger of x and y (rounded to prec bits). |
| `func Min(x, y *big.Float, prec uint) *big.Float` | Min returns the smaller of x and y (rounded to prec bits). |
| `func Mod(x, y *big.Float, prec uint) *big.Float` | Mod returns the floating-point remainder x - y*Trunc(x/y), to prec bits. |
| `func Mul(x, y *big.Float, prec uint) *big.Float` | Mul returns x * y rounded to prec bits. |
| `func MustParse(s string, prec uint) *big.Float` | MustParse is like Parse but panics on a parse error. |
| `func Neg(x *big.Float, prec uint) *big.Float` | Neg returns -x rounded to prec bits. |
| `func New(prec uint) *big.Float` | New returns a new zero value with the given precision in bits. |
| `func One(prec uint) *big.Float` | One returns 1 with the given precision. |
| `func Parse(s string, prec uint) (*big.Float, error)` | Parse parses s (in any base-10 float syntax accepted by big.Float.SetString) into a *big.Float with the given precision. |
| `func Pi(prec uint) *big.Float` | Pi returns the ratio of a circle's circumference to its diameter, pi, to prec bits. |
| `func PlasticNumber(prec uint) *big.Float` | PlasticNumber returns the plastic number, the real root of x^3 = x + 1, to prec bits. |
| `func Pochhammer(x *big.Float, n uint, prec uint) *big.Float` | Pochhammer is an alias for RisingFactorial. |
| `func Pow(x, y *big.Float, prec uint) (*big.Float, error)` | Pow returns x**y to prec bits, defined as exp(y*log(x)) for x > 0. |
| `func PowInt(x *big.Float, n int64, prec uint) *big.Float` | PowInt returns x**n for integer n to prec bits, via binary exponentiation. |
| `func QuarterPi(prec uint) *big.Float` | QuarterPi returns pi/4 to prec bits. |
| `func Rad2Deg(x *big.Float, prec uint) *big.Float` | Rad2Deg converts x radians to degrees to prec bits. |
| `func Recip(x *big.Float, prec uint) *big.Float` | Recip returns 1 / x rounded to prec bits. |
| `func Remainder(x, y *big.Float, prec uint) *big.Float` | Remainder returns the IEEE 754 remainder x - y*Round(x/y), to prec bits. |
| `func RisingFactorial(x *big.Float, n uint, prec uint) *big.Float` | RisingFactorial returns the Pochhammer symbol x*(x+1)*...*(x+n-1) to prec bits. |
| `func Root(x *big.Float, n int, prec uint) (*big.Float, error)` | Root returns the real n-th root of x to prec bits. |
| `func Round(x *big.Float, prec uint) *big.Float` | Round returns x rounded to the nearest integer, ties away from zero. |
| `func Rsqrt(x *big.Float, prec uint) (*big.Float, error)` | Rsqrt returns the reciprocal square root 1/sqrt(x) to prec bits. |
| `func Sec(x *big.Float, prec uint) *big.Float` | Sec returns the secant 1/cos(x) to prec bits. |
| `func Sech(x *big.Float, prec uint) *big.Float` | Sech returns the hyperbolic secant 1/cosh(x) to prec bits. |
| `func Sigmoid(x *big.Float, prec uint) *big.Float` | Sigmoid returns the logistic function 1/(1+e**-x) to prec bits. |
| `func Sign(x *big.Float) int` | Sign returns -1, 0, or +1 according to the sign of x. |
| `func Signbit(x *big.Float) bool` | Signbit reports whether x is negative or negative zero. |
| `func SilverRatio(prec uint) *big.Float` | SilverRatio returns the silver ratio 1+sqrt(2) to prec bits. |
| `func Sin(x *big.Float, prec uint) *big.Float` | Sin returns the sine of x (radians) to prec bits. |
| `func SinCos(x *big.Float, prec uint) (sin, cos *big.Float)` | SinCos returns both the sine and cosine of x (radians) to prec bits. |
| `func Sinc(x *big.Float, prec uint) *big.Float` | Sinc returns the unnormalised cardinal sine sin(x)/x to prec bits, with Sinc(0) == 1. |
| `func Sinh(x *big.Float, prec uint) *big.Float` | Sinh returns the hyperbolic sine of x to prec bits. |
| `func Sinpi(x *big.Float, prec uint) *big.Float` | Sinpi returns sin(pi*x) to prec bits, reducing x modulo 2 exactly so it stays accurate for large x. |
| `func Softplus(x *big.Float, prec uint) *big.Float` | Softplus returns log(1+e**x) to prec bits, computed in a numerically stable way for large \|x\|. |
| `func Sqrt(x *big.Float, prec uint) (*big.Float, error)` | Sqrt returns the square root of x to prec bits. |
| `func Sqrt2(prec uint) *big.Float` | Sqrt2 returns the square root of 2 to prec bits. |
| `func Sqrt3(prec uint) *big.Float` | Sqrt3 returns the square root of 3 to prec bits. |
| `func Sqrt5(prec uint) *big.Float` | Sqrt5 returns the square root of 5 to prec bits. |
| `func SqrtPi(prec uint) *big.Float` | SqrtPi returns sqrt(pi) to prec bits. |
| `func SqrtTwoPi(prec uint) *big.Float` | SqrtTwoPi returns sqrt(2*pi) to prec bits. |
| `func Square(x *big.Float, prec uint) *big.Float` | Square returns x*x rounded to prec bits. |
| `func String(x *big.Float) string` | String formats x using big.Float's default text representation. |
| `func Sub(x, y *big.Float, prec uint) *big.Float` | Sub returns x - y rounded to prec bits. |
| `func Tan(x *big.Float, prec uint) *big.Float` | Tan returns the tangent of x (radians) to prec bits. |
| `func Tanh(x *big.Float, prec uint) *big.Float` | Tanh returns the hyperbolic tangent of x to prec bits. |
| `func Tanpi(x *big.Float, prec uint) *big.Float` | Tanpi returns tan(pi*x) to prec bits. |
| `func Text(x *big.Float, format byte, digits int) string` | Text formats x with the given format byte and number of digits, matching big.Float.Text. |
| `func Trunc(x *big.Float, prec uint) *big.Float` | Trunc returns the integer part of x (rounding toward zero), to prec bits. |
| `func Two(prec uint) *big.Float` | Two returns 2 with the given precision. |
| `func TwoInvPi(prec uint) *big.Float` | TwoInvPi returns 2/pi to prec bits. |
| `func TwoPi(prec uint) *big.Float` | TwoPi returns 2*pi to prec bits. |
| `func Versin(x *big.Float, prec uint) *big.Float` | Versin returns the versed sine 1 - cos(x) to prec bits. |
| `func Zero(prec uint) *big.Float` | Zero returns 0 with the given precision. |
| `func Zeta3(prec uint) *big.Float` | Zeta3 is an alias for Apery; it returns zeta(3) to prec bits. |

### Variables

`ErrDomain`, `ErrNegative`, `ErrPole`, `ErrParse`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/bigfloat).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
