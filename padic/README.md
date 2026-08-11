# padic — self-contained toolkit for p-adic numbers and p-adic arithmetic

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/padic.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/padic)

Package padic is a self-contained toolkit for p-adic numbers and p-adic
arithmetic, built entirely on the Go standard library (math, math/big,
math/rand, sort, errors, fmt, strings).

It is a sibling subpackage of github.com/malcolmston/algebra but does not
depend on it. Every value is expressed in terms of the type `Padic`, a p-adic
number carried to a finite, tracked precision.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/padic"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func AbsValueBigRat(p *big.Int, r *big.Rat) *big.Rat` | AbsValueBigRat returns the p-adic absolute value of a big.Rat as an exact big.Rat. |
| `func AbsValueFloat(p, a, b *big.Int) float64` | AbsValueFloat returns the p-adic absolute value p^(-v) of the rational a/b as a float64. |
| `func AbsValueRat(p, a, b *big.Int) *big.Rat` | AbsValueRat returns the p-adic absolute value \|a/b\|_p = p^(-v) as an exact big.Rat, where v is the valuation of the rational a/b. |
| `func Binomial(n, k int) *big.Int` | Binomial returns the binomial coefficient C(n, k) as a big.Int. |
| `func CRTPair(r1, m1, r2, m2 *big.Int) (*big.Int, *big.Int, error)` | CRTPair returns the unique residue x modulo m1*m2 with x == r1 mod m1 and x == r2 mod m2, for coprime moduli m1, m2. |
| `func ConvergenceRadiusSlope(p *big.Int, coeffs []*big.Rat) (*big.Rat, bool)` | ConvergenceRadiusSlope returns the "slope" that bounds the region of convergence of the power series: the series sum a_n x^n converges for arguments… |
| `func ConvergesAt(p *big.Int, coeffs []*big.Rat, xval int) bool` | ConvergesAt reports whether sum a_n x^n converges for an argument of valuation xval, i.e. |
| `func ConvergesInUnitBall(p *big.Int, coeffs []*big.Rat) bool` | ConvergesInUnitBall reports whether the power series with the given coefficients converges everywhere on the closed unit ball, i.e. |
| `func Coprime(a, b *big.Int) bool` | Coprime reports whether a and b share no common factor greater than 1. |
| `func DigitsToRat(p *big.Int, val int, digits []*big.Int) *big.Rat` | DigitsToRat reconstructs the rational p^val * sum(digits`i`*p^i) from a valuation and a slice of digits. |
| `func ExpandRational(p *big.Int, a, b *big.Int, count int) (int, []*big.Int, error)` | ExpandRational returns the p-adic expansion of the rational a/b: the valuation v and the first count digits d[0], d[1], ... |
| `func Factorial(n int) *big.Int` | Factorial returns n! as a big.Int for n >= 0. |
| `func GCD(a, b *big.Int) *big.Int` | GCD returns the greatest common divisor of a and b as a fresh big.Int. |
| `func HenselLift(f []*big.Int, r0, p *big.Int, prec int) (*big.Int, error)` | HenselLift lifts a simple root r0 of the integer polynomial f modulo p to a root modulo p^prec, returning the lifted root in [0, p^prec). |
| `func InvMod(a, m *big.Int) (*big.Int, error)` | InvMod returns the multiplicative inverse of a modulo m, in the range [0, m), or an error if a is not invertible modulo m. |
| `func IsEisenstein(p *big.Int, coeffs []*big.Int) bool` | IsEisenstein reports whether the integer polynomial coeffs (low-to-high) is Eisenstein at p: the leading coefficient is a unit, every lower… |
| `func IsPrime(n *big.Int) bool` | IsPrime reports whether n is a prime, using a strong probabilistic test (Baillie-PSW plus Miller-Rabin rounds) that is deterministic in practice for… |
| `func IsPrimePower(n *big.Int) (p *big.Int, k int, ok bool)` | IsPrimePower reports whether n is a prime power p^k (k >= 1); if so it also returns p and k. |
| `func IsQuadraticResidue(a, p *big.Int) bool` | IsQuadraticResidue reports whether a is a non-zero quadratic residue modulo the odd prime p. |
| `func IsRootOfUnityInt(p, a *big.Int, prec int) bool` | IsRootOfUnityInt reports whether the integer a modulo p^prec is a root of unity in Z_p, i.e. |
| `func IsTeichmuller(p, a *big.Int, prec int) bool` | IsTeichmuller reports whether the integer a, taken modulo p^prec, is its own Teichmuller representative, i.e. |
| `func JacobiSymbol(a, n *big.Int) int` | JacobiSymbol returns the Jacobi symbol (a/n) for an odd positive n, one of -1, 0, or +1. |
| `func LegendreSymbol(a, p *big.Int) int` | LegendreSymbol returns the Legendre symbol (a/p) for an odd prime p, one of -1, 0, or +1. |
| `func MillerRabin(n *big.Int, rounds int, rng *rand.Rand) bool` | MillerRabin reports whether n is probably prime using the given number of Miller-Rabin rounds with witnesses drawn from rng. |
| `func Mod(a, m *big.Int) *big.Int` | Mod returns a reduced into the canonical range [0, m). |
| `func MultiplicativeOrderModP(a, p *big.Int) int` | MultiplicativeOrderModP returns the multiplicative order of a modulo the prime p (the least e > 0 with a^e == 1 mod p), or 0 if a is divisible by p. |
| `func NextPrime(n *big.Int) *big.Int` | NextPrime returns the smallest prime strictly greater than n. |
| `func PPow(p *big.Int, k int) *big.Int` | PPow returns p raised to the non-negative power k as a fresh big.Int. |
| `func PolyAdd(a, b []*big.Int) []*big.Int` | PolyAdd returns the sum of two integer polynomials given low-to-high. |
| `func PolyContent(p *big.Int, a []*big.Int) *big.Int` | PolyContent returns the p-adic content of a: p raised to the minimum coefficient valuation, i.e. |
| `func PolyDegree(coeffs []*big.Int) int` | PolyDegree returns the degree of the integer polynomial coeffs (low-to-high), or -1 for the zero polynomial. |
| `func PolyDerivative(coeffs []*big.Int) []*big.Int` | PolyDerivative returns the formal derivative of an integer polynomial given low-to-high. |
| `func PolyEval(coeffs []*big.Int, x *big.Int) *big.Int` | PolyEval evaluates the integer polynomial coeffs (low-to-high) at x exactly. |
| `func PolyEvalRat(coeffs []*big.Int, x *big.Rat) *big.Rat` | PolyEvalRat evaluates an integer polynomial at a rational point x = n/d. |
| `func PolyIsZero(coeffs []*big.Int) bool` | PolyIsZero reports whether every coefficient is zero. |
| `func PolyModPk(a []*big.Int, p *big.Int, k int) []*big.Int` | PolyModPk reduces every coefficient of a modulo p^k into [0, p^k). |
| `func PolyMul(a, b []*big.Int) []*big.Int` | PolyMul returns the product of two integer polynomials given low-to-high. |
| `func PolyNeg(a []*big.Int) []*big.Int` | PolyNeg returns the negation of an integer polynomial. |
| `func PolyScale(a []*big.Int, c *big.Int) []*big.Int` | PolyScale multiplies every coefficient of a by the scalar c. |
| `func PolySub(a, b []*big.Int) []*big.Int` | PolySub returns a - b for integer polynomials given low-to-high. |
| `func PolyTrim(coeffs []*big.Int) []*big.Int` | PolyTrim removes trailing (high-degree) zero coefficients, returning the polynomial in canonical low-to-high form. |
| `func PolyValuation(p *big.Int, a []*big.Int) (int, bool)` | PolyValuation returns the minimum p-adic valuation over the non-zero coefficients of a, and false if a is the zero polynomial. |
| `func PowMod(base, exp, m *big.Int) *big.Int` | PowMod returns base^exp mod m using fast modular exponentiation. |
| `func PrevPrime(n *big.Int) *big.Int` | PrevPrime returns the largest prime strictly less than n, or nil if there is none (n <= 2). |
| `func PrimitiveRootModP(p *big.Int) *big.Int` | PrimitiveRootModP returns the smallest primitive root modulo the odd prime p, or nil if p is 2 (where the group is trivial) or has none. |
| `func RootsModP(f []*big.Int, p *big.Int) []*big.Int` | RootsModP returns all residues r in [0, p) with f(r) == 0 mod p, including multiple roots. |
| `func SameValuation(p, m, n *big.Int) bool` | SameValuation reports whether two integers have the same p-adic valuation. |
| `func SeriesValuationAt(p *big.Int, coeffs []*big.Rat, xval int) (int, bool)` | SeriesValuationAt returns a lower bound on the p-adic valuation of the value f(x) = sum a_n x^n, namely min_n (val(a_n) + n*xval), where xval is the… |
| `func SimpleRootsModP(f []*big.Int, p *big.Int) []*big.Int` | SimpleRootsModP returns the residues r in [0, p) with f(r) == 0 mod p and f'(r) != 0 mod p, i.e. |
| `func SplitInt(p, n *big.Int) (int, *big.Int)` | SplitInt returns the valuation v and unit u such that n = p^v * u with u coprime to p. |
| `func SqrtInt(p, a *big.Int, prec int) (*big.Int, error)` | SqrtInt returns a square root of the integer a modulo p^prec (an integer r in [0, p^prec) with r^2 == a mod p^prec), or an error when none exists. |
| `func SqrtModP(a, p *big.Int, rng *rand.Rand) (*big.Int, error)` | SqrtModP returns a square root of a modulo the odd prime p using the Tonelli-Shanks algorithm. |
| `func SqrtUnitMod(p, u *big.Int, prec int) (*big.Int, error)` | SqrtUnitMod returns a square root of the unit u modulo p^prec, i.e. |
| `func StrassmannBound(p *big.Int, coeffs []*big.Rat) (int, error)` | StrassmannBound returns the number of zeros, counted with multiplicity, that a convergent power series sum a_n x^n has in the closed unit ball of Q_p… |
| `func StrassmannBoundInts(p *big.Int, coeffs []*big.Int) (int, error)` | StrassmannBoundInts is StrassmannBound for integer coefficients. |
| `func StrassmannBoundPadics(coeffs []*Padic) (int, error)` | StrassmannBoundPadics is StrassmannBound for p-adic coefficients, using their exact valuations. |
| `func TeichmullerCharacter(p, a *big.Int, prec int) (*big.Int, error)` | TeichmullerCharacter returns the value omega(a) mod p^prec of the Teichmuller character on the integer a, a synonym for TeichmullerRep that reads… |
| `func TeichmullerRep(p, a *big.Int, prec int) (*big.Int, error)` | TeichmullerRep returns the Teichmuller representative of the unit a modulo p^prec: the unique (p-1)-th root of unity (for p = 2, the unique element… |
| `func UnitPartInt(p, n *big.Int) *big.Int` | UnitPartInt returns n with all factors of p removed, i.e. |
| `func ValuationBigRat(p *big.Int, r *big.Rat) int` | ValuationBigRat returns the p-adic valuation of a big.Rat. |
| `func ValuationFactorial(p *big.Int, n int) int` | ValuationFactorial returns the p-adic valuation of n! using Legendre's formula, v_p(n!) = sum_{i>=1} floor(n / p^i) = (n - s_p(n)) / (p - 1). |
| `func ValuationInt(p, n *big.Int) int` | ValuationInt returns the p-adic valuation of the integer n: the exponent of the highest power of p dividing n. |
| `func ValuationInt64(p, n int64) int` | ValuationInt64 returns the p-adic valuation of the machine integer n for the machine prime p (p >= 2). |
| `func ValuationRat(p, a, b *big.Int) int` | ValuationRat returns the p-adic valuation of the rational a/b, defined as ValuationInt(a) - ValuationInt(b). |
| `func WeierstrassDegree(p *big.Int, coeffs []*big.Rat) (int, error)` | WeierstrassDegree returns the Weierstrass degree of a convergent power series: the number of zeros in the open unit ball, equal to the largest index… |

### Types

| Type | What it is |
| --- | --- |
| `NPPoint` | NPPoint is a lattice point (X, Y) of a Newton polygon, where X is a coefficient index and Y is the p-adic valuation of that coefficient. |
| `NewtonPolygon` | NewtonPolygon is the lower convex hull of the points (i, val_p(a_i)) of a polynomial sum a_i x^i, ignoring indices with a_i == 0. |
| `Padic` | Padic is a p-adic number carried to a finite, tracked precision. |

<details>
<summary><code>NewtonPolygon</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewtonPolygonFromInts(p *big.Int, coeffs []*big.Int) *NewtonPolygon` | NewtonPolygonFromInts builds the Newton polygon of the integer polynomial coeffs (low-to-high) for the prime p. |
| `func NewtonPolygonFromPadics(coeffs []*Padic) (*NewtonPolygon, error)` | NewtonPolygonFromPadics builds the Newton polygon from p-adic coefficients (low-to-high), using their valuations. |
| `func NewtonPolygonFromRats(p *big.Int, coeffs []*big.Rat) *NewtonPolygon` | NewtonPolygonFromRats builds the Newton polygon of the rational-coefficient polynomial coeffs (low-to-high) for the prime p. |
| `func (np *NewtonPolygon) DistinctSlopes() []*big.Rat` | DistinctSlopes returns the distinct slopes of the Newton polygon in increasing order (each root valuation is the negative of a slope). |
| `func (np *NewtonPolygon) IsPure() bool` | IsPure reports whether the Newton polygon consists of a single segment, i.e. |
| `func (np *NewtonPolygon) LowerVertices() []NPPoint` | LowerVertices returns a copy of the vertices of the lower hull. |
| `func (np *NewtonPolygon) NumRoots() int` | NumRoots returns the total number of roots (with multiplicity) accounted for by the Newton polygon, equal to the horizontal width of the lower hull. |
| `func (np *NewtonPolygon) PrincipalSlope() *big.Rat` | PrincipalSlope returns the smallest (most negative) slope of the polygon, corresponding to the roots of largest valuation, or nil if undefined. |
| `func (np *NewtonPolygon) RootValuations() []*big.Rat` | RootValuations returns the p-adic valuations of the roots of the polynomial, each repeated by the horizontal length of its segment. |
| `func (np *NewtonPolygon) SegmentLengths() []int` | SegmentLengths returns the horizontal lengths of the hull's segments, i.e. |
| `func (np *NewtonPolygon) Slopes() []*big.Rat` | Slopes returns the slopes of the segments of the lower hull, in increasing order, as exact big.Rat values. |
| `func (np *NewtonPolygon) Valuation(x int) *big.Rat` | Valuation returns the polygon's value (minimum height of the hull) at a given integer abscissa x, interpolating linearly along the hull segments. |
| `func (np *NewtonPolygon) Width() int` | Width returns the horizontal extent of the Newton polygon. |

</details>

<details>
<summary><code>Padic</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func DigitsToPadic(p *big.Int, val int, digits []*big.Int) *Padic` | DigitsToPadic builds a p-adic number from a valuation and a slice of digits (least significant first), to relative precision len(digits). |
| `func EvalPolyPadic(f []*big.Int, x *Padic) (*Padic, error)` | EvalPolyPadic evaluates the integer polynomial f (low-to-high) at the p-adic number x. |
| `func FromBigInt(p *big.Int, n *big.Int, prec int) *Padic` | FromBigInt returns the p-adic value of the integer n to relative precision prec. |
| `func FromInt(p *big.Int, n int, prec int) *Padic` | FromInt is an alias for FromInt64 for callers that pass an int. |
| `func FromInt64(p *big.Int, n int64, prec int) *Padic` | FromInt64 returns the p-adic value of the machine integer n to relative precision prec. |
| `func FromRat(p *big.Int, r *big.Rat, prec int) (*Padic, error)` | FromRat returns the p-adic value of the big.Rat r to relative precision prec. |
| `func FromRational(p *big.Int, a, b *big.Int, prec int) (*Padic, error)` | FromRational returns the p-adic value of the rational a/b to relative precision prec. |
| `func FromString(p *big.Int, s string, prec int) (*Padic, error)` | FromString parses a decimal integer or rational ("a" or "a/b") and returns its p-adic value to relative precision prec. |
| `func HenselLiftPadic(f []*big.Int, r0, p *big.Int, prec int) (*Padic, error)` | HenselLiftPadic lifts a simple root r0 of f modulo p to a p-adic number to relative precision prec. |
| `func Log1p(t *Padic) (*Padic, error)` | Log1p returns log(1 + t) for a p-adic t with valuation at least 1. |
| `func New(p *big.Int, val int, unit *big.Int, prec int) (*Padic, error)` | New constructs the p-adic number p^val * unit known to relative precision prec (so absolute precision val+prec). |
| `func NewUnit(p, unit *big.Int, prec int) (*Padic, error)` | NewUnit constructs a p-adic unit (valuation 0) from an integer unit value to relative precision prec. |
| `func NewtonStepPadic(f []*big.Int, r *Padic) (*Padic, error)` | NewtonStepPadic performs one full-precision Newton step r' = r - f(r)/f'(r) on p-adic numbers, given f and f' evaluated via the integer coefficient… |
| `func One(p *big.Int, prec int) *Padic` | One returns the p-adic number 1 to relative precision prec. |
| `func PadicRoots(f []*big.Int, p *big.Int, prec int) ([]*Padic, error)` | PadicRoots returns all p-adic roots (to relative precision prec) obtainable from the simple roots of the integer polynomial f modulo p by Hensel… |
| `func RootOfUnity(p *big.Int, prec int) (*Padic, error)` | RootOfUnity returns a primitive (p-1)-th root of unity in Z_p to relative precision prec: the Teichmuller representative of a primitive root modulo p. |
| `func Teichmuller(p, a *big.Int, prec int) (*Padic, error)` | Teichmuller returns the Teichmuller representative of the integer a as a p-adic unit to relative precision prec. |
| `func Uniformizer(p *big.Int, prec int) *Padic` | Uniformizer returns the p-adic number p (the standard uniformizer) to relative precision prec. |
| `func Zero(p *big.Int, absPrec int) *Padic` | Zero returns the p-adic zero known to the given absolute precision. |
| `func (x *Padic) AbsValue() *big.Rat` | AbsValue returns the p-adic absolute value \|x\|_p = p^(-val) as an exact big.Rat. |
| `func (x *Padic) AbsolutePrecision() int` | AbsolutePrecision returns the exponent v such that x is known modulo p^v. |
| `func (x *Padic) Add(y *Padic) (*Padic, error)` | Add returns x + y. |
| `func (x *Padic) AddInt(n int64) *Padic` | AddInt returns x + n for an integer n, tracked at x's precision. |
| `func (x *Padic) BigInt() (*big.Int, error)` | BigInt returns the canonical integer representative of x in [0, p^absprec) when x is integral (valuation >= 0), or an error otherwise. |
| `func (x *Padic) ConstantResidue() (*big.Int, error)` | ConstantResidue returns x reduced modulo p as an integer in [0, p), which requires x to be integral (valuation >= 0). |
| `func (x *Padic) Copy() *Padic` | Copy returns an independent deep copy of x. |
| `func (x *Padic) Digit(i int) (*big.Int, error)` | Digit returns the coefficient of p^i in the standard p-adic expansion of x, an integer in [0, p). |
| `func (x *Padic) Digits() ([]*big.Int, error)` | Digits returns the p-adic digits of the unit part of x: a slice of length RelativePrecision, each entry in [0, p), least significant first. |
| `func (x *Padic) Distance(y *Padic) *big.Rat` | Distance returns the p-adic distance \|x - y\|_p between x and y as an exact big.Rat. |
| `func (x *Padic) Div(y *Padic) (*Padic, error)` | Div returns x / y, or an error if y is zero or the primes differ. |
| `func (x *Padic) Equal(y *Padic) bool` | Equal reports whether x and y agree as p-adic numbers to the smaller of their two absolute precisions. |
| `func (x *Padic) Exp() (*Padic, error)` | Exp returns the p-adic exponential of x, defined by sum_{n>=0} x^n / n!, which converges when val(x) > 1/(p-1): for odd p this means val(x) >= 1, and… |
| `func (x *Padic) ExpConverges() bool` | ExpConverges reports whether Exp(x) converges, i.e. |
| `func (x *Padic) Format(f fmt.State, verb rune)` | Format implements fmt.Formatter for the %v, %s and %d verbs by delegating to String, so Padic values print cleanly. |
| `func (x *Padic) FullDigits() ([]*big.Int, error)` | FullDigits returns the p-adic digits of x from the constant term up to (but not including) the absolute precision. |
| `func (x *Padic) Inv() (*Padic, error)` | Inv returns the multiplicative inverse of x, or an error if x is zero. |
| `func (x *Padic) IsCloseTo(y *Padic, n int) bool` | IsCloseTo reports whether x and y agree modulo p^n, i.e. |
| `func (x *Padic) IsIntegral() bool` | IsIntegral reports whether x lies in the ring of p-adic integers, i.e. |
| `func (x *Padic) IsOne() bool` | IsOne reports whether x equals 1 to its tracked precision. |
| `func (x *Padic) IsRootOfUnity() bool` | IsRootOfUnity reports whether x is a root of unity in Z_p, i.e. |
| `func (x *Padic) IsSquare() bool` | IsSquare reports whether x is a square in Q_p, to the extent its tracked precision allows. |
| `func (x *Padic) IsUnit() bool` | IsUnit reports whether x is a p-adic unit, i.e. |
| `func (x *Padic) IsZero() bool` | IsZero reports whether x is the p-adic zero to its tracked precision. |
| `func (x *Padic) Log() (*Padic, error)` | Log returns the p-adic logarithm of x, defined by the series |
| `func (x *Padic) LogConverges() bool` | LogConverges reports whether Log(x) converges, i.e. |
| `func (x *Padic) Mul(y *Padic) (*Padic, error)` | Mul returns x * y. |
| `func (x *Padic) MulInt(n int64) *Padic` | MulInt returns x * n for an integer n. |
| `func (x *Padic) MulPow(k int) *Padic` | MulPow returns x multiplied by p^k, i.e. |
| `func (x *Padic) Neg() *Padic` | Neg returns -x with the same precision as x. |
| `func (x *Padic) NormFloat() float64` | NormFloat returns the p-adic absolute value of x as a float64. |
| `func (x *Padic) Pow(n int) (*Padic, error)` | Pow returns x raised to the integer power n, using binary exponentiation. |
| `func (x *Padic) Precision() int` | Precision returns the absolute precision of x, a synonym for AbsolutePrecision that reads naturally at call sites. |
| `func (x *Padic) Prime() *big.Int` | Prime returns a copy of the prime p of x. |
| `func (x *Padic) Rat() *big.Rat` | Rat returns an exact rational representative of x, namely unit * p^val with unit in [0, p^prec). |
| `func (x *Padic) RationalReconstruction() (*big.Int, *big.Int, error)` | RationalReconstruction attempts to recover a rational a/b congruent to x modulo p^absprec, with \|a\|, \|b\| bounded by sqrt(p^absprec / 2). |
| `func (x *Padic) ReduceTo(absPrec int) *Padic` | ReduceTo returns x truncated to the given absolute precision. |
| `func (x *Padic) RelativePrecision() int` | RelativePrecision returns the number of known p-adic digits of the unit part of x. |
| `func (x *Padic) ShiftDown(k int) *Padic` | ShiftDown returns x / p^k for k >= 0 (a synonym for MulPow with negative k). |
| `func (x *Padic) ShiftUp(k int) *Padic` | ShiftUp returns x * p^k for k >= 0 (a synonym for MulPow with non-negative k). |
| `func (x *Padic) Sqrt() (*Padic, error)` | Sqrt returns a square root of the p-adic number x to relative precision equal to that of x, or ErrNoRoot if x is not a square in Q_p. |
| `func (x *Padic) SqrtBothRoots() ([]*Padic, error)` | SqrtBothRoots returns both square roots (r and -r) of x, or an error if x is not a square. |
| `func (x *Padic) Square() *Padic` | Square returns x*x. |
| `func (x *Padic) String() string` | String renders x as "unit*p^val + O(p^absprec)" or "O(p^absprec)" for zero. |
| `func (x *Padic) Sub(y *Padic) (*Padic, error)` | Sub returns x - y to the smaller of the two absolute precisions. |
| `func (x *Padic) SubInt(n int64) *Padic` | SubInt returns x - n for an integer n. |
| `func (x *Padic) Teichmuller() (*Padic, error)` | Teichmuller returns the Teichmuller representative of a p-adic unit x as a new p-adic number of valuation 0 (a root of unity). |
| `func (x *Padic) TeichmullerDigits() ([]*big.Int, error)` | TeichmullerDigits returns the Teichmuller representatives of the standard p-adic digits of x, i.e. |
| `func (x *Padic) Unit() *big.Int` | Unit returns the unit part of x (an integer coprime to p) reduced modulo p^prec, or 0 for the zero element. |
| `func (x *Padic) Valuation() int` | Valuation returns the p-adic valuation of x. |

</details>

### Variables

`ErrNotPrime`, `ErrZeroDivision`, `ErrNotInvertible`, `ErrPrimeMismatch`, `ErrNoRoot`, `ErrDomain`, `ErrPrecision`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/padic).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
