# rootfind — implements polynomial and nonlinear root finding in pure Go

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/rootfind.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/rootfind)

Package rootfind implements polynomial and nonlinear root finding in pure Go.

The package is organized around two polynomial representations and a family of
solvers that operate on them or on arbitrary scalar functions.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/rootfind"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func AberthEhrlich(c CPoly, tol float64, maxIter int) ([]complex128, int, error)` | AberthEhrlich finds all complex roots of a polynomial simultaneously using the Aberth-Ehrlich iteration, a third-order-per-step refinement of… |
| `func AllRootsReal(p Poly, tol float64) bool` | AllRootsReal reports whether every root of p is real, i.e. |
| `func AnnulusBounds(p Poly) (lo, hi float64)` | AnnulusBounds returns a lower and an upper bound r0 <= \|z\| <= r1 that enclose every nonzero root z of p in a complex annulus. |
| `func Bairstow(p Poly, u0, v0, tol float64, maxIter int) (u, v float64, iters int, err error)` | Bairstow extracts a real quadratic factor x^2 + u*x + v from the real polynomial p, starting from initial guesses u0, v0, using Bairstow's method. |
| `func BairstowRoots(p Poly, tol float64, maxIter int) ([]complex128, error)` | BairstowRoots finds all roots of the real polynomial p by repeatedly applying `Bairstow` to peel off quadratic factors and deflating, dropping to a… |
| `func BracketOutward(f Func, a, b, factor float64, maxIter int) (lo, hi float64, err error)` | BracketOutward searches outward from the interval [a, b] by geometrically expanding it until f changes sign across the endpoints, returning the… |
| `func BudanFourierCount(p Poly, a, b float64) int` | BudanFourierCount returns V(a) - V(b), where V is the sign-variation count of the Fourier sequence of p. |
| `func BudanFourierUpperBound(p Poly) int` | BudanFourierUpperBound returns an upper bound on the number of real roots of p (with multiplicity) in the open positive axis (0, +inf), via the… |
| `func CPolyRoots(c CPoly) ([]complex128, error)` | CPolyRoots returns all complex roots of the complex polynomial c using the Aberth-Ehrlich method, sorted by real then imaginary part. |
| `func CauchyBound(p Poly) float64` | CauchyBound returns Cauchy's bound: every (real or complex) root z of p satisfies \|z\| <= 1 + max_{i<n} \|a_i / a_n\|, where a_n is the leading… |
| `func CoeffSignVariations(p Poly) int` | CoeffSignVariations returns the number of sign changes among the coefficients of p taken in order of ascending power (zeros skipped). |
| `func CompanionEigenvalues(p Poly) ([]complex128, error)` | CompanionEigenvalues returns all roots of the polynomial p as the eigenvalues of its companion matrix, computed with the Francis double-shift QR… |
| `func CompanionMatrix(p Poly) ([][]float64, error)` | CompanionMatrix returns the companion matrix of the polynomial p as a dense row-major n-by-n slice, where n is the degree of p. |
| `func CountComplexRoots(c CPoly) int` | CountComplexRoots returns the number of roots of c counted with multiplicity, which by the fundamental theorem of algebra equals its degree. |
| `func CountRealRoots(p Poly) int` | CountRealRoots returns the total number of distinct real roots of p (each counted once regardless of multiplicity) using Sturm's theorem over the… |
| `func CountRealRootsInInterval(p Poly, a, b float64) int` | CountRealRootsInInterval returns the number of distinct real roots of p in the half-open interval (a, b]. |
| `func DescartesNegativeBound(p Poly) int` | DescartesNegativeBound returns the Descartes upper bound on the number of negative real roots of p, obtained by applying the rule of signs to p(-x). |
| `func DescartesPositiveBound(p Poly) int` | DescartesPositiveBound returns the Descartes upper bound on the number of positive real roots of p (counted with multiplicity): the sign-variation… |
| `func DescartesRuleOfSigns(p Poly) (positive, negative int)` | DescartesRuleOfSigns returns both Descartes bounds at once: an upper bound on the number of positive real roots and an upper bound on the number of… |
| `func DiscriminantCubic(a, b, c, d float64) float64` | DiscriminantCubic returns the discriminant of the cubic a*x^3+b*x^2+c*x+d, 18abcd - 4b^3d + b^2c^2 - 4ac^3 - 27a^2d^2. |
| `func DiscriminantQuadratic(a, b, c float64) float64` | DiscriminantQuadratic returns the discriminant b^2 - 4ac of a*x^2 + b*x + c. |
| `func DistinctRealRoots(p Poly, tol float64) []float64` | DistinctRealRoots returns the distinct real roots of p, each listed once, refined to tolerance tol. |
| `func DurandKerner(c CPoly, tol float64, maxIter int) ([]complex128, int, error)` | DurandKerner finds all n complex roots of a degree-n polynomial simultaneously using the Durand-Kerner (Weierstrass) iteration. |
| `func DurandKernerWithInit(c CPoly, init []complex128, tol float64, maxIter int) ([]complex128, int, error)` | DurandKernerWithInit runs the Durand-Kerner iteration from a caller-supplied set of initial guesses, one per root. |
| `func FindBrackets(f Func, a, b float64, n int) [][2]float64` | FindBrackets scans the interval [a, b] on a uniform grid of n subintervals and returns every subinterval across which f changes sign. |
| `func FujiwaraBound(p Poly) float64` | FujiwaraBound returns Fujiwara's bound on the root moduli: \|z\| <= 2 * max_i ( \|a_{n-i}/a_n\| ^ (1/i) ), which is one of the sharpest simple upper… |
| `func HalleyComplex(c CPoly, x0 complex128, tol float64, maxIter int) (complex128, int, error)` | HalleyComplex finds a root of the complex polynomial c from x0 using Halley's cubically convergent method in complex arithmetic. |
| `func Horner(coeffs []float64, x float64) float64` | Horner evaluates the polynomial given by ascending-order real coefficients at x using Horner's method, without constructing a Poly. |
| `func HornerComplex(coeffs []complex128, x complex128) complex128` | HornerComplex evaluates the polynomial given by ascending-order complex coefficients at x using Horner's method. |
| `func IsolateRoots(p Poly) [][2]float64` | IsolateRoots returns a set of disjoint intervals, each containing exactly one distinct real root of p, covering all real roots. |
| `func IsolateRootsInterval(p Poly, a, b float64) [][2]float64` | IsolateRootsInterval returns disjoint isolating intervals for the distinct real roots of p that lie in (a, b], using recursive Sturm bisection. |
| `func KojimaBound(p Poly) float64` | KojimaBound returns Kojima's bound on the root moduli, formed from the ratios of consecutive coefficients: \|z\| <= 2 * max_i q_i, where q_1 =… |
| `func LagrangeBound(p Poly) float64` | LagrangeBound returns Lagrange's bound: \|z\| <= max(1, sum_{i<n} \|a_i/a_n\|) for every root z of p. |
| `func Laguerre(c CPoly, x0 complex128, tol float64, maxIter int) (complex128, int, error)` | Laguerre finds a single root of the polynomial c starting from x0 using Laguerre's method, which has cubic convergence for simple roots and is… |
| `func LaguerreRoots(c CPoly, tol float64, maxIter int) ([]complex128, error)` | LaguerreRoots finds all roots of the polynomial c by repeated application of `Laguerre` followed by deflation: each root found is divided out and the… |
| `func LowerRootBound(p Poly) float64` | LowerRootBound returns a positive lower bound on the moduli of the nonzero roots of p: no nonzero root satisfies \|z\| < LowerRootBound(p). |
| `func MaxResidual(c CPoly, roots []complex128) float64` | MaxResidual returns the largest modulus \|c(z)\| over the given candidate roots, a scale-free measure of how well a computed root set satisfies the… |
| `func Multiplicity(p Poly, r, tol float64) int` | Multiplicity returns the multiplicity of the real value r as a root of p, i.e. |
| `func MultiplicityComplex(c CPoly, r complex128, tol float64) int` | MultiplicityComplex returns the multiplicity of the complex value r as a root of the complex polynomial c, via the complex derivative test with… |
| `func NewtonComplex(c CPoly, x0 complex128, tol float64, maxIter int) (complex128, int, error)` | NewtonComplex finds a root of the complex polynomial c from the start x0 using Newton's method in complex arithmetic. |
| `func PolishComplexRoots(c CPoly, roots []complex128, tol float64, steps int) []complex128` | PolishComplexRoots refines each approximate root in roots with a few Newton steps on the complex polynomial c, improving accuracy after a global… |
| `func PolyRoots(p Poly) ([]complex128, error)` | PolyRoots returns all complex roots of the real polynomial p using the Aberth-Ehrlich method, the recommended general-purpose driver. |
| `func QuadraticRoots(a, b, c float64) (complex128, complex128)` | QuadraticRoots returns the two roots of the quadratic a*x^2 + b*x + c using a numerically stable formula that avoids catastrophic cancellation. |
| `func RealRootInterval(p Poly) (lo, hi float64)` | RealRootInterval returns a symmetric interval [-b, b] that contains all real roots of p, where b is the Lagrange bound. |
| `func RealRoots(p Poly, imagTol float64) ([]float64, error)` | RealRoots returns the real roots of p, that is the roots whose imaginary part is negligible relative to imagTol. |
| `func SeparateRoots(roots []complex128, imagTol float64) (reals []float64, complexes []complex128)` | SeparateRoots partitions a list of complex roots into those that are effectively real (\|imag\| <= imagTol) and those that are genuinely complex. |
| `func SignChange(f Func, a, b float64) bool` | SignChange reports whether f(a) and f(b) have strictly opposite signs, which guarantees a root of a continuous f in (a, b) by the intermediate value… |
| `func SignVariations(vals []float64) int` | SignVariations returns the number of sign changes in the sequence of values, skipping zeros. |
| `func SolveCubic(a, b, c, d float64) ([]complex128, error)` | SolveCubic returns the three roots of a*x^3 + b*x^2 + c*x + d as complex numbers, using Cardano's method with complex cube roots so the formula is… |
| `func SolveCubicReal(a, b, c, d float64) ([]float64, error)` | SolveCubicReal returns the real roots of the cubic a*x^3+b*x^2+c*x+d in ascending order, extracted from `SolveCubic` by keeping roots whose imaginary… |
| `func SolveLinear(a, b float64) (float64, error)` | SolveLinear returns the root of the linear equation a*x + b = 0. |
| `func SolveQuadraticReal(a, b, c float64) []float64` | SolveQuadraticReal returns the real roots of a*x^2 + b*x + c in ascending order: two values when the discriminant is positive, one (repeated) when it… |
| `func SortComplex(z []complex128) []complex128` | SortComplex returns a copy of z sorted by real part then imaginary part. |
| `func SturmCountRoots(seq []Poly, a, b float64) int` | SturmCountRoots returns the number of distinct real roots of the polynomial underlying the Sturm sequence seq that lie in the half-open interval (a,… |
| `func SturmRealRoots(p Poly, tol float64) []float64` | SturmRealRoots returns all distinct real roots of p refined to tolerance tol. |
| `func SturmRefine(seq []Poly, a, b, tol float64) float64` | SturmRefine refines an isolating interval [a, b] known to contain exactly one distinct real root of p down to width tol, using Sturm-count bisection. |
| `func SturmVariations(seq []Poly, x float64) int` | SturmVariations returns the number of sign changes in the Sturm sequence seq evaluated at x, ignoring terms that evaluate to zero. |
| `func SturmVariationsAtNegInf(seq []Poly) int` | SturmVariationsAtNegInf returns the sign-variation count of the Sturm sequence as x -> -infinity, determined from the leading coefficients and… |
| `func SturmVariationsAtPosInf(seq []Poly) int` | SturmVariationsAtPosInf returns the sign-variation count of the Sturm sequence as x -> +infinity, determined from the signs of the leading… |
| `func TotalRealRoots(p Poly, tol float64) int` | TotalRealRoots returns the number of real roots of p counted with multiplicity, summing the multiplicities reported by `RealRootsWithMultiplicity`. |

### Types

| Type | What it is |
| --- | --- |
| `CPoly` | CPoly is a complex polynomial stored as coefficients in ascending order of power: c`i` is the coefficient of x^i and c[len(c)-1] is the leading… |
| `ComplexRootMultiplicity` | ComplexRootMultiplicity pairs a complex root with its multiplicity. |
| `Func` | Func is a real scalar function of one real variable. |
| `Poly` | Poly is a real polynomial stored as coefficients in ascending order of power: p`i` is the coefficient of x^i, so p[len(p)-1] is the leading… |
| `Result` | Result records the outcome of a scalar root-finding iteration. |
| `RootMultiplicity` | RootMultiplicity pairs a real root value with its multiplicity in p. |
| `SquareFreeFactor` | SquareFreeFactor is one factor of a squarefree factorization: the polynomial Factor collects exactly the roots that occur with the given Multiplicity. |

<details>
<summary><code>CPoly</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func CFromRoots(roots ...complex128) CPoly` | CFromRoots builds the monic complex polynomial with the given roots, the product (x - r0)(x - r1)... |
| `func NewCPoly(coeffs ...complex128) CPoly` | NewCPoly returns a CPoly built from the given ascending-order coefficients. |
| `func (c CPoly) Add(q CPoly) CPoly` | Add returns the sum c+q. |
| `func (c CPoly) At(x complex128) complex128` | At is an alias for Eval, evaluating c(x). |
| `func (c CPoly) Clone() CPoly` | Clone returns an independent copy of c. |
| `func (c CPoly) Coeff(i int) complex128` | Coeff returns the coefficient of x^i, or 0 when i is out of range. |
| `func (c CPoly) Deflate(r complex128) (quo CPoly, remainder complex128)` | Deflate divides c by the linear factor (x - r), returning the quotient and the remainder c(r) via synthetic division. |
| `func (c CPoly) Degree() int` | Degree returns the degree of c, the largest index with a nonzero coefficient, or -1 for the zero polynomial. |
| `func (c CPoly) Derivative() CPoly` | Derivative returns the formal derivative c'(x). |
| `func (c CPoly) Eval(x complex128) complex128` | Eval evaluates c(x) using Horner's method. |
| `func (c CPoly) EvalDeriv(x complex128) (val, deriv complex128)` | EvalDeriv evaluates c(x) and c'(x) together in a single Horner sweep. |
| `func (c CPoly) EvalDeriv2(x complex128) (val, d1, d2 complex128)` | EvalDeriv2 evaluates c, c', and c” at x in one Horner sweep. |
| `func (c CPoly) IsMonic(tol float64) bool` | IsMonic reports whether the leading coefficient of c equals 1 within tol. |
| `func (c CPoly) IsZero() bool` | IsZero reports whether c is the zero polynomial. |
| `func (c CPoly) LeadingCoeff() complex128` | LeadingCoeff returns the leading coefficient, or 0 for the zero polynomial. |
| `func (c CPoly) Monic() (CPoly, error)` | Monic returns c divided by its leading coefficient. |
| `func (c CPoly) Mul(q CPoly) CPoly` | Mul returns the product c*q by convolution. |
| `func (c CPoly) Neg() CPoly` | Neg returns -c. |
| `func (c CPoly) NumTerms() int` | NumTerms returns the number of nonzero coefficients of c. |
| `func (c CPoly) Pow(k int) CPoly` | Pow returns c raised to the nonnegative integer power k by repeated squaring. |
| `func (c CPoly) Scale(s complex128) CPoly` | Scale returns the polynomial s*c. |
| `func (c CPoly) String() string` | String renders c in descending-power notation with parenthesized complex coefficients, for example "(1+2i)x^2 + (3+0i)". |
| `func (c CPoly) Sub(q CPoly) CPoly` | Sub returns the difference c-q. |
| `func (c CPoly) ToReal() Poly` | ToReal returns the real polynomial formed by taking the real part of each coefficient of c, discarding imaginary parts. |
| `func (c CPoly) Trim() CPoly` | Trim returns c with trailing zero coefficients removed; the underlying array is shared with c. |

</details>

<details>
<summary><code>ComplexRootMultiplicity</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func GroupComplexRoots(roots []complex128, tol float64) []ComplexRootMultiplicity` | GroupComplexRoots clusters a list of complex roots so that roots within distance tol of one another are merged into a single representative (their… |

</details>

<details>
<summary><code>Poly</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func DeflateRoots(p Poly, roots []float64) Poly` | DeflateRoots divides p successively by (x - r) for each r in roots, returning the final quotient. |
| `func FourierSequence(p Poly) []Poly` | FourierSequence returns the Fourier (derivative) sequence of p: p, p', p”, ... |
| `func FromRoots(roots ...float64) Poly` | FromRoots builds the monic real polynomial whose roots are exactly the given values, i.e. |
| `func FromRootsWithLead(lead float64, roots ...float64) Poly` | FromRootsWithLead builds the real polynomial lead*(x-r0)(x-r1)... |
| `func NewPoly(coeffs ...float64) Poly` | NewPoly returns a Poly built from the given ascending-order coefficients. |
| `func PolyFromDesc(coeffs ...float64) Poly` | PolyFromDesc returns a Poly from coefficients given in descending order of power, i.e. |
| `func SquareFree(p Poly) (Poly, error)` | SquareFree returns the squarefree part of p, that is p divided by gcd(p, p'). |
| `func SturmSequence(p Poly) []Poly` | SturmSequence returns the canonical Sturm chain of the polynomial p: |
| `func (p Poly) Add(q Poly) Poly` | Add returns the sum p+q. |
| `func (p Poly) At(x float64) float64` | At is an alias for Eval, evaluating p(x). |
| `func (p Poly) Clone() Poly` | Clone returns an independent copy of p. |
| `func (p Poly) Coeff(i int) float64` | Coeff returns the coefficient of x^i, or 0 when i is out of range. |
| `func (p Poly) CoeffsDesc() []float64` | CoeffsDesc returns the coefficients of p in descending order of power, highest-degree first, trimmed of leading zeros. |
| `func (p Poly) Compose(q Poly) Poly` | Compose returns the polynomial p(q(x)) formed by substituting q into p, using Horner's method over polynomial arithmetic. |
| `func (p Poly) DeflateReal(r float64) (quo Poly, remainder float64)` | DeflateReal divides p by the linear factor (x - r), returning the quotient and the remainder p(r). |
| `func (p Poly) Degree() int` | Degree returns the degree of p, that is the largest index whose coefficient is nonzero. |
| `func (p Poly) Derivative() Poly` | Derivative returns the formal derivative p'(x). |
| `func (p Poly) DivMod(d Poly) (q, r Poly, err error)` | DivMod divides p by d and returns the quotient q and remainder r satisfying p = q*d + r with deg(r) < deg(d). |
| `func (p Poly) Equal(q Poly, tol float64) bool` | Equal reports whether p and q are equal as polynomials, that is every coefficient agrees within tol after ignoring trailing zeros. |
| `func (p Poly) Eval(x float64) float64` | Eval evaluates p(x) using Horner's method, which is both fast and numerically stable. |
| `func (p Poly) EvalComplex(x complex128) complex128` | EvalComplex evaluates p at a complex argument using Horner's method. |
| `func (p Poly) EvalDeriv(x float64) (val, deriv float64)` | EvalDeriv evaluates p(x) and p'(x) simultaneously with a single Horner sweep, returning the value and the first derivative. |
| `func (p Poly) EvalDeriv2(x float64) (val, d1, d2 float64)` | EvalDeriv2 evaluates p, p', and p” at x in a single Horner sweep. |
| `func (p Poly) GCD(q Poly) Poly` | GCD returns a greatest common divisor of p and q using the Euclidean algorithm, normalized to be monic. |
| `func (p Poly) InfNorm() float64` | InfNorm returns the largest absolute coefficient of p (the sup norm of its coefficient vector). |
| `func (p Poly) Integral(c float64) Poly` | Integral returns an antiderivative of p whose constant term equals c. |
| `func (p Poly) IsConstant() bool` | IsConstant reports whether p has degree 0 or is the zero polynomial. |
| `func (p Poly) IsMonic(tol float64) bool` | IsMonic reports whether the leading coefficient of p equals 1 within tol. |
| `func (p Poly) IsZero() bool` | IsZero reports whether p is the zero polynomial, i.e. |
| `func (p Poly) L1Norm() float64` | L1Norm returns the sum of absolute values of the coefficients of p. |
| `func (p Poly) L2Norm() float64` | L2Norm returns the Euclidean norm of the coefficient vector of p. |
| `func (p Poly) LeadingCoeff() float64` | LeadingCoeff returns the coefficient of the highest-degree term, or 0 for the zero polynomial. |
| `func (p Poly) Monic() (Poly, error)` | Monic returns p divided by its leading coefficient, so the result has leading coefficient 1. |
| `func (p Poly) Mul(q Poly) Poly` | Mul returns the product p*q computed by convolution of the coefficient sequences. |
| `func (p Poly) Neg() Poly` | Neg returns -p. |
| `func (p Poly) NumTerms() int` | NumTerms returns the number of nonzero coefficients of p. |
| `func (p Poly) Pow(k int) Poly` | Pow returns p raised to the nonnegative integer power k by repeated squaring. |
| `func (p Poly) ProductOfRoots() float64` | ProductOfRoots returns the product of all roots of p counted with multiplicity, which by Vieta's formulas equals (-1)^n * a_0/a_n. |
| `func (p Poly) Quo(d Poly) (Poly, error)` | Quo returns the quotient of p divided by d, discarding the remainder. |
| `func (p Poly) ReflectX() Poly` | ReflectX returns the polynomial p(-x), whose positive roots are the negatives of the negative roots of p. |
| `func (p Poly) Rem(d Poly) (Poly, error)` | Rem returns the remainder of p divided by d. |
| `func (p Poly) Reverse() Poly` | Reverse returns the reversal (reciprocal) polynomial x^n * p(1/x), whose coefficient sequence is that of p reversed. |
| `func (p Poly) Scale(s float64) Poly` | Scale returns the polynomial s*p. |
| `func (p Poly) ShiftScale(a, b float64) Poly` | ShiftScale returns the polynomial p(a*x + b), the composition of p with the affine map x \|-> a*x + b. |
| `func (p Poly) String() string` | String renders p in conventional descending-power notation, for example "2x^3 - x + 5". |
| `func (p Poly) Sub(q Poly) Poly` | Sub returns the difference p-q. |
| `func (p Poly) SumOfRoots() float64` | SumOfRoots returns the sum of all roots of p counted with multiplicity, which by Vieta's formulas equals -a_{n-1}/a_n. |
| `func (p Poly) ToComplex() CPoly` | ToComplex converts a real polynomial to its complex counterpart. |
| `func (p Poly) Trim() Poly` | Trim returns p with trailing (high-order) zero coefficients removed. |

</details>

<details>
<summary><code>Result</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Bisection(f Func, a, b, tol float64, maxIter int) (Result, error)` | Bisection finds a root of f in the bracketing interval [a, b] by repeated halving. |
| `func Brent(f Func, a, b, tol float64, maxIter int) (Result, error)` | Brent finds a root of f in [a, b] using Brent's method, which combines bisection, the secant method, and inverse quadratic interpolation. |
| `func FalsePosition(f Func, a, b, tol float64, maxIter int) (Result, error)` | FalsePosition (regula falsi) finds a root in [a, b] by interpolating a secant line through the bracket endpoints and keeping the subinterval that… |
| `func FixedPoint(g Func, x0, tol float64, maxIter int) (Result, error)` | FixedPoint iterates x <- g(x) from x0 to find a fixed point of g, which is a root of f(x) = g(x) - x. |
| `func Halley(f, df, d2f Func, x0, tol float64, maxIter int) (Result, error)` | Halley finds a root of f using Halley's third-order method, given the first derivative df and second derivative d2f. |
| `func Illinois(f Func, a, b, tol float64, maxIter int) (Result, error)` | Illinois is the Illinois variant of the false-position method: when an endpoint is retained across consecutive iterations its function value is… |
| `func Newton(f, df Func, x0, tol float64, maxIter int) (Result, error)` | Newton finds a root of f using Newton's method, given the derivative df. |
| `func Ridders(f Func, a, b, tol float64, maxIter int) (Result, error)` | Ridders finds a root of f in the bracket [a, b] using Ridders' method, which applies an exponential correction to the false-position estimate and… |
| `func Secant(f Func, x0, x1, tol float64, maxIter int) (Result, error)` | Secant finds a root of f from two initial guesses x0 and x1 using the secant iteration, which approximates the derivative by a finite difference. |
| `func Steffensen(f Func, x0, tol float64, maxIter int) (Result, error)` | Steffensen finds a root of f using Steffensen's method, a derivative-free iteration that achieves quadratic convergence by estimating the derivative… |

</details>

<details>
<summary><code>RootMultiplicity</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func RealRootsWithMultiplicity(p Poly, tol float64) []RootMultiplicity` | RealRootsWithMultiplicity returns the distinct real roots of p together with their multiplicities. |

</details>

<details>
<summary><code>SquareFreeFactor</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func SquareFreeFactorization(p Poly) ([]SquareFreeFactor, error)` | SquareFreeFactorization returns the squarefree factorization of p using Yun's algorithm: it decomposes p (up to its leading constant) as a product… |

</details>

### Constants

`DefaultMaxIter`, `DefaultRootTol`, `DefaultTol`

### Variables

`ErrBadInput`, `ErrDegreeTooLow`, `ErrEmptyInterval`, `ErrNoBracket`, `ErrNoConvergence`, `ErrZeroDerivative`, `ErrZeroPolynomial`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/rootfind).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
