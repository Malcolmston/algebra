# diffalgebra — provides symbolic tools for differential algebra and linear ordinary differential

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/diffalgebra.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/diffalgebra)

Package diffalgebra provides symbolic tools for differential algebra and
linear ordinary differential equations over the rational numbers.

- Exact rational numbers via math/big and dense univariate polynomials `Poly`
over Q with the full complement of ring operations: addition, multiplication,
Euclidean division, greatest common divisors, extended GCD, resultants,
discriminants, square-free factorisation, rational-root isolation, composition
and formal differentiation and integration.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/diffalgebra"
```

## Usage

This is the package's own `ExampleIntegrateRationalReal`, so it compiles and its output is
asserted on every `go test ./diffalgebra/`.

```go
f, _ := diffalgebra.NewRatFunc(diffalgebra.OnePoly(), diffalgebra.PolyFromInts(1, 0, 1))

	rt, _ := diffalgebra.IntegrateRational(f)
	fmt.Println("Rothstein-Trager:", rt, "all residues rational:", rt.AllResiduesRational)

	real, _ := diffalgebra.IntegrateRationalReal(f)
	fmt.Println("real form:", real, "complete:", real.Complete)
	fmt.Printf("value at x=1: %.6f\n", real.EvalFloat(1))
```

```
Rothstein-Trager: 0 all residues rational: false
real form: 2/sqrt(4)*atan((2*x)/sqrt(4)) complete: true
value at x=1: 0.785398
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func CharacteristicComplexRoots(coeffs []float64, seed int64) []complex128` | CharacteristicComplexRoots returns the complex characteristic roots of the coefficient vector (coeffs`i` the coefficient of y^(i)). |
| `func HermiteReduce(f RatFunc) (rational RatFunc, remaining RatFunc)` | HermiteReduce performs Hermite reduction of the rational function f, returning the elementary rational part of its integral together with a remaining… |
| `func KovacicSolveSecondOrder(a2, a1, a0 Poly) (RatFunc, KovacicResult, error)` | KovacicSolveSecondOrder is a convenience wrapper that reduces the general second-order equation a2 y” + a1 y' + a0 y = 0 to normal form and runs… |
| `func LinearlyIndependentPoly(ps []Poly) bool` | LinearlyIndependentPoly reports whether the polynomials are linearly independent over Q, tested by a non-vanishing Wronskian. |
| `func LinearlyIndependentRatFunc(fs []RatFunc) bool` | LinearlyIndependentRatFunc reports whether the rational functions are linearly independent over the constants, tested by a non-vanishing Wronskian. |
| `func RatFromFloat(f float64) *big.Rat` | RatFromFloat returns the exact rational value of the IEEE-754 float f. |
| `func RatFromFrac(p, q int64) *big.Rat` | RatFromFrac returns the rational number p/q. |
| `func RatFromInt(n int64) *big.Rat` | RatFromInt returns the rational number n/1. |
| `func RatFromString(s string) (*big.Rat, bool)` | RatFromString parses a rational number written as "p/q", a decimal, or an integer. |
| `func RatString(r *big.Rat) string` | RatString formats r as "p/q" (or "p" when the denominator is one). |
| `func RatToFloat(r *big.Rat) float64` | RatToFloat returns the nearest float64 to r. |
| `func RecurrenceValues(coeffs []float64, seeds []float64, count int) ([]float64, error)` | RecurrenceValues iterates the recurrence directly to produce the first count terms from the given seed values, for cross-checking a closed-form… |
| `func StructureTheoremExp(f, g Poly) (string, bool)` | StructureTheoremExp applies the Risch structure-theorem heuristic to decide whether f e^{g} is elementary and returns a human-readable description of… |
| `func WronskianMatrixPoly(ps []Poly) [][]Poly` | WronskianMatrixPoly returns the Wronskian matrix of the polynomials ps. |
| `func WronskianMatrixRatFunc(fs []RatFunc) [][]RatFunc` | WronskianMatrixRatFunc returns the Wronskian matrix of the functions fs: row i holds the (i)-th derivatives of every function. |

### Types

| Type | What it is |
| --- | --- |
| `AtanTerm` | AtanTerm is an arctangent summand of a real antiderivative, arising from an irreducible quadratic factor of the denominator whose residues are a… |
| `Derivation` | Derivation is the standard derivation d/dx on the differential field Q(x). |
| `ExpIntegrand` | ExpIntegrand pairs a polynomial coefficient with an exponential argument, representing f(x) e^{g(x)}. |
| `KovacicPole` | KovacicPole describes a rational pole c of r together with its order and the leading Laurent coefficient b (the coefficient of 1/(x-c)^order). |
| `KovacicResult` | KovacicResult is the outcome of the Kovacic algorithm for y” = r y. |
| `LogTerm` | LogTerm is a summand Coeff * log(Arg) of the logarithmic part of an integral, where Coeff is a rational residue and Arg is a monic polynomial. |
| `ODESolution` | ODESolution is the general homogeneous solution of a constant-coefficient linear ODE, described by its fundamental system of basis terms and the… |
| `ODETerm` | ODETerm is one basis solution of a constant-coefficient linear ODE, of the form x^Power * e^{Alpha x} * trig(Beta x). |
| `ODETermKind` | ODETermKind classifies a basis solution of a constant-coefficient linear ODE. |
| `Operator` | Operator is a linear differential operator sum_i a_i(x) D^i with rational-function coefficients, an element of the non-commutative Ore ring Q(x)`D`… |
| `PartialFractionTerm` | PartialFractionTerm is a summand A/(F^k) of a partial-fraction decomposition, where F is a monic square-free polynomial appearing to power Power and… |
| `Poly` | Poly is a dense univariate polynomial over the rational numbers Q. The zero polynomial is represented by an empty coefficient slice. |
| `RatFunc` | RatFunc is a rational function p/q over Q, kept in reduced form with a monic denominator. |
| `RationalIntegral` | RationalIntegral is the result of integrating a rational function: an elementary rational part plus a sum of logarithm terms. |
| `RealRationalIntegral` | RealRationalIntegral is the antiderivative of a rational function expressed over the reals: a rational part, real logarithms and arctangents. |
| `RecurrenceSolution` | RecurrenceSolution is the general solution of a constant-coefficient linear recurrence, described by its fundamental basis terms and characteristic… |
| `RecurrenceTerm` | RecurrenceTerm is one basis solution of a linear recurrence, of the form n^Power * Rho^n * trig(n*Theta) (with Theta zero for real roots, in which… |
| `RecurrenceTermKind` | RecurrenceTermKind classifies a basis solution of a constant-coefficient linear recurrence. |
| `RootCluster` | RootCluster is a complex root together with the multiplicity inferred by grouping numerically coincident roots. |
| `SquareFreeFactor` | SquareFreeFactor is one factor of a square-free factorisation: the monic square-free polynomial Factor appears to the power Mult in the original. |

<details>
<summary><code>AtanTerm</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (t AtanTerm) EvalFloat(x float64) float64` | EvalFloat evaluates the arctangent term numerically at x. |
| `func (t AtanTerm) String() string` | String renders the term as "c/sqrt(d)*atan(lin/sqrt(d))". |

</details>

<details>
<summary><code>Derivation</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func StandardDerivation() Derivation` | StandardDerivation returns the derivation d/dx on Q(x). |
| `func (Derivation) Apply(f RatFunc) RatFunc` | Apply returns the derivative of a rational function. |
| `func (d Derivation) ApplyN(f RatFunc, n int) RatFunc` | ApplyN returns the n-th derivative of a rational function (n >= 0). |
| `func (d Derivation) ApplyNPoly(p Poly, n int) Poly` | ApplyNPoly returns the n-th derivative of a polynomial (n >= 0). |
| `func (Derivation) ApplyPoly(p Poly) Poly` | ApplyPoly returns the derivative of a polynomial. |
| `func (d Derivation) ConstantField(f RatFunc) bool` | ConstantField reports whether every element derived from f by the derivation eventually vanishes; for Q(x) this is equivalent to f being constant. |
| `func (d Derivation) IsConstant(f RatFunc) bool` | IsConstant reports whether f lies in the kernel of the derivation, i.e. |
| `func (d Derivation) Leibniz(f, g RatFunc) (direct, leibniz RatFunc)` | Leibniz verifies the product rule for two rational functions, returning D(fg) computed both directly and via the Leibniz expansion f'g + fg'; the two… |
| `func (d Derivation) LogarithmicDerivative(f RatFunc) (RatFunc, error)` | LogarithmicDerivative returns f'/f, the logarithmic derivative of f. |

</details>

<details>
<summary><code>ExpIntegrand</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (e ExpIntegrand) EvalFloat(x float64) (float64, bool)` | EvalFloat evaluates the elementary antiderivative R(x) e^{g(x)} at x, when it exists; the second return reports elementary integrability. |
| `func (e ExpIntegrand) Integrate() (Poly, bool)` | Integrate returns the polynomial R such that the antiderivative equals R(x) e^{Arg(x)}, together with whether the integrand is elementary. |

</details>

<details>
<summary><code>KovacicPole</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (p KovacicPole) String() string` | String renders a pole for diagnostics. |

</details>

<details>
<summary><code>KovacicResult</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Kovacic(r RatFunc) (KovacicResult, error)` | Kovacic runs Case 1 of the Kovacic algorithm on y” = r y, attempting to construct a Liouvillian solution y = P exp(integral omega). |
| `func (k KovacicResult) EvalFloat(x float64) float64` | EvalFloat evaluates the constructed solution at x, using \|x-c\| inside the logarithms; it returns 0 when no solution was found. |
| `func (k KovacicResult) SolutionString() string` | SolutionString renders the constructed solution y = P * exp(integral Omega). |

</details>

<details>
<summary><code>LogTerm</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (t LogTerm) String() string` | String renders the log term as "c*log(arg)". |

</details>

<details>
<summary><code>ODESolution</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func SolveLinearConstantODE(coeffs []float64, seed int64, tol float64) (ODESolution, error)` | SolveLinearConstantODE builds the general solution of the homogeneous constant-coefficient linear ODE sum_i coeffs`i` y^(i) = 0, where coeffs`i` is… |
| `func SolveODEIVP(coeffs []float64, x0 float64, ic []float64, seed int64, tol float64) ([]float64, ODESolution, error)` | SolveODEIVP solves the initial-value problem for the homogeneous constant-coefficient ODE sum_i coeffs`i` y^(i) = 0 with the initial data ic`d` =… |
| `func (s ODESolution) Basis() []ODETerm` | Basis returns the fundamental system of basis terms. |
| `func (s ODESolution) Dimension() int` | Dimension returns the number of basis solutions (the order of the ODE). |
| `func (s ODESolution) EvalBasis(x float64) []float64` | EvalBasis returns the values of every basis term at x. |
| `func (s ODESolution) Evaluate(consts []float64, x float64) float64` | Evaluate returns sum_i consts`i` * basis_i(x). |
| `func (s ODESolution) String() string` | String renders the general solution as a linear combination C1*..+C2*.. |

</details>

<details>
<summary><code>ODETerm</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (t ODETerm) Eval(x float64) float64` | Eval returns the value of the basis term at x. |
| `func (t ODETerm) String() string` | String renders the term in readable form. |

</details>

<details>
<summary><code>Operator</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func ConstOperator(a RatFunc) Operator` | ConstOperator returns the zeroth-order operator that multiplies by the rational function a. |
| `func DOperator() Operator` | DOperator returns the derivation operator D itself. |
| `func IdentityOperator() Operator` | IdentityOperator returns the identity operator (multiplication by 1). |
| `func NewOperator(coeffs ...RatFunc) Operator` | NewOperator builds an operator from coefficients in ascending order of D (coeffs[0] is the zeroth-order term). |
| `func OperatorFromPolys(coeffs ...Poly) Operator` | OperatorFromPolys builds an operator whose coefficients are the given polynomials, in ascending order of D. |
| `func ZeroOperator() Operator` | ZeroOperator returns the zero operator. |
| `func (o Operator) Add(p Operator) Operator` | Add returns o+p. |
| `func (o Operator) Adjoint() Operator` | Adjoint returns the formal adjoint operator L^* defined by L^*(y) = sum_i (-1)^i (a_i y)^(i), obtained by the standard alternating-sign transpose. |
| `func (o Operator) ApplyPoly(y Poly) RatFunc` | ApplyPoly applies the operator to a polynomial y. |
| `func (o Operator) ApplyRatFunc(y RatFunc) RatFunc` | ApplyRatFunc applies the operator to a rational function y, returning sum_i a_i y^(i). |
| `func (o Operator) Coeff(i int) RatFunc` | Coeff returns the coefficient of D^i, or zero when out of range. |
| `func (o Operator) Coeffs() []RatFunc` | Coeffs returns a copy of the coefficient slice in ascending order of D. |
| `func (o Operator) Equal(p Operator) bool` | Equal reports whether o and p are equal operators. |
| `func (o Operator) IndicialPolynomial() (Poly, bool)` | IndicialPolynomial returns the indicial polynomial at the ordinary point x=0 for an operator with polynomial coefficients, computed from the… |
| `func (o Operator) IsZero() bool` | IsZero reports whether o is the zero operator. |
| `func (o Operator) LeadingCoeff() RatFunc` | LeadingCoeff returns the coefficient of the highest-order term. |
| `func (o Operator) Mul(p Operator) Operator` | Mul returns the composition o∘p in the Ore ring, using D*a = a*D + a'. |
| `func (o Operator) Neg() Operator` | Neg returns -o. |
| `func (o Operator) Order() int` | Order returns the order (highest power of D) of the operator; the zero operator has order -1. |
| `func (o Operator) Pow(n int) Operator` | Pow raises the operator to the non-negative integer power n (composition). |
| `func (o Operator) ScalarMul(a RatFunc) Operator` | ScalarMul multiplies o on the left by the rational function a (a*o). |
| `func (o Operator) String() string` | String renders the operator in descending order of D. |
| `func (o Operator) Sub(p Operator) Operator` | Sub returns o-p. |
| `func (o Operator) SymbolPoly() (Poly, bool)` | SymbolPoly returns the symbol of the operator: the polynomial in a formal variable obtained by replacing D^i with the i-th power and freezing the… |

</details>

<details>
<summary><code>Poly</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func ConstPoly(r *big.Rat) Poly` | ConstPoly returns the constant polynomial equal to r. |
| `func ConstPolyInt(n int64) Poly` | ConstPolyInt returns the constant polynomial equal to the integer n. |
| `func DeterminantPoly(m [][]Poly) (Poly, error)` | DeterminantPoly returns the determinant of a square matrix of polynomials by evaluating over the field Q(x). |
| `func Monomial(coeff *big.Rat, deg int) Poly` | Monomial returns the polynomial coeff * x^deg. |
| `func NewPoly(coeffs ...*big.Rat) Poly` | NewPoly builds a polynomial from coefficients given in ascending degree order (coeffs[0] is the constant term). |
| `func OnePoly() Poly` | OnePoly returns the constant polynomial 1. |
| `func PolyFromInts(coeffs ...int64) Poly` | PolyFromInts builds a polynomial from integer coefficients in ascending degree order. |
| `func RischExpIntegrate(f, g Poly) (Poly, bool)` | RischExpIntegrate decides whether the integrand f(x) e^{g(x)} has an elementary antiderivative for polynomials f and g, and if so returns the… |
| `func RothsteinTragerResultant(f RatFunc) Poly` | RothsteinTragerResultant returns the Rothstein-Trager resultant R(z) = Res_x(D, C - z D') for the square-free proper integrand C/D, whose rational… |
| `func SolveRischDE(a, b Poly) (Poly, bool)` | SolveRischDE solves the Risch differential equation R' + a R = b for a polynomial R, given polynomials a and b. |
| `func WronskianPoly(ps []Poly) (Poly, error)` | WronskianPoly returns the Wronskian determinant of the polynomials ps. |
| `func XPoly() Poly` | XPoly returns the polynomial x. |
| `func ZeroPoly() Poly` | ZeroPoly returns the zero polynomial. |
| `func (p Poly) Add(q Poly) Poly` | Add returns p+q. |
| `func (p Poly) Clone() Poly` | Clone returns a deep copy of p. |
| `func (p Poly) Coeff(i int) *big.Rat` | Coeff returns the coefficient of x^i, or zero when i is out of range. |
| `func (p Poly) Coeffs() []*big.Rat` | Coeffs returns a fresh slice of the coefficients in ascending degree order. |
| `func (p Poly) ComplexRootsFloat(seed int64) []complex128` | ComplexRootsFloat returns the roots of p as complex numbers using the Durand-Kerner iteration seeded by the given value. |
| `func (p Poly) Compose(q Poly) Poly` | Compose returns p(q(x)). |
| `func (p Poly) ConstantTerm() *big.Rat` | ConstantTerm returns the coefficient of x^0. |
| `func (p Poly) Content() *big.Rat` | Content returns the (positive) rational content of p: the gcd of the integer numerators divided by the lcm of the denominators. |
| `func (p Poly) Degree() int` | Degree returns the degree of p. |
| `func (p Poly) Derivative() Poly` | Derivative returns the formal derivative dp/dx. |
| `func (p Poly) Discriminant() *big.Rat` | Discriminant returns the discriminant of p, defined via the resultant of p and its derivative. |
| `func (p Poly) DivMod(d Poly) (q, r Poly, err error)` | DivMod returns the quotient and remainder of Euclidean division p = q*d + r with deg r < deg d. |
| `func (p Poly) Equal(q Poly) bool` | Equal reports whether p and q are equal as polynomials. |
| `func (p Poly) EvalComplex(z complex128) complex128` | EvalComplex evaluates p at the complex point z. |
| `func (p Poly) EvalFloat(x float64) float64` | EvalFloat evaluates p at the floating-point point x. |
| `func (p Poly) EvalRat(x *big.Rat) *big.Rat` | EvalRat evaluates p at the rational point x using Horner's rule. |
| `func (p Poly) ExtendedGCD(q Poly) (g, s, t Poly)` | ExtendedGCD returns the monic gcd g of p and q together with cofactors s and t satisfying s*p + t*q = g. |
| `func (p Poly) GCD(q Poly) Poly` | GCD returns the monic greatest common divisor of p and q. |
| `func (p Poly) Integral() Poly` | Integral returns the antiderivative of p with zero constant of integration. |
| `func (p Poly) IsConstant() bool` | IsConstant reports whether p has degree at most zero. |
| `func (p Poly) IsMonic() bool` | IsMonic reports whether p is monic. |
| `func (p Poly) IsZero() bool` | IsZero reports whether p is the zero polynomial. |
| `func (p Poly) LeadingCoeff() *big.Rat` | LeadingCoeff returns the coefficient of the highest-degree term, or zero for the zero polynomial. |
| `func (p Poly) Monic() Poly` | Monic returns p scaled to be monic (leading coefficient 1). |
| `func (p Poly) Mul(q Poly) Poly` | Mul returns p*q. |
| `func (p Poly) Neg() Poly` | Neg returns -p. |
| `func (p Poly) PolyComplexRootClusters(seed int64, tol float64) []RootCluster` | PolyComplexRootClusters returns the clustered roots of p with inferred multiplicities, using tol to merge coincident roots. |
| `func (p Poly) Pow(n int) Poly` | Pow returns p raised to the non-negative integer power n. |
| `func (p Poly) PrimitivePart() Poly` | PrimitivePart returns p divided by its content, so the result has integer coprime coefficients and positive leading sign preserved. |
| `func (p Poly) Quo(d Poly) (Poly, error)` | Quo returns the quotient of p divided by d (see DivMod). |
| `func (p Poly) RationalRoots() []*big.Rat` | RationalRoots returns all distinct rational roots of p (without multiplicity) using the rational-root theorem. |
| `func (p Poly) Rem(d Poly) (Poly, error)` | Rem returns the remainder of p divided by d (see DivMod). |
| `func (p Poly) Resultant(q Poly) *big.Rat` | Resultant returns the resultant Res(p, q) computed with the Euclidean recurrence. |
| `func (p Poly) Reverse() Poly` | Reverse returns the reversal x^deg * p(1/x), i.e. |
| `func (p Poly) ScalarMul(r *big.Rat) Poly` | ScalarMul returns r*p for a rational scalar r. |
| `func (p Poly) Shift(a *big.Rat) Poly` | Shift returns p(x + a) for a rational shift a. |
| `func (p Poly) SquareFreeFactorization() []SquareFreeFactor` | SquareFreeFactorization returns the square-free factorisation of p using Yun's algorithm. |
| `func (p Poly) SquareFreePart() Poly` | SquareFreePart returns the square-free part of p, namely p / gcd(p, p'). |
| `func (p Poly) String() string` | String renders p in descending-degree human-readable form, e.g. |
| `func (p Poly) Sub(q Poly) Poly` | Sub returns p-q. |

</details>

<details>
<summary><code>RatFunc</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func ConstRatFunc(r *big.Rat) RatFunc` | ConstRatFunc returns the constant rational function equal to r. |
| `func DeterminantRatFunc(m [][]RatFunc) (RatFunc, error)` | DeterminantRatFunc returns the determinant of a square matrix of rational functions using fraction-free-free Gaussian elimination over the field Q(x). |
| `func LogarithmicDerivativeIsRational(f RatFunc) (RatFunc, bool)` | LogarithmicDerivativeIsRational reports whether the logarithmic derivative of a rational function f is again a proper rational function whose partial… |
| `func NewRatFunc(num, den Poly) (RatFunc, error)` | NewRatFunc builds the reduced rational function num/den. |
| `func OneRatFunc() RatFunc` | OneRatFunc returns the constant rational function 1. |
| `func RatFuncFromPoly(p Poly) RatFunc` | RatFuncFromPoly returns the rational function p/1. |
| `func ReduceToNormalForm(a2, a1, a0 Poly) (RatFunc, error)` | ReduceToNormalForm converts the second-order linear ODE a2 y” + a1 y' + a0 y = 0 into the reduced normal form z” = r z, returning r as a rational… |
| `func VariationOfParameters(ys []RatFunc, g RatFunc) (RatFunc, error)` | VariationOfParameters attempts to build a particular solution of the monic linear ODE L`y` = g from the fundamental system ys by integrating the… |
| `func VariationOfParametersIntegrands(ys []RatFunc, g RatFunc) ([]RatFunc, error)` | VariationOfParametersIntegrands returns the derivatives u_i'(x) of the variation-of-parameters coefficient functions for the monic linear ODE L`y` =… |
| `func WronskianRatFunc(fs []RatFunc) (RatFunc, error)` | WronskianRatFunc returns the Wronskian determinant of the rational functions fs. |
| `func XRatFunc() RatFunc` | XRatFunc returns the rational function x. |
| `func ZeroRatFunc() RatFunc` | ZeroRatFunc returns the zero rational function. |
| `func (f RatFunc) Add(g RatFunc) RatFunc` | Add returns f+g. |
| `func (f RatFunc) Den() Poly` | Den returns a copy of the denominator. |
| `func (f RatFunc) Derivative() RatFunc` | Derivative returns the derivative df/dx via the quotient rule. |
| `func (f RatFunc) Div(g RatFunc) (RatFunc, error)` | Div returns f/g. |
| `func (f RatFunc) Equal(g RatFunc) bool` | Equal reports whether f and g are equal as rational functions. |
| `func (f RatFunc) EvalComplex(z complex128) complex128` | EvalComplex evaluates f at the complex point z. |
| `func (f RatFunc) EvalFloat(x float64) float64` | EvalFloat evaluates f at the floating-point point x. |
| `func (f RatFunc) EvalRat(x *big.Rat) (*big.Rat, error)` | EvalRat evaluates f at the rational point x. |
| `func (f RatFunc) Inv() (RatFunc, error)` | Inv returns 1/f. |
| `func (f RatFunc) IsPolynomial() bool` | IsPolynomial reports whether the denominator is a constant. |
| `func (f RatFunc) IsProper() bool` | IsProper reports whether deg(num) < deg(den). |
| `func (f RatFunc) IsZero() bool` | IsZero reports whether f is the zero rational function. |
| `func (f RatFunc) LogDerivative() (RatFunc, error)` | LogDerivative returns the logarithmic derivative f'/f. |
| `func (f RatFunc) Mul(g RatFunc) RatFunc` | Mul returns f*g. |
| `func (f RatFunc) Neg() RatFunc` | Neg returns -f. |
| `func (f RatFunc) Num() Poly` | Num returns a copy of the numerator. |
| `func (f RatFunc) PartialFractions() (Poly, []PartialFractionTerm)` | PartialFractions returns the polynomial part together with the partial-fraction terms of f over Q, grouping the denominator by its square-free… |
| `func (f RatFunc) PolynomialPart() (quotient Poly, remainder RatFunc)` | PolynomialPart returns the polynomial quotient and proper remainder so that f = quotient + remainder with the remainder proper. |
| `func (f RatFunc) Pow(n int) (RatFunc, error)` | Pow raises f to the integer power n (n may be negative). |
| `func (f RatFunc) ScalarMul(r *big.Rat) RatFunc` | ScalarMul returns r*f for a rational scalar r. |
| `func (f RatFunc) String() string` | String renders f as "(num)/(den)" or just the numerator when polynomial. |
| `func (f RatFunc) Sub(g RatFunc) RatFunc` | Sub returns f-g. |

</details>

<details>
<summary><code>RationalIntegral</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func IntegrateRational(f RatFunc) (RationalIntegral, error)` | IntegrateRational returns the elementary integral of the rational function f, as a rational part plus a sum of logarithms with rational coefficients. |
| `func (r RationalIntegral) EvalFloat(x float64) float64` | EvalFloat evaluates the integral numerically at x, using log\|Arg\| for the logarithmic terms. |
| `func (r RationalIntegral) String() string` | String renders the integral as its rational part followed by its log terms. |

</details>

<details>
<summary><code>RealRationalIntegral</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func IntegrateRationalReal(f RatFunc) (RealRationalIntegral, error)` | IntegrateRationalReal returns the antiderivative of the rational function f over the reals, as a rational part plus real logarithms and arctangents. |
| `func (r RealRationalIntegral) EvalFloat(x float64) float64` | EvalFloat evaluates the resolved terms numerically at x, using log\|Arg\| for the logarithms. |
| `func (r RealRationalIntegral) String() string` | String renders the integral as the sum of its parts. |

</details>

<details>
<summary><code>RecurrenceSolution</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func SolveLinearRecurrence(coeffs []float64, seed int64, tol float64) (RecurrenceSolution, error)` | SolveLinearRecurrence builds the general solution of the homogeneous constant-coefficient linear recurrence sum_i coeffs`i` a_{n+i} = 0, where… |
| `func SolveRecurrenceIVP(coeffs []float64, initial []float64, seed int64, tol float64) ([]float64, RecurrenceSolution, error)` | SolveRecurrenceIVP fits the constants of the general solution to the initial data initial`i` = a_i for i = 0..order-1. |
| `func (s RecurrenceSolution) Basis() []RecurrenceTerm` | Basis returns the fundamental basis terms. |
| `func (s RecurrenceSolution) Dimension() int` | Dimension returns the number of basis terms (the order of the recurrence). |
| `func (s RecurrenceSolution) EvalBasis(n int) []float64` | EvalBasis returns the values of every basis term at index n. |
| `func (s RecurrenceSolution) Evaluate(consts []float64, n int) float64` | Evaluate returns sum_i consts`i` * basis_i(n). |
| `func (s RecurrenceSolution) String() string` | String renders the general solution as a linear combination. |

</details>

<details>
<summary><code>RecurrenceTerm</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (t RecurrenceTerm) Eval(n int) float64` | Eval returns the value of the basis term at integer index n. |
| `func (t RecurrenceTerm) String() string` | String renders the recurrence basis term. |

</details>

### Variables

`ErrDivByZero`, `ErrDim`, `ErrNotSquare`, `ErrSingular`, `ErrEmpty`, `ErrDegree`, `ErrNoSolution`, `ErrNonRational`, `ErrConverge`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/diffalgebra).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
