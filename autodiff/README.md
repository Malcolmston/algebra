# autodiff — implements automatic differentiation in pure Go using only the standard library

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/autodiff.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/autodiff)

Package autodiff implements automatic differentiation in pure Go using only
the standard library.

- Forward mode via `Dual` numbers. A dual number carries a value together with
its derivative along a chosen direction. Evaluating a function on a seeded
dual number yields the value and the directional derivative in a single pass.
Forward mode is efficient when the number of inputs is small relative to the
number of outputs and underlies `Derivative`, `Gradient`, `Jacobian`,
`DirectionalDerivative`, `PartialDerivative` and `Divergence`.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/autodiff"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func Derivative(f func(Dual) Dual, x float64) float64` | Derivative returns f'(x), the first derivative of the scalar function f at the point x, obtained by seeding a unit derivative and reading the… |
| `func Derivatives2(f func(HyperDual) HyperDual, x float64) (value, first, second float64)` | Derivatives2 returns f(x), f'(x) and f”(x) together from one hyper-dual evaluation of the scalar function f at x. |
| `func DirectionalDerivative(f func([]Dual) Dual, x, dir []float64) float64` | DirectionalDerivative returns the derivative of the scalar field f: Rⁿ → R at x along the direction dir, that is ∇f(x)·dir. |
| `func Divergence(f func([]Dual) []Dual, x []float64) float64` | Divergence returns ∇·F, the divergence of the vector field F: Rⁿ → Rⁿ at x, that is the trace Σᵢ ∂Fᵢ/∂xᵢ of its Jacobian. |
| `func Gradient(f func([]Dual) Dual, x []float64) []float64` | Gradient returns the gradient ∇f of the scalar field f: Rⁿ → R evaluated at x. |
| `func GradientHessian(f func([]HyperDual) HyperDual, x []float64) (grad []float64, hess [][]float64)` | GradientHessian returns both the gradient and the Hessian of the scalar field f: Rⁿ → R at x from the same set of hyper-dual passes. |
| `func GradientReverse(f func(t *Tape, x []Var) Var, x []float64) []float64` | GradientReverse returns the gradient of the scalar function f: Rⁿ → R at x using a single reverse-mode sweep. |
| `func Hessian(f func([]HyperDual) HyperDual, x []float64) [][]float64` | Hessian returns the Hessian matrix of second partial derivatives of the scalar field f: Rⁿ → R at x. |
| `func HessianDiagonal(f func([]HyperDual) HyperDual, x []float64) []float64` | HessianDiagonal returns the diagonal of the Hessian of the scalar field f: Rⁿ → R at x, that is the pure second partial derivatives… |
| `func HessianVectorProduct(f func([]HyperDual) HyperDual, x, v []float64) []float64` | HessianVectorProduct returns H(x)·v, the product of the Hessian of the scalar field f: Rⁿ → R at x with the vector v, without forming the full… |
| `func Jacobian(f func([]Dual) []Dual, x []float64) [][]float64` | Jacobian returns the Jacobian matrix of the vector field f: Rⁿ → Rᵐ at x. |
| `func JacobianVectorProduct(f func([]Dual) []Dual, x, v []float64) []float64` | JacobianVectorProduct returns J(x)·v, the product of the Jacobian of f: Rⁿ → Rᵐ at x with the vector v, without forming J. It seeds the inputs… |
| `func Laplacian(f func([]HyperDual) HyperDual, x []float64) float64` | Laplacian returns ∇²f, the Laplacian of the scalar field f: Rⁿ → R at x, i.e. |
| `func PartialDerivative(f func([]Dual) Dual, x []float64, i int) float64` | PartialDerivative returns ∂f/∂xᵢ, the partial derivative of the scalar field f: Rⁿ → R with respect to its i-th argument, evaluated at x. |
| `func SecondDerivative(f func(HyperDual) HyperDual, x float64) float64` | SecondDerivative returns f”(x), the second derivative of the scalar function f at x, from a single hyper-dual evaluation with both first-order… |
| `func ValueAndDerivative(f func(Dual) Dual, x float64) (value, derivative float64)` | ValueAndDerivative returns both f(x) and f'(x) from a single forward-mode evaluation of f at x. |
| `func ValueAndGradient(f func([]Dual) Dual, x []float64) (value float64, grad []float64)` | ValueAndGradient returns f(x) together with its gradient ∇f(x). |
| `func ValueAndGradientReverse(f func(t *Tape, x []Var) Var, x []float64) (value float64, grad []float64)` | ValueAndGradientReverse returns both f(x) and its gradient ∇f(x) from a single reverse-mode sweep. |

### Types

| Type | What it is |
| --- | --- |
| `Dual` | Dual is a forward-mode dual number of the form Val + Der*ε, where ε is a nilpotent infinitesimal satisfying ε² = 0. |
| `HyperDual` | HyperDual is a second-order forward-mode number of the form |
| `Tape` | Tape records the elementary operations of a computation in evaluation order so that adjoints can be propagated backward from an output to every… |
| `Var` | Var is a handle to a single node on a `Tape`. |

<details>
<summary><code>Dual</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Abs(x Dual) Dual` | Abs returns \|x\|; its derivative is the sign of x and is undefined at zero, where this implementation returns a zero derivative. |
| `func Acos(x Dual) Dual` | Acos returns arccos(x); its derivative is -1/√(1-x²). |
| `func Acosh(x Dual) Dual` | Acosh returns the inverse hyperbolic cosine arcosh(x); its derivative is 1/√(x²-1) and is defined for x > 1. |
| `func Asin(x Dual) Dual` | Asin returns arcsin(x); its derivative is 1/√(1-x²). |
| `func Asinh(x Dual) Dual` | Asinh returns the inverse hyperbolic sine arsinh(x); its derivative is 1/√(x²+1). |
| `func Atan(x Dual) Dual` | Atan returns arctan(x); its derivative is 1/(1+x²). |
| `func Atan2(y, x Dual) Dual` | Atan2 returns the two-argument arctangent atan2(y, x). |
| `func Atanh(x Dual) Dual` | Atanh returns the inverse hyperbolic tangent artanh(x); its derivative is 1/(1-x²) and is defined for \|x\| < 1. |
| `func Cbrt(x Dual) Dual` | Cbrt returns the cube root ∛x; its derivative is 1/(3·x^(2/3)). |
| `func Constant(val float64) Dual` | Constant returns the dual number val + 0*ε, representing a quantity whose derivative with respect to the independent variable is zero. |
| `func Cos(x Dual) Dual` | Cos returns cos(x); its derivative is -sin(x). |
| `func Cosh(x Dual) Dual` | Cosh returns cosh(x); its derivative is sinh(x). |
| `func Cot(x Dual) Dual` | Cot returns the cotangent cot(x) = cos(x)/sin(x); its derivative is -csc²(x) = -1/sin²(x). |
| `func Coth(x Dual) Dual` | Coth returns the hyperbolic cotangent coth(x) = cosh(x)/sinh(x); its derivative is -csch²(x) = -1/sinh²(x). |
| `func Csc(x Dual) Dual` | Csc returns the cosecant csc(x) = 1/sin(x); its derivative is -csc(x)·cot(x). |
| `func Csch(x Dual) Dual` | Csch returns the hyperbolic cosecant csch(x) = 1/sinh(x); its derivative is -csch(x)·coth(x). |
| `func Erf(x Dual) Dual` | Erf returns the error function erf(x); its derivative is (2/√π)·e^(-x²). |
| `func Erfc(x Dual) Dual` | Erfc returns the complementary error function erfc(x) = 1-erf(x); its derivative is -(2/√π)·e^(-x²). |
| `func Exp(x Dual) Dual` | Exp returns e^x; its derivative is e^x. |
| `func Exp2(x Dual) Dual` | Exp2 returns 2^x; its derivative is ln(2)·2^x. |
| `func Expm1(x Dual) Dual` | Expm1 returns e^x - 1 accurately for small x; its derivative is e^x. |
| `func Hypot(x, y Dual) Dual` | Hypot returns √(x²+y²) computed without spurious overflow; its differential is (x·dx + y·dy)/hypot(x,y). |
| `func Log(x Dual) Dual` | Log returns the natural logarithm ln(x); its derivative is 1/x. |
| `func Log10(x Dual) Dual` | Log10 returns the base-10 logarithm log₁₀(x); its derivative is 1/(x·ln10). |
| `func Log1p(x Dual) Dual` | Log1p returns ln(1+x) accurately for small x; its derivative is 1/(1+x). |
| `func Log2(x Dual) Dual` | Log2 returns the base-2 logarithm log₂(x); its derivative is 1/(x·ln2). |
| `func Max(x, y Dual) Dual` | Max returns the larger of x and y as a dual number, propagating the derivative of the selected argument. |
| `func Min(x, y Dual) Dual` | Min returns the smaller of x and y as a dual number, propagating the derivative of the selected argument. |
| `func NewDual(val, der float64) Dual` | NewDual returns the dual number val + der*ε. |
| `func Pow(x, y Dual) Dual` | Pow returns x^y where both the base and exponent are dual numbers. |
| `func PowReal(x Dual, p float64) Dual` | PowReal returns x^p for a constant real exponent p; its derivative is p·x^(p-1). |
| `func Powf(a float64, x Dual) Dual` | Powf returns a^x for a constant real base a and dual exponent x. |
| `func Relu(x Dual) Dual` | Relu returns the rectified linear unit max(0, x); its derivative is 1 for positive x and 0 otherwise, with the sub-gradient at zero taken as 0. |
| `func Sec(x Dual) Dual` | Sec returns the secant sec(x) = 1/cos(x); its derivative is sec(x)·tan(x). |
| `func Sech(x Dual) Dual` | Sech returns the hyperbolic secant sech(x) = 1/cosh(x); its derivative is -sech(x)·tanh(x). |
| `func Sigmoid(x Dual) Dual` | Sigmoid returns the logistic function σ(x) = 1/(1+e^(-x)); its derivative is σ(x)·(1-σ(x)). |
| `func Sin(x Dual) Dual` | Sin returns sin(x); its derivative is cos(x). |
| `func Sinh(x Dual) Dual` | Sinh returns sinh(x); its derivative is cosh(x). |
| `func Softplus(x Dual) Dual` | Softplus returns the smooth rectifier ln(1+e^x); its derivative is the logistic function σ(x). |
| `func Sqrt(x Dual) Dual` | Sqrt returns the square root √x; its derivative is 1/(2√x). |
| `func Tan(x Dual) Dual` | Tan returns tan(x); its derivative is sec²(x) = 1/cos²(x). |
| `func Tanh(x Dual) Dual` | Tanh returns tanh(x); its derivative is sech²(x) = 1-tanh²(x). |
| `func Variable(val float64) Dual` | Variable returns the dual number val + 1*ε, the seed for the independent variable of a forward-mode differentiation. |
| `func (x Dual) Add(y Dual) Dual` | Add returns the sum x + y. |
| `func (x Dual) AddReal(k float64) Dual` | AddReal returns x + k where k is a real constant, leaving the derivative unchanged. |
| `func (x Dual) Div(y Dual) Dual` | Div returns the quotient x / y using the quotient rule (a/b)' = (a'*b - a*b')/b². |
| `func (x Dual) Inv() Dual` | Inv returns the reciprocal 1/x, whose derivative is -x'/x². |
| `func (x Dual) Mul(y Dual) Dual` | Mul returns the product x * y using the product rule (a*b)' = a'*b + a*b'. |
| `func (x Dual) Neg() Dual` | Neg returns the negation -x. |
| `func (x Dual) Real() float64` | Real returns the primal (real) part of the dual number, discarding the derivative. |
| `func (x Dual) Scale(k float64) Dual` | Scale returns x multiplied by the real scalar k, differentiating k as a constant so that (k*x)' = k*x'. |
| `func (x Dual) String() string` | String formats the dual number in the canonical "a+bε" notation. |
| `func (x Dual) Sub(y Dual) Dual` | Sub returns the difference x - y. |

</details>

<details>
<summary><code>HyperDual</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func HyperConstant(val float64) HyperDual` | HyperConstant returns the hyper-dual number val with all infinitesimal parts zero, representing a quantity with vanishing first and second… |
| `func HyperVariable(val float64) HyperDual` | HyperVariable returns the hyper-dual seed val + ε₁ + ε₂ for a single independent variable. |
| `func NewHyperDual(val, e1, e2, e12 float64) HyperDual` | NewHyperDual returns the hyper-dual number val + e1·ε₁ + e2·ε₂ + e12·ε₁ε₂. |
| `func (x HyperDual) Abs() HyperDual` | Abs returns \|x\|; the first derivative is the sign of x and the second is zero away from the origin. |
| `func (x HyperDual) Acos() HyperDual` | Acos returns arccos(x) for \|x\| < 1; it uses (acos)' = -1/√(1-x²) and (acos)” = -x/(1-x²)^{3/2}. |
| `func (x HyperDual) Add(y HyperDual) HyperDual` | Add returns the sum x + y. |
| `func (x HyperDual) AddReal(k float64) HyperDual` | AddReal returns x + k where k is a real constant. |
| `func (x HyperDual) Asin() HyperDual` | Asin returns arcsin(x) for \|x\| < 1; it uses (asin)' = 1/√(1-x²) and (asin)” = x/(1-x²)^{3/2}. |
| `func (x HyperDual) Atan() HyperDual` | Atan returns arctan(x); it uses (atan)' = 1/(1+x²) and (atan)” = -2x/(1+x²)². |
| `func (x HyperDual) Cos() HyperDual` | Cos returns cos(x); it uses cos' = -sin and cos” = -cos. |
| `func (x HyperDual) Cosh() HyperDual` | Cosh returns cosh(x); it uses cosh' = sinh and cosh” = cosh. |
| `func (x HyperDual) Div(y HyperDual) HyperDual` | Div returns the quotient x / y as x · y⁻¹. |
| `func (x HyperDual) Erf() HyperDual` | Erf returns the error function erf(x); it uses erf' = (2/√π)·e^(-x²) and erf” = -2x·erf'. |
| `func (x HyperDual) Exp() HyperDual` | Exp returns e^x; every derivative of the exponential is itself. |
| `func (x HyperDual) Inv() HyperDual` | Inv returns the reciprocal 1/x, obtained by applying the reciprocal function and its derivatives r(a)=1/a, r'(a)=-1/a², r”(a)=2/a³. |
| `func (x HyperDual) Log() HyperDual` | Log returns ln(x); it uses (ln)' = 1/x and (ln)” = -1/x². |
| `func (x HyperDual) Mul(y HyperDual) HyperDual` | Mul returns the product x · y using the full second-order product rule, including the cross term E1·E2 + E2·E1 that feeds the mixed slot. |
| `func (x HyperDual) Neg() HyperDual` | Neg returns the negation -x. |
| `func (x HyperDual) PowReal(p float64) HyperDual` | PowReal returns x^p for a constant real exponent p, with first derivative p·x^(p-1) and second derivative p·(p-1)·x^(p-2). |
| `func (x HyperDual) Real() float64` | Real returns the primal part of the hyper-dual number. |
| `func (x HyperDual) Scale(k float64) HyperDual` | Scale returns x multiplied by the real scalar k. |
| `func (x HyperDual) Sigmoid() HyperDual` | Sigmoid returns the logistic function σ(x) = 1/(1+e^(-x)); it uses σ' = σ(1-σ) and σ” = σ(1-σ)(1-2σ). |
| `func (x HyperDual) Sin() HyperDual` | Sin returns sin(x); it uses sin' = cos and sin” = -sin. |
| `func (x HyperDual) Sinh() HyperDual` | Sinh returns sinh(x); it uses sinh' = cosh and sinh” = sinh. |
| `func (x HyperDual) Sqrt() HyperDual` | Sqrt returns √x; it uses derivatives 1/(2√x) and -1/(4·x^(3/2)). |
| `func (x HyperDual) Sub(y HyperDual) HyperDual` | Sub returns the difference x - y. |
| `func (x HyperDual) Tan() HyperDual` | Tan returns tan(x); it uses tan' = sec² and tan” = 2·sec²·tan. |
| `func (x HyperDual) Tanh() HyperDual` | Tanh returns tanh(x); it uses tanh' = 1-tanh² and tanh” = -2·tanh·(1-tanh²). |

</details>

<details>
<summary><code>Tape</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewTape() *Tape` | NewTape returns a new, empty operation tape ready to record a computation. |
| `func (t *Tape) Backward(y Var) []float64` | Backward propagates adjoints from the output node y backward through the tape and returns the gradient of y with respect to every registered… |
| `func (t *Tape) Constant(v float64) Var` | Constant records a constant with value v on the tape. |
| `func (t *Tape) NumNodes() int` | NumNodes returns the total number of nodes (constants, variables and intermediate results) currently recorded on the tape. |
| `func (t *Tape) NumVars() int` | NumVars returns the number of independent variables registered on the tape. |
| `func (t *Tape) Reset()` | Reset discards every recorded node so the tape can be reused for a fresh computation. |
| `func (t *Tape) Variable(v float64) Var` | Variable records an independent variable with value v on the tape and registers it so that its adjoint appears in the gradient returned by… |

</details>

<details>
<summary><code>Var</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (v Var) Abs() Var` | Abs records and returns \|v\| with derivative sign(v); the derivative at zero is taken as zero. |
| `func (v Var) Add(w Var) Var` | Add records and returns the sum v + w. |
| `func (v Var) AddReal(k float64) Var` | AddReal records and returns v + k for a real constant k. |
| `func (v Var) Cos() Var` | Cos records and returns cos(v) with derivative -sin(v). |
| `func (v Var) Cosh() Var` | Cosh records and returns cosh(v) with derivative sinh(v). |
| `func (v Var) Div(w Var) Var` | Div records and returns the quotient v / w with partials 1/w and -v/w². |
| `func (v Var) Erf() Var` | Erf records and returns the error function erf(v) with derivative (2/√π)·e^(-v²). |
| `func (v Var) Exp() Var` | Exp records and returns e^v. |
| `func (v Var) Log() Var` | Log records and returns ln(v) with derivative 1/v. |
| `func (v Var) Mul(w Var) Var` | Mul records and returns the product v · w, whose local partials are the opposite operand's value. |
| `func (v Var) Neg() Var` | Neg records and returns the negation -v. |
| `func (v Var) Pow(w Var) Var` | Pow records and returns v^w for two tape variables, using v^w = exp(w·ln v) so that the partials are w·v^(w-1) and v^w·ln v. |
| `func (v Var) PowReal(p float64) Var` | PowReal records and returns v^p for a constant real exponent p, with derivative p·v^(p-1). |
| `func (v Var) Recip() Var` | Recip records and returns the reciprocal 1/v with derivative -1/v². |
| `func (v Var) Relu() Var` | Relu records and returns the rectified linear unit max(0, v); the derivative is 1 for positive v and 0 otherwise, with the sub-gradient at zero taken… |
| `func (v Var) Scale(k float64) Var` | Scale records and returns k·v for a real constant k. |
| `func (v Var) Sigmoid() Var` | Sigmoid records and returns the logistic function σ(v) with derivative σ(v)·(1-σ(v)). |
| `func (v Var) Sin() Var` | Sin records and returns sin(v) with derivative cos(v). |
| `func (v Var) Sinh() Var` | Sinh records and returns sinh(v) with derivative cosh(v). |
| `func (v Var) Softplus() Var` | Softplus records and returns the smooth rectifier ln(1+e^v), evaluated stably; its derivative is the logistic function σ(v). |
| `func (v Var) Sqrt() Var` | Sqrt records and returns √v with derivative 1/(2√v). |
| `func (v Var) Sub(w Var) Var` | Sub records and returns the difference v - w. |
| `func (v Var) Tan() Var` | Tan records and returns tan(v) with derivative sec²(v). |
| `func (v Var) Tanh() Var` | Tanh records and returns tanh(v) with derivative 1-tanh²(v). |
| `func (v Var) Value() float64` | Value returns the primal value stored at this node. |

</details>

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/autodiff).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
