// Package algebra is a small computer-algebra system (CAS) written entirely
// with the Go standard library. It represents mathematical expressions as
// immutable trees, manipulates them symbolically, and prints them back as
// readable infix strings. It is a from-scratch, spiritual port of a tiny
// subset of Python's SymPy.
//
// # Expression model
//
// Every expression implements the [Expr] interface. The concrete node kinds
// are:
//
//   - [Symbol]   — a named variable such as x.
//   - [Integer]  — an arbitrary-precision integer (backed by math/big).
//   - [Rational] — an exact fraction with denominator > 1 (backed by math/big).
//   - [Float]    — an inexact float64 literal.
//   - [Constant] — a named mathematical constant; [Pi], [E], the imaginary
//     unit [I] and the infinities [Inf]/[NegInf] are provided.
//   - internal Add, Mul and Pow nodes for sums, products and powers.
//   - internal function nodes for the elementary, trigonometric, hyperbolic and
//     special functions listed below.
//
// Expressions are value-based and immutable: constructors never mutate their
// arguments and the big.Int/big.Rat payloads are never modified after
// construction. Structural equality is available through [Expr.Equal].
//
// # Building expressions
//
// Expressions can be built with the constructors [Sym], [Int], [Rat], [Flt],
// [Add], [Mul], [Pow], [Sin], [Cos], [Tan], [Exp], [Log] and [Sqrt], with the
// fluent builder methods [Expr.Add], [Expr.Mul] and [Expr.Pow], or by parsing
// a string with [Parse]:
//
//	x := Sym("x")
//	e := x.Pow(Int(2)).Add(x.Mul(Int(2)), Int(1)) // x^2 + 2*x + 1
//	e, _ = Parse("x^2 + 2*x + 1")                  // same thing
//
// The [Add], [Mul] and [Pow] constructors are canonicalizing: they flatten
// nested sums/products, fold numeric constants, apply the algebraic identities
// (x+0, x*1, x*0, x^1, x^0), combine like terms, combine repeated powers and
// sort arguments into a deterministic canonical order. As a result equal
// expressions compare equal with [Expr.Equal] and print identically.
//
// # Operations
//
// The package provides [Simplify], [Expand], [Diff] (symbolic
// differentiation), [Subs] (substitution), [Eval] and [Evalf] (numeric
// evaluation), [Integrate] (symbolic integration of a documented subset),
// [Solve] (polynomial equations) and the bonus [Factor] and [Collect] helpers
// for univariate polynomials.
//
// # Functions
//
// Beyond the original sin/cos/tan/exp/log/sqrt, the package now provides the
// reciprocal and inverse trigonometric functions ([Sec], [Csc], [Cot], [Asin],
// [Acos], [Atan], [Acot], [Asec], [Acsc], [Atan2]); the hyperbolic functions
// and their inverses ([Sinh], [Cosh], [Tanh], [Coth], [Sech], [Csch], [Asinh],
// [Acosh], [Atanh]); and the special functions [Abs], [Sign], [Floor], [Ceil],
// [Factorial], [Gamma], [Beta], [Erf] and [Erfc]. [Simplify] returns exact
// values at the standard rational multiples of Pi and applies the Pythagorean
// and double-angle identities.
//
// # Complex numbers
//
// The imaginary unit [I] participates in ordinary arithmetic (I^2 folds to -1
// and clean Euler identities such as exp(I*Pi) -> -1 are recognised).
// [Conjugate], [Re], [Im], [Abs] (modulus) and [Arg] decompose a complex value,
// and [Evalc] evaluates any expression to a complex128.
//
// # Calculus and solving
//
// [Limit] computes limits (including L'Hôpital's rule and limits at infinity),
// [Series] produces Taylor/Maclaurin expansions, and [Summation] and [Product]
// find closed forms for polynomial, geometric and finite sums and products.
// [Integrate] additionally handles the arctangent and arcsine forms, the
// remaining trigonometric and hyperbolic antiderivatives, integration by parts
// for polynomial×exp/sin/cos, and rational functions via partial fractions.
// [Solve] handles polynomials of any degree — exact rational, quadratic (with
// complex conjugate roots) and higher-degree factors, with a numeric fallback —
// and [SolveSystem] solves linear systems by Gaussian elimination.
//
// Correctness is favoured over coverage: an operation that does not know how
// to transform an expression returns it unchanged (or an explicit unevaluated
// node) rather than producing a wrong answer.
package algebra
