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
//   - [Constant] — a named mathematical constant; [Pi] and [E] are provided.
//   - internal Add, Mul and Pow nodes for sums, products and powers.
//   - internal function nodes for sin, cos, tan, exp, log and sqrt.
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
// [Solve] (linear and quadratic equations), and the bonus [Factor] and
// [Collect] helpers for univariate polynomials.
//
// Correctness is favoured over coverage: an operation that does not know how
// to transform an expression returns it unchanged (or, for [Integrate], an
// explicit unevaluated Integral node) rather than producing a wrong answer.
package algebra
