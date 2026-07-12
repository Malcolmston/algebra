# algebra

[![Go Test](https://github.com/Malcolmston/algebra/actions/workflows/go-test.yml/badge.svg)](https://github.com/Malcolmston/algebra/actions/workflows/go-test.yml)
[![Go Lint](https://github.com/Malcolmston/algebra/actions/workflows/go-lint.yml/badge.svg)](https://github.com/Malcolmston/algebra/actions/workflows/go-lint.yml)
[![Go Vuln](https://github.com/Malcolmston/algebra/actions/workflows/go-vuln.yml/badge.svg)](https://github.com/Malcolmston/algebra/actions/workflows/go-vuln.yml)
[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra.svg)](https://pkg.go.dev/github.com/malcolmston/algebra)
[![Go Report Card](https://goreportcard.com/badge/github.com/malcolmston/algebra)](https://goreportcard.com/report/github.com/malcolmston/algebra)
[![Go Version](https://img.shields.io/github/go-mod/go-version/malcolmston/algebra)](go.mod)
[![Release](https://img.shields.io/github/v/release/malcolmston/algebra?sort=semver)](https://github.com/malcolmston/algebra/releases)
[![Last Commit](https://img.shields.io/github/last-commit/malcolmston/algebra)](https://github.com/malcolmston/algebra/commits)
[![Code Size](https://img.shields.io/github/languages/code-size/malcolmston/algebra)](https://github.com/malcolmston/algebra)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg)](CONTRIBUTING.md)
[![Docs](https://img.shields.io/badge/docs-pages-2f9bff)](https://malcolmston.github.io/algebra/)

A standard-library-only Go port of SymPy — symbolic mathematics: build,
simplify, differentiate, integrate, and solve algebraic expressions.

## What it is

`algebra` is a small computer-algebra system (CAS) written entirely with the Go
standard library. It represents mathematical expressions as immutable trees,
manipulates them symbolically, and prints them back as readable infix strings.
It is a from-scratch, spiritual port of a tiny subset of Python's
[SymPy](https://github.com/sympy/sympy).

Every expression is a value implementing the `Expr` interface. The concrete
node kinds are symbols, arbitrary-precision integers and rationals (backed by
`math/big`), inexact floats, named constants (`Pi`, `E`), sums, products,
powers, and the elementary functions `sin`, `cos`, `tan`, `exp`, `log` and
`sqrt`. Expressions are immutable and value-based: the `Add`, `Mul` and `Pow`
constructors canonicalize as they build, so mathematically equal expressions
compare equal with `Expr.Equal` and print identically.

## Installation

```sh
go get github.com/malcolmston/algebra
```

Requires Go 1.24 or newer. There are no third-party dependencies.

## Quick start

Parse an expression, differentiate it, and simplify the result:

```go
package main

import (
	"fmt"

	"github.com/malcolmston/algebra"
)

func main() {
	// f(x) = x^2 + 2x + 1
	e := algebra.MustParse("x^2 + 2*x + 1")

	// f'(x) = 2x + 2
	d := algebra.Diff(e, algebra.Sym("x"))
	fmt.Println(algebra.Simplify(d)) // 2*x + 2
}
```

Solve a quadratic equation (`Solve` returns the real roots of `expr == 0`):

```go
package main

import (
	"fmt"

	"github.com/malcolmston/algebra"
)

func main() {
	x := algebra.Sym("x")

	// x^2 - 5x + 6 = 0  ->  x = 2, x = 3
	roots, err := algebra.Solve(algebra.MustParse("x^2 - 5*x + 6"), x)
	if err != nil {
		panic(err)
	}
	for _, r := range roots {
		fmt.Println(r) // 2, then 3
	}
}
```

See [`examples/main.go`](examples/main.go) for a full tour: parsing,
differentiation, simplification, expansion, factoring, substitution, numeric
evaluation, integration and solving.

## Features

- **Expression trees** — immutable, value-based `Expr` nodes with structural
  equality and a deterministic canonical ordering.
- **Parser** — `Parse` / `MustParse` read ordinary infix notation, including
  `+ - * / ^`, unary signs, parentheses, decimal and integer literals, the
  constants `pi` and `E`, function calls (`sin`, `cos`, `tan`, `exp`, `log`,
  `ln`, `sqrt`), and implicit multiplication (`2x`, `2(x+1)`, `x y`).
- **Simplify & expand** — `Simplify` folds numeric arithmetic, applies the
  algebraic identities and combines like terms; `Expand` distributes products
  over sums and expands integer powers (binomial/multinomial).
- **Differentiation** — `Diff` implements the sum, product, power, quotient and
  chain rules plus the derivatives of the elementary functions.
- **Integration** — `Integrate` returns a symbolic antiderivative for a
  documented subset, or an explicit unevaluated `Integral` node otherwise.
- **Equation solving** — `Solve` returns the exact real roots of linear and
  quadratic polynomials; `Factor` and `Collect` are univariate-polynomial
  helpers.
- **Substitution** — `Subs` replaces a symbol with any expression and rebuilds
  through the canonicalizing constructors.
- **Numeric evaluation** — `Eval` (with a symbol environment) and `Evalf`
  (fully numeric) reduce an expression to a `float64`.
- **Pretty-printing** — `String` renders each node as a readable infix string
  with correct precedence, parenthesization and textbook term ordering.
- **Zero dependencies** — pure Go standard library; nothing to audit but the
  toolchain.

## Coverage & limits

Correctness is favoured over coverage: an operation that does not know how to
transform an expression returns it unchanged rather than producing a wrong
answer. Concretely, as reported by the core:

- **Integration** covers constants; powers `x^n` for any integer `n` (including
  `1/x -> log x`); sums term by term; constant multiples; the elementary
  functions `exp`, `sin` and `cos`; and — via a linear substitution — any of
  the above whose argument or base is linear in the variable (`a*v + b`).
  Anything else comes back as an unevaluated `Integral` node.
- **Solving** covers polynomials of degree 1 (linear) and 2 (quadratic), with
  exact roots. Higher-degree and complex roots are deferred.

## Documentation

- Full API reference on [pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra).
- Docs site: <https://malcolmston.github.io/algebra/>.

## License

MIT
