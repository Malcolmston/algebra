# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.0] - 2026-07-12
### Added
- Initial release — a standard-library-only Go port of Python SymPy: a small
  computer-algebra system for symbolic mathematics.
- Immutable, value-based `Expr` expression trees (symbols, arbitrary-precision
  integers and rationals, floats, named constants `Pi`/`E`, sums, products,
  powers, and the elementary functions `sin`, `cos`, `tan`, `exp`, `log`,
  `sqrt`) with canonicalizing `Add`/`Mul`/`Pow` constructors, structural
  equality and a deterministic canonical ordering.
- Infix parser (`Parse` / `MustParse`) with unary signs, parentheses, integer
  and decimal literals, the constants `pi` and `E`, function calls and implicit
  multiplication.
- Symbolic operations: `Simplify`, `Expand`, `Diff` (differentiation), `Subs`
  (substitution), `Integrate` (a documented subset), `Solve` (linear and
  quadratic equations), and the `Factor` and `Collect` polynomial helpers.
- Numeric evaluation (`Eval` / `Evalf`) and readable infix pretty-printing.
- CI: gofmt · vet · build gate a `-race` + coverage run, plus golangci-lint v2,
  govulncheck, cross-compile, dependency review, and VERSION-driven releases.

[Unreleased]: https://github.com/malcolmston/algebra/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/malcolmston/algebra/releases/tag/v0.1.0
