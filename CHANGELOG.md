# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Added
- Trigonometry completion: `Sec`, `Csc`, `Cot`, the inverse functions `Asin`,
  `Acos`, `Atan`, `Acot`, `Asec`, `Acsc` and the two-argument `Atan2`, with
  exact values at the standard rational multiples of Pi.
- Hyperbolic functions `Sinh`, `Cosh`, `Tanh`, `Coth`, `Sech`, `Csch` and the
  inverses `Asinh`, `Acosh`, `Atanh`.
- Special functions `Abs`, `Sign`, `Floor`, `Ceil`, `Factorial`, `Gamma`,
  `Beta`, `Erf`, `Erfc` (with a numeric `digamma`).
- Complex numbers: the imaginary unit `I` (with `I^2 -> -1` and Euler folds such
  as `exp(I*Pi) -> -1`), plus `Conjugate`, `Re`, `Im`, `Abs` (modulus), `Arg`
  and a complex128 evaluator `Evalc`/`EvalComplex`.
- Differentiation of all the above (chain rule throughout), and simplification
  identities: `sin^2+cos^2 -> 1`, double-angle, `sin/cos -> tan`, `log(exp x)`
  and `exp(log x)` inverses.
- Integration extensions: `atan`/`asin`/`asinh`/`acosh` forms, the remaining
  trig/hyperbolic antiderivatives, integration by parts for polynomial×exp/sin/
  cos, and rational functions via partial fractions.
- Calculus additions `Limit` (with L'Hôpital and limits at infinity), `Series`
  (Taylor/Maclaurin), `Summation` and `Product` (arithmetic/geometric/finite).
- Solver extensions: cubic and quartic equations (exact rational/quadratic
  factors, complex conjugate roots, numeric fallback for irreducible factors)
  and `SolveSystem` for linear systems via Gaussian elimination.
- Parser support for all new function names, the `I`/`oo` constants and postfix
  factorial (`n!`).

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
