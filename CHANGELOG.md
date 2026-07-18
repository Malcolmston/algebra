# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.6.0] - 2026-07-18
### Added
- Seven new sub-packages, each covering a numerical-mathematics domain with the
  Go standard library only:
  - **special** — Bessel (J/Y/I/K, spherical, Kelvin, Struve), Airy, elliptic
    integrals, the error-function and Fresnel family, exponential/sine/cosine
    integrals, the zeta family, Lambert W, polylogarithms and the gamma family.
  - **orthopoly** — the classical orthogonal polynomials (Legendre, Chebyshev,
    Hermite, Laguerre, Gegenbauer, Jacobi) together with their derivatives,
    norms and the associated Gauss quadrature rules.
  - **combin** — combinatorial counting (factorials, binomials, Stirling,
    Bell, Catalan, partitions) plus enumerators for permutations, combinations
    and related sequences.
  - **numint** — one- and multi-dimensional numerical integration (adaptive,
    Gauss, Romberg, Monte Carlo and cubature rules).
  - **optimize** — unconstrained minimizers, numeric gradients/Hessians/
    Jacobians and a suite of one-dimensional root finders (Newton, Halley,
    Schröder, Steffensen and more).
  - **interp** — polynomial and spline interpolation.
  - **geom2d** — planar geometric primitives and shapes.

### Fixed
- Reconciled the duplicated exported symbols and shared helpers that the two
  parallel authors of each new package had defined twice, keeping a single
  complete implementation of each so the whole module builds, vets and tests
  cleanly.
- **special**: `specialEulerGamma` and the Bessel/Airy routines are now declared
  once (in `bessel.go`); `integrals.go` retains only the unique `BesselJ/Y/I/K`
  order-n wrappers and the four-value `Airy` function.
- **optimize**: `NumericHessian` now computes correct diagonal second
  derivatives (the mixed stencil previously collapsed to zero on the diagonal),
  and `Schroeder` recognises convergence at a multiple root instead of failing
  with a zero-derivative error.
- **orthopoly**: corrected the `GegenbauerDerivative` expectation (d/dx C₂¹ =
  8x, so 4 at x = ½) and added the missing Gauss–Laguerre weight assertions.
- **geom2d**: resolved the clashing `pointApprox` test helpers.

## [0.5.0] - 2026-07-18
### Changed
- Integrated the parallel `matrix`, `ntheory`, `stats` and `physics` additions
  into a single building, vetting, linting and testing tree.
- Completed API documentation: every exported symbol across all packages —
  including previously undocumented distribution and result struct fields
  (`Normal`, `Binomial`, `Beta`, `FDist`, `Weibull`, `Hypergeometric`,
  `TestResult`, `BootstrapResult`, `ODEError.Reason`, and the CAS expression
  fields) — now carries a doc comment (100% godoc).

### Fixed
- **physics**: corrected two new tests that encoded wrong expected values —
  `MagneticFieldStraightWire` is now checked against the package's CODATA
  vacuum permeability (μ0/2π ≈ 2.0000000011×10⁻⁷) rather than the pre-2019
  exact 2×10⁻⁷, and `VelocityFromDistance(0,2,8)` now expects √32 = 4√2 for
  v = √(v0²+2ad).
- **staticcheck** cleanups so `golangci-lint run ./...` passes with zero issues:
  renamed unexported error vars to the `errFoo` form (`errQRNoConverge`,
  `errEmptyNetwork`, `errZeroElement`), replaced the deprecated `doc.Synopsis`
  in the docs generator with `Package.Synopsis`, escaped the invisible MathML
  operators in `latex.go` as `⁡`/`⁢`, used `copy` in
  `Vector.RowMatrix`, converted `LogNormal` to `Normal` via a type conversion,
  and made the Pollard-Brent determinism test compare two distinct evaluations.

## [0.4.0] - 2026-07-18
### Added
- **matrix** sub-package: dense linear-algebra expansion — LU/QR/Cholesky and
  singular-value decompositions (`SVD`, `SingularValues`), symmetric and general
  eigensolvers, numerical `RankNumeric`, pseudoinverse and least-squares
  (`lstsq`), matrix functions (`matfunc`), a family of matrix/vector `norms`,
  and cache-friendly dense kernels.
- **ntheory** sub-package: overflow-safe uint64 fast paths — deterministic
  `IsPrimeU64` (Montgomery-accelerated Miller-Rabin with a proven 7-base witness
  set), `NextPrimeU64`/`PrevPrimeU64`/`PrevPrime` and Baillie-PSW
  `IsProbablePrimeBig`/`NextPrimeBig`/`PrevPrimeBig`; a segmented `SegmentedSieve`
  with a streaming `PrimeSieve`, `PrimePiRange` and `NthPrime`; `MulModU64` plus
  `Montgomery`/`Barrett` fast modular arithmetic; Pollard-Brent factorization
  (`PollardBrentU64`, `FactorizeU64`, `FactorizeBig`, `FactorListUint64`);
  modular square roots (`SqrtMod`, `SqrtModBig`, `AllSqrtModComposite`);
  primitive roots and `MultiplicativeOrder`; generalized divisor functions
  (`DivisorSigma`, `TotientSieve`, `MobiusSieve`, `MertensFunction`); and
  continued fractions with `PellFundamental`.
- **stats** sub-package: more distributions (`Beta`, `FDist`, `LogNormal`,
  `Weibull`, `Geometric`, `NegativeBinomial`, `GeometricFailures`,
  `NegativeBinomialInt`, `Hypergeometric`) with PDF/PMF/CDF/quantile and seeded
  sampling; streaming moments (`Accumulator`, `CovAccumulator`); hypothesis
  tests (t-tests, chi-square, one-way ANOVA, Mann-Whitney U); multiple linear
  and ridge regression with covariance/correlation matrices; and confidence
  intervals plus bootstrap resampling.
- **physics** sub-package: zero-allocation `Vec3` algebra; kinematics, orbital
  mechanics and 1D collisions; special relativity (`Beta`, Lorentz factor,
  relativistic energy/momentum, velocity addition); ideal-gas and heat-transfer
  thermodynamics; electromagnetism (Coulomb field/potential, resistor and
  capacitor networks); and measurement-uncertainty propagation.

## [0.3.0] - 2026-07-15
### Added
- Root CAS expansion: multivariate polynomial tools (`poly`), first- and
  second-order ODE solving (`odesolve`), and LaTeX rendering of expressions
  (`latex`).

## [0.2.0] - 2026-07-12
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
- New `matrix` subpackage: symbolic + numeric linear algebra — `Matrix`/`Vector`,
  exact `Det`/`Inverse`/`RREF`/`Solve`, `Adjugate`, `Rank`, `Transpose`, `Kron`,
  `Dot`/`Cross`/`Norm`, and `CharPoly`/`Eigenvalues`.
- New `ntheory` subpackage: number theory & combinatorics — `GCD`/`LCM`/`ExtendedGCD`,
  `IsPrime` (deterministic Miller-Rabin), `Factorize`, `PrimesUpTo`, `EulerPhi`,
  `MobiusMu`, modular arithmetic (`ModPow`/`ModInverse`/`CRT`/`Jacobi`/`DiscreteLog`),
  `Binomial`/`Factorial`/`Catalan`/`Partition`, and `Fibonacci`/`Lucas`/`Bernoulli`.
- New `stats` subpackage: descriptive statistics (`Mean`/`Median`/`Variance`/`StdDev`/
  `Correlation`/`LinearRegression`) and probability distributions (`Normal`, `Binomial`,
  `Poisson`, `Uniform`, `Exponential`, `StudentT`, `ChiSquared`, `Gamma`).
- New `physics` subpackage: SI-2019/CODATA physical constants (enumerable via
  `Constants()`/`Lookup`), unit conversions across length/mass/time/temperature/
  energy/angle, and common physics formulas (kinematics, relativity, waves, EM).

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

[0.2.0]: https://github.com/malcolmston/algebra/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/malcolmston/algebra/releases/tag/v0.1.0
