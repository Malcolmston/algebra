// Package analyticnt is a self-contained analytic number theory toolkit built
// entirely on the Go standard library.
//
// It is a sibling subpackage of github.com/malcolmston/algebra but does not
// depend on it. Everything here is implemented directly on machine floats,
// integers, and math/big / math/cmplx values.
//
// # Contents
//
// Prime enumeration and exact counting: [Sieve], [PrimesUpTo], [NthPrime],
// [IsPrime], [PrimePi], [PrimePiLegendre], [PrimePiMeissel], [LegendrePhi].
//
// Prime-counting approximations: [Li], [LiOffset], [Ei], [E1], [RiemannR],
// [PrimePiApprox], [NthPrimeApprox], [SoldnerConstant].
//
// Chebyshev and von Mangoldt functions: [ChebyshevTheta], [ChebyshevPsi],
// [VonMangoldt], [VonMangoldtLambda], [SecondChebyshev].
//
// The Riemann zeta function: [Zeta], [ZetaComplex], [DirichletEta],
// [DirichletEtaComplex], [RiemannSiegelTheta], [RiemannSiegelZ], [ZetaZero],
// [ZetaZeros], [HardyZ], [GramPoint].
//
// Dirichlet characters and L-functions: [DirichletCharacter],
// [PrincipalCharacter], [LegendreCharacter], [DirichletL], [LFunctionReal].
//
// Multiplicative functions and Mertens: [MobiusMu], [MertensFunction],
// [EulerPhi], [Liouville], [Omega], [BigOmega], [MertensConstant],
// [TwinPrimeConstant], [BrunConstant], [ArtinConstant].
//
// Prime-number-theorem error terms and prime-gap statistics: [PiError],
// [PsiError], [PrimeGaps], [MaximalPrimeGaps], [TwinPrimesUpTo],
// [CousinPrimesUpTo], [SexyPrimesUpTo], [PrimeCountingBias].
//
// # Conventions
//
// Real-valued approximation routines return float64 and document the input
// domain on which they are valid. Functions that violate a documented
// precondition (for example a negative count) panic rather than returning a
// silently wrong value. Any routine that consumes randomness takes an explicit
// int64 seed and uses a deterministic math/rand source, so results are fully
// reproducible.
package analyticnt
