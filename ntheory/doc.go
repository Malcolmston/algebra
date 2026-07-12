// Package ntheory is a self-contained number-theory and combinatorics toolkit
// built entirely on the Go standard library (only math and math/big).
//
// It is a sibling subpackage of github.com/malcolmston/algebra but does not
// depend on it: every routine here works directly on machine integers or on
// arbitrary-precision math/big values.
//
// # Design
//
// Two numeric styles are offered depending on what is natural for each
// problem:
//
//   - int64 conveniences for the small, fast integer functions where 64 bits
//     comfortably hold every intermediate value (divisibility, primality,
//     modular arithmetic). Internally these use math/big where an intermediate
//     product could otherwise overflow (for example [ModPow]), so the results
//     are always correct for any inputs that fit in an int64.
//   - *big.Int / *big.Rat results for functions whose outputs grow without
//     bound, such as [Factorial], [Binomial], [Fibonacci], [Partition] and
//     [Bernoulli].
//
// # Contents
//
// Divisibility: [GCD], [LCM], [ExtendedGCD], [Divisors], [SumDivisors],
// [CountDivisors], [IsPerfect].
//
// Primality & factorization: [IsPrime], [IsPrimeBig], [NextPrime],
// [PrimesUpTo], [PrimePi], [Factorize], [FactorList], [EulerPhi], [MobiusMu],
// [Radical].
//
// Modular arithmetic: [ModPow], [ModInverse], [CRT], [LegendreSymbol],
// [JacobiSymbol], [IsQuadraticResidue], [DiscreteLog], [Order].
//
// Combinatorics: [Factorial], [DoubleFactorial], [Binomial], [Multinomial],
// [Permutations], [CatalanNumber], [StirlingSecond], [Partition].
//
// Sequences: [Fibonacci], [Lucas], [Tribonacci], [IsSquare], [IsqrtBig],
// [Bernoulli].
//
// # Conventions
//
// Unless stated otherwise, functions taking a modulus require it to be
// positive, and functions taking a "count" argument (n for a factorial, an
// index into a sequence, and so on) require it to be non-negative. Violating a
// documented precondition panics rather than returning a silently wrong value.
package ntheory
