// Package quasirandom implements low-discrepancy (quasi-random) sequences and
// the numerical machinery that surrounds them.
//
// Quasi-random, or quasi-Monte-Carlo (QMC), point sets fill the unit cube
// [0,1)^d far more evenly than independent pseudo-random samples. Their
// defining quantity is the discrepancy, a measure of how far the empirical
// distribution of the points departs from the uniform distribution. The
// Koksma–Hlawka inequality bounds the error of an equal-weight cubature rule by
// the product of the star discrepancy of the nodes and the variation of the
// integrand, so sequences with asymptotically small discrepancy — of order
// (log N)^d / N rather than the 1/sqrt(N) of plain Monte Carlo — give faster,
// deterministic convergence.
//
// The package provides:
//
//   - The radical inverse and the van der Corput sequence in any base, with
//     digit scrambling and permutation-based (Faure–Tezuka style) scrambling.
//   - Halton sequences in arbitrary dimension (prime and custom bases, leaped
//     and scrambled variants) and the finite Hammersley point sets.
//   - Faure sequences built from the generalized Pascal generator matrices over
//     a prime base at least as large as the dimension.
//   - Sobol sequences generated from primitive polynomials over GF(2) and their
//     direction numbers, using the Gray-code recurrence, with a built-in table
//     of direction numbers and support for user-supplied ones.
//   - Closed-form L2-type discrepancies (Warnock, centered, wrap-around and
//     symmetric), the exact one-dimensional star discrepancy and an exact
//     grid-based star discrepancy for small point sets.
//   - QMC integration helpers over the unit cube together with the
//     Koksma–Hlawka error bound.
//
// Every generator is fully deterministic: no source of randomness is consulted
// and results depend only on the requested indices, bases and dimensions.
// Only the Go standard library is used.
package quasirandom
