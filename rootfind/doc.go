// Package rootfind implements polynomial and nonlinear root finding in pure Go.
//
// The package is organized around two polynomial representations and a family
// of solvers that operate on them or on arbitrary scalar functions.
//
// # Polynomial representation
//
// Polynomials are stored as coefficient slices in ascending order of power: the
// slice p represents the polynomial
//
//	p(x) = p[0] + p[1]*x + p[2]*x^2 + ... + p[n]*x^n
//
// so that p[i] is the coefficient of x^i and the last element is the leading
// coefficient. [Poly] holds real coefficients ([]float64) and [CPoly] holds
// complex coefficients ([]complex128). This ordering makes indexing by power
// trivial and matches Horner evaluation from the high-order end.
//
// # Global polynomial solvers
//
// Several methods find all roots of a polynomial at once. [DurandKerner]
// (the Weierstrass method) and [AberthEhrlich] are simultaneous iterations that
// converge to every complex root in parallel; [Bairstow] extracts real
// quadratic factors to recover complex-conjugate pairs using only real
// arithmetic; and [CompanionEigenvalues] finds roots as the eigenvalues of the
// companion matrix via the Francis double-shift QR algorithm. [PolyRoots] is a
// high-level convenience that returns all roots, and [RealRoots] filters those
// with negligible imaginary part.
//
// # Real-root theory
//
// [SturmSequence] and the counting routines built on it give the exact number
// of distinct real roots in any interval, from which [IsolateRoots] produces
// disjoint bracketing intervals that are then refined to full precision. The
// classical sign-based theorems are provided too: [DescartesRuleOfSigns],
// [BudanFourierCount], and a collection of a-priori root-magnitude bounds such
// as [CauchyBound], [LagrangeBound], and [FujiwaraBound].
//
// # Scalar solvers
//
// For a general continuous function of one variable the package offers bracketing
// methods ([Bisection], [FalsePosition], [Brent], [Ridders]), open methods
// ([Secant], [Steffensen], [FixedPoint]), and derivative-based methods
// ([Newton], [Halley]). [Laguerre] is a globally very reliable polynomial
// solver. Each returns a [Result] recording the located root, the residual, the
// iteration count, and whether convergence was achieved.
//
// # Multiplicity and deflation
//
// gcd(p, p') isolates repeated factors: [SquareFree] returns the squarefree
// part, [SquareFreeFactorization] performs Yun's algorithm to group roots by
// multiplicity, and [Multiplicity] measures the order of a specific root.
// Roots may be removed one at a time with [Poly.DeflateReal] or [CPoly.Deflate].
//
// All algorithms use only the Go standard library and perform no allocation
// beyond what their results require.
package rootfind
