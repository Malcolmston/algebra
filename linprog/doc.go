// Package linprog implements linear and convex optimization in pure Go using
// only the standard library.
//
// The package is organized around a small number of problem types and the
// classical algorithms that solve them:
//
//   - [LP] holds a general linear program (minimize or maximize a linear
//     objective subject to linear <=, =, >= constraints with x >= 0). It can be
//     reduced to [StandardLP] equality/nonnegative form and solved with the
//     two-phase primal [Simplex] method, which uses Bland's rule to guarantee
//     termination without cycling.
//
//   - Integer programs are solved by branch and bound on the LP relaxation via
//     [LP.SolveInteger].
//
//   - Duality utilities build the [DualCanonical] program and check
//     complementary slackness and the LP KKT conditions.
//
//   - The transportation and assignment problems are solved by [Transportation]
//     (reducing to a standard LP) and the O(n^3) [Hungarian] algorithm.
//
//   - Convex quadratic programs are solved either exactly for the
//     equality-constrained case ([SolveQPEquality], a single KKT linear system)
//     or by a primal [SolveQP] active-set method for inequality constraints,
//     with [ComputeKKTResidual] reporting the first-order optimality residuals.
//
// All routines are deterministic and depend only on the Go standard library.
// Vectors are []float64 and matrices are row-major [][]float64.
package linprog

// Eps is the default numerical tolerance used throughout the package for
// pivot selection, feasibility checks and reduced-cost sign tests.
const Eps = 1e-9
