// Package matrix provides exact symbolic and numeric linear algebra built on
// top of the github.com/malcolmston/algebra expression engine.
//
// Every entry of a [Matrix] (and every component of a [Vector]) is an
// algebra.Expr, so a matrix can hold arbitrary-precision integers, exact
// rationals, symbols and full symbolic expressions side by side. All
// arithmetic keeps results exact by canonicalizing every entry through
// algebra.Simplify, which means operations such as [Matrix.Det],
// [Matrix.Inverse], [Solve] and [Matrix.CharPoly] return closed-form,
// deterministic answers rather than floating-point approximations.
//
// # Design
//
//   - Values are immutable from the caller's perspective: every operation
//     returns a fresh matrix or vector and never mutates its receiver or
//     arguments.
//   - Determinants and adjugates use cofactor (Laplace) expansion, which never
//     divides and therefore stays exact for symbolic entries. This is O(n!) and
//     intended for the small matrices typical of symbolic work.
//   - Row reduction ([Matrix.RREF], [Matrix.Rank], [Solve]) treats an entry as
//     zero only when it structurally simplifies to zero. Pivots are assumed
//     invertible; symbolic pivots produce rational-function entries.
//   - A numeric float64 fast path is available via [FromFloats] and
//     [Matrix.Floats] and is used internally for the numeric eigenvalue
//     fallback.
//
// # Eigenvalues
//
// [Matrix.CharPoly] returns the characteristic polynomial det(A - λI) as an
// algebra.Expr for a square matrix of any size. [Matrix.Eigenvalues] returns
// exact eigenvalues for 1×1 and 2×2 matrices (via the parent package's Solve),
// and for 3×3 matrices when a rational eigenvalue can be peeled off to leave a
// quadratic. When no exact factorization is available the eigenvalues are
// computed numerically and returned as algebra.Flt values; only real roots are
// reported. Matrices of size 4×4 and larger return an error from Eigenvalues —
// callers can still obtain the characteristic polynomial and solve it
// themselves.
package matrix
