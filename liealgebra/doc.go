// Package liealgebra implements computational tools for Lie algebras and root
// systems using only the Go standard library.
//
// The package is organised around a handful of small, self-contained numeric
// types and a large collection of genuinely distinct operations on them:
//
//   - Dense real and complex matrices ([Matrix], [CMatrix]) with the arithmetic
//     needed for matrix Lie algebras: addition, multiplication, transpose,
//     conjugate transpose, trace, Frobenius norm, LU based determinant, inverse
//     and linear solve.
//
//   - The Lie bracket (commutator) [Bracket] and anticommutator, the adjoint
//     action [AdjointAction], the Jacobi identity residual and test
//     ([JacobiResidual], [SatisfiesJacobi]), structure constants
//     ([StructureConstants]), the Killing form ([KillingForm], [KillingFormValue])
//     and the trace form, together with semisimplicity tests.
//
//   - Concrete generators of the low dimensional Lie algebras: sl(2) with its
//     standard E, F, H basis ([SL2Generators]), so(3) rotation generators
//     ([SO3Generators]), su(2) via the Pauli matrices ([SU2Generators],
//     [PauliMatrices]), su(3) via the Gell-Mann matrices ([GellMannMatrices]),
//     and arbitrary spin-j representations ([SpinMatrices]).
//
//   - The exponential map for matrices ([MatExp], [CMatExp]) implemented with
//     scaling and squaring and a diagonal Padé approximant, plus integer powers
//     and the leading Baker-Campbell-Hausdorff terms.
//
//   - Root system combinatorics for the classical families A, B, C and D and
//     the exceptional types: Cartan matrices ([CartanMatrix]), simple and
//     positive roots ([SimpleRoots], [PositiveRoots]), fundamental weights,
//     Weyl group orders ([WeylGroupOrder]), Coxeter numbers, Dynkin diagram
//     adjacency, Weyl reflections and the Weyl dimension and Casimir formulas.
//
// All computation is deterministic. Operations that can fail (dimension
// mismatches, singular systems, unknown Dynkin types) return a sentinel error
// such as [ErrDim], [ErrSingular] or [ErrType]; a few pure accessors panic on
// out-of-range indices.
package liealgebra
