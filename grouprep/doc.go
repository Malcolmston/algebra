// Package grouprep implements the ordinary representation theory of finite
// groups over the complex numbers, built entirely on the Go standard library.
//
// A (matrix) representation of a finite group G is a homomorphism
// ρ: G → GL(d, ℂ) assigning to every group element an invertible d×d complex
// matrix such that ρ(gh) = ρ(g)ρ(h). This package provides the concrete data
// structures and algorithms of the classical theory: complex matrices
// ([Matrix]), finite groups given by their Cayley table ([Group]), matrix
// representations ([Rep]), characters ([Character]) and character tables
// ([CharacterTable]).
//
// # Groups
//
// A [Group] is stored as an explicit multiplication (Cayley) table together
// with an inverse map and identity index. Standard families are provided:
// [CyclicGroup], [DihedralGroup], [SymmetricGroup], [AlternatingGroup],
// [QuaternionGroup], [KleinFourGroup], [TrivialGroup] and the
// [DirectProduct] of two groups. Structural queries include
// [Group.ConjugacyClasses], [Group.Center], [Group.IsAbelian],
// [Group.ElementOrder] and subgroup closures via [Group.GeneratedBy].
//
// # Representations
//
// A [Rep] pins a group and a matrix for every element. Constructions include
// [TrivialRep], [RegularRep], [CyclicRep], [DihedralRep2D],
// [NaturalRepSymmetric], [SignRepSymmetric] and [QuaternionRep]. Reps combine
// under [Rep.DirectSum] and [Rep.TensorProduct] and can be dualised with
// [Rep.Dual]. Every rep is validated as a genuine homomorphism.
//
// # Characters
//
// The character of a rep is the trace class function χ(g) = tr ρ(g), obtained
// with [Rep.Character]. Characters live in the Hermitian inner product space of
// class functions with 〈χ, ψ〉 = (1/|G|) Σ χ(g)⁻ ψ(g), computed by
// [InnerProduct]. A character is irreducible exactly when its norm squared is
// one ([IsIrreducible]); an arbitrary character decomposes uniquely into
// irreducibles ([DecomposeCharacter]) — the numerical content of Maschke's
// theorem and complete reducibility. Restriction and induction between a group
// and a subgroup are provided by [RestrictCharacter] and [InduceCharacter],
// and satisfy Frobenius reciprocity.
//
// # Character tables
//
// A [CharacterTable] assembles a complete set of irreducible characters against
// the conjugacy classes of the group. Tables for the standard families are
// produced by [CharacterTableCyclic], [CharacterTableDihedral],
// [CharacterTableSymmetric], [CharacterTableQuaternion] and
// [CharacterTableKleinFour]. Each table satisfies the row and column
// orthogonality relations ([CharacterTable.RowOrthogonality],
// [CharacterTable.ColumnOrthogonality]) and the dimension identity
// Σ dᵢ² = |G| ([CharacterTable.SumOfSquaresOfDims]).
//
// All results are deterministic. Any randomness used by helper routines takes a
// caller-supplied seed through math/rand. Floating-point comparisons take an
// explicit tolerance; a value near 1e-9 is appropriate for the small groups the
// package targets.
package grouprep
