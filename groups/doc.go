// Package groups is a self-contained finite-group and modular-arithmetic
// toolkit built entirely on the Go standard library.
//
// It is a sibling subpackage of github.com/malcolmston/algebra but does not
// depend on it: every routine works directly on machine integers, small
// permutation slices, or explicit finite-group tables. All results are
// deterministic and validated in the accompanying tests against closed-form
// reference values.
//
// # Contents
//
// Integer gcd domain and the Euclidean algorithm: [Gcd], [ExtendedGcd], [Lcm],
// [GcdMany], [LcmMany], [Coprime].
//
// Modular ring Z/nZ: [Mod], [ModAdd], [ModSub], [ModMul], [ModNeg], [ModPow],
// [ModInverse], [CyclicOrder], [ElementOrderZn], [MultiplicativeOrder],
// [EulerTotient], [UnitsModN], [IsUnitModN], [IsPrimitiveRoot],
// [PrimitiveRoots].
//
// Prime field GF(p): [IsPrime], [GFAdd], [GFSub], [GFNeg], [GFMul], [GFInv],
// [GFDiv], [GFPow].
//
// Permutations and the symmetric/alternating groups: the [Perm] type with
// [Identity], [Compose], [Perm.Inverse], [Perm.Order], [Perm.Sign],
// [Perm.Cycles], [Perm.CycleType], [PermFromCycles], [Transposition],
// [NumInversions], [SymmetricGroup], [AlternatingGroup],
// [SymmetricGroupOrder], [AlternatingGroupOrder], [Factorial].
//
// Dihedral groups: [DihedralElement] with [DihedralRotation],
// [DihedralReflection], [DihedralIdentity], [DihedralCompose],
// [DihedralInverse], [DihedralElementOrder], [DihedralOrder], [DihedralGroup].
//
// Abstract finite groups (Cayley tables, subgroups, cosets): [FiniteGroup]
// with [NewFiniteGroup], [CyclicGroupZn], [SymmetricGroupTable],
// [FiniteGroup.CayleyTable], [FiniteGroup.Inverse], [FiniteGroup.ElementOrder],
// [FiniteGroup.IsAbelian], [FiniteGroup.IsCyclic], [FiniteGroup.Center],
// [FiniteGroup.GeneratedSubgroup], [FiniteGroup.IsSubgroup],
// [FiniteGroup.LeftCoset], [FiniteGroup.RightCoset], [FiniteGroup.Cosets],
// [FiniteGroup.Index], [FiniteGroup.IsValid].
//
// Polynomials over a field, gcd domain via the Euclidean algorithm: the [Poly]
// type over the rationals ([PolyAdd], [PolyMul], [PolyDivMod], [PolyGCD],
// [PolyEval], ...) and [][]int polynomials over GF(p) ([PolyModAdd],
// [PolyModMul], [PolyModDivMod], [PolyModGCD], [PolyModEval]).
//
// # Conventions
//
// Permutations are zero-based: a [Perm] p of length n is a bijection of
// {0,1,...,n-1} with p[i] the image of i. Modulus arguments must be positive
// and GF(p) routines require p to be prime; violating a documented
// precondition panics rather than returning a silently wrong value.
package groups
