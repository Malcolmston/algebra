// Package galois implements arithmetic and structure theory for finite fields
// GF(p^n), built entirely on the Go standard library.
//
// A finite field (Galois field) of order q = p^n exists for every prime p and
// positive integer n. This package models both the prime field GF(p) and the
// extension fields GF(p^n), providing element arithmetic (add, subtract,
// multiply, invert, divide, exponentiate) together with the deeper structural
// tools that make finite fields useful in coding theory and cryptography.
//
// # Representation
//
// Prime-field elements are ordinary residues handled by the modular-arithmetic
// helpers ([AddMod], [MulMod], [InvMod], [PowMod], …) operating on
// math/big integers reduced into the range [0, p).
//
// Polynomials over GF(p) are represented by [Poly], a little-endian slice of
// coefficients paired with the characteristic p. Extension fields are modelled
// by [Field], which pins a prime p, a degree n, and a monic irreducible
// modulus polynomial of degree n; elements are [Element] values, residues in
// GF(p)[x] modulo that polynomial.
//
// # Capabilities
//
//   - Prime-field modular arithmetic, quadratic residues and modular square
//     roots ([SqrtMod], Tonelli–Shanks), multiplicative order and primitive
//     roots.
//   - Dense polynomial arithmetic over GF(p): division, gcd, extended gcd,
//     modular exponentiation, derivatives, evaluation and composition.
//   - Irreducibility testing (Rabin), counting irreducibles, and searching for
//     irreducible and primitive polynomials, including a deterministic
//     Conway-style canonical choice.
//   - Full polynomial factorisation over GF(p): square-free decomposition and
//     Berlekamp's deterministic algorithm.
//   - Extension-field element arithmetic, the Frobenius automorphism, trace and
//     norm to the prime field, minimal polynomials, conjugates and cyclotomic
//     cosets.
//   - Multiplicative order of field elements, primitive elements, and discrete
//     logarithms by the baby-step giant-step method.
//   - The subfield lattice of GF(p^n).
//
// Everything is exact: computations run over math/big integers, so results are
// mathematically sound rather than floating-point approximations.
package galois
