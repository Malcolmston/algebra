// Package knottheory implements combinatorial knot and braid theory using only
// the Go standard library.
//
// The package provides self-contained, mathematically correct implementations
// of the classical objects and invariants of low-dimensional topology as they
// appear in a combinatorial (diagrammatic) setting:
//
//   - Laurent polynomials over the integers (the Laurent type), the coefficient
//     ring in which every polynomial knot invariant lives, with a full set of
//     ring operations, evaluation, substitution and normalisation helpers.
//   - Permutations (the Permutation type): composition, inversion, cycle
//     decomposition, order and sign, used to describe the underlying
//     permutation and closure of a braid.
//   - Braid groups (the Braid type): Artin generators and relations, free
//     reduction, exponent sum (writhe), the induced permutation, mirror and
//     reverse, torus braids, the full twist and the Garside half twist, and the
//     reduced Burau representation over the Laurent ring.
//   - Knot and link diagrams described by signed Gauss codes (the GaussCode and
//     Diagram types): validation, crossing number, writhe, self-writhe and the
//     Gauss linking number between components.
//   - Planar-diagram (PD) codes (the PDCode type): the Kauffman bracket via a
//     union-find state sum, the writhe-normalised bracket, and the Jones
//     polynomial for small diagrams, both in the standard variable t and in the
//     variable t^{1/2} for links.
//   - The Alexander polynomial, computed two independent ways: from the reduced
//     Burau matrix of a braid, and directly from a diagram via the Alexander
//     (Fox-derivative) matrix, together with the knot determinant and
//     p-colourability counts.
//   - Reidemeister-move recognisers (types I, II and III) that operate
//     combinatorially on Gauss codes.
//   - Torus-knot invariants T(p, q): crossing number, Seifert genus, bridge and
//     unknotting numbers, determinant, and the closed forms of the Alexander and
//     Jones polynomials.
//
// Unless documented otherwise polynomials are represented by the Laurent type,
// braids by words of non-zero integers (generator i is +i, its inverse is -i),
// and diagrams by Gauss or PD codes whose crossings are labelled by consecutive
// integers starting at one.
package knottheory
