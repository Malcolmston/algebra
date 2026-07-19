// Package ellipticcurves implements the arithmetic and number theory of
// elliptic curves in short Weierstrass form y^2 = x^3 + A*x + B, both over a
// prime finite field F_p and over the field of rational numbers Q. The focus is
// mathematical rather than cryptographic: the package is meant for exploring the
// group law, invariants, torsion, point counting and pairings of small curves.
//
// Over F_p the central types are CurveFp and its affine points PointFp. They
// support the geometric group law (addition, doubling, negation and scalar
// multiplication), the curve invariants discriminant and j-invariant, point
// counting (naive enumeration, baby-step/giant-step point order and a
// Mestre-style baby-step/giant-step curve order), division polynomials and the
// Weil pairing.
//
// Over Q the central types are CurveQ and PointQ, built on math/big rationals.
// They support the group law, invariants, the Nagell-Lutz determination of the
// torsion subgroup, naive and canonical heights, and heuristic rank lower
// bounds derived from the height pairing regulator of a set of rational points.
//
// Isomorphism and twist utilities relate curves with equal j-invariant over
// both fields.
//
// The implementation depends only on the Go standard library. Randomness, where
// used, is drawn from a caller-supplied math/rand source so that every routine
// is reproducible. Functions do not mutate their big.Int or big.Rat arguments;
// results are always freshly allocated.
package ellipticcurves
