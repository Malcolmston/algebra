// Package ecc implements elliptic-curve arithmetic in short Weierstrass form
// y^2 = x^3 + A*x + B, both over a prime finite field GF(p) and over the field
// of rational numbers Q.
//
// The package provides the geometric group law (point addition, doubling,
// negation and scalar multiplication), curve invariants (discriminant and
// j-invariant), point-counting utilities (naive enumeration and a
// baby-step/giant-step point-order routine), SEC1 point compression and
// decompression, and the finite-field number theory required to build higher
// level protocols on top: an Elliptic-Curve Diffie-Hellman shared-secret
// derivation and the ECDSA signing and verification equations over GF(p).
//
// Two concrete curve types are exposed:
//
//   - CurveFp models a curve over the prime field GF(p) using math/big
//     integers. Its points are affine PointFp values with an explicit
//     point-at-infinity flag.
//   - CurveQ models a curve over the rationals Q using math/big rationals. Its
//     points are affine PointQ values, likewise carrying an infinity flag.
//
// All routines are deterministic and depend only on the Go standard library.
// Functions never mutate their big.Int or big.Rat arguments; every result is a
// freshly allocated value.
package ecc
