// Package projectivegeom implements projective and inversive geometry in the
// real projective plane RP^2 and real projective space RP^3.
//
// The plane RP^2 is modelled with homogeneous coordinates: a [Point] is a
// non-zero triple [X Y Z] identified up to non-zero scale, and a [Line] is a
// non-zero triple [A B C] denoting the locus A*X+B*Y+C*Z=0. Points and lines
// are dual: joining two points ([Join]) and meeting two lines ([Meet]) are the
// same cross-product operation, and every construction has a dual obtained by
// swapping the two roles.
//
// On these primitives the package builds the classical toolbox of projective
// geometry: incidence, collinearity and concurrency predicates; the projective
// invariant [CrossRatio] and the [HarmonicConjugate] construction; homographies
// ([Homography], a 3x3 [Mat3]) with composition, inversion and the standard
// affine, similarity and Euclidean special cases; conics ([Conic], a symmetric
// 3x3 matrix) built from five points ([ConicFromFivePoints]) with tangents,
// pole/polar duality, intersection and affine classification; the Pappus and
// Desargues configuration theorems; and camera-style 3x4 projection matrices
// ([Camera]) mapping RP^3 to RP^2.
//
// Space RP^3 is modelled analogously with 4-vectors ([SPoint], [SPlane]) and
// Pluecker line coordinates ([LinePluecker]).
//
// All computation is performed in float64 using only the Go standard library.
// Results are deterministic. Predicates that must tolerate floating-point noise
// take or use an explicit tolerance so the caller controls the trade-off
// between robustness and strictness. Any routine that requires randomness takes
// a caller-supplied *math/rand.Rand so results are reproducible.
package projectivegeom
