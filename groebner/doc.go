// Package groebner implements exact computer algebra over sparse multivariate
// polynomials with rational coefficients (math/big.Rat).
//
// The package provides a full toolkit for ideal theory in the polynomial ring
// Q[x_1, ..., x_n]: monomials as integer exponent vectors, configurable
// monomial orders (lexicographic, graded lexicographic, graded reverse
// lexicographic, weighted and block/elimination orders), polynomial
// arithmetic (addition, subtraction, multiplication, powers, scalar and
// monomial multiplication), evaluation over the rationals and over the complex
// numbers, formal derivatives, multivariate division with remainder,
// S-polynomials, and Buchberger's algorithm (with the coprime-leading-term and
// chain criteria) for computing Gröbner bases.
//
// On top of Gröbner bases the package computes reduced and minimal Gröbner
// bases, decides ideal membership, and performs the standard ideal operations:
// sum, product, intersection, quotient (colon ideal), and elimination ideals.
// For zero-dimensional systems it can decide finiteness of the variety and
// numerically approximate all complex solutions by triangular back-solving
// combined with a self-contained Durand–Kerner univariate root finder.
//
// All arithmetic on coefficients is exact. Only the numerical variety solver
// works with floating point (complex128) and it takes a caller-supplied random
// seed so results are reproducible.
package groebner
