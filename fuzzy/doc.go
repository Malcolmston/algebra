// Package fuzzy implements fuzzy logic and fuzzy set theory using only the Go
// standard library.
//
// The package provides membership functions (triangular, trapezoidal,
// Gaussian, sigmoid, generalized bell, S/Z/Pi shapes and singletons), the
// standard family of t-norms and t-conorms (minimum/maximum, algebraic
// product/probabilistic sum, Lukasiewicz, drastic, Hamacher, Einstein and
// nilpotent variants), fuzzy complements (standard, Sugeno, Yager), discrete
// fuzzy sets with set-theoretic operations, alpha-cuts, linguistic hedges,
// scalar cardinality and convexity tests, fuzzy relations with max-min and
// max-product (and general max-t) composition, Mamdani and Sugeno (TSK)
// inference systems, and a full complement of defuzzification methods
// (centroid, bisector, mean/smallest/largest of maxima and weighted average).
//
// A membership function is modeled by the MF type, a plain func(float64)
// float64 that maps a point of the universe of discourse to a membership
// grade clamped to the closed interval [0, 1]. A discrete fuzzy set is modeled
// by Set, which stores a sorted universe together with the membership grade at
// each point. T-norms and t-conorms are modeled by the TNorm and TConorm
// function types so that inference and set operations can be parameterized by
// any binary aggregation operator.
//
// All computations are deterministic. The few routines that accept randomness
// take a caller supplied *math/rand.Rand so results are reproducible.
package fuzzy
