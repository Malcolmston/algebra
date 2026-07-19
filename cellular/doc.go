// Package cellular implements cellular automata in pure Go using only the
// standard library.
//
// The package covers the classical families of discrete dynamical systems
// studied by Wolfram, Conway and others:
//
//   - Elementary one-dimensional automata ([ElementaryRule]) with the standard
//     Wolfram numbering 0–255, single-step and spacetime evolution, the four
//     symmetry conjugates (mirror, complement and mirror-complement), the 88
//     inequivalent equivalence classes, additivity and totalistic tests, the
//     Langton lambda parameter and a Wolfram-class heuristic ([ElementaryRule.Class]).
//
//   - k-state, r-radius totalistic ([TotalisticRule]) and outer-totalistic
//     ([OuterTotalisticRule]) rules, together with fully general lookup-table
//     rules ([GeneralRule]). All of these satisfy the [Rule1D] interface and are
//     driven by the shared [Step1D] and [Evolve1D] evolvers under a selectable
//     [Boundary] condition (periodic, fixed-zero, fixed-one or reflecting).
//
//   - Two-dimensional life-like automata built from a [LifeRule] (a B/S
//     rulestring such as "B3/S23") on a toroidal or bounded [Grid], with the
//     Moore and von Neumann neighbourhoods and a library of named rules and
//     starting patterns.
//
//   - Reversible automata: second-order (Fredkin) elementary rules
//     ([SecondOrderCA]) that run forwards and backwards exactly, and
//     block/Margolus rules ([MargolusRule]) whose reversibility is guaranteed by
//     using a permutation of the block alphabet.
//
//   - Spacetime-diagram utilities: ASCII rendering, densities, Hamming distance
//     and damage spreading, and Shannon/block entropies used both for analysis
//     and for the class heuristics.
//
// States are represented as []int slices (values 0..k-1) for one-dimensional
// automata and by the [Grid] type in two dimensions. No third-party packages,
// randomness sources or cgo are used; where pseudo-random initial conditions are
// needed they are produced by a small deterministic splitmix64 generator seeded
// by an explicit uint64 so that every result in the package is reproducible.
package cellular
