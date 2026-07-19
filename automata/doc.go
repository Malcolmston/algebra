// Package automata implements finite automata, regular expressions and related
// machines of formal-language theory in pure Go using only the standard
// library.
//
// The package is organised around a small number of concrete machine types and
// a large collection of algorithms that transform and analyse them:
//
//   - [DFA] is a deterministic finite automaton with a possibly partial
//     transition function. It supports acceptance testing, completion,
//     complementation, trimming of unreachable and dead states, shortest-word
//     and word-enumeration queries, exact word counting via [DFA.CountAcceptedWords],
//     and emptiness/finiteness/universality decisions.
//
//   - [NFA] is a nondeterministic finite automaton with epsilon transitions. It
//     supports epsilon-closure computation ([NFA.EpsilonClosure]), simulation,
//     epsilon elimination, reversal and conversion to a [DFA] by the powerset
//     construction ([SubsetConstruction]).
//
//   - Regular expressions are parsed by [ParseRegex] into a [RegexNode] syntax
//     tree and compiled to automata by Thompson's construction ([Thompson]).
//     [Compile] returns a DFA-backed [Regexp] for fast repeated matching.
//
//   - DFA minimisation is available through both Hopcroft's ([Hopcroft]) and
//     Moore's ([Moore]) partition-refinement algorithms, which serve as mutual
//     cross-checks and expose the Myhill–Nerode equivalence classes.
//
//   - Boolean combinations of regular languages are built with the product
//     construction: [Union], [Intersection], [Difference] and
//     [SymmetricDifference]. Language equivalence ([Equivalent]), containment
//     ([Subset]) and disjointness ([Disjoint]) are decided by reducing to
//     emptiness, and [Witness] returns a shortest distinguishing string.
//
//   - Language operations on NFAs include [Concatenate], [UnionNFA], [StarNFA],
//     [PlusNFA], [OptionalNFA], [PowerNFA] and [ReverseNFA].
//
//   - The pumping lemma is made constructive by [PumpingLength] and
//     [PumpingDecomposition], which return a valid pumping length and an
//     explicit x·y·z factorisation of any sufficiently long accepted word.
//
//   - Beyond regular languages, [PDA] simulates nondeterministic pushdown
//     automata (accepting by final state or by empty stack) and [TM] simulates
//     deterministic single-tape Turing machines, both with bounded-resource
//     guarantees against non-termination.
//
// States are identified by the integers 0..NumStates-1; input, tape and stack
// symbols are runes, so the machines operate naturally on Go strings. All
// constructions are deterministic and free of global mutable state, and no
// algorithm depends on any package outside the Go standard library.
package automata
