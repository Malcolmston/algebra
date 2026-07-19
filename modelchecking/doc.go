// Package modelchecking implements temporal-logic model checking over finite
// Kripke structures in pure Go using only the standard library.
//
// The package covers the classical automata-theoretic and fixpoint approaches
// to verifying finite-state systems against linear-time and branching-time
// temporal specifications:
//
//   - [Kripke] is a finite labelled transition system over atomic propositions.
//     [StateSet] is a packed-bitset set of states supporting the Boolean
//     operations that drive the symbolic fixpoint computations.
//
//   - Linear temporal logic formulas are represented by the [LTL] syntax tree,
//     parsed from text by [ParseLTL], normalised by [LTL.NNF] and rewritten by
//     [LTL.Simplify]. Computation tree logic formulas use the [CTL] tree, parsed
//     by [ParseCTL] and reducible to an existential fragment by
//     [CTL.ExistentialNormalForm].
//
//   - CTL model checking is performed by the labelling algorithm [CTLCheck],
//     built from the predecessor images [PreExists] and [PreForall] and the
//     least/greatest fixpoint operators [SatEU], [SatEG] and their derived and
//     universal companions [SatAF], [SatAG], [SatAU], [SatAR] and [SatER].
//
//   - Linear temporal logic is handled automata-theoretically. [LTLToGenBuchi]
//     builds a generalized Büchi automaton from an LTL formula by the on-the-fly
//     tableau construction of Gerth, Peled, Vardi and Wolper; [Degeneralize]
//     converts it to an ordinary [Buchi] automaton via the counter construction.
//     Language emptiness is decided by [Buchi.IsEmpty] and [Buchi.IsEmptySCC]
//     with an accepting run recovered by [Buchi.AcceptingLasso].
//
//   - [LTLModelCheck] verifies a Kripke structure against an LTL property by
//     forming the product [ProductBuchiKripke] with the Büchi automaton of the
//     negated formula and testing emptiness, returning a [Counterexample] lasso
//     when the property fails. [LTLSatisfiable], [LTLValid] and [LTLEquivalent]
//     answer the corresponding decision problems for the logic itself.
//
//   - [Unroll] and the BMC helpers [BMCInvariant], [BoundedReachable],
//     [BMCFindLasso] and [BMCExistsGlobally] perform bounded model checking by
//     explicit depth-bounded unrolling of the transition relation.
//
//   - Fairness constraints are captured by [Fairness]; fair branching-time model
//     checking is provided by [FairCTLCheck] together with the fair operators
//     [FairEG], [FairEX], [FairEU] and [FairStates].
//
//   - Counterexamples and witnesses for individual operators are extracted by
//     [AGCounterexample], [AFCounterexample], [EFWitness], [EGWitness] and
//     [EUWitness].
//
// All algorithms are exact and operate on explicitly represented finite
// structures. The package performs no input/output and uses no concurrency.
package modelchecking
