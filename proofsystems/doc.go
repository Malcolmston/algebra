// Package proofsystems implements the classical machinery of formal logic proof
// systems using only the Go standard library.
//
// It provides a single term and formula representation shared across a family
// of decision and proof procedures for propositional and first-order logic:
//
//   - First-order terms and unification (term.go, subst.go, unify.go): the
//     [Term] type with variable, constant and function constructors, the
//     [Substitution] type with composition, and Robinson's [UnifyTerms] /
//     [MostGeneralUnifier] with the occurs-check, plus one-way [MatchTerm].
//   - Formulas and parsing (formula.go, parse.go): the [Formula] tree over the
//     usual connectives and quantifiers, structural utilities, and the
//     recursive-descent [ParseFormula] / [ParseTerm] readers.
//   - Propositional semantics (eval.go, semantics.go): [Eval], [NewTruthTable],
//     and the decision procedures [IsTautology], [IsSatisfiable], [Entails],
//     [Equivalent], [FindModel] and [CountModels].
//   - Normal forms (cnf.go, tseitin.go): [ToNNF], [ToCNF], [ToDNFFormula], the
//     clausal [PCNF] type, and the equisatisfiable [Tseitin] transformation.
//   - Refutation procedures (resolution.go, dpll.go, tableau.go): propositional
//     [ResolutionRefute], the [DPLL] satisfiability solver with
//     [UnitPropagate] and [PureLiteralElimination], and the analytic
//     [TableauClosed] method.
//   - First-order clauses and resolution (foclause.go, foresolution.go):
//     Skolemising [Clausify], the [FOClause] type, binary [FOResolve] with
//     [FOFactor], and the bounded [FOResolutionRefute] / [FOEntails] search.
//   - Horn-clause logic programming (sld.go): the [Program] type with SLD
//     resolution via [Program.Solve], [Program.Query] and [SLDResolve], a tiny
//     Prolog-style engine.
//   - Proof checking (naturaldeduction.go, sequent.go): the natural-deduction
//     [Derivation] with [CheckND] and [NDProves], and the sequent-calculus
//     [SequentDerivation] with [CheckSequent], the [ProveSequent] backward
//     prover and [ValidSequent].
//
// Every algorithm is a real, self-contained implementation depending only on
// the packages errors, fmt, sort, strings and math from the standard library.
package proofsystems
