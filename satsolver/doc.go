// Package satsolver implements Boolean logic and Boolean satisfiability using
// only the Go standard library.
//
// The package covers the classical toolbox of propositional reasoning and
// automated Boolean decision procedures:
//
//   - Literals and clauses (literal.go): the [Lit] type encoding a signed
//     variable reference, together with the [Clause] conjunctive/disjunctive
//     term type and its manipulation helpers.
//   - Conjunctive and disjunctive normal-form containers (cnf.go, dnf.go): the
//     [CNF] and [DNF] formula types with evaluation, simplification, DIMACS
//     serialisation ([CNF.DIMACS], [ParseDIMACS]) and unit/pure-literal
//     analysis.
//   - Boolean expression trees (expr.go): the sealed [Expr] interface with the
//     [Variable], [BoolConst], [Not], [And], [Or], [Xor], [Implies], [Iff],
//     [Nand], [Nor] and [Xnor] node types, the [Parse] recursive-descent
//     parser and structural utilities such as [Substitute] and [ExprEqual].
//   - Normal-form conversion (nnf.go, normalform.go): [ToNNF], [ToCNFExpr],
//     [ToDNFExpr], the Boolean [Simplify] rewriter and the [ToCNFFormula]
//     bridge into the clausal world.
//   - Tseitin encoding (tseitin.go): the [Tseitin] equisatisfiable
//     transformation and the incremental [TseitinEncoder].
//   - Truth tables and semantics (truthtable.go, semantic.go): [TruthTable],
//     [Minterms], [Maxterms] and the decision procedures [IsTautology],
//     [IsSatisfiable], [Equivalent], [Entails], [FindModel] and [CountModels].
//   - A DPLL/CDCL-lite SAT solver (solver.go): [DPLL], [SolveCNF],
//     [UnitPropagate], [PureLiteralElimination], resolution ([Resolve]) and
//     model enumeration ([AllSolutions], [CountSolutions]).
//   - Reduced ordered binary decision diagrams (bdd.go): the [BDD] manager with
//     [BDD.Apply], [BDD.Ite], [BDD.SatCount] and [BDD.FromExpr].
//   - Two-level minimisation (qm.go): the Quine-McCluskey engine
//     [QuineMcCluskey] with [PrimeImplicants] and [EssentialPrimeImplicants],
//     the [MinimizeSOP] cover routine and the [KarnaughMap] helper.
//
// Everything is implemented with real, self-contained algorithms and depends
// only on the packages math, math/big, sort, errors, fmt and strings.
package satsolver
