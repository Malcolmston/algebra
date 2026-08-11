# proofsystems — implements the classical machinery of formal logic proof systems using only the Go

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/proofsystems.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/proofsystems)

Package proofsystems implements the classical machinery of formal logic proof
systems using only the Go standard library.

- First-order terms and unification (term.go, subst.go, unify.go): the `Term`
type with variable, constant and function constructors, the `Substitution`
type with composition, and Robinson's `UnifyTerms` / `MostGeneralUnifier` with
the occurs-check, plus one-way `MatchTerm`. - Formulas and parsing
(formula.go, parse.go): the `Formula` tree over the usual connectives and
quantifiers, structural utilities, and the recursive-descent `ParseFormula` /
`ParseTerm` readers. - Propositional semantics (eval.go, semantics.go):
`Eval`, `NewTruthTable`, and the decision procedures `IsTautology`,
`IsSatisfiable`, `Entails`, `Equivalent`, `FindModel` and `CountModels`. -
Normal forms (cnf.go, tseitin.go): `ToNNF`, `ToCNF`, `ToDNFFormula`, the
clausal `PCNF` type, and the equisatisfiable `Tseitin` transformation. -
Refutation procedures (resolution.go, dpll.go, tableau.go): propositional
`ResolutionRefute`, the `DPLL` satisfiability solver with `UnitPropagate` and
`PureLiteralElimination`, and the analytic `TableauClosed` method. -
First-order clauses and resolution (foclause.go, foresolution.go): Skolemising
`Clausify`, the `FOClause` type, binary `FOResolve` with `FOFactor`, and the
bounded `FOResolutionRefute` / `FOEntails` search. - Horn-clause logic
programming (sld.go): the `Program` type with SLD resolution via
`Program.Solve`, `Program.Query` and `SLDResolve`, a tiny Prolog-style engine.
- Proof checking (naturaldeduction.go, sequent.go): the natural-deduction
`Derivation` with `CheckND` and `NDProves`, and the sequent-calculus
`SequentDerivation` with `CheckSequent`, the `ProveSequent` backward prover
and `ValidSequent`.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/proofsystems"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func CheckSequent(d SequentDerivation) error` | CheckSequent verifies that a sequent derivation is a valid LK proof: every leaf is an axiom and every internal node is a correct instance of its rule. |
| `func CollectVars(terms ...Term) []string` | CollectVars returns the sorted union of variable names occurring in all of the given terms. |
| `func Consistent(fs []Formula) bool` | Consistent reports whether a set of formulas is jointly satisfiable. |
| `func CountModels(f Formula) int` | CountModels returns the number of satisfying assignments of the formula. |
| `func DPLLSatisfiable(n PCNF) bool` | DPLLSatisfiable reports whether the clause set is satisfiable. |
| `func DPLLSatisfiableFormula(f Formula) bool` | DPLLSatisfiableFormula converts a propositional formula to CNF via the Tseitin encoding and decides its satisfiability with DPLL. |
| `func Entails(premises []Formula, conclusion Formula) bool` | Entails reports whether the set of premises semantically entails the conclusion: every assignment satisfying all premises also satisfies the… |
| `func EntailsSingle(premise, conclusion Formula) bool` | EntailsSingle reports whether a single premise entails the conclusion. |
| `func Equivalent(f, g Formula) bool` | Equivalent reports whether two propositional formulas have the same truth value under every assignment of their combined atoms. |
| `func Eval(f Formula, a Assignment) (bool, error)` | Eval evaluates a propositional formula under the given assignment, treating each nullary atom as a Boolean variable. |
| `func FOEntails(premises []Formula, conclusion Formula, maxClauses int) bool` | FOEntails reports whether the first-order premises entail the conclusion, searched by refuting the clause set of the premises together with the… |
| `func FOResolutionRefute(clauses []FOClause, maxClauses int) bool` | FOResolutionRefute attempts to derive the empty clause from a first-order clause set by saturation under binary resolution and factoring. |
| `func FOUnsatisfiable(clauses []FOClause, maxClauses int) bool` | FOUnsatisfiable reports whether a first-order clause set is refutable within the given bound. |
| `func FormulasEqual(a, b []Formula) bool` | FormulasEqual reports whether two formula slices are pointwise structurally equal. |
| `func Inconsistent(fs []Formula) bool` | Inconsistent reports whether a set of formulas is jointly unsatisfiable. |
| `func IsContingent(f Formula) bool` | IsContingent reports whether the formula is satisfiable but not a tautology, i.e. |
| `func IsContradiction(f Formula) bool` | IsContradiction reports whether the formula is false under every assignment. |
| `func IsHornCNF(n PCNF) bool` | IsHornCNF reports whether every clause of the CNF is a Horn clause. |
| `func IsHornClause(c PClause) bool` | IsHornClause reports whether a propositional clause is a Horn clause: it has at most one positive literal. |
| `func IsSatisfiable(f Formula) bool` | IsSatisfiable reports whether some assignment makes the formula true. |
| `func IsTautology(f Formula) bool` | IsTautology reports whether the propositional formula is true under every assignment of its atoms. |
| `func IsUnifiable(a, b Term) bool` | IsUnifiable reports whether two terms have a unifier. |
| `func IsUnsatisfiable(f Formula) bool` | IsUnsatisfiable is a synonym for IsContradiction. |
| `func IsValid(f Formula) bool` | IsValid is a synonym for IsTautology. |
| `func Matches(pattern, subject Term) bool` | Matches reports whether pattern one-way matches subject. |
| `func Models(a Assignment, f Formula) bool` | Models reports whether the assignment satisfies the formula. |
| `func MustEval(f Formula, a Assignment) bool` | MustEval evaluates f under a and panics on error. |
| `func NDProves(d Derivation) bool` | NDProves reports whether the derivation is a closed proof (no open assumptions) of its conclusion, i.e. |
| `func NDValid(d Derivation, premises []Formula) bool` | NDValid reports whether the derivation is well formed and its open assumptions are all contained in the given premise set, so it establishes premises… |
| `func PureLiteralElimination(n PCNF) (PCNF, Assignment)` | PureLiteralElimination repeatedly removes clauses satisfied by pure literals, returning the reduced clause set and the assignment fixing the pure… |
| `func ResolutionEntails(premises []Formula, conclusion Formula) bool` | ResolutionEntails reports whether the premises entail the conclusion by refuting the clause set of the premises together with the negated conclusion. |
| `func ResolutionRefute(n PCNF) bool` | ResolutionRefute attempts to derive the empty clause from a clause set by exhaustive propositional resolution. |
| `func ResolutionUnsatisfiable(n PCNF) bool` | ResolutionUnsatisfiable reports whether a CNF is unsatisfiable via resolution. |
| `func ResolutionValid(f Formula) bool` | ResolutionValid reports whether a propositional formula is valid by refuting its negation with resolution: f is valid iff ¬f is unsatisfiable. |
| `func SequentProvable(s Sequent) bool` | SequentProvable reports whether ProveSequent finds a proof of the sequent. |
| `func SkolemName(k int) string` | SkolemName returns the canonical name of the k-th Skolem symbol, matching the symbols introduced by Clausify. |
| `func TableauClosed(formulas []Formula) bool` | TableauClosed reports whether the analytic (semantic) tableau for the given set of propositional formulas closes on every branch, which holds exactly… |
| `func TableauEntails(premises []Formula, conclusion Formula) bool` | TableauEntails reports whether the premises entail the conclusion by closing the tableau for the premises together with the negated conclusion. |
| `func TableauProve(f Formula) bool` | TableauProve is a synonym for TableauValid: it reports whether f is a propositional tautology. |
| `func TableauSatisfiable(formulas []Formula) bool` | TableauSatisfiable reports whether a set of propositional formulas is satisfiable according to the tableau method (some branch stays open). |
| `func TableauValid(f Formula) bool` | TableauValid reports whether a propositional formula is valid by checking that the tableau for its negation closes. |
| `func TermsEqual(a, b []Term) bool` | TermsEqual reports whether two slices of terms are pointwise structurally equal. |
| `func UnitPropagate(n PCNF) (PCNF, Assignment, bool)` | UnitPropagate applies a single fixpoint of unit propagation to a clause set, returning the simplified clauses, the forced literal assignment, and… |
| `func ValidSequent(s Sequent) bool` | ValidSequent reports whether a propositional sequent is valid, i.e. |

### Types

| Type | What it is |
| --- | --- |
| `Assignment` | Assignment maps propositional atom names to Boolean truth values. |
| `Connective` | Connective identifies the top-level logical operator of a formula. |
| `DPLLResult` | DPLLResult reports the outcome of a DPLL search: whether the clause set is satisfiable and, when it is, a satisfying assignment. |
| `Derivation` | Derivation is a node of a natural-deduction proof tree. |
| `FOClause` | FOClause is a first-order clause: a disjunction of first-order literals with every variable implicitly universally quantified. |
| `FOLiteral` | FOLiteral is a first-order literal: a possibly negated atomic predicate applied to argument terms. |
| `Formula` | Formula is an immutable propositional or first-order formula tree. |
| `HornClause` | HornClause is a definite Horn clause head :- b1, ..., bn, where the head and body atoms are positive first-order literals. |
| `NDRule` | NDRule identifies a natural-deduction inference rule in the Gentzen/Prawitz intuitionistic-plus-classical propositional calculus implemented here. |
| `PCNF` | PCNF is a propositional formula in conjunctive normal form: a conjunction of clauses. |
| `PClause` | PClause is a propositional clause: a disjunction of literals. |
| `PLiteral` | PLiteral is a propositional literal: a Boolean variable name together with a polarity. |
| `Program` | Program is an ordered collection of Horn clauses forming a tiny Prolog-style logic program. |
| `SeqRule` | SeqRule identifies a rule of the propositional sequent calculus LK. |
| `Sequent` | Sequent is a two-sided propositional sequent Γ ⊢ Δ where Left is the antecedent (assumed formulas) and Right is the succedent (goal alternatives). |
| `SequentDerivation` | SequentDerivation is a node of a sequent-calculus proof tree: an inference of Concl from the conclusions of its premise sub-derivations by Rule. |
| `SolveOptions` | SolveOptions bounds an SLD search. |
| `Substitution` | Substitution maps variable names to replacement terms. |
| `Term` | Term is an immutable first-order term. |
| `TermKind` | TermKind classifies the three syntactic categories of a first-order term: logical variables, constant symbols, and compound function applications. |
| `TruthRow` | TruthRow is a single line of a truth table: the assignment and the resulting value of the formula. |
| `TruthTable` | TruthTable holds the complete truth table of a propositional formula. |
| `TseitinResult` | TseitinResult is the outcome of a Tseitin transformation: the equisatisfiable clause set together with the name of the auxiliary variable… |

<details>
<summary><code>Assignment</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func AllAssignments(atoms []string) []Assignment` | AllAssignments enumerates every total truth assignment over the given atom names in ascending binary order (the first atom is the most significant… |
| `func AllModels(f Formula) []Assignment` | AllModels returns every satisfying assignment of the formula over its atoms. |
| `func FindModel(f Formula) (Assignment, bool)` | FindModel returns a satisfying assignment for the formula and true, or the zero assignment and false when the formula is unsatisfiable. |
| `func (a Assignment) Clone() Assignment` | Clone returns an independent copy of the assignment. |
| `func (a Assignment) String() string` | String renders an assignment as {A=true, B=false} with names sorted. |

</details>

<details>
<summary><code>Connective</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (c Connective) String() string` | String returns a short mnemonic for the connective. |

</details>

<details>
<summary><code>DPLLResult</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func DPLL(n PCNF) DPLLResult` | DPLL decides the satisfiability of a propositional clause set using the Davis–Putnam–Logemann–Loveland procedure with unit propagation and… |

</details>

<details>
<summary><code>Derivation</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func AndElim1(d Derivation) Derivation` | AndElim1 builds a left conjunction-elimination node concluding the left conjunct of the premise. |
| `func AndElim2(d Derivation) Derivation` | AndElim2 builds a right conjunction-elimination node concluding the right conjunct of the premise. |
| `func AndIntro(left, right Derivation) Derivation` | AndIntro builds a conjunction-introduction node from proofs of the conjuncts. |
| `func Assume(f Formula) Derivation` | Assume constructs an assumption leaf whose conclusion is f. |
| `func ImpElim(a, imp Derivation) Derivation` | ImpElim builds a modus-ponens node from proofs of A and A → B, concluding B. |
| `func ImpIntro(assumption Formula, body Derivation) Derivation` | ImpIntro builds an implication-introduction node concluding assumption → body.Concl and discharging assumption from body. |
| `func NotElim(a, notA Derivation) Derivation` | NotElim builds a negation-elimination node concluding ⊥ from A and ¬A. |
| `func NotIntro(assumption Formula, bottom Derivation) Derivation` | NotIntro builds a negation-introduction node concluding ¬assumption from a derivation of ⊥, discharging assumption. |
| `func OrIntro1(d Derivation, other Formula) Derivation` | OrIntro1 builds a left disjunction-introduction node concluding d.Concl ∨ other. |
| `func OrIntro2(other Formula, d Derivation) Derivation` | OrIntro2 builds a right disjunction-introduction node concluding other ∨ d.Concl. |
| `func RAA(goal Formula, bottom Derivation) Derivation` | RAA builds a reductio node concluding goal from a derivation of ⊥, discharging the assumption ¬goal. |

</details>

<details>
<summary><code>FOClause</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Clausify(f Formula) []FOClause` | Clausify converts a first-order formula into an equisatisfiable set of clauses by NNF conversion, standardising variables apart, Skolemisation,… |
| `func ClausifyAll(fs ...Formula) []FOClause` | ClausifyAll clausifies a list of formulas into a single flat clause set. |
| `func FOFactor(c FOClause) []FOClause` | FOFactor computes the factors of a clause: for each pair of like-signed literals that unify, the clause with the unifier applied and the merged… |
| `func FOResolve(a, b FOClause) []FOClause` | FOResolve computes the binary resolvents of two first-order clauses. |
| `func NewFOClause(lits ...FOLiteral) FOClause` | NewFOClause builds a clause from the given literals. |
| `func (c FOClause) Canonical() FOClause` | Canonical returns a clause with literals de-duplicated and sorted by string form, giving a stable representative for set membership. |
| `func (c FOClause) IsEmpty() bool` | IsEmpty reports whether the clause has no literals. |
| `func (c FOClause) RenameApart(suffix string) FOClause` | RenameApart returns a copy of the clause with every variable renamed using the given suffix, so that two clauses can be resolved without variable… |
| `func (c FOClause) String() string` | String renders the clause as a bracketed disjunction, or [] for the empty clause. |
| `func (c FOClause) Vars() []string` | Vars returns the sorted variable names occurring in the clause. |

</details>

<details>
<summary><code>FOLiteral</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewFOLiteral(pred string, args ...Term) FOLiteral` | NewFOLiteral builds a positive first-order literal P(args...). |
| `func NewNegFOLiteral(pred string, args ...Term) FOLiteral` | NewNegFOLiteral builds a negative first-order literal ¬P(args...). |
| `func (l FOLiteral) Equal(o FOLiteral) bool` | Equal reports whether two first-order literals are structurally identical. |
| `func (l FOLiteral) IsComplementary(o FOLiteral) bool` | IsComplementary reports whether two literals share a predicate and arity but differ in sign, so they are candidates for resolution once their… |
| `func (l FOLiteral) Negated() FOLiteral` | Negated returns the complementary literal. |
| `func (l FOLiteral) String() string` | String renders the literal in P(a,b) or !P(a,b) form. |
| `func (l FOLiteral) Vars() []string` | Vars returns the sorted variable names occurring in the literal. |

</details>

<details>
<summary><code>Formula</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func And(a, b Formula) Formula` | And returns the conjunction of a and b. |
| `func Atom(pred string, args ...Term) Formula` | Atom builds a first-order atomic formula P(t1,...,tn). |
| `func Bot() Formula` | Bot returns the logical constant falsehood. |
| `func CheckND(d Derivation) ([]Formula, error)` | CheckND verifies a natural-deduction derivation and, on success, returns the multiset of undischarged (open) assumptions it depends on. |
| `func Conj(fs ...Formula) Formula` | Conj folds a slice of formulas into a right-nested conjunction. |
| `func Disj(fs ...Formula) Formula` | Disj folds a slice of formulas into a right-nested disjunction. |
| `func Exists(v string, body Formula) Formula` | Exists returns the existentially quantified formula ∃v. |
| `func FONNF(f Formula) Formula` | FONNF converts a first-order formula to negation normal form, eliminating implications and biconditionals and pushing negations through the… |
| `func Forall(v string, body Formula) Formula` | Forall returns the universally quantified formula ∀v. |
| `func Iff(a, b Formula) Formula` | Iff returns the biconditional a ↔ b. |
| `func Imp(a, b Formula) Formula` | Imp returns the implication a → b. |
| `func MustParseFormula(src string) Formula` | MustParseFormula parses src and panics on error. |
| `func Not(f Formula) Formula` | Not returns the negation of f. |
| `func Or(a, b Formula) Formula` | Or returns the disjunction of a and b. |
| `func ParseFormula(src string) (Formula, error)` | ParseFormula parses a propositional or first-order formula from its textual form. |
| `func Prop(name string) Formula` | Prop builds a propositional atom (a nullary predicate) with the given name. |
| `func SortFormulas(fs []Formula) []Formula` | SortFormulas returns a copy of the slice sorted by string form, giving a canonical ordering useful for set comparisons. |
| `func ToCNFFormula(f Formula) Formula` | ToCNFFormula converts a propositional formula into an equivalent formula in conjunctive normal form by NNF conversion followed by distribution of… |
| `func ToDNFFormula(f Formula) Formula` | ToDNFFormula converts a propositional formula into an equivalent formula in disjunctive normal form. |
| `func ToNNF(f Formula) Formula` | ToNNF converts a propositional formula to negation normal form: negations are pushed inward until they apply only to atoms, and implications and… |
| `func Top() Formula` | Top returns the logical constant truth. |
| `func (f Formula) Atoms() []string` | Atoms returns the sorted, de-duplicated list of atom string forms occurring in the formula. |
| `func (f Formula) Depth() int` | Depth returns the height of the formula tree; atoms and constants have depth 0. |
| `func (f Formula) Equal(g Formula) bool` | Equal reports whether two formulas are structurally identical (up to the literal names of bound variables — no alpha-conversion is performed). |
| `func (f Formula) FreeVars() []string` | FreeVars returns the sorted first-order variables occurring free in the formula (not bound by an enclosing quantifier). |
| `func (f Formula) IsAtom() bool` | IsAtom reports whether the formula is an atomic formula or propositional variable. |
| `func (f Formula) IsConstant() bool` | IsConstant reports whether the formula is one of the logical constants ⊤ or ⊥. |
| `func (f Formula) IsLiteral() bool` | IsLiteral reports whether the formula is an atom or the negation of an atom. |
| `func (f Formula) IsQuantified() bool` | IsQuantified reports whether the top connective is a quantifier. |
| `func (f Formula) PropVars() []string` | PropVars returns the sorted names of the propositional variables (nullary atoms) in the formula, ignoring first-order atoms with arguments. |
| `func (f Formula) Size() int` | Size returns the number of connective and atom nodes in the formula tree. |
| `func (f Formula) String() string` | String renders the formula with minimal parenthesisation using ASCII operators (! & \| -> <-> for the connectives and forall/exists for quantifiers). |
| `func (f Formula) Subformulas() []Formula` | Subformulas returns every subformula including f itself, in pre-order, with structural duplicates removed. |
| `func (f Formula) SubstituteProp(p string, r Formula) Formula` | SubstituteProp returns the formula with every propositional atom named p replaced by the formula r. |
| `func (f Formula) SubstituteTerm(v string, t Term) Formula` | SubstituteTerm returns the formula with every free occurrence of the first-order variable v replaced by the term t. |

</details>

<details>
<summary><code>HornClause</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Fact(head FOLiteral) HornClause` | Fact builds a Horn clause with the given head and no body. |
| `func Rule(head FOLiteral, body ...FOLiteral) HornClause` | Rule builds a Horn clause head :- body. |
| `func (h HornClause) IsFact() bool` | IsFact reports whether the clause has an empty body. |
| `func (h HornClause) String() string` | String renders the clause in Prolog-like syntax. |

</details>

<details>
<summary><code>NDRule</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (r NDRule) String() string` | String returns the conventional name of a natural-deduction rule. |

</details>

<details>
<summary><code>PCNF</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewPCNF(clauses ...PClause) PCNF` | NewPCNF builds a CNF from the given clauses. |
| `func ToCNF(f Formula) PCNF` | ToCNF converts a propositional formula to the clausal PCNF representation via full distribution. |
| `func TseitinCNF(f Formula, prefix string) PCNF` | TseitinCNF returns the equisatisfiable clause set of f with the top auxiliary variable asserted as a unit clause, so the returned CNF is satisfiable… |
| `func (n PCNF) Eval(a Assignment) bool` | Eval returns the truth value of the CNF (a conjunction of clauses) under the assignment. |
| `func (n PCNF) Simplify() PCNF` | Simplify removes tautological clauses and clauses subsumed by another, returning a logically equivalent but smaller CNF. |
| `func (n PCNF) String() string` | String renders the CNF as a conjunction of clause strings. |
| `func (n PCNF) Vars() []string` | Vars returns the sorted variable names occurring in the CNF. |

</details>

<details>
<summary><code>PClause</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func CNFClauses(f Formula) []PClause` | CNFClauses returns the list of clauses of the clausal form of f. |
| `func NewPClause(lits ...PLiteral) PClause` | NewPClause builds a clause from the given literals, removing exact duplicates and sorting for a canonical form. |
| `func Resolve(a, b PClause) []PClause` | Resolve computes all resolvents of two propositional clauses. |
| `func (c PClause) Contains(l PLiteral) bool` | Contains reports whether the clause contains the given literal. |
| `func (c PClause) Equal(o PClause) bool` | Equal reports whether two clauses contain the same set of literals. |
| `func (c PClause) Eval(a Assignment) bool` | Eval returns the truth value of the clause (a disjunction) under the assignment. |
| `func (c PClause) IsEmpty() bool` | IsEmpty reports whether the clause has no literals (the empty clause ⊥). |
| `func (c PClause) IsTautologyClause() bool` | IsTautologyClause reports whether the clause contains a variable and its negation and is therefore always true. |
| `func (c PClause) IsUnit() bool` | IsUnit reports whether the clause contains exactly one literal. |
| `func (c PClause) String() string` | String renders the clause as a bracketed disjunction, or [] for the empty clause. |
| `func (c PClause) Subsumes(o PClause) bool` | Subsumes reports whether clause c subsumes o, i.e. |
| `func (c PClause) Vars() []string` | Vars returns the sorted variable names occurring in the clause. |

</details>

<details>
<summary><code>PLiteral</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NegPLit(v string) PLiteral` | NegPLit returns the negative literal for the given variable. |
| `func PosPLit(v string) PLiteral` | PosPLit returns the positive literal for the given variable. |
| `func (l PLiteral) Equal(o PLiteral) bool` | Equal reports whether two literals are identical. |
| `func (l PLiteral) Eval(a Assignment) bool` | Eval returns the truth value of the literal under the assignment. |
| `func (l PLiteral) IsComplementary(o PLiteral) bool` | IsComplementary reports whether two literals are on the same variable with opposite polarity. |
| `func (l PLiteral) Negated() PLiteral` | Negated returns the complementary literal. |
| `func (l PLiteral) String() string` | String renders the literal as A or !A. |

</details>

<details>
<summary><code>Program</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewProgram(clauses ...HornClause) *Program` | NewProgram builds a program from the given clauses. |
| `func (p *Program) Add(c HornClause)` | Add appends a clause to the program. |
| `func (p *Program) AddFact(head FOLiteral)` | AddFact appends a fact to the program. |
| `func (p *Program) AddRule(head FOLiteral, body ...FOLiteral)` | AddRule appends a rule to the program. |
| `func (p *Program) Prove(goals ...FOLiteral) bool` | Prove is a synonym for Query used when the goals encode a yes/no question. |
| `func (p *Program) Query(goals ...FOLiteral) bool` | Query reports whether the goal list is derivable from the program (there is at least one SLD refutation) within the default search bounds. |
| `func (p *Program) Solve(goals []FOLiteral, opts SolveOptions) []Substitution` | Solve performs SLD resolution of the goal list against the program and returns the answer substitutions restricted to the variables of the query. |
| `func (p *Program) SolveOne(goals ...FOLiteral) (Substitution, bool)` | SolveOne returns the first answer substitution for the goal list and whether one was found. |
| `func (p *Program) String() string` | String renders the whole program, one clause per line. |

</details>

<details>
<summary><code>SeqRule</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (r SeqRule) String() string` | String returns the conventional name of the sequent rule. |

</details>

<details>
<summary><code>Sequent</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewSequent(left, right []Formula) Sequent` | NewSequent builds a sequent from antecedent and succedent slices. |
| `func (s Sequent) IsAxiom() bool` | IsAxiom reports whether the sequent is an axiom instance, i.e. |
| `func (s Sequent) String() string` | String renders the sequent using the ⊢ turnstile with comma-separated sides. |

</details>

<details>
<summary><code>SequentDerivation</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func ProveFormula(f Formula) (SequentDerivation, bool)` | ProveFormula attempts to prove the propositional formula f as the sequent ⊢ f, returning a sequent derivation when f is valid. |
| `func ProveSequent(s Sequent) (SequentDerivation, bool)` | ProveSequent searches for an LK proof of the sequent using backward application of the invertible propositional rules. |

</details>

<details>
<summary><code>SolveOptions</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func DefaultSolveOptions() SolveOptions` | DefaultSolveOptions returns reasonable search bounds. |

</details>

<details>
<summary><code>Substitution</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func MatchTerm(pattern, subject Term) (Substitution, error)` | MatchTerm computes a one-way match (pattern matching): a substitution s with domain among the variables of pattern such that s.ApplyTerm(pattern)… |
| `func MostGeneralUnifier(a, b Term) (Substitution, error)` | MostGeneralUnifier is a synonym for UnifyTerms returning the MGU of two terms. |
| `func NewSubstitution() Substitution` | NewSubstitution returns the empty (identity) substitution. |
| `func SLDResolve(clauses []HornClause, goals []FOLiteral) []Substitution` | SLDResolve is a package-level convenience that builds a program from clauses and solves the goals with the default options. |
| `func SingletonSubstitution(v string, t Term) Substitution` | SingletonSubstitution returns the substitution {v -> t}. |
| `func UnifyLiterals(a, b FOLiteral) (Substitution, error)` | UnifyLiterals unifies two first-order literals of the same predicate and arity (their signs are ignored, so a positive and a negative literal on the… |
| `func UnifyTermLists(a, b []Term) (Substitution, error)` | UnifyTermLists unifies two equal-length lists of terms position by position, threading the accumulated substitution. |
| `func UnifyTerms(a, b Term) (Substitution, error)` | UnifyTerms computes a most-general unifier of two terms using Robinson's algorithm with the occurs-check. |
| `func (s Substitution) ApplyClause(c FOClause) FOClause` | ApplyClause applies the substitution to every literal of a first-order clause. |
| `func (s Substitution) ApplyLiteral(l FOLiteral) FOLiteral` | ApplyLiteral applies the substitution to every argument term of a first-order literal. |
| `func (s Substitution) ApplyTerm(t Term) Term` | ApplyTerm applies the substitution to a term, recursively replacing bound variables by their images. |
| `func (s Substitution) ApplyTerms(ts []Term) []Term` | ApplyTerms applies the substitution to each term in a slice. |
| `func (s Substitution) Bind(v string, t Term) Substitution` | Bind returns a copy of the substitution with v additionally mapped to t. |
| `func (s Substitution) Compose(o Substitution) Substitution` | Compose returns the composition s then o (written o∘s): applying the result to a term is the same as applying s and then o. |
| `func (s Substitution) Domain() []string` | Domain returns the sorted list of variable names bound by the substitution. |
| `func (s Substitution) Equal(o Substitution) bool` | Equal reports whether two substitutions have identical bindings. |
| `func (s Substitution) Get(v string) (Term, bool)` | Get returns the image of variable v and whether it is bound. |
| `func (s Substitution) IsEmpty() bool` | IsEmpty reports whether the substitution has no bindings (is the identity). |
| `func (s Substitution) Len() int` | Len returns the number of variable bindings. |
| `func (s Substitution) Restrict(vars []string) Substitution` | Restrict returns the substitution limited to the given set of variable names. |
| `func (s Substitution) String() string` | String renders a substitution as {x -> t, y -> u} with the domain sorted. |

</details>

<details>
<summary><code>Term</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func MustParseTerm(src string) Term` | MustParseTerm parses src and panics on error. |
| `func NewConst(name string) Term` | NewConst builds a constant term (a nullary function symbol) with the given name. |
| `func NewFunc(name string, args ...Term) Term` | NewFunc builds a compound term applying the function symbol name to the given argument terms. |
| `func NewVar(name string) Term` | NewVar builds a first-order variable term with the given name. |
| `func ParseTerm(src string) (Term, error)` | ParseTerm parses a single first-order term from text, such as "f(x,a)". |
| `func (t Term) Arity() int` | Arity returns the number of argument sub-terms; it is zero for variables and constants. |
| `func (t Term) Depth() int` | Depth returns the nesting depth of the term: 0 for a variable or constant, and 1 plus the maximum argument depth for a function application. |
| `func (t Term) Equal(u Term) bool` | Equal reports whether two terms are structurally identical. |
| `func (t Term) IsConst() bool` | IsConst reports whether the term is a constant symbol. |
| `func (t Term) IsFunc() bool` | IsFunc reports whether the term is a compound function application. |
| `func (t Term) IsGround() bool` | IsGround reports whether the term contains no variables. |
| `func (t Term) IsVar() bool` | IsVar reports whether the term is a variable. |
| `func (t Term) Occurs(v string) bool` | Occurs reports whether the variable named v occurs anywhere within the term. |
| `func (t Term) Rename(r map[string]string) Term` | Rename returns a copy of the term with every variable name x replaced by the name r(x) if r has an entry for it. |
| `func (t Term) Size() int` | Size returns the total number of symbol occurrences (nodes) in the term. |
| `func (t Term) String() string` | String renders a term in standard f(a,g(x)) notation. |
| `func (t Term) Subterms() []Term` | Subterms returns every subterm of t including t itself, in pre-order. |
| `func (t Term) Vars() []string` | Vars returns the sorted, de-duplicated list of variable names occurring in the term. |

</details>

<details>
<summary><code>TermKind</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (k TermKind) String() string` | String returns the human-readable name of a term kind. |

</details>

<details>
<summary><code>TruthTable</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewTruthTable(f Formula) (TruthTable, error)` | NewTruthTable builds the full truth table of a propositional formula over its atoms. |
| `func (tt TruthTable) CountTrue() int` | CountTrue returns the number of assignments under which the formula is true. |
| `func (tt TruthTable) String() string` | String renders the truth table as a fixed-width grid. |
| `func (tt TruthTable) TrueRows() []Assignment` | TrueRows returns the assignments under which the formula is true. |

</details>

<details>
<summary><code>TseitinResult</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Tseitin(f Formula, prefix string) TseitinResult` | Tseitin performs the Tseitin transformation of a propositional formula, introducing a fresh auxiliary variable for each compound subformula and… |

</details>

### Variables

`ErrEval`, `ErrParse`, `ErrUnify`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/proofsystems).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
