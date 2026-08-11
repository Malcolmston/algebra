# satsolver — implements Boolean logic and Boolean satisfiability using only the Go standard library

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/satsolver.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/satsolver)

Package satsolver implements Boolean logic and Boolean satisfiability using
only the Go standard library.

- Literals and clauses (literal.go): the `Lit` type encoding a signed variable
reference, together with the `Clause` conjunctive/disjunctive term type and
its manipulation helpers. - Conjunctive and disjunctive normal-form containers
(cnf.go, dnf.go): the `CNF` and `DNF` formula types with evaluation,
simplification, DIMACS serialisation (`CNF.DIMACS`, `ParseDIMACS`) and
unit/pure-literal analysis. - Boolean expression trees (expr.go): the sealed
`Expr` interface with the `Variable`, `BoolConst`, `Not`, `And`, `Or`, `Xor`,
`Implies`, `Iff`, `Nand`, `Nor` and `Xnor` node types, the `Parse`
recursive-descent parser and structural utilities such as `Substitute` and
`ExprEqual`. - Normal-form conversion (nnf.go, normalform.go): `ToNNF`,
`ToCNFExpr`, `ToDNFExpr`, the Boolean `Simplify` rewriter and the
`ToCNFFormula` bridge into the clausal world. - Tseitin encoding (tseitin.go):
the `Tseitin` equisatisfiable transformation and the incremental
`TseitinEncoder`. - Truth tables and semantics (truthtable.go, semantic.go):
`TruthTable`, `Minterms`, `Maxterms` and the decision procedures
`IsTautology`, `IsSatisfiable`, `Equivalent`, `Entails`, `FindModel` and
`CountModels`. - A DPLL/CDCL-lite SAT solver (solver.go): `DPLL`, `SolveCNF`,
`UnitPropagate`, `PureLiteralElimination`, resolution (`Resolve`) and model
enumeration (`AllSolutions`, `CountSolutions`). - Reduced ordered binary
decision diagrams (bdd.go): the `BDD` manager with `BDD.Apply`, `BDD.Ite`,
`BDD.SatCount` and `BDD.FromExpr`. - Two-level minimisation (qm.go): the
Quine-McCluskey engine `QuineMcCluskey` with `PrimeImplicants` and
`EssentialPrimeImplicants`, the `MinimizeSOP` cover routine and the
`KarnaughMap` helper.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/satsolver"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func AllModels(e Expr) []map[string]bool` | AllModels returns every satisfying assignment of e over its variables. |
| `func AllSolutions(f CNF) []map[int]bool` | AllSolutions returns every complete satisfying assignment of f over its variables, enumerated in ascending mask order. |
| `func BitString(x, numVars int) string` | BitString returns the binary representation of x in exactly numVars digits, most-significant bit first. |
| `func CountModels(e Expr) int` | CountModels returns the number of satisfying assignments of e over its variables. |
| `func CountNodesByKind(e Expr) map[string]int` | CountNodesByKind returns a map from a short kind name ("var", "const", "and", ...) to the number of nodes of that kind in e. |
| `func CountSolutions(f CNF) int` | CountSolutions returns the number of complete satisfying assignments of f over its variables. |
| `func DPLL(f CNF) (map[int]bool, bool)` | DPLL solves the satisfiability of a CNF formula with the classic Davis-Putnam-Logemann-Loveland procedure: unit propagation, pure-literal elimination… |
| `func Depth(e Expr) int` | Depth returns the height of the expression tree; a leaf has depth 1. |
| `func Entails(a, b Expr) bool` | Entails reports whether a semantically entails b (a \|= b): every assignment satisfying a also satisfies b. |
| `func EnvToIndex(vars []string, env map[string]bool) int` | EnvToIndex converts a variable environment to its row index for the given sorted variable list, first variable most significant. |
| `func Equivalent(a, b Expr) bool` | Equivalent reports whether two expressions have identical truth values on every assignment to their combined variables. |
| `func EquivalentCNF(f, g CNF) bool` | EquivalentCNF reports whether two clausal CNF formulas are logically equivalent by comparing their truth tables over the union of their variables. |
| `func Eval(e Expr, env map[string]bool) bool` | Eval evaluates e under env; a convenience wrapper over the `Expr.Eval` method. |
| `func EvalString(s string, env map[string]bool) (bool, error)` | EvalString parses s and evaluates it under env, returning the value and any parse error. |
| `func ExprEqual(a, b Expr) bool` | ExprEqual reports whether two expression trees are structurally identical. |
| `func FindModel(e Expr) (map[string]bool, bool)` | FindModel returns a satisfying assignment for e and true, or nil and false if e is unsatisfiable. |
| `func GrayCode(n int) []int` | GrayCode returns the reflected binary Gray-code sequence for n bits. |
| `func GrayDecode(g int) int` | GrayDecode inverts `GrayEncode`. |
| `func GrayEncode(n int) int` | GrayEncode returns the reflected binary Gray code of n, equal to n XOR (n>>1). |
| `func Implies2(a, b Expr) bool` | Implies reports the same relation as `Entails`; provided as a readable alias for semantic implication. |
| `func IndexToEnv(vars []string, i int) map[string]bool` | IndexToEnv converts a row index to a variable environment for the given sorted variable list, first variable most significant. |
| `func IsCNFExpr(e Expr) bool` | IsCNFExpr reports whether e is syntactically in conjunctive normal form: a conjunction of disjunctions of literals. |
| `func IsContingent(e Expr) bool` | IsContingent reports whether e is neither a tautology nor a contradiction: true for some assignments and false for others. |
| `func IsContradiction(e Expr) bool` | IsContradiction reports whether e evaluates to false under every assignment. |
| `func IsDNFExpr(e Expr) bool` | IsDNFExpr reports whether e is syntactically in disjunctive normal form: a disjunction of conjunctions of literals. |
| `func IsNNF(e Expr) bool` | IsNNF reports whether e is already in negation normal form: it contains only And, Or, Not, Variable and BoolConst nodes and every Not wraps a… |
| `func IsSatisfiable(e Expr) bool` | IsSatisfiable reports whether e evaluates to true under at least one assignment. |
| `func IsSatisfiableCNF(f CNF) bool` | IsSatisfiableCNF reports whether the CNF formula has at least one satisfying assignment. |
| `func IsTautology(e Expr) bool` | IsTautology reports whether e evaluates to true under every assignment. |
| `func IsTautologyCNF(f CNF) bool` | IsTautologyCNF reports whether every assignment satisfies f, i.e. |
| `func IsUnsatisfiableCNF(f CNF) bool` | IsUnsatisfiableCNF reports whether the CNF formula has no satisfying assignment. |
| `func Leaves(e Expr) []string` | Leaves returns the sorted distinct variable names, identical to `Vars`; kept for readability at call sites that view the tree as a term. |
| `func Maxterms(e Expr) []int` | Maxterms returns the maxterm indices of e: the assignments where e is false. |
| `func Minterms(e Expr) []int` | Minterms returns the set of minterm indices of e: the assignments (as integer row indices with the first sorted variable most significant) where e is… |
| `func NewBDDFromExpr(e Expr) (*BDD, Node)` | NewBDDFromExpr is a convenience that builds a manager whose ordering is the sorted variables of e and returns both the manager and the root node. |
| `func PopCount(x uint) int` | PopCount returns the number of one bits (the Hamming weight) of x. |
| `func PrettyString(e Expr) string` | PrettyString is an alias for the expression's String rendering, provided for symmetry with other renderers. |
| `func SOPString(cover []Implicant, names []string) string` | SOPString renders a list of implicants as a sum-of-products Boolean expression using the given variable names, most significant first. |
| `func Size(e Expr) int` | Size returns the number of nodes in the expression tree. |
| `func SolveCNF(f CNF) (map[int]bool, bool)` | SolveCNF solves f and, on success, returns a complete assignment covering every variable of f (variables left free by the search are set to false). |
| `func ToCNFFormula(e Expr) (CNF, *VarMap)` | ToCNFFormula converts an expression to a clausal `CNF` together with the `VarMap` recording the name-to-index correspondence. |
| `func ToDNFFormula(e Expr) (DNF, *VarMap)` | ToDNFFormula converts an expression to a clausal `DNF` together with the `VarMap` recording the name-to-index correspondence. |
| `func Tseitin(e Expr) (CNF, *VarMap)` | Tseitin returns an equisatisfiable CNF encoding of e that asserts e is true, together with the variable map used. |
| `func Vars(e Expr) []string` | Vars returns the sorted set of distinct variable names occurring in e. |

### Types

| Type | What it is |
| --- | --- |
| `And` | And is logical conjunction of two operands. |
| `BDD` | BDD is a reduced ordered binary decision diagram manager. |
| `BoolConst` | BoolConst is a Boolean literal constant, either true or false. |
| `CNF` | CNF is a formula in conjunctive normal form: a conjunction (AND) of clauses, each of which is a disjunction (OR) of literals. |
| `Clause` | Clause is a disjunction of literals (an OR term) when used inside a `CNF` and a conjunction of literals (an AND term) when used inside a `DNF`. |
| `DNF` | DNF is a formula in disjunctive normal form: a disjunction (OR) of terms, each of which is a conjunction (AND) of literals. |
| `Expr` | Expr is the sealed interface implemented by every node of a Boolean expression tree. |
| `Iff` | Iff is the biconditional X <-> Y (logical equivalence). |
| `Implicant` | Implicant is a Boolean cube: a product term covering a set of minterms. |
| `Implies` | Implies is the material conditional X -> Y. |
| `KarnaughMap` | KarnaughMap is a two-dimensional Gray-coded arrangement of the truth table of a Boolean function, the layout used for manual minimisation. |
| `Lit` | Lit is a Boolean literal: a signed reference to a variable. |
| `Nand` | Nand is the negated conjunction of two operands. |
| `Node` | Node identifies a vertex of a reduced ordered binary decision diagram managed by a `BDD`. |
| `Nor` | Nor is the negated disjunction of two operands. |
| `Not` | Not is logical negation of its operand. |
| `Or` | Or is logical inclusive disjunction of two operands. |
| `Solver` | Solver is a small stateful wrapper around the DPLL search that remembers the formula and its most recent model. |
| `TruthRow` | TruthRow is a single line of a `TruthTable`: an assignment to the variables and the resulting value of the expression. |
| `TruthTable` | TruthTable is the exhaustive truth table of a Boolean expression over its variables, ordered so that the first variable is the most significant bit… |
| `TseitinEncoder` | TseitinEncoder builds an equisatisfiable CNF encoding of a Boolean expression by introducing a fresh auxiliary variable for every compound… |
| `VarMap` | VarMap maps variable names to positive integer indices (1-based), the numeric encoding used by the clausal solver. |
| `Variable` | Variable is a named Boolean variable. |
| `Xnor` | Xnor is the negated exclusive disjunction (equivalence) of two operands. |
| `Xor` | Xor is exclusive disjunction of two operands. |

<details>
<summary><code>And</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (a And) Eval(env map[string]bool) bool` | Eval evaluates the conjunction. |
| `func (a And) String() string` | String renders the conjunction as "(x & y)". |

</details>

<details>
<summary><code>BDD</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewBDD(order []string) *BDD` | NewBDD returns a manager whose variables are tested in the given order, most significant (closest to the root) first. |
| `func (b *BDD) And(f, g Node) Node` | And returns the conjunction of f and g. |
| `func (b *BDD) AnySat(n Node) (map[string]bool, bool)` | AnySat returns one satisfying assignment of the function n as a name->value map and true, or nil and false when n is unsatisfiable. |
| `func (b *BDD) Apply(op string, f, g Node) Node` | Apply combines f and g with the named binary operator ("and", "or", "xor", "nand", "nor", "iff", "implies"). |
| `func (b *BDD) Equal(f, g Node) bool` | Equal reports whether two nodes denote the same Boolean function. |
| `func (b *BDD) Eval(n Node, env map[string]bool) bool` | Eval evaluates the function n under the given environment. |
| `func (b *BDD) FromExpr(e Expr) Node` | FromExpr builds the BDD node for the expression e, testing its variables in the manager's ordering. |
| `func (b *BDD) Iff(f, g Node) Node` | Iff returns the biconditional (equivalence) of f and g. |
| `func (b *BDD) Implies(f, g Node) Node` | Implies returns the conditional f -> g. |
| `func (b *BDD) IsContradiction(n Node) bool` | IsContradiction reports whether n is the false terminal. |
| `func (b *BDD) IsSatisfiable(n Node) bool` | IsSatisfiable reports whether n is not the false terminal. |
| `func (b *BDD) IsTautology(n Node) bool` | IsTautology reports whether n is the true terminal, i.e. |
| `func (b *BDD) Ite(f, g, h Node) Node` | Ite computes the if-then-else operator ite(f, g, h) = (f AND g) OR (NOT f AND h), the universal ternary connective from which all binary Boolean… |
| `func (b *BDD) Nand(f, g Node) Node` | Nand returns the negated conjunction of f and g. |
| `func (b *BDD) Nor(f, g Node) Node` | Nor returns the negated disjunction of f and g. |
| `func (b *BDD) Not(n Node) Node` | Not returns the negation of n. |
| `func (b *BDD) NumNodes() int` | NumNodes returns the total number of distinct nodes allocated, including the two terminals. |
| `func (b *BDD) Or(f, g Node) Node` | Or returns the disjunction of f and g. |
| `func (b *BDD) Order() []string` | Order returns the variable ordering used by the manager. |
| `func (b *BDD) Restrict(n Node, name string, val bool) Node` | Restrict returns the cofactor of n obtained by fixing the variable name to val (Shannon restriction). |
| `func (b *BDD) SatCount(n Node) int` | SatCount returns the number of assignments over all ordering variables under which the function represented by n is true. |
| `func (b *BDD) ToExpr(n Node) Expr` | ToExpr reconstructs a Boolean `Expr` from the node n by Shannon expansion, using the manager's variable names. |
| `func (b *BDD) Var(name string) Node` | Var returns the node representing the single variable name. |
| `func (b *BDD) Xor(f, g Node) Node` | Xor returns the exclusive-or of f and g. |

</details>

<details>
<summary><code>BoolConst</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (c BoolConst) Eval(map[string]bool) bool` | Eval returns the constant's value. |
| `func (c BoolConst) String() string` | String returns "T" or "F". |

</details>

<details>
<summary><code>CNF</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewCNF(clauses ...Clause) CNF` | NewCNF returns a CNF formula built from the given clauses. |
| `func ParseDIMACS(s string) (CNF, error)` | ParseDIMACS parses a formula from DIMACS CNF text. |
| `func PureLiteralElimination(f CNF, assign map[int]bool) (CNF, map[int]bool)` | PureLiteralElimination assigns every pure literal (one occurring with a single polarity) to the value that satisfies its clauses, repeating until no… |
| `func TseitinCNF(e Expr) CNF` | TseitinCNF returns just the equisatisfiable CNF from `Tseitin`. |
| `func UnitPropagate(f CNF, assign map[int]bool) (CNF, map[int]bool, bool)` | UnitPropagate repeatedly assigns the literal of every unit clause and simplifies the formula accordingly. |
| `func (f CNF) AddClause(c Clause) CNF` | AddClause returns a new formula with c appended to the conjunction. |
| `func (f CNF) Assign(l Lit) CNF` | Assign returns the formula conditioned on literal l being true: clauses satisfied by l are removed and the complement of l is deleted from the… |
| `func (f CNF) Clone() CNF` | Clone returns a deep copy of the formula. |
| `func (f CNF) Concat(g CNF) CNF` | Concat returns the conjunction of two CNF formulas. |
| `func (f CNF) CountLiterals() int` | CountLiterals returns the total number of literal occurrences in the formula. |
| `func (f CNF) DIMACS() string` | DIMACS renders the formula in the standard DIMACS CNF text format, with a "p cnf <vars> <clauses>" header followed by one zero-terminated clause per… |
| `func (f CNF) Eval(assign map[int]bool) bool` | Eval evaluates the formula under the given assignment: true when every clause is satisfied. |
| `func (f CNF) HasEmptyClause() bool` | HasEmptyClause reports whether the formula contains an empty clause, which makes it unsatisfiable. |
| `func (f CNF) IsEmpty() bool` | IsEmpty reports whether the formula has no clauses, in which case it is trivially satisfied. |
| `func (f CNF) Literals() []Lit` | Literals returns every literal occurrence across all clauses, with repetition. |
| `func (f CNF) MaxVar() int` | MaxVar returns the largest variable index occurring in the formula, or 0 when there are none. |
| `func (f CNF) Negate() DNF` | Negate on a CNF returns a DNF equivalent to its negation by De Morgan's law. |
| `func (f CNF) NumClauses() int` | NumClauses returns the number of clauses in the formula. |
| `func (f CNF) NumVars() int` | NumVars returns the count of distinct variables appearing in the formula. |
| `func (f CNF) PureLiterals() []Lit` | PureLiterals returns the literals that occur with only one polarity across the whole formula. |
| `func (f CNF) RemoveDuplicateClauses() CNF` | RemoveDuplicateClauses returns an equivalent formula with syntactically duplicate clauses collapsed to a single copy. |
| `func (f CNF) Simplify() CNF` | Simplify returns an equivalent formula with tautological clauses dropped and duplicate literals removed from each remaining clause. |
| `func (f CNF) String() string` | String renders the formula as a conjunction of parenthesised clauses joined by " & ". |
| `func (f CNF) UnitClauses() []Lit` | UnitClauses returns the literals of every unit (single-literal) clause in the formula. |
| `func (f CNF) Vars() []int` | Vars returns the sorted set of distinct variable indices in the formula. |

</details>

<details>
<summary><code>Clause</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewClause(lits ...Lit) Clause` | NewClause returns a clause built from the given literals. |
| `func Resolve(c1, c2 Clause, v int) (Clause, bool)` | Resolve computes the resolvent of two clauses on the pivot variable v: the disjunction of all their other literals, provided one clause contains the… |
| `func (c Clause) Clone() Clause` | Clone returns an independent copy of the clause. |
| `func (c Clause) Contains(l Lit) bool` | Contains reports whether the literal l appears in the clause. |
| `func (c Clause) ContainsVar(v int) bool` | ContainsVar reports whether variable v appears in the clause with either polarity. |
| `func (c Clause) Dedup() Clause` | Dedup returns a copy of the clause with duplicate literals removed, order preserved. |
| `func (c Clause) EvalAnd(assign map[int]bool) bool` | EvalAnd evaluates the clause as a conjunction under the assignment: true when every literal is satisfied. |
| `func (c Clause) EvalOr(assign map[int]bool) bool` | EvalOr evaluates the clause as a disjunction under the assignment: true when at least one literal is satisfied. |
| `func (c Clause) IsEmpty() bool` | IsEmpty reports whether the clause has no literals. |
| `func (c Clause) IsTautology() bool` | IsTautology reports whether the disjunctive clause contains a literal and its negation, making it trivially true. |
| `func (c Clause) IsUnit() bool` | IsUnit reports whether the clause consists of a single literal. |
| `func (c Clause) Len() int` | Len returns the number of literals in the clause. |
| `func (c Clause) MaxVar() int` | MaxVar returns the largest variable index occurring in the clause, or 0 for the empty clause. |
| `func (c Clause) Negate() Clause` | Negate returns the literal-wise complement of the clause. |
| `func (c Clause) Sorted() Clause` | Sorted returns a copy of the clause sorted by variable index then sign. |
| `func (c Clause) String() string` | String renders the clause as a parenthesised disjunction such as "(x1 \| ~x2 \| x3)". |
| `func (c Clause) Vars() []int` | Vars returns the sorted set of distinct variable indices occurring in the clause. |

</details>

<details>
<summary><code>DNF</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewDNF(terms ...Clause) DNF` | NewDNF returns a DNF formula built from the given terms. |
| `func (d DNF) AddTerm(t Clause) DNF` | AddTerm returns a new formula with t appended to the disjunction. |
| `func (d DNF) Clone() DNF` | Clone returns a deep copy of the formula. |
| `func (d DNF) Eval(assign map[int]bool) bool` | Eval evaluates the formula under the assignment: true when at least one term is fully satisfied. |
| `func (d DNF) IsEmpty() bool` | IsEmpty reports whether the formula has no terms, in which case it is unsatisfiable. |
| `func (d DNF) MaxVar() int` | MaxVar returns the largest variable index occurring in the formula, or 0. |
| `func (d DNF) Negate() CNF` | Negate returns a CNF formula logically equivalent to the negation of the DNF, applying De Morgan's law term-by-term. |
| `func (d DNF) NumTerms() int` | NumTerms returns the number of terms in the formula. |
| `func (d DNF) String() string` | String renders the formula as a disjunction of parenthesised conjunctive terms joined by " \| ". |
| `func (d DNF) Vars() []int` | Vars returns the sorted set of distinct variable indices in the formula. |

</details>

<details>
<summary><code>Expr</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func AndAll(es ...Expr) Expr` | AndAll returns the conjunction of all expressions. |
| `func AndE(x, y Expr) Expr` | AndE returns the conjunction of x and y. |
| `func CNFToExpr(f CNF, vm *VarMap) Expr` | CNFToExpr renders a clausal CNF as an `Expr` tree using the variable names in vm. |
| `func Constant(b bool) Expr` | Constant returns the constant expression for b. |
| `func DNFToExpr(d DNF, vm *VarMap) Expr` | DNFToExpr renders a clausal DNF as an `Expr` tree using the names in vm. |
| `func EliminateImplications(e Expr) Expr` | EliminateImplications rewrites an expression so that it uses only the And, Or and Not connectives, expanding conditionals, biconditionals and the… |
| `func IffE(x, y Expr) Expr` | IffE returns the biconditional x <-> y. |
| `func ImpliesE(x, y Expr) Expr` | ImpliesE returns the conditional x -> y. |
| `func MinimizeSOPExpr(minterms, dontcares []int, names []string) Expr` | MinimizeSOPExpr returns the minimised sum-of-products cover of the function as a Boolean `Expr` built from the given variable names. |
| `func MustParse(s string) Expr` | MustParse is like `Parse` but panics on error. |
| `func NandE(x, y Expr) Expr` | NandE returns the nand of x and y. |
| `func NorE(x, y Expr) Expr` | NorE returns the nor of x and y. |
| `func NotE(x Expr) Expr` | NotE returns the negation of x, folding double negations and constants. |
| `func OrAll(es ...Expr) Expr` | OrAll returns the disjunction of all expressions. |
| `func OrE(x, y Expr) Expr` | OrE returns the disjunction of x and y. |
| `func Parse(s string) (Expr, error)` | Parse parses a Boolean expression from its textual form and returns the corresponding `Expr` tree. |
| `func Simplify(e Expr) Expr` | Simplify applies Boolean simplification rules — constant folding, negation, idempotence, complementation, identity, domination and absorption —… |
| `func Substitute(e Expr, name string, repl Expr) Expr` | Substitute returns a copy of e with every occurrence of the variable named name replaced by the expression repl. |
| `func ToCNFExpr(e Expr) Expr` | ToCNFExpr converts an arbitrary Boolean expression into an equivalent expression in conjunctive normal form using De Morgan's laws and distribution… |
| `func ToDNFExpr(e Expr) Expr` | ToDNFExpr converts an arbitrary Boolean expression into an equivalent expression in disjunctive normal form. |
| `func ToNNF(e Expr) Expr` | ToNNF converts an expression to negation normal form, in which negations apply only to variables and the tree uses solely And, Or and Not (over… |
| `func V(name string) Expr` | V returns a variable expression with the given name. |
| `func XnorE(x, y Expr) Expr` | XnorE returns the xnor of x and y. |
| `func XorAll(es ...Expr) Expr` | XorAll returns the exclusive-or of all expressions. |
| `func XorE(x, y Expr) Expr` | XorE returns the exclusive-or of x and y. |

</details>

<details>
<summary><code>Iff</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (i Iff) Eval(env map[string]bool) bool` | Eval evaluates the biconditional. |
| `func (i Iff) String() string` | String renders the biconditional as "(x <-> y)". |

</details>

<details>
<summary><code>Implicant</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func EssentialPrimeImplicants(minterms, dontcares []int, numVars int) []Implicant` | EssentialPrimeImplicants returns the prime implicants that are the sole cover of at least one required minterm (a minterm in the on-set, ignoring… |
| `func MinimizeSOP(minterms, dontcares []int, numVars int) []Implicant` | MinimizeSOP returns a minimal (greedy) sum-of-products cover of the function as a list of implicants: all essential prime implicants plus a greedy… |
| `func PrimeImplicants(minterms, dontcares []int, numVars int) []Implicant` | PrimeImplicants returns the prime implicants of the Boolean function defined by the given on-set minterms and don't-care terms over numVars… |
| `func QuineMcCluskey(minterms, dontcares []int, numVars int) []Implicant` | QuineMcCluskey is an alias for `PrimeImplicants`, returning all prime implicants of the function. |
| `func (im Implicant) Covers(m int) bool` | Covers reports whether the cube covers the minterm m. |
| `func (im Implicant) Minterms() []int` | Minterms returns the sorted list of minterm indices covered by the cube. |
| `func (im Implicant) Size() int` | Size returns the number of minterms the cube covers, which is 2 raised to the number of dash positions. |
| `func (im Implicant) String() string` | String renders the cube as a pattern of '0', '1' and '-' characters, most-significant variable first. |
| `func (im Implicant) Term(names []string) string` | Term renders the cube as a product term using the given variable names, most significant first, e.g. |

</details>

<details>
<summary><code>Implies</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (i Implies) Eval(env map[string]bool) bool` | Eval evaluates the conditional. |
| `func (i Implies) String() string` | String renders the conditional as "(x -> y)". |

</details>

<details>
<summary><code>KarnaughMap</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewKarnaughMap(minterms, dontcares []int, numVars int) KarnaughMap` | NewKarnaughMap builds the Karnaugh map of the function defined by minterms and dontcares over numVars variables. |
| `func (k KarnaughMap) String() string` | String renders the Karnaugh map as a grid of '0', '1' and 'X' (don't-care) characters with Gray-coded axis labels. |

</details>

<details>
<summary><code>Lit</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func MakeLit(v int, negated bool) Lit` | MakeLit builds the literal for variable v (which must be >= 1). |
| `func NegLit(v int) Lit` | NegLit returns the negative literal for variable v. |
| `func PosLit(v int) Lit` | PosLit returns the positive literal for variable v. |
| `func (l Lit) IsNeg() bool` | IsNeg reports whether the literal is negated (its variable asserted false). |
| `func (l Lit) IsPos() bool` | IsPos reports whether the literal is positive (its variable asserted true). |
| `func (l Lit) Negate() Lit` | Negate returns the complementary literal, flipping its polarity. |
| `func (l Lit) Satisfied(assign map[int]bool) bool` | Satisfied reports whether the literal evaluates to true under the given variable assignment. |
| `func (l Lit) Sign() int` | Sign returns +1 for a positive literal and -1 for a negative literal. |
| `func (l Lit) String() string` | String renders the literal in the form "x3" for a positive literal and "~x3" for a negative literal. |
| `func (l Lit) Valid() bool` | Valid reports whether the literal refers to a real variable (index >= 1). |
| `func (l Lit) Var() int` | Var returns the underlying variable index of the literal, always a non-negative number, discarding the sign. |

</details>

<details>
<summary><code>Nand</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (n Nand) Eval(env map[string]bool) bool` | Eval evaluates the nand. |
| `func (n Nand) String() string` | String renders the nand as "(x !& y)". |

</details>

<details>
<summary><code>Nor</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (n Nor) Eval(env map[string]bool) bool` | Eval evaluates the nor. |
| `func (n Nor) String() string` | String renders the nor as "(x !\| y)". |

</details>

<details>
<summary><code>Not</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (n Not) Eval(env map[string]bool) bool` | Eval evaluates the negation. |
| `func (n Not) String() string` | String renders the negation as "~(x)". |

</details>

<details>
<summary><code>Or</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (o Or) Eval(env map[string]bool) bool` | Eval evaluates the disjunction. |
| `func (o Or) String() string` | String renders the disjunction as "(x \| y)". |

</details>

<details>
<summary><code>Solver</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewSolver(f CNF) *Solver` | NewSolver returns a solver for the given formula. |
| `func (s *Solver) AddClause(c Clause)` | AddClause adds a clause to the solver's formula and invalidates any cached result. |
| `func (s *Solver) Model() map[int]bool` | Model returns the satisfying assignment found by the most recent successful `Solver.Solve`. |
| `func (s *Solver) Solve() bool` | Solve runs the search and returns whether the formula is satisfiable, caching the outcome and model. |

</details>

<details>
<summary><code>TruthTable</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewTruthTable(e Expr) TruthTable` | NewTruthTable builds the complete truth table of e. |
| `func (t TruthTable) Column() []bool` | Column returns the output column of the table as a slice of booleans in row order. |
| `func (t TruthTable) Maxterms() []int` | Maxterms returns the row indices at which the expression is false. |
| `func (t TruthTable) Minterms() []int` | Minterms returns the row indices at which the expression is true. |
| `func (t TruthTable) NumRows() int` | NumRows returns the number of rows, which is 2 raised to the number of variables. |
| `func (t TruthTable) String() string` | String renders the truth table as fixed-width text with a header row. |

</details>

<details>
<summary><code>TseitinEncoder</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewTseitinEncoder() *TseitinEncoder` | NewTseitinEncoder returns an empty encoder ready to encode expressions. |
| `func (t *TseitinEncoder) Assert(e Expr)` | Assert records that the encoded expression must be true, by adding a unit clause on its top-level literal. |
| `func (t *TseitinEncoder) CNF() CNF` | CNF returns the accumulated clauses as a `CNF` formula. |
| `func (t *TseitinEncoder) Clauses() []Clause` | Clauses returns the clauses accumulated so far. |
| `func (t *TseitinEncoder) Encode(e Expr) Lit` | Encode registers the clauses that define e and returns a literal that is true exactly when e is true. |
| `func (t *TseitinEncoder) VarMap() *VarMap` | VarMap returns the variable map that records both original and auxiliary variable names introduced during encoding. |

</details>

<details>
<summary><code>VarMap</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewVarMap(names []string) *VarMap` | NewVarMap builds a VarMap assigning indices 1..n to the given names in order. |
| `func (m *VarMap) Add(name string) int` | Add returns the index of name, allocating a fresh one if it is new. |
| `func (m *VarMap) Index(name string) (int, bool)` | Index returns the index of name and whether it is present. |
| `func (m *VarMap) Len() int` | Len returns the number of variables in the map. |
| `func (m *VarMap) Name(i int) string` | Name returns the variable name for index i (1-based), or the empty string if out of range. |

</details>

<details>
<summary><code>Variable</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (v Variable) Eval(env map[string]bool) bool` | Eval looks the variable up in env, defaulting to false. |
| `func (v Variable) String() string` | String returns the variable name. |

</details>

<details>
<summary><code>Xnor</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (n Xnor) Eval(env map[string]bool) bool` | Eval evaluates the xnor. |
| `func (n Xnor) String() string` | String renders the xnor as "(x ~^ y)". |

</details>

<details>
<summary><code>Xor</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (x Xor) Eval(env map[string]bool) bool` | Eval evaluates the exclusive-or. |
| `func (x Xor) String() string` | String renders the exclusive-or as "(x ^ y)". |

</details>

### Variables

`ErrEmptyFormula`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/satsolver).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
