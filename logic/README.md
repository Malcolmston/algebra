# logic — implements propositional logic and Boolean algebra using only the Go standard library

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/logic.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/logic)

Package logic implements propositional logic and Boolean algebra using only
the Go standard library.

- Boolean primitives (ops.go): the elementary connectives `And`, `Or`, `Not`,
`Xor`, `Nand`, `Nor`, `Xnor`, `Implies` and `Iff`, together with the variadic
reductions `AndAll`, `OrAll`, `XorAll` and the derived gates `Majority` and
`Mux`. - Expression trees, parsing and evaluation (expr.go): the sealed `Expr`
interface implemented by `Var`, `Const`, `UnaryExpr` and `BinaryExpr`; the
constructors `NewVar`, `NewConst`, `NewNot`, `NewAnd`, `NewOr`, `NewXor`,
`NewNand`, `NewNor`, `NewXnor`, `NewImplies` and `NewIff`; the
`Parse`/`MustParse` recursive-descent parser; and the helpers `Vars` and
`EvalString`. - Truth tables (truthtable.go): the `TruthTable` and `TruthRow`
types built by `NewTruthTable`, plus `Assignments`, `Minterms`, `Maxterms`,
`IndexToAssignment` and `AssignmentToIndex`. - Semantic analysis (analyze.go):
`IsTautology`, `IsContradiction`, `IsSatisfiable`, `IsContingency`,
`FindModel`, `AllModels`, `CountModels`, `Equivalent` and `Entails`. - Normal
forms (normalform.go): `ToNNF`, `ToCNF`, `ToDNF`, `Simplify`, `CNFString`,
`DNFString`, `IsCNF` and `IsDNF`. - Two-level minimisation (minimize.go): the
`Implicant` cube type, the `QuineMcCluskey` engine with `PrimeImplicants` and
`EssentialPrimeImplicants`, the cover routines `MinimizeSOP` and `SOPString`,
the Gray-code helpers `GrayCode`, `GrayEncode` and `GrayDecode`, the
`KarnaughMap` grouping type built by `NewKarnaughMap`, and the bit utilities
`PopCount` and `BitString`.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/logic"
```

## Usage

This is the package's own `ExampleAndAll`, so it compiles and its output is
asserted on every `go test ./logic/`.

```go
fmt.Println(logic.AndAll(), logic.OrAll(), logic.XorAll())
	fmt.Println(logic.AndAll(true, true), logic.OrAll(false, true), logic.XorAll(true, true, true))
```

```
true false false
true true true
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func AllModels(e Expr) []map[string]bool` | AllModels returns every satisfying assignment of e, in ascending assignment-index order. |
| `func And(a, b bool) bool` | And returns the logical conjunction of a and b: true when both operands are true. |
| `func AndAll(vals ...bool) bool` | AndAll returns the conjunction of every value in vals. |
| `func AssignmentToIndex(env map[string]bool, vars []string) int` | AssignmentToIndex encodes an assignment over vars back into its row index, the inverse of `IndexToAssignment`. |
| `func Assignments(vars []string) []map[string]bool` | Assignments enumerates all 2^len(vars) truth assignments over vars, in ascending index order. |
| `func BitString(x, numVars int) string` | BitString returns the binary representation of x in exactly numVars digits, most-significant bit first. |
| `func CNFString(e Expr) string` | CNFString renders the canonical conjunctive normal form of e as a product-of- sums string. |
| `func CountModels(e Expr) int` | CountModels returns the number of assignments that satisfy e. |
| `func DNFString(e Expr) string` | DNFString renders the canonical disjunctive normal form of e as a sum-of- products string, using juxtaposition-free "&"/"\|" notation with negation… |
| `func Entails(a, b Expr) bool` | Entails reports whether a logically entails b (written a ⊨ b): every assignment that satisfies a also satisfies b. |
| `func Equivalent(a, b Expr) bool` | Equivalent reports whether a and b are logically equivalent: they yield the same truth value under every assignment of their combined variables. |
| `func EvalString(s string, env map[string]bool) (bool, error)` | EvalString parses s and evaluates it under env in a single call. |
| `func FindModel(e Expr) (map[string]bool, bool)` | FindModel returns a satisfying assignment for e and true, or a nil map and false when e is unsatisfiable. |
| `func GrayCode(n int) []int` | GrayCode returns the reflected binary Gray code sequence for n bits: the 2^n values GrayEncode(0), GrayEncode(1), ..., GrayEncode(2^n-1), in which… |
| `func GrayDecode(g int) int` | GrayDecode inverts `GrayEncode`, recovering the integer whose Gray code is g. |
| `func GrayEncode(n int) int` | GrayEncode returns the reflected binary Gray code of the non-negative integer n, computed as n XOR (n>>1). |
| `func Iff(a, b bool) bool` | Iff returns the biconditional a <-> b (logical equivalence): true when the operands share the same truth value. |
| `func Implies(a, b bool) bool` | Implies returns the material conditional a -> b, which is false in exactly one case: a true and b false. |
| `func IndexToAssignment(i int, vars []string) map[string]bool` | IndexToAssignment decodes the row index i into a variable assignment over vars. |
| `func IsCNF(e Expr) bool` | IsCNF reports whether e is in conjunctive normal form: a conjunction of one or more disjunctions of literals. |
| `func IsContingency(e Expr) bool` | IsContingency reports whether e is neither a tautology nor a contradiction, that is, whether its value depends on the assignment. |
| `func IsContradiction(e Expr) bool` | IsContradiction reports whether e is false under every assignment of its variables. |
| `func IsDNF(e Expr) bool` | IsDNF reports whether e is in disjunctive normal form: a disjunction of one or more conjunctions of literals, where a literal is a variable, a… |
| `func IsSatisfiable(e Expr) bool` | IsSatisfiable reports whether at least one assignment makes e true. |
| `func IsTautology(e Expr) bool` | IsTautology reports whether e is true under every assignment of its variables. |
| `func Majority(a, b, c bool) bool` | Majority returns the value taken by the majority of its three inputs: true when at least two of a, b and c are true. |
| `func Maxterms(e Expr) []int` | Maxterms returns the sorted indices of the assignments for which e is false, using the same indexing convention as `Minterms`. |
| `func MinimizeString(e Expr) string` | MinimizeString minimises e and returns its minimal sum-of-products form as a string over the variables of e. |
| `func Minterms(e Expr) []int` | Minterms returns the sorted indices of the assignments for which e is true. |
| `func Mux(sel, a, b bool) bool` | Mux returns a two-to-one multiplexer output: b when sel is true, otherwise a. |
| `func Nand(a, b bool) bool` | Nand returns the negated conjunction of a and b, equivalent to Not(And(a, b)). |
| `func Nor(a, b bool) bool` | Nor returns the negated disjunction of a and b, equivalent to Not(Or(a, b)). |
| `func Not(a bool) bool` | Not returns the logical negation of a. |
| `func Or(a, b bool) bool` | Or returns the logical (inclusive) disjunction of a and b: true when at least one operand is true. |
| `func OrAll(vals ...bool) bool` | OrAll returns the disjunction of every value in vals. |
| `func PopCount(x uint) int` | PopCount returns the number of one bits (the Hamming weight) of x. |
| `func SOPString(cover []Implicant, vars []string) string` | SOPString renders a cover produced by `QuineMcCluskey` or `MinimizeSOP` as a sum-of-products string over vars. |
| `func Vars(e Expr) []string` | Vars returns the sorted, de-duplicated set of variable names appearing in e. |
| `func Xnor(a, b bool) bool` | Xnor returns the negated exclusive disjunction (logical equivalence) of a and b: true when the operands are equal. |
| `func Xor(a, b bool) bool` | Xor returns the exclusive disjunction of a and b: true when the operands differ. |
| `func XorAll(vals ...bool) bool` | XorAll returns the parity of vals: true when an odd number of the values are true. |

### Types

| Type | What it is |
| --- | --- |
| `BinaryExpr` | BinaryExpr applies a `BinaryOp` to a left and right operand. |
| `BinaryOp` | BinaryOp identifies a binary connective. |
| `Const` | Const is a Boolean constant, either true or false. |
| `Expr` | Expr is a node in a propositional-logic expression tree. |
| `Implicant` | Implicant is a product-term cube over numVars variables. |
| `KarnaughMap` | KarnaughMap is a Karnaugh map: a Gray-coded grid view of a Boolean function used to visualise adjacency and read off groupings. |
| `TruthRow` | TruthRow is a single line of a `TruthTable`: an assignment of truth values to the variables together with the resulting value of the expression. |
| `TruthTable` | TruthTable is the exhaustive evaluation of an expression over all 2^n assignments of its n variables. |
| `UnaryExpr` | UnaryExpr applies a `UnaryOp` to a single operand. |
| `UnaryOp` | UnaryOp identifies a unary connective. |
| `Var` | Var is a propositional variable identified by its name. |

<details>
<summary><code>BinaryExpr</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewAnd(l, r Expr) *BinaryExpr` | NewAnd returns the conjunction l & r. |
| `func NewIff(l, r Expr) *BinaryExpr` | NewIff returns the biconditional l <-> r. |
| `func NewImplies(l, r Expr) *BinaryExpr` | NewImplies returns the material conditional l -> r. |
| `func NewNand(l, r Expr) *BinaryExpr` | NewNand returns the negated conjunction of l and r. |
| `func NewNor(l, r Expr) *BinaryExpr` | NewNor returns the negated disjunction of l and r. |
| `func NewOr(l, r Expr) *BinaryExpr` | NewOr returns the disjunction l \| r. |
| `func NewXnor(l, r Expr) *BinaryExpr` | NewXnor returns the negated exclusive disjunction (equivalence) of l and r. |
| `func NewXor(l, r Expr) *BinaryExpr` | NewXor returns the exclusive disjunction l ^ r. |
| `func (e *BinaryExpr) Eval(env map[string]bool) (bool, error)` | Eval evaluates both operands and applies the binary connective. |
| `func (e *BinaryExpr) String() string` | String renders the binary expression with explicit parentheses. |

</details>

<details>
<summary><code>BinaryOp</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (op BinaryOp) String() string` | String returns the canonical source symbol for the binary operator. |

</details>

<details>
<summary><code>Const</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewConst(b bool) Const` | NewConst returns the constant expression for b. |
| `func (c Const) Eval(map[string]bool) (bool, error)` | Eval returns the constant's value; it never fails. |
| `func (c Const) String() string` | String renders the constant as "T" or "F". |

</details>

<details>
<summary><code>Expr</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func MustParse(s string) Expr` | MustParse is like `Parse` but panics if s does not parse. |
| `func Parse(s string) (Expr, error)` | Parse parses a propositional-logic expression. |
| `func Simplify(e Expr) Expr` | Simplify returns an equivalent expression reduced by constant folding and the standard Boolean identities (identity, domination, idempotence and… |
| `func ToCNF(e Expr) Expr` | ToCNF returns the canonical (full) conjunctive normal form of e: an And of maxterms, one disjunction of literals for each falsifying assignment. |
| `func ToDNF(e Expr) Expr` | ToDNF returns the canonical (full) disjunctive normal form of e: an Or of minterms, one conjunction of literals for each satisfying assignment. |
| `func ToNNF(e Expr) Expr` | ToNNF converts e to negation normal form: an equivalent expression built only from And, Or and negations applied directly to variables or constants. |

</details>

<details>
<summary><code>Implicant</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func EssentialPrimeImplicants(minterms, dontCares []int, numVars int) []Implicant` | EssentialPrimeImplicants returns the essential prime implicants of the function defined by minterms and dontCares over numVars variables: those prime… |
| `func MinimizeSOP(e Expr) []Implicant` | MinimizeSOP minimises the expression e and returns a minimal sum-of-products cover as a set of implicants over the variables of `Vars`(e). |
| `func PrimeImplicants(minterms, dontCares []int, numVars int) []Implicant` | PrimeImplicants computes the prime implicants of the Boolean function defined by the given minterms and don't-care terms over numVars variables,… |
| `func QuineMcCluskey(minterms, dontCares []int, numVars int) []Implicant` | QuineMcCluskey minimises the Boolean function given by minterms and dontCares over numVars variables and returns a sum-of-products cover as a set of… |
| `func (im Implicant) Covers(minterm int) bool` | Covers reports whether the cube contains the given minterm index. |
| `func (im Implicant) Expression(vars []string) string` | Expression renders the cube as a product term over the given variable names, where vars[0] is the most-significant position. |
| `func (im Implicant) LiteralCount() int` | LiteralCount returns the number of fixed variable positions in the cube, i.e. |
| `func (im Implicant) Minterms() []int` | Minterms returns the original minterms that the cube covers, as recorded when it was produced by the minimisation routines. |
| `func (im Implicant) String() string` | String renders the cube as a pattern of numVars characters, most-significant position first, using '1', '0' and '-' for a dash. |

</details>

<details>
<summary><code>KarnaughMap</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewKarnaughMap(minterms []int, numVars int) *KarnaughMap` | NewKarnaughMap builds a Karnaugh map for numVars variables in which the given minterms are set to true. |
| `func (k *KarnaughMap) At(row, col int) bool` | At returns the cell value at the given grid row and column. |
| `func (k *KarnaughMap) Cols() int` | Cols returns the number of columns in the map (2 raised to the column-variable count). |
| `func (k *KarnaughMap) Groups() []Implicant` | Groups returns the prime implicants of the mapped function, which correspond exactly to the maximal rectangular groupings of adjacent true cells. |
| `func (k *KarnaughMap) Rows() int` | Rows returns the number of rows in the map (2 raised to the row-variable count). |
| `func (k *KarnaughMap) String() string` | String renders the map as a Gray-coded grid of 1s and 0s with binary row and column headers. |

</details>

<details>
<summary><code>TruthTable</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewTruthTable(e Expr) *TruthTable` | NewTruthTable evaluates e over every assignment of its variables and returns the resulting truth table. |
| `func (t *TruthTable) String() string` | String renders the truth table as an aligned ASCII grid with a header row of variable names, a "=" result column, and one line per assignment using 1… |

</details>

<details>
<summary><code>UnaryExpr</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewNot(x Expr) *UnaryExpr` | NewNot returns the negation of x. |
| `func (e *UnaryExpr) Eval(env map[string]bool) (bool, error)` | Eval evaluates the operand and applies the unary connective. |
| `func (e *UnaryExpr) String() string` | String renders the negation, parenthesising a compound operand. |

</details>

<details>
<summary><code>UnaryOp</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (op UnaryOp) String() string` | String returns the source symbol for the unary operator. |

</details>

<details>
<summary><code>Var</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewVar(name string) Var` | NewVar returns the variable with the given name. |
| `func (v Var) Eval(env map[string]bool) (bool, error)` | Eval looks the variable up in env, returning an error when it is unbound. |
| `func (v Var) String() string` | String returns the variable's name. |

</details>

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/logic).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
