# modelchecking — implements temporal-logic model checking over finite Kripke structures in pure Go using

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/modelchecking.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/modelchecking)

Package modelchecking implements temporal-logic model checking over finite
Kripke structures in pure Go using only the standard library.

- `Kripke` is a finite labelled transition system over atomic propositions.
`StateSet` is a packed-bitset set of states supporting the Boolean operations
that drive the symbolic fixpoint computations.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/modelchecking"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func BMCReachInvariant(k *Kripke, p StateSet, bound int) bool` | BMCReachInvariant is the bounded analogue of proving AG p: it reports whether no violation of p is reachable within bound steps. |
| `func BoundedReachable(k *Kripke, target StateSet, bound int) ([]int, bool)` | BoundedReachable searches for a path of length at most bound from an initial state to some state in target. |
| `func CTLModelCheck(k *Kripke, f *CTL) (bool, error)` | CTLModelCheck reports whether every initial state of k satisfies f. |
| `func CTLSatisfies(k *Kripke, s int, f *CTL) (bool, error)` | CTLSatisfies reports whether state s of k satisfies f. |
| `func EFWitness(k *Kripke, s int, target StateSet) ([]int, bool)` | EFWitness returns a finite path from state s to a state satisfying target, witnessing EF target at s. |
| `func EUWitness(k *Kripke, s int, a, b StateSet) ([]int, bool)` | EUWitness returns a finite path from s along which a holds until b, witnessing E[a U b] at s. |
| `func EXWitness(k *Kripke, s int, target StateSet) (int, bool)` | EXWitness returns a single successor of s that satisfies target, witnessing EX target. |
| `func ExistsPathTo(k *Kripke, sources []int, target StateSet) []int` | ExistsPathTo returns a shortest path (inclusive of endpoints) from any of the sources to some state of target, following transitions of k without a… |
| `func FairCTLModelCheck(k *Kripke, f *CTL, fc Fairness) (bool, error)` | FairCTLModelCheck reports whether every initial state of k satisfies f under the fairness condition fc. |
| `func LTLEquivalent(f, g *LTL) bool` | LTLEquivalent reports whether two LTL formulas are satisfied by exactly the same infinite words, decided by checking validity of their biconditional. |
| `func LTLValid(f *LTL) bool` | LTLValid reports whether the LTL formula f holds on every infinite word, which is the case exactly when its negation is unsatisfiable. |

### Types

| Type | What it is |
| --- | --- |
| `Buchi` | Buchi is a nondeterministic Büchi automaton over the alphabet 2^AP with a single acceptance set. |
| `BuchiEdge` | BuchiEdge is a labelled transition of a `Buchi` or `GenBuchi` automaton. |
| `CTL` | CTL is a node of a computation tree logic formula tree. |
| `CTLKind` | CTLKind enumerates the node kinds of a computation tree logic syntax tree. |
| `Counterexample` | Counterexample records a violating infinite path of a Kripke structure, found while model checking an LTL (or path-form CTL) property. |
| `Fairness` | Fairness is a set of (unconditional / justice) fairness constraints for a Kripke structure. |
| `GenBuchi` | GenBuchi is a generalized Büchi automaton: like `Buchi` but with a family of acceptance sets. |
| `Guard` | Guard is a conjunction of literals labelling a Büchi automaton edge. |
| `Kripke` | Kripke is a finite Kripke structure: a labelled transition system over a set of atomic propositions. |
| `LTL` | LTL is a node of a linear temporal logic formula tree. |
| `LTLKind` | LTLKind enumerates the node kinds of a linear temporal logic syntax tree. |
| `Lasso` | Lasso is an ultimately periodic accepting run of a Büchi automaton or the witness path of a model-checking counterexample. |
| `StateSet` | StateSet is a fixed-universe set of states drawn from {0, 1, ..., n-1}. |
| `Unrolling` | Unrolling records the breadth-first unrolling of a Kripke structure's transition relation from its initial states up to a fixed bound. |

<details>
<summary><code>Buchi</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Degeneralize(g *GenBuchi) *Buchi` | Degeneralize converts a generalized Büchi automaton with k acceptance sets into an equivalent ordinary `Buchi` automaton using the standard counter… |
| `func LTLToBuchi(f *LTL) *Buchi` | LTLToBuchi builds an ordinary Büchi automaton equivalent to the LTL formula f by degeneralizing the output of `LTLToGenBuchi`. |
| `func NewBuchi(n int) *Buchi` | NewBuchi returns a Büchi automaton with n states, no edges, no initial states and no accepting states. |
| `func ProductBuchiKripke(k *Kripke, b *Buchi) *Buchi` | ProductBuchiKripke builds the synchronous product of a Kripke structure k and a Büchi automaton b whose edges are source-labelled by literal guards. |
| `func (b *Buchi) AcceptingLasso() (Lasso, bool)` | AcceptingLasso searches for an accepting run and returns it as a `Lasso`. |
| `func (b *Buchi) AcceptingStates() []int` | AcceptingStates returns the accepting states in ascending order. |
| `func (b *Buchi) AddEdge(from, to int, g Guard)` | AddEdge adds a guarded transition from -> to. |
| `func (b *Buchi) Edges(s int) []BuchiEdge` | Edges returns the outgoing edges of state s. |
| `func (b *Buchi) InitialStates() []int` | InitialStates returns the initial states in ascending order. |
| `func (b *Buchi) IsAccepting(s int) bool` | IsAccepting reports whether s is an accepting state. |
| `func (b *Buchi) IsEmpty() bool` | IsEmpty reports whether the automaton accepts no infinite word, i.e. |
| `func (b *Buchi) IsEmptySCC() bool` | IsEmptySCC decides emptiness via strongly connected components: the language is non-empty iff some non-trivial SCC reachable from an initial state… |
| `func (b *Buchi) IsInitial(s int) bool` | IsInitial reports whether s is an initial state. |
| `func (b *Buchi) NumEdges() int` | NumEdges returns the total number of edges. |
| `func (b *Buchi) NumStates() int` | NumStates returns the number of states. |
| `func (b *Buchi) Propositions() []string` | Propositions returns the recorded atomic proposition alphabet. |
| `func (b *Buchi) ReachableStates() StateSet` | ReachableStates returns the set of states reachable from the initial states via non-contradictory edges. |
| `func (b *Buchi) SCCs() [][]int` | SCCs returns the strongly connected components of the automaton. |
| `func (b *Buchi) SetAccepting(states ...int)` | SetAccepting marks the given states as accepting. |
| `func (b *Buchi) SetInitial(states ...int)` | SetInitial marks the given states as initial. |
| `func (b *Buchi) SetPropositions(props []string)` | SetPropositions records the atomic proposition alphabet, used for reporting. |
| `func (b *Buchi) String() string` | String renders a multi-line description of the automaton. |
| `func (b *Buchi) Successors(s int) []int` | Successors returns the distinct successor states of s in ascending order. |

</details>

<details>
<summary><code>CTL</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func AF(f *CTL) *CTL` | AF returns the formula AF f: every path eventually reaches f. |
| `func AG(f *CTL) *CTL` | AG returns the formula AG f: every path satisfies f globally. |
| `func AR(f, g *CTL) *CTL` | AR returns the formula A[f R g]: every path satisfies g released by f. |
| `func AU(f, g *CTL) *CTL` | AU returns the formula A[f U g]: every path satisfies f until g. |
| `func AX(f *CTL) *CTL` | AX returns the formula AX f: every successor satisfies f. |
| `func CTLAnd(f, g *CTL) *CTL` | CTLAnd returns the conjunction f ∧ g. |
| `func CTLBot() *CTL` | CTLBot returns the constant false. |
| `func CTLIff(f, g *CTL) *CTL` | CTLIff returns the biconditional f ↔ g. |
| `func CTLImplies(f, g *CTL) *CTL` | CTLImplies returns the implication f → g. |
| `func CTLNot(f *CTL) *CTL` | CTLNot returns the negation ¬f. |
| `func CTLOr(f, g *CTL) *CTL` | CTLOr returns the disjunction f ∨ g. |
| `func CTLTop() *CTL` | CTLTop returns the constant true. |
| `func CTLVar(name string) *CTL` | CTLVar returns an atomic proposition node for name. |
| `func EF(f *CTL) *CTL` | EF returns the formula EF f: some path eventually reaches f. |
| `func EG(f *CTL) *CTL` | EG returns the formula EG f: some path satisfies f globally. |
| `func ER(f, g *CTL) *CTL` | ER returns the formula E[f R g]: some path satisfies g released by f. |
| `func EU(f, g *CTL) *CTL` | EU returns the formula E[f U g]: some path satisfies f until g. |
| `func EX(f *CTL) *CTL` | EX returns the formula EX f: some successor satisfies f. |
| `func MustParseCTL(s string) *CTL` | MustParseCTL parses s and panics on error. |
| `func ParseCTL(s string) (*CTL, error)` | ParseCTL parses a CTL formula from text. |
| `func (f *CTL) Atoms() []string` | Atoms returns the atomic proposition names occurring in f, sorted. |
| `func (f *CTL) Clone() *CTL` | Clone returns a deep copy of f. |
| `func (f *CTL) Equal(g *CTL) bool` | Equal reports structural equality of two formulas. |
| `func (f *CTL) ExistentialNormalForm() *CTL` | ExistentialNormalForm rewrites f into the existential fragment {¬, ∧, EX, EU, EG} plus atoms, which is the minimal set of operators the fixpoint… |
| `func (f *CTL) Height() int` | Height returns the height of the formula tree. |
| `func (f *CTL) IsTemporal() bool` | IsTemporal reports whether f's root is a temporal (path-quantified) operator. |
| `func (f *CTL) Size() int` | Size returns the number of nodes in the formula tree. |
| `func (f *CTL) String() string` | String renders f with ASCII operators, e.g. |
| `func (f *CTL) Subformulas() []*CTL` | Subformulas returns every distinct subformula of f (including f) ordered by increasing size, the order in which the labelling model-checking… |

</details>

<details>
<summary><code>Counterexample</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func AFCounterexample(k *Kripke, p StateSet) (*Counterexample, bool)` | AFCounterexample returns an infinite path from an initial state that never reaches p, witnessing the failure of AF p (equivalently a witness to EG… |
| `func AGCounterexample(k *Kripke, p StateSet) (*Counterexample, bool)` | AGCounterexample returns a finite path from an initial state to a state violating p, witnessing the failure of AG p. |
| `func BMCInvariant(k *Kripke, inv StateSet, bound int) (*Counterexample, bool)` | BMCInvariant performs bounded model checking of the invariant "always inv": it searches for a path of length at most bound from an initial state to a… |
| `func LTLModelCheck(k *Kripke, f *LTL) (bool, *Counterexample, error)` | LTLModelCheck checks whether every initial path of k satisfies the LTL formula f. |
| `func (c *Counterexample) StatePath() []int` | StatePath returns one unrolling of the counterexample path (prefix followed by one copy of the loop). |
| `func (c *Counterexample) String() string` | String renders the counterexample as a sequence of state names, marking the beginning of the periodic part with "(loop)". |

</details>

<details>
<summary><code>Fairness</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewFairness(sets ...StateSet) Fairness` | NewFairness returns a fairness condition consisting of the given constraint sets, each copied. |
| `func (f Fairness) Add(set StateSet) Fairness` | Add appends a constraint set and returns the updated condition. |
| `func (f Fairness) NumConstraints() int` | NumConstraints returns the number of fairness constraints. |

</details>

<details>
<summary><code>GenBuchi</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func LTLToGenBuchi(f *LTL) *GenBuchi` | LTLToGenBuchi builds a generalized Büchi automaton whose language is exactly the set of infinite words satisfying the LTL formula f, using the… |
| `func NewGenBuchi(n int) *GenBuchi` | NewGenBuchi returns a generalized Büchi automaton with n states and no acceptance sets. |
| `func (g *GenBuchi) AcceptanceSet(i int) StateSet` | AcceptanceSet returns the i-th acceptance set. |
| `func (g *GenBuchi) AddAcceptanceSet(set StateSet)` | AddAcceptanceSet appends an acceptance condition (a set of states that must be visited infinitely often). |
| `func (g *GenBuchi) AddEdge(from, to int, guard Guard)` | AddEdge adds a guarded transition from -> to. |
| `func (g *GenBuchi) Edges(s int) []BuchiEdge` | Edges returns the outgoing edges of state s. |
| `func (g *GenBuchi) InitialStates() []int` | InitialStates returns the initial states in ascending order. |
| `func (g *GenBuchi) IsInitial(s int) bool` | IsInitial reports whether s is an initial state. |
| `func (g *GenBuchi) NumAcceptanceSets() int` | NumAcceptanceSets returns the number of acceptance sets. |
| `func (g *GenBuchi) NumEdges() int` | NumEdges returns the total number of edges. |
| `func (g *GenBuchi) NumStates() int` | NumStates returns the number of states. |
| `func (g *GenBuchi) Propositions() []string` | Propositions returns the recorded atomic proposition alphabet. |
| `func (g *GenBuchi) SetInitial(states ...int)` | SetInitial marks the given states as initial. |
| `func (g *GenBuchi) SetPropositions(props []string)` | SetPropositions records the atomic proposition alphabet. |

</details>

<details>
<summary><code>Guard</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewGuard(pos, neg []string) Guard` | NewGuard returns a guard requiring every proposition in pos to hold and every proposition in neg to fail. |
| `func TrueGuard() Guard` | TrueGuard returns the guard that matches every letter. |
| `func (g Guard) Contradictory() bool` | Contradictory reports whether the guard requires a proposition to be both true and false, so that it matches no letter. |
| `func (g Guard) IsTrue() bool` | IsTrue reports whether the guard imposes no constraints. |
| `func (g Guard) Satisfies(label map[string]bool) bool` | Satisfies reports whether a state label (as a set of true propositions) satisfies the guard: all positive literals present, all negative literals… |
| `func (g Guard) String() string` | String renders the guard as a conjunction such as "a & !b", or "true" when unconstrained. |

</details>

<details>
<summary><code>Kripke</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewKripke(n int) *Kripke` | NewKripke returns a Kripke structure with n states, no transitions, empty labels and no initial states. |
| `func (k *Kripke) AddLabel(s int, p string) error` | AddLabel adds proposition p to the label of state s. |
| `func (k *Kripke) AddSelfLoop(s int) error` | AddSelfLoop adds the edge s -> s. |
| `func (k *Kripke) AddTransition(from, to int) error` | AddTransition adds a directed edge from state from to state to. |
| `func (k *Kripke) ClearInitial(s int)` | ClearInitial unmarks state s as initial. |
| `func (k *Kripke) Clone() *Kripke` | Clone returns a deep copy of the structure. |
| `func (k *Kripke) Deadlocks() []int` | Deadlocks returns the states that have no successors. |
| `func (k *Kripke) HasTransition(from, to int) bool` | HasTransition reports whether the edge from -> to exists. |
| `func (k *Kripke) Holds(s int, p string) bool` | Holds reports whether atomic proposition p holds in state s. |
| `func (k *Kripke) InitialStates() []int` | InitialStates returns the initial states in ascending order. |
| `func (k *Kripke) IsInitial(s int) bool` | IsInitial reports whether state s is an initial state. |
| `func (k *Kripke) IsTotal() bool` | IsTotal reports whether every state has at least one successor. |
| `func (k *Kripke) Label(s int) []string` | Label returns the atomic propositions holding in state s, sorted. |
| `func (k *Kripke) LabelSet(s int) map[string]bool` | LabelSet returns the label of state s as a boolean map. |
| `func (k *Kripke) LabelStateSet(p string) StateSet` | LabelStateSet returns the set of states whose label contains proposition p. |
| `func (k *Kripke) MakeTotal() int` | MakeTotal adds a self-loop to every deadlock state so that the transition relation becomes total, as required by the usual CTL/LTL semantics over… |
| `func (k *Kripke) NumStates() int` | NumStates returns the number of states in the structure. |
| `func (k *Kripke) NumTransitions() int` | NumTransitions returns the total number of edges in the structure. |
| `func (k *Kripke) OutDegree(s int) int` | OutDegree returns the number of successors of state s. |
| `func (k *Kripke) Predecessors(s int) []int` | Predecessors returns the states with an edge into s, in ascending order. |
| `func (k *Kripke) Propositions() []string` | Propositions returns every atomic proposition mentioned in any state label, sorted and de-duplicated. |
| `func (k *Kripke) Reachable() StateSet` | Reachable returns the set of states reachable from the initial states by following transitions. |
| `func (k *Kripke) ReachableFrom(sources ...int) StateSet` | ReachableFrom returns the set of states reachable from the given sources. |
| `func (k *Kripke) RemoveLabel(s int, p string) error` | RemoveLabel removes proposition p from the label of state s. |
| `func (k *Kripke) RemoveTransition(from, to int) error` | RemoveTransition removes the edge from -> to if present. |
| `func (k *Kripke) Reverse() *Kripke` | Reverse returns a new Kripke structure with the same states, labels and initial marks but with every transition reversed. |
| `func (k *Kripke) SetInitial(states ...int)` | SetInitial marks the given states as initial (in addition to any already marked). |
| `func (k *Kripke) SetLabel(s int, props ...string) error` | SetLabel replaces the label of state s with exactly the given propositions. |
| `func (k *Kripke) SetStateName(s int, name string) error` | SetStateName attaches a human-readable name to state s, used by `Kripke.String` and trace rendering. |
| `func (k *Kripke) StateName(s int) string` | StateName returns the name of state s, or a default "sN" if none was set. |
| `func (k *Kripke) String() string` | String renders a compact multi-line description of the structure. |
| `func (k *Kripke) Successors(s int) []int` | Successors returns the successors of state s in ascending order. |
| `func (k *Kripke) Validate() error` | Validate checks structural invariants: every successor index is in range. |

</details>

<details>
<summary><code>LTL</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func LTLAnd(f, g *LTL) *LTL` | LTLAnd returns the conjunction f ∧ g. |
| `func LTLBot() *LTL` | LTLBot returns the constant false. |
| `func LTLEventually(f *LTL) *LTL` | LTLEventually returns the future formula F f (eventually f). |
| `func LTLGlobally(f *LTL) *LTL` | LTLGlobally returns the always formula G f (globally f). |
| `func LTLIff(f, g *LTL) *LTL` | LTLIff returns the biconditional f ↔ g. |
| `func LTLImplies(f, g *LTL) *LTL` | LTLImplies returns the implication f → g. |
| `func LTLNext(f *LTL) *LTL` | LTLNext returns the next-time formula X f. |
| `func LTLNot(f *LTL) *LTL` | LTLNot returns the negation of f. |
| `func LTLOr(f, g *LTL) *LTL` | LTLOr returns the disjunction f ∨ g. |
| `func LTLRelease(f, g *LTL) *LTL` | LTLRelease returns the release f R g (the dual of Until). |
| `func LTLTop() *LTL` | LTLTop returns the constant true. |
| `func LTLUntil(f, g *LTL) *LTL` | LTLUntil returns the strong until f U g. |
| `func LTLVar(name string) *LTL` | LTLVar returns an atomic proposition node for name. |
| `func LTLWeakUntil(f, g *LTL) *LTL` | LTLWeakUntil returns the weak until f W g, equivalent to (f U g) ∨ G f. |
| `func MustParseLTL(s string) *LTL` | MustParseLTL parses s and panics on error. |
| `func ParseLTL(s string) (*LTL, error)` | ParseLTL parses an LTL formula from its textual representation and returns the syntax tree. |
| `func (f *LTL) Atoms() []string` | Atoms returns the atomic proposition names occurring in f, sorted and de-duplicated. |
| `func (f *LTL) Clone() *LTL` | Clone returns a deep copy of f. |
| `func (f *LTL) Equal(g *LTL) bool` | Equal reports structural equality of two formulas. |
| `func (f *LTL) Height() int` | Height returns the height of the formula tree (a single node has height 1). |
| `func (f *LTL) IsLiteral() bool` | IsLiteral reports whether f is true, false, an atom or the negation of an atom. |
| `func (f *LTL) IsTemporal() bool` | IsTemporal reports whether f's root is a temporal operator (X, F, G, U, R, W). |
| `func (f *LTL) NNF() *LTL` | NNF returns an equivalent formula in negation normal form: negations apply only to atoms, implications and biconditionals are eliminated, and the… |
| `func (f *LTL) Negate() *LTL` | Negate returns a formula equivalent to ¬f without pushing the negation inward. |
| `func (f *LTL) Simplify() *LTL` | Simplify applies a set of validity-preserving Boolean and temporal simplifications (constant folding, idempotence, absorption of true/false) and… |
| `func (f *LTL) Size() int` | Size returns the number of nodes in the formula tree. |
| `func (f *LTL) String() string` | String renders f using ASCII operators: ! & \| -> <-> X F G U R W. |
| `func (f *LTL) Subformulas() []*LTL` | Subformulas returns every distinct subformula of f (including f itself) in a deterministic order by increasing size. |

</details>

<details>
<summary><code>Lasso</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func BMCExistsGlobally(k *Kripke, p StateSet, bound int) (Lasso, bool)` | BMCExistsGlobally searches within the bound for a lasso all of whose states satisfy p, i.e. |
| `func BMCFindLasso(k *Kripke, bound int) (Lasso, bool)` | BMCFindLasso searches for an ultimately periodic path (a lasso) reachable from an initial state whose total length (prefix plus loop) is at most… |
| `func EGWitness(k *Kripke, s int, invariant StateSet) (Lasso, bool)` | EGWitness returns an infinite (ultimately periodic) path starting at s all of whose states satisfy invariant, witnessing EG invariant at s. |
| `func LTLSatisfiable(f *LTL) (bool, Lasso)` | LTLSatisfiable reports whether the LTL formula f is satisfiable, i.e. |
| `func (l Lasso) Length() int` | Length returns the number of distinct states listed in the lasso (prefix plus loop). |
| `func (l Lasso) LoopStart() int` | LoopStart returns the first state of the periodic part, or -1 if the loop is empty. |
| `func (l Lasso) States() []int` | States returns prefix followed by loop, i.e. |

</details>

<details>
<summary><code>StateSet</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func BMCUnrollReachSet(k *Kripke, bound int) StateSet` | BMCUnrollReachSet returns the exact set of states reachable within bound steps, computed from the unrolling. |
| `func CTLCheck(k *Kripke, f *CTL) (StateSet, error)` | CTLCheck returns the set of states of k that satisfy the CTL formula f. |
| `func FairCTLCheck(k *Kripke, f *CTL, fc Fairness) (StateSet, error)` | FairCTLCheck returns the set of states of k satisfying the CTL formula f under the fairness condition fc. |
| `func FairEF(k *Kripke, a StateSet, fc Fairness) StateSet` | FairEF returns the set of states satisfying the fair EF operator: a fair path that eventually reaches a. |
| `func FairEG(k *Kripke, a StateSet, fc Fairness) StateSet` | FairEG returns the set of states from which there is a fair path all of whose states lie in a. |
| `func FairEU(k *Kripke, a, b StateSet, fc Fairness) StateSet` | FairEU returns the set of states satisfying the fair existential until: a fair path along which a holds until b holds. |
| `func FairEX(k *Kripke, phi StateSet, fc Fairness) StateSet` | FairEX returns the set of states satisfying the fair existential next operator: some successor satisfies phi and has a fair continuation. |
| `func FairStates(k *Kripke, fc Fairness) StateSet` | FairStates returns the set of states of k from which at least one fair path begins. |
| `func FullStateSet(n int) StateSet` | FullStateSet returns the set containing every state of the universe {0, ..., n-1}. |
| `func NewStateSet(n int) StateSet` | NewStateSet returns an empty set over the universe {0, ..., n-1}. |
| `func PreExists(k *Kripke, target StateSet) StateSet` | PreExists returns the existential predecessor image of target: the set of states having at least one successor in target. |
| `func PreForall(k *Kripke, target StateSet) StateSet` | PreForall returns the universal predecessor image of target: the set of states all of whose successors lie in target. |
| `func SatAF(k *Kripke, a StateSet) StateSet` | SatAF returns the set of states satisfying AF a ≡ ¬EG ¬a. |
| `func SatAG(k *Kripke, a StateSet) StateSet` | SatAG returns the set of states satisfying AG a ≡ ¬EF ¬a. |
| `func SatAR(k *Kripke, a, b StateSet) StateSet` | SatAR returns the set of states satisfying A[a R b] by the greatest fixpoint Z = b ∩ (a ∪ AX Z). |
| `func SatAU(k *Kripke, a, b StateSet) StateSet` | SatAU returns the set of states satisfying A[a U b] by the least fixpoint Z = b ∪ (a ∩ AX Z). |
| `func SatAX(k *Kripke, phi StateSet) StateSet` | SatAX returns the set of states satisfying AX phi. |
| `func SatEF(k *Kripke, a StateSet) StateSet` | SatEF returns the set of states satisfying EF a ≡ E[true U a]. |
| `func SatEG(k *Kripke, a StateSet) StateSet` | SatEG returns the set of states satisfying EG a by computing the greatest fixpoint Z = a ∩ EX Z starting from a. |
| `func SatER(k *Kripke, a, b StateSet) StateSet` | SatER returns the set of states satisfying E[a R b] by the greatest fixpoint Z = b ∩ (a ∪ EX Z). |
| `func SatEU(k *Kripke, a, b StateSet) StateSet` | SatEU returns the set of states satisfying E[a U b] by computing the least fixpoint Z = b ∪ (a ∩ EX Z). |
| `func SatEX(k *Kripke, phi StateSet) StateSet` | SatEX returns the set of states satisfying EX phi, given the set phi of states satisfying the argument. |
| `func StateSetFromSlice(n int, elems []int) StateSet` | StateSetFromSlice builds a set over the universe {0, ..., n-1} containing exactly the listed elements. |
| `func (s StateSet) Add(i int)` | Add inserts state i into the set. |
| `func (s StateSet) Any() (int, bool)` | Any returns some element of the set and true, or (0, false) if the set is empty. |
| `func (s StateSet) Clone() StateSet` | Clone returns an independent copy of the set. |
| `func (s StateSet) Complement() StateSet` | Complement returns the set of states of the universe not in s. |
| `func (s StateSet) Contains(i int) bool` | Contains reports whether state i is a member of the set. |
| `func (s StateSet) Difference(t StateSet) StateSet` | Difference returns the set of states in s but not in t. |
| `func (s StateSet) Elements() []int` | Elements returns the members of the set in ascending order. |
| `func (s StateSet) Equal(t StateSet) bool` | Equal reports whether s and t contain the same elements over the same universe. |
| `func (s StateSet) Intersect(t StateSet) StateSet` | Intersect returns the set of states in both s and t. |
| `func (s StateSet) Intersects(t StateSet) bool` | Intersects reports whether s and t share at least one element. |
| `func (s StateSet) IsEmpty() bool` | IsEmpty reports whether the set has no members. |
| `func (s StateSet) IsSubset(t StateSet) bool` | IsSubset reports whether every element of s is also in t. |
| `func (s StateSet) Len() int` | Len returns the number of states in the set. |
| `func (s StateSet) Remove(i int)` | Remove deletes state i from the set. |
| `func (s StateSet) String() string` | String renders the set as {a, b, c} with elements in ascending order. |
| `func (s StateSet) Union(t StateSet) StateSet` | Union returns the set of states in s or t. |
| `func (s StateSet) Universe() int` | Universe returns the size n of the universe over which the set is defined. |

</details>

<details>
<summary><code>Unrolling</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Unroll(k *Kripke, bound int) Unrolling` | Unroll computes the bounded breadth-first unrolling of k from its initial states to the given depth bound. |
| `func (u Unrolling) Depth() int` | Depth returns the number of populated layers minus one, i.e. |
| `func (u Unrolling) ReachableWithin() StateSet` | ReachableWithin returns the union of all layers of the unrolling: the states reachable within Bound steps of an initial state. |

</details>

### Variables

`ErrState`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/modelchecking).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
