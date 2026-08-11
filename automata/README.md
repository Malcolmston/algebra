# automata — implements finite automata, regular expressions and related machines of formal-language

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/automata.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/automata)

Package automata implements finite automata, regular expressions and related
machines of formal-language theory in pure Go using only the standard library.

- `DFA` is a deterministic finite automaton with a possibly partial transition
function. It supports acceptance testing, completion, complementation,
trimming of unreachable and dead states, shortest-word and word-enumeration
queries, exact word counting via `DFA.CountAcceptedWords`, and
emptiness/finiteness/universality decisions.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/automata"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func CollectAlphabet(node RegexNode) []rune` | CollectAlphabet returns the sorted set of literal symbols appearing in the syntax tree (wildcards contribute nothing). |
| `func DFAToRegexString(d *DFA) string` | DFAToRegexString returns the state-elimination regular expression for d rendered as a string. |
| `func Disjoint(a, b *DFA) bool` | Disjoint reports whether two DFAs accept no common string. |
| `func Equivalent(a, b *DFA) bool` | Equivalent reports whether two DFAs accept exactly the same language. |
| `func Match(pattern, input string) (bool, error)` | Match is a one-shot helper that reports whether pattern matches the entire input string. |
| `func NFAEquivalent(a, b *NFA) bool` | NFAEquivalent reports whether two NFAs accept the same language, by determinising both and comparing the resulting DFAs. |
| `func PumpingDecomposition(d *DFA) func(w string) (x, y, z string, ok bool)` | PumpingDecomposition splits an accepted string w of length at least PumpingLength into parts x, y, z with w = xyz, \|xy\| ≤ p, \|y\| ≥ 1, such… |
| `func PumpingLength(d *DFA) int` | PumpingLength returns a valid pumping length for the DFA's language: the number of states of the trimmed, minimal automaton. |
| `func Subset(a, b *DFA) bool` | Subset reports whether L(a) is a subset of L(b), i.e. |
| `func Witness(a, b *DFA) (string, bool)` | Witness returns a shortest string on which two DFAs disagree (accepted by exactly one of them) together with true, or the empty string and false if… |

### Types

| Type | What it is |
| --- | --- |
| `AcceptMode` | AcceptMode selects a pushdown automaton's acceptance condition. |
| `Alternate` | Alternate denotes the alternation (union) of two sub-expressions L or R. |
| `AnyChar` | AnyChar denotes the wildcard "." matching any single alphabet symbol. |
| `Concat` | Concat denotes the concatenation of two sub-expressions L then R. |
| `DFA` | DFA is a deterministic finite automaton over an alphabet of runes. |
| `Direction` | Direction is the head-movement direction of a Turing machine or the tape motion of a two-way automaton. |
| `EmptySet` | EmptySet denotes the regular expression matching no strings at all (∅). |
| `EmptyString` | EmptyString denotes the regular expression matching only the empty string (ε). |
| `NFA` | NFA is a nondeterministic finite automaton with epsilon transitions over an alphabet of runes. |
| `Optional` | Optional denotes zero-or-one occurrence of a sub-expression. |
| `PDA` | PDA is a nondeterministic pushdown automaton with epsilon moves over rune input and rune stack alphabets. |
| `PDATransition` | PDATransition is one nondeterministic move of a pushdown automaton: enter NextState and replace the popped stack top with the symbols in Push, where… |
| `Plus` | Plus denotes one-or-more repetitions of a sub-expression. |
| `RegexNode` | RegexNode is a node in a parsed regular-expression syntax tree. |
| `Regexp` | Regexp is a compiled regular expression backed by a DFA, offering fast repeated matching. |
| `Star` | Star denotes the Kleene star of a sub-expression (zero or more repetitions). |
| `StateSet` | StateSet is a set of automaton states represented as a map from state index to presence. |
| `Symbol` | Symbol denotes a regular expression matching the single character R. |
| `TM` | TM is a deterministic single-tape Turing machine. |
| `TMConfig` | TMConfig is an instantaneous description of a Turing-machine computation: the current state, the tape contents and the head position. |
| `TMResult` | TMResult reports the outcome of running a Turing machine. |
| `TMRule` | TMRule is the action of a deterministic Turing-machine transition: the state to enter, the symbol to write, and the direction to move the head. |

<details>
<summary><code>Alternate</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (a Alternate) String() string` | String renders an alternation with parentheses. |

</details>

<details>
<summary><code>AnyChar</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (AnyChar) String() string` | String renders the wildcard as ".". |

</details>

<details>
<summary><code>Concat</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (c Concat) String() string` | String renders a concatenation. |

</details>

<details>
<summary><code>DFA</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Difference(a, b *DFA) *DFA` | Difference returns a DFA accepting L(a) \ L(b). |
| `func Hopcroft(d *DFA) *DFA` | Hopcroft returns the minimal DFA equivalent to d using Hopcroft's partition refinement algorithm, which runs in O(n·\|Σ\|·log n) time. |
| `func Intersection(a, b *DFA) *DFA` | Intersection returns a DFA accepting L(a) ∩ L(b). |
| `func Moore(d *DFA) *DFA` | Moore returns the minimal DFA equivalent to d using Moore's classic O(n²·\|Σ\|) partition-refinement algorithm. |
| `func NFAToDFA(n *NFA) *DFA` | NFAToDFA is a package-level convenience wrapper around SubsetConstruction. |
| `func NewDFA(numStates int, alphabet []rune, start int) *DFA` | NewDFA constructs a DFA with the given number of states, alphabet and start state, and no transitions or accepting states. |
| `func Product(a, b *DFA, accept func(aAcc, bAcc bool) bool) *DFA` | Product builds the synchronous product of two DFAs over the combined alphabet. |
| `func RegexToDFA(pattern string) (*DFA, error)` | RegexToDFA parses pattern, builds an NFA and determinises it to a minimal DFA. |
| `func ReverseDFA(d *DFA) *DFA` | ReverseDFA returns a DFA accepting the reversal of L(d), obtained by reversing the automaton and determinising. |
| `func SubsetConstruction(n *NFA) *DFA` | SubsetConstruction converts an NFA (possibly with epsilon transitions) into an equivalent DFA using the classic powerset construction. |
| `func SymmetricDifference(a, b *DFA) *DFA` | SymmetricDifference returns a DFA accepting the symmetric difference (L(a) \ L(b)) ∪ (L(b) \ L(a)), i.e. |
| `func Union(a, b *DFA) *DFA` | Union returns a DFA accepting L(a) ∪ L(b). |
| `func (d *DFA) AcceptedWords(maxLen int) []string` | AcceptedWords returns all accepted strings of length at most maxLen, in shortlex order (shorter first, ties broken by symbol order). |
| `func (d *DFA) AcceptingStates() []int` | AcceptingStates returns the accepting states in ascending order. |
| `func (d *DFA) Accepts(input string) bool` | Accepts reports whether the DFA accepts the input string. |
| `func (d *DFA) AcceptsAll() bool` | AcceptsAll reports whether the DFA accepts every string over its alphabet, i.e. |
| `func (d *DFA) AddAccept(states ...int)` | AddAccept marks each of the given states as accepting. |
| `func (d *DFA) Clone() *DFA` | Clone returns a deep copy of the DFA. |
| `func (d *DFA) Complement() *DFA` | Complement returns a DFA accepting exactly the strings over the alphabet that this DFA rejects. |
| `func (d *DFA) Complete() *DFA` | Complete returns an equivalent DFA whose transition function is total. |
| `func (d *DFA) Concatenate(other *DFA) *DFA` | Concatenate returns a DFA accepting L(d)·L(other). |
| `func (d *DFA) CountAcceptedWords(length int) *big.Int` | CountAcceptedWords returns the exact number of accepted strings of the given length as an arbitrary-precision integer, using dynamic programming over… |
| `func (d *DFA) Difference(other *DFA) *DFA` | Difference returns a DFA accepting the set difference L(d) \ L(other). |
| `func (d *DFA) EquivalenceClasses() [][]int` | EquivalenceClasses returns the blocks of the Myhill–Nerode equivalence over the reachable, completed DFA: each block is a set of mutually… |
| `func (d *DFA) Intersection(other *DFA) *DFA` | Intersection returns a DFA accepting the intersection of the two languages. |
| `func (d *DFA) IsAccepting(state int) bool` | IsAccepting reports whether the given state is accepting. |
| `func (d *DFA) IsComplete() bool` | IsComplete reports whether the transition function is total, i.e. |
| `func (d *DFA) IsEmptyLanguage() bool` | IsEmptyLanguage reports whether the DFA accepts no strings at all. |
| `func (d *DFA) IsEquivalentTo(other *DFA) bool` | IsEquivalentTo reports whether this DFA is language-equivalent to other. |
| `func (d *DFA) IsFiniteLanguage() bool` | IsFiniteLanguage reports whether the language accepted by the DFA is finite. |
| `func (d *DFA) IsSubsetOf(other *DFA) bool` | IsSubsetOf reports whether this DFA's language is contained in other's. |
| `func (d *DFA) IsUniversal() bool` | IsUniversal is an alias for AcceptsAll. |
| `func (d *DFA) Minimize() *DFA` | Minimize returns the unique (up to renumbering) minimal DFA equivalent to the receiver. |
| `func (d *DFA) MinimizeMoore() *DFA` | MinimizeMoore is a method form of Moore's minimisation algorithm. |
| `func (d *DFA) NumTransitions() int` | NumTransitions returns the number of defined edges in the DFA. |
| `func (d *DFA) Plus() *DFA` | Plus returns a DFA accepting L(d)+. |
| `func (d *DFA) Productive() StateSet` | Productive returns the set of states from which some accepting state is reachable. |
| `func (d *DFA) Reachable() StateSet` | Reachable returns the set of states reachable from the start state by any sequence of defined transitions. |
| `func (d *DFA) RemoveDead() *DFA` | RemoveDead returns an equivalent DFA with all non-productive (dead) states removed. |
| `func (d *DFA) RemoveUnreachable() *DFA` | RemoveUnreachable returns an equivalent DFA with all states unreachable from the start removed and the remaining states renumbered compactly. |
| `func (d *DFA) Reverse() *DFA` | Reverse returns a DFA accepting the reversal of this DFA's language. |
| `func (d *DFA) Run(input string) (int, bool)` | Run feeds input to the DFA starting from the start state and returns the final state reached together with a boolean that is false if the run got… |
| `func (d *DFA) SetAccept(state int, accepting bool)` | SetAccept marks state as accepting when accepting is true, or non-accepting when it is false. |
| `func (d *DFA) SetTransition(from int, sym rune, to int)` | SetTransition sets the edge from state from on symbol sym to state to. |
| `func (d *DFA) ShortestAccepted() (string, bool)` | ShortestAccepted returns a shortest string accepted by the DFA and true, or the empty string and false if the language is empty. |
| `func (d *DFA) Star() *DFA` | Star returns a DFA accepting L(d)*. |
| `func (d *DFA) States() []int` | States returns all state indices in ascending order. |
| `func (d *DFA) Step(from int, sym rune) (int, bool)` | Step is an alias for Transition, advancing one symbol from a single state. |
| `func (d *DFA) Symbols() []rune` | Symbols returns the alphabet in ascending order. |
| `func (d *DFA) SymmetricDifference(other *DFA) *DFA` | SymmetricDifference returns a DFA accepting the symmetric difference of the two languages. |
| `func (d *DFA) ToNFA() *NFA` | ToNFA returns an NFA that accepts the same language as the DFA. |
| `func (d *DFA) Transition(from int, sym rune) (int, bool)` | Transition returns the target state of the edge from state from on symbol sym, and whether such an edge exists. |
| `func (d *DFA) Trim() *DFA` | Trim returns an equivalent DFA containing only states that are both reachable from the start and productive (able to reach an accepting state). |
| `func (d *DFA) Union(other *DFA) *DFA` | Union returns a DFA accepting the union of this DFA's language with other's. |
| `func (d *DFA) Validate() error` | Validate checks the structural well-formedness of the DFA: valid start state, transition endpoints and accepting states all within range. |

</details>

<details>
<summary><code>Direction</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (d Direction) String() string` | String renders a Direction as "L", "S" or "R". |

</details>

<details>
<summary><code>EmptySet</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (EmptySet) String() string` | String renders the empty set as the symbol ∅. |

</details>

<details>
<summary><code>EmptyString</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (EmptyString) String() string` | String renders the empty string as the symbol ε. |

</details>

<details>
<summary><code>NFA</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Concatenate(a, b *NFA) *NFA` | Concatenate returns an NFA accepting the concatenation L(a)·L(b): every string formed by an a-word followed by a b-word. |
| `func EmptyLanguageNFA(alphabet []rune) *NFA` | EmptyLanguageNFA returns an NFA over alphabet accepting no strings. |
| `func EpsilonNFA(alphabet []rune) *NFA` | EpsilonNFA returns an NFA over alphabet accepting only the empty string. |
| `func LiteralNFA(s string) *NFA` | LiteralNFA returns an NFA accepting exactly the string s. |
| `func NewNFA(numStates int, alphabet []rune, start int) *NFA` | NewNFA constructs an NFA with the given number of states, alphabet and start state, and no transitions or accepting states. |
| `func OptionalNFA(a *NFA) *NFA` | OptionalNFA returns an NFA accepting L(a) ∪ {ε}. |
| `func PlusNFA(a *NFA) *NFA` | PlusNFA returns an NFA accepting L(a)+ = L(a)·L(a)*, i.e. |
| `func PowerNFA(a *NFA, n int) *NFA` | PowerNFA returns an NFA accepting exactly n consecutive copies of L(a). |
| `func RegexToNFA(pattern string) (*NFA, error)` | RegexToNFA parses pattern and returns an equivalent NFA via Thompson's construction. |
| `func RegexToNFAWithAlphabet(pattern string, alphabet []rune) (*NFA, error)` | RegexToNFAWithAlphabet is like RegexToNFA but additionally includes the given alphabet symbols, which is required when the pattern uses '.' and the… |
| `func ReverseNFA(a *NFA) *NFA` | ReverseNFA returns an NFA accepting the reversal of L(a). |
| `func StarNFA(a *NFA) *NFA` | StarNFA returns an NFA accepting the Kleene star L(a)*, including the empty string. |
| `func SymbolNFA(sym rune) *NFA` | SymbolNFA returns an NFA accepting exactly the single-character string sym. |
| `func Thompson(node RegexNode, alphabet []rune) *NFA` | Thompson builds an NFA from a regex syntax tree using Thompson's construction. |
| `func UnionNFA(a, b *NFA) *NFA` | UnionNFA returns an NFA accepting L(a) ∪ L(b) via a fresh start state with epsilon edges into both operands. |
| `func (n *NFA) AcceptingStates() []int` | AcceptingStates returns the accepting states in ascending order. |
| `func (n *NFA) Accepts(input string) bool` | Accepts reports whether the NFA accepts input, i.e. |
| `func (n *NFA) AddAccept(states ...int)` | AddAccept marks each of the given states as accepting. |
| `func (n *NFA) AddEpsilon(from, to int)` | AddEpsilon adds an epsilon (empty-string) edge from state from to state to. |
| `func (n *NFA) AddTransition(from int, sym rune, to int)` | AddTransition adds an edge from state from on symbol sym to state to. |
| `func (n *NFA) Clone() *NFA` | Clone returns a deep copy of the NFA. |
| `func (n *NFA) Concat(other *NFA) *NFA` | Concat is the method form of Concatenate. |
| `func (n *NFA) Determinize() *DFA` | Determinize is an alias for ToDFA that additionally removes unreachable states (the powerset construction already avoids them) and is provided for… |
| `func (n *NFA) EpsilonClosure(states StateSet) StateSet` | EpsilonClosure returns the set of states reachable from any state in the input set using zero or more epsilon transitions. |
| `func (n *NFA) EpsilonClosureState(q int) StateSet` | EpsilonClosureState returns the epsilon-closure of a single state. |
| `func (n *NFA) EpsilonSuccessors(q int) []int` | EpsilonSuccessors returns the direct epsilon-successors of state q. |
| `func (n *NFA) HasEpsilon() bool` | HasEpsilon reports whether the NFA contains any epsilon transition. |
| `func (n *NFA) Move(states StateSet, sym rune) StateSet` | Move returns the set of states reachable from the input set by a single symbol edge on sym, without taking epsilon closures. |
| `func (n *NFA) NumEpsilon() int` | NumEpsilon returns the number of epsilon edges. |
| `func (n *NFA) NumTransitions() int` | NumTransitions returns the number of symbol edges (epsilon edges excluded). |
| `func (n *NFA) Optional() *NFA` | Optional is the method form of OptionalNFA. |
| `func (n *NFA) Plus() *NFA` | Plus is the method form of PlusNFA. |
| `func (n *NFA) Power(k int) *NFA` | Power is the method form of PowerNFA. |
| `func (n *NFA) Reachable() StateSet` | Reachable returns the set of states reachable from the start state using symbol and epsilon transitions. |
| `func (n *NFA) RemoveEpsilon() *NFA` | RemoveEpsilon returns an equivalent NFA without epsilon transitions. |
| `func (n *NFA) Reverse() *NFA` | Reverse returns an NFA accepting the reversal of this NFA's language: every edge is reversed, a fresh start state epsilon-connects to all former… |
| `func (n *NFA) Run(input string) StateSet` | Run simulates the NFA on input and returns the set of active states after the whole string has been consumed. |
| `func (n *NFA) SetAccept(state int, accepting bool)` | SetAccept marks state as accepting when accepting is true, else clears it. |
| `func (n *NFA) Star() *NFA` | Star is the method form of StarNFA. |
| `func (n *NFA) StartSet() StateSet` | StartSet returns the epsilon-closure of the start state, i.e. |
| `func (n *NFA) States() []int` | States returns all state indices in ascending order. |
| `func (n *NFA) Step(states StateSet, sym rune) StateSet` | Step advances the NFA one symbol from the current set of states, applying the epsilon closure after the symbol move. |
| `func (n *NFA) Symbols() []rune` | Symbols returns the alphabet in ascending order. |
| `func (n *NFA) Targets(from int, sym rune) []int` | Targets returns the set of states reachable from state from on symbol sym in one step (not including epsilon moves). |
| `func (n *NFA) ToDFA() *DFA` | ToDFA converts the NFA to an equivalent DFA via subset construction. |
| `func (n *NFA) Union(other *NFA) *NFA` | Union is the method form of UnionNFA. |
| `func (n *NFA) Validate() error` | Validate checks the structural well-formedness of the NFA. |

</details>

<details>
<summary><code>Optional</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (o Optional) String() string` | String renders an optional, parenthesising compound operands. |

</details>

<details>
<summary><code>PDA</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewPDA(start string, initialStack rune, mode AcceptMode) *PDA` | NewPDA constructs a pushdown automaton with the given start state, initial stack symbol and acceptance mode. |
| `func (p *PDA) Accepts(input string) bool` | Accepts reports whether the PDA accepts input under its acceptance mode. |
| `func (p *PDA) AddAccept(states ...string)` | AddAccept marks each of the given states as accepting. |
| `func (p *PDA) AddEpsilonTransition(state string, top rune, next string, push ...rune)` | AddEpsilonTransition adds an epsilon move (consuming no input) from state, popping top, entering next and pushing the given symbols. |
| `func (p *PDA) AddTransition(state string, in rune, top rune, next string, push ...rune)` | AddTransition adds a move on input symbol in (popping top): from state, read in and pop top, go to next and push the given symbols. |

</details>

<details>
<summary><code>Plus</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (p Plus) String() string` | String renders a plus, parenthesising compound operands. |

</details>

<details>
<summary><code>RegexNode</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Alt(nodes ...RegexNode) RegexNode` | Alt builds the left-associative alternation of the given nodes. |
| `func DFAToRegex(d *DFA) RegexNode` | DFAToRegex converts a DFA into an equivalent regular expression using the state-elimination method. |
| `func Dot() RegexNode` | Dot builds an AnyChar (.) wildcard node. |
| `func Empty() RegexNode` | Empty builds an EmptySet (∅) node. |
| `func Eps() RegexNode` | Eps builds an EmptyString (ε) node. |
| `func Kleene(x RegexNode) RegexNode` | Kleene builds a Star node over x. |
| `func Lit(r rune) RegexNode` | Lit builds a Symbol node matching the single rune r. |
| `func LiteralNode(s string) RegexNode` | LiteralNode builds a syntax tree matching exactly the string s (a concatenation of its runes; ε when s is empty). |
| `func MustParseRegex(pattern string) RegexNode` | MustParseRegex is like ParseRegex but panics if the pattern is invalid. |
| `func NFAToRegex(n *NFA) RegexNode` | NFAToRegex converts an NFA to an equivalent regular expression by first determinising it and then applying state elimination. |
| `func OneOrMore(x RegexNode) RegexNode` | OneOrMore builds a Plus node over x. |
| `func ParseRegex(pattern string) (RegexNode, error)` | ParseRegex parses a regular expression written in the supported syntax and returns its syntax tree. |
| `func Seq(nodes ...RegexNode) RegexNode` | Seq builds the left-associative concatenation of the given nodes. |
| `func ZeroOrOne(x RegexNode) RegexNode` | ZeroOrOne builds an Optional node over x. |

</details>

<details>
<summary><code>Regexp</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Compile(pattern string) (*Regexp, error)` | Compile parses pattern and returns a compiled Regexp whose alphabet is the set of literal symbols in the pattern. |
| `func CompileWithAlphabet(pattern string, alphabet []rune) (*Regexp, error)` | CompileWithAlphabet parses pattern using the supplied alphabet for wildcard expansion and returns a compiled Regexp. |
| `func MustCompile(pattern string) *Regexp` | MustCompile is like Compile but panics on a malformed pattern. |
| `func (re *Regexp) DFA() *DFA` | DFA returns the underlying deterministic automaton (a clone, so callers may mutate it freely). |
| `func (re *Regexp) Matches(input string) bool` | Matches reports whether the whole input string is matched by the regular expression. |
| `func (re *Regexp) Pattern() string` | Pattern returns the source pattern the Regexp was compiled from. |

</details>

<details>
<summary><code>Star</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (s Star) String() string` | String renders a Kleene star, parenthesising compound operands. |

</details>

<details>
<summary><code>StateSet</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewStateSet(states ...int) StateSet` | NewStateSet returns a StateSet containing exactly the given states. |
| `func (s StateSet) Add(q int) StateSet` | Add inserts state q into the set and returns the set for chaining. |
| `func (s StateSet) Clone() StateSet` | Clone returns an independent copy of the set. |
| `func (s StateSet) Contains(q int) bool` | Contains reports whether state q is a member of the set. |
| `func (s StateSet) Equal(other StateSet) bool` | Equal reports whether two StateSets contain exactly the same states. |
| `func (s StateSet) Intersect(other StateSet) StateSet` | Intersect returns a new StateSet containing states present in both operands. |
| `func (s StateSet) IsEmpty() bool` | IsEmpty reports whether the set contains no states. |
| `func (s StateSet) Key() string` | Key returns a canonical string key for the set, suitable for use as a map index. |
| `func (s StateSet) Len() int` | Len returns the number of states in the set. |
| `func (s StateSet) Sorted() []int` | Sorted returns the members of the set as an ascending slice of ints. |
| `func (s StateSet) Union(other StateSet) StateSet` | Union returns a new StateSet containing every state in either operand. |

</details>

<details>
<summary><code>Symbol</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (s Symbol) String() string` | String renders a single symbol, escaping regex metacharacters. |

</details>

<details>
<summary><code>TM</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewTM(start, accept, reject string, blank rune) *TM` | NewTM constructs a Turing machine with the given start, accept and reject states and blank symbol. |
| `func (m *TM) Accepts(input string, maxSteps int) (bool, error)` | Accepts reports whether the machine accepts input within maxSteps transitions and returns an error if the step budget was exhausted before halting. |
| `func (m *TM) Compute(input string, maxSteps int) (string, error)` | Compute runs the machine on input and returns the trimmed tape contents when it halts in the accept state, modelling a transducer. |
| `func (m *TM) InitialConfig(input string) TMConfig` | InitialConfig returns the starting configuration for the given input, with the head at the leftmost input cell. |
| `func (m *TM) Rule(state string, read rune) (TMRule, bool)` | Rule returns the transition for (state, read) and whether one is defined. |
| `func (m *TM) Run(input string, maxSteps int) TMResult` | Run simulates the machine on input for up to maxSteps transitions and reports the outcome. |
| `func (m *TM) SetRule(state string, read rune, next string, write rune, move Direction)` | SetRule installs the transition δ(state, read) = (next, write, move). |
| `func (m *TM) Step(c TMConfig) (TMConfig, bool)` | Step applies one transition to c and returns the successor configuration and true, or c unchanged and false if the machine has halted (no applicable… |

</details>

<details>
<summary><code>TMConfig</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (c *TMConfig) Read() rune` | Read returns the symbol under the head, extending the tape view with blanks as needed (without mutating the configuration). |
| `func (c *TMConfig) TapeString() string` | TapeString returns the tape contents as a string with leading and trailing blank symbols trimmed. |

</details>

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/automata).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
