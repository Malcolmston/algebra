package modelchecking

import (
	"fmt"
	"sort"
	"strings"
)

// Guard is a conjunction of literals labelling a Büchi automaton edge. Pos lists
// atomic propositions that must hold and Neg lists propositions that must not
// hold for the edge to be taken. An empty guard (no positive and no negative
// literals) is the constant true and matches every letter.
type Guard struct {
	Pos []string
	Neg []string
}

// TrueGuard returns the guard that matches every letter.
func TrueGuard() Guard { return Guard{} }

// NewGuard returns a guard requiring every proposition in pos to hold and every
// proposition in neg to fail. The inputs are copied, sorted and de-duplicated.
func NewGuard(pos, neg []string) Guard {
	return Guard{Pos: sortedUnique(pos), Neg: sortedUnique(neg)}
}

// IsTrue reports whether the guard imposes no constraints.
func (g Guard) IsTrue() bool { return len(g.Pos) == 0 && len(g.Neg) == 0 }

// Satisfies reports whether a state label (as a set of true propositions)
// satisfies the guard: all positive literals present, all negative literals
// absent.
func (g Guard) Satisfies(label map[string]bool) bool {
	for _, p := range g.Pos {
		if !label[p] {
			return false
		}
	}
	for _, p := range g.Neg {
		if label[p] {
			return false
		}
	}
	return true
}

// Contradictory reports whether the guard requires a proposition to be both true
// and false, so that it matches no letter.
func (g Guard) Contradictory() bool {
	set := map[string]bool{}
	for _, p := range g.Pos {
		set[p] = true
	}
	for _, p := range g.Neg {
		if set[p] {
			return true
		}
	}
	return false
}

// String renders the guard as a conjunction such as "a & !b", or "true" when
// unconstrained.
func (g Guard) String() string {
	if g.IsTrue() {
		return "true"
	}
	parts := make([]string, 0, len(g.Pos)+len(g.Neg))
	for _, p := range g.Pos {
		parts = append(parts, p)
	}
	for _, p := range g.Neg {
		parts = append(parts, "!"+p)
	}
	return strings.Join(parts, " & ")
}

// BuchiEdge is a labelled transition of a [Buchi] or [GenBuchi] automaton.
type BuchiEdge struct {
	To    int
	Guard Guard
}

// Buchi is a nondeterministic Büchi automaton over the alphabet 2^AP with a
// single acceptance set. A run is accepting if it visits some accepting state
// infinitely often. Edges are labelled by [Guard] conjunctions of literals.
type Buchi struct {
	n       int
	edges   [][]BuchiEdge
	initial []bool
	accept  []bool
	props   []string
}

// NewBuchi returns a Büchi automaton with n states, no edges, no initial states
// and no accepting states.
func NewBuchi(n int) *Buchi {
	if n < 0 {
		n = 0
	}
	return &Buchi{
		n:       n,
		edges:   make([][]BuchiEdge, n),
		initial: make([]bool, n),
		accept:  make([]bool, n),
	}
}

// NumStates returns the number of states.
func (b *Buchi) NumStates() int { return b.n }

// AddEdge adds a guarded transition from -> to. Out-of-range endpoints are a
// no-op.
func (b *Buchi) AddEdge(from, to int, g Guard) {
	if from < 0 || from >= b.n || to < 0 || to >= b.n {
		return
	}
	b.edges[from] = append(b.edges[from], BuchiEdge{To: to, Guard: g})
}

// Edges returns the outgoing edges of state s.
func (b *Buchi) Edges(s int) []BuchiEdge {
	if s < 0 || s >= b.n {
		return nil
	}
	return append([]BuchiEdge(nil), b.edges[s]...)
}

// Successors returns the distinct successor states of s in ascending order.
func (b *Buchi) Successors(s int) []int {
	if s < 0 || s >= b.n {
		return nil
	}
	seen := map[int]bool{}
	var out []int
	for _, e := range b.edges[s] {
		if !seen[e.To] {
			seen[e.To] = true
			out = append(out, e.To)
		}
	}
	sort.Ints(out)
	return out
}

// SetInitial marks the given states as initial.
func (b *Buchi) SetInitial(states ...int) {
	for _, s := range states {
		if s >= 0 && s < b.n {
			b.initial[s] = true
		}
	}
}

// SetAccepting marks the given states as accepting.
func (b *Buchi) SetAccepting(states ...int) {
	for _, s := range states {
		if s >= 0 && s < b.n {
			b.accept[s] = true
		}
	}
}

// IsInitial reports whether s is an initial state.
func (b *Buchi) IsInitial(s int) bool { return s >= 0 && s < b.n && b.initial[s] }

// IsAccepting reports whether s is an accepting state.
func (b *Buchi) IsAccepting(s int) bool { return s >= 0 && s < b.n && b.accept[s] }

// InitialStates returns the initial states in ascending order.
func (b *Buchi) InitialStates() []int {
	var out []int
	for s := 0; s < b.n; s++ {
		if b.initial[s] {
			out = append(out, s)
		}
	}
	return out
}

// AcceptingStates returns the accepting states in ascending order.
func (b *Buchi) AcceptingStates() []int {
	var out []int
	for s := 0; s < b.n; s++ {
		if b.accept[s] {
			out = append(out, s)
		}
	}
	return out
}

// NumEdges returns the total number of edges.
func (b *Buchi) NumEdges() int {
	c := 0
	for s := 0; s < b.n; s++ {
		c += len(b.edges[s])
	}
	return c
}

// SetPropositions records the atomic proposition alphabet, used for reporting.
func (b *Buchi) SetPropositions(props []string) { b.props = sortedUnique(props) }

// Propositions returns the recorded atomic proposition alphabet.
func (b *Buchi) Propositions() []string { return append([]string(nil), b.props...) }

// String renders a multi-line description of the automaton.
func (b *Buchi) String() string {
	var sb strings.Builder
	fmt.Fprintf(&sb, "Buchi(%d states):\n", b.n)
	for s := 0; s < b.n; s++ {
		mark := " "
		if b.initial[s] {
			mark = ">"
		}
		acc := " "
		if b.accept[s] {
			acc = "@"
		}
		fmt.Fprintf(&sb, " %s%sq%d:", mark, acc, s)
		for _, e := range b.edges[s] {
			fmt.Fprintf(&sb, " -%s->q%d", e.Guard.String(), e.To)
		}
		sb.WriteByte('\n')
	}
	return sb.String()
}

// GenBuchi is a generalized Büchi automaton: like [Buchi] but with a family of
// acceptance sets. A run is accepting if it visits each acceptance set
// infinitely often. Generalized acceptance is the natural output of the LTL
// tableau construction and is converted to ordinary Büchi acceptance by
// [Degeneralize].
type GenBuchi struct {
	n       int
	edges   [][]BuchiEdge
	initial []bool
	accept  []StateSet
	props   []string
}

// NewGenBuchi returns a generalized Büchi automaton with n states and no
// acceptance sets.
func NewGenBuchi(n int) *GenBuchi {
	if n < 0 {
		n = 0
	}
	return &GenBuchi{
		n:       n,
		edges:   make([][]BuchiEdge, n),
		initial: make([]bool, n),
	}
}

// NumStates returns the number of states.
func (g *GenBuchi) NumStates() int { return g.n }

// AddEdge adds a guarded transition from -> to.
func (g *GenBuchi) AddEdge(from, to int, guard Guard) {
	if from < 0 || from >= g.n || to < 0 || to >= g.n {
		return
	}
	g.edges[from] = append(g.edges[from], BuchiEdge{To: to, Guard: guard})
}

// Edges returns the outgoing edges of state s.
func (g *GenBuchi) Edges(s int) []BuchiEdge {
	if s < 0 || s >= g.n {
		return nil
	}
	return append([]BuchiEdge(nil), g.edges[s]...)
}

// SetInitial marks the given states as initial.
func (g *GenBuchi) SetInitial(states ...int) {
	for _, s := range states {
		if s >= 0 && s < g.n {
			g.initial[s] = true
		}
	}
}

// IsInitial reports whether s is an initial state.
func (g *GenBuchi) IsInitial(s int) bool { return s >= 0 && s < g.n && g.initial[s] }

// InitialStates returns the initial states in ascending order.
func (g *GenBuchi) InitialStates() []int {
	var out []int
	for s := 0; s < g.n; s++ {
		if g.initial[s] {
			out = append(out, s)
		}
	}
	return out
}

// AddAcceptanceSet appends an acceptance condition (a set of states that must be
// visited infinitely often).
func (g *GenBuchi) AddAcceptanceSet(set StateSet) {
	g.accept = append(g.accept, set.Clone())
}

// NumAcceptanceSets returns the number of acceptance sets.
func (g *GenBuchi) NumAcceptanceSets() int { return len(g.accept) }

// AcceptanceSet returns the i-th acceptance set.
func (g *GenBuchi) AcceptanceSet(i int) StateSet {
	if i < 0 || i >= len(g.accept) {
		return NewStateSet(g.n)
	}
	return g.accept[i].Clone()
}

// SetPropositions records the atomic proposition alphabet.
func (g *GenBuchi) SetPropositions(props []string) { g.props = sortedUnique(props) }

// Propositions returns the recorded atomic proposition alphabet.
func (g *GenBuchi) Propositions() []string { return append([]string(nil), g.props...) }

// NumEdges returns the total number of edges.
func (g *GenBuchi) NumEdges() int {
	c := 0
	for s := 0; s < g.n; s++ {
		c += len(g.edges[s])
	}
	return c
}

// Degeneralize converts a generalized Büchi automaton with k acceptance sets
// into an equivalent ordinary [Buchi] automaton using the standard counter
// construction: the state space is duplicated k times and a run advances the
// counter from layer i to layer i+1 (mod k) when it visits the i-th acceptance
// set, with the accepting states placed in layer 0. When k is zero every run is
// accepting, modelled by a single acceptance set equal to all states.
func Degeneralize(g *GenBuchi) *Buchi {
	k := len(g.accept)
	if k == 0 {
		b := NewBuchi(g.n)
		b.props = append([]string(nil), g.props...)
		for s := 0; s < g.n; s++ {
			for _, e := range g.edges[s] {
				b.AddEdge(s, e.To, e.Guard)
			}
			if g.initial[s] {
				b.SetInitial(s)
			}
			b.SetAccepting(s)
		}
		return b
	}
	idx := func(state, layer int) int { return layer*g.n + state }
	b := NewBuchi(g.n * k)
	b.props = append([]string(nil), g.props...)
	for layer := 0; layer < k; layer++ {
		for s := 0; s < g.n; s++ {
			nextLayer := layer
			if g.accept[layer].Contains(s) {
				nextLayer = (layer + 1) % k
			}
			for _, e := range g.edges[s] {
				b.AddEdge(idx(s, layer), idx(e.To, nextLayer), e.Guard)
			}
		}
	}
	for s := 0; s < g.n; s++ {
		if g.initial[s] {
			b.SetInitial(idx(s, 0))
		}
		if g.accept[0].Contains(s) {
			b.SetAccepting(idx(s, 0))
		}
	}
	return b
}
