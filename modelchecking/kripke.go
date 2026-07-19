package modelchecking

import (
	"errors"
	"fmt"
	"sort"
	"strings"
)

// ErrState is returned when a state index lies outside the range of a Kripke
// structure.
var ErrState = errors.New("modelchecking: state index out of range")

// Kripke is a finite Kripke structure: a labelled transition system over a set
// of atomic propositions. States are the integers {0, ..., NumStates()-1}. Each
// state carries a set of atomic propositions that hold there (its label), a set
// of successor states, and a flag marking it as initial. The transition
// relation need not be total; use [Kripke.MakeTotal] to add self-loops on
// deadlock states as required by the standard branching-time semantics.
type Kripke struct {
	n       int
	labels  []map[string]bool
	succ    [][]int
	initial []bool
	names   []string
}

// NewKripke returns a Kripke structure with n states, no transitions, empty
// labels and no initial states.
func NewKripke(n int) *Kripke {
	if n < 0 {
		n = 0
	}
	k := &Kripke{
		n:       n,
		labels:  make([]map[string]bool, n),
		succ:    make([][]int, n),
		initial: make([]bool, n),
		names:   make([]string, n),
	}
	for i := 0; i < n; i++ {
		k.labels[i] = map[string]bool{}
	}
	return k
}

// NumStates returns the number of states in the structure.
func (k *Kripke) NumStates() int { return k.n }

func (k *Kripke) valid(s int) bool { return s >= 0 && s < k.n }

// AddTransition adds a directed edge from state from to state to. It returns
// [ErrState] if either index is out of range and is idempotent for duplicate
// edges.
func (k *Kripke) AddTransition(from, to int) error {
	if !k.valid(from) || !k.valid(to) {
		return ErrState
	}
	for _, t := range k.succ[from] {
		if t == to {
			return nil
		}
	}
	k.succ[from] = append(k.succ[from], to)
	sort.Ints(k.succ[from])
	return nil
}

// RemoveTransition removes the edge from -> to if present.
func (k *Kripke) RemoveTransition(from, to int) error {
	if !k.valid(from) || !k.valid(to) {
		return ErrState
	}
	out := k.succ[from][:0]
	for _, t := range k.succ[from] {
		if t != to {
			out = append(out, t)
		}
	}
	k.succ[from] = out
	return nil
}

// AddSelfLoop adds the edge s -> s.
func (k *Kripke) AddSelfLoop(s int) error { return k.AddTransition(s, s) }

// HasTransition reports whether the edge from -> to exists.
func (k *Kripke) HasTransition(from, to int) bool {
	if !k.valid(from) {
		return false
	}
	for _, t := range k.succ[from] {
		if t == to {
			return true
		}
	}
	return false
}

// Successors returns the successors of state s in ascending order. The returned
// slice is a copy and may be modified freely by the caller.
func (k *Kripke) Successors(s int) []int {
	if !k.valid(s) {
		return nil
	}
	return append([]int(nil), k.succ[s]...)
}

// Predecessors returns the states with an edge into s, in ascending order.
func (k *Kripke) Predecessors(s int) []int {
	if !k.valid(s) {
		return nil
	}
	var out []int
	for u := 0; u < k.n; u++ {
		if k.HasTransition(u, s) {
			out = append(out, u)
		}
	}
	return out
}

// OutDegree returns the number of successors of state s.
func (k *Kripke) OutDegree(s int) int {
	if !k.valid(s) {
		return 0
	}
	return len(k.succ[s])
}

// NumTransitions returns the total number of edges in the structure.
func (k *Kripke) NumTransitions() int {
	c := 0
	for s := 0; s < k.n; s++ {
		c += len(k.succ[s])
	}
	return c
}

// SetLabel replaces the label of state s with exactly the given propositions.
func (k *Kripke) SetLabel(s int, props ...string) error {
	if !k.valid(s) {
		return ErrState
	}
	m := make(map[string]bool, len(props))
	for _, p := range props {
		m[p] = true
	}
	k.labels[s] = m
	return nil
}

// AddLabel adds proposition p to the label of state s.
func (k *Kripke) AddLabel(s int, p string) error {
	if !k.valid(s) {
		return ErrState
	}
	k.labels[s][p] = true
	return nil
}

// RemoveLabel removes proposition p from the label of state s.
func (k *Kripke) RemoveLabel(s int, p string) error {
	if !k.valid(s) {
		return ErrState
	}
	delete(k.labels[s], p)
	return nil
}

// Holds reports whether atomic proposition p holds in state s.
func (k *Kripke) Holds(s int, p string) bool {
	if !k.valid(s) {
		return false
	}
	return k.labels[s][p]
}

// Label returns the atomic propositions holding in state s, sorted.
func (k *Kripke) Label(s int) []string {
	if !k.valid(s) {
		return nil
	}
	out := make([]string, 0, len(k.labels[s]))
	for p := range k.labels[s] {
		out = append(out, p)
	}
	sort.Strings(out)
	return out
}

// LabelSet returns the label of state s as a boolean map. The returned map is a
// copy.
func (k *Kripke) LabelSet(s int) map[string]bool {
	m := map[string]bool{}
	if !k.valid(s) {
		return m
	}
	for p := range k.labels[s] {
		m[p] = true
	}
	return m
}

// SetInitial marks the given states as initial (in addition to any already
// marked). Out-of-range indices are ignored.
func (k *Kripke) SetInitial(states ...int) {
	for _, s := range states {
		if k.valid(s) {
			k.initial[s] = true
		}
	}
}

// ClearInitial unmarks state s as initial.
func (k *Kripke) ClearInitial(s int) {
	if k.valid(s) {
		k.initial[s] = false
	}
}

// IsInitial reports whether state s is an initial state.
func (k *Kripke) IsInitial(s int) bool {
	return k.valid(s) && k.initial[s]
}

// InitialStates returns the initial states in ascending order.
func (k *Kripke) InitialStates() []int {
	var out []int
	for s := 0; s < k.n; s++ {
		if k.initial[s] {
			out = append(out, s)
		}
	}
	return out
}

// SetStateName attaches a human-readable name to state s, used by [Kripke.String]
// and trace rendering.
func (k *Kripke) SetStateName(s int, name string) error {
	if !k.valid(s) {
		return ErrState
	}
	k.names[s] = name
	return nil
}

// StateName returns the name of state s, or a default "sN" if none was set.
func (k *Kripke) StateName(s int) string {
	if !k.valid(s) {
		return "?"
	}
	if k.names[s] != "" {
		return k.names[s]
	}
	return fmt.Sprintf("s%d", s)
}

// Propositions returns every atomic proposition mentioned in any state label,
// sorted and de-duplicated.
func (k *Kripke) Propositions() []string {
	set := map[string]bool{}
	for s := 0; s < k.n; s++ {
		for p := range k.labels[s] {
			set[p] = true
		}
	}
	out := make([]string, 0, len(set))
	for p := range set {
		out = append(out, p)
	}
	sort.Strings(out)
	return out
}

// IsTotal reports whether every state has at least one successor.
func (k *Kripke) IsTotal() bool {
	for s := 0; s < k.n; s++ {
		if len(k.succ[s]) == 0 {
			return false
		}
	}
	return true
}

// Deadlocks returns the states that have no successors.
func (k *Kripke) Deadlocks() []int {
	var out []int
	for s := 0; s < k.n; s++ {
		if len(k.succ[s]) == 0 {
			out = append(out, s)
		}
	}
	return out
}

// MakeTotal adds a self-loop to every deadlock state so that the transition
// relation becomes total, as required by the usual CTL/LTL semantics over
// infinite paths. It returns the number of self-loops added.
func (k *Kripke) MakeTotal() int {
	c := 0
	for s := 0; s < k.n; s++ {
		if len(k.succ[s]) == 0 {
			k.succ[s] = []int{s}
			c++
		}
	}
	return c
}

// LabelStateSet returns the set of states whose label contains proposition p.
func (k *Kripke) LabelStateSet(p string) StateSet {
	set := NewStateSet(k.n)
	for s := 0; s < k.n; s++ {
		if k.labels[s][p] {
			set.Add(s)
		}
	}
	return set
}

// Reachable returns the set of states reachable from the initial states by
// following transitions.
func (k *Kripke) Reachable() StateSet {
	return k.ReachableFrom(k.InitialStates()...)
}

// ReachableFrom returns the set of states reachable from the given sources.
func (k *Kripke) ReachableFrom(sources ...int) StateSet {
	seen := NewStateSet(k.n)
	stack := []int{}
	for _, s := range sources {
		if k.valid(s) && !seen.Contains(s) {
			seen.Add(s)
			stack = append(stack, s)
		}
	}
	for len(stack) > 0 {
		u := stack[len(stack)-1]
		stack = stack[:len(stack)-1]
		for _, v := range k.succ[u] {
			if !seen.Contains(v) {
				seen.Add(v)
				stack = append(stack, v)
			}
		}
	}
	return seen
}

// Reverse returns a new Kripke structure with the same states, labels and
// initial marks but with every transition reversed. It is used by the CTL
// backward-image computations conceptually; the labels are copied.
func (k *Kripke) Reverse() *Kripke {
	r := NewKripke(k.n)
	for s := 0; s < k.n; s++ {
		r.labels[s] = map[string]bool{}
		for p := range k.labels[s] {
			r.labels[s][p] = true
		}
		r.initial[s] = k.initial[s]
		r.names[s] = k.names[s]
		for _, t := range k.succ[s] {
			r.succ[t] = append(r.succ[t], s)
		}
	}
	for s := 0; s < r.n; s++ {
		sort.Ints(r.succ[s])
	}
	return r
}

// Clone returns a deep copy of the structure.
func (k *Kripke) Clone() *Kripke {
	c := NewKripke(k.n)
	for s := 0; s < k.n; s++ {
		for p := range k.labels[s] {
			c.labels[s][p] = true
		}
		c.succ[s] = append([]int(nil), k.succ[s]...)
		c.initial[s] = k.initial[s]
		c.names[s] = k.names[s]
	}
	return c
}

// Validate checks structural invariants: every successor index is in range.
// It returns nil for a well-formed structure.
func (k *Kripke) Validate() error {
	for s := 0; s < k.n; s++ {
		for _, t := range k.succ[s] {
			if !k.valid(t) {
				return fmt.Errorf("modelchecking: state %d has out-of-range successor %d", s, t)
			}
		}
	}
	return nil
}

// String renders a compact multi-line description of the structure.
func (k *Kripke) String() string {
	var b strings.Builder
	fmt.Fprintf(&b, "Kripke(%d states):\n", k.n)
	for s := 0; s < k.n; s++ {
		init := " "
		if k.initial[s] {
			init = "*"
		}
		fmt.Fprintf(&b, " %s%s {%s} ->", init, k.StateName(s), strings.Join(k.Label(s), ","))
		names := make([]string, len(k.succ[s]))
		for i, t := range k.succ[s] {
			names[i] = k.StateName(t)
		}
		fmt.Fprintf(&b, " [%s]\n", strings.Join(names, ","))
	}
	return b.String()
}
