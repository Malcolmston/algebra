package automata

import (
	"errors"
	"fmt"
	"math/big"
	"sort"
)

// DFA is a deterministic finite automaton over an alphabet of runes. States are
// identified by the integers 0..NumStates-1. The transition function may be
// partial: a symbol with no outgoing edge from a state is treated as a
// transition to an implicit non-accepting dead state, so that a run simply
// rejects. Use Complete to materialise the dead state explicitly when a total
// transition function is required.
type DFA struct {
	// NumStates is the number of states; valid states are 0..NumStates-1.
	NumStates int
	// Alphabet is the sorted, de-duplicated set of input symbols.
	Alphabet []rune
	// Start is the initial state.
	Start int
	// Accept holds the set of accepting (final) states.
	Accept StateSet
	// trans[q][a] gives the target of the edge from q on symbol a.
	trans map[int]map[rune]int
}

// NewDFA constructs a DFA with the given number of states, alphabet and start
// state, and no transitions or accepting states.
func NewDFA(numStates int, alphabet []rune, start int) *DFA {
	return &DFA{
		NumStates: numStates,
		Alphabet:  sortedRunes(alphabet),
		Start:     start,
		Accept:    NewStateSet(),
		trans:     make(map[int]map[rune]int),
	}
}

// SetTransition sets the edge from state from on symbol sym to state to. The
// symbol is added to the alphabet if not already present.
func (d *DFA) SetTransition(from int, sym rune, to int) {
	if d.trans[from] == nil {
		d.trans[from] = make(map[rune]int)
	}
	d.trans[from][sym] = to
	if !containsRune(d.Alphabet, sym) {
		d.Alphabet = mergeRunes(d.Alphabet, []rune{sym})
	}
}

// SetAccept marks state as accepting when accepting is true, or non-accepting
// when it is false.
func (d *DFA) SetAccept(state int, accepting bool) {
	if accepting {
		d.Accept[state] = true
	} else {
		delete(d.Accept, state)
	}
}

// AddAccept marks each of the given states as accepting.
func (d *DFA) AddAccept(states ...int) {
	for _, q := range states {
		d.Accept[q] = true
	}
}

// Transition returns the target state of the edge from state from on symbol
// sym, and whether such an edge exists.
func (d *DFA) Transition(from int, sym rune) (int, bool) {
	if m := d.trans[from]; m != nil {
		to, ok := m[sym]
		return to, ok
	}
	return 0, false
}

// Step is an alias for Transition, advancing one symbol from a single state.
func (d *DFA) Step(from int, sym rune) (int, bool) {
	return d.Transition(from, sym)
}

// Run feeds input to the DFA starting from the start state and returns the
// final state reached together with a boolean that is false if the run got
// stuck on a missing transition.
func (d *DFA) Run(input string) (int, bool) {
	q := d.Start
	for _, r := range input {
		next, ok := d.Transition(q, r)
		if !ok {
			return q, false
		}
		q = next
	}
	return q, true
}

// Accepts reports whether the DFA accepts the input string.
func (d *DFA) Accepts(input string) bool {
	q, ok := d.Run(input)
	return ok && d.Accept[q]
}

// IsAccepting reports whether the given state is accepting.
func (d *DFA) IsAccepting(state int) bool {
	return d.Accept[state]
}

// States returns all state indices in ascending order.
func (d *DFA) States() []int {
	out := make([]int, d.NumStates)
	for i := range out {
		out[i] = i
	}
	return out
}

// AcceptingStates returns the accepting states in ascending order.
func (d *DFA) AcceptingStates() []int {
	return d.Accept.Sorted()
}

// Symbols returns the alphabet in ascending order.
func (d *DFA) Symbols() []rune {
	return append([]rune{}, d.Alphabet...)
}

// NumTransitions returns the number of defined edges in the DFA.
func (d *DFA) NumTransitions() int {
	n := 0
	for _, m := range d.trans {
		n += len(m)
	}
	return n
}

// IsComplete reports whether the transition function is total, i.e. every state
// has an outgoing edge for every alphabet symbol.
func (d *DFA) IsComplete() bool {
	for q := 0; q < d.NumStates; q++ {
		for _, a := range d.Alphabet {
			if _, ok := d.Transition(q, a); !ok {
				return false
			}
		}
	}
	return true
}

// Complete returns an equivalent DFA whose transition function is total. If the
// original is already complete a clone is returned; otherwise a single dead
// (non-accepting, self-looping) state is added to absorb missing transitions.
func (d *DFA) Complete() *DFA {
	if d.IsComplete() {
		return d.Clone()
	}
	out := d.Clone()
	dead := out.NumStates
	out.NumStates++
	for q := 0; q < out.NumStates; q++ {
		for _, a := range out.Alphabet {
			if _, ok := out.Transition(q, a); !ok {
				out.SetTransition(q, a, dead)
			}
		}
	}
	return out
}

// Clone returns a deep copy of the DFA.
func (d *DFA) Clone() *DFA {
	c := NewDFA(d.NumStates, d.Alphabet, d.Start)
	for q, m := range d.trans {
		for a, to := range m {
			c.SetTransition(q, a, to)
		}
	}
	for q := range d.Accept {
		c.Accept[q] = true
	}
	return c
}

// Reachable returns the set of states reachable from the start state by any
// sequence of defined transitions.
func (d *DFA) Reachable() StateSet {
	seen := NewStateSet(d.Start)
	stack := []int{d.Start}
	for len(stack) > 0 {
		q := stack[len(stack)-1]
		stack = stack[:len(stack)-1]
		for _, a := range d.Alphabet {
			if to, ok := d.Transition(q, a); ok && !seen[to] {
				seen[to] = true
				stack = append(stack, to)
			}
		}
	}
	return seen
}

// Productive returns the set of states from which some accepting state is
// reachable. Such states are also called live or co-reachable states.
func (d *DFA) Productive() StateSet {
	// Build reverse adjacency.
	rev := make(map[int][]int)
	for q, m := range d.trans {
		for _, to := range m {
			rev[to] = append(rev[to], q)
		}
	}
	live := d.Accept.Clone()
	stack := d.Accept.Sorted()
	for len(stack) > 0 {
		q := stack[len(stack)-1]
		stack = stack[:len(stack)-1]
		for _, p := range rev[q] {
			if !live[p] {
				live[p] = true
				stack = append(stack, p)
			}
		}
	}
	return live
}

// renumber builds a DFA restricted to the given keep set, renumbering states to
// a compact 0-based range. If the start state is not kept, an empty-language
// DFA with a single dead state is returned.
func (d *DFA) renumber(keep StateSet) *DFA {
	if !keep[d.Start] {
		out := NewDFA(1, d.Alphabet, 0)
		return out
	}
	order := keep.Sorted()
	idx := make(map[int]int, len(order))
	for i, q := range order {
		idx[q] = i
	}
	out := NewDFA(len(order), d.Alphabet, idx[d.Start])
	for _, q := range order {
		for _, a := range d.Alphabet {
			if to, ok := d.Transition(q, a); ok && keep[to] {
				out.SetTransition(idx[q], a, idx[to])
			}
		}
		if d.Accept[q] {
			out.Accept[idx[q]] = true
		}
	}
	return out
}

// RemoveUnreachable returns an equivalent DFA with all states unreachable from
// the start removed and the remaining states renumbered compactly.
func (d *DFA) RemoveUnreachable() *DFA {
	return d.renumber(d.Reachable())
}

// RemoveDead returns an equivalent DFA with all non-productive (dead) states
// removed. The start state is always retained.
func (d *DFA) RemoveDead() *DFA {
	live := d.Productive()
	live[d.Start] = true
	return d.renumber(live)
}

// Trim returns an equivalent DFA containing only states that are both reachable
// from the start and productive (able to reach an accepting state).
func (d *DFA) Trim() *DFA {
	keep := d.Reachable().Intersect(d.Productive())
	keep[d.Start] = true
	return d.renumber(keep)
}

// Complement returns a DFA accepting exactly the strings over the alphabet that
// this DFA rejects. The automaton is completed first so that stuck runs become
// accepting runs in the complement.
func (d *DFA) Complement() *DFA {
	c := d.Complete()
	out := c.Clone()
	out.Accept = NewStateSet()
	for q := 0; q < out.NumStates; q++ {
		if !c.Accept[q] {
			out.Accept[q] = true
		}
	}
	return out
}

// ToNFA returns an NFA that accepts the same language as the DFA.
func (d *DFA) ToNFA() *NFA {
	n := NewNFA(d.NumStates, d.Alphabet, d.Start)
	for q, m := range d.trans {
		for a, to := range m {
			n.AddTransition(q, a, to)
		}
	}
	for q := range d.Accept {
		n.Accept[q] = true
	}
	return n
}

// IsEmptyLanguage reports whether the DFA accepts no strings at all.
func (d *DFA) IsEmptyLanguage() bool {
	reach := d.Reachable()
	for q := range reach {
		if d.Accept[q] {
			return false
		}
	}
	return true
}

// IsFiniteLanguage reports whether the language accepted by the DFA is finite.
// The language is infinite exactly when a cycle exists among states that are
// both reachable and productive.
func (d *DFA) IsFiniteLanguage() bool {
	keep := d.Reachable().Intersect(d.Productive())
	// Detect a cycle within the kept subgraph via DFS colouring.
	const (
		white = 0
		gray  = 1
		black = 2
	)
	color := make(map[int]int)
	var dfs func(q int) bool
	dfs = func(q int) bool {
		color[q] = gray
		for _, a := range d.Alphabet {
			if to, ok := d.Transition(q, a); ok && keep[to] {
				switch color[to] {
				case gray:
					return true
				case white:
					if dfs(to) {
						return true
					}
				}
			}
		}
		color[q] = black
		return false
	}
	for q := range keep {
		if color[q] == white {
			if dfs(q) {
				return false
			}
		}
	}
	return true
}

// ShortestAccepted returns a shortest string accepted by the DFA and true, or
// the empty string and false if the language is empty. When several shortest
// words exist the lexicographically smallest by symbol order is returned.
func (d *DFA) ShortestAccepted() (string, bool) {
	type node struct {
		state int
		word  []rune
	}
	visited := NewStateSet(d.Start)
	queue := []node{{d.Start, nil}}
	for len(queue) > 0 {
		cur := queue[0]
		queue = queue[1:]
		if d.Accept[cur.state] {
			return string(cur.word), true
		}
		for _, a := range d.Alphabet {
			if to, ok := d.Transition(cur.state, a); ok && !visited[to] {
				visited[to] = true
				w := append(append([]rune{}, cur.word...), a)
				queue = append(queue, node{to, w})
			}
		}
	}
	return "", false
}

// AcceptedWords returns all accepted strings of length at most maxLen, in
// shortlex order (shorter first, ties broken by symbol order).
func (d *DFA) AcceptedWords(maxLen int) []string {
	var out []string
	var rec func(q int, depth int, word []rune)
	rec = func(q int, depth int, word []rune) {
		if d.Accept[q] {
			out = append(out, string(word))
		}
		if depth == maxLen {
			return
		}
		for _, a := range d.Alphabet {
			if to, ok := d.Transition(q, a); ok {
				rec(to, depth+1, append(word, a))
			}
		}
	}
	rec(d.Start, 0, nil)
	sort.SliceStable(out, func(i, j int) bool {
		if len(out[i]) != len(out[j]) {
			return len(out[i]) < len(out[j])
		}
		return out[i] < out[j]
	})
	return out
}

// CountAcceptedWords returns the exact number of accepted strings of the given
// length as an arbitrary-precision integer, using dynamic programming over path
// counts.
func (d *DFA) CountAcceptedWords(length int) *big.Int {
	counts := make([]*big.Int, d.NumStates)
	for i := range counts {
		counts[i] = big.NewInt(0)
	}
	counts[d.Start] = big.NewInt(1)
	for step := 0; step < length; step++ {
		next := make([]*big.Int, d.NumStates)
		for i := range next {
			next[i] = big.NewInt(0)
		}
		for q := 0; q < d.NumStates; q++ {
			if counts[q].Sign() == 0 {
				continue
			}
			for _, a := range d.Alphabet {
				if to, ok := d.Transition(q, a); ok {
					next[to].Add(next[to], counts[q])
				}
			}
		}
		counts = next
	}
	total := big.NewInt(0)
	for q := range d.Accept {
		total.Add(total, counts[q])
	}
	return total
}

// Validate checks the structural well-formedness of the DFA: valid start state,
// transition endpoints and accepting states all within range.
func (d *DFA) Validate() error {
	if d.NumStates <= 0 {
		return errors.New("automata: DFA must have at least one state")
	}
	if d.Start < 0 || d.Start >= d.NumStates {
		return fmt.Errorf("automata: start state %d out of range", d.Start)
	}
	for q, m := range d.trans {
		if q < 0 || q >= d.NumStates {
			return fmt.Errorf("automata: transition source %d out of range", q)
		}
		for a, to := range m {
			if to < 0 || to >= d.NumStates {
				return fmt.Errorf("automata: transition target %d out of range on %q", to, a)
			}
		}
	}
	for q := range d.Accept {
		if q < 0 || q >= d.NumStates {
			return fmt.Errorf("automata: accepting state %d out of range", q)
		}
	}
	return nil
}

// containsRune reports whether rs contains r.
func containsRune(rs []rune, r rune) bool {
	for _, x := range rs {
		if x == r {
			return true
		}
	}
	return false
}
