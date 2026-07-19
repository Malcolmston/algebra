package automata

import (
	"errors"
	"fmt"
)

// NFA is a nondeterministic finite automaton with epsilon transitions over an
// alphabet of runes. States are the integers 0..NumStates-1. Symbol edges and
// epsilon (empty-string) edges are stored separately.
type NFA struct {
	// NumStates is the number of states; valid states are 0..NumStates-1.
	NumStates int
	// Alphabet is the sorted, de-duplicated set of input symbols.
	Alphabet []rune
	// Start is the initial state.
	Start int
	// Accept holds the set of accepting states.
	Accept StateSet
	// trans[q][a] gives the set of targets of edges from q on symbol a.
	trans map[int]map[rune][]int
	// eps[q] gives the epsilon-successors of q.
	eps map[int][]int
}

// NewNFA constructs an NFA with the given number of states, alphabet and start
// state, and no transitions or accepting states.
func NewNFA(numStates int, alphabet []rune, start int) *NFA {
	return &NFA{
		NumStates: numStates,
		Alphabet:  sortedRunes(alphabet),
		Start:     start,
		Accept:    NewStateSet(),
		trans:     make(map[int]map[rune][]int),
		eps:       make(map[int][]int),
	}
}

// AddTransition adds an edge from state from on symbol sym to state to. The
// symbol is added to the alphabet if not already present.
func (n *NFA) AddTransition(from int, sym rune, to int) {
	if n.trans[from] == nil {
		n.trans[from] = make(map[rune][]int)
	}
	if !containsInt(n.trans[from][sym], to) {
		n.trans[from][sym] = append(n.trans[from][sym], to)
	}
	if !containsRune(n.Alphabet, sym) {
		n.Alphabet = mergeRunes(n.Alphabet, []rune{sym})
	}
}

// AddEpsilon adds an epsilon (empty-string) edge from state from to state to.
func (n *NFA) AddEpsilon(from, to int) {
	if !containsInt(n.eps[from], to) {
		n.eps[from] = append(n.eps[from], to)
	}
}

// SetAccept marks state as accepting when accepting is true, else clears it.
func (n *NFA) SetAccept(state int, accepting bool) {
	if accepting {
		n.Accept[state] = true
	} else {
		delete(n.Accept, state)
	}
}

// AddAccept marks each of the given states as accepting.
func (n *NFA) AddAccept(states ...int) {
	for _, q := range states {
		n.Accept[q] = true
	}
}

// Targets returns the set of states reachable from state from on symbol sym in
// one step (not including epsilon moves).
func (n *NFA) Targets(from int, sym rune) []int {
	if m := n.trans[from]; m != nil {
		return append([]int{}, m[sym]...)
	}
	return nil
}

// EpsilonSuccessors returns the direct epsilon-successors of state q.
func (n *NFA) EpsilonSuccessors(q int) []int {
	return append([]int{}, n.eps[q]...)
}

// HasEpsilon reports whether the NFA contains any epsilon transition.
func (n *NFA) HasEpsilon() bool {
	for _, tos := range n.eps {
		if len(tos) > 0 {
			return true
		}
	}
	return false
}

// EpsilonClosure returns the set of states reachable from any state in the
// input set using zero or more epsilon transitions.
func (n *NFA) EpsilonClosure(states StateSet) StateSet {
	closure := states.Clone()
	stack := states.Sorted()
	for len(stack) > 0 {
		q := stack[len(stack)-1]
		stack = stack[:len(stack)-1]
		for _, to := range n.eps[q] {
			if !closure[to] {
				closure[to] = true
				stack = append(stack, to)
			}
		}
	}
	return closure
}

// EpsilonClosureState returns the epsilon-closure of a single state.
func (n *NFA) EpsilonClosureState(q int) StateSet {
	return n.EpsilonClosure(NewStateSet(q))
}

// Move returns the set of states reachable from the input set by a single
// symbol edge on sym, without taking epsilon closures.
func (n *NFA) Move(states StateSet, sym rune) StateSet {
	out := NewStateSet()
	for q := range states {
		if m := n.trans[q]; m != nil {
			for _, to := range m[sym] {
				out[to] = true
			}
		}
	}
	return out
}

// Step advances the NFA one symbol from the current set of states, applying the
// epsilon closure after the symbol move. The input set is assumed already
// epsilon-closed.
func (n *NFA) Step(states StateSet, sym rune) StateSet {
	return n.EpsilonClosure(n.Move(states, sym))
}

// StartSet returns the epsilon-closure of the start state, i.e. the set of
// states active before any input is read.
func (n *NFA) StartSet() StateSet {
	return n.EpsilonClosureState(n.Start)
}

// Run simulates the NFA on input and returns the set of active states after the
// whole string has been consumed.
func (n *NFA) Run(input string) StateSet {
	cur := n.StartSet()
	for _, r := range input {
		cur = n.Step(cur, r)
	}
	return cur
}

// Accepts reports whether the NFA accepts input, i.e. whether some accepting
// state is active after consuming it.
func (n *NFA) Accepts(input string) bool {
	final := n.Run(input)
	for q := range final {
		if n.Accept[q] {
			return true
		}
	}
	return false
}

// States returns all state indices in ascending order.
func (n *NFA) States() []int {
	out := make([]int, n.NumStates)
	for i := range out {
		out[i] = i
	}
	return out
}

// AcceptingStates returns the accepting states in ascending order.
func (n *NFA) AcceptingStates() []int {
	return n.Accept.Sorted()
}

// Symbols returns the alphabet in ascending order.
func (n *NFA) Symbols() []rune {
	return append([]rune{}, n.Alphabet...)
}

// Clone returns a deep copy of the NFA.
func (n *NFA) Clone() *NFA {
	c := NewNFA(n.NumStates, n.Alphabet, n.Start)
	for q, m := range n.trans {
		for a, tos := range m {
			for _, to := range tos {
				c.AddTransition(q, a, to)
			}
		}
	}
	for q, tos := range n.eps {
		for _, to := range tos {
			c.AddEpsilon(q, to)
		}
	}
	for q := range n.Accept {
		c.Accept[q] = true
	}
	return c
}

// Reachable returns the set of states reachable from the start state using
// symbol and epsilon transitions.
func (n *NFA) Reachable() StateSet {
	seen := NewStateSet(n.Start)
	stack := []int{n.Start}
	for len(stack) > 0 {
		q := stack[len(stack)-1]
		stack = stack[:len(stack)-1]
		for _, to := range n.eps[q] {
			if !seen[to] {
				seen[to] = true
				stack = append(stack, to)
			}
		}
		if m := n.trans[q]; m != nil {
			for _, tos := range m {
				for _, to := range tos {
					if !seen[to] {
						seen[to] = true
						stack = append(stack, to)
					}
				}
			}
		}
	}
	return seen
}

// NumTransitions returns the number of symbol edges (epsilon edges excluded).
func (n *NFA) NumTransitions() int {
	total := 0
	for _, m := range n.trans {
		for _, tos := range m {
			total += len(tos)
		}
	}
	return total
}

// NumEpsilon returns the number of epsilon edges.
func (n *NFA) NumEpsilon() int {
	total := 0
	for _, tos := range n.eps {
		total += len(tos)
	}
	return total
}

// RemoveEpsilon returns an equivalent NFA without epsilon transitions. Each
// state's epsilon-closure is folded into direct symbol edges, and a state
// becomes accepting when its closure contains an accepting state.
func (n *NFA) RemoveEpsilon() *NFA {
	out := NewNFA(n.NumStates, n.Alphabet, n.Start)
	for q := 0; q < n.NumStates; q++ {
		closure := n.EpsilonClosureState(q)
		// Accepting if any state in closure is accepting.
		for c := range closure {
			if n.Accept[c] {
				out.Accept[q] = true
			}
			// For each symbol edge out of a closure member, target's closure
			// becomes reachable from q on that symbol.
			if m := n.trans[c]; m != nil {
				for a, tos := range m {
					for _, to := range tos {
						out.AddTransition(q, a, to)
					}
				}
			}
		}
	}
	return out
}

// Reverse returns an NFA accepting the reversal of this NFA's language: every
// edge is reversed, a fresh start state epsilon-connects to all former
// accepting states, and the former start state becomes the sole accepting
// state.
func (n *NFA) Reverse() *NFA {
	newStart := n.NumStates
	out := NewNFA(n.NumStates+1, n.Alphabet, newStart)
	for q, m := range n.trans {
		for a, tos := range m {
			for _, to := range tos {
				out.AddTransition(to, a, q)
			}
		}
	}
	for q, tos := range n.eps {
		for _, to := range tos {
			out.AddEpsilon(to, q)
		}
	}
	for q := range n.Accept {
		out.AddEpsilon(newStart, q)
	}
	out.Accept[n.Start] = true
	return out
}

// Validate checks the structural well-formedness of the NFA.
func (n *NFA) Validate() error {
	if n.NumStates <= 0 {
		return errors.New("automata: NFA must have at least one state")
	}
	if n.Start < 0 || n.Start >= n.NumStates {
		return fmt.Errorf("automata: start state %d out of range", n.Start)
	}
	for q, m := range n.trans {
		if q < 0 || q >= n.NumStates {
			return fmt.Errorf("automata: transition source %d out of range", q)
		}
		for a, tos := range m {
			for _, to := range tos {
				if to < 0 || to >= n.NumStates {
					return fmt.Errorf("automata: transition target %d out of range on %q", to, a)
				}
			}
		}
	}
	for q, tos := range n.eps {
		for _, to := range tos {
			if to < 0 || to >= n.NumStates {
				return fmt.Errorf("automata: epsilon target %d out of range from %d", to, q)
			}
		}
	}
	for q := range n.Accept {
		if q < 0 || q >= n.NumStates {
			return fmt.Errorf("automata: accepting state %d out of range", q)
		}
	}
	return nil
}

// containsInt reports whether xs contains v.
func containsInt(xs []int, v int) bool {
	for _, x := range xs {
		if x == v {
			return true
		}
	}
	return false
}
