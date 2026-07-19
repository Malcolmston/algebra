package modelchecking

// Lasso is an ultimately periodic accepting run of a Büchi automaton or the
// witness path of a model-checking counterexample. The run visits the states of
// Prefix once and then repeats the states of Loop forever. There is an edge from
// the last state of Prefix (or from an initial state, if Prefix is empty) to
// Loop[0], and an edge from the last state of Loop back to Loop[0].
type Lasso struct {
	Prefix []int
	Loop   []int
}

// Length returns the number of distinct states listed in the lasso (prefix plus
// loop).
func (l Lasso) Length() int { return len(l.Prefix) + len(l.Loop) }

// LoopStart returns the first state of the periodic part, or -1 if the loop is
// empty.
func (l Lasso) LoopStart() int {
	if len(l.Loop) == 0 {
		return -1
	}
	return l.Loop[0]
}

// States returns prefix followed by loop, i.e. one unrolling of the run.
func (l Lasso) States() []int {
	out := make([]int, 0, len(l.Prefix)+len(l.Loop))
	out = append(out, l.Prefix...)
	out = append(out, l.Loop...)
	return out
}

// buchiSucc returns the successors of s reachable by a non-contradictory edge.
func buchiSucc(b *Buchi, s int) []int {
	var out []int
	seen := map[int]bool{}
	for _, e := range b.edges[s] {
		if e.Guard.Contradictory() {
			continue
		}
		if !seen[e.To] {
			seen[e.To] = true
			out = append(out, e.To)
		}
	}
	return out
}

// ReachableStates returns the set of states reachable from the initial states
// via non-contradictory edges.
func (b *Buchi) ReachableStates() StateSet {
	seen := NewStateSet(b.n)
	var stack []int
	for _, s := range b.InitialStates() {
		if !seen.Contains(s) {
			seen.Add(s)
			stack = append(stack, s)
		}
	}
	for len(stack) > 0 {
		u := stack[len(stack)-1]
		stack = stack[:len(stack)-1]
		for _, v := range buchiSucc(b, u) {
			if !seen.Contains(v) {
				seen.Add(v)
				stack = append(stack, v)
			}
		}
	}
	return seen
}

// bfsPath returns a shortest path (as a slice of states, inclusive of endpoints)
// from any state in sources to target, using non-contradictory edges. It
// returns nil if target is unreachable. A source equal to target yields the
// singleton path.
func bfsPath(b *Buchi, sources []int, target int) []int {
	parent := make([]int, b.n)
	for i := range parent {
		parent[i] = -2
	}
	var queue []int
	for _, s := range sources {
		if s >= 0 && s < b.n && parent[s] == -2 {
			parent[s] = -1
			queue = append(queue, s)
		}
	}
	for _, s := range sources {
		if s == target {
			return []int{target}
		}
	}
	for len(queue) > 0 {
		u := queue[0]
		queue = queue[1:]
		for _, v := range buchiSucc(b, u) {
			if parent[v] == -2 {
				parent[v] = u
				if v == target {
					return reconstruct(parent, target)
				}
				queue = append(queue, v)
			}
		}
	}
	return nil
}

func reconstruct(parent []int, target int) []int {
	var rev []int
	for x := target; x != -1; x = parent[x] {
		rev = append(rev, x)
	}
	for i, j := 0, len(rev)-1; i < j; i, j = i+1, j-1 {
		rev[i], rev[j] = rev[j], rev[i]
	}
	return rev
}

// bfsCycle returns a shortest cycle through a of length at least one, as the
// loop states [a, q1, ..., q_{L-1}] with an implicit edge from the last state
// back to a. It returns nil if a is not on any cycle.
func bfsCycle(b *Buchi, a int) []int {
	parent := make([]int, b.n)
	for i := range parent {
		parent[i] = -2
	}
	parent[a] = -1 // sentinel so reconstruct stops at the cycle head
	var queue []int
	for _, v := range buchiSucc(b, a) {
		if v == a {
			return []int{a}
		}
		if parent[v] == -2 {
			parent[v] = a
			queue = append(queue, v)
		}
	}
	for len(queue) > 0 {
		u := queue[0]
		queue = queue[1:]
		for _, v := range buchiSucc(b, u) {
			if v == a {
				path := reconstruct(parent, u) // a ... u
				return path
			}
			if parent[v] == -2 {
				parent[v] = u
				queue = append(queue, v)
			}
		}
	}
	return nil
}

// AcceptingLasso searches for an accepting run and returns it as a [Lasso]. The
// second result is false when the automaton's language is empty. The witness is
// a shortest-prefix, shortest-loop run through some accepting state that lies on
// a cycle and is reachable from an initial state.
func (b *Buchi) AcceptingLasso() (Lasso, bool) {
	reach := b.ReachableStates()
	init := b.InitialStates()
	// Deterministic order over accepting, reachable states.
	for a := 0; a < b.n; a++ {
		if !b.accept[a] || !reach.Contains(a) {
			continue
		}
		loop := bfsCycle(b, a)
		if loop == nil {
			continue
		}
		stem := bfsPath(b, init, a)
		if stem == nil {
			continue
		}
		// stem ends at a; drop the last element (a) so the loop begins at a.
		prefix := stem[:len(stem)-1]
		return Lasso{Prefix: prefix, Loop: loop}, true
	}
	return Lasso{}, false
}

// IsEmpty reports whether the automaton accepts no infinite word, i.e. whether
// no accepting state that lies on a cycle is reachable from an initial state.
func (b *Buchi) IsEmpty() bool {
	_, ok := b.AcceptingLasso()
	return !ok
}

// tarjanSCC returns the strongly connected components of the automaton's
// successor graph as a slice of state slices. Each component is returned in the
// order discovered; singleton components without a self-loop represent trivial
// SCCs.
func tarjanSCC(b *Buchi) [][]int {
	index := make([]int, b.n)
	low := make([]int, b.n)
	onStack := make([]bool, b.n)
	for i := range index {
		index[i] = -1
	}
	var stack []int
	counter := 0
	var comps [][]int
	var strongconnect func(v int)
	strongconnect = func(v int) {
		index[v] = counter
		low[v] = counter
		counter++
		stack = append(stack, v)
		onStack[v] = true
		for _, w := range buchiSucc(b, v) {
			if index[w] == -1 {
				strongconnect(w)
				if low[w] < low[v] {
					low[v] = low[w]
				}
			} else if onStack[w] {
				if index[w] < low[v] {
					low[v] = index[w]
				}
			}
		}
		if low[v] == index[v] {
			var comp []int
			for {
				w := stack[len(stack)-1]
				stack = stack[:len(stack)-1]
				onStack[w] = false
				comp = append(comp, w)
				if w == v {
					break
				}
			}
			comps = append(comps, comp)
		}
	}
	for v := 0; v < b.n; v++ {
		if index[v] == -1 {
			strongconnect(v)
		}
	}
	return comps
}

// SCCs returns the strongly connected components of the automaton.
func (b *Buchi) SCCs() [][]int { return tarjanSCC(b) }

// hasSelfLoop reports whether s has a non-contradictory edge to itself.
func (b *Buchi) hasSelfLoop(s int) bool {
	for _, v := range buchiSucc(b, s) {
		if v == s {
			return true
		}
	}
	return false
}

// IsEmptySCC decides emptiness via strongly connected components: the language
// is non-empty iff some non-trivial SCC reachable from an initial state contains
// an accepting state. It is an independent cross-check of [Buchi.IsEmpty].
func (b *Buchi) IsEmptySCC() bool {
	reach := b.ReachableStates()
	for _, comp := range tarjanSCC(b) {
		nontrivial := len(comp) > 1
		if !nontrivial && len(comp) == 1 {
			nontrivial = b.hasSelfLoop(comp[0])
		}
		if !nontrivial {
			continue
		}
		for _, s := range comp {
			if b.accept[s] && reach.Contains(s) {
				return false
			}
		}
	}
	return true
}
