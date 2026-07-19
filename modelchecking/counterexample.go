package modelchecking

// ExistsPathTo returns a shortest path (inclusive of endpoints) from any of the
// sources to some state of target, following transitions of k without a length
// bound. It returns nil when target is unreachable from the sources.
func ExistsPathTo(k *Kripke, sources []int, target StateSet) []int {
	return bfsKripkePath(k, sources, target, k.n)
}

// EFWitness returns a finite path from state s to a state satisfying target,
// witnessing EF target at s. The second result is false when no such path
// exists (s does not satisfy EF target).
func EFWitness(k *Kripke, s int, target StateSet) ([]int, bool) {
	p := bfsKripkePath(k, []int{s}, target, k.n)
	return p, p != nil
}

// EGWitness returns an infinite (ultimately periodic) path starting at s all of
// whose states satisfy invariant, witnessing EG invariant at s. The second
// result is false when s does not satisfy EG invariant.
func EGWitness(k *Kripke, s int, invariant StateSet) (Lasso, bool) {
	eg := SatEG(k, invariant)
	if !eg.Contains(s) {
		return Lasso{}, false
	}
	// Walk within eg until a state repeats, forming the loop. Every state of eg
	// has a successor in eg, so the walk cannot get stuck.
	var path []int
	pos := map[int]int{}
	cur := s
	for {
		if p, ok := pos[cur]; ok {
			return Lasso{Prefix: append([]int(nil), path[:p]...), Loop: append([]int(nil), path[p:]...)}, true
		}
		pos[cur] = len(path)
		path = append(path, cur)
		next := -1
		for _, v := range k.succ[cur] {
			if eg.Contains(v) {
				next = v
				break
			}
		}
		if next == -1 {
			// Should not happen for states in EG; guard defensively.
			return Lasso{Prefix: path, Loop: []int{cur}}, true
		}
		cur = next
	}
}

// EUWitness returns a finite path from s along which a holds until b, witnessing
// E[a U b] at s. The second result is false when s does not satisfy E[a U b].
func EUWitness(k *Kripke, s int, a, b StateSet) ([]int, bool) {
	eu := SatEU(k, a, b)
	if !eu.Contains(s) {
		return nil, false
	}
	// Greedily descend towards b staying inside eu; a holds at every non-final
	// state by construction of the fixpoint.
	var path []int
	cur := s
	visited := NewStateSet(k.n)
	for {
		path = append(path, cur)
		if b.Contains(cur) {
			return path, true
		}
		visited.Add(cur)
		next := -1
		for _, v := range k.succ[cur] {
			if eu.Contains(v) && !visited.Contains(v) {
				next = v
				break
			}
		}
		if next == -1 {
			return path, true
		}
		cur = next
	}
}

// AGCounterexample returns a finite path from an initial state to a state
// violating p, witnessing the failure of AG p. The second result is false when
// AG p holds at every initial state.
func AGCounterexample(k *Kripke, p StateSet) (*Counterexample, bool) {
	sat := SatAG(k, p)
	var bad []int
	for _, s := range k.InitialStates() {
		if !sat.Contains(s) {
			bad = append(bad, s)
		}
	}
	if len(bad) == 0 {
		return nil, false
	}
	path := bfsKripkePath(k, bad, p.Complement(), k.n)
	if path == nil {
		return nil, false
	}
	return &Counterexample{Kripke: k, Lasso: Lasso{Prefix: path}, Formula: "AG(p)"}, true
}

// AFCounterexample returns an infinite path from an initial state that never
// reaches p, witnessing the failure of AF p (equivalently a witness to EG ¬p).
// The second result is false when AF p holds at every initial state.
func AFCounterexample(k *Kripke, p StateSet) (*Counterexample, bool) {
	notP := p.Complement()
	egNotP := SatEG(k, notP)
	var start int = -1
	for _, s := range k.InitialStates() {
		if egNotP.Contains(s) {
			start = s
			break
		}
	}
	if start == -1 {
		return nil, false
	}
	lasso, ok := EGWitness(k, start, notP)
	if !ok {
		return nil, false
	}
	return &Counterexample{Kripke: k, Lasso: lasso, Formula: "AF(p)"}, true
}

// EXWitness returns a single successor of s that satisfies target, witnessing
// EX target. The second result is false when s has no such successor.
func EXWitness(k *Kripke, s int, target StateSet) (int, bool) {
	for _, v := range k.succ[s] {
		if target.Contains(v) {
			return v, true
		}
	}
	return -1, false
}
