package modelchecking

// Fairness is a set of (unconditional / justice) fairness constraints for a
// Kripke structure. Each constraint is a set of states; a path is fair when it
// visits every constraint set infinitely often. With no constraints every
// infinite path is fair and the fair operators coincide with the ordinary ones.
type Fairness struct {
	Constraints []StateSet
}

// NewFairness returns a fairness condition consisting of the given constraint
// sets, each copied.
func NewFairness(sets ...StateSet) Fairness {
	cs := make([]StateSet, len(sets))
	for i, s := range sets {
		cs[i] = s.Clone()
	}
	return Fairness{Constraints: cs}
}

// NumConstraints returns the number of fairness constraints.
func (f Fairness) NumConstraints() int { return len(f.Constraints) }

// Add appends a constraint set and returns the updated condition.
func (f Fairness) Add(set StateSet) Fairness {
	f.Constraints = append(f.Constraints, set.Clone())
	return f
}

// sccKripkeInduced computes the strongly connected components of the subgraph of
// k induced by the state set within, considering only edges whose endpoints are
// both in within.
func sccKripkeInduced(k *Kripke, within StateSet) [][]int {
	index := make([]int, k.n)
	low := make([]int, k.n)
	onStack := make([]bool, k.n)
	for i := range index {
		index[i] = -1
	}
	var stack []int
	counter := 0
	var comps [][]int
	var sc func(v int)
	sc = func(v int) {
		index[v] = counter
		low[v] = counter
		counter++
		stack = append(stack, v)
		onStack[v] = true
		for _, w := range k.succ[v] {
			if !within.Contains(w) {
				continue
			}
			if index[w] == -1 {
				sc(w)
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
	for _, v := range within.Elements() {
		if index[v] == -1 {
			sc(v)
		}
	}
	return comps
}

// kripkeInternalEdge reports whether some edge exists from u into the set within
// (used to decide whether a singleton SCC is non-trivial via a self-loop).
func kripkeSelfLoopIn(k *Kripke, u int) bool {
	for _, w := range k.succ[u] {
		if w == u {
			return true
		}
	}
	return false
}

// fairSCCStates returns the union of the states of every non-trivial strongly
// connected component of the subgraph induced by within that meets each
// fairness constraint. From any such component an infinite fair path exists.
func fairSCCStates(k *Kripke, within StateSet, fc Fairness) StateSet {
	out := NewStateSet(k.n)
	for _, comp := range sccKripkeInduced(k, within) {
		nontrivial := len(comp) > 1
		if !nontrivial && len(comp) == 1 {
			nontrivial = kripkeSelfLoopIn(k, comp[0])
		}
		if !nontrivial {
			continue
		}
		compSet := StateSetFromSlice(k.n, comp)
		fair := true
		for _, c := range fc.Constraints {
			if !compSet.Intersects(c) {
				fair = false
				break
			}
		}
		if fair {
			out = out.Union(compSet)
		}
	}
	return out
}

// FairEG returns the set of states from which there is a fair path all of whose
// states lie in a. It is computed by locating the non-trivial fair strongly
// connected components inside a and taking the states of a that can reach one
// while remaining in a.
func FairEG(k *Kripke, a StateSet, fc Fairness) StateSet {
	good := fairSCCStates(k, a, fc)
	return SatEU(k, a, good)
}

// FairStates returns the set of states of k from which at least one fair path
// begins. It equals FairEG over the full state set.
func FairStates(k *Kripke, fc Fairness) StateSet {
	return FairEG(k, FullStateSet(k.n), fc)
}

// FairEX returns the set of states satisfying the fair existential next operator:
// some successor satisfies phi and has a fair continuation.
func FairEX(k *Kripke, phi StateSet, fc Fairness) StateSet {
	return PreExists(k, phi.Intersect(FairStates(k, fc)))
}

// FairEU returns the set of states satisfying the fair existential until: a fair
// path along which a holds until b holds.
func FairEU(k *Kripke, a, b StateSet, fc Fairness) StateSet {
	return SatEU(k, a, b.Intersect(FairStates(k, fc)))
}

// FairEF returns the set of states satisfying the fair EF operator: a fair path
// that eventually reaches a.
func FairEF(k *Kripke, a StateSet, fc Fairness) StateSet {
	return FairEU(k, FullStateSet(k.n), a, fc)
}

// FairCTLCheck returns the set of states of k satisfying the CTL formula f under
// the fairness condition fc. The formula is first reduced to the existential
// fragment {¬, ∧, EX, EU, EG} by [CTL.ExistentialNormalForm], and the temporal
// operators are then evaluated with their fair variants.
func FairCTLCheck(k *Kripke, f *CTL, fc Fairness) (StateSet, error) {
	enf := f.ExistentialNormalForm()
	fairSet := FairStates(k, fc)
	return fairEval(k, enf, fc, fairSet)
}

func fairEval(k *Kripke, f *CTL, fc Fairness, fairSet StateSet) (StateSet, error) {
	switch f.Kind {
	case CTLTrueKind:
		return FullStateSet(k.n), nil
	case CTLFalseKind:
		return NewStateSet(k.n), nil
	case CTLAtomKind:
		return k.LabelStateSet(f.Atom), nil
	case CTLNotKind:
		a, err := fairEval(k, f.L, fc, fairSet)
		if err != nil {
			return a, err
		}
		return a.Complement(), nil
	case CTLAndKind:
		a, err := fairEval(k, f.L, fc, fairSet)
		if err != nil {
			return a, err
		}
		b, err := fairEval(k, f.R, fc, fairSet)
		if err != nil {
			return b, err
		}
		return a.Intersect(b), nil
	case CTLEXKind:
		a, err := fairEval(k, f.L, fc, fairSet)
		if err != nil {
			return a, err
		}
		return PreExists(k, a.Intersect(fairSet)), nil
	case CTLEUKind:
		a, err := fairEval(k, f.L, fc, fairSet)
		if err != nil {
			return a, err
		}
		b, err := fairEval(k, f.R, fc, fairSet)
		if err != nil {
			return b, err
		}
		return SatEU(k, a, b.Intersect(fairSet)), nil
	case CTLEGKind:
		a, err := fairEval(k, f.L, fc, fairSet)
		if err != nil {
			return a, err
		}
		return FairEG(k, a, fc), nil
	}
	// Any remaining operator should have been eliminated by the normal form;
	// fall back to the ordinary checker for robustness.
	return CTLCheck(k, f)
}

// FairCTLModelCheck reports whether every initial state of k satisfies f under
// the fairness condition fc.
func FairCTLModelCheck(k *Kripke, f *CTL, fc Fairness) (bool, error) {
	set, err := FairCTLCheck(k, f, fc)
	if err != nil {
		return false, err
	}
	for _, s := range k.InitialStates() {
		if !set.Contains(s) {
			return false, nil
		}
	}
	return true, nil
}
