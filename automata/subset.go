package automata

// SubsetConstruction converts an NFA (possibly with epsilon transitions) into
// an equivalent DFA using the classic powerset construction. Each DFA state
// corresponds to an epsilon-closed set of NFA states reachable from the start.
// Only subsets actually reachable during the construction are materialised.
func SubsetConstruction(n *NFA) *DFA {
	alphabet := n.Alphabet
	start := n.StartSet()

	index := map[string]int{}
	var sets []StateSet

	getIndex := func(s StateSet) int {
		key := s.Key()
		if i, ok := index[key]; ok {
			return i
		}
		i := len(sets)
		index[key] = i
		sets = append(sets, s)
		return i
	}

	startIdx := getIndex(start)
	dfa := NewDFA(1, alphabet, startIdx)

	for i := 0; i < len(sets); i++ {
		cur := sets[i]
		for _, a := range alphabet {
			next := n.Step(cur, a)
			if next.IsEmpty() {
				continue
			}
			before := len(sets)
			j := getIndex(next)
			if j >= before {
				dfa.NumStates = len(sets)
			}
			dfa.SetTransition(i, a, j)
		}
	}
	dfa.NumStates = len(sets)

	// Accepting DFA states are those whose subset contains an NFA accepting
	// state.
	for i, s := range sets {
		for q := range s {
			if n.Accept[q] {
				dfa.Accept[i] = true
				break
			}
		}
	}
	return dfa
}

// ToDFA converts the NFA to an equivalent DFA via subset construction.
func (n *NFA) ToDFA() *DFA {
	return SubsetConstruction(n)
}

// Determinize is an alias for ToDFA that additionally removes unreachable
// states (the powerset construction already avoids them) and is provided for
// readability at call sites.
func (n *NFA) Determinize() *DFA {
	return SubsetConstruction(n).RemoveUnreachable()
}

// NFAToDFA is a package-level convenience wrapper around SubsetConstruction.
func NFAToDFA(n *NFA) *DFA {
	return SubsetConstruction(n)
}
