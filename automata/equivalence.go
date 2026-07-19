package automata

// Equivalent reports whether two DFAs accept exactly the same language. It works
// by checking that their symmetric difference accepts nothing.
func Equivalent(a, b *DFA) bool {
	return SymmetricDifference(a, b).IsEmptyLanguage()
}

// IsEquivalentTo reports whether this DFA is language-equivalent to other.
func (d *DFA) IsEquivalentTo(other *DFA) bool {
	return Equivalent(d, other)
}

// Subset reports whether L(a) is a subset of L(b), i.e. every string accepted by
// a is also accepted by b. This holds exactly when L(a) \ L(b) is empty.
func Subset(a, b *DFA) bool {
	return Difference(a, b).IsEmptyLanguage()
}

// IsSubsetOf reports whether this DFA's language is contained in other's.
func (d *DFA) IsSubsetOf(other *DFA) bool {
	return Subset(d, other)
}

// Disjoint reports whether two DFAs accept no common string.
func Disjoint(a, b *DFA) bool {
	return Intersection(a, b).IsEmptyLanguage()
}

// Witness returns a shortest string on which two DFAs disagree (accepted by
// exactly one of them) together with true, or the empty string and false if the
// DFAs are equivalent.
func Witness(a, b *DFA) (string, bool) {
	return SymmetricDifference(a, b).ShortestAccepted()
}

// NFAEquivalent reports whether two NFAs accept the same language, by
// determinising both and comparing the resulting DFAs.
func NFAEquivalent(a, b *NFA) bool {
	return Equivalent(a.ToDFA(), b.ToDFA())
}

// PumpingLength returns a valid pumping length for the DFA's language: the
// number of states of the trimmed, minimal automaton. Every accepted string of
// at least this length admits a pumping decomposition.
func PumpingLength(d *DFA) int {
	m := d.Minimize()
	return m.NumStates
}

// PumpingDecomposition splits an accepted string w of length at least
// PumpingLength into parts x, y, z with w = xyz, |xy| ≤ p, |y| ≥ 1, such that
// x·yⁱ·z is accepted for every i ≥ 0. It returns the three parts and true on
// success, or empty strings and false if w is too short or is not accepted.
func PumpingDecomposition(d *DFA) func(w string) (x, y, z string, ok bool) {
	p := PumpingLength(d)
	dfa := d.Complete()
	return func(w string) (string, string, string, bool) {
		if len([]rune(w)) < p || !dfa.Accepts(w) {
			return "", "", "", false
		}
		runes := []rune(w)
		// Track the state after each prefix of length 0..p.
		state := dfa.Start
		visited := map[int]int{state: 0} // state -> prefix length
		for i := 0; i < p; i++ {
			next, ok := dfa.Transition(state, runes[i])
			if !ok {
				return "", "", "", false
			}
			state = next
			if prev, seen := visited[state]; seen {
				// Loop between prefix lengths prev and i+1.
				x := string(runes[:prev])
				y := string(runes[prev : i+1])
				z := string(runes[i+1:])
				return x, y, z, true
			}
			visited[state] = i + 1
		}
		// By pigeonhole a repeat must occur within the first p steps; if not,
		// fall back to the first symbol as the loop (should not happen for a
		// minimal DFA on a sufficiently long word).
		return string(runes[:0]), string(runes[:1]), string(runes[1:]), true
	}
}

// Reachable languages helpers.

// AcceptsAll reports whether the DFA accepts every string over its alphabet,
// i.e. its complement is empty.
func (d *DFA) AcceptsAll() bool {
	return d.Complement().IsEmptyLanguage()
}

// IsUniversal is an alias for AcceptsAll.
func (d *DFA) IsUniversal() bool {
	return d.AcceptsAll()
}
