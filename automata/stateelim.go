package automata

// unionR returns a regex node for L(a) ∪ L(b), simplifying away ∅ and nil
// operands (nil denotes the empty language).
func unionR(a, b RegexNode) RegexNode {
	if a == nil {
		return b
	}
	if b == nil {
		return a
	}
	if _, ok := a.(EmptySet); ok {
		return b
	}
	if _, ok := b.(EmptySet); ok {
		return a
	}
	return Alternate{L: a, R: b}
}

// concatR returns a regex node for L(a)·L(b), simplifying with ∅ and ε.
func concatR(a, b RegexNode) RegexNode {
	if a == nil || b == nil {
		return nil
	}
	if _, ok := a.(EmptySet); ok {
		return EmptySet{}
	}
	if _, ok := b.(EmptySet); ok {
		return EmptySet{}
	}
	if _, ok := a.(EmptyString); ok {
		return b
	}
	if _, ok := b.(EmptyString); ok {
		return a
	}
	return Concat{L: a, R: b}
}

// starR returns a regex node for L(a)*, simplifying the empty language and empty
// string to ε.
func starR(a RegexNode) RegexNode {
	if a == nil {
		return EmptyString{}
	}
	switch a.(type) {
	case EmptySet:
		return EmptyString{}
	case EmptyString:
		return EmptyString{}
	}
	return Star{X: a}
}

// DFAToRegex converts a DFA into an equivalent regular expression using the
// state-elimination method. Fresh start and final states are added, then the
// original states are eliminated one by one, accumulating regex labels on the
// remaining edges. The returned node denotes exactly the DFA's language; it is
// EmptySet when the language is empty.
func DFAToRegex(d *DFA) RegexNode {
	n := d.NumStates
	newStart := n
	newFinal := n + 1
	total := n + 2

	// R[i][j] is the current edge label from i to j (nil == empty language).
	R := make([][]RegexNode, total)
	for i := range R {
		R[i] = make([]RegexNode, total)
	}

	for q, m := range d.trans {
		for a, to := range m {
			R[q][to] = unionR(R[q][to], Symbol{R: a})
		}
	}
	R[newStart][d.Start] = unionR(R[newStart][d.Start], EmptyString{})
	for q := range d.Accept {
		R[q][newFinal] = unionR(R[q][newFinal], EmptyString{})
	}

	// Eliminate the original states 0..n-1.
	for k := 0; k < n; k++ {
		loop := starR(R[k][k])
		for i := 0; i < total; i++ {
			if i == k || R[i][k] == nil {
				continue
			}
			for j := 0; j < total; j++ {
				if j == k || R[k][j] == nil {
					continue
				}
				add := concatR(concatR(R[i][k], loop), R[k][j])
				R[i][j] = unionR(R[i][j], add)
			}
		}
		// Remove state k's incident edges.
		for x := 0; x < total; x++ {
			R[x][k] = nil
			R[k][x] = nil
		}
	}

	result := R[newStart][newFinal]
	if result == nil {
		return EmptySet{}
	}
	return result
}

// DFAToRegexString returns the state-elimination regular expression for d
// rendered as a string.
func DFAToRegexString(d *DFA) string {
	return DFAToRegex(d).String()
}

// NFAToRegex converts an NFA to an equivalent regular expression by first
// determinising it and then applying state elimination.
func NFAToRegex(n *NFA) RegexNode {
	return DFAToRegex(n.ToDFA())
}
