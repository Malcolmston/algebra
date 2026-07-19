package automata

// Product builds the synchronous product of two DFAs over the combined
// alphabet. A product state is accepting exactly when accept(aAcc, bAcc) is
// true for the accepting status of its two components, allowing union,
// intersection, difference and symmetric difference to be expressed uniformly.
// Both inputs are completed over the shared alphabet before the product is
// formed, and the result is trimmed of unreachable states.
func Product(a, b *DFA, accept func(aAcc, bAcc bool) bool) *DFA {
	alphabet := mergeRunes(a.Alphabet, b.Alphabet)
	ca := a.withAlphabet(alphabet).Complete()
	cb := b.withAlphabet(alphabet).Complete()

	nb := cb.NumStates
	idx := func(qa, qb int) int { return qa*nb + qb }

	out := NewDFA(ca.NumStates*nb, alphabet, idx(ca.Start, cb.Start))
	for qa := 0; qa < ca.NumStates; qa++ {
		for qb := 0; qb < nb; qb++ {
			for _, sym := range alphabet {
				ta, oka := ca.Transition(qa, sym)
				tb, okb := cb.Transition(qb, sym)
				if oka && okb {
					out.SetTransition(idx(qa, qb), sym, idx(ta, tb))
				}
			}
			if accept(ca.Accept[qa], cb.Accept[qb]) {
				out.Accept[idx(qa, qb)] = true
			}
		}
	}
	return out.RemoveUnreachable()
}

// withAlphabet returns a copy of the DFA whose alphabet is extended to include
// every symbol in alphabet (transitions are unchanged; new symbols are simply
// undefined until Complete adds a dead state).
func (d *DFA) withAlphabet(alphabet []rune) *DFA {
	out := d.Clone()
	out.Alphabet = mergeRunes(out.Alphabet, alphabet)
	return out
}

// Union returns a DFA accepting L(a) ∪ L(b).
func Union(a, b *DFA) *DFA {
	return Product(a, b, func(x, y bool) bool { return x || y })
}

// Intersection returns a DFA accepting L(a) ∩ L(b).
func Intersection(a, b *DFA) *DFA {
	return Product(a, b, func(x, y bool) bool { return x && y })
}

// Difference returns a DFA accepting L(a) \ L(b).
func Difference(a, b *DFA) *DFA {
	return Product(a, b, func(x, y bool) bool { return x && !y })
}

// SymmetricDifference returns a DFA accepting the symmetric difference
// (L(a) \ L(b)) ∪ (L(b) \ L(a)), i.e. strings accepted by exactly one of a, b.
func SymmetricDifference(a, b *DFA) *DFA {
	return Product(a, b, func(x, y bool) bool { return x != y })
}

// Union returns a DFA accepting the union of this DFA's language with other's.
func (d *DFA) Union(other *DFA) *DFA {
	return Union(d, other)
}

// Intersection returns a DFA accepting the intersection of the two languages.
func (d *DFA) Intersection(other *DFA) *DFA {
	return Intersection(d, other)
}

// Difference returns a DFA accepting the set difference L(d) \ L(other).
func (d *DFA) Difference(other *DFA) *DFA {
	return Difference(d, other)
}

// SymmetricDifference returns a DFA accepting the symmetric difference of the
// two languages.
func (d *DFA) SymmetricDifference(other *DFA) *DFA {
	return SymmetricDifference(d, other)
}
