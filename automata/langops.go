package automata

// copyInto copies every edge and accepting mark of src into dst, shifting all
// state indices by offset. It does not touch dst's start state.
func copyInto(dst, src *NFA, offset int) {
	for q, m := range src.trans {
		for a, tos := range m {
			for _, to := range tos {
				dst.AddTransition(q+offset, a, to+offset)
			}
		}
	}
	for q, tos := range src.eps {
		for _, to := range tos {
			dst.AddEpsilon(q+offset, to+offset)
		}
	}
	for q := range src.Accept {
		dst.Accept[q+offset] = true
	}
}

// EmptyLanguageNFA returns an NFA over alphabet accepting no strings.
func EmptyLanguageNFA(alphabet []rune) *NFA {
	n := NewNFA(1, alphabet, 0)
	return n
}

// EpsilonNFA returns an NFA over alphabet accepting only the empty string.
func EpsilonNFA(alphabet []rune) *NFA {
	n := NewNFA(1, alphabet, 0)
	n.Accept[0] = true
	return n
}

// SymbolNFA returns an NFA accepting exactly the single-character string sym.
func SymbolNFA(sym rune) *NFA {
	n := NewNFA(2, []rune{sym}, 0)
	n.AddTransition(0, sym, 1)
	n.Accept[1] = true
	return n
}

// LiteralNFA returns an NFA accepting exactly the string s.
func LiteralNFA(s string) *NFA {
	rs := []rune(s)
	n := NewNFA(len(rs)+1, rs, 0)
	for i, r := range rs {
		n.AddTransition(i, r, i+1)
	}
	n.Accept[len(rs)] = true
	return n
}

// Concatenate returns an NFA accepting the concatenation L(a)·L(b): every
// string formed by an a-word followed by a b-word.
func Concatenate(a, b *NFA) *NFA {
	na := a.NumStates
	out := NewNFA(na+b.NumStates, mergeRunes(a.Alphabet, b.Alphabet), a.Start)
	copyInto(out, a, 0)
	copyInto(out, b, na)
	// a's accepting states are no longer accepting; connect them to b's start.
	out.Accept = NewStateSet()
	for q := range a.Accept {
		out.AddEpsilon(q, b.Start+na)
	}
	for q := range b.Accept {
		out.Accept[q+na] = true
	}
	return out
}

// UnionNFA returns an NFA accepting L(a) ∪ L(b) via a fresh start state with
// epsilon edges into both operands.
func UnionNFA(a, b *NFA) *NFA {
	na := a.NumStates
	nb := b.NumStates
	newStart := na + nb
	out := NewNFA(na+nb+1, mergeRunes(a.Alphabet, b.Alphabet), newStart)
	out.Accept = NewStateSet()
	copyInto(out, a, 0)
	copyInto(out, b, na)
	out.AddEpsilon(newStart, a.Start)
	out.AddEpsilon(newStart, b.Start+na)
	return out
}

// StarNFA returns an NFA accepting the Kleene star L(a)*, including the empty
// string.
func StarNFA(a *NFA) *NFA {
	na := a.NumStates
	newStart := na
	out := NewNFA(na+1, a.Alphabet, newStart)
	out.Accept = NewStateSet()
	copyInto(out, a, 0)
	out.AddEpsilon(newStart, a.Start)
	out.Accept[newStart] = true
	for q := range a.Accept {
		out.AddEpsilon(q, a.Start)
		out.Accept[q] = true
	}
	return out
}

// PlusNFA returns an NFA accepting L(a)+ = L(a)·L(a)*, i.e. one or more
// repetitions.
func PlusNFA(a *NFA) *NFA {
	out := a.Clone()
	for q := range a.Accept {
		out.AddEpsilon(q, a.Start)
	}
	return out
}

// OptionalNFA returns an NFA accepting L(a) ∪ {ε}.
func OptionalNFA(a *NFA) *NFA {
	na := a.NumStates
	newStart := na
	out := NewNFA(na+1, a.Alphabet, newStart)
	out.Accept = NewStateSet()
	copyInto(out, a, 0)
	out.AddEpsilon(newStart, a.Start)
	out.Accept[newStart] = true
	return out
}

// PowerNFA returns an NFA accepting exactly n consecutive copies of L(a). For
// n = 0 it accepts only the empty string.
func PowerNFA(a *NFA, n int) *NFA {
	if n <= 0 {
		return EpsilonNFA(a.Alphabet)
	}
	out := a.Clone()
	for i := 1; i < n; i++ {
		out = Concatenate(out, a)
	}
	return out
}

// ReverseNFA returns an NFA accepting the reversal of L(a).
func ReverseNFA(a *NFA) *NFA {
	return a.Reverse()
}

// Concat is the method form of Concatenate.
func (n *NFA) Concat(other *NFA) *NFA { return Concatenate(n, other) }

// Union is the method form of UnionNFA.
func (n *NFA) Union(other *NFA) *NFA { return UnionNFA(n, other) }

// Star is the method form of StarNFA.
func (n *NFA) Star() *NFA { return StarNFA(n) }

// Plus is the method form of PlusNFA.
func (n *NFA) Plus() *NFA { return PlusNFA(n) }

// Optional is the method form of OptionalNFA.
func (n *NFA) Optional() *NFA { return OptionalNFA(n) }

// Power is the method form of PowerNFA.
func (n *NFA) Power(k int) *NFA { return PowerNFA(n, k) }

// ReverseDFA returns a DFA accepting the reversal of L(d), obtained by reversing
// the automaton and determinising.
func ReverseDFA(d *DFA) *DFA {
	return d.ToNFA().Reverse().ToDFA()
}

// Reverse returns a DFA accepting the reversal of this DFA's language.
func (d *DFA) Reverse() *DFA {
	return ReverseDFA(d)
}

// Concatenate returns a DFA accepting L(d)·L(other).
func (d *DFA) Concatenate(other *DFA) *DFA {
	return Concatenate(d.ToNFA(), other.ToNFA()).ToDFA()
}

// Star returns a DFA accepting L(d)*.
func (d *DFA) Star() *DFA {
	return StarNFA(d.ToNFA()).ToDFA()
}

// Plus returns a DFA accepting L(d)+.
func (d *DFA) Plus() *DFA {
	return PlusNFA(d.ToNFA()).ToDFA()
}
