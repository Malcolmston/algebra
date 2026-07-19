package automata

// thompsonBuilder incrementally allocates states and edges while translating a
// regex syntax tree into an NFA fragment.
type thompsonBuilder struct {
	nfa      *NFA
	alphabet []rune
}

// newState allocates a fresh state and returns its index.
func (b *thompsonBuilder) newState() int {
	q := b.nfa.NumStates
	b.nfa.NumStates++
	return q
}

// fragment holds the start and single accepting state of a sub-NFA.
type fragment struct {
	start, accept int
}

// Thompson builds an NFA from a regex syntax tree using Thompson's
// construction. Wildcards ('.') are expanded over alphabet together with any
// literal symbols in the tree; if alphabet is nil only the literals are used.
// The resulting NFA has exactly one accepting state and uses epsilon edges to
// glue fragments together.
func Thompson(node RegexNode, alphabet []rune) *NFA {
	full := mergeRunes(CollectAlphabet(node), alphabet)
	b := &thompsonBuilder{
		nfa:      NewNFA(0, full, 0),
		alphabet: full,
	}
	frag := b.build(node)
	b.nfa.Start = frag.start
	b.nfa.Accept = NewStateSet(frag.accept)
	return b.nfa
}

// build recursively constructs the fragment for node.
func (b *thompsonBuilder) build(node RegexNode) fragment {
	switch t := node.(type) {
	case EmptySet:
		// No path from start to accept.
		s := b.newState()
		a := b.newState()
		return fragment{s, a}
	case EmptyString:
		s := b.newState()
		a := b.newState()
		b.nfa.AddEpsilon(s, a)
		return fragment{s, a}
	case Symbol:
		s := b.newState()
		a := b.newState()
		b.nfa.AddTransition(s, t.R, a)
		return fragment{s, a}
	case AnyChar:
		s := b.newState()
		a := b.newState()
		for _, r := range b.alphabet {
			b.nfa.AddTransition(s, r, a)
		}
		return fragment{s, a}
	case Concat:
		l := b.build(t.L)
		r := b.build(t.R)
		b.nfa.AddEpsilon(l.accept, r.start)
		return fragment{l.start, r.accept}
	case Alternate:
		s := b.newState()
		a := b.newState()
		l := b.build(t.L)
		r := b.build(t.R)
		b.nfa.AddEpsilon(s, l.start)
		b.nfa.AddEpsilon(s, r.start)
		b.nfa.AddEpsilon(l.accept, a)
		b.nfa.AddEpsilon(r.accept, a)
		return fragment{s, a}
	case Star:
		s := b.newState()
		a := b.newState()
		inner := b.build(t.X)
		b.nfa.AddEpsilon(s, inner.start)
		b.nfa.AddEpsilon(s, a)
		b.nfa.AddEpsilon(inner.accept, inner.start)
		b.nfa.AddEpsilon(inner.accept, a)
		return fragment{s, a}
	case Plus:
		s := b.newState()
		a := b.newState()
		inner := b.build(t.X)
		b.nfa.AddEpsilon(s, inner.start)
		b.nfa.AddEpsilon(inner.accept, inner.start)
		b.nfa.AddEpsilon(inner.accept, a)
		return fragment{s, a}
	case Optional:
		s := b.newState()
		a := b.newState()
		inner := b.build(t.X)
		b.nfa.AddEpsilon(s, inner.start)
		b.nfa.AddEpsilon(s, a)
		b.nfa.AddEpsilon(inner.accept, a)
		return fragment{s, a}
	default:
		// Unknown node: treat as empty set.
		s := b.newState()
		a := b.newState()
		return fragment{s, a}
	}
}

// RegexToNFA parses pattern and returns an equivalent NFA via Thompson's
// construction. The NFA's alphabet is the set of literal symbols in the
// pattern; wildcards expand over that same set.
func RegexToNFA(pattern string) (*NFA, error) {
	node, err := ParseRegex(pattern)
	if err != nil {
		return nil, err
	}
	return Thompson(node, nil), nil
}

// RegexToNFAWithAlphabet is like RegexToNFA but additionally includes the given
// alphabet symbols, which is required when the pattern uses '.' and the
// intended alphabet is larger than the literals present.
func RegexToNFAWithAlphabet(pattern string, alphabet []rune) (*NFA, error) {
	node, err := ParseRegex(pattern)
	if err != nil {
		return nil, err
	}
	return Thompson(node, alphabet), nil
}

// RegexToDFA parses pattern, builds an NFA and determinises it to a minimal
// DFA.
func RegexToDFA(pattern string) (*DFA, error) {
	n, err := RegexToNFA(pattern)
	if err != nil {
		return nil, err
	}
	return SubsetConstruction(n).Minimize(), nil
}

// Regexp is a compiled regular expression backed by a DFA, offering fast
// repeated matching.
type Regexp struct {
	pattern string
	dfa     *DFA
}

// Compile parses pattern and returns a compiled Regexp whose alphabet is the
// set of literal symbols in the pattern.
func Compile(pattern string) (*Regexp, error) {
	dfa, err := RegexToDFA(pattern)
	if err != nil {
		return nil, err
	}
	return &Regexp{pattern: pattern, dfa: dfa}, nil
}

// CompileWithAlphabet parses pattern using the supplied alphabet for wildcard
// expansion and returns a compiled Regexp.
func CompileWithAlphabet(pattern string, alphabet []rune) (*Regexp, error) {
	n, err := RegexToNFAWithAlphabet(pattern, alphabet)
	if err != nil {
		return nil, err
	}
	return &Regexp{pattern: pattern, dfa: SubsetConstruction(n).Minimize()}, nil
}

// MustCompile is like Compile but panics on a malformed pattern.
func MustCompile(pattern string) *Regexp {
	re, err := Compile(pattern)
	if err != nil {
		panic(err)
	}
	return re
}

// Matches reports whether the whole input string is matched by the regular
// expression.
func (re *Regexp) Matches(input string) bool {
	return re.dfa.Accepts(input)
}

// Pattern returns the source pattern the Regexp was compiled from.
func (re *Regexp) Pattern() string {
	return re.pattern
}

// DFA returns the underlying deterministic automaton (a clone, so callers may
// mutate it freely).
func (re *Regexp) DFA() *DFA {
	return re.dfa.Clone()
}

// Match is a one-shot helper that reports whether pattern matches the entire
// input string.
func Match(pattern, input string) (bool, error) {
	re, err := Compile(pattern)
	if err != nil {
		return false, err
	}
	return re.Matches(input), nil
}
