package automata

import (
	"errors"
	"fmt"
	"strings"
)

// RegexNode is a node in a parsed regular-expression syntax tree. The concrete
// implementations are EmptySet, EmptyString, Symbol, Concat, Alternate, Star,
// Plus, Optional and AnyChar.
type RegexNode interface {
	// String renders the node back to regular-expression syntax.
	String() string
	// isRegexNode is an unexported marker method.
	isRegexNode()
}

// EmptySet denotes the regular expression matching no strings at all (∅).
type EmptySet struct{}

// EmptyString denotes the regular expression matching only the empty string
// (ε).
type EmptyString struct{}

// Symbol denotes a regular expression matching the single character R.
type Symbol struct{ R rune }

// AnyChar denotes the wildcard "." matching any single alphabet symbol.
type AnyChar struct{}

// Concat denotes the concatenation of two sub-expressions L then R.
type Concat struct{ L, R RegexNode }

// Alternate denotes the alternation (union) of two sub-expressions L or R.
type Alternate struct{ L, R RegexNode }

// Star denotes the Kleene star of a sub-expression (zero or more repetitions).
type Star struct{ X RegexNode }

// Plus denotes one-or-more repetitions of a sub-expression.
type Plus struct{ X RegexNode }

// Optional denotes zero-or-one occurrence of a sub-expression.
type Optional struct{ X RegexNode }

func (EmptySet) isRegexNode()    {}
func (EmptyString) isRegexNode() {}
func (Symbol) isRegexNode()      {}
func (AnyChar) isRegexNode()     {}
func (Concat) isRegexNode()      {}
func (Alternate) isRegexNode()   {}
func (Star) isRegexNode()        {}
func (Plus) isRegexNode()        {}
func (Optional) isRegexNode()    {}

// String renders the empty set as the symbol ∅.
func (EmptySet) String() string { return "∅" }

// String renders the empty string as the symbol ε.
func (EmptyString) String() string { return "ε" }

// String renders a single symbol, escaping regex metacharacters.
func (s Symbol) String() string {
	if strings.ContainsRune("()|*+?.\\", s.R) {
		return "\\" + string(s.R)
	}
	return string(s.R)
}

// String renders the wildcard as ".".
func (AnyChar) String() string { return "." }

// String renders a concatenation.
func (c Concat) String() string { return c.L.String() + c.R.String() }

// String renders an alternation with parentheses.
func (a Alternate) String() string { return "(" + a.L.String() + "|" + a.R.String() + ")" }

// String renders a Kleene star, parenthesising compound operands.
func (s Star) String() string { return starWrap(s.X) + "*" }

// String renders a plus, parenthesising compound operands.
func (p Plus) String() string { return starWrap(p.X) + "+" }

// String renders an optional, parenthesising compound operands.
func (o Optional) String() string { return starWrap(o.X) + "?" }

func starWrap(x RegexNode) string {
	switch x.(type) {
	case Symbol, AnyChar, EmptyString, EmptySet:
		return x.String()
	default:
		return "(" + x.String() + ")"
	}
}

// Constructor helpers for building syntax trees programmatically.

// Lit builds a Symbol node matching the single rune r.
func Lit(r rune) RegexNode { return Symbol{R: r} }

// Eps builds an EmptyString (ε) node.
func Eps() RegexNode { return EmptyString{} }

// Empty builds an EmptySet (∅) node.
func Empty() RegexNode { return EmptySet{} }

// Dot builds an AnyChar (.) wildcard node.
func Dot() RegexNode { return AnyChar{} }

// Seq builds the left-associative concatenation of the given nodes. With no
// arguments it returns ε; with one it returns that node unchanged.
func Seq(nodes ...RegexNode) RegexNode {
	if len(nodes) == 0 {
		return EmptyString{}
	}
	acc := nodes[0]
	for _, n := range nodes[1:] {
		acc = Concat{L: acc, R: n}
	}
	return acc
}

// Alt builds the left-associative alternation of the given nodes. With no
// arguments it returns ∅; with one it returns that node unchanged.
func Alt(nodes ...RegexNode) RegexNode {
	if len(nodes) == 0 {
		return EmptySet{}
	}
	acc := nodes[0]
	for _, n := range nodes[1:] {
		acc = Alternate{L: acc, R: n}
	}
	return acc
}

// Kleene builds a Star node over x.
func Kleene(x RegexNode) RegexNode { return Star{X: x} }

// OneOrMore builds a Plus node over x.
func OneOrMore(x RegexNode) RegexNode { return Plus{X: x} }

// ZeroOrOne builds an Optional node over x.
func ZeroOrOne(x RegexNode) RegexNode { return Optional{X: x} }

// LiteralNode builds a syntax tree matching exactly the string s (a
// concatenation of its runes; ε when s is empty).
func LiteralNode(s string) RegexNode {
	rs := []rune(s)
	nodes := make([]RegexNode, len(rs))
	for i, r := range rs {
		nodes[i] = Symbol{R: r}
	}
	return Seq(nodes...)
}

// regexParser is a recursive-descent parser for the supported grammar:
//
//	regex   := alt
//	alt     := concat ('|' concat)*
//	concat  := repeat*
//	repeat  := atom ('*' | '+' | '?')*
//	atom    := '(' regex ')' | '.' | '\' char | char
type regexParser struct {
	src []rune
	pos int
}

// ParseRegex parses a regular expression written in the supported syntax and
// returns its syntax tree. Supported constructs are concatenation, alternation
// with '|', Kleene star '*', plus '+', optional '?', grouping with parentheses,
// the wildcard '.', and backslash escaping of metacharacters. An empty pattern
// parses to ε.
func ParseRegex(pattern string) (RegexNode, error) {
	p := &regexParser{src: []rune(pattern)}
	node, err := p.parseAlt()
	if err != nil {
		return nil, err
	}
	if p.pos != len(p.src) {
		return nil, fmt.Errorf("automata: unexpected character %q at position %d", string(p.src[p.pos]), p.pos)
	}
	return node, nil
}

// MustParseRegex is like ParseRegex but panics if the pattern is invalid. It is
// intended for package-level expressions and tests.
func MustParseRegex(pattern string) RegexNode {
	node, err := ParseRegex(pattern)
	if err != nil {
		panic(err)
	}
	return node
}

func (p *regexParser) peek() (rune, bool) {
	if p.pos < len(p.src) {
		return p.src[p.pos], true
	}
	return 0, false
}

func (p *regexParser) parseAlt() (RegexNode, error) {
	left, err := p.parseConcat()
	if err != nil {
		return nil, err
	}
	for {
		r, ok := p.peek()
		if !ok || r != '|' {
			break
		}
		p.pos++ // consume '|'
		right, err := p.parseConcat()
		if err != nil {
			return nil, err
		}
		left = Alternate{L: left, R: right}
	}
	return left, nil
}

func (p *regexParser) parseConcat() (RegexNode, error) {
	var factors []RegexNode
	for {
		r, ok := p.peek()
		if !ok || r == '|' || r == ')' {
			break
		}
		f, err := p.parseRepeat()
		if err != nil {
			return nil, err
		}
		factors = append(factors, f)
	}
	if len(factors) == 0 {
		return EmptyString{}, nil
	}
	return Seq(factors...), nil
}

func (p *regexParser) parseRepeat() (RegexNode, error) {
	atom, err := p.parseAtom()
	if err != nil {
		return nil, err
	}
	for {
		r, ok := p.peek()
		if !ok {
			break
		}
		switch r {
		case '*':
			p.pos++
			atom = Star{X: atom}
		case '+':
			p.pos++
			atom = Plus{X: atom}
		case '?':
			p.pos++
			atom = Optional{X: atom}
		default:
			return atom, nil
		}
	}
	return atom, nil
}

func (p *regexParser) parseAtom() (RegexNode, error) {
	r, ok := p.peek()
	if !ok {
		return nil, errors.New("automata: unexpected end of pattern")
	}
	switch r {
	case '(':
		p.pos++
		inner, err := p.parseAlt()
		if err != nil {
			return nil, err
		}
		if cr, ok := p.peek(); !ok || cr != ')' {
			return nil, errors.New("automata: missing closing parenthesis")
		}
		p.pos++ // consume ')'
		return inner, nil
	case '.':
		p.pos++
		return AnyChar{}, nil
	case '\\':
		p.pos++
		esc, ok := p.peek()
		if !ok {
			return nil, errors.New("automata: dangling escape at end of pattern")
		}
		p.pos++
		return Symbol{R: esc}, nil
	case '*', '+', '?':
		return nil, fmt.Errorf("automata: unexpected quantifier %q at position %d", string(r), p.pos)
	case ')':
		return nil, errors.New("automata: unexpected closing parenthesis")
	default:
		p.pos++
		return Symbol{R: r}, nil
	}
}

// CollectAlphabet returns the sorted set of literal symbols appearing in the
// syntax tree (wildcards contribute nothing).
func CollectAlphabet(node RegexNode) []rune {
	seen := map[rune]bool{}
	var walk func(RegexNode)
	walk = func(n RegexNode) {
		switch t := n.(type) {
		case Symbol:
			seen[t.R] = true
		case Concat:
			walk(t.L)
			walk(t.R)
		case Alternate:
			walk(t.L)
			walk(t.R)
		case Star:
			walk(t.X)
		case Plus:
			walk(t.X)
		case Optional:
			walk(t.X)
		}
	}
	walk(node)
	out := make([]rune, 0, len(seen))
	for r := range seen {
		out = append(out, r)
	}
	return sortedRunes(out)
}
