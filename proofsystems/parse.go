package proofsystems

import (
	"errors"
	"fmt"
)

// ErrParse indicates a syntax error while reading a formula or term from text.
var ErrParse = errors.New("proofsystems: parse error")

// ParseFormula parses a propositional or first-order formula from its textual
// form. The grammar uses ! for negation, & for conjunction, | for disjunction,
// -> for implication, <-> for the biconditional, T and F for the logical
// constants, and the keywords forall/exists for quantifiers written as
// "forall x. body". Implication and the biconditional are right-associative;
// conjunction and disjunction are left-associative; negation binds tightest.
// Identifiers may be predicate symbols applied to comma-separated argument
// terms, as in "P(x,f(a))".
func ParseFormula(src string) (Formula, error) {
	p := &parser{toks: tokenize(src)}
	f, err := p.parseIff()
	if err != nil {
		return Formula{}, err
	}
	if p.pos != len(p.toks) {
		return Formula{}, fmt.Errorf("%w: unexpected trailing input %q", ErrParse, p.peek())
	}
	return f, nil
}

// MustParseFormula parses src and panics on error. It is intended for tests and
// static formula literals.
func MustParseFormula(src string) Formula {
	f, err := ParseFormula(src)
	if err != nil {
		panic(err)
	}
	return f
}

// ParseTerm parses a single first-order term from text, such as "f(x,a)".
// Lowercase-initial bare identifiers are read as variables; identifiers applied
// to arguments are function symbols; other bare identifiers are constants by
// the convention that a single leading uppercase or digit denotes a constant.
func ParseTerm(src string) (Term, error) {
	p := &parser{toks: tokenize(src)}
	t, err := p.parseTerm()
	if err != nil {
		return Term{}, err
	}
	if p.pos != len(p.toks) {
		return Term{}, fmt.Errorf("%w: trailing input %q", ErrParse, p.peek())
	}
	return t, nil
}

// MustParseTerm parses src and panics on error.
func MustParseTerm(src string) Term {
	t, err := ParseTerm(src)
	if err != nil {
		panic(err)
	}
	return t
}

type parser struct {
	toks []string
	pos  int
}

func (p *parser) peek() string {
	if p.pos < len(p.toks) {
		return p.toks[p.pos]
	}
	return ""
}

func (p *parser) next() string {
	t := p.peek()
	p.pos++
	return t
}

func (p *parser) parseIff() (Formula, error) {
	left, err := p.parseImp()
	if err != nil {
		return Formula{}, err
	}
	if p.peek() == "<->" {
		p.next()
		right, err := p.parseIff()
		if err != nil {
			return Formula{}, err
		}
		return Iff(left, right), nil
	}
	return left, nil
}

func (p *parser) parseImp() (Formula, error) {
	left, err := p.parseOr()
	if err != nil {
		return Formula{}, err
	}
	if p.peek() == "->" {
		p.next()
		right, err := p.parseImp()
		if err != nil {
			return Formula{}, err
		}
		return Imp(left, right), nil
	}
	return left, nil
}

func (p *parser) parseOr() (Formula, error) {
	left, err := p.parseAnd()
	if err != nil {
		return Formula{}, err
	}
	for p.peek() == "|" {
		p.next()
		right, err := p.parseAnd()
		if err != nil {
			return Formula{}, err
		}
		left = Or(left, right)
	}
	return left, nil
}

func (p *parser) parseAnd() (Formula, error) {
	left, err := p.parseUnary()
	if err != nil {
		return Formula{}, err
	}
	for p.peek() == "&" {
		p.next()
		right, err := p.parseUnary()
		if err != nil {
			return Formula{}, err
		}
		left = And(left, right)
	}
	return left, nil
}

func (p *parser) parseUnary() (Formula, error) {
	switch p.peek() {
	case "!":
		p.next()
		sub, err := p.parseUnary()
		if err != nil {
			return Formula{}, err
		}
		return Not(sub), nil
	case "forall", "exists":
		kw := p.next()
		v := p.next()
		if !isIdent(v) {
			return Formula{}, fmt.Errorf("%w: expected variable after %s", ErrParse, kw)
		}
		if p.peek() == "." {
			p.next()
		}
		body, err := p.parseUnary()
		if err != nil {
			return Formula{}, err
		}
		if kw == "forall" {
			return Forall(v, body), nil
		}
		return Exists(v, body), nil
	default:
		return p.parseAtom()
	}
}

func (p *parser) parseAtom() (Formula, error) {
	tok := p.peek()
	switch {
	case tok == "(":
		p.next()
		f, err := p.parseIff()
		if err != nil {
			return Formula{}, err
		}
		if p.next() != ")" {
			return Formula{}, fmt.Errorf("%w: missing closing parenthesis", ErrParse)
		}
		return f, nil
	case tok == "T":
		p.next()
		return Top(), nil
	case tok == "F":
		p.next()
		return Bot(), nil
	case isIdent(tok):
		p.next()
		if p.peek() == "(" {
			p.next()
			args, err := p.parseTermArgs()
			if err != nil {
				return Formula{}, err
			}
			return Atom(tok, args...), nil
		}
		return Prop(tok), nil
	default:
		return Formula{}, fmt.Errorf("%w: unexpected token %q", ErrParse, tok)
	}
}

func (p *parser) parseTermArgs() ([]Term, error) {
	var args []Term
	for {
		t, err := p.parseTerm()
		if err != nil {
			return nil, err
		}
		args = append(args, t)
		switch p.peek() {
		case ",":
			p.next()
		case ")":
			p.next()
			return args, nil
		default:
			return nil, fmt.Errorf("%w: expected , or ) in argument list", ErrParse)
		}
	}
}

func (p *parser) parseTerm() (Term, error) {
	tok := p.peek()
	if !isIdent(tok) {
		return Term{}, fmt.Errorf("%w: expected term, got %q", ErrParse, tok)
	}
	p.next()
	if p.peek() == "(" {
		p.next()
		args, err := p.parseTermArgs()
		if err != nil {
			return Term{}, err
		}
		return NewFunc(tok, args...), nil
	}
	if isVarName(tok) {
		return NewVar(tok), nil
	}
	return NewConst(tok), nil
}

func tokenize(src string) []string {
	var toks []string
	i := 0
	for i < len(src) {
		c := src[i]
		switch {
		case c == ' ' || c == '\t' || c == '\n' || c == '\r':
			i++
		case c == '(' || c == ')' || c == '!' || c == '&' || c == '|' || c == ',' || c == '.':
			toks = append(toks, string(c))
			i++
		case c == '-' && i+1 < len(src) && src[i+1] == '>':
			toks = append(toks, "->")
			i += 2
		case c == '<' && i+2 < len(src) && src[i+1] == '-' && src[i+2] == '>':
			toks = append(toks, "<->")
			i += 3
		case isIdentByte(c):
			j := i
			for j < len(src) && isIdentByte(src[j]) {
				j++
			}
			toks = append(toks, src[i:j])
			i = j
		default:
			// Skip unknown byte; the parser will report a structural error.
			i++
		}
	}
	return toks
}

func isIdentByte(c byte) bool {
	return c >= 'a' && c <= 'z' || c >= 'A' && c <= 'Z' || c >= '0' && c <= '9' || c == '_'
}

func isIdent(s string) bool {
	if s == "" {
		return false
	}
	for i := 0; i < len(s); i++ {
		if !isIdentByte(s[i]) {
			return false
		}
	}
	return true
}

// isVarName treats an identifier starting with one of the lowercase letters
// u,v,w,x,y,z as a first-order variable, following the usual mathematical
// convention; every other unapplied identifier is read as a constant symbol.
func isVarName(s string) bool {
	if s == "" {
		return false
	}
	c := s[0]
	return c >= 'u' && c <= 'z'
}
