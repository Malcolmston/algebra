package satsolver

import (
	"fmt"
	"strings"
)

// tokenKind enumerates the lexical categories of the expression parser.
type tokenKind int

const (
	tokEOF tokenKind = iota
	tokIdent
	tokConst
	tokNot
	tokAnd
	tokOr
	tokXor
	tokImplies
	tokIff
	tokLParen
	tokRParen
)

type token struct {
	kind tokenKind
	text string
	val  bool
}

func lex(s string) ([]token, error) {
	var toks []token
	i := 0
	for i < len(s) {
		c := s[i]
		switch {
		case c == ' ' || c == '\t' || c == '\n' || c == '\r':
			i++
		case c == '(':
			toks = append(toks, token{kind: tokLParen})
			i++
		case c == ')':
			toks = append(toks, token{kind: tokRParen})
			i++
		case c == '~' || c == '!':
			// "!&" and "!|" are treated as their bare form is not accepted;
			// a lone '!' or '~' is negation.
			toks = append(toks, token{kind: tokNot})
			i++
		case c == '&':
			toks = append(toks, token{kind: tokAnd})
			i++
		case c == '|':
			toks = append(toks, token{kind: tokOr})
			i++
		case c == '^':
			toks = append(toks, token{kind: tokXor})
			i++
		case c == '-' && i+1 < len(s) && s[i+1] == '>':
			toks = append(toks, token{kind: tokImplies})
			i += 2
		case c == '<' && strings.HasPrefix(s[i:], "<->"):
			toks = append(toks, token{kind: tokIff})
			i += 3
		case c == '=' && i+1 < len(s) && s[i+1] == '=':
			toks = append(toks, token{kind: tokIff})
			i += 2
		case isIdentStart(c):
			j := i + 1
			for j < len(s) && isIdentPart(s[j]) {
				j++
			}
			word := s[i:j]
			switch strings.ToLower(word) {
			case "and":
				toks = append(toks, token{kind: tokAnd})
			case "or":
				toks = append(toks, token{kind: tokOr})
			case "xor":
				toks = append(toks, token{kind: tokXor})
			case "not":
				toks = append(toks, token{kind: tokNot})
			case "iff":
				toks = append(toks, token{kind: tokIff})
			case "true", "t", "1":
				toks = append(toks, token{kind: tokConst, val: true})
			case "false", "f", "0":
				toks = append(toks, token{kind: tokConst, val: false})
			default:
				toks = append(toks, token{kind: tokIdent, text: word})
			}
			i = j
		default:
			return nil, fmt.Errorf("satsolver: unexpected character %q at position %d", c, i)
		}
	}
	toks = append(toks, token{kind: tokEOF})
	return toks, nil
}

func isIdentStart(c byte) bool {
	return c == '_' || (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z')
}

func isIdentPart(c byte) bool {
	return isIdentStart(c) || (c >= '0' && c <= '9')
}

type parser struct {
	toks []token
	pos  int
}

func (p *parser) peek() token         { return p.toks[p.pos] }
func (p *parser) next() token         { t := p.toks[p.pos]; p.pos++; return t }
func (p *parser) at(k tokenKind) bool { return p.toks[p.pos].kind == k }

// Parse parses a Boolean expression from its textual form and returns the
// corresponding [Expr] tree. Supported operators, from lowest to highest
// precedence, are the biconditional ("<->" or "=="), the conditional ("->"),
// disjunction ("|" or "or"), exclusive-or ("^" or "xor"), conjunction ("&" or
// "and") and negation ("~", "!" or "not"). Identifiers are variables and the
// words true/false (also T/F, 1/0) are constants.
func Parse(s string) (Expr, error) {
	toks, err := lex(s)
	if err != nil {
		return nil, err
	}
	p := &parser{toks: toks}
	e, err := p.parseIff()
	if err != nil {
		return nil, err
	}
	if !p.at(tokEOF) {
		return nil, fmt.Errorf("satsolver: unexpected trailing input near token %d", p.pos)
	}
	return e, nil
}

// MustParse is like [Parse] but panics on error. It is intended for use with
// constant expressions in tests and package-level initialisers.
func MustParse(s string) Expr {
	e, err := Parse(s)
	if err != nil {
		panic(err)
	}
	return e
}

func (p *parser) parseIff() (Expr, error) {
	left, err := p.parseImplies()
	if err != nil {
		return nil, err
	}
	for p.at(tokIff) {
		p.next()
		right, err := p.parseImplies()
		if err != nil {
			return nil, err
		}
		left = Iff{X: left, Y: right}
	}
	return left, nil
}

func (p *parser) parseImplies() (Expr, error) {
	left, err := p.parseOr()
	if err != nil {
		return nil, err
	}
	if p.at(tokImplies) {
		p.next()
		// right associative
		right, err := p.parseImplies()
		if err != nil {
			return nil, err
		}
		return Implies{X: left, Y: right}, nil
	}
	return left, nil
}

func (p *parser) parseOr() (Expr, error) {
	left, err := p.parseXor()
	if err != nil {
		return nil, err
	}
	for p.at(tokOr) {
		p.next()
		right, err := p.parseXor()
		if err != nil {
			return nil, err
		}
		left = Or{X: left, Y: right}
	}
	return left, nil
}

func (p *parser) parseXor() (Expr, error) {
	left, err := p.parseAnd()
	if err != nil {
		return nil, err
	}
	for p.at(tokXor) {
		p.next()
		right, err := p.parseAnd()
		if err != nil {
			return nil, err
		}
		left = Xor{X: left, Y: right}
	}
	return left, nil
}

func (p *parser) parseAnd() (Expr, error) {
	left, err := p.parseUnary()
	if err != nil {
		return nil, err
	}
	for p.at(tokAnd) {
		p.next()
		right, err := p.parseUnary()
		if err != nil {
			return nil, err
		}
		left = And{X: left, Y: right}
	}
	return left, nil
}

func (p *parser) parseUnary() (Expr, error) {
	if p.at(tokNot) {
		p.next()
		x, err := p.parseUnary()
		if err != nil {
			return nil, err
		}
		return Not{X: x}, nil
	}
	return p.parseAtom()
}

func (p *parser) parseAtom() (Expr, error) {
	t := p.next()
	switch t.kind {
	case tokIdent:
		return Variable(t.text), nil
	case tokConst:
		return BoolConst(t.val), nil
	case tokLParen:
		e, err := p.parseIff()
		if err != nil {
			return nil, err
		}
		if !p.at(tokRParen) {
			return nil, fmt.Errorf("satsolver: missing closing parenthesis")
		}
		p.next()
		return e, nil
	}
	return nil, fmt.Errorf("satsolver: unexpected token in expression")
}
