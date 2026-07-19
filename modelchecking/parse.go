package modelchecking

import (
	"fmt"
	"strings"
)

type tokType int

const (
	tEOF tokType = iota
	tAtom
	tTrue
	tFalse
	tNot
	tAnd
	tOr
	tImplies
	tIff
	tLParen
	tRParen
	tLBrack
	tRBrack
	tComma
	tX
	tF
	tG
	tU
	tR
	tW
	tE
	tA
	tEX
	tEF
	tEG
	tAX
	tAF
	tAG
)

type token struct {
	typ tokType
	txt string
	pos int
}

var reservedWords = map[string]tokType{
	"true":  tTrue,
	"false": tFalse,
	"X":     tX,
	"F":     tF,
	"G":     tG,
	"U":     tU,
	"R":     tR,
	"W":     tW,
	"E":     tE,
	"A":     tA,
	"EX":    tEX,
	"EF":    tEF,
	"EG":    tEG,
	"AX":    tAX,
	"AF":    tAF,
	"AG":    tAG,
}

// tokenize converts s into a slice of tokens terminated by an EOF token.
func tokenize(s string) ([]token, error) {
	var toks []token
	i := 0
	for i < len(s) {
		c := s[i]
		switch {
		case c == ' ' || c == '\t' || c == '\n' || c == '\r':
			i++
		case c == '(':
			toks = append(toks, token{tLParen, "(", i})
			i++
		case c == ')':
			toks = append(toks, token{tRParen, ")", i})
			i++
		case c == '[':
			toks = append(toks, token{tLBrack, "[", i})
			i++
		case c == ']':
			toks = append(toks, token{tRBrack, "]", i})
			i++
		case c == ',':
			toks = append(toks, token{tComma, ",", i})
			i++
		case c == '!' || c == '~':
			toks = append(toks, token{tNot, "!", i})
			i++
		case c == '&':
			i++
			if i < len(s) && s[i] == '&' {
				i++
			}
			toks = append(toks, token{tAnd, "&", i})
		case c == '|':
			i++
			if i < len(s) && s[i] == '|' {
				i++
			}
			toks = append(toks, token{tOr, "|", i})
		case c == '-':
			if i+1 < len(s) && s[i+1] == '>' {
				toks = append(toks, token{tImplies, "->", i})
				i += 2
			} else {
				return nil, fmt.Errorf("modelchecking: unexpected '-' at position %d", i)
			}
		case c == '<':
			if i+2 < len(s) && s[i+1] == '-' && s[i+2] == '>' {
				toks = append(toks, token{tIff, "<->", i})
				i += 3
			} else {
				return nil, fmt.Errorf("modelchecking: unexpected '<' at position %d", i)
			}
		case isIdentStart(c):
			j := i + 1
			for j < len(s) && isIdentPart(s[j]) {
				j++
			}
			word := s[i:j]
			if kw, ok := reservedWords[word]; ok {
				toks = append(toks, token{kw, word, i})
			} else {
				toks = append(toks, token{tAtom, word, i})
			}
			i = j
		default:
			return nil, fmt.Errorf("modelchecking: unexpected character %q at position %d", string(c), i)
		}
	}
	toks = append(toks, token{tEOF, "", len(s)})
	return toks, nil
}

func isIdentStart(c byte) bool {
	return c == '_' || (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z')
}

func isIdentPart(c byte) bool {
	return isIdentStart(c) || (c >= '0' && c <= '9') || c == '.' || c == '\''
}

// ltlParser is a recursive-descent parser over a token stream.
type ltlParser struct {
	toks []token
	pos  int
}

func (p *ltlParser) peek() token { return p.toks[p.pos] }
func (p *ltlParser) next() token { t := p.toks[p.pos]; p.pos++; return t }
func (p *ltlParser) atEnd() bool { return p.toks[p.pos].typ == tEOF }
func (p *ltlParser) accept(t tokType) bool {
	if p.toks[p.pos].typ == t {
		p.pos++
		return true
	}
	return false
}

// ParseLTL parses an LTL formula from its textual representation and returns the
// syntax tree. The accepted grammar uses the operators ! (not), & (and),
// | (or), -> (implies), <-> (iff), X (next), F (eventually), G (globally),
// U (until), R (release) and W (weak until), with the usual precedences and
// parentheses for grouping.
func ParseLTL(s string) (*LTL, error) {
	toks, err := tokenize(s)
	if err != nil {
		return nil, err
	}
	p := &ltlParser{toks: toks}
	f, err := p.parseIff()
	if err != nil {
		return nil, err
	}
	if !p.atEnd() {
		return nil, fmt.Errorf("modelchecking: trailing tokens starting at position %d (%q)", p.peek().pos, p.peek().txt)
	}
	return f, nil
}

// MustParseLTL parses s and panics on error. It is intended for constant
// formulas in tests and examples.
func MustParseLTL(s string) *LTL {
	f, err := ParseLTL(s)
	if err != nil {
		panic(err)
	}
	return f
}

func (p *ltlParser) parseIff() (*LTL, error) {
	l, err := p.parseImplies()
	if err != nil {
		return nil, err
	}
	for p.accept(tIff) {
		r, err := p.parseImplies()
		if err != nil {
			return nil, err
		}
		l = LTLIff(l, r)
	}
	return l, nil
}

func (p *ltlParser) parseImplies() (*LTL, error) {
	l, err := p.parseOr()
	if err != nil {
		return nil, err
	}
	if p.accept(tImplies) {
		r, err := p.parseImplies()
		if err != nil {
			return nil, err
		}
		return LTLImplies(l, r), nil
	}
	return l, nil
}

func (p *ltlParser) parseOr() (*LTL, error) {
	l, err := p.parseAnd()
	if err != nil {
		return nil, err
	}
	for p.accept(tOr) {
		r, err := p.parseAnd()
		if err != nil {
			return nil, err
		}
		l = LTLOr(l, r)
	}
	return l, nil
}

func (p *ltlParser) parseAnd() (*LTL, error) {
	l, err := p.parseUntil()
	if err != nil {
		return nil, err
	}
	for p.accept(tAnd) {
		r, err := p.parseUntil()
		if err != nil {
			return nil, err
		}
		l = LTLAnd(l, r)
	}
	return l, nil
}

func (p *ltlParser) parseUntil() (*LTL, error) {
	l, err := p.parseUnary()
	if err != nil {
		return nil, err
	}
	switch p.peek().typ {
	case tU:
		p.next()
		r, err := p.parseUntil()
		if err != nil {
			return nil, err
		}
		return LTLUntil(l, r), nil
	case tR:
		p.next()
		r, err := p.parseUntil()
		if err != nil {
			return nil, err
		}
		return LTLRelease(l, r), nil
	case tW:
		p.next()
		r, err := p.parseUntil()
		if err != nil {
			return nil, err
		}
		return LTLWeakUntil(l, r), nil
	}
	return l, nil
}

func (p *ltlParser) parseUnary() (*LTL, error) {
	switch p.peek().typ {
	case tNot:
		p.next()
		f, err := p.parseUnary()
		if err != nil {
			return nil, err
		}
		return LTLNot(f), nil
	case tX:
		p.next()
		f, err := p.parseUnary()
		if err != nil {
			return nil, err
		}
		return LTLNext(f), nil
	case tF:
		p.next()
		f, err := p.parseUnary()
		if err != nil {
			return nil, err
		}
		return LTLEventually(f), nil
	case tG:
		p.next()
		f, err := p.parseUnary()
		if err != nil {
			return nil, err
		}
		return LTLGlobally(f), nil
	}
	return p.parsePrimary()
}

func (p *ltlParser) parsePrimary() (*LTL, error) {
	t := p.next()
	switch t.typ {
	case tTrue:
		return LTLTop(), nil
	case tFalse:
		return LTLBot(), nil
	case tAtom:
		return LTLVar(t.txt), nil
	case tLParen:
		f, err := p.parseIff()
		if err != nil {
			return nil, err
		}
		if !p.accept(tRParen) {
			return nil, fmt.Errorf("modelchecking: expected ')' at position %d", p.peek().pos)
		}
		return f, nil
	}
	return nil, fmt.Errorf("modelchecking: unexpected token %q at position %d", t.txt, t.pos)
}

// ctlParser is a recursive-descent parser for CTL sharing the lexer.
type ctlParser struct {
	toks []token
	pos  int
}

func (p *ctlParser) peek() token { return p.toks[p.pos] }
func (p *ctlParser) next() token { t := p.toks[p.pos]; p.pos++; return t }
func (p *ctlParser) atEnd() bool { return p.toks[p.pos].typ == tEOF }
func (p *ctlParser) accept(t tokType) bool {
	if p.toks[p.pos].typ == t {
		p.pos++
		return true
	}
	return false
}

// ParseCTL parses a CTL formula from text. The grammar uses the Boolean
// operators ! & | -> <->, the unary temporal operators EX EF EG AX AF AG, and
// the binary forms E[a U b], A[a U b], E[a R b] and A[a R b].
func ParseCTL(s string) (*CTL, error) {
	toks, err := tokenize(s)
	if err != nil {
		return nil, err
	}
	p := &ctlParser{toks: toks}
	f, err := p.parseIff()
	if err != nil {
		return nil, err
	}
	if !p.atEnd() {
		return nil, fmt.Errorf("modelchecking: trailing tokens starting at position %d (%q)", p.peek().pos, p.peek().txt)
	}
	return f, nil
}

// MustParseCTL parses s and panics on error.
func MustParseCTL(s string) *CTL {
	f, err := ParseCTL(s)
	if err != nil {
		panic(err)
	}
	return f
}

func (p *ctlParser) parseIff() (*CTL, error) {
	l, err := p.parseImplies()
	if err != nil {
		return nil, err
	}
	for p.accept(tIff) {
		r, err := p.parseImplies()
		if err != nil {
			return nil, err
		}
		l = CTLIff(l, r)
	}
	return l, nil
}

func (p *ctlParser) parseImplies() (*CTL, error) {
	l, err := p.parseOr()
	if err != nil {
		return nil, err
	}
	if p.accept(tImplies) {
		r, err := p.parseImplies()
		if err != nil {
			return nil, err
		}
		return CTLImplies(l, r), nil
	}
	return l, nil
}

func (p *ctlParser) parseOr() (*CTL, error) {
	l, err := p.parseAnd()
	if err != nil {
		return nil, err
	}
	for p.accept(tOr) {
		r, err := p.parseAnd()
		if err != nil {
			return nil, err
		}
		l = CTLOr(l, r)
	}
	return l, nil
}

func (p *ctlParser) parseAnd() (*CTL, error) {
	l, err := p.parseUnary()
	if err != nil {
		return nil, err
	}
	for p.accept(tAnd) {
		r, err := p.parseUnary()
		if err != nil {
			return nil, err
		}
		l = CTLAnd(l, r)
	}
	return l, nil
}

func (p *ctlParser) parseUnary() (*CTL, error) {
	switch p.peek().typ {
	case tNot:
		p.next()
		f, err := p.parseUnary()
		if err != nil {
			return nil, err
		}
		return CTLNot(f), nil
	case tEX:
		p.next()
		f, err := p.parseUnary()
		if err != nil {
			return nil, err
		}
		return EX(f), nil
	case tEF:
		p.next()
		f, err := p.parseUnary()
		if err != nil {
			return nil, err
		}
		return EF(f), nil
	case tEG:
		p.next()
		f, err := p.parseUnary()
		if err != nil {
			return nil, err
		}
		return EG(f), nil
	case tAX:
		p.next()
		f, err := p.parseUnary()
		if err != nil {
			return nil, err
		}
		return AX(f), nil
	case tAF:
		p.next()
		f, err := p.parseUnary()
		if err != nil {
			return nil, err
		}
		return AF(f), nil
	case tAG:
		p.next()
		f, err := p.parseUnary()
		if err != nil {
			return nil, err
		}
		return AG(f), nil
	case tE:
		p.next()
		return p.parseBracket(true)
	case tA:
		p.next()
		return p.parseBracket(false)
	}
	return p.parsePrimary()
}

// parseBracket parses the body of E[..] or A[..]; exists selects the path
// quantifier.
func (p *ctlParser) parseBracket(exists bool) (*CTL, error) {
	if !p.accept(tLBrack) {
		return nil, fmt.Errorf("modelchecking: expected '[' after path quantifier at position %d", p.peek().pos)
	}
	a, err := p.parseIff()
	if err != nil {
		return nil, err
	}
	var isUntil bool
	switch p.peek().typ {
	case tU:
		isUntil = true
		p.next()
	case tR:
		isUntil = false
		p.next()
	default:
		return nil, fmt.Errorf("modelchecking: expected 'U' or 'R' inside path formula at position %d", p.peek().pos)
	}
	b, err := p.parseIff()
	if err != nil {
		return nil, err
	}
	if !p.accept(tRBrack) {
		return nil, fmt.Errorf("modelchecking: expected ']' at position %d", p.peek().pos)
	}
	switch {
	case exists && isUntil:
		return EU(a, b), nil
	case exists && !isUntil:
		return ER(a, b), nil
	case !exists && isUntil:
		return AU(a, b), nil
	default:
		return AR(a, b), nil
	}
}

func (p *ctlParser) parsePrimary() (*CTL, error) {
	t := p.next()
	switch t.typ {
	case tTrue:
		return CTLTop(), nil
	case tFalse:
		return CTLBot(), nil
	case tAtom:
		return CTLVar(t.txt), nil
	case tLParen:
		f, err := p.parseIff()
		if err != nil {
			return nil, err
		}
		if !p.accept(tRParen) {
			return nil, fmt.Errorf("modelchecking: expected ')' at position %d", p.peek().pos)
		}
		return f, nil
	}
	return nil, fmt.Errorf("modelchecking: unexpected token %q at position %d", t.txt, t.pos)
}

// normalizeSpaces collapses runs of whitespace for stable formula printing in
// tests.
func normalizeSpaces(s string) string {
	return strings.Join(strings.Fields(s), " ")
}
