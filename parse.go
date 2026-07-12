package algebra

import (
	"fmt"
	"math/big"
	"strings"
	"unicode"
)

// Parse parses a mathematical expression written in ordinary infix notation
// and returns its expression tree. It understands:
//
//   - the binary operators + - * / ^ (^ is right-associative),
//   - unary plus and minus,
//   - parentheses for grouping,
//   - integer and decimal-float literals,
//   - symbols (identifiers), the constants pi and E,
//   - function calls sin, cos, tan, exp, log, ln and sqrt, and
//   - implicit multiplication, so "2x", "2(x+1)" and "x y" all parse as
//     products.
//
// The returned expression is canonicalized by the constructors as it is built.
func Parse(s string) (Expr, error) {
	toks, err := tokenize(s)
	if err != nil {
		return nil, err
	}
	p := &parser{toks: toks}
	e, err := p.parseExpr(0)
	if err != nil {
		return nil, err
	}
	if p.peek().kind != tEOF {
		return nil, fmt.Errorf("algebra: unexpected %q at position %d", p.peek().text, p.peek().pos)
	}
	return e, nil
}

// MustParse is like [Parse] but panics on error. It is convenient for tests
// and package-level variables built from constant strings.
func MustParse(s string) Expr {
	e, err := Parse(s)
	if err != nil {
		panic(err)
	}
	return e
}

type tokenKind int

const (
	tEOF tokenKind = iota
	tNum
	tIdent
	tOp
	tLParen
	tRParen
)

type token struct {
	kind tokenKind
	text string
	pos  int
}

func tokenize(s string) ([]token, error) {
	var toks []token
	runes := []rune(s)
	i := 0
	for i < len(runes) {
		c := runes[i]
		switch {
		case unicode.IsSpace(c):
			i++
		case c >= '0' && c <= '9' || c == '.':
			start := i
			for i < len(runes) && (runes[i] >= '0' && runes[i] <= '9' || runes[i] == '.') {
				i++
			}
			toks = append(toks, token{tNum, string(runes[start:i]), start})
		case unicode.IsLetter(c) || c == '_':
			start := i
			for i < len(runes) && (unicode.IsLetter(runes[i]) || unicode.IsDigit(runes[i]) || runes[i] == '_') {
				i++
			}
			toks = append(toks, token{tIdent, string(runes[start:i]), start})
		case c == '+' || c == '-' || c == '*' || c == '/' || c == '^':
			toks = append(toks, token{tOp, string(c), i})
			i++
		case c == '(':
			toks = append(toks, token{tLParen, "(", i})
			i++
		case c == ')':
			toks = append(toks, token{tRParen, ")", i})
			i++
		default:
			return nil, fmt.Errorf("algebra: unexpected character %q at position %d", string(c), i)
		}
	}
	toks = append(toks, token{tEOF, "", len(runes)})
	return toks, nil
}

type parser struct {
	toks []token
	pos  int
}

func (p *parser) peek() token { return p.toks[p.pos] }
func (p *parser) next() token { t := p.toks[p.pos]; p.pos++; return t }

// Operator precedences.
const (
	precAdd   = 10
	precMul   = 20
	precUnary = 30
	precPow   = 40
)

// parseExpr is a precedence-climbing (Pratt) parser.
func (p *parser) parseExpr(minPrec int) (Expr, error) {
	var left Expr
	var err error

	switch t := p.peek(); {
	case t.kind == tOp && t.text == "-":
		p.next()
		operand, e := p.parseExpr(precUnary)
		if e != nil {
			return nil, e
		}
		left = neg(operand)
	case t.kind == tOp && t.text == "+":
		p.next()
		left, err = p.parseExpr(precUnary)
		if err != nil {
			return nil, err
		}
	default:
		left, err = p.parsePrimary()
		if err != nil {
			return nil, err
		}
	}

	for {
		op, prec, rightAssoc, implicit := p.peekBinary()
		if op == "" || prec < minPrec {
			break
		}
		if !implicit {
			p.next()
		}
		nextMin := prec + 1
		if rightAssoc {
			nextMin = prec
		}
		right, e := p.parseExpr(nextMin)
		if e != nil {
			return nil, e
		}
		left = applyOp(op, left, right)
	}
	return left, nil
}

// peekBinary classifies the current token as a binary operator, reporting an
// implicit multiplication when a primary directly follows another primary.
func (p *parser) peekBinary() (op string, prec int, rightAssoc, implicit bool) {
	t := p.peek()
	switch {
	case t.kind == tOp && (t.text == "+" || t.text == "-"):
		return t.text, precAdd, false, false
	case t.kind == tOp && (t.text == "*" || t.text == "/"):
		return t.text, precMul, false, false
	case t.kind == tOp && t.text == "^":
		return t.text, precPow, true, false
	case t.kind == tNum || t.kind == tIdent || t.kind == tLParen:
		return "*", precMul, false, true
	}
	return "", 0, false, false
}

func (p *parser) parsePrimary() (Expr, error) {
	t := p.peek()
	switch t.kind {
	case tNum:
		p.next()
		return parseNumber(t.text, t.pos)
	case tLParen:
		p.next()
		e, err := p.parseExpr(0)
		if err != nil {
			return nil, err
		}
		if p.peek().kind != tRParen {
			return nil, fmt.Errorf("algebra: expected ')' at position %d", p.peek().pos)
		}
		p.next()
		return e, nil
	case tIdent:
		p.next()
		if isFunctionName(t.text) && p.peek().kind == tLParen {
			p.next()
			arg, err := p.parseExpr(0)
			if err != nil {
				return nil, err
			}
			if p.peek().kind != tRParen {
				return nil, fmt.Errorf("algebra: expected ')' after %s( at position %d", t.text, p.peek().pos)
			}
			p.next()
			name := t.text
			if name == "ln" {
				name = "log"
			}
			return applyFn(name, arg), nil
		}
		switch t.text {
		case "pi":
			return Pi, nil
		case "E":
			return E, nil
		}
		return Sym(t.text), nil
	}
	return nil, fmt.Errorf("algebra: unexpected %q at position %d", t.text, t.pos)
}

func parseNumber(text string, pos int) (Expr, error) {
	if strings.Contains(text, ".") {
		f := new(big.Float)
		if _, _, err := f.Parse(text, 10); err != nil {
			return nil, fmt.Errorf("algebra: invalid number %q at position %d", text, pos)
		}
		v, _ := f.Float64()
		return Flt(v), nil
	}
	n, ok := new(big.Int).SetString(text, 10)
	if !ok {
		return nil, fmt.Errorf("algebra: invalid integer %q at position %d", text, pos)
	}
	return newInteger(n), nil
}

func isFunctionName(name string) bool {
	switch name {
	case "sin", "cos", "tan", "exp", "log", "ln", "sqrt":
		return true
	}
	return false
}

func applyOp(op string, l, r Expr) Expr {
	switch op {
	case "+":
		return Add(l, r)
	case "-":
		return Add(l, neg(r))
	case "*":
		return Mul(l, r)
	case "/":
		return Mul(l, Pow(r, Int(-1)))
	case "^":
		return Pow(l, r)
	}
	return l
}
