package logic

import (
	"fmt"
	"sort"
	"strings"
)

// Expr is a node in a propositional-logic expression tree. The interface is
// sealed: only the [Var], [Const], [UnaryExpr] and [BinaryExpr] types defined
// in this package implement it, which lets the package reason exhaustively over
// the concrete node kinds.
type Expr interface {
	// Eval evaluates the expression under the given variable assignment. It
	// returns an error if the expression references a variable that env does
	// not define.
	Eval(env map[string]bool) (bool, error)
	// String renders the expression as a re-parseable string.
	String() string

	collectVars(m map[string]bool)
}

// UnaryOp identifies a unary connective. The only unary connective is [NotOp].
type UnaryOp int

// NotOp is the logical-negation connective.
const NotOp UnaryOp = iota

// String returns the source symbol for the unary operator.
func (op UnaryOp) String() string {
	if op == NotOp {
		return "!"
	}
	return fmt.Sprintf("UnaryOp(%d)", int(op))
}

// BinaryOp identifies a binary connective.
type BinaryOp int

// The binary connectives, ordered from loosest to tightest binding as parsed by
// [Parse].
const (
	IffOp     BinaryOp = iota // <-> biconditional (equivalence)
	ImpliesOp                 // -> material conditional
	OrOp                      // | inclusive disjunction
	NorOp                     // negated disjunction
	XorOp                     // ^ exclusive disjunction
	XnorOp                    // negated exclusive disjunction
	AndOp                     // & conjunction
	NandOp                    // negated conjunction
)

// String returns the canonical source symbol for the binary operator.
func (op BinaryOp) String() string {
	switch op {
	case IffOp:
		return "<->"
	case ImpliesOp:
		return "->"
	case OrOp:
		return "|"
	case NorOp:
		return "nor"
	case XorOp:
		return "^"
	case XnorOp:
		return "xnor"
	case AndOp:
		return "&"
	case NandOp:
		return "nand"
	default:
		return fmt.Sprintf("BinaryOp(%d)", int(op))
	}
}

// apply computes the operator's truth function on two operand values.
func (op BinaryOp) apply(a, b bool) bool {
	switch op {
	case IffOp:
		return Iff(a, b)
	case ImpliesOp:
		return Implies(a, b)
	case OrOp:
		return Or(a, b)
	case NorOp:
		return Nor(a, b)
	case XorOp:
		return Xor(a, b)
	case XnorOp:
		return Xnor(a, b)
	case AndOp:
		return And(a, b)
	case NandOp:
		return Nand(a, b)
	default:
		return false
	}
}

// Var is a propositional variable identified by its name.
type Var string

// NewVar returns the variable with the given name.
func NewVar(name string) Var { return Var(name) }

// Eval looks the variable up in env, returning an error when it is unbound.
func (v Var) Eval(env map[string]bool) (bool, error) {
	b, ok := env[string(v)]
	if !ok {
		return false, fmt.Errorf("logic: unbound variable %q", string(v))
	}
	return b, nil
}

// String returns the variable's name.
func (v Var) String() string { return string(v) }

func (v Var) collectVars(m map[string]bool) { m[string(v)] = true }

// Const is a Boolean constant, either true or false.
type Const bool

// True and False are the Boolean constant expressions.
const (
	True  Const = true
	False Const = false
)

// NewConst returns the constant expression for b.
func NewConst(b bool) Const { return Const(b) }

// Eval returns the constant's value; it never fails.
func (c Const) Eval(map[string]bool) (bool, error) { return bool(c), nil }

// String renders the constant as "T" or "F".
func (c Const) String() string {
	if bool(c) {
		return "T"
	}
	return "F"
}

func (c Const) collectVars(map[string]bool) {}

// UnaryExpr applies a [UnaryOp] to a single operand.
type UnaryExpr struct {
	Op UnaryOp // the connective, always NotOp
	X  Expr    // the operand
}

// NewNot returns the negation of x.
func NewNot(x Expr) *UnaryExpr { return &UnaryExpr{Op: NotOp, X: x} }

// Eval evaluates the operand and applies the unary connective.
func (e *UnaryExpr) Eval(env map[string]bool) (bool, error) {
	x, err := e.X.Eval(env)
	if err != nil {
		return false, err
	}
	return Not(x), nil
}

// String renders the negation, parenthesising a compound operand.
func (e *UnaryExpr) String() string {
	if _, ok := e.X.(*BinaryExpr); ok {
		return e.Op.String() + "(" + e.X.String() + ")"
	}
	return e.Op.String() + e.X.String()
}

func (e *UnaryExpr) collectVars(m map[string]bool) { e.X.collectVars(m) }

// BinaryExpr applies a [BinaryOp] to a left and right operand.
type BinaryExpr struct {
	Op   BinaryOp // the connective
	L, R Expr     // the operands
}

// newBinary builds a binary node.
func newBinary(op BinaryOp, l, r Expr) *BinaryExpr {
	return &BinaryExpr{Op: op, L: l, R: r}
}

// NewAnd returns the conjunction l & r.
func NewAnd(l, r Expr) *BinaryExpr { return newBinary(AndOp, l, r) }

// NewOr returns the disjunction l | r.
func NewOr(l, r Expr) *BinaryExpr { return newBinary(OrOp, l, r) }

// NewXor returns the exclusive disjunction l ^ r.
func NewXor(l, r Expr) *BinaryExpr { return newBinary(XorOp, l, r) }

// NewNand returns the negated conjunction of l and r.
func NewNand(l, r Expr) *BinaryExpr { return newBinary(NandOp, l, r) }

// NewNor returns the negated disjunction of l and r.
func NewNor(l, r Expr) *BinaryExpr { return newBinary(NorOp, l, r) }

// NewXnor returns the negated exclusive disjunction (equivalence) of l and r.
func NewXnor(l, r Expr) *BinaryExpr { return newBinary(XnorOp, l, r) }

// NewImplies returns the material conditional l -> r.
func NewImplies(l, r Expr) *BinaryExpr { return newBinary(ImpliesOp, l, r) }

// NewIff returns the biconditional l <-> r.
func NewIff(l, r Expr) *BinaryExpr { return newBinary(IffOp, l, r) }

// Eval evaluates both operands and applies the binary connective.
func (e *BinaryExpr) Eval(env map[string]bool) (bool, error) {
	l, err := e.L.Eval(env)
	if err != nil {
		return false, err
	}
	r, err := e.R.Eval(env)
	if err != nil {
		return false, err
	}
	return e.Op.apply(l, r), nil
}

// String renders the binary expression with explicit parentheses.
func (e *BinaryExpr) String() string {
	return "(" + e.L.String() + " " + e.Op.String() + " " + e.R.String() + ")"
}

func (e *BinaryExpr) collectVars(m map[string]bool) {
	e.L.collectVars(m)
	e.R.collectVars(m)
}

// Vars returns the sorted, de-duplicated set of variable names appearing in e.
func Vars(e Expr) []string {
	m := map[string]bool{}
	e.collectVars(m)
	out := make([]string, 0, len(m))
	for k := range m {
		out = append(out, k)
	}
	sort.Strings(out)
	return out
}

// EvalString parses s and evaluates it under env in a single call.
func EvalString(s string, env map[string]bool) (bool, error) {
	e, err := Parse(s)
	if err != nil {
		return false, err
	}
	return e.Eval(env)
}

// MustParse is like [Parse] but panics if s does not parse. It is intended for
// package-level expressions known to be valid and for tests.
func MustParse(s string) Expr {
	e, err := Parse(s)
	if err != nil {
		panic(err)
	}
	return e
}

// tokenKind enumerates the lexical categories produced by the tokenizer.
type tokenKind int

const (
	tkEOF tokenKind = iota
	tkVar
	tkConst
	tkNot
	tkAnd
	tkOr
	tkXor
	tkNand
	tkNor
	tkXnor
	tkImplies
	tkIff
	tkLParen
	tkRParen
)

// token is a single lexical unit.
type token struct {
	kind tokenKind
	text string
	val  bool // for tkConst
}

// keywordKinds maps word operators and constants to their token kind.
var keywordKinds = map[string]tokenKind{
	"not": tkNot, "and": tkAnd, "or": tkOr, "xor": tkXor,
	"nand": tkNand, "nor": tkNor, "xnor": tkXnor,
	"implies": tkImplies, "iff": tkIff, "eq": tkIff,
}

// logicTokenize converts the source string into a token slice.
func logicTokenize(s string) ([]token, error) {
	var toks []token
	i := 0
	for i < len(s) {
		c := s[i]
		switch {
		case c == ' ' || c == '\t' || c == '\n' || c == '\r':
			i++
		case c == '(':
			toks = append(toks, token{kind: tkLParen})
			i++
		case c == ')':
			toks = append(toks, token{kind: tkRParen})
			i++
		case strings.HasPrefix(s[i:], "<->"):
			toks = append(toks, token{kind: tkIff})
			i += 3
		case strings.HasPrefix(s[i:], "<=>"):
			toks = append(toks, token{kind: tkIff})
			i += 3
		case strings.HasPrefix(s[i:], "->"):
			toks = append(toks, token{kind: tkImplies})
			i += 2
		case strings.HasPrefix(s[i:], "=>"):
			toks = append(toks, token{kind: tkImplies})
			i += 2
		case strings.HasPrefix(s[i:], "=="):
			toks = append(toks, token{kind: tkIff})
			i += 2
		case strings.HasPrefix(s[i:], "&&"):
			toks = append(toks, token{kind: tkAnd})
			i += 2
		case strings.HasPrefix(s[i:], "||"):
			toks = append(toks, token{kind: tkOr})
			i += 2
		case c == '!' || c == '~':
			toks = append(toks, token{kind: tkNot})
			i++
		case c == '&' || c == '*' || c == '.':
			toks = append(toks, token{kind: tkAnd})
			i++
		case c == '|' || c == '+':
			toks = append(toks, token{kind: tkOr})
			i++
		case c == '^':
			toks = append(toks, token{kind: tkXor})
			i++
		case c == '=':
			toks = append(toks, token{kind: tkIff})
			i++
		case logicIsIdentStart(c):
			j := i + 1
			for j < len(s) && logicIsIdentPart(s[j]) {
				j++
			}
			word := s[i:j]
			i = j
			lower := strings.ToLower(word)
			if k, ok := keywordKinds[lower]; ok {
				toks = append(toks, token{kind: k, text: word})
				continue
			}
			switch word {
			case "true", "True", "TRUE", "T", "1":
				toks = append(toks, token{kind: tkConst, val: true, text: word})
			case "false", "False", "FALSE", "F", "0":
				toks = append(toks, token{kind: tkConst, val: false, text: word})
			default:
				toks = append(toks, token{kind: tkVar, text: word})
			}
		case c == '0' || c == '1':
			toks = append(toks, token{kind: tkConst, val: c == '1', text: string(c)})
			i++
		default:
			return nil, fmt.Errorf("logic: unexpected character %q at position %d", string(c), i)
		}
	}
	toks = append(toks, token{kind: tkEOF})
	return toks, nil
}

// logicIsIdentStart reports whether c may begin an identifier.
func logicIsIdentStart(c byte) bool {
	return c == '_' || (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z')
}

// logicIsIdentPart reports whether c may continue an identifier.
func logicIsIdentPart(c byte) bool {
	return logicIsIdentStart(c) || (c >= '0' && c <= '9')
}

// logicParser is a recursive-descent parser over a token slice.
type logicParser struct {
	toks []token
	pos  int
}

// peek returns the current token without consuming it.
func (p *logicParser) peek() token { return p.toks[p.pos] }

// next consumes and returns the current token.
func (p *logicParser) next() token {
	t := p.toks[p.pos]
	p.pos++
	return t
}

// Parse parses a propositional-logic expression. See the package documentation
// for the accepted grammar. Binding tightens from <-> through ->, |, ^, & to
// the unary !; the word forms and, or, not, xor, nand, nor, xnor, implies and
// iff are accepted alongside the symbols.
func Parse(s string) (Expr, error) {
	toks, err := logicTokenize(s)
	if err != nil {
		return nil, err
	}
	p := &logicParser{toks: toks}
	e, err := p.parseIff()
	if err != nil {
		return nil, err
	}
	if p.peek().kind != tkEOF {
		return nil, fmt.Errorf("logic: unexpected trailing input near %q", p.peek().text)
	}
	return e, nil
}

// parseIff handles the loosest-binding biconditional (left associative).
func (p *logicParser) parseIff() (Expr, error) {
	left, err := p.parseImplies()
	if err != nil {
		return nil, err
	}
	for p.peek().kind == tkIff {
		p.next()
		right, err := p.parseImplies()
		if err != nil {
			return nil, err
		}
		left = NewIff(left, right)
	}
	return left, nil
}

// parseImplies handles the material conditional (right associative).
func (p *logicParser) parseImplies() (Expr, error) {
	left, err := p.parseOr()
	if err != nil {
		return nil, err
	}
	if p.peek().kind == tkImplies {
		p.next()
		right, err := p.parseImplies()
		if err != nil {
			return nil, err
		}
		return NewImplies(left, right), nil
	}
	return left, nil
}

// parseOr handles inclusive and negated disjunction (left associative).
func (p *logicParser) parseOr() (Expr, error) {
	left, err := p.parseXor()
	if err != nil {
		return nil, err
	}
	for {
		switch p.peek().kind {
		case tkOr:
			p.next()
			right, err := p.parseXor()
			if err != nil {
				return nil, err
			}
			left = NewOr(left, right)
		case tkNor:
			p.next()
			right, err := p.parseXor()
			if err != nil {
				return nil, err
			}
			left = NewNor(left, right)
		default:
			return left, nil
		}
	}
}

// parseXor handles exclusive and negated-exclusive disjunction.
func (p *logicParser) parseXor() (Expr, error) {
	left, err := p.parseAnd()
	if err != nil {
		return nil, err
	}
	for {
		switch p.peek().kind {
		case tkXor:
			p.next()
			right, err := p.parseAnd()
			if err != nil {
				return nil, err
			}
			left = NewXor(left, right)
		case tkXnor:
			p.next()
			right, err := p.parseAnd()
			if err != nil {
				return nil, err
			}
			left = NewXnor(left, right)
		default:
			return left, nil
		}
	}
}

// parseAnd handles conjunction and negated conjunction (tightest binary level).
func (p *logicParser) parseAnd() (Expr, error) {
	left, err := p.parseNot()
	if err != nil {
		return nil, err
	}
	for {
		switch p.peek().kind {
		case tkAnd:
			p.next()
			right, err := p.parseNot()
			if err != nil {
				return nil, err
			}
			left = NewAnd(left, right)
		case tkNand:
			p.next()
			right, err := p.parseNot()
			if err != nil {
				return nil, err
			}
			left = NewNand(left, right)
		default:
			return left, nil
		}
	}
}

// parseNot handles unary negation.
func (p *logicParser) parseNot() (Expr, error) {
	if p.peek().kind == tkNot {
		p.next()
		x, err := p.parseNot()
		if err != nil {
			return nil, err
		}
		return NewNot(x), nil
	}
	return p.parsePrimary()
}

// parsePrimary handles constants, variables and parenthesised subexpressions.
func (p *logicParser) parsePrimary() (Expr, error) {
	t := p.next()
	switch t.kind {
	case tkConst:
		return NewConst(t.val), nil
	case tkVar:
		return NewVar(t.text), nil
	case tkLParen:
		e, err := p.parseIff()
		if err != nil {
			return nil, err
		}
		if p.peek().kind != tkRParen {
			return nil, fmt.Errorf("logic: missing closing parenthesis")
		}
		p.next()
		return e, nil
	case tkEOF:
		return nil, fmt.Errorf("logic: unexpected end of input")
	default:
		return nil, fmt.Errorf("logic: unexpected token %q", t.text)
	}
}
