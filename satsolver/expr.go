package satsolver

import (
	"sort"
	"strings"
)

// Expr is the sealed interface implemented by every node of a Boolean
// expression tree. Concrete implementations are [BoolConst], [Variable],
// [Not], [And], [Or], [Xor], [Implies], [Iff], [Nand], [Nor] and [Xnor].
type Expr interface {
	// Eval evaluates the expression under a variable environment. Variables
	// absent from env are treated as false.
	Eval(env map[string]bool) bool
	// String renders the expression in fully parenthesised infix form.
	String() string
	isExpr()
}

// BoolConst is a Boolean literal constant, either true or false.
type BoolConst bool

// True is the constant true expression.
var True Expr = BoolConst(true)

// False is the constant false expression.
var False Expr = BoolConst(false)

// Constant returns the constant expression for b.
func Constant(b bool) Expr { return BoolConst(b) }

// Eval returns the constant's value.
func (c BoolConst) Eval(map[string]bool) bool { return bool(c) }

// String returns "T" or "F".
func (c BoolConst) String() string {
	if bool(c) {
		return "T"
	}
	return "F"
}
func (c BoolConst) isExpr() {}

// Variable is a named Boolean variable.
type Variable string

// V returns a variable expression with the given name.
func V(name string) Expr { return Variable(name) }

// Eval looks the variable up in env, defaulting to false.
func (v Variable) Eval(env map[string]bool) bool { return env[string(v)] }

// String returns the variable name.
func (v Variable) String() string { return string(v) }
func (v Variable) isExpr()        {}

// Not is logical negation of its operand.
type Not struct{ X Expr }

// NotE returns the negation of x, folding double negations and constants.
func NotE(x Expr) Expr {
	switch t := x.(type) {
	case BoolConst:
		return BoolConst(!bool(t))
	case Not:
		return t.X
	}
	return Not{X: x}
}

// Eval evaluates the negation.
func (n Not) Eval(env map[string]bool) bool { return !n.X.Eval(env) }

// String renders the negation as "~(x)".
func (n Not) String() string { return "~" + paren(n.X) }
func (n Not) isExpr()        {}

// And is logical conjunction of two operands.
type And struct{ X, Y Expr }

// AndE returns the conjunction of x and y.
func AndE(x, y Expr) Expr { return And{X: x, Y: y} }

// Eval evaluates the conjunction.
func (a And) Eval(env map[string]bool) bool { return a.X.Eval(env) && a.Y.Eval(env) }

// String renders the conjunction as "(x & y)".
func (a And) String() string { return "(" + a.X.String() + " & " + a.Y.String() + ")" }
func (a And) isExpr()        {}

// Or is logical inclusive disjunction of two operands.
type Or struct{ X, Y Expr }

// OrE returns the disjunction of x and y.
func OrE(x, y Expr) Expr { return Or{X: x, Y: y} }

// Eval evaluates the disjunction.
func (o Or) Eval(env map[string]bool) bool { return o.X.Eval(env) || o.Y.Eval(env) }

// String renders the disjunction as "(x | y)".
func (o Or) String() string { return "(" + o.X.String() + " | " + o.Y.String() + ")" }
func (o Or) isExpr()        {}

// Xor is exclusive disjunction of two operands.
type Xor struct{ X, Y Expr }

// XorE returns the exclusive-or of x and y.
func XorE(x, y Expr) Expr { return Xor{X: x, Y: y} }

// Eval evaluates the exclusive-or.
func (x Xor) Eval(env map[string]bool) bool { return x.X.Eval(env) != x.Y.Eval(env) }

// String renders the exclusive-or as "(x ^ y)".
func (x Xor) String() string { return "(" + x.X.String() + " ^ " + x.Y.String() + ")" }
func (x Xor) isExpr()        {}

// Implies is the material conditional X -> Y.
type Implies struct{ X, Y Expr }

// ImpliesE returns the conditional x -> y.
func ImpliesE(x, y Expr) Expr { return Implies{X: x, Y: y} }

// Eval evaluates the conditional.
func (i Implies) Eval(env map[string]bool) bool { return !i.X.Eval(env) || i.Y.Eval(env) }

// String renders the conditional as "(x -> y)".
func (i Implies) String() string { return "(" + i.X.String() + " -> " + i.Y.String() + ")" }
func (i Implies) isExpr()        {}

// Iff is the biconditional X <-> Y (logical equivalence).
type Iff struct{ X, Y Expr }

// IffE returns the biconditional x <-> y.
func IffE(x, y Expr) Expr { return Iff{X: x, Y: y} }

// Eval evaluates the biconditional.
func (i Iff) Eval(env map[string]bool) bool { return i.X.Eval(env) == i.Y.Eval(env) }

// String renders the biconditional as "(x <-> y)".
func (i Iff) String() string { return "(" + i.X.String() + " <-> " + i.Y.String() + ")" }
func (i Iff) isExpr()        {}

// Nand is the negated conjunction of two operands.
type Nand struct{ X, Y Expr }

// NandE returns the nand of x and y.
func NandE(x, y Expr) Expr { return Nand{X: x, Y: y} }

// Eval evaluates the nand.
func (n Nand) Eval(env map[string]bool) bool { return !(n.X.Eval(env) && n.Y.Eval(env)) }

// String renders the nand as "(x !& y)".
func (n Nand) String() string { return "(" + n.X.String() + " !& " + n.Y.String() + ")" }
func (n Nand) isExpr()        {}

// Nor is the negated disjunction of two operands.
type Nor struct{ X, Y Expr }

// NorE returns the nor of x and y.
func NorE(x, y Expr) Expr { return Nor{X: x, Y: y} }

// Eval evaluates the nor.
func (n Nor) Eval(env map[string]bool) bool { return !(n.X.Eval(env) || n.Y.Eval(env)) }

// String renders the nor as "(x !| y)".
func (n Nor) String() string { return "(" + n.X.String() + " !| " + n.Y.String() + ")" }
func (n Nor) isExpr()        {}

// Xnor is the negated exclusive disjunction (equivalence) of two operands. It
// is logically identical to [Iff] but retained as a distinct gate.
type Xnor struct{ X, Y Expr }

// XnorE returns the xnor of x and y.
func XnorE(x, y Expr) Expr { return Xnor{X: x, Y: y} }

// Eval evaluates the xnor.
func (n Xnor) Eval(env map[string]bool) bool { return n.X.Eval(env) == n.Y.Eval(env) }

// String renders the xnor as "(x ~^ y)".
func (n Xnor) String() string { return "(" + n.X.String() + " ~^ " + n.Y.String() + ")" }
func (n Xnor) isExpr()        {}

func paren(e Expr) string {
	switch e.(type) {
	case Variable, BoolConst:
		return e.String()
	}
	return e.String()
}

// AndAll returns the conjunction of all expressions. The empty conjunction is
// [True].
func AndAll(es ...Expr) Expr {
	if len(es) == 0 {
		return True
	}
	acc := es[0]
	for _, e := range es[1:] {
		acc = And{X: acc, Y: e}
	}
	return acc
}

// OrAll returns the disjunction of all expressions. The empty disjunction is
// [False].
func OrAll(es ...Expr) Expr {
	if len(es) == 0 {
		return False
	}
	acc := es[0]
	for _, e := range es[1:] {
		acc = Or{X: acc, Y: e}
	}
	return acc
}

// XorAll returns the exclusive-or of all expressions. The empty case is
// [False].
func XorAll(es ...Expr) Expr {
	if len(es) == 0 {
		return False
	}
	acc := es[0]
	for _, e := range es[1:] {
		acc = Xor{X: acc, Y: e}
	}
	return acc
}

// Vars returns the sorted set of distinct variable names occurring in e.
func Vars(e Expr) []string {
	seen := map[string]bool{}
	collectVars(e, seen)
	out := make([]string, 0, len(seen))
	for v := range seen {
		out = append(out, v)
	}
	sort.Strings(out)
	return out
}

func collectVars(e Expr, seen map[string]bool) {
	switch t := e.(type) {
	case Variable:
		seen[string(t)] = true
	case BoolConst:
	case Not:
		collectVars(t.X, seen)
	case And:
		collectVars(t.X, seen)
		collectVars(t.Y, seen)
	case Or:
		collectVars(t.X, seen)
		collectVars(t.Y, seen)
	case Xor:
		collectVars(t.X, seen)
		collectVars(t.Y, seen)
	case Implies:
		collectVars(t.X, seen)
		collectVars(t.Y, seen)
	case Iff:
		collectVars(t.X, seen)
		collectVars(t.Y, seen)
	case Nand:
		collectVars(t.X, seen)
		collectVars(t.Y, seen)
	case Nor:
		collectVars(t.X, seen)
		collectVars(t.Y, seen)
	case Xnor:
		collectVars(t.X, seen)
		collectVars(t.Y, seen)
	}
}

// Eval evaluates e under env; a convenience wrapper over the [Expr.Eval]
// method.
func Eval(e Expr, env map[string]bool) bool { return e.Eval(env) }

// Size returns the number of nodes in the expression tree.
func Size(e Expr) int {
	switch t := e.(type) {
	case Variable, BoolConst:
		return 1
	case Not:
		return 1 + Size(t.X)
	case And:
		return 1 + Size(t.X) + Size(t.Y)
	case Or:
		return 1 + Size(t.X) + Size(t.Y)
	case Xor:
		return 1 + Size(t.X) + Size(t.Y)
	case Implies:
		return 1 + Size(t.X) + Size(t.Y)
	case Iff:
		return 1 + Size(t.X) + Size(t.Y)
	case Nand:
		return 1 + Size(t.X) + Size(t.Y)
	case Nor:
		return 1 + Size(t.X) + Size(t.Y)
	case Xnor:
		return 1 + Size(t.X) + Size(t.Y)
	}
	return 0
}

// Depth returns the height of the expression tree; a leaf has depth 1.
func Depth(e Expr) int {
	switch t := e.(type) {
	case Variable, BoolConst:
		return 1
	case Not:
		return 1 + Depth(t.X)
	case And:
		return 1 + max2(Depth(t.X), Depth(t.Y))
	case Or:
		return 1 + max2(Depth(t.X), Depth(t.Y))
	case Xor:
		return 1 + max2(Depth(t.X), Depth(t.Y))
	case Implies:
		return 1 + max2(Depth(t.X), Depth(t.Y))
	case Iff:
		return 1 + max2(Depth(t.X), Depth(t.Y))
	case Nand:
		return 1 + max2(Depth(t.X), Depth(t.Y))
	case Nor:
		return 1 + max2(Depth(t.X), Depth(t.Y))
	case Xnor:
		return 1 + max2(Depth(t.X), Depth(t.Y))
	}
	return 0
}

func max2(a, b int) int {
	if a > b {
		return a
	}
	return b
}

// Substitute returns a copy of e with every occurrence of the variable named
// name replaced by the expression repl.
func Substitute(e Expr, name string, repl Expr) Expr {
	switch t := e.(type) {
	case Variable:
		if string(t) == name {
			return repl
		}
		return t
	case BoolConst:
		return t
	case Not:
		return Not{X: Substitute(t.X, name, repl)}
	case And:
		return And{X: Substitute(t.X, name, repl), Y: Substitute(t.Y, name, repl)}
	case Or:
		return Or{X: Substitute(t.X, name, repl), Y: Substitute(t.Y, name, repl)}
	case Xor:
		return Xor{X: Substitute(t.X, name, repl), Y: Substitute(t.Y, name, repl)}
	case Implies:
		return Implies{X: Substitute(t.X, name, repl), Y: Substitute(t.Y, name, repl)}
	case Iff:
		return Iff{X: Substitute(t.X, name, repl), Y: Substitute(t.Y, name, repl)}
	case Nand:
		return Nand{X: Substitute(t.X, name, repl), Y: Substitute(t.Y, name, repl)}
	case Nor:
		return Nor{X: Substitute(t.X, name, repl), Y: Substitute(t.Y, name, repl)}
	case Xnor:
		return Xnor{X: Substitute(t.X, name, repl), Y: Substitute(t.Y, name, repl)}
	}
	return e
}

// ExprEqual reports whether two expression trees are structurally identical.
func ExprEqual(a, b Expr) bool {
	switch ta := a.(type) {
	case Variable:
		tb, ok := b.(Variable)
		return ok && ta == tb
	case BoolConst:
		tb, ok := b.(BoolConst)
		return ok && ta == tb
	case Not:
		tb, ok := b.(Not)
		return ok && ExprEqual(ta.X, tb.X)
	case And:
		tb, ok := b.(And)
		return ok && ExprEqual(ta.X, tb.X) && ExprEqual(ta.Y, tb.Y)
	case Or:
		tb, ok := b.(Or)
		return ok && ExprEqual(ta.X, tb.X) && ExprEqual(ta.Y, tb.Y)
	case Xor:
		tb, ok := b.(Xor)
		return ok && ExprEqual(ta.X, tb.X) && ExprEqual(ta.Y, tb.Y)
	case Implies:
		tb, ok := b.(Implies)
		return ok && ExprEqual(ta.X, tb.X) && ExprEqual(ta.Y, tb.Y)
	case Iff:
		tb, ok := b.(Iff)
		return ok && ExprEqual(ta.X, tb.X) && ExprEqual(ta.Y, tb.Y)
	case Nand:
		tb, ok := b.(Nand)
		return ok && ExprEqual(ta.X, tb.X) && ExprEqual(ta.Y, tb.Y)
	case Nor:
		tb, ok := b.(Nor)
		return ok && ExprEqual(ta.X, tb.X) && ExprEqual(ta.Y, tb.Y)
	case Xnor:
		tb, ok := b.(Xnor)
		return ok && ExprEqual(ta.X, tb.X) && ExprEqual(ta.Y, tb.Y)
	}
	return false
}

// EvalString parses s and evaluates it under env, returning the value and any
// parse error.
func EvalString(s string, env map[string]bool) (bool, error) {
	e, err := Parse(s)
	if err != nil {
		return false, err
	}
	return e.Eval(env), nil
}

// CountNodesByKind returns a map from a short kind name ("var", "const", "and",
// ...) to the number of nodes of that kind in e.
func CountNodesByKind(e Expr) map[string]int {
	out := map[string]int{}
	var walk func(Expr)
	walk = func(x Expr) {
		switch t := x.(type) {
		case Variable:
			out["var"]++
		case BoolConst:
			out["const"]++
		case Not:
			out["not"]++
			walk(t.X)
		case And:
			out["and"]++
			walk(t.X)
			walk(t.Y)
		case Or:
			out["or"]++
			walk(t.X)
			walk(t.Y)
		case Xor:
			out["xor"]++
			walk(t.X)
			walk(t.Y)
		case Implies:
			out["implies"]++
			walk(t.X)
			walk(t.Y)
		case Iff:
			out["iff"]++
			walk(t.X)
			walk(t.Y)
		case Nand:
			out["nand"]++
			walk(t.X)
			walk(t.Y)
		case Nor:
			out["nor"]++
			walk(t.X)
			walk(t.Y)
		case Xnor:
			out["xnor"]++
			walk(t.X)
			walk(t.Y)
		}
	}
	walk(e)
	return out
}

// binaryChildren returns the two children of a binary node and true, or nil,
// nil, false for non-binary nodes.
func binaryChildren(e Expr) (Expr, Expr, bool) {
	switch t := e.(type) {
	case And:
		return t.X, t.Y, true
	case Or:
		return t.X, t.Y, true
	case Xor:
		return t.X, t.Y, true
	case Implies:
		return t.X, t.Y, true
	case Iff:
		return t.X, t.Y, true
	case Nand:
		return t.X, t.Y, true
	case Nor:
		return t.X, t.Y, true
	case Xnor:
		return t.X, t.Y, true
	}
	return nil, nil, false
}

// Leaves returns the sorted distinct variable names, identical to [Vars]; kept
// for readability at call sites that view the tree as a term.
func Leaves(e Expr) []string { return Vars(e) }

// PrettyString is an alias for the expression's String rendering, provided for
// symmetry with other renderers.
func PrettyString(e Expr) string { return strings.TrimSpace(e.String()) }
