package logic

import "strings"

// ToNNF converts e to negation normal form: an equivalent expression built only
// from And, Or and negations applied directly to variables or constants. All
// derived connectives (implication, equivalence, exclusive-or and the negated
// gates) are expanded and negations are driven inward with the De Morgan and
// double-negation laws.
func ToNNF(e Expr) Expr {
	return logicNNF(e, false)
}

// logicNNF returns the negation normal form of e, negated when neg is true.
func logicNNF(e Expr, neg bool) Expr {
	switch t := e.(type) {
	case Const:
		if neg {
			return NewConst(!bool(t))
		}
		return t
	case Var:
		if neg {
			return NewNot(t)
		}
		return t
	case *UnaryExpr: // NotOp
		return logicNNF(t.X, !neg)
	case *BinaryExpr:
		return logicNNFBinary(t, neg)
	default:
		return e
	}
}

// logicNNFBinary handles the binary-operator cases of the NNF transform.
func logicNNFBinary(t *BinaryExpr, neg bool) Expr {
	l, r := t.L, t.R
	switch t.Op {
	case AndOp:
		if neg {
			return NewOr(logicNNF(l, true), logicNNF(r, true))
		}
		return NewAnd(logicNNF(l, false), logicNNF(r, false))
	case OrOp:
		if neg {
			return NewAnd(logicNNF(l, true), logicNNF(r, true))
		}
		return NewOr(logicNNF(l, false), logicNNF(r, false))
	case NandOp: // !(l & r)
		return logicNNFBinary(&BinaryExpr{Op: AndOp, L: l, R: r}, !neg)
	case NorOp: // !(l | r)
		return logicNNFBinary(&BinaryExpr{Op: OrOp, L: l, R: r}, !neg)
	case ImpliesOp: // !l | r
		if neg {
			return NewAnd(logicNNF(l, false), logicNNF(r, true))
		}
		return NewOr(logicNNF(l, true), logicNNF(r, false))
	case IffOp:
		return logicNNFEquiv(l, r, neg)
	case XnorOp:
		return logicNNFEquiv(l, r, neg)
	case XorOp:
		return logicNNFEquiv(l, r, !neg)
	default:
		return t
	}
}

// logicNNFEquiv expands the equivalence l<->r (or its negation) into NNF.
func logicNNFEquiv(l, r Expr, neg bool) Expr {
	if neg { // exclusive-or: (l & !r) | (!l & r)
		return NewOr(
			NewAnd(logicNNF(l, false), logicNNF(r, true)),
			NewAnd(logicNNF(l, true), logicNNF(r, false)),
		)
	}
	// equivalence: (l & r) | (!l & !r)
	return NewOr(
		NewAnd(logicNNF(l, false), logicNNF(r, false)),
		NewAnd(logicNNF(l, true), logicNNF(r, true)),
	)
}

// logicLiteral builds the literal for variable v that is true when the assigned
// bit equals want.
func logicLiteral(v string, want bool) Expr {
	if want {
		return NewVar(v)
	}
	return NewNot(NewVar(v))
}

// ToDNF returns the canonical (full) disjunctive normal form of e: an Or of
// minterms, one conjunction of literals for each satisfying assignment. A
// contradiction yields the constant F; a tautology over no variables yields T.
func ToDNF(e Expr) Expr {
	vars := Vars(e)
	tt := NewTruthTable(e)
	var terms []Expr
	for i, row := range tt.Rows {
		if !row.Result {
			continue
		}
		terms = append(terms, logicMintermExpr(i, vars))
	}
	if len(terms) == 0 {
		return NewConst(false)
	}
	return logicOrChain(terms)
}

// ToCNF returns the canonical (full) conjunctive normal form of e: an And of
// maxterms, one disjunction of literals for each falsifying assignment. A
// tautology yields the constant T; a contradiction over no variables yields F.
func ToCNF(e Expr) Expr {
	vars := Vars(e)
	tt := NewTruthTable(e)
	var clauses []Expr
	for i, row := range tt.Rows {
		if row.Result {
			continue
		}
		clauses = append(clauses, logicMaxtermExpr(i, vars))
	}
	if len(clauses) == 0 {
		return NewConst(true)
	}
	return logicAndChain(clauses)
}

// logicMintermExpr builds the conjunction of literals true only at assignment i.
func logicMintermExpr(i int, vars []string) Expr {
	env := IndexToAssignment(i, vars)
	if len(vars) == 0 {
		return NewConst(true)
	}
	var lits []Expr
	for _, v := range vars {
		lits = append(lits, logicLiteral(v, env[v]))
	}
	return logicAndChain(lits)
}

// logicMaxtermExpr builds the disjunction of literals false only at assignment i.
func logicMaxtermExpr(i int, vars []string) Expr {
	env := IndexToAssignment(i, vars)
	if len(vars) == 0 {
		return NewConst(false)
	}
	var lits []Expr
	for _, v := range vars {
		// The clause must be false at env, so each literal is false there.
		lits = append(lits, logicLiteral(v, !env[v]))
	}
	return logicOrChain(lits)
}

// logicOrChain folds terms into a right-nested Or; it requires len(terms) >= 1.
func logicOrChain(terms []Expr) Expr {
	acc := terms[len(terms)-1]
	for i := len(terms) - 2; i >= 0; i-- {
		acc = NewOr(terms[i], acc)
	}
	return acc
}

// logicAndChain folds terms into a right-nested And; it requires len(terms) >= 1.
func logicAndChain(terms []Expr) Expr {
	acc := terms[len(terms)-1]
	for i := len(terms) - 2; i >= 0; i-- {
		acc = NewAnd(terms[i], acc)
	}
	return acc
}

// DNFString renders the canonical disjunctive normal form of e as a sum-of-
// products string, using juxtaposition-free "&"/"|" notation with negation
// written as a leading "!". A contradiction renders as "F".
func DNFString(e Expr) string {
	vars := Vars(e)
	tt := NewTruthTable(e)
	var terms []string
	for i, row := range tt.Rows {
		if !row.Result {
			continue
		}
		terms = append(terms, logicMintermString(i, vars))
	}
	if len(terms) == 0 {
		return "F"
	}
	return strings.Join(terms, " | ")
}

// CNFString renders the canonical conjunctive normal form of e as a product-of-
// sums string. A tautology renders as "T".
func CNFString(e Expr) string {
	vars := Vars(e)
	tt := NewTruthTable(e)
	var clauses []string
	for i, row := range tt.Rows {
		if row.Result {
			continue
		}
		clauses = append(clauses, logicMaxtermString(i, vars))
	}
	if len(clauses) == 0 {
		return "T"
	}
	return strings.Join(clauses, " & ")
}

// logicMintermString renders one product term of the DNF.
func logicMintermString(i int, vars []string) string {
	if len(vars) == 0 {
		return "T"
	}
	env := IndexToAssignment(i, vars)
	var b strings.Builder
	for j, v := range vars {
		if j > 0 {
			b.WriteByte('&')
		}
		if !env[v] {
			b.WriteByte('!')
		}
		b.WriteString(v)
	}
	return b.String()
}

// logicMaxtermString renders one sum clause of the CNF.
func logicMaxtermString(i int, vars []string) string {
	if len(vars) == 0 {
		return "F"
	}
	env := IndexToAssignment(i, vars)
	var b strings.Builder
	b.WriteByte('(')
	for j, v := range vars {
		if j > 0 {
			b.WriteByte('|')
		}
		if env[v] {
			b.WriteByte('!')
		}
		b.WriteString(v)
	}
	b.WriteByte(')')
	return b.String()
}

// Simplify returns an equivalent expression reduced by constant folding and the
// standard Boolean identities (identity, domination, idempotence and
// complementation). It performs local, semantics-preserving rewrites rather
// than full minimisation; use [MinimizeSOP] for optimal two-level results.
func Simplify(e Expr) Expr {
	switch t := e.(type) {
	case Const, Var:
		return e
	case *UnaryExpr:
		x := Simplify(t.X)
		if c, ok := x.(Const); ok {
			return NewConst(!bool(c))
		}
		if inner, ok := x.(*UnaryExpr); ok { // double negation
			return inner.X
		}
		return NewNot(x)
	case *BinaryExpr:
		return logicSimplifyBinary(t.Op, Simplify(t.L), Simplify(t.R))
	default:
		return e
	}
}

// logicSimplifyBinary applies identity rules to a binary node whose operands are
// already simplified.
func logicSimplifyBinary(op BinaryOp, l, r Expr) Expr {
	lc, lok := l.(Const)
	rc, rok := r.(Const)
	if lok && rok { // constant fold
		return NewConst(op.apply(bool(lc), bool(rc)))
	}
	eq := logicEqual(l, r)
	comp := logicIsNegation(l, r)
	switch op {
	case AndOp:
		if (lok && !bool(lc)) || (rok && !bool(rc)) {
			return NewConst(false)
		}
		if lok {
			return r
		}
		if rok {
			return l
		}
		if eq {
			return l
		}
		if comp {
			return NewConst(false)
		}
	case OrOp:
		if (lok && bool(lc)) || (rok && bool(rc)) {
			return NewConst(true)
		}
		if lok {
			return r
		}
		if rok {
			return l
		}
		if eq {
			return l
		}
		if comp {
			return NewConst(true)
		}
	case NandOp:
		return Simplify(NewNot(logicSimplifyBinary(AndOp, l, r)))
	case NorOp:
		return Simplify(NewNot(logicSimplifyBinary(OrOp, l, r)))
	case XorOp:
		if lok {
			return logicXorConst(bool(lc), r)
		}
		if rok {
			return logicXorConst(bool(rc), l)
		}
		if eq {
			return NewConst(false)
		}
		if comp {
			return NewConst(true)
		}
	case XnorOp:
		return Simplify(NewNot(logicSimplifyBinary(XorOp, l, r)))
	case ImpliesOp:
		if lok {
			if !bool(lc) {
				return NewConst(true)
			}
			return r
		}
		if rok {
			if bool(rc) {
				return NewConst(true)
			}
			return Simplify(NewNot(l))
		}
		if eq {
			return NewConst(true)
		}
	case IffOp:
		if lok {
			if bool(lc) {
				return r
			}
			return Simplify(NewNot(r))
		}
		if rok {
			if bool(rc) {
				return l
			}
			return Simplify(NewNot(l))
		}
		if eq {
			return NewConst(true)
		}
		if comp {
			return NewConst(false)
		}
	}
	return newBinary(op, l, r)
}

// logicXorConst simplifies c ^ x for a constant c.
func logicXorConst(c bool, x Expr) Expr {
	if c {
		return Simplify(NewNot(x))
	}
	return x
}

// logicEqual reports structural equality of two expressions.
func logicEqual(a, b Expr) bool {
	switch at := a.(type) {
	case Const:
		bt, ok := b.(Const)
		return ok && at == bt
	case Var:
		bt, ok := b.(Var)
		return ok && at == bt
	case *UnaryExpr:
		bt, ok := b.(*UnaryExpr)
		return ok && at.Op == bt.Op && logicEqual(at.X, bt.X)
	case *BinaryExpr:
		bt, ok := b.(*BinaryExpr)
		return ok && at.Op == bt.Op && logicEqual(at.L, bt.L) && logicEqual(at.R, bt.R)
	default:
		return false
	}
}

// logicIsNegation reports whether a and b are structural negations of each
// other, i.e. one is Not applied to the other.
func logicIsNegation(a, b Expr) bool {
	if u, ok := a.(*UnaryExpr); ok && u.Op == NotOp && logicEqual(u.X, b) {
		return true
	}
	if u, ok := b.(*UnaryExpr); ok && u.Op == NotOp && logicEqual(u.X, a) {
		return true
	}
	return false
}

// IsDNF reports whether e is in disjunctive normal form: a disjunction of one or
// more conjunctions of literals, where a literal is a variable, a negated
// variable or a constant.
func IsDNF(e Expr) bool {
	if b, ok := e.(*BinaryExpr); ok && b.Op == OrOp {
		return IsDNF(b.L) && IsDNF(b.R)
	}
	return logicIsProduct(e)
}

// IsCNF reports whether e is in conjunctive normal form: a conjunction of one or
// more disjunctions of literals.
func IsCNF(e Expr) bool {
	if b, ok := e.(*BinaryExpr); ok && b.Op == AndOp {
		return IsCNF(b.L) && IsCNF(b.R)
	}
	return logicIsSum(e)
}

// logicIsProduct reports whether e is a conjunction of literals.
func logicIsProduct(e Expr) bool {
	if b, ok := e.(*BinaryExpr); ok && b.Op == AndOp {
		return logicIsProduct(b.L) && logicIsProduct(b.R)
	}
	return logicIsLiteral(e)
}

// logicIsSum reports whether e is a disjunction of literals.
func logicIsSum(e Expr) bool {
	if b, ok := e.(*BinaryExpr); ok && b.Op == OrOp {
		return logicIsSum(b.L) && logicIsSum(b.R)
	}
	return logicIsLiteral(e)
}

// logicIsLiteral reports whether e is a variable, a negated variable or a
// constant.
func logicIsLiteral(e Expr) bool {
	switch t := e.(type) {
	case Var, Const:
		return true
	case *UnaryExpr:
		if t.Op != NotOp {
			return false
		}
		_, ok := t.X.(Var)
		return ok
	default:
		return false
	}
}
