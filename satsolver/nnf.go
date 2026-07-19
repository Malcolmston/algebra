package satsolver

// EliminateImplications rewrites an expression so that it uses only the And,
// Or and Not connectives, expanding conditionals, biconditionals and the
// derived gates Xor, Nand, Nor and Xnor into that basis.
func EliminateImplications(e Expr) Expr {
	switch t := e.(type) {
	case Variable, BoolConst:
		return t
	case Not:
		return Not{X: EliminateImplications(t.X)}
	case And:
		return And{X: EliminateImplications(t.X), Y: EliminateImplications(t.Y)}
	case Or:
		return Or{X: EliminateImplications(t.X), Y: EliminateImplications(t.Y)}
	case Implies:
		return Or{X: Not{X: EliminateImplications(t.X)}, Y: EliminateImplications(t.Y)}
	case Iff:
		a := EliminateImplications(t.X)
		b := EliminateImplications(t.Y)
		return And{
			X: Or{X: Not{X: a}, Y: b},
			Y: Or{X: Not{X: b}, Y: a},
		}
	case Xnor:
		a := EliminateImplications(t.X)
		b := EliminateImplications(t.Y)
		return And{
			X: Or{X: Not{X: a}, Y: b},
			Y: Or{X: Not{X: b}, Y: a},
		}
	case Xor:
		a := EliminateImplications(t.X)
		b := EliminateImplications(t.Y)
		return Or{
			X: And{X: a, Y: Not{X: b}},
			Y: And{X: Not{X: a}, Y: b},
		}
	case Nand:
		a := EliminateImplications(t.X)
		b := EliminateImplications(t.Y)
		return Not{X: And{X: a, Y: b}}
	case Nor:
		a := EliminateImplications(t.X)
		b := EliminateImplications(t.Y)
		return Not{X: Or{X: a, Y: b}}
	}
	return e
}

// ToNNF converts an expression to negation normal form, in which negations
// apply only to variables and the tree uses solely And, Or and Not (over
// literals). Implications and biconditionals are first eliminated.
func ToNNF(e Expr) Expr {
	return pushNeg(EliminateImplications(e), false)
}

// pushNeg drives De Morgan's laws inward; neg records whether the current
// subtree is under an odd number of negations.
func pushNeg(e Expr, neg bool) Expr {
	switch t := e.(type) {
	case BoolConst:
		if neg {
			return BoolConst(!bool(t))
		}
		return t
	case Variable:
		if neg {
			return Not{X: t}
		}
		return t
	case Not:
		return pushNeg(t.X, !neg)
	case And:
		if neg {
			return Or{X: pushNeg(t.X, true), Y: pushNeg(t.Y, true)}
		}
		return And{X: pushNeg(t.X, false), Y: pushNeg(t.Y, false)}
	case Or:
		if neg {
			return And{X: pushNeg(t.X, true), Y: pushNeg(t.Y, true)}
		}
		return Or{X: pushNeg(t.X, false), Y: pushNeg(t.Y, false)}
	}
	// Should not happen once implications are eliminated.
	if neg {
		return Not{X: e}
	}
	return e
}

// IsNNF reports whether e is already in negation normal form: it contains only
// And, Or, Not, Variable and BoolConst nodes and every Not wraps a variable.
func IsNNF(e Expr) bool {
	switch t := e.(type) {
	case Variable, BoolConst:
		return true
	case Not:
		switch t.X.(type) {
		case Variable, BoolConst:
			return true
		}
		return false
	case And:
		return IsNNF(t.X) && IsNNF(t.Y)
	case Or:
		return IsNNF(t.X) && IsNNF(t.Y)
	}
	return false
}

// Simplify applies Boolean simplification rules — constant folding, negation,
// idempotence, complementation, identity, domination and absorption — repeatedly
// until the expression reaches a fixed point.
func Simplify(e Expr) Expr {
	for {
		s := simplifyOnce(e)
		if ExprEqual(s, e) {
			return s
		}
		e = s
	}
}

func simplifyOnce(e Expr) Expr {
	switch t := e.(type) {
	case Variable, BoolConst:
		return t
	case Not:
		x := simplifyOnce(t.X)
		if c, ok := x.(BoolConst); ok {
			return BoolConst(!bool(c))
		}
		if n, ok := x.(Not); ok {
			return n.X
		}
		return Not{X: x}
	case And:
		x := simplifyOnce(t.X)
		y := simplifyOnce(t.Y)
		return simplifyAnd(x, y)
	case Or:
		x := simplifyOnce(t.X)
		y := simplifyOnce(t.Y)
		return simplifyOr(x, y)
	case Xor:
		x := simplifyOnce(t.X)
		y := simplifyOnce(t.Y)
		return simplifyXor(x, y)
	case Implies:
		x := simplifyOnce(t.X)
		y := simplifyOnce(t.Y)
		return simplifyOr(Not{X: x}, y)
	case Iff:
		x := simplifyOnce(t.X)
		y := simplifyOnce(t.Y)
		return simplifyIff(x, y)
	case Xnor:
		x := simplifyOnce(t.X)
		y := simplifyOnce(t.Y)
		return simplifyIff(x, y)
	case Nand:
		x := simplifyOnce(t.X)
		y := simplifyOnce(t.Y)
		return simplifyOnce(Not{X: And{X: x, Y: y}})
	case Nor:
		x := simplifyOnce(t.X)
		y := simplifyOnce(t.Y)
		return simplifyOnce(Not{X: Or{X: x, Y: y}})
	}
	return e
}

func simplifyAnd(x, y Expr) Expr {
	if c, ok := x.(BoolConst); ok {
		if !bool(c) {
			return False
		}
		return y
	}
	if c, ok := y.(BoolConst); ok {
		if !bool(c) {
			return False
		}
		return x
	}
	if ExprEqual(x, y) {
		return x
	}
	if isComplement(x, y) {
		return False
	}
	return And{X: x, Y: y}
}

func simplifyOr(x, y Expr) Expr {
	if c, ok := x.(BoolConst); ok {
		if bool(c) {
			return True
		}
		return y
	}
	if c, ok := y.(BoolConst); ok {
		if bool(c) {
			return True
		}
		return x
	}
	if ExprEqual(x, y) {
		return x
	}
	if isComplement(x, y) {
		return True
	}
	return Or{X: x, Y: y}
}

func simplifyXor(x, y Expr) Expr {
	if c, ok := x.(BoolConst); ok {
		if bool(c) {
			return simplifyOnce(Not{X: y})
		}
		return y
	}
	if c, ok := y.(BoolConst); ok {
		if bool(c) {
			return simplifyOnce(Not{X: x})
		}
		return x
	}
	if ExprEqual(x, y) {
		return False
	}
	if isComplement(x, y) {
		return True
	}
	return Xor{X: x, Y: y}
}

func simplifyIff(x, y Expr) Expr {
	if c, ok := x.(BoolConst); ok {
		if bool(c) {
			return y
		}
		return simplifyOnce(Not{X: y})
	}
	if c, ok := y.(BoolConst); ok {
		if bool(c) {
			return x
		}
		return simplifyOnce(Not{X: x})
	}
	if ExprEqual(x, y) {
		return True
	}
	if isComplement(x, y) {
		return False
	}
	return Iff{X: x, Y: y}
}

// isComplement reports whether a and b are syntactic complements (one is the
// negation of the other).
func isComplement(a, b Expr) bool {
	if n, ok := a.(Not); ok && ExprEqual(n.X, b) {
		return true
	}
	if n, ok := b.(Not); ok && ExprEqual(n.X, a) {
		return true
	}
	return false
}
