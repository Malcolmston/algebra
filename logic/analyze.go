package logic

// IsTautology reports whether e is true under every assignment of its
// variables. A variable-free expression is a tautology exactly when it
// evaluates to true.
func IsTautology(e Expr) bool {
	tt := NewTruthTable(e)
	for _, r := range tt.Rows {
		if !r.Result {
			return false
		}
	}
	return true
}

// IsContradiction reports whether e is false under every assignment of its
// variables.
func IsContradiction(e Expr) bool {
	tt := NewTruthTable(e)
	for _, r := range tt.Rows {
		if r.Result {
			return false
		}
	}
	return true
}

// IsSatisfiable reports whether at least one assignment makes e true. It is the
// negation of [IsContradiction].
func IsSatisfiable(e Expr) bool {
	return !IsContradiction(e)
}

// IsContingency reports whether e is neither a tautology nor a contradiction,
// that is, whether its value depends on the assignment.
func IsContingency(e Expr) bool {
	tt := NewTruthTable(e)
	sawTrue, sawFalse := false, false
	for _, r := range tt.Rows {
		if r.Result {
			sawTrue = true
		} else {
			sawFalse = true
		}
		if sawTrue && sawFalse {
			return true
		}
	}
	return false
}

// FindModel returns a satisfying assignment for e and true, or a nil map and
// false when e is unsatisfiable. The search proceeds in ascending
// assignment-index order, so the returned model is deterministic.
func FindModel(e Expr) (map[string]bool, bool) {
	tt := NewTruthTable(e)
	for _, r := range tt.Rows {
		if r.Result {
			return r.Values, true
		}
	}
	return nil, false
}

// AllModels returns every satisfying assignment of e, in ascending
// assignment-index order.
func AllModels(e Expr) []map[string]bool {
	tt := NewTruthTable(e)
	var out []map[string]bool
	for _, r := range tt.Rows {
		if r.Result {
			out = append(out, r.Values)
		}
	}
	return out
}

// CountModels returns the number of assignments that satisfy e.
func CountModels(e Expr) int {
	tt := NewTruthTable(e)
	n := 0
	for _, r := range tt.Rows {
		if r.Result {
			n++
		}
	}
	return n
}

// Equivalent reports whether a and b are logically equivalent: they yield the
// same truth value under every assignment of their combined variables.
func Equivalent(a, b Expr) bool {
	return IsTautology(NewIff(a, b))
}

// Entails reports whether a logically entails b (written a ⊨ b): every
// assignment that satisfies a also satisfies b. Equivalently, a -> b is a
// tautology.
func Entails(a, b Expr) bool {
	return IsTautology(NewImplies(a, b))
}
