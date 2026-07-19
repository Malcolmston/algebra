package proofsystems

// TableauClosed reports whether the analytic (semantic) tableau for the given
// set of propositional formulas closes on every branch, which holds exactly
// when the set is unsatisfiable. The method decomposes formulas with the
// standard alpha (non-branching) and beta (branching) rules until every branch
// either contains a formula together with its negation or is saturated with
// literals.
func TableauClosed(formulas []Formula) bool {
	branch := make([]Formula, len(formulas))
	copy(branch, formulas)
	return closes(branch)
}

// TableauSatisfiable reports whether a set of propositional formulas is
// satisfiable according to the tableau method (some branch stays open).
func TableauSatisfiable(formulas []Formula) bool {
	return !TableauClosed(formulas)
}

// TableauValid reports whether a propositional formula is valid by checking
// that the tableau for its negation closes.
func TableauValid(f Formula) bool {
	return TableauClosed([]Formula{Not(f)})
}

// TableauProve is a synonym for TableauValid: it reports whether f is a
// propositional tautology.
func TableauProve(f Formula) bool { return TableauValid(f) }

// TableauEntails reports whether the premises entail the conclusion by closing
// the tableau for the premises together with the negated conclusion.
func TableauEntails(premises []Formula, conclusion Formula) bool {
	set := append([]Formula{}, premises...)
	set = append(set, Not(conclusion))
	return TableauClosed(set)
}

func closes(branch []Formula) bool {
	if hasContradiction(branch) {
		return true
	}
	// Find a non-literal formula to expand; prefer alpha rules.
	idx := -1
	for i, f := range branch {
		if !isTableauLiteral(f) {
			idx = i
			break
		}
	}
	if idx == -1 {
		// Saturated branch with only literals and no contradiction: open.
		return false
	}
	f := branch[idx]
	rest := make([]Formula, 0, len(branch)-1)
	rest = append(rest, branch[:idx]...)
	rest = append(rest, branch[idx+1:]...)

	alpha, beta := tableauRule(f)
	if alpha != nil {
		next := append(rest, alpha...)
		return closes(next)
	}
	// Beta rule: branch into each option; the node closes only if all do.
	for _, opt := range beta {
		next := append(append([]Formula{}, rest...), opt...)
		if !closes(next) {
			return false
		}
	}
	return true
}

// tableauRule returns either the alpha (conjunctive) expansion of f or, when
// alpha is nil, its beta (disjunctive) expansion as a list of alternative
// branches.
func tableauRule(f Formula) (alpha []Formula, beta [][]Formula) {
	switch f.Conn {
	case ConnNot:
		g := f.Subs[0]
		switch g.Conn {
		case ConnNot:
			return []Formula{g.Subs[0]}, nil
		case ConnAnd:
			return nil, [][]Formula{{Not(g.Subs[0])}, {Not(g.Subs[1])}}
		case ConnOr:
			return []Formula{Not(g.Subs[0]), Not(g.Subs[1])}, nil
		case ConnImp:
			return []Formula{g.Subs[0], Not(g.Subs[1])}, nil
		case ConnIff:
			return nil, [][]Formula{
				{g.Subs[0], Not(g.Subs[1])},
				{Not(g.Subs[0]), g.Subs[1]},
			}
		case ConnTrue:
			return []Formula{Bot()}, nil
		case ConnFalse:
			return []Formula{Top()}, nil
		}
	case ConnAnd:
		return []Formula{f.Subs[0], f.Subs[1]}, nil
	case ConnOr:
		return nil, [][]Formula{{f.Subs[0]}, {f.Subs[1]}}
	case ConnImp:
		return nil, [][]Formula{{Not(f.Subs[0])}, {f.Subs[1]}}
	case ConnIff:
		return nil, [][]Formula{
			{f.Subs[0], f.Subs[1]},
			{Not(f.Subs[0]), Not(f.Subs[1])},
		}
	}
	return nil, nil
}

func isTableauLiteral(f Formula) bool {
	switch f.Conn {
	case ConnAtom, ConnTrue, ConnFalse:
		return true
	case ConnNot:
		g := f.Subs[0]
		return g.Conn == ConnAtom
	default:
		return false
	}
}

func hasContradiction(branch []Formula) bool {
	for _, f := range branch {
		if f.Conn == ConnFalse {
			return true
		}
		if f.Conn == ConnNot && f.Subs[0].Conn == ConnTrue {
			return true
		}
	}
	for i := range branch {
		for j := range branch {
			if i == j {
				continue
			}
			if branch[i].Conn == ConnNot && branch[i].Subs[0].Equal(branch[j]) {
				return true
			}
		}
	}
	return false
}
