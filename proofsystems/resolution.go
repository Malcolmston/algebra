package proofsystems

// Resolve computes all resolvents of two propositional clauses. For each pair
// of complementary literals (one clause has p, the other ¬p) it produces the
// clause obtained by union of the remaining literals. Tautological resolvents
// are discarded. Most clause pairs share at most one complementary variable, so
// the result is usually a single clause.
func Resolve(a, b PClause) []PClause {
	var out []PClause
	for _, la := range a.Lits {
		for _, lb := range b.Lits {
			if la.IsComplementary(lb) {
				var lits []PLiteral
				for _, x := range a.Lits {
					if !x.Equal(la) {
						lits = append(lits, x)
					}
				}
				for _, y := range b.Lits {
					if !y.Equal(lb) {
						lits = append(lits, y)
					}
				}
				res := NewPClause(lits...)
				if !res.IsTautologyClause() {
					out = append(out, res)
				}
			}
		}
	}
	return out
}

// ResolutionRefute attempts to derive the empty clause from a clause set by
// exhaustive propositional resolution. It returns true when the set is
// unsatisfiable (a refutation exists) and false when saturation is reached
// without producing the empty clause.
func ResolutionRefute(n PCNF) bool {
	clauses := make([]PClause, 0, len(n.Clauses))
	seen := map[string]bool{}
	add := func(c PClause) bool {
		if c.IsTautologyClause() {
			return false
		}
		k := clauseKey(c)
		if seen[k] {
			return false
		}
		seen[k] = true
		clauses = append(clauses, c)
		return true
	}
	for _, c := range n.Clauses {
		add(c)
		if c.IsEmpty() {
			return true
		}
	}
	for {
		var newClauses []PClause
		for i := 0; i < len(clauses); i++ {
			for j := i + 1; j < len(clauses); j++ {
				for _, r := range Resolve(clauses[i], clauses[j]) {
					if r.IsEmpty() {
						return true
					}
					k := clauseKey(r)
					if !seen[k] && !r.IsTautologyClause() {
						newClauses = append(newClauses, r)
					}
				}
			}
		}
		progress := false
		for _, c := range newClauses {
			if add(c) {
				progress = true
			}
		}
		if !progress {
			return false
		}
	}
}

// ResolutionUnsatisfiable reports whether a CNF is unsatisfiable via resolution.
func ResolutionUnsatisfiable(n PCNF) bool { return ResolutionRefute(n) }

// ResolutionValid reports whether a propositional formula is valid by refuting
// its negation with resolution: f is valid iff ¬f is unsatisfiable.
func ResolutionValid(f Formula) bool {
	return ResolutionRefute(ToCNF(Not(f)))
}

// ResolutionEntails reports whether the premises entail the conclusion by
// refuting the clause set of the premises together with the negated conclusion.
func ResolutionEntails(premises []Formula, conclusion Formula) bool {
	var clauses []PClause
	for _, p := range premises {
		clauses = append(clauses, ToCNF(p).Clauses...)
	}
	clauses = append(clauses, ToCNF(Not(conclusion)).Clauses...)
	return ResolutionRefute(PCNF{Clauses: clauses})
}

func clauseKey(c PClause) string {
	return NewPClause(c.Lits...).String()
}
