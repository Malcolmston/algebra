package proofsystems

import "fmt"

// FOResolve computes the binary resolvents of two first-order clauses. The
// clauses are first renamed apart, then for every pair of complementary
// literals whose argument lists unify, the resolvent is the union of the
// remaining literals with the most-general unifier applied. Each returned
// clause is in canonical form.
func FOResolve(a, b FOClause) []FOClause {
	a = a.RenameApart("_a")
	b = b.RenameApart("_b")
	var out []FOClause
	for i, la := range a.Lits {
		for j, lb := range b.Lits {
			if !la.IsComplementary(lb) {
				continue
			}
			sub, err := UnifyTermLists(la.Args, lb.Args)
			if err != nil {
				continue
			}
			var lits []FOLiteral
			for k, l := range a.Lits {
				if k == i {
					continue
				}
				lits = append(lits, sub.ApplyLiteral(l))
			}
			for k, l := range b.Lits {
				if k == j {
					continue
				}
				lits = append(lits, sub.ApplyLiteral(l))
			}
			out = append(out, FOClause{Lits: lits}.Canonical())
		}
	}
	return out
}

// FOFactor computes the factors of a clause: for each pair of like-signed
// literals that unify, the clause with the unifier applied and the merged
// literal collapsed. Factoring is required for the completeness of binary
// resolution.
func FOFactor(c FOClause) []FOClause {
	var out []FOClause
	for i := 0; i < len(c.Lits); i++ {
		for j := i + 1; j < len(c.Lits); j++ {
			li, lj := c.Lits[i], c.Lits[j]
			if li.Neg != lj.Neg || li.Pred != lj.Pred || len(li.Args) != len(lj.Args) {
				continue
			}
			sub, err := UnifyTermLists(li.Args, lj.Args)
			if err != nil {
				continue
			}
			out = append(out, sub.ApplyClause(c).Canonical())
		}
	}
	return out
}

// FOResolutionRefute attempts to derive the empty clause from a first-order
// clause set by saturation under binary resolution and factoring. Because
// first-order validity is only semi-decidable the search is bounded by
// maxClauses; the function returns true if a refutation is found, and false if
// the bound is reached or the set saturates without the empty clause.
func FOResolutionRefute(clauses []FOClause, maxClauses int) bool {
	work := make([]FOClause, 0, len(clauses))
	seen := map[string]bool{}
	add := func(c FOClause) bool {
		cc := c.Canonical()
		k := cc.String()
		if seen[k] {
			return false
		}
		seen[k] = true
		work = append(work, cc)
		return true
	}
	for _, c := range clauses {
		if c.IsEmpty() {
			return true
		}
		add(c)
		for _, fc := range FOFactor(c) {
			add(fc)
		}
	}
	for len(work) < maxClauses {
		var generated []FOClause
		for i := 0; i < len(work); i++ {
			for j := i; j < len(work); j++ {
				for _, r := range FOResolve(work[i], work[j]) {
					if r.IsEmpty() {
						return true
					}
					generated = append(generated, r)
					generated = append(generated, FOFactor(r)...)
				}
			}
		}
		progress := false
		for _, c := range generated {
			if c.IsEmpty() {
				return true
			}
			if add(c) {
				progress = true
			}
			if len(work) >= maxClauses {
				break
			}
		}
		if !progress {
			return false
		}
	}
	return false
}

// FOEntails reports whether the first-order premises entail the conclusion,
// searched by refuting the clause set of the premises together with the negated
// conclusion under a resolution bound. A true result is sound; a false result
// means no refutation was found within the bound.
func FOEntails(premises []Formula, conclusion Formula, maxClauses int) bool {
	var clauses []FOClause
	for _, p := range premises {
		clauses = append(clauses, Clausify(p)...)
	}
	clauses = append(clauses, Clausify(Not(conclusion))...)
	return FOResolutionRefute(clauses, maxClauses)
}

// FOUnsatisfiable reports whether a first-order clause set is refutable within
// the given bound.
func FOUnsatisfiable(clauses []FOClause, maxClauses int) bool {
	return FOResolutionRefute(clauses, maxClauses)
}

// SkolemName returns the canonical name of the k-th Skolem symbol, matching the
// symbols introduced by Clausify. It is exported to help callers recognise
// Skolem terms in resolved clauses.
func SkolemName(k int) string { return fmt.Sprintf("sk%d", k) }
