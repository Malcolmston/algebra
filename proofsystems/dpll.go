package proofsystems

import "sort"

// DPLLResult reports the outcome of a DPLL search: whether the clause set is
// satisfiable and, when it is, a satisfying assignment.
type DPLLResult struct {
	Sat   bool
	Model Assignment
}

// DPLL decides the satisfiability of a propositional clause set using the
// Davis–Putnam–Logemann–Loveland procedure with unit propagation and
// pure-literal elimination followed by splitting on the first unassigned
// variable. The returned model, when satisfiable, assigns every variable of the
// clause set.
func DPLL(n PCNF) DPLLResult {
	vars := n.Vars()
	model, ok := dpllSolve(n.Clauses, Assignment{})
	if ok {
		// Fill in any variable left free (from pure elimination) with false.
		for _, v := range vars {
			if _, present := model[v]; !present {
				model[v] = false
			}
		}
		return DPLLResult{Sat: true, Model: model}
	}
	return DPLLResult{Sat: false}
}

// DPLLSatisfiable reports whether the clause set is satisfiable.
func DPLLSatisfiable(n PCNF) bool { return DPLL(n).Sat }

func dpllSolve(clauses []PClause, assign Assignment) (Assignment, bool) {
	clauses, assign, ok := propagate(clauses, assign)
	if !ok {
		return nil, false
	}
	if len(clauses) == 0 {
		return assign, true
	}
	for _, c := range clauses {
		if c.IsEmpty() {
			return nil, false
		}
	}
	// Pure-literal elimination.
	if lit, ok := findPureLiteral(clauses); ok {
		na := assign.Clone()
		na[lit.Var] = !lit.Neg
		return dpllSolve(simplify(clauses, lit), na)
	}
	// Choose a branching variable and split.
	v := branchVar(clauses)
	for _, val := range []bool{true, false} {
		na := assign.Clone()
		na[v] = val
		lit := PLiteral{Var: v, Neg: !val}
		if model, ok := dpllSolve(simplify(clauses, lit), na); ok {
			return model, true
		}
	}
	return nil, false
}

// propagate performs unit propagation to a fixpoint, extending assign and
// simplifying clauses. It returns ok=false on a conflict (an empty clause).
func propagate(clauses []PClause, assign Assignment) ([]PClause, Assignment, bool) {
	assign = assign.Clone()
	for {
		var unit *PLiteral
		for _, c := range clauses {
			if c.IsEmpty() {
				return clauses, assign, false
			}
			if c.IsUnit() {
				l := c.Lits[0]
				unit = &l
				break
			}
		}
		if unit == nil {
			return clauses, assign, true
		}
		assign[unit.Var] = !unit.Neg
		clauses = simplify(clauses, *unit)
		for _, c := range clauses {
			if c.IsEmpty() {
				return clauses, assign, false
			}
		}
	}
}

// simplify returns the clause set conditioned on the literal being true:
// clauses containing the literal are removed and the complementary literal is
// deleted from the remaining clauses.
func simplify(clauses []PClause, lit PLiteral) []PClause {
	comp := lit.Negated()
	var out []PClause
	for _, c := range clauses {
		if c.Contains(lit) {
			continue
		}
		var lits []PLiteral
		for _, l := range c.Lits {
			if l.Equal(comp) {
				continue
			}
			lits = append(lits, l)
		}
		out = append(out, PClause{Lits: lits})
	}
	return out
}

// findPureLiteral returns a literal whose variable appears with only one
// polarity across all clauses, if any.
func findPureLiteral(clauses []PClause) (PLiteral, bool) {
	pos := map[string]bool{}
	neg := map[string]bool{}
	for _, c := range clauses {
		for _, l := range c.Lits {
			if l.Neg {
				neg[l.Var] = true
			} else {
				pos[l.Var] = true
			}
		}
	}
	vars := map[string]bool{}
	for v := range pos {
		vars[v] = true
	}
	for v := range neg {
		vars[v] = true
	}
	for _, v := range sortedKeys(vars) {
		if pos[v] && !neg[v] {
			return PosPLit(v), true
		}
		if neg[v] && !pos[v] {
			return NegPLit(v), true
		}
	}
	return PLiteral{}, false
}

func branchVar(clauses []PClause) string {
	set := map[string]bool{}
	for _, c := range clauses {
		for _, l := range c.Lits {
			set[l.Var] = true
		}
	}
	keys := sortedKeys(set)
	sort.Strings(keys)
	return keys[0]
}

// UnitPropagate applies a single fixpoint of unit propagation to a clause set,
// returning the simplified clauses, the forced literal assignment, and whether
// a conflict (an empty clause) was reached.
func UnitPropagate(n PCNF) (PCNF, Assignment, bool) {
	clauses, assign, ok := propagate(n.Clauses, Assignment{})
	return PCNF{Clauses: clauses}, assign, ok
}

// PureLiteralElimination repeatedly removes clauses satisfied by pure literals,
// returning the reduced clause set and the assignment fixing the pure literals.
func PureLiteralElimination(n PCNF) (PCNF, Assignment) {
	clauses := n.Clauses
	assign := Assignment{}
	for {
		lit, ok := findPureLiteral(clauses)
		if !ok {
			break
		}
		assign[lit.Var] = !lit.Neg
		clauses = simplify(clauses, lit)
	}
	return PCNF{Clauses: clauses}, assign
}

// DPLLSatisfiableFormula converts a propositional formula to CNF via the
// Tseitin encoding and decides its satisfiability with DPLL. This avoids the
// exponential blow-up of naive CNF conversion.
func DPLLSatisfiableFormula(f Formula) bool {
	return DPLLSatisfiable(TseitinCNF(f, "_t"))
}
