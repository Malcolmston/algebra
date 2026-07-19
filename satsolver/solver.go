package satsolver

import "sort"

// copyAssign returns an independent copy of a variable assignment.
func copyAssign(a map[int]bool) map[int]bool {
	out := make(map[int]bool, len(a))
	for k, v := range a {
		out[k] = v
	}
	return out
}

// UnitPropagate repeatedly assigns the literal of every unit clause and
// simplifies the formula accordingly. It returns the simplified formula, the
// extended assignment, and a boolean that is false when propagation derived a
// conflict (an empty clause).
func UnitPropagate(f CNF, assign map[int]bool) (CNF, map[int]bool, bool) {
	a := copyAssign(assign)
	for {
		if f.HasEmptyClause() {
			return f, a, false
		}
		units := f.UnitClauses()
		if len(units) == 0 {
			return f, a, true
		}
		l := units[0]
		a[l.Var()] = l.IsPos()
		f = f.Assign(l)
	}
}

// PureLiteralElimination assigns every pure literal (one occurring with a single
// polarity) to the value that satisfies its clauses, repeating until no pure
// literals remain. It returns the simplified formula and the extended
// assignment.
func PureLiteralElimination(f CNF, assign map[int]bool) (CNF, map[int]bool) {
	a := copyAssign(assign)
	for {
		pures := f.PureLiterals()
		if len(pures) == 0 {
			return f, a
		}
		for _, l := range pures {
			a[l.Var()] = l.IsPos()
			f = f.Assign(l)
		}
	}
}

// chooseLiteral selects a branching literal from the formula: the first literal
// of the shortest clause.
func chooseLiteral(f CNF) (Lit, bool) {
	best := -1
	var lit Lit
	for _, c := range f.Clauses {
		if len(c) == 0 {
			continue
		}
		if best == -1 || len(c) < best {
			best = len(c)
			lit = c[0]
		}
	}
	if best == -1 {
		return 0, false
	}
	return lit, true
}

// DPLL solves the satisfiability of a CNF formula with the classic
// Davis-Putnam-Logemann-Loveland procedure: unit propagation, pure-literal
// elimination and backtracking search. It returns a satisfying partial
// assignment and true, or nil and false when the formula is unsatisfiable.
func DPLL(f CNF) (map[int]bool, bool) {
	return dpll(f.Simplify(), map[int]bool{})
}

func dpll(f CNF, assign map[int]bool) (map[int]bool, bool) {
	f, assign, ok := UnitPropagate(f, assign)
	if !ok {
		return nil, false
	}
	f, assign = PureLiteralElimination(f, assign)
	if f.IsEmpty() {
		return assign, true
	}
	if f.HasEmptyClause() {
		return nil, false
	}
	l, ok := chooseLiteral(f)
	if !ok {
		return assign, true
	}
	// Branch on l = true.
	a1 := copyAssign(assign)
	a1[l.Var()] = l.IsPos()
	if res, ok := dpll(f.Assign(l), a1); ok {
		return res, true
	}
	// Branch on l = false.
	a2 := copyAssign(assign)
	a2[l.Var()] = !l.IsPos()
	return dpll(f.Assign(l.Negate()), a2)
}

// SolveCNF solves f and, on success, returns a complete assignment covering
// every variable of f (variables left free by the search are set to false).
func SolveCNF(f CNF) (map[int]bool, bool) {
	partial, ok := DPLL(f)
	if !ok {
		return nil, false
	}
	full := copyAssign(partial)
	for _, v := range f.Vars() {
		if _, seen := full[v]; !seen {
			full[v] = false
		}
	}
	return full, true
}

// IsSatisfiableCNF reports whether the CNF formula has at least one satisfying
// assignment.
func IsSatisfiableCNF(f CNF) bool {
	_, ok := DPLL(f)
	return ok
}

// IsUnsatisfiableCNF reports whether the CNF formula has no satisfying
// assignment.
func IsUnsatisfiableCNF(f CNF) bool { return !IsSatisfiableCNF(f) }

// IsTautologyCNF reports whether every assignment satisfies f, i.e. the
// negation of f is unsatisfiable. It is decided by enumeration over the
// variables of f.
func IsTautologyCNF(f CNF) bool {
	vars := f.Vars()
	n := len(vars)
	for mask := 0; mask < (1 << n); mask++ {
		assign := maskToAssign(vars, mask)
		if !f.Eval(assign) {
			return false
		}
	}
	return true
}

// AllSolutions returns every complete satisfying assignment of f over its
// variables, enumerated in ascending mask order.
func AllSolutions(f CNF) []map[int]bool {
	vars := f.Vars()
	n := len(vars)
	var out []map[int]bool
	for mask := 0; mask < (1 << n); mask++ {
		assign := maskToAssign(vars, mask)
		if f.Eval(assign) {
			out = append(out, assign)
		}
	}
	return out
}

// CountSolutions returns the number of complete satisfying assignments of f
// over its variables.
func CountSolutions(f CNF) int {
	vars := f.Vars()
	n := len(vars)
	count := 0
	for mask := 0; mask < (1 << n); mask++ {
		if f.Eval(maskToAssign(vars, mask)) {
			count++
		}
	}
	return count
}

func maskToAssign(vars []int, mask int) map[int]bool {
	assign := make(map[int]bool, len(vars))
	for i, v := range vars {
		assign[v] = (mask>>i)&1 == 1
	}
	return assign
}

// Resolve computes the resolvent of two clauses on the pivot variable v: the
// disjunction of all their other literals, provided one clause contains the
// positive and the other the negative literal of v. It returns the resolvent
// and true, or nil and false when the clauses cannot be resolved on v.
func Resolve(c1, c2 Clause, v int) (Clause, bool) {
	has1pos := c1.Contains(PosLit(v))
	has1neg := c1.Contains(NegLit(v))
	has2pos := c2.Contains(PosLit(v))
	has2neg := c2.Contains(NegLit(v))
	if !((has1pos && has2neg) || (has1neg && has2pos)) {
		return nil, false
	}
	seen := map[Lit]bool{}
	var out Clause
	add := func(c Clause) {
		for _, l := range c {
			if l.Var() == v {
				continue
			}
			if !seen[l] {
				seen[l] = true
				out = append(out, l)
			}
		}
	}
	add(c1)
	add(c2)
	// A resolvent containing a literal and its complement is a tautology.
	for _, l := range out {
		if seen[l.Negate()] {
			// still return it; caller may discard tautologies
			break
		}
	}
	sort.Slice(out, func(i, j int) bool {
		if out[i].Var() != out[j].Var() {
			return out[i].Var() < out[j].Var()
		}
		return out[i] > out[j]
	})
	return out, true
}

// Solver is a small stateful wrapper around the DPLL search that remembers the
// formula and its most recent model.
type Solver struct {
	Formula CNF
	model   map[int]bool
	solved  bool
	sat     bool
}

// NewSolver returns a solver for the given formula.
func NewSolver(f CNF) *Solver { return &Solver{Formula: f.Clone()} }

// AddClause adds a clause to the solver's formula and invalidates any cached
// result.
func (s *Solver) AddClause(c Clause) {
	s.Formula = s.Formula.AddClause(c)
	s.solved = false
}

// Solve runs the search and returns whether the formula is satisfiable, caching
// the outcome and model.
func (s *Solver) Solve() bool {
	if s.solved {
		return s.sat
	}
	m, ok := SolveCNF(s.Formula)
	s.model = m
	s.sat = ok
	s.solved = true
	return ok
}

// Model returns the satisfying assignment found by the most recent successful
// [Solver.Solve]. It is nil when the formula is unsatisfiable or unsolved.
func (s *Solver) Model() map[int]bool {
	if !s.solved {
		s.Solve()
	}
	return s.model
}
