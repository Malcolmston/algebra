package satsolver

import "sort"

// enumerate calls fn for every assignment over vars, stopping early if fn
// returns false.
func enumerate(vars []string, fn func(env map[string]bool) bool) bool {
	n := len(vars)
	for i := 0; i < (1 << n); i++ {
		if !fn(indexEnv(vars, i)) {
			return false
		}
	}
	return true
}

// IsTautology reports whether e evaluates to true under every assignment.
func IsTautology(e Expr) bool {
	return enumerate(Vars(e), func(env map[string]bool) bool {
		return e.Eval(env)
	})
}

// IsContradiction reports whether e evaluates to false under every assignment.
func IsContradiction(e Expr) bool {
	return enumerate(Vars(e), func(env map[string]bool) bool {
		return !e.Eval(env)
	})
}

// IsSatisfiable reports whether e evaluates to true under at least one
// assignment.
func IsSatisfiable(e Expr) bool { return !IsContradiction(e) }

// IsContingent reports whether e is neither a tautology nor a contradiction:
// true for some assignments and false for others.
func IsContingent(e Expr) bool {
	return !IsTautology(e) && !IsContradiction(e)
}

// Equivalent reports whether two expressions have identical truth values on
// every assignment to their combined variables.
func Equivalent(a, b Expr) bool {
	vars := unionVars(a, b)
	return enumerate(vars, func(env map[string]bool) bool {
		return a.Eval(env) == b.Eval(env)
	})
}

// Entails reports whether a semantically entails b (a |= b): every assignment
// satisfying a also satisfies b.
func Entails(a, b Expr) bool {
	vars := unionVars(a, b)
	return enumerate(vars, func(env map[string]bool) bool {
		return !a.Eval(env) || b.Eval(env)
	})
}

// Implies reports the same relation as [Entails]; provided as a readable alias
// for semantic implication.
func Implies2(a, b Expr) bool { return Entails(a, b) }

// FindModel returns a satisfying assignment for e and true, or nil and false if
// e is unsatisfiable.
func FindModel(e Expr) (map[string]bool, bool) {
	vars := Vars(e)
	var found map[string]bool
	ok := !enumerate(vars, func(env map[string]bool) bool {
		if e.Eval(env) {
			found = copyEnv(env)
			return false
		}
		return true
	})
	return found, ok
}

// AllModels returns every satisfying assignment of e over its variables.
func AllModels(e Expr) []map[string]bool {
	vars := Vars(e)
	var out []map[string]bool
	enumerate(vars, func(env map[string]bool) bool {
		if e.Eval(env) {
			out = append(out, copyEnv(env))
		}
		return true
	})
	return out
}

// CountModels returns the number of satisfying assignments of e over its
// variables.
func CountModels(e Expr) int {
	vars := Vars(e)
	count := 0
	enumerate(vars, func(env map[string]bool) bool {
		if e.Eval(env) {
			count++
		}
		return true
	})
	return count
}

// unionVars returns the sorted union of the variable names of a and b.
func unionVars(a, b Expr) []string {
	seen := map[string]bool{}
	for _, v := range Vars(a) {
		seen[v] = true
	}
	for _, v := range Vars(b) {
		seen[v] = true
	}
	out := make([]string, 0, len(seen))
	for v := range seen {
		out = append(out, v)
	}
	sort.Strings(out)
	return out
}

func copyEnv(env map[string]bool) map[string]bool {
	out := make(map[string]bool, len(env))
	for k, v := range env {
		out[k] = v
	}
	return out
}

// EquivalentCNF reports whether two clausal CNF formulas are logically
// equivalent by comparing their truth tables over the union of their variables.
func EquivalentCNF(f, g CNF) bool {
	seen := map[int]bool{}
	for _, v := range f.Vars() {
		seen[v] = true
	}
	for _, v := range g.Vars() {
		seen[v] = true
	}
	vars := make([]int, 0, len(seen))
	for v := range seen {
		vars = append(vars, v)
	}
	sort.Ints(vars)
	n := len(vars)
	for mask := 0; mask < (1 << n); mask++ {
		a := maskToAssign(vars, mask)
		if f.Eval(a) != g.Eval(a) {
			return false
		}
	}
	return true
}
