package proofsystems

// IsTautology reports whether the propositional formula is true under every
// assignment of its atoms. A quantifier-free variable-free formula is a
// tautology exactly when it evaluates to true.
func IsTautology(f Formula) bool {
	for _, a := range AllAssignments(f.PropVars()) {
		v, err := Eval(f, a)
		if err != nil || !v {
			return false
		}
	}
	return true
}

// IsValid is a synonym for IsTautology.
func IsValid(f Formula) bool { return IsTautology(f) }

// IsContradiction reports whether the formula is false under every assignment.
func IsContradiction(f Formula) bool {
	for _, a := range AllAssignments(f.PropVars()) {
		v, err := Eval(f, a)
		if err != nil || v {
			return false
		}
	}
	return true
}

// IsUnsatisfiable is a synonym for IsContradiction.
func IsUnsatisfiable(f Formula) bool { return IsContradiction(f) }

// IsSatisfiable reports whether some assignment makes the formula true.
func IsSatisfiable(f Formula) bool {
	for _, a := range AllAssignments(f.PropVars()) {
		if v, err := Eval(f, a); err == nil && v {
			return true
		}
	}
	return false
}

// IsContingent reports whether the formula is satisfiable but not a tautology,
// i.e. true under some assignments and false under others.
func IsContingent(f Formula) bool {
	return IsSatisfiable(f) && !IsTautology(f)
}

// FindModel returns a satisfying assignment for the formula and true, or the
// zero assignment and false when the formula is unsatisfiable.
func FindModel(f Formula) (Assignment, bool) {
	for _, a := range AllAssignments(f.PropVars()) {
		if v, err := Eval(f, a); err == nil && v {
			return a, true
		}
	}
	return nil, false
}

// AllModels returns every satisfying assignment of the formula over its atoms.
func AllModels(f Formula) []Assignment {
	var out []Assignment
	for _, a := range AllAssignments(f.PropVars()) {
		if v, err := Eval(f, a); err == nil && v {
			out = append(out, a)
		}
	}
	return out
}

// CountModels returns the number of satisfying assignments of the formula.
func CountModels(f Formula) int {
	n := 0
	for _, a := range AllAssignments(f.PropVars()) {
		if v, err := Eval(f, a); err == nil && v {
			n++
		}
	}
	return n
}

// Models reports whether the assignment satisfies the formula. Atoms of the
// formula not present in the assignment default to false.
func Models(a Assignment, f Formula) bool {
	v, err := Eval(f, a)
	return err == nil && v
}

// Equivalent reports whether two propositional formulas have the same truth
// value under every assignment of their combined atoms.
func Equivalent(f, g Formula) bool {
	atoms := unionAtoms(f.PropVars(), g.PropVars())
	for _, a := range AllAssignments(atoms) {
		vf, ef := Eval(f, a)
		vg, eg := Eval(g, a)
		if ef != nil || eg != nil || vf != vg {
			return false
		}
	}
	return true
}

// Entails reports whether the set of premises semantically entails the
// conclusion: every assignment satisfying all premises also satisfies the
// conclusion.
func Entails(premises []Formula, conclusion Formula) bool {
	names := conclusion.PropVars()
	for _, p := range premises {
		names = unionAtoms(names, p.PropVars())
	}
	for _, a := range AllAssignments(names) {
		if allTrue(premises, a) {
			if v, err := Eval(conclusion, a); err != nil || !v {
				return false
			}
		}
	}
	return true
}

// EntailsSingle reports whether a single premise entails the conclusion.
func EntailsSingle(premise, conclusion Formula) bool {
	return Entails([]Formula{premise}, conclusion)
}

// Consistent reports whether a set of formulas is jointly satisfiable.
func Consistent(fs []Formula) bool {
	var names []string
	for _, f := range fs {
		names = unionAtoms(names, f.PropVars())
	}
	for _, a := range AllAssignments(names) {
		if allTrue(fs, a) {
			return true
		}
	}
	return false
}

// Inconsistent reports whether a set of formulas is jointly unsatisfiable.
func Inconsistent(fs []Formula) bool { return !Consistent(fs) }

func allTrue(fs []Formula, a Assignment) bool {
	for _, f := range fs {
		v, err := Eval(f, a)
		if err != nil || !v {
			return false
		}
	}
	return true
}

func unionAtoms(a, b []string) []string {
	set := map[string]bool{}
	for _, x := range a {
		set[x] = true
	}
	for _, x := range b {
		set[x] = true
	}
	return sortedKeys(set)
}
