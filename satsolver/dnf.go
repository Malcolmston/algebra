package satsolver

import (
	"sort"
	"strings"
)

// DNF is a formula in disjunctive normal form: a disjunction (OR) of terms,
// each of which is a conjunction (AND) of literals. An empty term list is the
// unsatisfiable formula (false); a formula containing an empty term is valid.
type DNF struct {
	// Terms are the disjuncts of the formula.
	Terms []Clause
}

// NewDNF returns a DNF formula built from the given terms.
func NewDNF(terms ...Clause) DNF {
	ts := make([]Clause, len(terms))
	for i, t := range terms {
		ts[i] = t.Clone()
	}
	return DNF{Terms: ts}
}

// AddTerm returns a new formula with t appended to the disjunction.
func (d DNF) AddTerm(t Clause) DNF {
	out := d.Clone()
	out.Terms = append(out.Terms, t.Clone())
	return out
}

// NumTerms returns the number of terms in the formula.
func (d DNF) NumTerms() int { return len(d.Terms) }

// IsEmpty reports whether the formula has no terms, in which case it is
// unsatisfiable.
func (d DNF) IsEmpty() bool { return len(d.Terms) == 0 }

// Clone returns a deep copy of the formula.
func (d DNF) Clone() DNF {
	ts := make([]Clause, len(d.Terms))
	for i, t := range d.Terms {
		ts[i] = t.Clone()
	}
	return DNF{Terms: ts}
}

// MaxVar returns the largest variable index occurring in the formula, or 0.
func (d DNF) MaxVar() int {
	m := 0
	for _, t := range d.Terms {
		if v := t.MaxVar(); v > m {
			m = v
		}
	}
	return m
}

// Vars returns the sorted set of distinct variable indices in the formula.
func (d DNF) Vars() []int {
	seen := map[int]bool{}
	for _, t := range d.Terms {
		for _, l := range t {
			seen[l.Var()] = true
		}
	}
	out := make([]int, 0, len(seen))
	for v := range seen {
		out = append(out, v)
	}
	sort.Ints(out)
	return out
}

// Eval evaluates the formula under the assignment: true when at least one term
// is fully satisfied.
func (d DNF) Eval(assign map[int]bool) bool {
	for _, t := range d.Terms {
		if t.EvalAnd(assign) {
			return true
		}
	}
	return false
}

// Negate returns a CNF formula logically equivalent to the negation of the DNF,
// applying De Morgan's law term-by-term.
func (d DNF) Negate() CNF {
	cs := make([]Clause, len(d.Terms))
	for i, t := range d.Terms {
		cs[i] = t.Negate()
	}
	return CNF{Clauses: cs}
}

// String renders the formula as a disjunction of parenthesised conjunctive
// terms joined by " | ".
func (d DNF) String() string {
	if len(d.Terms) == 0 {
		return "false"
	}
	parts := make([]string, len(d.Terms))
	for i, t := range d.Terms {
		lits := make([]string, len(t))
		for j, l := range t {
			lits[j] = l.String()
		}
		parts[i] = "(" + strings.Join(lits, " & ") + ")"
	}
	return strings.Join(parts, " | ")
}

// Negate on a CNF returns a DNF equivalent to its negation by De Morgan's law.
func (f CNF) Negate() DNF {
	ts := make([]Clause, len(f.Clauses))
	for i, c := range f.Clauses {
		ts[i] = c.Negate()
	}
	return DNF{Terms: ts}
}
