package satsolver

import (
	"errors"
	"fmt"
	"sort"
	"strconv"
	"strings"
)

// CNF is a formula in conjunctive normal form: a conjunction (AND) of clauses,
// each of which is a disjunction (OR) of literals. An empty clause list is the
// valid formula (true); a formula containing an empty clause is unsatisfiable.
type CNF struct {
	// Clauses are the conjuncts of the formula.
	Clauses []Clause
}

// NewCNF returns a CNF formula built from the given clauses.
func NewCNF(clauses ...Clause) CNF {
	cs := make([]Clause, len(clauses))
	for i, c := range clauses {
		cs[i] = c.Clone()
	}
	return CNF{Clauses: cs}
}

// AddClause returns a new formula with c appended to the conjunction.
func (f CNF) AddClause(c Clause) CNF {
	out := f.Clone()
	out.Clauses = append(out.Clauses, c.Clone())
	return out
}

// NumClauses returns the number of clauses in the formula.
func (f CNF) NumClauses() int { return len(f.Clauses) }

// IsEmpty reports whether the formula has no clauses, in which case it is
// trivially satisfied.
func (f CNF) IsEmpty() bool { return len(f.Clauses) == 0 }

// HasEmptyClause reports whether the formula contains an empty clause, which
// makes it unsatisfiable.
func (f CNF) HasEmptyClause() bool {
	for _, c := range f.Clauses {
		if len(c) == 0 {
			return true
		}
	}
	return false
}

// MaxVar returns the largest variable index occurring in the formula, or 0 when
// there are none.
func (f CNF) MaxVar() int {
	m := 0
	for _, c := range f.Clauses {
		if v := c.MaxVar(); v > m {
			m = v
		}
	}
	return m
}

// Vars returns the sorted set of distinct variable indices in the formula.
func (f CNF) Vars() []int {
	seen := map[int]bool{}
	for _, c := range f.Clauses {
		for _, l := range c {
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

// NumVars returns the count of distinct variables appearing in the formula.
func (f CNF) NumVars() int { return len(f.Vars()) }

// Eval evaluates the formula under the given assignment: true when every clause
// is satisfied.
func (f CNF) Eval(assign map[int]bool) bool {
	for _, c := range f.Clauses {
		if !c.EvalOr(assign) {
			return false
		}
	}
	return true
}

// Clone returns a deep copy of the formula.
func (f CNF) Clone() CNF {
	cs := make([]Clause, len(f.Clauses))
	for i, c := range f.Clauses {
		cs[i] = c.Clone()
	}
	return CNF{Clauses: cs}
}

// Literals returns every literal occurrence across all clauses, with
// repetition.
func (f CNF) Literals() []Lit {
	var out []Lit
	for _, c := range f.Clauses {
		out = append(out, c...)
	}
	return out
}

// CountLiterals returns the total number of literal occurrences in the formula.
func (f CNF) CountLiterals() int {
	n := 0
	for _, c := range f.Clauses {
		n += len(c)
	}
	return n
}

// UnitClauses returns the literals of every unit (single-literal) clause in the
// formula.
func (f CNF) UnitClauses() []Lit {
	var out []Lit
	for _, c := range f.Clauses {
		if len(c) == 1 {
			out = append(out, c[0])
		}
	}
	return out
}

// PureLiterals returns the literals that occur with only one polarity across
// the whole formula. Such literals can always be set true without falsifying a
// clause.
func (f CNF) PureLiterals() []Lit {
	pos := map[int]bool{}
	neg := map[int]bool{}
	for _, c := range f.Clauses {
		for _, l := range c {
			if l.IsNeg() {
				neg[l.Var()] = true
			} else {
				pos[l.Var()] = true
			}
		}
	}
	var out []Lit
	vars := make([]int, 0, len(pos)+len(neg))
	seen := map[int]bool{}
	for v := range pos {
		if !seen[v] {
			seen[v] = true
			vars = append(vars, v)
		}
	}
	for v := range neg {
		if !seen[v] {
			seen[v] = true
			vars = append(vars, v)
		}
	}
	sort.Ints(vars)
	for _, v := range vars {
		switch {
		case pos[v] && !neg[v]:
			out = append(out, PosLit(v))
		case neg[v] && !pos[v]:
			out = append(out, NegLit(v))
		}
	}
	return out
}

// Simplify returns an equivalent formula with tautological clauses dropped and
// duplicate literals removed from each remaining clause.
func (f CNF) Simplify() CNF {
	var cs []Clause
	for _, c := range f.Clauses {
		if c.IsTautology() {
			continue
		}
		cs = append(cs, c.Dedup())
	}
	return CNF{Clauses: cs}
}

// Assign returns the formula conditioned on literal l being true: clauses
// satisfied by l are removed and the complement of l is deleted from the
// remaining clauses (unit propagation of a single literal).
func (f CNF) Assign(l Lit) CNF {
	comp := l.Negate()
	var cs []Clause
	for _, c := range f.Clauses {
		if c.Contains(l) {
			continue
		}
		if c.Contains(comp) {
			nc := make(Clause, 0, len(c))
			for _, x := range c {
				if x != comp {
					nc = append(nc, x)
				}
			}
			cs = append(cs, nc)
		} else {
			cs = append(cs, c.Clone())
		}
	}
	return CNF{Clauses: cs}
}

// String renders the formula as a conjunction of parenthesised clauses joined
// by " & ".
func (f CNF) String() string {
	if len(f.Clauses) == 0 {
		return "true"
	}
	parts := make([]string, len(f.Clauses))
	for i, c := range f.Clauses {
		parts[i] = c.String()
	}
	return strings.Join(parts, " & ")
}

// DIMACS renders the formula in the standard DIMACS CNF text format, with a
// "p cnf <vars> <clauses>" header followed by one zero-terminated clause per
// line.
func (f CNF) DIMACS() string {
	var b strings.Builder
	fmt.Fprintf(&b, "p cnf %d %d\n", f.MaxVar(), len(f.Clauses))
	for _, c := range f.Clauses {
		for _, l := range c {
			b.WriteString(strconv.Itoa(int(l)))
			b.WriteByte(' ')
		}
		b.WriteString("0\n")
	}
	return b.String()
}

// ParseDIMACS parses a formula from DIMACS CNF text. Comment lines beginning
// with 'c' and the optional 'p cnf' header are ignored; literals are read as
// signed integers with 0 terminating each clause.
func ParseDIMACS(s string) (CNF, error) {
	var f CNF
	var cur Clause
	lines := strings.Split(s, "\n")
	for _, line := range lines {
		line = strings.TrimSpace(line)
		if line == "" || strings.HasPrefix(line, "c") || strings.HasPrefix(line, "p") {
			continue
		}
		for _, tok := range strings.Fields(line) {
			n, err := strconv.Atoi(tok)
			if err != nil {
				return CNF{}, fmt.Errorf("satsolver: bad DIMACS token %q: %w", tok, err)
			}
			if n == 0 {
				f.Clauses = append(f.Clauses, cur)
				cur = nil
			} else {
				cur = append(cur, Lit(n))
			}
		}
	}
	if len(cur) > 0 {
		f.Clauses = append(f.Clauses, cur)
	}
	return f, nil
}

// ErrEmptyFormula is returned by routines that require at least one clause.
var ErrEmptyFormula = errors.New("satsolver: empty formula")

// Concat returns the conjunction of two CNF formulas.
func (f CNF) Concat(g CNF) CNF {
	out := f.Clone()
	for _, c := range g.Clauses {
		out.Clauses = append(out.Clauses, c.Clone())
	}
	return out
}

// RemoveDuplicateClauses returns an equivalent formula with syntactically
// duplicate clauses collapsed to a single copy.
func (f CNF) RemoveDuplicateClauses() CNF {
	seen := map[string]bool{}
	var cs []Clause
	for _, c := range f.Clauses {
		key := clauseKey(c)
		if seen[key] {
			continue
		}
		seen[key] = true
		cs = append(cs, c.Clone())
	}
	return CNF{Clauses: cs}
}

func clauseKey(c Clause) string {
	s := c.Sorted()
	parts := make([]string, len(s))
	for i, l := range s {
		parts[i] = strconv.Itoa(int(l))
	}
	return strings.Join(parts, ",")
}
