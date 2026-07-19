package satsolver

import (
	"sort"
	"strconv"
	"strings"
)

// Lit is a Boolean literal: a signed reference to a variable. A positive value
// v denotes the variable v asserted true; a negative value -v denotes the
// variable v asserted false. Variables are numbered from 1 upward, following
// the DIMACS convention, so the zero value is not a valid literal.
type Lit int

// MakeLit builds the literal for variable v (which must be >= 1). If negated is
// true the returned literal is the negation of v, otherwise it is v asserted
// positively.
func MakeLit(v int, negated bool) Lit {
	if negated {
		return Lit(-v)
	}
	return Lit(v)
}

// PosLit returns the positive literal for variable v.
func PosLit(v int) Lit { return Lit(v) }

// NegLit returns the negative literal for variable v.
func NegLit(v int) Lit { return Lit(-v) }

// Var returns the underlying variable index of the literal, always a
// non-negative number, discarding the sign.
func (l Lit) Var() int {
	if l < 0 {
		return int(-l)
	}
	return int(l)
}

// IsNeg reports whether the literal is negated (its variable asserted false).
func (l Lit) IsNeg() bool { return l < 0 }

// IsPos reports whether the literal is positive (its variable asserted true).
func (l Lit) IsPos() bool { return l > 0 }

// Negate returns the complementary literal, flipping its polarity.
func (l Lit) Negate() Lit { return -l }

// Valid reports whether the literal refers to a real variable (index >= 1).
func (l Lit) Valid() bool { return l != 0 }

// Sign returns +1 for a positive literal and -1 for a negative literal.
func (l Lit) Sign() int {
	if l < 0 {
		return -1
	}
	return 1
}

// Satisfied reports whether the literal evaluates to true under the given
// variable assignment. A variable missing from the assignment is treated as
// false.
func (l Lit) Satisfied(assign map[int]bool) bool {
	val := assign[l.Var()]
	if l.IsNeg() {
		return !val
	}
	return val
}

// String renders the literal in the form "x3" for a positive literal and
// "~x3" for a negative literal.
func (l Lit) String() string {
	if l.IsNeg() {
		return "~x" + strconv.Itoa(l.Var())
	}
	return "x" + strconv.Itoa(l.Var())
}

// Clause is a disjunction of literals (an OR term) when used inside a [CNF] and
// a conjunction of literals (an AND term) when used inside a [DNF]. The slice
// order carries no logical meaning.
type Clause []Lit

// NewClause returns a clause built from the given literals.
func NewClause(lits ...Lit) Clause {
	c := make(Clause, len(lits))
	copy(c, lits)
	return c
}

// Len returns the number of literals in the clause.
func (c Clause) Len() int { return len(c) }

// IsEmpty reports whether the clause has no literals. An empty disjunctive
// clause is unsatisfiable (false); an empty conjunctive term is valid (true).
func (c Clause) IsEmpty() bool { return len(c) == 0 }

// IsUnit reports whether the clause consists of a single literal.
func (c Clause) IsUnit() bool { return len(c) == 1 }

// Contains reports whether the literal l appears in the clause.
func (c Clause) Contains(l Lit) bool {
	for _, x := range c {
		if x == l {
			return true
		}
	}
	return false
}

// ContainsVar reports whether variable v appears in the clause with either
// polarity.
func (c Clause) ContainsVar(v int) bool {
	for _, x := range c {
		if x.Var() == v {
			return true
		}
	}
	return false
}

// Vars returns the sorted set of distinct variable indices occurring in the
// clause.
func (c Clause) Vars() []int {
	seen := map[int]bool{}
	for _, l := range c {
		seen[l.Var()] = true
	}
	out := make([]int, 0, len(seen))
	for v := range seen {
		out = append(out, v)
	}
	sort.Ints(out)
	return out
}

// IsTautology reports whether the disjunctive clause contains a literal and its
// negation, making it trivially true.
func (c Clause) IsTautology() bool {
	seen := map[Lit]bool{}
	for _, l := range c {
		if seen[l.Negate()] {
			return true
		}
		seen[l] = true
	}
	return false
}

// Dedup returns a copy of the clause with duplicate literals removed, order
// preserved.
func (c Clause) Dedup() Clause {
	seen := map[Lit]bool{}
	out := make(Clause, 0, len(c))
	for _, l := range c {
		if !seen[l] {
			seen[l] = true
			out = append(out, l)
		}
	}
	return out
}

// Sorted returns a copy of the clause sorted by variable index then sign.
func (c Clause) Sorted() Clause {
	out := make(Clause, len(c))
	copy(out, c)
	sort.Slice(out, func(i, j int) bool {
		if out[i].Var() != out[j].Var() {
			return out[i].Var() < out[j].Var()
		}
		return out[i] > out[j]
	})
	return out
}

// Negate returns the literal-wise complement of the clause. Interpreting the
// input as a disjunction, the result is the conjunction (as a slice of unit
// clauses) equivalent to its negation by De Morgan's law.
func (c Clause) Negate() Clause {
	out := make(Clause, len(c))
	for i, l := range c {
		out[i] = l.Negate()
	}
	return out
}

// EvalOr evaluates the clause as a disjunction under the assignment: true when
// at least one literal is satisfied.
func (c Clause) EvalOr(assign map[int]bool) bool {
	for _, l := range c {
		if l.Satisfied(assign) {
			return true
		}
	}
	return false
}

// EvalAnd evaluates the clause as a conjunction under the assignment: true when
// every literal is satisfied.
func (c Clause) EvalAnd(assign map[int]bool) bool {
	for _, l := range c {
		if !l.Satisfied(assign) {
			return false
		}
	}
	return true
}

// String renders the clause as a parenthesised disjunction such as
// "(x1 | ~x2 | x3)".
func (c Clause) String() string {
	if len(c) == 0 {
		return "()"
	}
	parts := make([]string, len(c))
	for i, l := range c {
		parts[i] = l.String()
	}
	return "(" + strings.Join(parts, " | ") + ")"
}

// Clone returns an independent copy of the clause.
func (c Clause) Clone() Clause {
	out := make(Clause, len(c))
	copy(out, c)
	return out
}

// MaxVar returns the largest variable index occurring in the clause, or 0 for
// the empty clause.
func (c Clause) MaxVar() int {
	m := 0
	for _, l := range c {
		if v := l.Var(); v > m {
			m = v
		}
	}
	return m
}
