package proofsystems

import (
	"sort"
	"strings"
)

// PLiteral is a propositional literal: a Boolean variable name together with a
// polarity. Neg true denotes the negated literal ¬Var.
type PLiteral struct {
	Var string
	Neg bool
}

// PosPLit returns the positive literal for the given variable.
func PosPLit(v string) PLiteral { return PLiteral{Var: v} }

// NegPLit returns the negative literal for the given variable.
func NegPLit(v string) PLiteral { return PLiteral{Var: v, Neg: true} }

// Negated returns the complementary literal.
func (l PLiteral) Negated() PLiteral { return PLiteral{Var: l.Var, Neg: !l.Neg} }

// String renders the literal as A or !A.
func (l PLiteral) String() string {
	if l.Neg {
		return "!" + l.Var
	}
	return l.Var
}

// Equal reports whether two literals are identical.
func (l PLiteral) Equal(o PLiteral) bool { return l.Var == o.Var && l.Neg == o.Neg }

// IsComplementary reports whether two literals are on the same variable with
// opposite polarity.
func (l PLiteral) IsComplementary(o PLiteral) bool { return l.Var == o.Var && l.Neg != o.Neg }

// Eval returns the truth value of the literal under the assignment.
func (l PLiteral) Eval(a Assignment) bool {
	if l.Neg {
		return !a[l.Var]
	}
	return a[l.Var]
}

// PClause is a propositional clause: a disjunction of literals. The empty
// clause represents falsehood and is the goal of a resolution refutation.
type PClause struct {
	Lits []PLiteral
}

// NewPClause builds a clause from the given literals, removing exact
// duplicates and sorting for a canonical form.
func NewPClause(lits ...PLiteral) PClause {
	return PClause{Lits: dedupeLits(lits)}
}

// IsEmpty reports whether the clause has no literals (the empty clause ⊥).
func (c PClause) IsEmpty() bool { return len(c.Lits) == 0 }

// IsUnit reports whether the clause contains exactly one literal.
func (c PClause) IsUnit() bool { return len(c.Lits) == 1 }

// IsTautologyClause reports whether the clause contains a variable and its
// negation and is therefore always true.
func (c PClause) IsTautologyClause() bool {
	for i := range c.Lits {
		for j := i + 1; j < len(c.Lits); j++ {
			if c.Lits[i].IsComplementary(c.Lits[j]) {
				return true
			}
		}
	}
	return false
}

// Contains reports whether the clause contains the given literal.
func (c PClause) Contains(l PLiteral) bool {
	for _, x := range c.Lits {
		if x.Equal(l) {
			return true
		}
	}
	return false
}

// Vars returns the sorted variable names occurring in the clause.
func (c PClause) Vars() []string {
	set := map[string]bool{}
	for _, l := range c.Lits {
		set[l.Var] = true
	}
	return sortedKeys(set)
}

// Eval returns the truth value of the clause (a disjunction) under the
// assignment.
func (c PClause) Eval(a Assignment) bool {
	for _, l := range c.Lits {
		if l.Eval(a) {
			return true
		}
	}
	return false
}

// Equal reports whether two clauses contain the same set of literals.
func (c PClause) Equal(o PClause) bool {
	if len(c.Lits) != len(o.Lits) {
		return false
	}
	a := dedupeLits(c.Lits)
	b := dedupeLits(o.Lits)
	for i := range a {
		if !a[i].Equal(b[i]) {
			return false
		}
	}
	return true
}

// Subsumes reports whether clause c subsumes o, i.e. every literal of c also
// occurs in o. A subsumed clause is redundant.
func (c PClause) Subsumes(o PClause) bool {
	for _, l := range c.Lits {
		if !o.Contains(l) {
			return false
		}
	}
	return true
}

// String renders the clause as a bracketed disjunction, or [] for the empty
// clause.
func (c PClause) String() string {
	if len(c.Lits) == 0 {
		return "[]"
	}
	parts := make([]string, len(c.Lits))
	for i, l := range c.Lits {
		parts[i] = l.String()
	}
	return "(" + strings.Join(parts, " | ") + ")"
}

// PCNF is a propositional formula in conjunctive normal form: a conjunction of
// clauses.
type PCNF struct {
	Clauses []PClause
}

// NewPCNF builds a CNF from the given clauses.
func NewPCNF(clauses ...PClause) PCNF { return PCNF{Clauses: clauses} }

// Eval returns the truth value of the CNF (a conjunction of clauses) under the
// assignment.
func (n PCNF) Eval(a Assignment) bool {
	for _, c := range n.Clauses {
		if !c.Eval(a) {
			return false
		}
	}
	return true
}

// Vars returns the sorted variable names occurring in the CNF.
func (n PCNF) Vars() []string {
	set := map[string]bool{}
	for _, c := range n.Clauses {
		for _, l := range c.Lits {
			set[l.Var] = true
		}
	}
	return sortedKeys(set)
}

// String renders the CNF as a conjunction of clause strings.
func (n PCNF) String() string {
	if len(n.Clauses) == 0 {
		return "T"
	}
	parts := make([]string, len(n.Clauses))
	for i, c := range n.Clauses {
		parts[i] = c.String()
	}
	return strings.Join(parts, " & ")
}

// Simplify removes tautological clauses and clauses subsumed by another,
// returning a logically equivalent but smaller CNF.
func (n PCNF) Simplify() PCNF {
	var kept []PClause
	for _, c := range n.Clauses {
		if c.IsTautologyClause() {
			continue
		}
		kept = append(kept, c)
	}
	var out []PClause
	for i, c := range kept {
		redundant := false
		for j, d := range kept {
			if i == j {
				continue
			}
			if d.Subsumes(c) && (len(d.Lits) < len(c.Lits) || j < i) {
				redundant = true
				break
			}
		}
		if !redundant {
			out = append(out, c)
		}
	}
	return PCNF{Clauses: out}
}

func dedupeLits(lits []PLiteral) []PLiteral {
	seen := map[string]bool{}
	var out []PLiteral
	for _, l := range lits {
		k := l.String()
		if !seen[k] {
			seen[k] = true
			out = append(out, l)
		}
	}
	sort.Slice(out, func(i, j int) bool {
		if out[i].Var != out[j].Var {
			return out[i].Var < out[j].Var
		}
		return !out[i].Neg && out[j].Neg
	})
	return out
}

// ToNNF converts a propositional formula to negation normal form: negations are
// pushed inward until they apply only to atoms, and implications and
// biconditionals are eliminated. Quantifier-free input is required.
func ToNNF(f Formula) Formula {
	return nnf(f, false)
}

func nnf(f Formula, neg bool) Formula {
	switch f.Conn {
	case ConnAtom:
		if neg {
			return Not(f)
		}
		return f
	case ConnTrue:
		if neg {
			return Bot()
		}
		return Top()
	case ConnFalse:
		if neg {
			return Top()
		}
		return Bot()
	case ConnNot:
		return nnf(f.Subs[0], !neg)
	case ConnAnd:
		if neg {
			return Or(nnf(f.Subs[0], true), nnf(f.Subs[1], true))
		}
		return And(nnf(f.Subs[0], false), nnf(f.Subs[1], false))
	case ConnOr:
		if neg {
			return And(nnf(f.Subs[0], true), nnf(f.Subs[1], true))
		}
		return Or(nnf(f.Subs[0], false), nnf(f.Subs[1], false))
	case ConnImp:
		// a -> b == !a | b
		return nnf(Or(Not(f.Subs[0]), f.Subs[1]), neg)
	case ConnIff:
		// a <-> b == (a -> b) & (b -> a)
		return nnf(And(Imp(f.Subs[0], f.Subs[1]), Imp(f.Subs[1], f.Subs[0])), neg)
	default:
		return f
	}
}

// ToCNFFormula converts a propositional formula into an equivalent formula in
// conjunctive normal form by NNF conversion followed by distribution of
// disjunction over conjunction. The result may be exponentially larger than the
// input; use Tseitin for an equisatisfiable linear-size encoding.
func ToCNFFormula(f Formula) Formula {
	return distributeCNF(ToNNF(f))
}

func distributeCNF(f Formula) Formula {
	switch f.Conn {
	case ConnAnd:
		return And(distributeCNF(f.Subs[0]), distributeCNF(f.Subs[1]))
	case ConnOr:
		l := distributeCNF(f.Subs[0])
		r := distributeCNF(f.Subs[1])
		return distributeOr(l, r)
	default:
		return f
	}
}

func distributeOr(l, r Formula) Formula {
	if l.Conn == ConnAnd {
		return And(distributeOr(l.Subs[0], r), distributeOr(l.Subs[1], r))
	}
	if r.Conn == ConnAnd {
		return And(distributeOr(l, r.Subs[0]), distributeOr(l, r.Subs[1]))
	}
	return Or(l, r)
}

// ToDNFFormula converts a propositional formula into an equivalent formula in
// disjunctive normal form.
func ToDNFFormula(f Formula) Formula {
	return distributeDNF(ToNNF(f))
}

func distributeDNF(f Formula) Formula {
	switch f.Conn {
	case ConnOr:
		return Or(distributeDNF(f.Subs[0]), distributeDNF(f.Subs[1]))
	case ConnAnd:
		l := distributeDNF(f.Subs[0])
		r := distributeDNF(f.Subs[1])
		return distributeAnd(l, r)
	default:
		return f
	}
}

func distributeAnd(l, r Formula) Formula {
	if l.Conn == ConnOr {
		return Or(distributeAnd(l.Subs[0], r), distributeAnd(l.Subs[1], r))
	}
	if r.Conn == ConnOr {
		return Or(distributeAnd(l, r.Subs[0]), distributeAnd(l, r.Subs[1]))
	}
	return And(l, r)
}

// ToCNF converts a propositional formula to the clausal PCNF representation via
// full distribution. Tautological clauses are dropped. The input must be
// quantifier-free.
func ToCNF(f Formula) PCNF {
	cnf := ToCNFFormula(f)
	var clauses []PClause
	for _, conj := range splitConj(cnf) {
		c := clauseOfDisj(conj)
		if !c.IsTautologyClause() {
			clauses = append(clauses, c)
		}
	}
	return PCNF{Clauses: clauses}
}

func splitConj(f Formula) []Formula {
	if f.Conn == ConnAnd {
		return append(splitConj(f.Subs[0]), splitConj(f.Subs[1])...)
	}
	return []Formula{f}
}

func splitDisj(f Formula) []Formula {
	if f.Conn == ConnOr {
		return append(splitDisj(f.Subs[0]), splitDisj(f.Subs[1])...)
	}
	return []Formula{f}
}

func clauseOfDisj(f Formula) PClause {
	var lits []PLiteral
	for _, d := range splitDisj(f) {
		switch d.Conn {
		case ConnAtom:
			lits = append(lits, PosPLit(d.Pred))
		case ConnNot:
			lits = append(lits, NegPLit(d.Subs[0].Pred))
		case ConnTrue:
			// A true disjunct makes the whole clause a tautology; encode it
			// with complementary literals so IsTautologyClause detects it.
			lits = append(lits, PosPLit("__true__"), NegPLit("__true__"))
		case ConnFalse:
			// A false disjunct contributes nothing.
		}
	}
	return NewPClause(lits...)
}

// CNFClauses returns the list of clauses of the clausal form of f.
func CNFClauses(f Formula) []PClause { return ToCNF(f).Clauses }

// IsHornClause reports whether a propositional clause is a Horn clause: it has
// at most one positive literal.
func IsHornClause(c PClause) bool {
	pos := 0
	for _, l := range c.Lits {
		if !l.Neg {
			pos++
		}
	}
	return pos <= 1
}

// IsHornCNF reports whether every clause of the CNF is a Horn clause.
func IsHornCNF(n PCNF) bool {
	for _, c := range n.Clauses {
		if !IsHornClause(c) {
			return false
		}
	}
	return true
}
