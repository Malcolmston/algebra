package proofsystems

import (
	"fmt"
	"sort"
	"strings"
)

// FOLiteral is a first-order literal: a possibly negated atomic predicate
// applied to argument terms.
type FOLiteral struct {
	Neg  bool
	Pred string
	Args []Term
}

// NewFOLiteral builds a positive first-order literal P(args...).
func NewFOLiteral(pred string, args ...Term) FOLiteral {
	cp := make([]Term, len(args))
	copy(cp, args)
	return FOLiteral{Pred: pred, Args: cp}
}

// NewNegFOLiteral builds a negative first-order literal ¬P(args...).
func NewNegFOLiteral(pred string, args ...Term) FOLiteral {
	l := NewFOLiteral(pred, args...)
	l.Neg = true
	return l
}

// Negated returns the complementary literal.
func (l FOLiteral) Negated() FOLiteral {
	return FOLiteral{Neg: !l.Neg, Pred: l.Pred, Args: l.Args}
}

// Equal reports whether two first-order literals are structurally identical.
func (l FOLiteral) Equal(o FOLiteral) bool {
	return l.Neg == o.Neg && l.Pred == o.Pred && TermsEqual(l.Args, o.Args)
}

// IsComplementary reports whether two literals share a predicate and arity but
// differ in sign, so they are candidates for resolution once their arguments
// unify.
func (l FOLiteral) IsComplementary(o FOLiteral) bool {
	return l.Neg != o.Neg && l.Pred == o.Pred && len(l.Args) == len(o.Args)
}

// Vars returns the sorted variable names occurring in the literal.
func (l FOLiteral) Vars() []string { return CollectVars(l.Args...) }

// String renders the literal in P(a,b) or !P(a,b) form.
func (l FOLiteral) String() string {
	var b strings.Builder
	if l.Neg {
		b.WriteString("!")
	}
	b.WriteString(l.Pred)
	if len(l.Args) > 0 {
		parts := make([]string, len(l.Args))
		for i, a := range l.Args {
			parts[i] = a.String()
		}
		b.WriteString("(" + strings.Join(parts, ",") + ")")
	}
	return b.String()
}

// FOClause is a first-order clause: a disjunction of first-order literals with
// every variable implicitly universally quantified. The empty clause denotes
// falsehood.
type FOClause struct {
	Lits []FOLiteral
}

// NewFOClause builds a clause from the given literals.
func NewFOClause(lits ...FOLiteral) FOClause {
	cp := make([]FOLiteral, len(lits))
	copy(cp, lits)
	return FOClause{Lits: cp}
}

// IsEmpty reports whether the clause has no literals.
func (c FOClause) IsEmpty() bool { return len(c.Lits) == 0 }

// Vars returns the sorted variable names occurring in the clause.
func (c FOClause) Vars() []string {
	set := map[string]bool{}
	for _, l := range c.Lits {
		for _, v := range l.Vars() {
			set[v] = true
		}
	}
	return sortedKeys(set)
}

// String renders the clause as a bracketed disjunction, or [] for the empty
// clause.
func (c FOClause) String() string {
	if len(c.Lits) == 0 {
		return "[]"
	}
	parts := make([]string, len(c.Lits))
	for i, l := range c.Lits {
		parts[i] = l.String()
	}
	return "{" + strings.Join(parts, " | ") + "}"
}

// RenameApart returns a copy of the clause with every variable renamed using
// the given suffix, so that two clauses can be resolved without variable
// capture.
func (c FOClause) RenameApart(suffix string) FOClause {
	r := map[string]string{}
	for _, v := range c.Vars() {
		r[v] = v + suffix
	}
	lits := make([]FOLiteral, len(c.Lits))
	for i, l := range c.Lits {
		args := make([]Term, len(l.Args))
		for j, a := range l.Args {
			args[j] = a.Rename(r)
		}
		lits[i] = FOLiteral{Neg: l.Neg, Pred: l.Pred, Args: args}
	}
	return FOClause{Lits: lits}
}

// Canonical returns a clause with literals de-duplicated and sorted by string
// form, giving a stable representative for set membership.
func (c FOClause) Canonical() FOClause {
	seen := map[string]bool{}
	var lits []FOLiteral
	for _, l := range c.Lits {
		k := l.String()
		if !seen[k] {
			seen[k] = true
			lits = append(lits, l)
		}
	}
	sort.Slice(lits, func(i, j int) bool { return lits[i].String() < lits[j].String() })
	return FOClause{Lits: lits}
}

// literalToFO converts a propositional/first-order literal formula (an atom or
// negated atom) into an FOLiteral.
func literalToFO(f Formula) (FOLiteral, bool) {
	switch f.Conn {
	case ConnAtom:
		return FOLiteral{Pred: f.Pred, Args: f.Args}, true
	case ConnNot:
		if f.Subs[0].Conn == ConnAtom {
			g := f.Subs[0]
			return FOLiteral{Neg: true, Pred: g.Pred, Args: g.Args}, true
		}
	}
	return FOLiteral{}, false
}

// FONNF converts a first-order formula to negation normal form, eliminating
// implications and biconditionals and pushing negations through the
// quantifiers (¬∀ becomes ∃¬ and ¬∃ becomes ∀¬).
func FONNF(f Formula) Formula { return fonnf(f, false) }

func fonnf(f Formula, neg bool) Formula {
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
		return fonnf(f.Subs[0], !neg)
	case ConnAnd:
		if neg {
			return Or(fonnf(f.Subs[0], true), fonnf(f.Subs[1], true))
		}
		return And(fonnf(f.Subs[0], false), fonnf(f.Subs[1], false))
	case ConnOr:
		if neg {
			return And(fonnf(f.Subs[0], true), fonnf(f.Subs[1], true))
		}
		return Or(fonnf(f.Subs[0], false), fonnf(f.Subs[1], false))
	case ConnImp:
		return fonnf(Or(Not(f.Subs[0]), f.Subs[1]), neg)
	case ConnIff:
		return fonnf(And(Imp(f.Subs[0], f.Subs[1]), Imp(f.Subs[1], f.Subs[0])), neg)
	case ConnForall:
		if neg {
			return Exists(f.Bound, fonnf(f.Subs[0], true))
		}
		return Forall(f.Bound, fonnf(f.Subs[0], false))
	case ConnExists:
		if neg {
			return Forall(f.Bound, fonnf(f.Subs[0], true))
		}
		return Exists(f.Bound, fonnf(f.Subs[0], false))
	default:
		return f
	}
}

type clausifier struct {
	varCtr int
	skCtr  int
}

func (c *clausifier) freshVar() string {
	c.varCtr++
	return fmt.Sprintf("_v%d", c.varCtr)
}

func (c *clausifier) freshSkolem() string {
	c.skCtr++
	return fmt.Sprintf("sk%d", c.skCtr)
}

// standardize renames every bound variable to a fresh name to avoid clashes.
func (c *clausifier) standardize(f Formula, ren map[string]string) Formula {
	switch f.Conn {
	case ConnAtom:
		sub := NewSubstitution()
		for old, nw := range ren {
			sub = sub.Bind(old, NewVar(nw))
		}
		return Formula{Conn: ConnAtom, Pred: f.Pred, Args: sub.ApplyTerms(f.Args)}
	case ConnForall, ConnExists:
		nv := c.freshVar()
		nr := map[string]string{}
		for k, v := range ren {
			nr[k] = v
		}
		nr[f.Bound] = nv
		return Formula{Conn: f.Conn, Bound: nv, Subs: []Formula{c.standardize(f.Subs[0], nr)}}
	default:
		subs := make([]Formula, len(f.Subs))
		for i, s := range f.Subs {
			subs[i] = c.standardize(s, ren)
		}
		return Formula{Conn: f.Conn, Subs: subs}
	}
}

// skolemize removes existential quantifiers, replacing each existential
// variable by a Skolem term over the universally quantified variables in scope,
// and drops universal quantifiers.
func (c *clausifier) skolemize(f Formula, univ []string) Formula {
	switch f.Conn {
	case ConnForall:
		return c.skolemize(f.Subs[0], append(append([]string{}, univ...), f.Bound))
	case ConnExists:
		var sk Term
		if len(univ) == 0 {
			sk = NewConst(c.freshSkolem())
		} else {
			args := make([]Term, len(univ))
			for i, v := range univ {
				args[i] = NewVar(v)
			}
			sk = NewFunc(c.freshSkolem(), args...)
		}
		body := f.Subs[0].SubstituteTerm(f.Bound, sk)
		return c.skolemize(body, univ)
	case ConnAnd:
		return And(c.skolemize(f.Subs[0], univ), c.skolemize(f.Subs[1], univ))
	case ConnOr:
		return Or(c.skolemize(f.Subs[0], univ), c.skolemize(f.Subs[1], univ))
	default:
		return f
	}
}

// Clausify converts a first-order formula into an equisatisfiable set of clauses
// by NNF conversion, standardising variables apart, Skolemisation, dropping the
// (now universal) quantifiers and distributing to conjunctive normal form.
func Clausify(f Formula) []FOClause {
	c := &clausifier{}
	nnf := FONNF(f)
	std := c.standardize(nnf, map[string]string{})
	sk := c.skolemize(std, nil)
	matrix := distributeCNF(sk)
	var clauses []FOClause
	for _, conj := range splitConj(matrix) {
		var lits []FOLiteral
		taut := false
		for _, d := range splitDisj(conj) {
			l, ok := literalToFO(d)
			if !ok {
				continue
			}
			lits = append(lits, l)
		}
		// Drop clauses containing complementary literals (always true).
		for i := 0; i < len(lits) && !taut; i++ {
			for j := i + 1; j < len(lits); j++ {
				if lits[i].Neg != lits[j].Neg && lits[i].Pred == lits[j].Pred && TermsEqual(lits[i].Args, lits[j].Args) {
					taut = true
					break
				}
			}
		}
		if !taut {
			clauses = append(clauses, FOClause{Lits: lits})
		}
	}
	return clauses
}

// ClausifyAll clausifies a list of formulas into a single flat clause set.
func ClausifyAll(fs ...Formula) []FOClause {
	var out []FOClause
	for _, f := range fs {
		out = append(out, Clausify(f)...)
	}
	return out
}
