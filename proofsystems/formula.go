package proofsystems

import (
	"sort"
	"strings"
)

// Connective identifies the top-level logical operator of a formula. Atoms and
// the logical constants are treated uniformly as connectives with zero
// subformulas.
type Connective int

const (
	// ConnAtom marks an atomic formula P(t1,...,tn); a propositional variable
	// is the nullary case.
	ConnAtom Connective = iota
	// ConnTrue is the logical constant truth (⊤).
	ConnTrue
	// ConnFalse is the logical constant falsehood (⊥).
	ConnFalse
	// ConnNot is negation (¬).
	ConnNot
	// ConnAnd is conjunction (∧).
	ConnAnd
	// ConnOr is disjunction (∨).
	ConnOr
	// ConnImp is material implication (→).
	ConnImp
	// ConnIff is the biconditional (↔).
	ConnIff
	// ConnForall is the universal quantifier (∀).
	ConnForall
	// ConnExists is the existential quantifier (∃).
	ConnExists
)

// String returns a short mnemonic for the connective.
func (c Connective) String() string {
	switch c {
	case ConnAtom:
		return "atom"
	case ConnTrue:
		return "true"
	case ConnFalse:
		return "false"
	case ConnNot:
		return "not"
	case ConnAnd:
		return "and"
	case ConnOr:
		return "or"
	case ConnImp:
		return "imp"
	case ConnIff:
		return "iff"
	case ConnForall:
		return "forall"
	case ConnExists:
		return "exists"
	default:
		return "?"
	}
}

// Formula is an immutable propositional or first-order formula tree. The Conn
// field selects the interpretation of the remaining fields: for ConnAtom the
// Pred and Args fields hold a predicate symbol and its argument terms; for the
// unary and binary connectives the Subs field holds the operands; for the
// quantifiers the Bound field holds the bound variable name and Subs[0] the
// body.
type Formula struct {
	Conn  Connective
	Pred  string
	Args  []Term
	Subs  []Formula
	Bound string
}

// Prop builds a propositional atom (a nullary predicate) with the given name.
func Prop(name string) Formula {
	return Formula{Conn: ConnAtom, Pred: name}
}

// Atom builds a first-order atomic formula P(t1,...,tn).
func Atom(pred string, args ...Term) Formula {
	cp := make([]Term, len(args))
	copy(cp, args)
	return Formula{Conn: ConnAtom, Pred: pred, Args: cp}
}

// Top returns the logical constant truth.
func Top() Formula { return Formula{Conn: ConnTrue} }

// Bot returns the logical constant falsehood.
func Bot() Formula { return Formula{Conn: ConnFalse} }

// Not returns the negation of f.
func Not(f Formula) Formula { return Formula{Conn: ConnNot, Subs: []Formula{f}} }

// And returns the conjunction of a and b.
func And(a, b Formula) Formula { return Formula{Conn: ConnAnd, Subs: []Formula{a, b}} }

// Or returns the disjunction of a and b.
func Or(a, b Formula) Formula { return Formula{Conn: ConnOr, Subs: []Formula{a, b}} }

// Imp returns the implication a → b.
func Imp(a, b Formula) Formula { return Formula{Conn: ConnImp, Subs: []Formula{a, b}} }

// Iff returns the biconditional a ↔ b.
func Iff(a, b Formula) Formula { return Formula{Conn: ConnIff, Subs: []Formula{a, b}} }

// Forall returns the universally quantified formula ∀v. body.
func Forall(v string, body Formula) Formula {
	return Formula{Conn: ConnForall, Bound: v, Subs: []Formula{body}}
}

// Exists returns the existentially quantified formula ∃v. body.
func Exists(v string, body Formula) Formula {
	return Formula{Conn: ConnExists, Bound: v, Subs: []Formula{body}}
}

// Conj folds a slice of formulas into a right-nested conjunction. The empty
// slice yields Top and a singleton yields its element.
func Conj(fs ...Formula) Formula {
	if len(fs) == 0 {
		return Top()
	}
	acc := fs[len(fs)-1]
	for i := len(fs) - 2; i >= 0; i-- {
		acc = And(fs[i], acc)
	}
	return acc
}

// Disj folds a slice of formulas into a right-nested disjunction. The empty
// slice yields Bot and a singleton yields its element.
func Disj(fs ...Formula) Formula {
	if len(fs) == 0 {
		return Bot()
	}
	acc := fs[len(fs)-1]
	for i := len(fs) - 2; i >= 0; i-- {
		acc = Or(fs[i], acc)
	}
	return acc
}

// IsAtom reports whether the formula is an atomic formula or propositional
// variable.
func (f Formula) IsAtom() bool { return f.Conn == ConnAtom }

// IsLiteral reports whether the formula is an atom or the negation of an atom.
func (f Formula) IsLiteral() bool {
	if f.Conn == ConnAtom {
		return true
	}
	return f.Conn == ConnNot && f.Subs[0].Conn == ConnAtom
}

// IsConstant reports whether the formula is one of the logical constants ⊤ or
// ⊥.
func (f Formula) IsConstant() bool { return f.Conn == ConnTrue || f.Conn == ConnFalse }

// IsQuantified reports whether the top connective is a quantifier.
func (f Formula) IsQuantified() bool { return f.Conn == ConnForall || f.Conn == ConnExists }

// Equal reports whether two formulas are structurally identical (up to the
// literal names of bound variables — no alpha-conversion is performed).
func (f Formula) Equal(g Formula) bool {
	if f.Conn != g.Conn {
		return false
	}
	switch f.Conn {
	case ConnAtom:
		return f.Pred == g.Pred && TermsEqual(f.Args, g.Args)
	case ConnTrue, ConnFalse:
		return true
	case ConnForall, ConnExists:
		return f.Bound == g.Bound && f.Subs[0].Equal(g.Subs[0])
	default:
		if len(f.Subs) != len(g.Subs) {
			return false
		}
		for i := range f.Subs {
			if !f.Subs[i].Equal(g.Subs[i]) {
				return false
			}
		}
		return true
	}
}

// String renders the formula with minimal parenthesisation using ASCII
// operators (! & | -> <-> for the connectives and forall/exists for
// quantifiers).
func (f Formula) String() string {
	switch f.Conn {
	case ConnAtom:
		if len(f.Args) == 0 {
			return f.Pred
		}
		parts := make([]string, len(f.Args))
		for i, a := range f.Args {
			parts[i] = a.String()
		}
		return f.Pred + "(" + strings.Join(parts, ",") + ")"
	case ConnTrue:
		return "T"
	case ConnFalse:
		return "F"
	case ConnNot:
		return "!" + f.Subs[0].wrapped()
	case ConnAnd:
		return f.Subs[0].wrapped() + " & " + f.Subs[1].wrapped()
	case ConnOr:
		return f.Subs[0].wrapped() + " | " + f.Subs[1].wrapped()
	case ConnImp:
		return f.Subs[0].wrapped() + " -> " + f.Subs[1].wrapped()
	case ConnIff:
		return f.Subs[0].wrapped() + " <-> " + f.Subs[1].wrapped()
	case ConnForall:
		return "forall " + f.Bound + ". " + f.Subs[0].wrapped()
	case ConnExists:
		return "exists " + f.Bound + ". " + f.Subs[0].wrapped()
	default:
		return "?"
	}
}

func (f Formula) wrapped() string {
	switch f.Conn {
	case ConnAtom, ConnTrue, ConnFalse, ConnNot:
		return f.String()
	default:
		return "(" + f.String() + ")"
	}
}

// Atoms returns the sorted, de-duplicated list of atom string forms occurring
// in the formula. For propositional formulas this is the set of variable
// names.
func (f Formula) Atoms() []string {
	set := map[string]bool{}
	f.collectAtoms(set)
	return sortedKeys(set)
}

func (f Formula) collectAtoms(set map[string]bool) {
	if f.Conn == ConnAtom {
		set[f.atomKey()] = true
		return
	}
	for _, s := range f.Subs {
		s.collectAtoms(set)
	}
}

func (f Formula) atomKey() string {
	if len(f.Args) == 0 {
		return f.Pred
	}
	parts := make([]string, len(f.Args))
	for i, a := range f.Args {
		parts[i] = a.String()
	}
	return f.Pred + "(" + strings.Join(parts, ",") + ")"
}

// PropVars returns the sorted names of the propositional variables (nullary
// atoms) in the formula, ignoring first-order atoms with arguments.
func (f Formula) PropVars() []string {
	set := map[string]bool{}
	f.collectPropVars(set)
	return sortedKeys(set)
}

func (f Formula) collectPropVars(set map[string]bool) {
	if f.Conn == ConnAtom {
		if len(f.Args) == 0 {
			set[f.Pred] = true
		}
		return
	}
	for _, s := range f.Subs {
		s.collectPropVars(set)
	}
}

// Size returns the number of connective and atom nodes in the formula tree.
func (f Formula) Size() int {
	n := 1
	for _, s := range f.Subs {
		n += s.Size()
	}
	return n
}

// Depth returns the height of the formula tree; atoms and constants have depth
// 0.
func (f Formula) Depth() int {
	if len(f.Subs) == 0 {
		return 0
	}
	max := 0
	for _, s := range f.Subs {
		if d := s.Depth(); d > max {
			max = d
		}
	}
	return max + 1
}

// Subformulas returns every subformula including f itself, in pre-order, with
// structural duplicates removed.
func (f Formula) Subformulas() []Formula {
	var out []Formula
	seen := map[string]bool{}
	var walk func(g Formula)
	walk = func(g Formula) {
		key := g.String()
		if !seen[key] {
			seen[key] = true
			out = append(out, g)
		}
		for _, s := range g.Subs {
			walk(s)
		}
	}
	walk(f)
	return out
}

// FreeVars returns the sorted first-order variables occurring free in the
// formula (not bound by an enclosing quantifier).
func (f Formula) FreeVars() []string {
	set := map[string]bool{}
	f.collectFreeVars(map[string]bool{}, set)
	return sortedKeys(set)
}

func (f Formula) collectFreeVars(bound, free map[string]bool) {
	switch f.Conn {
	case ConnAtom:
		for _, a := range f.Args {
			for _, v := range a.Vars() {
				if !bound[v] {
					free[v] = true
				}
			}
		}
	case ConnForall, ConnExists:
		nb := map[string]bool{}
		for k := range bound {
			nb[k] = true
		}
		nb[f.Bound] = true
		f.Subs[0].collectFreeVars(nb, free)
	default:
		for _, s := range f.Subs {
			s.collectFreeVars(bound, free)
		}
	}
}

// SubstituteTerm returns the formula with every free occurrence of the
// first-order variable v replaced by the term t. Substitution stops at a
// quantifier that rebinds v, preserving capture avoidance for that variable.
func (f Formula) SubstituteTerm(v string, t Term) Formula {
	switch f.Conn {
	case ConnAtom:
		sub := SingletonSubstitution(v, t)
		return Formula{Conn: ConnAtom, Pred: f.Pred, Args: sub.ApplyTerms(f.Args)}
	case ConnForall, ConnExists:
		if f.Bound == v {
			return f
		}
		return Formula{Conn: f.Conn, Bound: f.Bound, Subs: []Formula{f.Subs[0].SubstituteTerm(v, t)}}
	default:
		subs := make([]Formula, len(f.Subs))
		for i, s := range f.Subs {
			subs[i] = s.SubstituteTerm(v, t)
		}
		return Formula{Conn: f.Conn, Pred: f.Pred, Args: f.Args, Subs: subs, Bound: f.Bound}
	}
}

// SubstituteProp returns the formula with every propositional atom named p
// replaced by the formula r.
func (f Formula) SubstituteProp(p string, r Formula) Formula {
	if f.Conn == ConnAtom {
		if len(f.Args) == 0 && f.Pred == p {
			return r
		}
		return f
	}
	subs := make([]Formula, len(f.Subs))
	for i, s := range f.Subs {
		subs[i] = s.SubstituteProp(p, r)
	}
	return Formula{Conn: f.Conn, Pred: f.Pred, Args: f.Args, Subs: subs, Bound: f.Bound}
}

// FormulasEqual reports whether two formula slices are pointwise structurally
// equal.
func FormulasEqual(a, b []Formula) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if !a[i].Equal(b[i]) {
			return false
		}
	}
	return true
}

// SortFormulas returns a copy of the slice sorted by string form, giving a
// canonical ordering useful for set comparisons.
func SortFormulas(fs []Formula) []Formula {
	out := make([]Formula, len(fs))
	copy(out, fs)
	sort.Slice(out, func(i, j int) bool { return out[i].String() < out[j].String() })
	return out
}
