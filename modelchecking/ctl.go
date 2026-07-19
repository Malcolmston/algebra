package modelchecking

import "sort"

// CTLKind enumerates the node kinds of a computation tree logic syntax tree.
type CTLKind int

// The CTL node kinds. Every temporal operator pairs a path quantifier (A: for
// all paths, E: there exists a path) with a temporal modality (X, F, G, U, R).
const (
	CTLTrueKind CTLKind = iota
	CTLFalseKind
	CTLAtomKind
	CTLNotKind
	CTLAndKind
	CTLOrKind
	CTLImpliesKind
	CTLIffKind
	CTLEXKind
	CTLEFKind
	CTLEGKind
	CTLEUKind
	CTLERKind
	CTLAXKind
	CTLAFKind
	CTLAGKind
	CTLAUKind
	CTLARKind
)

// CTL is a node of a computation tree logic formula tree. Atom nodes carry a
// proposition name; unary operators use L; binary operators (EU, ER, AU, AR and
// the Boolean binaries) use L and R.
type CTL struct {
	Kind CTLKind
	Atom string
	L    *CTL
	R    *CTL
}

// CTLTop returns the constant true.
func CTLTop() *CTL { return &CTL{Kind: CTLTrueKind} }

// CTLBot returns the constant false.
func CTLBot() *CTL { return &CTL{Kind: CTLFalseKind} }

// CTLVar returns an atomic proposition node for name.
func CTLVar(name string) *CTL { return &CTL{Kind: CTLAtomKind, Atom: name} }

// CTLNot returns the negation ¬f.
func CTLNot(f *CTL) *CTL { return &CTL{Kind: CTLNotKind, L: f} }

// CTLAnd returns the conjunction f ∧ g.
func CTLAnd(f, g *CTL) *CTL { return &CTL{Kind: CTLAndKind, L: f, R: g} }

// CTLOr returns the disjunction f ∨ g.
func CTLOr(f, g *CTL) *CTL { return &CTL{Kind: CTLOrKind, L: f, R: g} }

// CTLImplies returns the implication f → g.
func CTLImplies(f, g *CTL) *CTL { return &CTL{Kind: CTLImpliesKind, L: f, R: g} }

// CTLIff returns the biconditional f ↔ g.
func CTLIff(f, g *CTL) *CTL { return &CTL{Kind: CTLIffKind, L: f, R: g} }

// EX returns the formula EX f: some successor satisfies f.
func EX(f *CTL) *CTL { return &CTL{Kind: CTLEXKind, L: f} }

// EF returns the formula EF f: some path eventually reaches f.
func EF(f *CTL) *CTL { return &CTL{Kind: CTLEFKind, L: f} }

// EG returns the formula EG f: some path satisfies f globally.
func EG(f *CTL) *CTL { return &CTL{Kind: CTLEGKind, L: f} }

// EU returns the formula E[f U g]: some path satisfies f until g.
func EU(f, g *CTL) *CTL { return &CTL{Kind: CTLEUKind, L: f, R: g} }

// ER returns the formula E[f R g]: some path satisfies g released by f.
func ER(f, g *CTL) *CTL { return &CTL{Kind: CTLERKind, L: f, R: g} }

// AX returns the formula AX f: every successor satisfies f.
func AX(f *CTL) *CTL { return &CTL{Kind: CTLAXKind, L: f} }

// AF returns the formula AF f: every path eventually reaches f.
func AF(f *CTL) *CTL { return &CTL{Kind: CTLAFKind, L: f} }

// AG returns the formula AG f: every path satisfies f globally.
func AG(f *CTL) *CTL { return &CTL{Kind: CTLAGKind, L: f} }

// AU returns the formula A[f U g]: every path satisfies f until g.
func AU(f, g *CTL) *CTL { return &CTL{Kind: CTLAUKind, L: f, R: g} }

// AR returns the formula A[f R g]: every path satisfies g released by f.
func AR(f, g *CTL) *CTL { return &CTL{Kind: CTLARKind, L: f, R: g} }

// Size returns the number of nodes in the formula tree.
func (f *CTL) Size() int {
	if f == nil {
		return 0
	}
	return 1 + f.L.Size() + f.R.Size()
}

// Height returns the height of the formula tree.
func (f *CTL) Height() int {
	if f == nil {
		return 0
	}
	lh, rh := f.L.Height(), f.R.Height()
	if lh < rh {
		lh = rh
	}
	return 1 + lh
}

// Atoms returns the atomic proposition names occurring in f, sorted.
func (f *CTL) Atoms() []string {
	set := map[string]bool{}
	var rec func(*CTL)
	rec = func(n *CTL) {
		if n == nil {
			return
		}
		if n.Kind == CTLAtomKind {
			set[n.Atom] = true
		}
		rec(n.L)
		rec(n.R)
	}
	rec(f)
	out := make([]string, 0, len(set))
	for a := range set {
		out = append(out, a)
	}
	sort.Strings(out)
	return out
}

// Equal reports structural equality of two formulas.
func (f *CTL) Equal(g *CTL) bool {
	if f == nil || g == nil {
		return f == nil && g == nil
	}
	if f.Kind != g.Kind || f.Atom != g.Atom {
		return false
	}
	return f.L.Equal(g.L) && f.R.Equal(g.R)
}

// Clone returns a deep copy of f.
func (f *CTL) Clone() *CTL {
	if f == nil {
		return nil
	}
	return &CTL{Kind: f.Kind, Atom: f.Atom, L: f.L.Clone(), R: f.R.Clone()}
}

// IsTemporal reports whether f's root is a temporal (path-quantified) operator.
func (f *CTL) IsTemporal() bool {
	switch f.Kind {
	case CTLEXKind, CTLEFKind, CTLEGKind, CTLEUKind, CTLERKind,
		CTLAXKind, CTLAFKind, CTLAGKind, CTLAUKind, CTLARKind:
		return true
	}
	return false
}

// Subformulas returns every distinct subformula of f (including f) ordered by
// increasing size, the order in which the labelling model-checking algorithm
// processes them.
func (f *CTL) Subformulas() []*CTL {
	var all []*CTL
	seen := map[string]bool{}
	var rec func(*CTL)
	rec = func(n *CTL) {
		if n == nil {
			return
		}
		rec(n.L)
		rec(n.R)
		key := n.String()
		if !seen[key] {
			seen[key] = true
			all = append(all, n)
		}
	}
	rec(f)
	sort.SliceStable(all, func(i, j int) bool { return all[i].Size() < all[j].Size() })
	return all
}

// ExistentialNormalForm rewrites f into the existential fragment {¬, ∧, EX, EU,
// EG} plus atoms, which is the minimal set of operators the fixpoint checker
// implements directly. All A-quantified and derived operators are expressed via
// their existential duals. The result is logically equivalent to f.
func (f *CTL) ExistentialNormalForm() *CTL {
	if f == nil {
		return nil
	}
	switch f.Kind {
	case CTLTrueKind, CTLFalseKind, CTLAtomKind:
		return f.Clone()
	case CTLNotKind:
		return CTLNot(f.L.ExistentialNormalForm())
	case CTLAndKind:
		return CTLAnd(f.L.ExistentialNormalForm(), f.R.ExistentialNormalForm())
	case CTLOrKind:
		// a ∨ b ≡ ¬(¬a ∧ ¬b)
		return CTLNot(CTLAnd(CTLNot(f.L.ExistentialNormalForm()), CTLNot(f.R.ExistentialNormalForm())))
	case CTLImpliesKind:
		return CTLNot(CTLAnd(f.L.ExistentialNormalForm(), CTLNot(f.R.ExistentialNormalForm())))
	case CTLIffKind:
		a := f.L.ExistentialNormalForm()
		b := f.R.ExistentialNormalForm()
		// (a∧b) ∨ (¬a∧¬b)
		return CTLNot(CTLAnd(
			CTLNot(CTLAnd(a, b)),
			CTLNot(CTLAnd(CTLNot(a.Clone()), CTLNot(b.Clone())))))
	case CTLEXKind:
		return EX(f.L.ExistentialNormalForm())
	case CTLEUKind:
		return EU(f.L.ExistentialNormalForm(), f.R.ExistentialNormalForm())
	case CTLEGKind:
		return EG(f.L.ExistentialNormalForm())
	case CTLEFKind:
		// EF a ≡ E[true U a]
		return EU(CTLTop(), f.L.ExistentialNormalForm())
	case CTLERKind:
		// E[a R b] ≡ ¬A[¬a U ¬b] ≡ EG b ∧ ... use E[a R b] = ¬(A[¬a U ¬b])
		a := f.L.ExistentialNormalForm()
		b := f.R.ExistentialNormalForm()
		return existentialRelease(a, b)
	case CTLAXKind:
		// AX a ≡ ¬EX ¬a
		return CTLNot(EX(CTLNot(f.L.ExistentialNormalForm())))
	case CTLAFKind:
		// AF a ≡ ¬EG ¬a
		return CTLNot(EG(CTLNot(f.L.ExistentialNormalForm())))
	case CTLAGKind:
		// AG a ≡ ¬EF ¬a ≡ ¬E[true U ¬a]
		return CTLNot(EU(CTLTop(), CTLNot(f.L.ExistentialNormalForm())))
	case CTLAUKind:
		// A[a U b] ≡ ¬( E[¬b U (¬a∧¬b)] ∨ EG ¬b )
		a := f.L.ExistentialNormalForm()
		b := f.R.ExistentialNormalForm()
		notA := CTLNot(a)
		notB := CTLNot(b)
		left := EU(notB, CTLAnd(notA, CTLNot(b.Clone())))
		right := EG(CTLNot(b.Clone()))
		return CTLAnd(CTLNot(left), CTLNot(right))
	case CTLARKind:
		// A[a R b] ≡ ¬E[¬a U ¬b]
		a := f.L.ExistentialNormalForm()
		b := f.R.ExistentialNormalForm()
		return CTLNot(EU(CTLNot(a), CTLNot(b)))
	}
	return f.Clone()
}

// existentialRelease expresses E[a R b] using EU and EG:
// E[a R b] ≡ b ∧ ( a ∨ EX E[a R b] ). Equivalently, using the identity
// E[a R b] ≡ EG b ∨ E[b U (a ∧ b)].
func existentialRelease(a, b *CTL) *CTL {
	return CTLOr(EG(b.Clone()), EU(b.Clone(), CTLAnd(a.Clone(), b.Clone())))
}

// String renders f with ASCII operators, e.g. "AG (req -> AF ack)".
func (f *CTL) String() string {
	if f == nil {
		return ""
	}
	switch f.Kind {
	case CTLTrueKind:
		return "true"
	case CTLFalseKind:
		return "false"
	case CTLAtomKind:
		return f.Atom
	case CTLNotKind:
		return "!" + f.L.cparen()
	case CTLAndKind:
		return f.L.cparen() + " & " + f.R.cparen()
	case CTLOrKind:
		return f.L.cparen() + " | " + f.R.cparen()
	case CTLImpliesKind:
		return f.L.cparen() + " -> " + f.R.cparen()
	case CTLIffKind:
		return f.L.cparen() + " <-> " + f.R.cparen()
	case CTLEXKind:
		return "EX " + f.L.cparen()
	case CTLEFKind:
		return "EF " + f.L.cparen()
	case CTLEGKind:
		return "EG " + f.L.cparen()
	case CTLAXKind:
		return "AX " + f.L.cparen()
	case CTLAFKind:
		return "AF " + f.L.cparen()
	case CTLAGKind:
		return "AG " + f.L.cparen()
	case CTLEUKind:
		return "E[" + f.L.String() + " U " + f.R.String() + "]"
	case CTLERKind:
		return "E[" + f.L.String() + " R " + f.R.String() + "]"
	case CTLAUKind:
		return "A[" + f.L.String() + " U " + f.R.String() + "]"
	case CTLARKind:
		return "A[" + f.L.String() + " R " + f.R.String() + "]"
	}
	return "?"
}

func (f *CTL) cparen() string {
	if f == nil {
		return ""
	}
	switch f.Kind {
	case CTLTrueKind, CTLFalseKind, CTLAtomKind:
		return f.String()
	case CTLNotKind, CTLEXKind, CTLEFKind, CTLEGKind, CTLAXKind, CTLAFKind, CTLAGKind,
		CTLEUKind, CTLERKind, CTLAUKind, CTLARKind:
		return f.String()
	}
	return "(" + f.String() + ")"
}
