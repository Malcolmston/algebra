package modelchecking

import (
	"sort"
	"strings"
)

// LTLKind enumerates the node kinds of a linear temporal logic syntax tree.
type LTLKind int

// The LTL node kinds. Boolean connectives, the unary future modalities and the
// binary temporal operators Until, Release and Weak-Until are all represented.
const (
	LTLTrueKind LTLKind = iota
	LTLFalseKind
	LTLAtomKind
	LTLNotKind
	LTLAndKind
	LTLOrKind
	LTLImpliesKind
	LTLIffKind
	LTLNextKind
	LTLUntilKind
	LTLReleaseKind
	LTLEventuallyKind
	LTLGloballyKind
	LTLWeakUntilKind
)

// LTL is a node of a linear temporal logic formula tree. Atom nodes carry a
// proposition name in Atom; unary nodes use L; binary nodes use L and R.
type LTL struct {
	Kind LTLKind
	Atom string
	L    *LTL
	R    *LTL
}

// LTLTop returns the constant true.
func LTLTop() *LTL { return &LTL{Kind: LTLTrueKind} }

// LTLBot returns the constant false.
func LTLBot() *LTL { return &LTL{Kind: LTLFalseKind} }

// LTLVar returns an atomic proposition node for name.
func LTLVar(name string) *LTL { return &LTL{Kind: LTLAtomKind, Atom: name} }

// LTLNot returns the negation of f.
func LTLNot(f *LTL) *LTL { return &LTL{Kind: LTLNotKind, L: f} }

// LTLAnd returns the conjunction f ∧ g.
func LTLAnd(f, g *LTL) *LTL { return &LTL{Kind: LTLAndKind, L: f, R: g} }

// LTLOr returns the disjunction f ∨ g.
func LTLOr(f, g *LTL) *LTL { return &LTL{Kind: LTLOrKind, L: f, R: g} }

// LTLImplies returns the implication f → g.
func LTLImplies(f, g *LTL) *LTL { return &LTL{Kind: LTLImpliesKind, L: f, R: g} }

// LTLIff returns the biconditional f ↔ g.
func LTLIff(f, g *LTL) *LTL { return &LTL{Kind: LTLIffKind, L: f, R: g} }

// LTLNext returns the next-time formula X f.
func LTLNext(f *LTL) *LTL { return &LTL{Kind: LTLNextKind, L: f} }

// LTLEventually returns the future formula F f (eventually f).
func LTLEventually(f *LTL) *LTL { return &LTL{Kind: LTLEventuallyKind, L: f} }

// LTLGlobally returns the always formula G f (globally f).
func LTLGlobally(f *LTL) *LTL { return &LTL{Kind: LTLGloballyKind, L: f} }

// LTLUntil returns the strong until f U g.
func LTLUntil(f, g *LTL) *LTL { return &LTL{Kind: LTLUntilKind, L: f, R: g} }

// LTLRelease returns the release f R g (the dual of Until).
func LTLRelease(f, g *LTL) *LTL { return &LTL{Kind: LTLReleaseKind, L: f, R: g} }

// LTLWeakUntil returns the weak until f W g, equivalent to (f U g) ∨ G f.
func LTLWeakUntil(f, g *LTL) *LTL { return &LTL{Kind: LTLWeakUntilKind, L: f, R: g} }

// IsLiteral reports whether f is true, false, an atom or the negation of an
// atom.
func (f *LTL) IsLiteral() bool {
	switch f.Kind {
	case LTLTrueKind, LTLFalseKind, LTLAtomKind:
		return true
	case LTLNotKind:
		return f.L != nil && f.L.Kind == LTLAtomKind
	}
	return false
}

// IsTemporal reports whether f's root is a temporal operator (X, F, G, U, R, W).
func (f *LTL) IsTemporal() bool {
	switch f.Kind {
	case LTLNextKind, LTLEventuallyKind, LTLGloballyKind,
		LTLUntilKind, LTLReleaseKind, LTLWeakUntilKind:
		return true
	}
	return false
}

// Size returns the number of nodes in the formula tree.
func (f *LTL) Size() int {
	if f == nil {
		return 0
	}
	return 1 + f.L.Size() + f.R.Size()
}

// Height returns the height of the formula tree (a single node has height 1).
func (f *LTL) Height() int {
	if f == nil {
		return 0
	}
	lh, rh := f.L.Height(), f.R.Height()
	if lh < rh {
		lh = rh
	}
	return 1 + lh
}

// Atoms returns the atomic proposition names occurring in f, sorted and
// de-duplicated.
func (f *LTL) Atoms() []string {
	set := map[string]bool{}
	var rec func(*LTL)
	rec = func(n *LTL) {
		if n == nil {
			return
		}
		if n.Kind == LTLAtomKind {
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
func (f *LTL) Equal(g *LTL) bool {
	if f == nil || g == nil {
		return f == nil && g == nil
	}
	if f.Kind != g.Kind || f.Atom != g.Atom {
		return false
	}
	return f.L.Equal(g.L) && f.R.Equal(g.R)
}

// Clone returns a deep copy of f.
func (f *LTL) Clone() *LTL {
	if f == nil {
		return nil
	}
	return &LTL{Kind: f.Kind, Atom: f.Atom, L: f.L.Clone(), R: f.R.Clone()}
}

// Negate returns a formula equivalent to ¬f without pushing the negation inward.
func (f *LTL) Negate() *LTL { return LTLNot(f) }

// NNF returns an equivalent formula in negation normal form: negations apply
// only to atoms, implications and biconditionals are eliminated, and the
// derived operators F, G and W are rewritten in terms of U and R. The result
// uses only the connectives ∧, ∨, X, U and R over literals, which is the input
// form expected by [LTLToGenBuchi].
func (f *LTL) NNF() *LTL { return nnf(f, false) }

// nnf converts f to negation normal form; neg indicates that the caller wants
// the negation of f.
func nnf(f *LTL, neg bool) *LTL {
	if f == nil {
		return nil
	}
	switch f.Kind {
	case LTLTrueKind:
		if neg {
			return LTLBot()
		}
		return LTLTop()
	case LTLFalseKind:
		if neg {
			return LTLTop()
		}
		return LTLBot()
	case LTLAtomKind:
		if neg {
			return LTLNot(LTLVar(f.Atom))
		}
		return LTLVar(f.Atom)
	case LTLNotKind:
		return nnf(f.L, !neg)
	case LTLAndKind:
		if neg {
			return LTLOr(nnf(f.L, true), nnf(f.R, true))
		}
		return LTLAnd(nnf(f.L, false), nnf(f.R, false))
	case LTLOrKind:
		if neg {
			return LTLAnd(nnf(f.L, true), nnf(f.R, true))
		}
		return LTLOr(nnf(f.L, false), nnf(f.R, false))
	case LTLImpliesKind:
		// f -> g  ≡  ¬f ∨ g
		if neg {
			return LTLAnd(nnf(f.L, false), nnf(f.R, true))
		}
		return LTLOr(nnf(f.L, true), nnf(f.R, false))
	case LTLIffKind:
		a, b := f.L, f.R
		if neg {
			// ¬(a↔b) ≡ (a ∧ ¬b) ∨ (¬a ∧ b)
			return LTLOr(
				LTLAnd(nnf(a, false), nnf(b, true)),
				LTLAnd(nnf(a, true), nnf(b, false)))
		}
		return LTLOr(
			LTLAnd(nnf(a, false), nnf(b, false)),
			LTLAnd(nnf(a, true), nnf(b, true)))
	case LTLNextKind:
		// X distributes through negation: ¬X f ≡ X ¬f
		return LTLNext(nnf(f.L, neg))
	case LTLEventuallyKind:
		// F g ≡ true U g ; ¬F g ≡ G ¬g ≡ false R ¬g
		if neg {
			return LTLRelease(LTLBot(), nnf(f.L, true))
		}
		return LTLUntil(LTLTop(), nnf(f.L, false))
	case LTLGloballyKind:
		// G g ≡ false R g ; ¬G g ≡ true U ¬g
		if neg {
			return LTLUntil(LTLTop(), nnf(f.L, true))
		}
		return LTLRelease(LTLBot(), nnf(f.L, false))
	case LTLUntilKind:
		// ¬(a U b) ≡ ¬a R ¬b
		if neg {
			return LTLRelease(nnf(f.L, true), nnf(f.R, true))
		}
		return LTLUntil(nnf(f.L, false), nnf(f.R, false))
	case LTLReleaseKind:
		// ¬(a R b) ≡ ¬a U ¬b
		if neg {
			return LTLUntil(nnf(f.L, true), nnf(f.R, true))
		}
		return LTLRelease(nnf(f.L, false), nnf(f.R, false))
	case LTLWeakUntilKind:
		// a W b ≡ b R (a ∨ b)
		eq := LTLRelease(f.R, LTLOr(f.L, f.R))
		return nnf(eq, neg)
	}
	return LTLTop()
}

// Subformulas returns every distinct subformula of f (including f itself) in a
// deterministic order by increasing size.
func (f *LTL) Subformulas() []*LTL {
	var all []*LTL
	seen := map[string]bool{}
	var rec func(*LTL)
	rec = func(n *LTL) {
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

// String renders f using ASCII operators: ! & | -> <-> X F G U R W.
func (f *LTL) String() string {
	if f == nil {
		return ""
	}
	switch f.Kind {
	case LTLTrueKind:
		return "true"
	case LTLFalseKind:
		return "false"
	case LTLAtomKind:
		return f.Atom
	case LTLNotKind:
		return "!" + f.L.paren()
	case LTLNextKind:
		return "X " + f.L.paren()
	case LTLEventuallyKind:
		return "F " + f.L.paren()
	case LTLGloballyKind:
		return "G " + f.L.paren()
	case LTLAndKind:
		return f.L.paren() + " & " + f.R.paren()
	case LTLOrKind:
		return f.L.paren() + " | " + f.R.paren()
	case LTLImpliesKind:
		return f.L.paren() + " -> " + f.R.paren()
	case LTLIffKind:
		return f.L.paren() + " <-> " + f.R.paren()
	case LTLUntilKind:
		return f.L.paren() + " U " + f.R.paren()
	case LTLReleaseKind:
		return f.L.paren() + " R " + f.R.paren()
	case LTLWeakUntilKind:
		return f.L.paren() + " W " + f.R.paren()
	}
	return "?"
}

func (f *LTL) paren() string {
	if f == nil {
		return ""
	}
	switch f.Kind {
	case LTLTrueKind, LTLFalseKind, LTLAtomKind:
		return f.String()
	case LTLNotKind, LTLNextKind, LTLEventuallyKind, LTLGloballyKind:
		if f.L.IsLiteral() || f.L.Kind == LTLAtomKind {
			return f.String()
		}
		return f.String()
	}
	return "(" + f.String() + ")"
}

// Simplify applies a set of validity-preserving Boolean and temporal
// simplifications (constant folding, idempotence, absorption of true/false)
// and returns the simplified formula. The original is not modified.
func (f *LTL) Simplify() *LTL {
	if f == nil {
		return nil
	}
	l := f.L.Simplify()
	r := f.R.Simplify()
	switch f.Kind {
	case LTLNotKind:
		if l.Kind == LTLTrueKind {
			return LTLBot()
		}
		if l.Kind == LTLFalseKind {
			return LTLTop()
		}
		if l.Kind == LTLNotKind {
			return l.L
		}
		return LTLNot(l)
	case LTLAndKind:
		if l.Kind == LTLFalseKind || r.Kind == LTLFalseKind {
			return LTLBot()
		}
		if l.Kind == LTLTrueKind {
			return r
		}
		if r.Kind == LTLTrueKind {
			return l
		}
		if l.Equal(r) {
			return l
		}
		return LTLAnd(l, r)
	case LTLOrKind:
		if l.Kind == LTLTrueKind || r.Kind == LTLTrueKind {
			return LTLTop()
		}
		if l.Kind == LTLFalseKind {
			return r
		}
		if r.Kind == LTLFalseKind {
			return l
		}
		if l.Equal(r) {
			return l
		}
		return LTLOr(l, r)
	case LTLImpliesKind:
		if l.Kind == LTLFalseKind || r.Kind == LTLTrueKind {
			return LTLTop()
		}
		if l.Kind == LTLTrueKind {
			return r
		}
		return LTLImplies(l, r)
	case LTLNextKind:
		if l.Kind == LTLTrueKind {
			return LTLTop()
		}
		if l.Kind == LTLFalseKind {
			return LTLBot()
		}
		return LTLNext(l)
	case LTLEventuallyKind:
		if l.Kind == LTLTrueKind {
			return LTLTop()
		}
		if l.Kind == LTLFalseKind {
			return LTLBot()
		}
		if l.Kind == LTLEventuallyKind {
			return l
		}
		return LTLEventually(l)
	case LTLGloballyKind:
		if l.Kind == LTLTrueKind {
			return LTLTop()
		}
		if l.Kind == LTLFalseKind {
			return LTLBot()
		}
		if l.Kind == LTLGloballyKind {
			return l
		}
		return LTLGlobally(l)
	case LTLUntilKind:
		if r.Kind == LTLTrueKind {
			return LTLTop()
		}
		if r.Kind == LTLFalseKind {
			return LTLBot()
		}
		if l.Kind == LTLFalseKind {
			return r
		}
		return LTLUntil(l, r)
	case LTLReleaseKind:
		if r.Kind == LTLFalseKind {
			return LTLBot()
		}
		if l.Kind == LTLTrueKind {
			return r
		}
		return LTLRelease(l, r)
	case LTLIffKind:
		return LTLIff(l, r)
	case LTLWeakUntilKind:
		return LTLWeakUntil(l, r)
	}
	return f.Clone()
}

// ltlKey returns a canonical string key for use as a map index.
func (f *LTL) ltlKey() string {
	if f == nil {
		return "∅"
	}
	var b strings.Builder
	var rec func(*LTL)
	rec = func(n *LTL) {
		if n == nil {
			b.WriteByte('_')
			return
		}
		b.WriteByte(byte('a' + int(n.Kind)))
		b.WriteString(n.Atom)
		b.WriteByte(':')
		rec(n.L)
		rec(n.R)
	}
	rec(f)
	return b.String()
}
