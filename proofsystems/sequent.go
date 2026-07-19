package proofsystems

import (
	"fmt"
	"sort"
	"strings"
)

// Sequent is a two-sided propositional sequent Γ ⊢ Δ where Left is the
// antecedent (assumed formulas) and Right is the succedent (goal
// alternatives). Its intended meaning is that the conjunction of Left entails
// the disjunction of Right.
type Sequent struct {
	Left  []Formula
	Right []Formula
}

// NewSequent builds a sequent from antecedent and succedent slices.
func NewSequent(left, right []Formula) Sequent {
	l := make([]Formula, len(left))
	copy(l, left)
	r := make([]Formula, len(right))
	copy(r, right)
	return Sequent{Left: l, Right: r}
}

// String renders the sequent using the ⊢ turnstile with comma-separated sides.
func (s Sequent) String() string {
	return strings.Join(formulaStrings(s.Left), ", ") + " |- " + strings.Join(formulaStrings(s.Right), ", ")
}

func formulaStrings(fs []Formula) []string {
	out := make([]string, len(fs))
	for i, f := range fs {
		out[i] = f.String()
	}
	return out
}

// IsAxiom reports whether the sequent is an axiom instance, i.e. some formula
// occurs on both sides.
func (s Sequent) IsAxiom() bool {
	for _, l := range s.Left {
		for _, r := range s.Right {
			if l.Equal(r) {
				return true
			}
		}
	}
	// ⊥ on the left or ⊤ on the right also closes the sequent.
	for _, l := range s.Left {
		if l.Conn == ConnFalse {
			return true
		}
	}
	for _, r := range s.Right {
		if r.Conn == ConnTrue {
			return true
		}
	}
	return false
}

// SeqRule identifies a rule of the propositional sequent calculus LK.
type SeqRule int

const (
	// SeqAxiom is the identity axiom Γ, A ⊢ A, Δ.
	SeqAxiom SeqRule = iota
	// SeqNotL is the left negation rule.
	SeqNotL
	// SeqNotR is the right negation rule.
	SeqNotR
	// SeqAndL is the left conjunction rule.
	SeqAndL
	// SeqAndR is the right conjunction rule.
	SeqAndR
	// SeqOrL is the left disjunction rule.
	SeqOrL
	// SeqOrR is the right disjunction rule.
	SeqOrR
	// SeqImpL is the left implication rule.
	SeqImpL
	// SeqImpR is the right implication rule.
	SeqImpR
	// SeqIffL is the left biconditional rule.
	SeqIffL
	// SeqIffR is the right biconditional rule.
	SeqIffR
)

// String returns the conventional name of the sequent rule.
func (r SeqRule) String() string {
	switch r {
	case SeqAxiom:
		return "Ax"
	case SeqNotL:
		return "¬L"
	case SeqNotR:
		return "¬R"
	case SeqAndL:
		return "∧L"
	case SeqAndR:
		return "∧R"
	case SeqOrL:
		return "∨L"
	case SeqOrR:
		return "∨R"
	case SeqImpL:
		return "→L"
	case SeqImpR:
		return "→R"
	case SeqIffL:
		return "↔L"
	case SeqIffR:
		return "↔R"
	default:
		return "?"
	}
}

// SequentDerivation is a node of a sequent-calculus proof tree: an inference of
// Concl from the conclusions of its premise sub-derivations by Rule.
type SequentDerivation struct {
	Rule     SeqRule
	Concl    Sequent
	Premises []SequentDerivation
}

// CheckSequent verifies that a sequent derivation is a valid LK proof: every
// leaf is an axiom and every internal node is a correct instance of its rule.
// It returns nil on success or a descriptive error.
func CheckSequent(d SequentDerivation) error {
	if d.Rule == SeqAxiom {
		if !d.Concl.IsAxiom() {
			return fmt.Errorf("proofsystems: %s is not an axiom", d.Concl)
		}
		return nil
	}
	for i := range d.Premises {
		if err := CheckSequent(d.Premises[i]); err != nil {
			return err
		}
	}
	prem := make([]Sequent, len(d.Premises))
	for i, p := range d.Premises {
		prem[i] = p.Concl
	}
	if !validSeqStep(d.Rule, d.Concl, prem) {
		return fmt.Errorf("proofsystems: invalid %s step at %s", d.Rule, d.Concl)
	}
	return nil
}

// validSeqStep reports whether applying rule to some principal formula of concl
// yields exactly the given premises (as multisets). It searches candidate
// principal formulas so the check is independent of formula ordering.
func validSeqStep(rule SeqRule, concl Sequent, prem []Sequent) bool {
	switch rule {
	case SeqNotL:
		for i, f := range concl.Left {
			if f.Conn != ConnNot {
				continue
			}
			gamma := removeAt(concl.Left, i)
			// premise: gamma ⊢ f.Subs[0], Δ
			want := Sequent{Left: gamma, Right: append([]Formula{f.Subs[0]}, concl.Right...)}
			if len(prem) == 1 && seqEqual(prem[0], want) {
				return true
			}
		}
	case SeqNotR:
		for i, f := range concl.Right {
			if f.Conn != ConnNot {
				continue
			}
			delta := removeAt(concl.Right, i)
			want := Sequent{Left: append([]Formula{f.Subs[0]}, concl.Left...), Right: delta}
			if len(prem) == 1 && seqEqual(prem[0], want) {
				return true
			}
		}
	case SeqAndL:
		for i, f := range concl.Left {
			if f.Conn != ConnAnd {
				continue
			}
			gamma := removeAt(concl.Left, i)
			want := Sequent{Left: append([]Formula{f.Subs[0], f.Subs[1]}, gamma...), Right: concl.Right}
			if len(prem) == 1 && seqEqual(prem[0], want) {
				return true
			}
		}
	case SeqAndR:
		for i, f := range concl.Right {
			if f.Conn != ConnAnd {
				continue
			}
			delta := removeAt(concl.Right, i)
			w1 := Sequent{Left: concl.Left, Right: append([]Formula{f.Subs[0]}, delta...)}
			w2 := Sequent{Left: concl.Left, Right: append([]Formula{f.Subs[1]}, delta...)}
			if len(prem) == 2 && matchTwo(prem, w1, w2) {
				return true
			}
		}
	case SeqOrL:
		for i, f := range concl.Left {
			if f.Conn != ConnOr {
				continue
			}
			gamma := removeAt(concl.Left, i)
			w1 := Sequent{Left: append([]Formula{f.Subs[0]}, gamma...), Right: concl.Right}
			w2 := Sequent{Left: append([]Formula{f.Subs[1]}, gamma...), Right: concl.Right}
			if len(prem) == 2 && matchTwo(prem, w1, w2) {
				return true
			}
		}
	case SeqOrR:
		for i, f := range concl.Right {
			if f.Conn != ConnOr {
				continue
			}
			delta := removeAt(concl.Right, i)
			want := Sequent{Left: concl.Left, Right: append([]Formula{f.Subs[0], f.Subs[1]}, delta...)}
			if len(prem) == 1 && seqEqual(prem[0], want) {
				return true
			}
		}
	case SeqImpL:
		for i, f := range concl.Left {
			if f.Conn != ConnImp {
				continue
			}
			gamma := removeAt(concl.Left, i)
			w1 := Sequent{Left: gamma, Right: append([]Formula{f.Subs[0]}, concl.Right...)}
			w2 := Sequent{Left: append([]Formula{f.Subs[1]}, gamma...), Right: concl.Right}
			if len(prem) == 2 && matchTwo(prem, w1, w2) {
				return true
			}
		}
	case SeqImpR:
		for i, f := range concl.Right {
			if f.Conn != ConnImp {
				continue
			}
			delta := removeAt(concl.Right, i)
			want := Sequent{Left: append([]Formula{f.Subs[0]}, concl.Left...), Right: append([]Formula{f.Subs[1]}, delta...)}
			if len(prem) == 1 && seqEqual(prem[0], want) {
				return true
			}
		}
	case SeqIffL:
		for i, f := range concl.Left {
			if f.Conn != ConnIff {
				continue
			}
			gamma := removeAt(concl.Left, i)
			a, b := f.Subs[0], f.Subs[1]
			w1 := Sequent{Left: append([]Formula{a, b}, gamma...), Right: concl.Right}
			w2 := Sequent{Left: gamma, Right: append([]Formula{a, b}, concl.Right...)}
			if len(prem) == 2 && matchTwo(prem, w1, w2) {
				return true
			}
		}
	case SeqIffR:
		for i, f := range concl.Right {
			if f.Conn != ConnIff {
				continue
			}
			delta := removeAt(concl.Right, i)
			a, b := f.Subs[0], f.Subs[1]
			w1 := Sequent{Left: append([]Formula{a}, concl.Left...), Right: append([]Formula{b}, delta...)}
			w2 := Sequent{Left: append([]Formula{b}, concl.Left...), Right: append([]Formula{a}, delta...)}
			if len(prem) == 2 && matchTwo(prem, w1, w2) {
				return true
			}
		}
	}
	return false
}

// ValidSequent reports whether a propositional sequent is valid, i.e. the
// conjunction of its antecedent entails the disjunction of its succedent.
func ValidSequent(s Sequent) bool {
	ant := Conj(s.Left...)
	if len(s.Left) == 0 {
		ant = Top()
	}
	suc := Disj(s.Right...)
	if len(s.Right) == 0 {
		suc = Bot()
	}
	return IsTautology(Imp(ant, suc))
}

// ProveSequent searches for an LK proof of the sequent using backward
// application of the invertible propositional rules. Because those rules are
// invertible and strictly reduce the connective count, the search always
// terminates; it returns a checkable derivation and true when the sequent is
// valid, or the zero derivation and false otherwise.
func ProveSequent(s Sequent) (SequentDerivation, bool) {
	if s.IsAxiom() {
		return SequentDerivation{Rule: SeqAxiom, Concl: s}, true
	}
	// Try to decompose a non-atomic formula on either side.
	if i, ok := firstCompound(s.Left); ok {
		f := s.Left[i]
		gamma := removeAt(s.Left, i)
		switch f.Conn {
		case ConnNot:
			sub := Sequent{Left: gamma, Right: append([]Formula{f.Subs[0]}, s.Right...)}
			return buildUnary(SeqNotL, s, sub)
		case ConnAnd:
			sub := Sequent{Left: append([]Formula{f.Subs[0], f.Subs[1]}, gamma...), Right: s.Right}
			return buildUnary(SeqAndL, s, sub)
		case ConnOr:
			s1 := Sequent{Left: append([]Formula{f.Subs[0]}, gamma...), Right: s.Right}
			s2 := Sequent{Left: append([]Formula{f.Subs[1]}, gamma...), Right: s.Right}
			return buildBinary(SeqOrL, s, s1, s2)
		case ConnImp:
			s1 := Sequent{Left: gamma, Right: append([]Formula{f.Subs[0]}, s.Right...)}
			s2 := Sequent{Left: append([]Formula{f.Subs[1]}, gamma...), Right: s.Right}
			return buildBinary(SeqImpL, s, s1, s2)
		case ConnIff:
			s1 := Sequent{Left: append([]Formula{f.Subs[0], f.Subs[1]}, gamma...), Right: s.Right}
			s2 := Sequent{Left: gamma, Right: append([]Formula{f.Subs[0], f.Subs[1]}, s.Right...)}
			return buildBinary(SeqIffL, s, s1, s2)
		}
	}
	if i, ok := firstCompound(s.Right); ok {
		f := s.Right[i]
		delta := removeAt(s.Right, i)
		switch f.Conn {
		case ConnNot:
			sub := Sequent{Left: append([]Formula{f.Subs[0]}, s.Left...), Right: delta}
			return buildUnary(SeqNotR, s, sub)
		case ConnOr:
			sub := Sequent{Left: s.Left, Right: append([]Formula{f.Subs[0], f.Subs[1]}, delta...)}
			return buildUnary(SeqOrR, s, sub)
		case ConnImp:
			sub := Sequent{Left: append([]Formula{f.Subs[0]}, s.Left...), Right: append([]Formula{f.Subs[1]}, delta...)}
			return buildUnary(SeqImpR, s, sub)
		case ConnAnd:
			s1 := Sequent{Left: s.Left, Right: append([]Formula{f.Subs[0]}, delta...)}
			s2 := Sequent{Left: s.Left, Right: append([]Formula{f.Subs[1]}, delta...)}
			return buildBinary(SeqAndR, s, s1, s2)
		case ConnIff:
			s1 := Sequent{Left: append([]Formula{f.Subs[0]}, s.Left...), Right: append([]Formula{f.Subs[1]}, delta...)}
			s2 := Sequent{Left: append([]Formula{f.Subs[1]}, s.Left...), Right: append([]Formula{f.Subs[0]}, delta...)}
			return buildBinary(SeqIffR, s, s1, s2)
		}
	}
	return SequentDerivation{}, false
}

// SequentProvable reports whether ProveSequent finds a proof of the sequent.
func SequentProvable(s Sequent) bool {
	_, ok := ProveSequent(s)
	return ok
}

// ProveFormula attempts to prove the propositional formula f as the sequent
// ⊢ f, returning a sequent derivation when f is valid.
func ProveFormula(f Formula) (SequentDerivation, bool) {
	return ProveSequent(Sequent{Right: []Formula{f}})
}

func buildUnary(rule SeqRule, concl, sub Sequent) (SequentDerivation, bool) {
	d, ok := ProveSequent(sub)
	if !ok {
		return SequentDerivation{}, false
	}
	return SequentDerivation{Rule: rule, Concl: concl, Premises: []SequentDerivation{d}}, true
}

func buildBinary(rule SeqRule, concl, s1, s2 Sequent) (SequentDerivation, bool) {
	d1, ok := ProveSequent(s1)
	if !ok {
		return SequentDerivation{}, false
	}
	d2, ok := ProveSequent(s2)
	if !ok {
		return SequentDerivation{}, false
	}
	return SequentDerivation{Rule: rule, Concl: concl, Premises: []SequentDerivation{d1, d2}}, true
}

func firstCompound(fs []Formula) (int, bool) {
	for i, f := range fs {
		switch f.Conn {
		case ConnAtom, ConnTrue, ConnFalse:
			continue
		default:
			return i, true
		}
	}
	return 0, false
}

func removeAt(fs []Formula, i int) []Formula {
	out := make([]Formula, 0, len(fs)-1)
	out = append(out, fs[:i]...)
	out = append(out, fs[i+1:]...)
	return out
}

func seqEqual(a, b Sequent) bool {
	return multisetEqual(a.Left, b.Left) && multisetEqual(a.Right, b.Right)
}

func matchTwo(prem []Sequent, w1, w2 Sequent) bool {
	return (seqEqual(prem[0], w1) && seqEqual(prem[1], w2)) ||
		(seqEqual(prem[0], w2) && seqEqual(prem[1], w1))
}

func multisetEqual(a, b []Formula) bool {
	if len(a) != len(b) {
		return false
	}
	as := formulaStrings(a)
	bs := formulaStrings(b)
	sort.Strings(as)
	sort.Strings(bs)
	for i := range as {
		if as[i] != bs[i] {
			return false
		}
	}
	return true
}
