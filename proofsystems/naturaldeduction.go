package proofsystems

import (
	"errors"
	"fmt"
)

// NDRule identifies a natural-deduction inference rule in the Gentzen/Prawitz
// intuitionistic-plus-classical propositional calculus implemented here.
type NDRule int

const (
	// NDAssume introduces an assumption; it has no premises and its conclusion
	// becomes an open assumption.
	NDAssume NDRule = iota
	// NDAndI is conjunction introduction: from A and B infer A ∧ B.
	NDAndI
	// NDAndE1 is left conjunction elimination: from A ∧ B infer A.
	NDAndE1
	// NDAndE2 is right conjunction elimination: from A ∧ B infer B.
	NDAndE2
	// NDOrI1 is left disjunction introduction: from A infer A ∨ B.
	NDOrI1
	// NDOrI2 is right disjunction introduction: from B infer A ∨ B.
	NDOrI2
	// NDOrE is disjunction elimination (proof by cases).
	NDOrE
	// NDImpI is implication introduction, discharging an assumption.
	NDImpI
	// NDImpE is implication elimination (modus ponens).
	NDImpE
	// NDNotI is negation introduction: assume A, derive ⊥, conclude ¬A.
	NDNotI
	// NDNotE is negation elimination: from A and ¬A infer ⊥.
	NDNotE
	// NDBotE is ex falso quodlibet: from ⊥ infer any formula.
	NDBotE
	// NDRAA is reductio ad absurdum (classical): assume ¬A, derive ⊥,
	// conclude A.
	NDRAA
	// NDIffI is biconditional introduction from the two implications.
	NDIffI
	// NDIffE1 extracts the left-to-right implication from a biconditional.
	NDIffE1
	// NDIffE2 extracts the right-to-left implication from a biconditional.
	NDIffE2
	// NDTopI introduces the constant ⊤ with no premises.
	NDTopI
)

// String returns the conventional name of a natural-deduction rule.
func (r NDRule) String() string {
	switch r {
	case NDAssume:
		return "assume"
	case NDAndI:
		return "∧I"
	case NDAndE1:
		return "∧E1"
	case NDAndE2:
		return "∧E2"
	case NDOrI1:
		return "∨I1"
	case NDOrI2:
		return "∨I2"
	case NDOrE:
		return "∨E"
	case NDImpI:
		return "→I"
	case NDImpE:
		return "→E"
	case NDNotI:
		return "¬I"
	case NDNotE:
		return "¬E"
	case NDBotE:
		return "⊥E"
	case NDRAA:
		return "RAA"
	case NDIffI:
		return "↔I"
	case NDIffE1:
		return "↔E1"
	case NDIffE2:
		return "↔E2"
	case NDTopI:
		return "⊤I"
	default:
		return "?"
	}
}

// Derivation is a node of a natural-deduction proof tree. Rule names the
// inference, Concl is the formula it establishes, Premises are the
// sub-derivations of its immediate premises, and Discharge names the assumption
// formula released by a discharging rule (→I, ¬I, RAA and the two branches of
// ∨E).
type Derivation struct {
	Rule      NDRule
	Concl     Formula
	Premises  []Derivation
	Discharge Formula
}

// Assume constructs an assumption leaf whose conclusion is f.
func Assume(f Formula) Derivation { return Derivation{Rule: NDAssume, Concl: f} }

// AndIntro builds a conjunction-introduction node from proofs of the conjuncts.
func AndIntro(left, right Derivation) Derivation {
	return Derivation{Rule: NDAndI, Concl: And(left.Concl, right.Concl), Premises: []Derivation{left, right}}
}

// AndElim1 builds a left conjunction-elimination node concluding the left
// conjunct of the premise.
func AndElim1(d Derivation) Derivation {
	return Derivation{Rule: NDAndE1, Concl: d.Concl.Subs[0], Premises: []Derivation{d}}
}

// AndElim2 builds a right conjunction-elimination node concluding the right
// conjunct of the premise.
func AndElim2(d Derivation) Derivation {
	return Derivation{Rule: NDAndE2, Concl: d.Concl.Subs[1], Premises: []Derivation{d}}
}

// OrIntro1 builds a left disjunction-introduction node concluding d.Concl ∨
// other.
func OrIntro1(d Derivation, other Formula) Derivation {
	return Derivation{Rule: NDOrI1, Concl: Or(d.Concl, other), Premises: []Derivation{d}}
}

// OrIntro2 builds a right disjunction-introduction node concluding other ∨
// d.Concl.
func OrIntro2(other Formula, d Derivation) Derivation {
	return Derivation{Rule: NDOrI2, Concl: Or(other, d.Concl), Premises: []Derivation{d}}
}

// ImpIntro builds an implication-introduction node concluding assumption →
// body.Concl and discharging assumption from body.
func ImpIntro(assumption Formula, body Derivation) Derivation {
	return Derivation{Rule: NDImpI, Concl: Imp(assumption, body.Concl), Premises: []Derivation{body}, Discharge: assumption}
}

// ImpElim builds a modus-ponens node from proofs of A and A → B, concluding B.
func ImpElim(a, imp Derivation) Derivation {
	return Derivation{Rule: NDImpE, Concl: imp.Concl.Subs[1], Premises: []Derivation{a, imp}}
}

// NotIntro builds a negation-introduction node concluding ¬assumption from a
// derivation of ⊥, discharging assumption.
func NotIntro(assumption Formula, bottom Derivation) Derivation {
	return Derivation{Rule: NDNotI, Concl: Not(assumption), Premises: []Derivation{bottom}, Discharge: assumption}
}

// NotElim builds a negation-elimination node concluding ⊥ from A and ¬A.
func NotElim(a, notA Derivation) Derivation {
	return Derivation{Rule: NDNotE, Concl: Bot(), Premises: []Derivation{a, notA}}
}

// RAA builds a reductio node concluding goal from a derivation of ⊥,
// discharging the assumption ¬goal.
func RAA(goal Formula, bottom Derivation) Derivation {
	return Derivation{Rule: NDRAA, Concl: goal, Premises: []Derivation{bottom}, Discharge: Not(goal)}
}

// CheckND verifies a natural-deduction derivation and, on success, returns the
// multiset of undischarged (open) assumptions it depends on. An error is
// returned when any inference is not a valid instance of its rule.
func CheckND(d Derivation) ([]Formula, error) {
	switch d.Rule {
	case NDAssume:
		return []Formula{d.Concl}, nil
	case NDTopI:
		if d.Concl.Conn != ConnTrue {
			return nil, errors.New("proofsystems: ⊤I must conclude ⊤")
		}
		return nil, nil
	case NDAndI:
		return checkBinary(d, func(a, b Formula) error {
			if d.Concl.Conn != ConnAnd || !d.Concl.Subs[0].Equal(a) || !d.Concl.Subs[1].Equal(b) {
				return errUnsound("∧I", d)
			}
			return nil
		})
	case NDAndE1:
		return checkUnary(d, func(p Formula) error {
			if p.Conn != ConnAnd || !p.Subs[0].Equal(d.Concl) {
				return errUnsound("∧E1", d)
			}
			return nil
		})
	case NDAndE2:
		return checkUnary(d, func(p Formula) error {
			if p.Conn != ConnAnd || !p.Subs[1].Equal(d.Concl) {
				return errUnsound("∧E2", d)
			}
			return nil
		})
	case NDOrI1:
		return checkUnary(d, func(p Formula) error {
			if d.Concl.Conn != ConnOr || !d.Concl.Subs[0].Equal(p) {
				return errUnsound("∨I1", d)
			}
			return nil
		})
	case NDOrI2:
		return checkUnary(d, func(p Formula) error {
			if d.Concl.Conn != ConnOr || !d.Concl.Subs[1].Equal(p) {
				return errUnsound("∨I2", d)
			}
			return nil
		})
	case NDImpE:
		return checkBinary(d, func(a, imp Formula) error {
			if imp.Conn != ConnImp || !imp.Subs[0].Equal(a) || !imp.Subs[1].Equal(d.Concl) {
				return errUnsound("→E", d)
			}
			return nil
		})
	case NDNotE:
		return checkBinary(d, func(a, notA Formula) error {
			if notA.Conn != ConnNot || !notA.Subs[0].Equal(a) || d.Concl.Conn != ConnFalse {
				return errUnsound("¬E", d)
			}
			return nil
		})
	case NDBotE:
		return checkUnary(d, func(p Formula) error {
			if p.Conn != ConnFalse {
				return errUnsound("⊥E", d)
			}
			return nil
		})
	case NDImpI:
		if d.Concl.Conn != ConnImp || !d.Concl.Subs[0].Equal(d.Discharge) {
			return nil, errUnsound("→I", d)
		}
		open, err := CheckND(d.Premises[0])
		if err != nil {
			return nil, err
		}
		if !d.Premises[0].Concl.Equal(d.Concl.Subs[1]) {
			return nil, errUnsound("→I", d)
		}
		return discharge(open, d.Discharge), nil
	case NDNotI:
		if d.Concl.Conn != ConnNot || !d.Concl.Subs[0].Equal(d.Discharge) {
			return nil, errUnsound("¬I", d)
		}
		open, err := CheckND(d.Premises[0])
		if err != nil {
			return nil, err
		}
		if d.Premises[0].Concl.Conn != ConnFalse {
			return nil, errUnsound("¬I", d)
		}
		return discharge(open, d.Discharge), nil
	case NDRAA:
		open, err := CheckND(d.Premises[0])
		if err != nil {
			return nil, err
		}
		if d.Premises[0].Concl.Conn != ConnFalse || !d.Discharge.Equal(Not(d.Concl)) {
			return nil, errUnsound("RAA", d)
		}
		return discharge(open, d.Discharge), nil
	case NDIffI:
		return checkBinary(d, func(l, r Formula) error {
			if l.Conn != ConnImp || r.Conn != ConnImp || d.Concl.Conn != ConnIff {
				return errUnsound("↔I", d)
			}
			if !l.Subs[0].Equal(d.Concl.Subs[0]) || !l.Subs[1].Equal(d.Concl.Subs[1]) {
				return errUnsound("↔I", d)
			}
			if !r.Subs[0].Equal(d.Concl.Subs[1]) || !r.Subs[1].Equal(d.Concl.Subs[0]) {
				return errUnsound("↔I", d)
			}
			return nil
		})
	case NDIffE1:
		return checkUnary(d, func(p Formula) error {
			if p.Conn != ConnIff || d.Concl.Conn != ConnImp ||
				!d.Concl.Subs[0].Equal(p.Subs[0]) || !d.Concl.Subs[1].Equal(p.Subs[1]) {
				return errUnsound("↔E1", d)
			}
			return nil
		})
	case NDIffE2:
		return checkUnary(d, func(p Formula) error {
			if p.Conn != ConnIff || d.Concl.Conn != ConnImp ||
				!d.Concl.Subs[0].Equal(p.Subs[1]) || !d.Concl.Subs[1].Equal(p.Subs[0]) {
				return errUnsound("↔E2", d)
			}
			return nil
		})
	case NDOrE:
		return checkOrE(d)
	default:
		return nil, fmt.Errorf("proofsystems: unknown ND rule %v", d.Rule)
	}
}

// NDValid reports whether the derivation is well formed and its open assumptions
// are all contained in the given premise set, so it establishes premises ⊢
// Concl.
func NDValid(d Derivation, premises []Formula) bool {
	open, err := CheckND(d)
	if err != nil {
		return false
	}
	for _, a := range open {
		if !containsFormula(premises, a) {
			return false
		}
	}
	return true
}

// NDProves reports whether the derivation is a closed proof (no open
// assumptions) of its conclusion, i.e. a theorem.
func NDProves(d Derivation) bool {
	open, err := CheckND(d)
	return err == nil && len(open) == 0
}

func checkUnary(d Derivation, ok func(premise Formula) error) ([]Formula, error) {
	if len(d.Premises) != 1 {
		return nil, errUnsound(d.Rule.String(), d)
	}
	open, err := CheckND(d.Premises[0])
	if err != nil {
		return nil, err
	}
	if err := ok(d.Premises[0].Concl); err != nil {
		return nil, err
	}
	return open, nil
}

func checkBinary(d Derivation, ok func(a, b Formula) error) ([]Formula, error) {
	if len(d.Premises) != 2 {
		return nil, errUnsound(d.Rule.String(), d)
	}
	o1, err := CheckND(d.Premises[0])
	if err != nil {
		return nil, err
	}
	o2, err := CheckND(d.Premises[1])
	if err != nil {
		return nil, err
	}
	if err := ok(d.Premises[0].Concl, d.Premises[1].Concl); err != nil {
		return nil, err
	}
	return append(o1, o2...), nil
}

// checkOrE validates disjunction elimination: premise 0 proves A ∨ B, premise 1
// proves C from assumption A, premise 2 proves C from assumption B; the
// conclusion is C, and A and B are discharged from their respective branches.
func checkOrE(d Derivation) ([]Formula, error) {
	if len(d.Premises) != 3 {
		return nil, errUnsound("∨E", d)
	}
	disj, err := CheckND(d.Premises[0])
	if err != nil {
		return nil, err
	}
	left, err := CheckND(d.Premises[1])
	if err != nil {
		return nil, err
	}
	right, err := CheckND(d.Premises[2])
	if err != nil {
		return nil, err
	}
	dj := d.Premises[0].Concl
	if dj.Conn != ConnOr {
		return nil, errUnsound("∨E", d)
	}
	if !d.Premises[1].Concl.Equal(d.Concl) || !d.Premises[2].Concl.Equal(d.Concl) {
		return nil, errUnsound("∨E", d)
	}
	open := disj
	open = append(open, discharge(left, dj.Subs[0])...)
	open = append(open, discharge(right, dj.Subs[1])...)
	return open, nil
}

func discharge(open []Formula, a Formula) []Formula {
	var out []Formula
	for _, f := range open {
		if !f.Equal(a) {
			out = append(out, f)
		}
	}
	return out
}

func containsFormula(fs []Formula, f Formula) bool {
	for _, g := range fs {
		if g.Equal(f) {
			return true
		}
	}
	return false
}

func errUnsound(rule string, d Derivation) error {
	return fmt.Errorf("proofsystems: invalid %s inference concluding %s", rule, d.Concl.String())
}
