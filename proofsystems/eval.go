package proofsystems

import (
	"fmt"
	"sort"
	"strings"
)

// Assignment maps propositional atom names to Boolean truth values. Atoms
// absent from the map are taken to be false by Eval unless a total assignment
// is required.
type Assignment map[string]bool

// Clone returns an independent copy of the assignment.
func (a Assignment) Clone() Assignment {
	out := make(Assignment, len(a))
	for k, v := range a {
		out[k] = v
	}
	return out
}

// String renders an assignment as {A=true, B=false} with names sorted.
func (a Assignment) String() string {
	keys := make([]string, 0, len(a))
	for k := range a {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	parts := make([]string, len(keys))
	for i, k := range keys {
		parts[i] = fmt.Sprintf("%s=%v", k, a[k])
	}
	return "{" + strings.Join(parts, ", ") + "}"
}

// ErrEval is returned when a formula cannot be evaluated as a propositional
// formula (for example because it contains quantifiers).
var ErrEval = fmt.Errorf("proofsystems: cannot evaluate as propositional formula")

// Eval evaluates a propositional formula under the given assignment, treating
// each nullary atom as a Boolean variable. Missing atoms default to false. It
// panics-free and returns ErrEval for quantified or first-order atoms.
func Eval(f Formula, a Assignment) (bool, error) {
	switch f.Conn {
	case ConnAtom:
		if len(f.Args) != 0 {
			return false, ErrEval
		}
		return a[f.Pred], nil
	case ConnTrue:
		return true, nil
	case ConnFalse:
		return false, nil
	case ConnNot:
		v, err := Eval(f.Subs[0], a)
		return !v, err
	case ConnAnd:
		l, err := Eval(f.Subs[0], a)
		if err != nil {
			return false, err
		}
		r, err := Eval(f.Subs[1], a)
		return l && r, err
	case ConnOr:
		l, err := Eval(f.Subs[0], a)
		if err != nil {
			return false, err
		}
		r, err := Eval(f.Subs[1], a)
		return l || r, err
	case ConnImp:
		l, err := Eval(f.Subs[0], a)
		if err != nil {
			return false, err
		}
		r, err := Eval(f.Subs[1], a)
		return !l || r, err
	case ConnIff:
		l, err := Eval(f.Subs[0], a)
		if err != nil {
			return false, err
		}
		r, err := Eval(f.Subs[1], a)
		return l == r, err
	default:
		return false, ErrEval
	}
}

// MustEval evaluates f under a and panics on error. It is convenient for
// formulas known to be propositional.
func MustEval(f Formula, a Assignment) bool {
	v, err := Eval(f, a)
	if err != nil {
		panic(err)
	}
	return v
}

// AllAssignments enumerates every total truth assignment over the given atom
// names in ascending binary order (the first atom is the most significant bit).
func AllAssignments(atoms []string) []Assignment {
	names := make([]string, len(atoms))
	copy(names, atoms)
	sort.Strings(names)
	n := len(names)
	total := 1 << n
	out := make([]Assignment, 0, total)
	for mask := 0; mask < total; mask++ {
		a := make(Assignment, n)
		for i, name := range names {
			bit := (mask >> (n - 1 - i)) & 1
			a[name] = bit == 1
		}
		out = append(out, a)
	}
	return out
}

// TruthRow is a single line of a truth table: the assignment and the resulting
// value of the formula.
type TruthRow struct {
	Assignment Assignment
	Value      bool
}

// TruthTable holds the complete truth table of a propositional formula.
type TruthTable struct {
	Atoms []string
	Rows  []TruthRow
}

// NewTruthTable builds the full truth table of a propositional formula over its
// atoms. It returns ErrEval if the formula is not propositional.
func NewTruthTable(f Formula) (TruthTable, error) {
	atoms := f.PropVars()
	tt := TruthTable{Atoms: atoms}
	for _, a := range AllAssignments(atoms) {
		v, err := Eval(f, a)
		if err != nil {
			return TruthTable{}, err
		}
		tt.Rows = append(tt.Rows, TruthRow{Assignment: a, Value: v})
	}
	return tt, nil
}

// String renders the truth table as a fixed-width grid.
func (tt TruthTable) String() string {
	var b strings.Builder
	for _, a := range tt.Atoms {
		b.WriteString(a)
		b.WriteString(" ")
	}
	b.WriteString("| value\n")
	for _, r := range tt.Rows {
		for _, a := range tt.Atoms {
			if r.Assignment[a] {
				b.WriteString("T")
			} else {
				b.WriteString("F")
			}
			b.WriteString(strings.Repeat(" ", len(a)))
		}
		b.WriteString("| ")
		if r.Value {
			b.WriteString("T")
		} else {
			b.WriteString("F")
		}
		b.WriteString("\n")
	}
	return b.String()
}

// TrueRows returns the assignments under which the formula is true.
func (tt TruthTable) TrueRows() []Assignment {
	var out []Assignment
	for _, r := range tt.Rows {
		if r.Value {
			out = append(out, r.Assignment)
		}
	}
	return out
}

// CountTrue returns the number of assignments under which the formula is true.
func (tt TruthTable) CountTrue() int {
	n := 0
	for _, r := range tt.Rows {
		if r.Value {
			n++
		}
	}
	return n
}
