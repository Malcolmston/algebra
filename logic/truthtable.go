package logic

import (
	"fmt"
	"sort"
	"strings"
)

// TruthRow is a single line of a [TruthTable]: an assignment of truth values to
// the variables together with the resulting value of the expression.
type TruthRow struct {
	// Values maps each variable name to its assigned truth value.
	Values map[string]bool
	// Result is the value of the expression under Values.
	Result bool
}

// TruthTable is the exhaustive evaluation of an expression over all 2^n
// assignments of its n variables.
type TruthTable struct {
	// Vars holds the variable names in the column order used by the rows,
	// sorted lexicographically.
	Vars []string
	// Rows holds the 2^len(Vars) rows in ascending assignment-index order,
	// where the first variable is the most-significant bit.
	Rows []TruthRow
}

// NewTruthTable evaluates e over every assignment of its variables and returns
// the resulting truth table. Variables are ordered lexicographically and the
// rows run from the all-false assignment to the all-true one, treating the
// first variable as the most-significant bit. Eval errors cannot occur because
// every variable is bound.
func NewTruthTable(e Expr) *TruthTable {
	vars := Vars(e)
	n := len(vars)
	rows := make([]TruthRow, 0, 1<<n)
	for i := 0; i < (1 << n); i++ {
		env := IndexToAssignment(i, vars)
		res, _ := e.Eval(env)
		rows = append(rows, TruthRow{Values: env, Result: res})
	}
	return &TruthTable{Vars: vars, Rows: rows}
}

// IndexToAssignment decodes the row index i into a variable assignment over
// vars. Bit position len(vars)-1 (the least-significant bit) corresponds to the
// last variable, so vars[0] is the most-significant bit.
func IndexToAssignment(i int, vars []string) map[string]bool {
	n := len(vars)
	env := make(map[string]bool, n)
	for b := 0; b < n; b++ {
		// vars[b] occupies bit (n-1-b).
		env[vars[b]] = (i>>(n-1-b))&1 == 1
	}
	return env
}

// AssignmentToIndex encodes an assignment over vars back into its row index,
// the inverse of [IndexToAssignment]. A variable missing from env is treated as
// false.
func AssignmentToIndex(env map[string]bool, vars []string) int {
	n := len(vars)
	idx := 0
	for b := 0; b < n; b++ {
		if env[vars[b]] {
			idx |= 1 << (n - 1 - b)
		}
	}
	return idx
}

// Assignments enumerates all 2^len(vars) truth assignments over vars, in
// ascending index order. The input slice is sorted before enumeration so the
// ordering matches [NewTruthTable].
func Assignments(vars []string) []map[string]bool {
	sorted := append([]string(nil), vars...)
	sort.Strings(sorted)
	n := len(sorted)
	out := make([]map[string]bool, 0, 1<<n)
	for i := 0; i < (1 << n); i++ {
		out = append(out, IndexToAssignment(i, sorted))
	}
	return out
}

// Minterms returns the sorted indices of the assignments for which e is true.
// Each index is the minterm number using the lexicographic variable order of
// [Vars], with the first variable as the most-significant bit.
func Minterms(e Expr) []int {
	tt := NewTruthTable(e)
	var out []int
	for i, r := range tt.Rows {
		if r.Result {
			out = append(out, i)
		}
	}
	return out
}

// Maxterms returns the sorted indices of the assignments for which e is false,
// using the same indexing convention as [Minterms].
func Maxterms(e Expr) []int {
	tt := NewTruthTable(e)
	var out []int
	for i, r := range tt.Rows {
		if !r.Result {
			out = append(out, i)
		}
	}
	return out
}

// String renders the truth table as an aligned ASCII grid with a header row of
// variable names, a "=" result column, and one line per assignment using 1 and
// 0 for true and false.
func (t *TruthTable) String() string {
	var b strings.Builder
	widths := make([]int, len(t.Vars))
	for i, v := range t.Vars {
		widths[i] = len(v)
		if widths[i] < 1 {
			widths[i] = 1
		}
	}
	for i, v := range t.Vars {
		if i > 0 {
			b.WriteByte(' ')
		}
		fmt.Fprintf(&b, "%*s", widths[i], v)
	}
	b.WriteString(" | =\n")
	for _, r := range t.Rows {
		for i, v := range t.Vars {
			if i > 0 {
				b.WriteByte(' ')
			}
			fmt.Fprintf(&b, "%*s", widths[i], logicBit(r.Values[v]))
		}
		b.WriteString(" | ")
		b.WriteString(logicBit(r.Result))
		b.WriteByte('\n')
	}
	return b.String()
}

// logicBit renders a bool as "1" or "0".
func logicBit(b bool) string {
	if b {
		return "1"
	}
	return "0"
}
