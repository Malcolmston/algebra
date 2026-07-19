package satsolver

import (
	"sort"
	"strings"
)

// TruthRow is a single line of a [TruthTable]: an assignment to the variables
// and the resulting value of the expression.
type TruthRow struct {
	// Assignment maps each variable name to its truth value in this row.
	Assignment map[string]bool
	// Value is the value of the expression under the assignment.
	Value bool
}

// TruthTable is the exhaustive truth table of a Boolean expression over its
// variables, ordered so that the first variable is the most significant bit of
// the row index.
type TruthTable struct {
	// Vars are the variable names in most-significant-first order.
	Vars []string
	// Rows are the 2^len(Vars) rows in ascending index order.
	Rows []TruthRow
}

// NewTruthTable builds the complete truth table of e. Variables are taken in
// sorted order.
func NewTruthTable(e Expr) TruthTable {
	vars := Vars(e)
	sort.Strings(vars)
	n := len(vars)
	rows := make([]TruthRow, 0, 1<<n)
	for i := 0; i < (1 << n); i++ {
		env := indexEnv(vars, i)
		rows = append(rows, TruthRow{Assignment: env, Value: e.Eval(env)})
	}
	return TruthTable{Vars: vars, Rows: rows}
}

// indexEnv maps row index i to an environment where the first variable is the
// most significant bit.
func indexEnv(vars []string, i int) map[string]bool {
	n := len(vars)
	env := make(map[string]bool, n)
	for k, v := range vars {
		bit := (i >> (n - 1 - k)) & 1
		env[v] = bit == 1
	}
	return env
}

// NumRows returns the number of rows, which is 2 raised to the number of
// variables.
func (t TruthTable) NumRows() int { return len(t.Rows) }

// Minterms returns the row indices at which the expression is true.
func (t TruthTable) Minterms() []int {
	var out []int
	for i, r := range t.Rows {
		if r.Value {
			out = append(out, i)
		}
	}
	return out
}

// Maxterms returns the row indices at which the expression is false.
func (t TruthTable) Maxterms() []int {
	var out []int
	for i, r := range t.Rows {
		if !r.Value {
			out = append(out, i)
		}
	}
	return out
}

// Column returns the output column of the table as a slice of booleans in row
// order.
func (t TruthTable) Column() []bool {
	out := make([]bool, len(t.Rows))
	for i, r := range t.Rows {
		out[i] = r.Value
	}
	return out
}

// String renders the truth table as fixed-width text with a header row.
func (t TruthTable) String() string {
	var b strings.Builder
	for _, v := range t.Vars {
		b.WriteString(v)
		b.WriteByte(' ')
	}
	b.WriteString("| f\n")
	for i, r := range t.Rows {
		env := indexEnv(t.Vars, i)
		for _, v := range t.Vars {
			if env[v] {
				b.WriteString("1")
			} else {
				b.WriteString("0")
			}
			b.WriteString(strings.Repeat(" ", len(v)))
		}
		b.WriteString("| ")
		if r.Value {
			b.WriteString("1")
		} else {
			b.WriteString("0")
		}
		b.WriteByte('\n')
	}
	return b.String()
}

// Minterms returns the set of minterm indices of e: the assignments (as integer
// row indices with the first sorted variable most significant) where e is true.
func Minterms(e Expr) []int { return NewTruthTable(e).Minterms() }

// Maxterms returns the maxterm indices of e: the assignments where e is false.
func Maxterms(e Expr) []int { return NewTruthTable(e).Maxterms() }

// IndexToEnv converts a row index to a variable environment for the given
// sorted variable list, first variable most significant.
func IndexToEnv(vars []string, i int) map[string]bool { return indexEnv(vars, i) }

// EnvToIndex converts a variable environment to its row index for the given
// sorted variable list, first variable most significant.
func EnvToIndex(vars []string, env map[string]bool) int {
	n := len(vars)
	idx := 0
	for k, v := range vars {
		if env[v] {
			idx |= 1 << (n - 1 - k)
		}
	}
	return idx
}
