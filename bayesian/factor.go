package bayesian

import (
	"errors"
	"sort"
)

// ErrFactor is returned for malformed factor operations, such as a table whose
// length does not match the product of the variable cardinalities.
var ErrFactor = errors.New("bayesian: invalid factor")

// Factor is a discrete potential over a set of named categorical variables. The
// Table holds one non-negative value per joint assignment, laid out in
// row-major order with the first variable most significant and the last least
// significant.
type Factor struct {
	Vars  []string
	Card  []int
	Table []float64
}

// NewFactor constructs a Factor from parallel variable-name and cardinality
// slices and a value table. It returns ErrFactor when the shapes are
// inconsistent.
func NewFactor(vars []string, card []int, table []float64) (Factor, error) {
	if len(vars) != len(card) {
		return Factor{}, ErrFactor
	}
	size := 1
	for _, c := range card {
		if c <= 0 {
			return Factor{}, ErrFactor
		}
		size *= c
	}
	if len(table) != size {
		return Factor{}, ErrFactor
	}
	v := make([]string, len(vars))
	copy(v, vars)
	cc := make([]int, len(card))
	copy(cc, card)
	t := make([]float64, len(table))
	copy(t, table)
	return Factor{Vars: v, Card: cc, Table: t}, nil
}

// Size returns the number of entries in the factor's table.
func (f Factor) Size() int { return len(f.Table) }

// strides returns the row-major strides for each variable.
func (f Factor) strides() []int {
	n := len(f.Card)
	s := make([]int, n)
	acc := 1
	for i := n - 1; i >= 0; i-- {
		s[i] = acc
		acc *= f.Card[i]
	}
	return s
}

// indexToAssignment converts a flat table index into a per-variable assignment.
func (f Factor) indexToAssignment(idx int) []int {
	s := f.strides()
	a := make([]int, len(f.Card))
	for i := range f.Card {
		a[i] = (idx / s[i]) % f.Card[i]
	}
	return a
}

// assignmentToIndex converts a per-variable assignment into a flat table index.
func (f Factor) assignmentToIndex(a []int) int {
	s := f.strides()
	idx := 0
	for i := range a {
		idx += a[i] * s[i]
	}
	return idx
}

// varIndex returns the position of variable name within f.Vars, or −1.
func (f Factor) varIndex(name string) int {
	for i, v := range f.Vars {
		if v == name {
			return i
		}
	}
	return -1
}

// Get returns the factor value for a full or partial assignment given as a map
// from variable name to category index. Every variable of the factor must be
// present in assign; a missing variable yields 0.
func (f Factor) Get(assign map[string]int) float64 {
	a := make([]int, len(f.Vars))
	for i, v := range f.Vars {
		val, ok := assign[v]
		if !ok || val < 0 || val >= f.Card[i] {
			return 0
		}
		a[i] = val
	}
	return f.Table[f.assignmentToIndex(a)]
}

// Sum returns the total mass of the factor.
func (f Factor) Sum() float64 {
	var s float64
	for _, v := range f.Table {
		s += v
	}
	return s
}

// Normalize returns a copy of the factor scaled so its entries sum to one. If
// the total mass is zero the original factor is returned unchanged.
func (f Factor) Normalize() Factor {
	s := f.Sum()
	if s == 0 {
		return f
	}
	t := make([]float64, len(f.Table))
	for i, v := range f.Table {
		t[i] = v / s
	}
	return Factor{Vars: append([]string(nil), f.Vars...), Card: append([]int(nil), f.Card...), Table: t}
}

// Multiply returns the factor product of a and b, whose scope is the union of
// their variables. Shared variables are matched by name.
func (a Factor) Multiply(b Factor) Factor {
	// Build union scope.
	vars := append([]string(nil), a.Vars...)
	card := append([]int(nil), a.Card...)
	for i, v := range b.Vars {
		if a.varIndex(v) == -1 {
			vars = append(vars, v)
			card = append(card, b.Card[i])
		}
	}
	size := 1
	for _, c := range card {
		size *= c
	}
	res := Factor{Vars: vars, Card: card, Table: make([]float64, size)}
	assign := make(map[string]int, len(vars))
	for idx := 0; idx < size; idx++ {
		a2 := res.indexToAssignment(idx)
		for i, v := range vars {
			assign[v] = a2[i]
		}
		res.Table[idx] = a.Get(assign) * b.Get(assign)
	}
	return res
}

// Marginalize returns the factor obtained by summing out the named variable. If
// the variable is not in scope the factor is returned unchanged.
func (f Factor) Marginalize(name string) Factor {
	vi := f.varIndex(name)
	if vi == -1 {
		return f
	}
	vars := make([]string, 0, len(f.Vars)-1)
	card := make([]int, 0, len(f.Card)-1)
	for i, v := range f.Vars {
		if i == vi {
			continue
		}
		vars = append(vars, v)
		card = append(card, f.Card[i])
	}
	size := 1
	for _, c := range card {
		size *= c
	}
	res := Factor{Vars: vars, Card: card, Table: make([]float64, size)}
	assign := make(map[string]int, len(f.Vars))
	for idx := 0; idx < len(f.Table); idx++ {
		a := f.indexToAssignment(idx)
		for i, v := range f.Vars {
			assign[v] = a[i]
		}
		delete(assign, name)
		res.Table[res.assignmentToIndex(mapToAssign(res.Vars, assign))] += f.Table[idx]
	}
	return res
}

// mapToAssign extracts the per-variable assignment for vars from a map.
func mapToAssign(vars []string, assign map[string]int) []int {
	a := make([]int, len(vars))
	for i, v := range vars {
		a[i] = assign[v]
	}
	return a
}

// Reduce returns the factor obtained by fixing the variables in evidence to the
// given category indices and dropping them from scope. Variables not in scope
// are ignored.
func (f Factor) Reduce(evidence map[string]int) Factor {
	// Determine remaining variables.
	keep := make([]bool, len(f.Vars))
	vars := make([]string, 0, len(f.Vars))
	card := make([]int, 0, len(f.Card))
	for i, v := range f.Vars {
		if _, ok := evidence[v]; ok {
			keep[i] = false
		} else {
			keep[i] = true
			vars = append(vars, v)
			card = append(card, f.Card[i])
		}
	}
	size := 1
	for _, c := range card {
		size *= c
	}
	res := Factor{Vars: vars, Card: card, Table: make([]float64, size)}
	assign := make(map[string]int, len(f.Vars))
	for idx := 0; idx < len(f.Table); idx++ {
		a := f.indexToAssignment(idx)
		consistent := true
		for i, v := range f.Vars {
			if val, ok := evidence[v]; ok && val != a[i] {
				consistent = false
				break
			}
		}
		if !consistent {
			continue
		}
		for i, v := range f.Vars {
			if keep[i] {
				assign[v] = a[i]
			}
		}
		res.Table[res.assignmentToIndex(mapToAssign(res.Vars, assign))] = f.Table[idx]
	}
	return res
}

// Scope returns a sorted copy of the factor's variable names.
func (f Factor) Scope() []string {
	out := append([]string(nil), f.Vars...)
	sort.Strings(out)
	return out
}
