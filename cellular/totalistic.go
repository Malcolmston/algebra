package cellular

import (
	"fmt"
	"strings"
)

// TotalisticRule is a k-state, radius-r totalistic one-dimensional rule: the new
// value of a cell depends only on the sum of the values in its neighbourhood
// (including the centre). The rule is identified by a base-k code whose digit s
// gives the output for neighbourhood sum s.
type TotalisticRule struct {
	k     int
	r     int
	code  int
	table []int // index = neighbourhood sum, value = next state
}

// NewTotalisticRule builds a k-state, radius-r totalistic rule from its code. It
// returns an error for k < 2, r < 0 or a code outside the valid range.
func NewTotalisticRule(k, r, code int) (*TotalisticRule, error) {
	if k < 2 {
		return nil, fmt.Errorf("cellular: NewTotalisticRule needs k >= 2, got %d", k)
	}
	if r < 0 {
		return nil, fmt.Errorf("cellular: NewTotalisticRule needs r >= 0, got %d", r)
	}
	entries := (2*r+1)*(k-1) + 1
	table, err := DecodeBaseK(code, entries, k)
	if err != nil {
		return nil, fmt.Errorf("cellular: totalistic code %d invalid for k=%d r=%d: %w", code, k, r, err)
	}
	return &TotalisticRule{k: k, r: r, code: code, table: table}, nil
}

// States reports the number of cell states k.
func (t *TotalisticRule) States() int { return t.k }

// Radius reports the neighbourhood radius r.
func (t *TotalisticRule) Radius() int { return t.r }

// Code returns the base-k code identifying the rule.
func (t *TotalisticRule) Code() int { return t.code }

// MaxSum returns the largest possible neighbourhood sum, (2r+1)(k-1).
func (t *TotalisticRule) MaxSum() int { return (2*t.r + 1) * (t.k - 1) }

// Table returns a copy of the output table indexed by neighbourhood sum.
func (t *TotalisticRule) Table() []int { return append([]int(nil), t.table...) }

// Apply implements Rule1D by summing the neighbourhood and looking up the sum in
// the table.
func (t *TotalisticRule) Apply(neighbourhood []int) int {
	s := 0
	for _, v := range neighbourhood {
		s += v
	}
	if s < 0 || s >= len(t.table) {
		return 0
	}
	return t.table[s]
}

// String returns a description such as "totalistic k=2 r=1 code=90".
func (t *TotalisticRule) String() string {
	return fmt.Sprintf("totalistic k=%d r=%d code=%d", t.k, t.r, t.code)
}

// OuterTotalisticRule is a k-state, radius-r outer-totalistic rule: the new value
// of a cell depends on the centre cell together with the sum of the surrounding
// cells (excluding the centre). The rule is identified by a base-k code indexed
// by centre*(outerMax+1) + outerSum.
type OuterTotalisticRule struct {
	k        int
	r        int
	code     int
	outerMax int
	table    []int
}

// NewOuterTotalisticRule builds a k-state, radius-r outer-totalistic rule from
// its code. It returns an error for invalid parameters or an out-of-range code.
func NewOuterTotalisticRule(k, r, code int) (*OuterTotalisticRule, error) {
	if k < 2 {
		return nil, fmt.Errorf("cellular: NewOuterTotalisticRule needs k >= 2, got %d", k)
	}
	if r < 1 {
		return nil, fmt.Errorf("cellular: NewOuterTotalisticRule needs r >= 1, got %d", r)
	}
	outerMax := 2 * r * (k - 1)
	entries := k * (outerMax + 1)
	table, err := DecodeBaseK(code, entries, k)
	if err != nil {
		return nil, fmt.Errorf("cellular: outer-totalistic code %d invalid for k=%d r=%d: %w", code, k, r, err)
	}
	return &OuterTotalisticRule{k: k, r: r, code: code, outerMax: outerMax, table: table}, nil
}

// States reports the number of cell states k.
func (o *OuterTotalisticRule) States() int { return o.k }

// Radius reports the neighbourhood radius r.
func (o *OuterTotalisticRule) Radius() int { return o.r }

// Code returns the base-k code identifying the rule.
func (o *OuterTotalisticRule) Code() int { return o.code }

// OuterMax returns the largest possible outer-neighbour sum, 2r(k-1).
func (o *OuterTotalisticRule) OuterMax() int { return o.outerMax }

// Table returns a copy of the output table.
func (o *OuterTotalisticRule) Table() []int { return append([]int(nil), o.table...) }

// Apply implements Rule1D. neighbourhood must have length 2r+1; the centre cell
// is the middle element.
func (o *OuterTotalisticRule) Apply(neighbourhood []int) int {
	mid := len(neighbourhood) / 2
	centre := neighbourhood[mid]
	outer := 0
	for i, v := range neighbourhood {
		if i != mid {
			outer += v
		}
	}
	idx := centre*(o.outerMax+1) + outer
	if idx < 0 || idx >= len(o.table) {
		return 0
	}
	return o.table[idx]
}

// String returns a description such as "outer-totalistic k=2 r=1 code=6".
func (o *OuterTotalisticRule) String() string {
	return fmt.Sprintf("outer-totalistic k=%d r=%d code=%d", o.k, o.r, o.code)
}

// GeneralRule is a fully general k-state, radius-r rule specified by a complete
// lookup table over all k^(2r+1) neighbourhoods. Neighbourhoods are indexed in
// base k with the leftmost cell as the most significant digit.
type GeneralRule struct {
	k     int
	r     int
	table []int
}

// NewGeneralRule builds a general rule from an explicit table whose length must
// equal k^(2r+1). Table entry i is the output for the neighbourhood whose base-k
// encoding (leftmost cell most significant) is i.
func NewGeneralRule(k, r int, table []int) (*GeneralRule, error) {
	if k < 2 {
		return nil, fmt.Errorf("cellular: NewGeneralRule needs k >= 2, got %d", k)
	}
	if r < 0 {
		return nil, fmt.Errorf("cellular: NewGeneralRule needs r >= 0, got %d", r)
	}
	want := 1
	for i := 0; i < 2*r+1; i++ {
		want *= k
	}
	if len(table) != want {
		return nil, fmt.Errorf("cellular: NewGeneralRule table length %d, want %d", len(table), want)
	}
	for i, v := range table {
		if v < 0 || v >= k {
			return nil, fmt.Errorf("cellular: NewGeneralRule table[%d]=%d out of range [0,%d)", i, v, k)
		}
	}
	return &GeneralRule{k: k, r: r, table: append([]int(nil), table...)}, nil
}

// GeneralRuleFromCode builds a general rule whose table is the base-k
// representation of code, least significant digit indexing neighbourhood 0.
func GeneralRuleFromCode(k, r, code int) (*GeneralRule, error) {
	if k < 2 || r < 0 {
		return nil, fmt.Errorf("cellular: GeneralRuleFromCode needs k >= 2 and r >= 0")
	}
	entries := 1
	for i := 0; i < 2*r+1; i++ {
		entries *= k
	}
	table := make([]int, entries)
	c := code
	for i := 0; i < entries; i++ {
		table[i] = c % k
		c /= k
	}
	if c != 0 {
		return nil, fmt.Errorf("cellular: GeneralRuleFromCode code %d too large for k=%d r=%d", code, k, r)
	}
	return &GeneralRule{k: k, r: r, table: table}, nil
}

// States reports the number of cell states k.
func (g *GeneralRule) States() int { return g.k }

// Radius reports the neighbourhood radius r.
func (g *GeneralRule) Radius() int { return g.r }

// Table returns a copy of the full lookup table.
func (g *GeneralRule) Table() []int { return append([]int(nil), g.table...) }

// Index returns the base-k index of a neighbourhood slice (leftmost cell most
// significant), the position used to look it up in the table.
func (g *GeneralRule) Index(neighbourhood []int) int {
	idx := 0
	for _, v := range neighbourhood {
		idx = idx*g.k + v
	}
	return idx
}

// Apply implements Rule1D by indexing the full lookup table.
func (g *GeneralRule) Apply(neighbourhood []int) int {
	idx := g.Index(neighbourhood)
	if idx < 0 || idx >= len(g.table) {
		return 0
	}
	return g.table[idx]
}

// String returns a description such as "general k=2 r=1".
func (g *GeneralRule) String() string {
	return fmt.Sprintf("general k=%d r=%d", g.k, g.r)
}

// ElementaryAsGeneral converts an elementary rule to the equivalent GeneralRule
// with k = 2 and r = 1, useful for exercising the general machinery.
func ElementaryAsGeneral(r ElementaryRule) *GeneralRule {
	t := r.Table()
	table := make([]int, 8)
	// GeneralRule indexes leftmost cell most significant, matching the
	// elementary p = 4*left+2*centre+right convention exactly.
	copy(table, t[:])
	g, _ := NewGeneralRule(2, 1, table)
	return g
}

// LangtonLambda returns Langton's lambda parameter for a general rule: the
// fraction of neighbourhoods that map to a non-quiescent state, where state 0 is
// taken to be quiescent.
func LangtonLambda(rule *GeneralRule) float64 {
	if len(rule.table) == 0 {
		return 0
	}
	active := 0
	for _, v := range rule.table {
		if v != 0 {
			active++
		}
	}
	return float64(active) / float64(len(rule.table))
}

// TotalisticRuleString renders the sum-to-output table of a totalistic rule as a
// single line such as "0->0 1->1 2->0 3->1".
func TotalisticRuleString(t *TotalisticRule) string {
	var b strings.Builder
	for s, out := range t.table {
		if s > 0 {
			b.WriteByte(' ')
		}
		fmt.Fprintf(&b, "%d->%d", s, out)
	}
	return b.String()
}
