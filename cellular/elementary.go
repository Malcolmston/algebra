package cellular

import (
	"fmt"
	"math/bits"
	"strings"
)

// ElementaryRule is an elementary (k=2, r=1) one-dimensional cellular-automaton
// rule identified by its Wolfram number 0..255. Bit p of the number (for a
// neighbourhood value p = 4*left + 2*centre + right) gives the next value of the
// centre cell.
type ElementaryRule uint8

// NewElementaryRule returns the ElementaryRule with Wolfram number n. It returns
// an error if n is not in the range 0..255.
func NewElementaryRule(n int) (ElementaryRule, error) {
	if n < 0 || n > 255 {
		return 0, fmt.Errorf("cellular: elementary rule %d out of range [0,255]", n)
	}
	return ElementaryRule(n), nil
}

// MustElementaryRule is like NewElementaryRule but panics on an out-of-range
// number. It is intended for package-level rule constants and tests.
func MustElementaryRule(n int) ElementaryRule {
	r, err := NewElementaryRule(n)
	if err != nil {
		panic(err)
	}
	return r
}

// Number returns the Wolfram rule number as an int.
func (r ElementaryRule) Number() int { return int(r) }

// States reports the number of cell states, which is always 2 for an elementary
// rule.
func (r ElementaryRule) States() int { return 2 }

// Radius reports the neighbourhood radius, which is always 1 for an elementary
// rule.
func (r ElementaryRule) Radius() int { return 1 }

// Bit returns the output bit (0 or 1) for the neighbourhood value p in 0..7,
// where p = 4*left + 2*centre + right.
func (r ElementaryRule) Bit(p int) int {
	return int(r>>uint(p&7)) & 1
}

// ApplyLCR returns the next centre value for the ordered neighbourhood
// (left, centre, right), each 0 or 1.
func (r ElementaryRule) ApplyLCR(left, centre, right int) int {
	p := (left&1)<<2 | (centre&1)<<1 | (right & 1)
	return r.Bit(p)
}

// Apply implements the Rule1D interface for a length-3 neighbourhood slice
// (left, centre, right).
func (r ElementaryRule) Apply(neighbourhood []int) int {
	return r.ApplyLCR(neighbourhood[0], neighbourhood[1], neighbourhood[2])
}

// Table returns the 8-entry output table indexed by neighbourhood value 0..7.
func (r ElementaryRule) Table() [8]int {
	var t [8]int
	for p := 0; p < 8; p++ {
		t[p] = r.Bit(p)
	}
	return t
}

// String returns the rule in the form "rule 110".
func (r ElementaryRule) String() string {
	return fmt.Sprintf("rule %d", int(r))
}

// ActiveCount returns the number of neighbourhoods that map to 1 (the population
// count of the rule number).
func (r ElementaryRule) ActiveCount() int {
	return bits.OnesCount8(uint8(r))
}

// LambdaParameter returns Langton's lambda for the rule: the fraction of
// neighbourhoods mapping to a non-quiescent state, i.e. ActiveCount/8.
func (r ElementaryRule) LambdaParameter() float64 {
	return float64(r.ActiveCount()) / 8.0
}

// IsQuiescent reports whether the all-zero neighbourhood maps to 0, so that a
// blank background stays blank.
func (r ElementaryRule) IsQuiescent() bool {
	return r.Bit(0) == 0
}

// IsSymmetric reports whether the rule is invariant under left-right reflection,
// i.e. it equals its own mirror rule.
func (r ElementaryRule) IsSymmetric() bool {
	return r == r.MirrorRule()
}

// IsLegal reports whether the rule is "legal" in Wolfram's original sense:
// quiescent (blank stays blank) and left-right symmetric.
func (r ElementaryRule) IsLegal() bool {
	return r.IsQuiescent() && r.IsSymmetric()
}

// MirrorRule returns the rule obtained by reflecting every neighbourhood
// left-to-right. Applying MirrorRule to a configuration is equivalent to
// applying r to the spatially reversed configuration.
func (r ElementaryRule) MirrorRule() ElementaryRule {
	var out ElementaryRule
	for p := 0; p < 8; p++ {
		l := (p >> 2) & 1
		c := (p >> 1) & 1
		rt := p & 1
		q := (rt << 2) | (c << 1) | l
		if r.Bit(p) == 1 {
			out |= 1 << uint(q)
		}
	}
	return out
}

// ComplementRule returns the rule obtained by exchanging the roles of 0 and 1 in
// both inputs and outputs. Its evolution is the bitwise complement of r's
// evolution from the complemented initial condition.
func (r ElementaryRule) ComplementRule() ElementaryRule {
	var out ElementaryRule
	for p := 0; p < 8; p++ {
		q := (^p) & 7
		if r.Bit(p) == 0 {
			out |= 1 << uint(q)
		}
	}
	return out
}

// MirrorComplementRule returns the composition of MirrorRule and ComplementRule.
func (r ElementaryRule) MirrorComplementRule() ElementaryRule {
	return r.MirrorRule().ComplementRule()
}

// Equivalents returns the sorted set of rules in the same equivalence class as r
// under left-right reflection and 0/1 complementation (1, 2 or 4 members).
func (r ElementaryRule) Equivalents() []ElementaryRule {
	set := map[ElementaryRule]bool{
		r:                        true,
		r.MirrorRule():           true,
		r.ComplementRule():       true,
		r.MirrorComplementRule(): true,
	}
	out := make([]ElementaryRule, 0, len(set))
	for k := range set {
		out = append(out, k)
	}
	// simple insertion sort to avoid importing sort for such a tiny slice
	for i := 1; i < len(out); i++ {
		for j := i; j > 0 && out[j-1] > out[j]; j-- {
			out[j-1], out[j] = out[j], out[j-1]
		}
	}
	return out
}

// MinimalEquivalent returns the smallest rule number in r's equivalence class,
// the canonical representative used to enumerate the 88 inequivalent elementary
// rules.
func (r ElementaryRule) MinimalEquivalent() ElementaryRule {
	eq := r.Equivalents()
	return eq[0]
}

// EquivalenceClassSize returns the number of distinct rules in r's equivalence
// class (1, 2 or 4).
func (r ElementaryRule) EquivalenceClassSize() int {
	return len(r.Equivalents())
}

// IsCanonical reports whether r is the minimal representative of its equivalence
// class.
func (r ElementaryRule) IsCanonical() bool {
	return r == r.MinimalEquivalent()
}

// IsAdditive reports whether the rule is additive (linear over GF(2)): its
// output equals a fixed XOR combination of the three neighbourhood cells. The
// additive elementary rules are 0, 60, 90, 102, 150, 170, 204 and 240.
func (r ElementaryRule) IsAdditive() bool {
	// Find the linear coefficients from the images of the unit neighbourhoods.
	// f(l,c,r) = a*l XOR b*c XOR d*r XOR e, with e = f(0,0,0).
	e := r.ApplyLCR(0, 0, 0)
	if e != 0 {
		// A strictly additive rule maps the blank neighbourhood to 0; a non-zero
		// constant term makes the rule affine but not additive.
		return false
	}
	a := r.ApplyLCR(1, 0, 0) ^ e
	b := r.ApplyLCR(0, 1, 0) ^ e
	d := r.ApplyLCR(0, 0, 1) ^ e
	for p := 0; p < 8; p++ {
		l := (p >> 2) & 1
		c := (p >> 1) & 1
		rt := p & 1
		want := (a*l ^ b*c ^ d*rt ^ e) & 1
		if r.Bit(p) != want {
			return false
		}
	}
	return true
}

// IsOuterTotalistic reports whether the rule depends only on the centre cell and
// the sum of its two neighbours (left+right), rather than on their individual
// values.
func (r ElementaryRule) IsOuterTotalistic() bool {
	for c := 0; c < 2; c++ {
		// neighbours summing to 1 can be (0,1) or (1,0); outputs must agree.
		o01 := r.ApplyLCR(0, c, 1)
		o10 := r.ApplyLCR(1, c, 0)
		if o01 != o10 {
			return false
		}
	}
	return true
}

// IsTotalistic reports whether the rule depends only on the total number of 1s
// in the whole three-cell neighbourhood.
func (r ElementaryRule) IsTotalistic() bool {
	byCount := map[int]int{}
	for p := 0; p < 8; p++ {
		cnt := bits.OnesCount8(uint8(p) & 7)
		out := r.Bit(p)
		if prev, ok := byCount[cnt]; ok {
			if prev != out {
				return false
			}
		} else {
			byCount[cnt] = out
		}
	}
	return true
}

// ElementaryStep advances a binary configuration by one time step under the
// elementary rule and boundary condition, a thin wrapper over Step1D.
func ElementaryStep(r ElementaryRule, state []int, bc Boundary) []int {
	return Step1D(r, state, bc)
}

// ElementaryEvolve returns the steps+1 row spacetime diagram of the elementary
// rule starting from initial.
func ElementaryEvolve(r ElementaryRule, initial []int, steps int, bc Boundary) [][]int {
	return Evolve1D(r, initial, steps, bc)
}

// ElementaryEvolveSeed evolves the rule for steps steps from a single central
// seed on a background of width cells, using a fixed-zero boundary so the light
// cone is fully contained. It returns the spacetime diagram.
func ElementaryEvolveSeed(r ElementaryRule, width, steps int) [][]int {
	return Evolve1D(r, SingleSeedState(width), steps, FixedZero)
}

// AllElementaryRules returns all 256 elementary rules in numeric order.
func AllElementaryRules() []ElementaryRule {
	out := make([]ElementaryRule, 256)
	for i := range out {
		out[i] = ElementaryRule(i)
	}
	return out
}

// InequivalentElementaryRules returns the 88 canonical elementary rules, one per
// equivalence class under reflection and complementation, in numeric order.
func InequivalentElementaryRules() []ElementaryRule {
	var out []ElementaryRule
	for i := 0; i < 256; i++ {
		r := ElementaryRule(i)
		if r.IsCanonical() {
			out = append(out, r)
		}
	}
	return out
}

// RuleToBits returns the 8 output bits of an elementary rule, indexed by
// neighbourhood value 0..7 (identical to Table but returned as a slice).
func RuleToBits(r ElementaryRule) []int {
	t := r.Table()
	return t[:]
}

// RuleFromBits builds an elementary rule from an 8-entry output table indexed by
// neighbourhood value 0..7. It returns an error if bits has the wrong length or
// contains values other than 0 and 1.
func RuleFromBits(bitsTable []int) (ElementaryRule, error) {
	if len(bitsTable) != 8 {
		return 0, fmt.Errorf("cellular: RuleFromBits needs 8 entries, got %d", len(bitsTable))
	}
	var r ElementaryRule
	for p, v := range bitsTable {
		switch v {
		case 0:
		case 1:
			r |= 1 << uint(p)
		default:
			return 0, fmt.Errorf("cellular: RuleFromBits entry %d is %d, want 0 or 1", p, v)
		}
	}
	return r, nil
}

// RuleTableString returns a multi-line human-readable transition table for an
// elementary rule, listing each neighbourhood pattern and its output.
func RuleTableString(r ElementaryRule) string {
	var b strings.Builder
	for p := 7; p >= 0; p-- {
		l := (p >> 2) & 1
		c := (p >> 1) & 1
		rt := p & 1
		fmt.Fprintf(&b, "%d%d%d -> %d\n", l, c, rt, r.Bit(p))
	}
	return b.String()
}
