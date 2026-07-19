package cellular

import "fmt"

// SecondOrderStep computes one forward step of the second-order (Fredkin)
// construction of an elementary rule: the next configuration is the rule applied
// to the current configuration, exclusive-ORed cell by cell with the previous
// configuration. The construction is exactly reversible. prev and cur must have
// equal length.
func SecondOrderStep(r ElementaryRule, prev, cur []int, bc Boundary) []int {
	f := Step1D(r, cur, bc)
	out := make([]int, len(cur))
	for i := range cur {
		out[i] = f[i] ^ prev[i]
	}
	return out
}

// SecondOrderReverse computes one backward step of the second-order
// construction: given the current and next configurations it recovers the
// previous one. Because XOR is its own inverse the formula is identical to the
// forward step.
func SecondOrderReverse(r ElementaryRule, next, cur []int, bc Boundary) []int {
	return SecondOrderStep(r, next, cur, bc)
}

// SecondOrderCA holds the state of a running second-order elementary automaton:
// the underlying rule, the previous and current generations and the boundary
// condition. It can be stepped forward and backward without loss of
// information.
type SecondOrderCA struct {
	Rule ElementaryRule
	Prev []int
	Cur  []int
	BC   Boundary
}

// NewSecondOrderCA creates a second-order automaton from a rule and two initial
// generations. It returns an error if the generations differ in length.
func NewSecondOrderCA(r ElementaryRule, prev, cur []int, bc Boundary) (*SecondOrderCA, error) {
	if len(prev) != len(cur) {
		return nil, fmt.Errorf("cellular: NewSecondOrderCA generations differ in length %d != %d", len(prev), len(cur))
	}
	return &SecondOrderCA{
		Rule: r,
		Prev: CloneState(prev),
		Cur:  CloneState(cur),
		BC:   bc,
	}, nil
}

// Step advances the automaton one generation forward and returns the new current
// configuration.
func (s *SecondOrderCA) Step() []int {
	next := SecondOrderStep(s.Rule, s.Prev, s.Cur, s.BC)
	s.Prev, s.Cur = s.Cur, next
	return s.Cur
}

// StepBack reverses the automaton one generation and returns the recovered
// configuration, exactly undoing a call to Step.
func (s *SecondOrderCA) StepBack() []int {
	prev := SecondOrderReverse(s.Rule, s.Cur, s.Prev, s.BC)
	s.Cur, s.Prev = s.Prev, prev
	return s.Cur
}

// State returns a copy of the current configuration.
func (s *SecondOrderCA) State() []int { return CloneState(s.Cur) }

// Clone returns an independent copy of the automaton.
func (s *SecondOrderCA) Clone() *SecondOrderCA {
	return &SecondOrderCA{
		Rule: s.Rule,
		Prev: CloneState(s.Prev),
		Cur:  CloneState(s.Cur),
		BC:   s.BC,
	}
}

// SecondOrderEvolve runs the second-order construction of a rule for steps
// generations from the two initial configurations and returns the spacetime
// diagram, which has steps+2 rows (the two seeds followed by each new
// generation).
func SecondOrderEvolve(r ElementaryRule, prev, cur []int, steps int, bc Boundary) [][]int {
	rows := make([][]int, 0, steps+2)
	rows = append(rows, CloneState(prev), CloneState(cur))
	p, c := CloneState(prev), CloneState(cur)
	for t := 0; t < steps; t++ {
		next := SecondOrderStep(r, p, c, bc)
		rows = append(rows, next)
		p, c = c, next
	}
	return rows
}

// MargolusRule is a reversible two-dimensional block rule acting on 2x2
// Margolus-neighbourhood blocks of a binary grid. Entry i of the array gives the
// output block for input block i, where a block encodes its four cells as bits
// nw=8, ne=4, sw=2, se=1. The rule is reversible if and only if the array is a
// permutation of 0..15, which IsPermutation verifies.
type MargolusRule [16]int

// IsPermutation reports whether the rule is a bijection of the 16 block states
// and therefore reversible.
func (m MargolusRule) IsPermutation() bool {
	var seen [16]bool
	for _, v := range m {
		if v < 0 || v > 15 || seen[v] {
			return false
		}
		seen[v] = true
	}
	return true
}

// Inverse returns the inverse Margolus rule. It returns an error if the rule is
// not a permutation.
func (m MargolusRule) Inverse() (MargolusRule, error) {
	if !m.IsPermutation() {
		return MargolusRule{}, fmt.Errorf("cellular: MargolusRule is not a permutation")
	}
	var inv MargolusRule
	for i, v := range m {
		inv[v] = i
	}
	return inv, nil
}

// Rotate180Rule returns the Margolus rule that rotates every 2x2 block by 180
// degrees (nw<->se, ne<->sw). It is an involutive, reversible block rule.
func Rotate180Rule() MargolusRule {
	var m MargolusRule
	for i := 0; i < 16; i++ {
		nw := (i >> 3) & 1
		ne := (i >> 2) & 1
		sw := (i >> 1) & 1
		se := i & 1
		// after 180 rotation
		m[i] = (se << 3) | (sw << 2) | (ne << 1) | nw
	}
	return m
}

// CrittersRule returns an involutive, Critters-style reversible block rule:
// blocks with exactly two live cells are left unchanged, and every other block
// is complemented. Because complementation maps a population count p to 4-p, it
// preserves the "population equals two" partition, so the map is a bijection.
func CrittersRule() MargolusRule {
	var m MargolusRule
	for i := 0; i < 16; i++ {
		pop := bitCount4(i)
		if pop == 2 {
			m[i] = i
		} else {
			m[i] = 15 - i
		}
	}
	return m
}

// bitCount4 returns the number of set bits among the low four bits of v.
func bitCount4(v int) int {
	c := 0
	for j := 0; j < 4; j++ {
		c += (v >> uint(j)) & 1
	}
	return c
}

// MargolusStep applies a Margolus block rule to a binary grid for one time step.
// When even is true the 2x2 blocks are aligned to the origin; when false they are
// offset by (1,1), realising the alternating Margolus partition. The grid
// dimensions must both be even. The input grid is not modified.
func MargolusStep(m MargolusRule, g *Grid, even bool) (*Grid, error) {
	if g.rows%2 != 0 || g.cols%2 != 0 {
		return nil, fmt.Errorf("cellular: MargolusStep needs even dimensions, got %dx%d", g.rows, g.cols)
	}
	out := g.Clone()
	off := 0
	if !even {
		off = 1
	}
	for br := 0; br < g.rows; br += 2 {
		for bc := 0; bc < g.cols; bc += 2 {
			r0 := (br + off) % g.rows
			c0 := (bc + off) % g.cols
			r1 := (r0 + 1) % g.rows
			c1 := (c0 + 1) % g.cols
			nw := g.at(r0, c0)
			ne := g.at(r0, c1)
			sw := g.at(r1, c0)
			se := g.at(r1, c1)
			in := (nw << 3) | (ne << 2) | (sw << 1) | se
			o := m[in]
			out.set(r0, c0, (o>>3)&1)
			out.set(r0, c1, (o>>2)&1)
			out.set(r1, c0, (o>>1)&1)
			out.set(r1, c1, o&1)
		}
	}
	return out, nil
}

// MargolusEvolve runs a Margolus block rule for steps time steps, alternating the
// block partition each step (starting with the even partition), and returns the
// sequence of steps+1 grids beginning with a copy of g.
func MargolusEvolve(m MargolusRule, g *Grid, steps int) ([]*Grid, error) {
	frames := make([]*Grid, 0, steps+1)
	frames = append(frames, g.Clone())
	cur := g.Clone()
	for t := 0; t < steps; t++ {
		next, err := MargolusStep(m, cur, t%2 == 0)
		if err != nil {
			return nil, err
		}
		frames = append(frames, next)
		cur = next
	}
	return frames, nil
}
