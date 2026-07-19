package cellular

import (
	"errors"
	"fmt"
	"strings"
)

// Boundary selects how cells beyond the ends of a one-dimensional lattice (or
// the borders of a two-dimensional grid) are read when applying a rule.
type Boundary int

const (
	// Periodic wraps the lattice into a ring (or torus in 2-D): the cell to the
	// left of index 0 is the last cell.
	Periodic Boundary = iota
	// FixedZero treats every out-of-range cell as the quiescent state 0.
	FixedZero
	// FixedOne treats every out-of-range cell as state 1.
	FixedOne
	// Reflect mirrors the lattice at each end, duplicating the edge cell.
	Reflect
)

// String returns a human-readable name for the boundary condition.
func (b Boundary) String() string {
	switch b {
	case Periodic:
		return "periodic"
	case FixedZero:
		return "fixed-zero"
	case FixedOne:
		return "fixed-one"
	case Reflect:
		return "reflect"
	default:
		return fmt.Sprintf("Boundary(%d)", int(b))
	}
}

// Rule1D is the common interface implemented by every one-dimensional cellular
// automaton rule in the package. States() reports the number of distinct cell
// values k, Radius() the neighbourhood radius r, and Apply maps a neighbourhood
// slice of length 2r+1 (leftmost cell first) to the next value of the centre
// cell.
type Rule1D interface {
	States() int
	Radius() int
	Apply(neighbourhood []int) int
}

// ErrBadState is returned when a cell value falls outside the range 0..k-1.
var ErrBadState = errors.New("cellular: cell value out of range")

// cellAt reads the cell at (possibly out-of-range) index i under the given
// boundary condition.
func cellAt(state []int, i int, bc Boundary) int {
	n := len(state)
	if i >= 0 && i < n {
		return state[i]
	}
	switch bc {
	case Periodic:
		return state[((i%n)+n)%n]
	case FixedOne:
		return 1
	case Reflect:
		if n == 1 {
			return state[0]
		}
		m := ((i % (2 * n)) + 2*n) % (2 * n)
		if m >= n {
			m = 2*n - 1 - m
		}
		return state[m]
	default: // FixedZero
		return 0
	}
}

// Step1D advances a one-dimensional configuration by one time step under rule
// and the given boundary condition, returning a freshly allocated slice of the
// same length. The input is not modified.
func Step1D(rule Rule1D, state []int, bc Boundary) []int {
	n := len(state)
	r := rule.Radius()
	out := make([]int, n)
	nb := make([]int, 2*r+1)
	for i := 0; i < n; i++ {
		for j := -r; j <= r; j++ {
			nb[j+r] = cellAt(state, i+j, bc)
		}
		out[i] = rule.Apply(nb)
	}
	return out
}

// Evolve1D returns the spacetime diagram produced by iterating rule for steps
// time steps starting from initial. The result has steps+1 rows; row 0 is a copy
// of initial and row t is the configuration after t applications of the rule.
func Evolve1D(rule Rule1D, initial []int, steps int, bc Boundary) [][]int {
	if steps < 0 {
		steps = 0
	}
	rows := make([][]int, steps+1)
	rows[0] = append([]int(nil), initial...)
	cur := rows[0]
	for t := 1; t <= steps; t++ {
		cur = Step1D(rule, cur, bc)
		rows[t] = cur
	}
	return rows
}

// IterateState applies rule steps times and returns only the final
// configuration, allocating one working buffer instead of a full diagram.
func IterateState(rule Rule1D, initial []int, steps int, bc Boundary) []int {
	cur := append([]int(nil), initial...)
	for t := 0; t < steps; t++ {
		cur = Step1D(rule, cur, bc)
	}
	return cur
}

// Population returns the number of non-zero cells in state.
func Population(state []int) int {
	c := 0
	for _, v := range state {
		if v != 0 {
			c++
		}
	}
	return c
}

// Density returns the fraction of non-zero cells in state, or 0 for an empty
// slice.
func Density(state []int) float64 {
	if len(state) == 0 {
		return 0
	}
	return float64(Population(state)) / float64(len(state))
}

// Sum returns the arithmetic sum of the cell values in state.
func Sum(state []int) int {
	s := 0
	for _, v := range state {
		s += v
	}
	return s
}

// HammingDistance returns the number of positions at which a and b differ. It
// returns an error if the slices have different lengths.
func HammingDistance(a, b []int) (int, error) {
	if len(a) != len(b) {
		return 0, fmt.Errorf("cellular: HammingDistance length mismatch %d != %d", len(a), len(b))
	}
	d := 0
	for i := range a {
		if a[i] != b[i] {
			d++
		}
	}
	return d, nil
}

// EqualState reports whether a and b have identical length and contents.
func EqualState(a, b []int) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

// Uniform reports whether every cell of state holds the same value. The empty
// slice is considered uniform.
func Uniform(state []int) bool {
	for i := 1; i < len(state); i++ {
		if state[i] != state[0] {
			return false
		}
	}
	return true
}

// CloneState returns an independent copy of state.
func CloneState(state []int) []int {
	return append([]int(nil), state...)
}

// SingleSeedState returns a length-n configuration of zeros with a single 1 in
// the middle cell (index n/2). It returns nil for n <= 0.
func SingleSeedState(n int) []int {
	if n <= 0 {
		return nil
	}
	s := make([]int, n)
	s[n/2] = 1
	return s
}

// SingleSeedStateAt returns a length-n configuration of zeros with a single 1 at
// index pos. It returns an error if pos is out of range.
func SingleSeedStateAt(n, pos int) ([]int, error) {
	if n <= 0 {
		return nil, errors.New("cellular: SingleSeedStateAt needs n > 0")
	}
	if pos < 0 || pos >= n {
		return nil, fmt.Errorf("cellular: SingleSeedStateAt position %d out of range [0,%d)", pos, n)
	}
	s := make([]int, n)
	s[pos] = 1
	return s, nil
}

// splitmix64 is a small deterministic pseudo-random generator used to build
// reproducible initial conditions without importing math/rand.
type splitmix64 struct{ s uint64 }

func (p *splitmix64) next() uint64 {
	p.s += 0x9E3779B97F4A7C15
	z := p.s
	z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9
	z = (z ^ (z >> 27)) * 0x94D049BB133111EB
	return z ^ (z >> 31)
}

// RandomState returns a reproducible length-n configuration whose cells are
// uniformly distributed over 0..k-1, generated by a deterministic splitmix64
// stream seeded with seed. It returns an error if n < 0 or k < 1.
func RandomState(n, k int, seed uint64) ([]int, error) {
	if n < 0 {
		return nil, errors.New("cellular: RandomState needs n >= 0")
	}
	if k < 1 {
		return nil, errors.New("cellular: RandomState needs k >= 1")
	}
	g := splitmix64{s: seed}
	s := make([]int, n)
	for i := range s {
		s[i] = int(g.next() % uint64(k))
	}
	return s, nil
}

// RandomBinaryState returns a reproducible length-n binary configuration, a
// convenience wrapper around RandomState with k = 2. It returns nil for n <= 0.
func RandomBinaryState(n int, seed uint64) []int {
	if n <= 0 {
		return nil
	}
	s, _ := RandomState(n, 2, seed)
	return s
}

// StateFromString parses a string into a binary configuration, mapping each rune
// equal to on to 1 and every other rune to 0.
func StateFromString(s string, on rune) []int {
	out := make([]int, 0, len(s))
	for _, r := range s {
		if r == on {
			out = append(out, 1)
		} else {
			out = append(out, 0)
		}
	}
	return out
}

// StateToString renders a configuration as a string using off for 0 and on for
// every non-zero value.
func StateToString(state []int, off, on rune) string {
	var b strings.Builder
	for _, v := range state {
		if v == 0 {
			b.WriteRune(off)
		} else {
			b.WriteRune(on)
		}
	}
	return b.String()
}

// ValidateState reports an error if any cell of state lies outside 0..k-1.
func ValidateState(state []int, k int) error {
	for i, v := range state {
		if v < 0 || v >= k {
			return fmt.Errorf("cellular: %w: state[%d]=%d not in [0,%d)", ErrBadState, i, v, k)
		}
	}
	return nil
}
