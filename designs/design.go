package designs

import (
	"errors"
	"sort"
)

// Design is an incidence structure: a finite set of points {0,1,...,Points-1}
// together with a list of blocks, each block being a sorted set of distinct
// point indices. It is the common representation returned by the constructions
// in this package and carries methods to recover design parameters.
type Design struct {
	Points int
	Blocks [][]int
}

// NewDesign builds a Design on the given number of points from the supplied
// blocks. Each block is copied, de-duplicated and sorted. It reports an error
// when points<=0 or a block references a point outside [0,points).
func NewDesign(points int, blocks [][]int) (*Design, error) {
	if points <= 0 {
		return nil, errors.New("designs: number of points must be positive")
	}
	out := make([][]int, 0, len(blocks))
	for _, b := range blocks {
		seen := make(map[int]bool, len(b))
		nb := make([]int, 0, len(b))
		for _, x := range b {
			if x < 0 || x >= points {
				return nil, errors.New("designs: point index out of range")
			}
			if !seen[x] {
				seen[x] = true
				nb = append(nb, x)
			}
		}
		sort.Ints(nb)
		out = append(out, nb)
	}
	return &Design{Points: points, Blocks: out}, nil
}

// NumPoints returns the number of points v of the design.
func (d *Design) NumPoints() int { return d.Points }

// NumBlocks returns the number of blocks b of the design.
func (d *Design) NumBlocks() int { return len(d.Blocks) }

// IncidenceMatrix returns the points-by-blocks 0/1 incidence matrix M, where
// M[i][j] is 1 exactly when point i lies in block j.
func (d *Design) IncidenceMatrix() [][]int {
	m := make([][]int, d.Points)
	for i := range m {
		m[i] = make([]int, len(d.Blocks))
	}
	for j, b := range d.Blocks {
		for _, x := range b {
			m[x][j] = 1
		}
	}
	return m
}

// PointDegree returns the number of blocks containing point i (its replication
// number).
func (d *Design) PointDegree(i int) int {
	c := 0
	for _, b := range d.Blocks {
		for _, x := range b {
			if x == i {
				c++
				break
			}
		}
	}
	return c
}

// ReplicationNumbers returns the replication number of every point, indexed by
// point.
func (d *Design) ReplicationNumbers() []int {
	r := make([]int, d.Points)
	for _, b := range d.Blocks {
		for _, x := range b {
			r[x]++
		}
	}
	return r
}

// IsPointRegular reports whether every point has the same replication number,
// returning that common value r and true, or (0,false) otherwise.
func (d *Design) IsPointRegular() (int, bool) {
	r := d.ReplicationNumbers()
	if len(r) == 0 {
		return 0, false
	}
	for _, x := range r {
		if x != r[0] {
			return 0, false
		}
	}
	return r[0], true
}

// BlockSizes returns the size of every block, indexed by block.
func (d *Design) BlockSizes() []int {
	s := make([]int, len(d.Blocks))
	for j, b := range d.Blocks {
		s[j] = len(b)
	}
	return s
}

// IsUniform reports whether all blocks have the same size k, returning that
// common size and true, or (0,false) otherwise.
func (d *Design) IsUniform() (int, bool) {
	if len(d.Blocks) == 0 {
		return 0, false
	}
	k := len(d.Blocks[0])
	for _, b := range d.Blocks {
		if len(b) != k {
			return 0, false
		}
	}
	return k, true
}

// Concurrence returns the number of blocks that contain both point i and point
// j. For i==j this equals the replication number of point i.
func (d *Design) Concurrence(i, j int) int {
	c := 0
	for _, b := range d.Blocks {
		hi, hj := false, false
		for _, x := range b {
			if x == i {
				hi = true
			}
			if x == j {
				hj = true
			}
		}
		if hi && hj {
			c++
		}
	}
	return c
}

// ConcurrenceMatrix returns the symmetric points-by-points matrix whose (i,j)
// entry is the number of blocks containing both points i and j; the diagonal
// holds the replication numbers.
func (d *Design) ConcurrenceMatrix() [][]int {
	m := make([][]int, d.Points)
	for i := range m {
		m[i] = make([]int, d.Points)
	}
	for _, b := range d.Blocks {
		for _, x := range b {
			for _, y := range b {
				m[x][y]++
			}
		}
	}
	return m
}

// IsPairBalanced reports whether every pair of distinct points lies in the same
// number of blocks lambda, returning lambda and true, or (0,false) otherwise. A
// design with fewer than two points is not pair balanced.
func (d *Design) IsPairBalanced() (int, bool) {
	if d.Points < 2 {
		return 0, false
	}
	m := d.ConcurrenceMatrix()
	lambda := m[0][1]
	for i := 0; i < d.Points; i++ {
		for j := i + 1; j < d.Points; j++ {
			if m[i][j] != lambda {
				return 0, false
			}
		}
	}
	return lambda, true
}

// Dual returns the dual design, in which the roles of points and blocks are
// exchanged: the dual has one point per original block, and one block per
// original point listing the original blocks through that point.
func (d *Design) Dual() *Design {
	inc := make([][]int, d.Points)
	for i := range inc {
		inc[i] = nil
	}
	for j, b := range d.Blocks {
		for _, x := range b {
			inc[x] = append(inc[x], j)
		}
	}
	nd, _ := NewDesign(len(d.Blocks), inc)
	return nd
}

// Complement returns the complementary design on the same points, whose blocks
// are the set complements of the original blocks.
func (d *Design) Complement() *Design {
	blocks := make([][]int, len(d.Blocks))
	for j, b := range d.Blocks {
		inBlock := make([]bool, d.Points)
		for _, x := range b {
			inBlock[x] = true
		}
		var nb []int
		for x := 0; x < d.Points; x++ {
			if !inBlock[x] {
				nb = append(nb, x)
			}
		}
		blocks[j] = nb
	}
	nd, _ := NewDesign(d.Points, blocks)
	return nd
}

// IsSimple reports whether the design has no repeated blocks.
func (d *Design) IsSimple() bool {
	seen := make(map[string]bool, len(d.Blocks))
	for _, b := range d.Blocks {
		key := intsKey(b)
		if seen[key] {
			return false
		}
		seen[key] = true
	}
	return true
}

// BlockContains reports whether block index j contains point i.
func (d *Design) BlockContains(j, i int) bool {
	if j < 0 || j >= len(d.Blocks) {
		return false
	}
	for _, x := range d.Blocks[j] {
		if x == i {
			return true
		}
	}
	return false
}

// Parameters attempts to recover the balanced incomplete block design
// parameters (v,b,r,k,lambda) of the design. It reports an error when the
// design is not point-regular, block-uniform and pair-balanced.
func (d *Design) Parameters() (BIBDParams, error) {
	r, ok := d.IsPointRegular()
	if !ok {
		return BIBDParams{}, errors.New("designs: design is not point-regular")
	}
	k, ok := d.IsUniform()
	if !ok {
		return BIBDParams{}, errors.New("designs: design is not block-uniform")
	}
	lambda, ok := d.IsPairBalanced()
	if !ok {
		return BIBDParams{}, errors.New("designs: design is not pair-balanced")
	}
	return BIBDParams{V: d.Points, B: len(d.Blocks), R: r, K: k, Lambda: lambda}, nil
}

// IsBIBD reports whether the design is a balanced incomplete block design and,
// when it is, returns its parameters.
func (d *Design) IsBIBD() (BIBDParams, bool) {
	p, err := d.Parameters()
	if err != nil {
		return BIBDParams{}, false
	}
	return p, true
}

func intsKey(b []int) string {
	c := append([]int(nil), b...)
	sort.Ints(c)
	// build a compact key
	buf := make([]byte, 0, len(c)*3)
	for _, x := range c {
		for x >= 0 {
			buf = append(buf, byte('0'+x%10))
			x /= 10
			if x == 0 {
				break
			}
		}
		buf = append(buf, ',')
	}
	return string(buf)
}
