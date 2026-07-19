package designs

import "errors"

// IsParallelClass reports whether the given block indices of the design form a
// parallel class: the referenced blocks are pairwise disjoint and together
// cover every point of the design exactly once.
func IsParallelClass(d *Design, blockIdx []int) bool {
	seen := make([]bool, d.Points)
	count := 0
	for _, j := range blockIdx {
		if j < 0 || j >= len(d.Blocks) {
			return false
		}
		for _, x := range d.Blocks[j] {
			if seen[x] {
				return false
			}
			seen[x] = true
			count++
		}
	}
	return count == d.Points
}

// Resolution pairs a design with a partition of its blocks into parallel
// classes. A design admitting such a partition is called resolvable.
type Resolution struct {
	Design  *Design
	Classes [][]int // each entry is a list of block indices
}

// NewResolution validates that the given classes partition all blocks of the
// design into parallel classes and returns the resolution. It reports an error
// otherwise.
func NewResolution(d *Design, classes [][]int) (*Resolution, error) {
	usedBlock := make([]bool, len(d.Blocks))
	total := 0
	for _, c := range classes {
		if !IsParallelClass(d, c) {
			return nil, errors.New("designs: class is not a parallel class")
		}
		for _, j := range c {
			if usedBlock[j] {
				return nil, errors.New("designs: block appears in more than one class")
			}
			usedBlock[j] = true
			total++
		}
	}
	if total != len(d.Blocks) {
		return nil, errors.New("designs: classes do not cover every block")
	}
	return &Resolution{Design: d, Classes: classes}, nil
}

// NumClasses returns the number of parallel classes in the resolution.
func (r *Resolution) NumClasses() int { return len(r.Classes) }

// Class returns the block indices making up parallel class i.
func (r *Resolution) Class(i int) []int { return append([]int(nil), r.Classes[i]...) }

// IsValid reports whether the resolution is well formed: every class is a
// parallel class and together they partition the blocks.
func (r *Resolution) IsValid() bool {
	_, err := NewResolution(r.Design, r.Classes)
	return err == nil
}

// AffineResolution returns the affine plane AG(2,q) as a resolvable design
// together with its natural resolution into q+1 parallel classes (the parallel
// classes of lines). It reports an error when q is not a prime power.
func AffineResolution(q int) (*Resolution, error) {
	ap, err := NewAffinePlane(q)
	if err != nil {
		return nil, err
	}
	return NewResolution(ap.IncidenceDesign(), ap.ParallelClasses())
}

// IsResolvableDesign reports whether the design admits the supplied resolution
// into parallel classes.
func IsResolvableDesign(d *Design, classes [][]int) bool {
	_, err := NewResolution(d, classes)
	return err == nil
}

// OneFactorization returns a one-factorization of the complete graph K_m on m
// vertices (m even): a partition of its edges into m-1 perfect matchings
// (one-factors), computed by the classical circle (round-robin) method. Each
// factor is a slice of {u,v} vertex pairs. It reports an error when m is odd or
// less than 2.
func OneFactorization(m int) ([][][2]int, error) {
	if m < 2 || m%2 != 0 {
		return nil, errors.New("designs: one-factorization needs an even number of vertices")
	}
	fixed := m - 1
	rounds := m - 1
	half := m / 2
	out := make([][][2]int, 0, rounds)
	for r := 0; r < rounds; r++ {
		var factor [][2]int
		factor = append(factor, edgeSorted(fixed, r%(m-1)))
		for i := 1; i < half; i++ {
			a := (r + i) % (m - 1)
			b := ((r-i)%(m-1) + (m - 1)) % (m - 1)
			factor = append(factor, edgeSorted(a, b))
		}
		out = append(out, factor)
	}
	return out, nil
}

func edgeSorted(a, b int) [2]int {
	if a < b {
		return [2]int{a, b}
	}
	return [2]int{b, a}
}

// IsOneFactorization reports whether the given factors form a one-factorization
// of K_m: every factor is a perfect matching and every edge of K_m appears in
// exactly one factor.
func IsOneFactorization(m int, factors [][][2]int) bool {
	if m < 2 || m%2 != 0 {
		return false
	}
	if len(factors) != m-1 {
		return false
	}
	edgeSeen := make(map[[2]int]bool)
	for _, f := range factors {
		covered := make([]bool, m)
		if len(f) != m/2 {
			return false
		}
		for _, e := range f {
			u, v := e[0], e[1]
			if u < 0 || v < 0 || u >= m || v >= m || u == v {
				return false
			}
			if covered[u] || covered[v] {
				return false
			}
			covered[u] = true
			covered[v] = true
			key := edgeSorted(u, v)
			if edgeSeen[key] {
				return false
			}
			edgeSeen[key] = true
		}
	}
	return len(edgeSeen) == m*(m-1)/2
}

// RoundRobinSchedule returns a round-robin tournament schedule for the given
// number of teams: a list of rounds, each round a set of {home,away} pairings,
// such that every pair of teams meets exactly once. When the number of teams is
// odd a bye is handled internally and simply omitted from the pairings. It
// reports an error when teams<2.
func RoundRobinSchedule(teams int) ([][][2]int, error) {
	if teams < 2 {
		return nil, errors.New("designs: need at least two teams")
	}
	m := teams
	bye := -1
	if m%2 == 1 {
		bye = m
		m = m + 1
	}
	factors, err := OneFactorization(m)
	if err != nil {
		return nil, err
	}
	out := make([][][2]int, 0, len(factors))
	for _, f := range factors {
		var round [][2]int
		for _, e := range f {
			if e[0] == bye || e[1] == bye {
				continue
			}
			round = append(round, e)
		}
		out = append(out, round)
	}
	return out, nil
}

// NearOneFactorization returns a near-one-factorization of the complete graph
// K_m on an odd number m of vertices: m near-perfect matchings, each missing a
// distinct vertex, together covering every edge exactly once. It reports an
// error when m is even or less than 1.
func NearOneFactorization(m int) ([][][2]int, error) {
	if m < 1 || m%2 == 0 {
		return nil, errors.New("designs: near-one-factorization needs an odd number of vertices")
	}
	out := make([][][2]int, 0, m)
	for r := 0; r < m; r++ {
		var factor [][2]int
		for i := 1; i <= (m-1)/2; i++ {
			a := (r + i) % m
			b := ((r-i)%m + m) % m
			factor = append(factor, edgeSorted(a, b))
		}
		out = append(out, factor)
	}
	return out, nil
}
