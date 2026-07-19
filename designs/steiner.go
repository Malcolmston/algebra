package designs

import (
	"errors"
	"math/rand"
	"sort"
)

// SteinerTripleParameters reports the number of blocks v(v-1)/6 of a Steiner
// triple system on v points. It reports an error unless v is congruent to 1 or
// 3 modulo 6.
func SteinerTripleParameters(v int) (blocks int, err error) {
	if v%6 != 1 && v%6 != 3 {
		return 0, errors.New("designs: STS requires v = 1 or 3 (mod 6)")
	}
	return v * (v - 1) / 6, nil
}

// SteinerTripleAdmissible reports whether v is an admissible order for a Steiner
// triple system, i.e. v is congruent to 1 or 3 modulo 6 (or v is 0 or 1).
func SteinerTripleAdmissible(v int) bool {
	return v == 0 || v == 1 || v%6 == 1 || v%6 == 3
}

// BoseSteinerTripleSystem constructs a Steiner triple system on n = 6t+3 points
// by the Bose construction, using the idempotent commutative quasigroup
// a*b = ((m+1)/2)(a+b) mod m on Z_m with m = 2t+1. It reports an error when n is
// not congruent to 3 modulo 6.
func BoseSteinerTripleSystem(n int) (*Design, error) {
	if n%6 != 3 {
		return nil, errors.New("designs: Bose construction requires n = 3 (mod 6)")
	}
	m := n / 3 // m = 2t+1, odd
	half := (m + 1) / 2
	op := func(a, b int) int { return half * (a + b) % m }
	idx := func(i, c int) int { return 3*i + c }
	var blocks [][]int
	// Vertical triples.
	for i := 0; i < m; i++ {
		blocks = append(blocks, []int{idx(i, 0), idx(i, 1), idx(i, 2)})
	}
	// Cross triples.
	for i := 0; i < m; i++ {
		for j := i + 1; j < m; j++ {
			ij := op(i, j)
			for c := 0; c < 3; c++ {
				blocks = append(blocks, []int{idx(i, c), idx(j, c), idx(ij, (c+1)%3)})
			}
		}
	}
	return NewDesign(n, blocks)
}

// SteinerTripleSystem constructs a Steiner triple system on v points using
// Stinson's hill-climbing algorithm, driven by the supplied random source so
// that the result is reproducible from a seed. It reports an error when v is not
// admissible (not congruent to 1 or 3 modulo 6).
func SteinerTripleSystem(v int, rng *rand.Rand) (*Design, error) {
	if !SteinerTripleAdmissible(v) {
		return nil, errors.New("designs: STS requires v = 1 or 3 (mod 6)")
	}
	if v <= 1 {
		return NewDesign(maxInt(v, 1), nil)
	}
	if v == 3 {
		return NewDesign(3, [][]int{{0, 1, 2}})
	}
	if rng == nil {
		rng = rand.New(rand.NewSource(1))
	}
	blocks, ok := hillClimbSTS(v, rng)
	if !ok {
		return nil, errors.New("designs: hill-climbing failed to converge")
	}
	return NewDesign(v, blocks)
}

// hillClimbSTS runs Stinson's hill-climbing algorithm for a Steiner triple
// system of order n, returning the triples and whether it converged.
func hillClimbSTS(n int, rng *rand.Rand) ([][]int, bool) {
	r := (n - 1) / 2
	target := n * (n - 1) / 6
	covered := make([][]bool, n)
	third := make([][]int, n)
	for i := 0; i < n; i++ {
		covered[i] = make([]bool, n)
		third[i] = make([]int, n)
	}
	degree := make([]int, n)
	type tri [3]int
	triples := make(map[tri]bool)

	setPair := func(a, b, c int) {
		covered[a][b], covered[b][a] = true, true
		third[a][b], third[b][a] = c, c
	}
	clearPair := func(a, b int) {
		covered[a][b], covered[b][a] = false, false
	}
	addTriple := func(x, y, z int) {
		s := sortTriple(x, y, z)
		triples[tri{s[0], s[1], s[2]}] = true
		setPair(x, y, z)
		setPair(x, z, y)
		setPair(y, z, x)
		degree[x]++
		degree[y]++
		degree[z]++
	}
	removeTriple := func(x, y, z int) {
		s := sortTriple(x, y, z)
		delete(triples, tri{s[0], s[1], s[2]})
		clearPair(x, y)
		clearPair(x, z)
		clearPair(y, z)
		degree[x]--
		degree[y]--
		degree[z]--
	}

	maxIter := 500*n*n + 100000
	for it := 0; it < maxIter; it++ {
		if len(triples) == target {
			out := make([][]int, 0, target)
			for t := range triples {
				out = append(out, []int{t[0], t[1], t[2]})
			}
			sort.Slice(out, func(a, b int) bool {
				if out[a][0] != out[b][0] {
					return out[a][0] < out[b][0]
				}
				if out[a][1] != out[b][1] {
					return out[a][1] < out[b][1]
				}
				return out[a][2] < out[b][2]
			})
			return out, true
		}
		// Choose a live point x.
		var live []int
		for i := 0; i < n; i++ {
			if degree[i] < r {
				live = append(live, i)
			}
		}
		x := live[rng.Intn(len(live))]
		// Uncovered partners of x.
		var part []int
		for y := 0; y < n; y++ {
			if y != x && !covered[x][y] {
				part = append(part, y)
			}
		}
		// Need at least two; guaranteed for a live point.
		if len(part) < 2 {
			continue
		}
		i1 := rng.Intn(len(part))
		i2 := rng.Intn(len(part))
		for i2 == i1 {
			i2 = rng.Intn(len(part))
		}
		y, z := part[i1], part[i2]
		if !covered[y][z] {
			addTriple(x, y, z)
		} else {
			w := third[y][z]
			removeTriple(y, z, w)
			addTriple(x, y, z)
		}
	}
	return nil, false
}

func sortTriple(a, b, c int) [3]int {
	s := [3]int{a, b, c}
	if s[0] > s[1] {
		s[0], s[1] = s[1], s[0]
	}
	if s[1] > s[2] {
		s[1], s[2] = s[2], s[1]
	}
	if s[0] > s[1] {
		s[0], s[1] = s[1], s[0]
	}
	return s
}

func maxInt(a, b int) int {
	if a > b {
		return a
	}
	return b
}

// IsSteinerTripleSystem reports whether the design is a Steiner triple system: a
// 2-(v,3,1) design in which every pair of points lies in exactly one block of
// size three.
func IsSteinerTripleSystem(d *Design) bool {
	p, err := d.Parameters()
	if err != nil {
		return false
	}
	return p.K == 3 && p.Lambda == 1
}

// SteinerSystemBlockCount returns the number of blocks C(v,t)/C(k,t) of a
// Steiner system S(t,k,v), the number of k-subsets needed so every t-subset is
// covered once. It reports an error when the divisibility fails or arguments are
// out of range.
func SteinerSystemBlockCount(t, k, v int) (int, error) {
	if t < 1 || k < t || v < k {
		return 0, errors.New("designs: require 1<=t<=k<=v")
	}
	num := Binomial(v, t)
	den := Binomial(k, t)
	if den == 0 || num%den != 0 {
		return 0, errors.New("designs: parameters do not yield an integer block count")
	}
	return num / den, nil
}

// SteinerQuadrupleAdmissible reports whether v is an admissible order for a
// Steiner quadruple system S(3,4,v), i.e. v is congruent to 2 or 4 modulo 6 (or
// v is 0 or 1).
func SteinerQuadrupleAdmissible(v int) bool {
	return v == 0 || v == 1 || v%6 == 2 || v%6 == 4
}

// SteinerQuadrupleParameters returns the number of blocks v(v-1)(v-2)/24 of a
// Steiner quadruple system on v points. It reports an error when v is not
// admissible.
func SteinerQuadrupleParameters(v int) (int, error) {
	if !SteinerQuadrupleAdmissible(v) {
		return 0, errors.New("designs: SQS requires v = 2 or 4 (mod 6)")
	}
	return v * (v - 1) * (v - 2) / 24, nil
}

// BooleanQuadrupleSystem constructs the Steiner quadruple system on the 2**k
// points of the vector space GF(2)^k, whose blocks are the 4-subsets
// {a,b,c,d} with a XOR b XOR c XOR d = 0. It reports an error when k<2.
func BooleanQuadrupleSystem(k int) (*Design, error) {
	if k < 2 {
		return nil, errors.New("designs: require k>=2")
	}
	n := 1 << k
	seen := make(map[[4]int]bool)
	var blocks [][]int
	for a := 0; a < n; a++ {
		for b := a + 1; b < n; b++ {
			for c := b + 1; c < n; c++ {
				d := a ^ b ^ c
				if d == a || d == b || d == c {
					continue
				}
				q := sortQuad(a, b, c, d)
				if !seen[q] {
					seen[q] = true
					blocks = append(blocks, []int{q[0], q[1], q[2], q[3]})
				}
			}
		}
	}
	return NewDesign(n, blocks)
}

// SteinerQuadrupleSystem constructs a Steiner quadruple system on v points. It
// is implemented for v a power of two (the Boolean construction) and for the
// trivial orders 1 and 2; it reports an error for other admissible orders,
// which require constructions beyond the scope of this package.
func SteinerQuadrupleSystem(v int) (*Design, error) {
	if !SteinerQuadrupleAdmissible(v) {
		return nil, errors.New("designs: SQS requires v = 2 or 4 (mod 6)")
	}
	if v <= 1 {
		return NewDesign(maxInt(v, 1), nil)
	}
	if v == 2 {
		return NewDesign(2, nil)
	}
	// Power of two?
	k := 0
	for t := v; t > 1; t >>= 1 {
		if t&1 == 1 {
			return nil, errors.New("designs: SQS construction only implemented for powers of two")
		}
		k++
	}
	return BooleanQuadrupleSystem(k)
}

// IsSteinerQuadrupleSystem reports whether the design is a Steiner quadruple
// system: every 3-subset of points lies in exactly one block of size four.
func IsSteinerQuadrupleSystem(d *Design) bool {
	for _, b := range d.Blocks {
		if len(b) != 4 {
			return false
		}
	}
	v := d.Points
	type tri [3]int
	count := make(map[tri]int)
	for _, b := range d.Blocks {
		for i := 0; i < 4; i++ {
			for j := i + 1; j < 4; j++ {
				for l := j + 1; l < 4; l++ {
					count[tri{b[i], b[j], b[l]}]++
				}
			}
		}
	}
	if len(count) != Binomial(v, 3) {
		return false
	}
	for _, c := range count {
		if c != 1 {
			return false
		}
	}
	return true
}

func sortQuad(a, b, c, d int) [4]int {
	s := []int{a, b, c, d}
	sort.Ints(s)
	return [4]int{s[0], s[1], s[2], s[3]}
}
