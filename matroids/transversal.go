package matroids

// TransversalMatroid is the transversal matroid of a family of sets. The
// ground set is {0, ..., n-1}. Each ground element e is associated with a set
// of "targets" (right vertices). A subset S of the ground set is independent
// exactly when it admits a system of distinct representatives — a matching that
// assigns each element of S a distinct target — so the rank of S is the size of
// a maximum matching saturating a subset of S. Ranks are computed by bipartite
// matching (Kuhn's augmenting-path algorithm).
type TransversalMatroid struct {
	n       int
	targets [][]int // targets[e] lists the right vertices reachable from e
	nr      int     // number of distinct right vertices used
}

// NewTransversalMatroid builds a transversal matroid on n ground elements. adj
// must have length n; adj[e] lists the target vertices associated with element
// e (arbitrary non-negative integers). Targets are re-indexed internally. It
// panics on a length mismatch or a negative target.
func NewTransversalMatroid(n int, adj [][]int) *TransversalMatroid {
	if len(adj) != n {
		panic("matroids: adj length must equal n")
	}
	remap := make(map[int]int)
	targets := make([][]int, n)
	for e := 0; e < n; e++ {
		seen := make(map[int]bool)
		for _, t := range adj[e] {
			if t < 0 {
				panic("matroids: negative target vertex")
			}
			if seen[t] {
				continue
			}
			seen[t] = true
			idx, ok := remap[t]
			if !ok {
				idx = len(remap)
				remap[t] = idx
			}
			targets[e] = append(targets[e], idx)
		}
	}
	return &TransversalMatroid{n: n, targets: targets, nr: len(remap)}
}

// Size returns the number of ground-set elements.
func (m *TransversalMatroid) Size() int { return m.n }

// NumTargets returns the number of distinct target (right) vertices.
func (m *TransversalMatroid) NumTargets() int { return m.nr }

// Targets returns the internal target indices associated with element e.
func (m *TransversalMatroid) Targets(e int) []int {
	out := make([]int, len(m.targets[e]))
	copy(out, m.targets[e])
	return out
}

// Rank returns the size of a maximum matching that saturates a subset of the
// distinct in-range elements of set, i.e. the largest partial transversal
// available from those elements.
func (m *TransversalMatroid) Rank(set []int) int {
	left := distinctInRange(set, m.n)
	matchTo := make([]int, m.nr) // right vertex -> matched left element, or -1
	for i := range matchTo {
		matchTo[i] = -1
	}
	count := 0
	for _, e := range left {
		visited := make([]bool, m.nr)
		if m.augment(e, matchTo, visited) {
			count++
		}
	}
	return count
}

// augment tries to find an augmenting path from left element e using Kuhn's
// algorithm; matchTo maps each right vertex to its matched left element.
func (m *TransversalMatroid) augment(e int, matchTo []int, visited []bool) bool {
	for _, t := range m.targets[e] {
		if visited[t] {
			continue
		}
		visited[t] = true
		if matchTo[t] == -1 || m.augment(matchTo[t], matchTo, visited) {
			matchTo[t] = e
			return true
		}
	}
	return false
}

// SystemOfDistinctRepresentatives returns, for an independent set S (given as
// set), an assignment mapping each element of S to a distinct target index, or
// nil together with false if S is dependent (no such system exists). The
// mapping keys are the original ground elements; the values are internal target
// indices.
func (m *TransversalMatroid) SystemOfDistinctRepresentatives(set []int) (map[int]int, bool) {
	left := distinctInRange(set, m.n)
	matchTo := make([]int, m.nr)
	for i := range matchTo {
		matchTo[i] = -1
	}
	count := 0
	for _, e := range left {
		visited := make([]bool, m.nr)
		if m.augment(e, matchTo, visited) {
			count++
		}
	}
	if count != len(left) {
		return nil, false
	}
	assign := make(map[int]int, len(left))
	for t, e := range matchTo {
		if e != -1 {
			assign[e] = t
		}
	}
	return assign, true
}
