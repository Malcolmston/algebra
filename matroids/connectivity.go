package matroids

import "sort"

// ConnectivityFunction returns the matroid connectivity function
// λ(S) = r(S) + r(E\S) - r(E), where r is the rank function of m and E the
// ground set. It measures how strongly S is linked to its complement.
func ConnectivityFunction(m Matroid, set []int) int {
	s := distinctInRange(set, m.Size())
	comp := Complement(m.Size(), s)
	return m.Rank(s) + m.Rank(comp) - FullRank(m)
}

// IsSeparator reports whether set is a separator of m, i.e. λ(set) = 0. For a
// separator S, m is the direct sum of its restrictions to S and to E\S.
func IsSeparator(m Matroid, set []int) bool {
	return ConnectivityFunction(m, set) == 0
}

// Components returns the connected components of m as a partition of the ground
// set. Two elements lie in the same component exactly when some circuit of m
// contains both; loops and coloops form singleton components. Components are
// sorted, ordered by their smallest element.
func Components(m Matroid) [][]int {
	n := m.Size()
	uf := newUnionFind(n)
	for _, c := range Circuits(m) {
		for i := 1; i < len(c); i++ {
			uf.union(c[0], c[i])
		}
	}
	groups := make(map[int][]int)
	for e := 0; e < n; e++ {
		r := uf.find(e)
		groups[r] = append(groups[r], e)
	}
	var out [][]int
	for _, g := range groups {
		sort.Ints(g)
		out = append(out, g)
	}
	sort.SliceStable(out, func(i, j int) bool { return out[i][0] < out[j][0] })
	return out
}

// NumComponents returns the number of connected components of m.
func NumComponents(m Matroid) int { return len(Components(m)) }

// IsConnected reports whether m is connected: it has at most one component
// (every pair of elements lies in a common circuit). Matroids with fewer than
// two elements are connected by convention.
func IsConnected(m Matroid) bool {
	if m.Size() <= 1 {
		return true
	}
	return NumComponents(m) == 1
}

// ComponentOf returns the connected component of m that contains element e.
func ComponentOf(m Matroid, e int) []int {
	for _, comp := range Components(m) {
		if SetContains(comp, e) {
			return comp
		}
	}
	return nil
}
