package matroids

import (
	"math"
	"sort"
)

// appendElem returns a copy of base with e appended.
func appendElem(base []int, e int) []int {
	out := make([]int, len(base)+1)
	copy(out, base)
	out[len(base)] = e
	return out
}

// removeAndAdd returns the set base with x removed and y added.
func removeAndAdd(base []int, x, y int) []int {
	out := make([]int, 0, len(base)+1)
	for _, e := range base {
		if e != x {
			out = append(out, e)
		}
	}
	out = append(out, y)
	return out
}

// Intersection returns a maximum-cardinality set that is independent in both
// m1 and m2, computed by repeatedly finding a shortest augmenting path in the
// exchange graph. m1 and m2 must share a ground set (equal Size); it panics
// otherwise.
func Intersection(m1, m2 Matroid) []int {
	n := m1.Size()
	if m2.Size() != n {
		panic(ErrDimensionMismatch)
	}
	inI := make([]bool, n)
	for {
		path, ok := findAugment(m1, m2, inI)
		if !ok {
			break
		}
		for _, e := range path {
			inI[e] = !inI[e]
		}
	}
	var out []int
	for e := 0; e < n; e++ {
		if inI[e] {
			out = append(out, e)
		}
	}
	return out
}

// IntersectionSize returns the size of a maximum common independent set of m1
// and m2.
func IntersectionSize(m1, m2 Matroid) int { return len(Intersection(m1, m2)) }

// findAugment finds a shortest augmenting path in the exchange graph of the
// current common independent set (the elements marked in inI), or reports that
// none exists.
func findAugment(m1, m2 Matroid, inI []bool) ([]int, bool) {
	n := len(inI)
	var I, notI []int
	for e := 0; e < n; e++ {
		if inI[e] {
			I = append(I, e)
		} else {
			notI = append(notI, e)
		}
	}
	x1 := make([]bool, n)
	x2 := make([]bool, n)
	for _, y := range notI {
		if Independent(m1, appendElem(I, y)) {
			x1[y] = true
		}
		if Independent(m2, appendElem(I, y)) {
			x2[y] = true
		}
	}
	// direct augmentation: an element addable to both matroids.
	for y := range x1 {
		if x1[y] && x2[y] {
			return []int{y}, true
		}
	}
	adj := make([][]int, n)
	for _, x := range I {
		for _, y := range notI {
			s := removeAndAdd(I, x, y)
			if Independent(m1, s) {
				adj[x] = append(adj[x], y) // x -> y
			}
			if Independent(m2, s) {
				adj[y] = append(adj[y], x) // y -> x
			}
		}
	}
	prev := make([]int, n)
	for i := range prev {
		prev[i] = -2 // unvisited
	}
	var queue []int
	for y := 0; y < n; y++ {
		if x1[y] {
			prev[y] = -1
			queue = append(queue, y)
		}
	}
	target := -1
	for len(queue) > 0 && target == -1 {
		u := queue[0]
		queue = queue[1:]
		for _, v := range adj[u] {
			if prev[v] == -2 {
				prev[v] = u
				if x2[v] {
					target = v
					break
				}
				queue = append(queue, v)
			}
		}
	}
	if target == -1 {
		return nil, false
	}
	var path []int
	for v := target; v != -1; v = prev[v] {
		path = append(path, v)
	}
	return path, true
}

// WeightedIntersection returns a maximum-weight set that is independent in both
// m1 and m2. It runs the weighted matroid-intersection algorithm, augmenting
// along a maximum-weight (fewest-arcs among those) path in the exchange graph
// at each step and returning the best set seen over all cardinalities. weights
// must have length equal to the shared ground-set size. m1 and m2 must share a
// ground set; it panics otherwise.
func WeightedIntersection(m1, m2 Matroid, weights []float64) []int {
	n := m1.Size()
	if m2.Size() != n {
		panic(ErrDimensionMismatch)
	}
	if len(weights) != n {
		panic("matroids: weights length must equal Size()")
	}
	inI := make([]bool, n)
	var best []int
	bestW := 0.0 // empty set has weight 0
	for {
		path, ok := findWeightedAugment(m1, m2, inI, weights)
		if !ok {
			break
		}
		for _, e := range path {
			inI[e] = !inI[e]
		}
		var cur []int
		w := 0.0
		for e := 0; e < n; e++ {
			if inI[e] {
				cur = append(cur, e)
				w += weights[e]
			}
		}
		if w > bestW+1e-12 {
			bestW = w
			best = append([]int(nil), cur...)
		}
	}
	sort.Ints(best)
	return best
}

// WeightedIntersectionWeight returns the total weight of the set returned by
// [WeightedIntersection].
func WeightedIntersectionWeight(m1, m2 Matroid, weights []float64) float64 {
	set := WeightedIntersection(m1, m2, weights)
	total := 0.0
	for _, e := range set {
		total += weights[e]
	}
	return total
}

// findWeightedAugment finds a maximum-weight augmenting path (fewest arcs among
// maximum-weight ones) in the exchange graph, using Bellman-Ford over node
// lengths. It returns the path (as a set of elements to toggle) or false.
func findWeightedAugment(m1, m2 Matroid, inI []bool, weights []float64) ([]int, bool) {
	n := len(inI)
	var I, notI []int
	for e := 0; e < n; e++ {
		if inI[e] {
			I = append(I, e)
		} else {
			notI = append(notI, e)
		}
	}
	x1 := make([]bool, n)
	x2 := make([]bool, n)
	anyX1, anyX2 := false, false
	for _, y := range notI {
		if Independent(m1, appendElem(I, y)) {
			x1[y] = true
			anyX1 = true
		}
		if Independent(m2, appendElem(I, y)) {
			x2[y] = true
			anyX2 = true
		}
	}
	if !anyX1 || !anyX2 {
		return nil, false
	}
	adj := make([][]int, n)
	for _, x := range I {
		for _, y := range notI {
			s := removeAndAdd(I, x, y)
			if Independent(m1, s) {
				adj[x] = append(adj[x], y)
			}
			if Independent(m2, s) {
				adj[y] = append(adj[y], x)
			}
		}
	}
	nodeLen := func(v int) float64 {
		if inI[v] {
			return weights[v]
		}
		return -weights[v]
	}
	const eps = 1e-12
	inf := math.Inf(1)
	dist := make([]float64, n)
	hops := make([]int, n)
	prev := make([]int, n)
	for i := range dist {
		dist[i] = inf
		hops[i] = 1 << 30
		prev[i] = -1
	}
	for y := 0; y < n; y++ {
		if x1[y] {
			dist[y] = nodeLen(y)
			hops[y] = 0
			prev[y] = -1
		}
	}
	for iter := 0; iter < n+1; iter++ {
		changed := false
		for u := 0; u < n; u++ {
			if math.IsInf(dist[u], 1) {
				continue
			}
			for _, v := range adj[u] {
				nd := dist[u] + nodeLen(v)
				if nd < dist[v]-eps || (math.Abs(nd-dist[v]) <= eps && hops[u]+1 < hops[v]) {
					dist[v] = nd
					hops[v] = hops[u] + 1
					prev[v] = u
					changed = true
				}
			}
		}
		if !changed {
			break
		}
	}
	target := -1
	for v := 0; v < n; v++ {
		if !x2[v] || math.IsInf(dist[v], 1) {
			continue
		}
		if target == -1 || dist[v] < dist[target]-eps ||
			(math.Abs(dist[v]-dist[target]) <= eps && hops[v] < hops[target]) {
			target = v
		}
	}
	if target == -1 {
		return nil, false
	}
	var path []int
	for v := target; v != -1; v = prev[v] {
		path = append(path, v)
	}
	return path, true
}
