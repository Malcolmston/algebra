package matroids

import "sort"

// GreedyResult holds the outcome of a greedy matroid optimisation: the selected
// set (sorted) and its total weight.
type GreedyResult struct {
	Set    []int
	Weight float64
}

// Greedy runs the matroid greedy algorithm for a maximum-weight independent
// set. Elements are considered in order of non-increasing weight; an element is
// added whenever it keeps the current set independent and its weight is
// strictly positive. For any matroid this yields a maximum-weight independent
// set (an optimal solution of the weighted matroid problem). weights must have
// length m.Size(); it panics otherwise.
func Greedy(m Matroid, weights []float64) GreedyResult {
	if len(weights) != m.Size() {
		panic("matroids: weights length must equal Size()")
	}
	order := weightOrder(weights)
	var chosen []int
	total := 0.0
	for _, e := range order {
		if weights[e] <= 0 {
			break
		}
		trial := append(append([]int(nil), chosen...), e)
		if Independent(m, trial) {
			chosen = trial
			total += weights[e]
		}
	}
	sort.Ints(chosen)
	return GreedyResult{Set: chosen, Weight: total}
}

// GreedyMaxWeightBasis returns a maximum-weight basis of m using the greedy
// algorithm: elements are considered in order of non-increasing weight and each
// is added when it keeps the set independent, regardless of sign, until a basis
// is reached. weights must have length m.Size().
func GreedyMaxWeightBasis(m Matroid, weights []float64) GreedyResult {
	if len(weights) != m.Size() {
		panic("matroids: weights length must equal Size()")
	}
	order := weightOrder(weights)
	var chosen []int
	total := 0.0
	for _, e := range order {
		trial := append(append([]int(nil), chosen...), e)
		if Independent(m, trial) {
			chosen = trial
			total += weights[e]
		}
	}
	sort.Ints(chosen)
	return GreedyResult{Set: chosen, Weight: total}
}

// GreedyMinWeightBasis returns a minimum-weight basis of m. It is the greedy
// algorithm applied to negated weights; equivalently elements are considered in
// order of non-decreasing weight. This generalises Kruskal's minimum spanning
// tree algorithm (which is the graphic-matroid special case). weights must have
// length m.Size().
func GreedyMinWeightBasis(m Matroid, weights []float64) GreedyResult {
	if len(weights) != m.Size() {
		panic("matroids: weights length must equal Size()")
	}
	neg := make([]float64, len(weights))
	for i, w := range weights {
		neg[i] = -w
	}
	res := GreedyMaxWeightBasis(m, neg)
	res.Weight = -res.Weight
	return res
}

// MaxWeightIndependentSet returns the set found by [Greedy]; it is a convenience
// wrapper returning only the element list.
func MaxWeightIndependentSet(m Matroid, weights []float64) []int {
	return Greedy(m, weights).Set
}

// weightOrder returns the ground-element indices sorted by non-increasing
// weight, ties broken by index for determinism.
func weightOrder(weights []float64) []int {
	order := make([]int, len(weights))
	for i := range order {
		order[i] = i
	}
	sort.SliceStable(order, func(i, j int) bool {
		if weights[order[i]] != weights[order[j]] {
			return weights[order[i]] > weights[order[j]]
		}
		return order[i] < order[j]
	})
	return order
}
