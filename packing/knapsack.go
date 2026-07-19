package packing

import (
	"errors"
	"math"
	"sort"
)

// ErrDimensionMismatch is returned when parallel slices (values and weights, or
// sets and costs) have different lengths.
var ErrDimensionMismatch = errors.New("packing: slice length mismatch")

// ErrNegativeWeight is returned when a knapsack weight or capacity is negative.
var ErrNegativeWeight = errors.New("packing: negative weight or capacity")

// ErrNotCoverable is returned when the union of the given sets does not cover
// the whole universe, so no set cover exists.
var ErrNotCoverable = errors.New("packing: universe is not coverable by the given sets")

// ----------------------------------------------------------------------------
// Knapsack.
// ----------------------------------------------------------------------------

// ratioOrder returns item indices sorted by non-increasing value-to-weight
// ratio, with zero-weight positive-value items first (infinite ratio) and ties
// broken by index for determinism.
func ratioOrder(values, weights []float64) []int {
	idx := make([]int, len(values))
	for i := range idx {
		idx[i] = i
	}
	ratio := func(i int) float64 {
		if weights[i] <= 0 {
			if values[i] > 0 {
				return math.Inf(1)
			}
			return 0
		}
		return values[i] / weights[i]
	}
	sort.SliceStable(idx, func(a, b int) bool {
		return ratio(idx[a]) > ratio(idx[b])
	})
	return idx
}

// FractionalKnapsack solves the fractional (continuous) knapsack problem
// exactly by the greedy value-to-weight rule: items are taken whole in order of
// decreasing density until the capacity is exhausted, and the final item is
// taken fractionally. It returns the optimal value and the fraction (in [0,1])
// taken of each item in the original order.
func FractionalKnapsack(values, weights []float64, capacity float64) (value float64, fraction []float64, err error) {
	if len(values) != len(weights) {
		return 0, nil, ErrDimensionMismatch
	}
	if capacity < 0 {
		return 0, nil, ErrNegativeWeight
	}
	for _, w := range weights {
		if w < 0 {
			return 0, nil, ErrNegativeWeight
		}
	}
	fraction = make([]float64, len(values))
	remaining := capacity
	for _, i := range ratioOrder(values, weights) {
		if weights[i] <= 0 {
			// Zero-weight item: take fully if it has positive value.
			if values[i] > 0 {
				fraction[i] = 1
				value += values[i]
			}
			continue
		}
		if remaining <= 0 {
			break
		}
		if weights[i] <= remaining {
			fraction[i] = 1
			value += values[i]
			remaining -= weights[i]
		} else {
			f := remaining / weights[i]
			fraction[i] = f
			value += f * values[i]
			remaining = 0
		}
	}
	return value, fraction, nil
}

// GreedyKnapsack applies the greedy value-to-weight heuristic to the 0/1
// knapsack problem: it considers items in order of decreasing density and takes
// each whole item that still fits. It returns the achieved value, the sorted
// indices of the chosen items and their total weight. The result is a feasible
// lower bound on the optimum, not necessarily optimal.
func GreedyKnapsack(values, weights []float64, capacity float64) (value float64, chosen []int, weight float64, err error) {
	if len(values) != len(weights) {
		return 0, nil, 0, ErrDimensionMismatch
	}
	if capacity < 0 {
		return 0, nil, 0, ErrNegativeWeight
	}
	for _, w := range weights {
		if w < 0 {
			return 0, nil, 0, ErrNegativeWeight
		}
	}
	remaining := capacity
	for _, i := range ratioOrder(values, weights) {
		if weights[i] <= remaining+eps {
			chosen = append(chosen, i)
			value += values[i]
			weight += weights[i]
			remaining -= weights[i]
		}
	}
	sort.Ints(chosen)
	return value, chosen, weight, nil
}

// GreedyKnapsackBest runs the greedy heuristic but also compares it against the
// single most valuable item that fits, returning whichever gives more value.
// This modified greedy achieves at least half of the optimum for the 0/1
// knapsack problem.
func GreedyKnapsackBest(values, weights []float64, capacity float64) (value float64, chosen []int, weight float64, err error) {
	gv, gc, gw, err := GreedyKnapsack(values, weights, capacity)
	if err != nil {
		return 0, nil, 0, err
	}
	// Best single item that fits on its own.
	bestI, bestV := -1, 0.0
	for i := range values {
		if weights[i] <= capacity+eps && values[i] > bestV {
			bestV, bestI = values[i], i
		}
	}
	if bestI >= 0 && bestV > gv {
		return bestV, []int{bestI}, weights[bestI], nil
	}
	return gv, gc, gw, nil
}

// KnapsackUpperBound returns the linear-programming (fractional) relaxation
// value of the 0/1 knapsack instance, an upper bound on the integer optimum.
func KnapsackUpperBound(values, weights []float64, capacity float64) (float64, error) {
	v, _, err := FractionalKnapsack(values, weights, capacity)
	return v, err
}

// KnapsackValue returns the total value of the items whose indices are given.
func KnapsackValue(values []float64, chosen []int) float64 {
	s := 0.0
	for _, i := range chosen {
		s += values[i]
	}
	return s
}

// KnapsackWeight returns the total weight of the items whose indices are given.
func KnapsackWeight(weights []float64, chosen []int) float64 {
	s := 0.0
	for _, i := range chosen {
		s += weights[i]
	}
	return s
}

// KnapsackFeasible reports whether the chosen items fit within the capacity.
func KnapsackFeasible(weights []float64, chosen []int, capacity float64) bool {
	return KnapsackWeight(weights, chosen) <= capacity+eps
}

// ----------------------------------------------------------------------------
// Set cover.
// ----------------------------------------------------------------------------

// coverable reports whether the union of the sets contains every element of the
// universe {0,...,universe-1} and returns the set of covered elements.
func coverable(universe int, sets [][]int) (bool, []bool) {
	covered := make([]bool, universe)
	for _, s := range sets {
		for _, e := range s {
			if e >= 0 && e < universe {
				covered[e] = true
			}
		}
	}
	for _, c := range covered {
		if !c {
			return false, covered
		}
	}
	return true, covered
}

// SetCoverGreedy applies the greedy set-cover heuristic to cover the universe
// {0,...,universe-1}: repeatedly choose the set covering the most still
// uncovered elements. It returns the chosen set indices in the order selected,
// or [ErrNotCoverable] if the sets do not cover the universe. The number of
// sets chosen is at most H(k) times the optimum, where k is the size of the
// largest set (see [SetCoverGreedyBound]).
func SetCoverGreedy(universe int, sets [][]int) ([]int, error) {
	if universe <= 0 {
		return nil, nil
	}
	if ok, _ := coverable(universe, sets); !ok {
		return nil, ErrNotCoverable
	}
	covered := make([]bool, universe)
	remaining := universe
	var chosen []int
	used := make([]bool, len(sets))
	for remaining > 0 {
		best, bestGain := -1, 0
		for j, s := range sets {
			if used[j] {
				continue
			}
			gain := 0
			for _, e := range s {
				if e >= 0 && e < universe && !covered[e] {
					gain++
				}
			}
			if gain > bestGain {
				bestGain, best = gain, j
			}
		}
		if best < 0 {
			// Should not happen given coverable check.
			return nil, ErrNotCoverable
		}
		used[best] = true
		chosen = append(chosen, best)
		for _, e := range sets[best] {
			if e >= 0 && e < universe && !covered[e] {
				covered[e] = true
				remaining--
			}
		}
	}
	return chosen, nil
}

// WeightedSetCoverGreedy applies the greedy heuristic to the weighted set-cover
// problem: repeatedly choose the set minimizing cost divided by the number of
// newly covered elements. costs[j] is the cost of sets[j]. It returns the
// chosen set indices in selection order and their total cost, or
// [ErrNotCoverable] if the sets do not cover the universe.
func WeightedSetCoverGreedy(universe int, sets [][]int, costs []float64) (chosen []int, totalCost float64, err error) {
	if len(sets) != len(costs) {
		return nil, 0, ErrDimensionMismatch
	}
	if universe <= 0 {
		return nil, 0, nil
	}
	for _, c := range costs {
		if c < 0 {
			return nil, 0, ErrNegativeWeight
		}
	}
	if ok, _ := coverable(universe, sets); !ok {
		return nil, 0, ErrNotCoverable
	}
	covered := make([]bool, universe)
	remaining := universe
	used := make([]bool, len(sets))
	for remaining > 0 {
		best := -1
		bestEff := math.Inf(1)
		bestGain := 0
		for j, s := range sets {
			if used[j] {
				continue
			}
			gain := 0
			for _, e := range s {
				if e >= 0 && e < universe && !covered[e] {
					gain++
				}
			}
			if gain == 0 {
				continue
			}
			eff := costs[j] / float64(gain)
			if eff < bestEff-1e-12 || (math.Abs(eff-bestEff) <= 1e-12 && gain > bestGain) {
				bestEff, best, bestGain = eff, j, gain
			}
		}
		if best < 0 {
			return nil, 0, ErrNotCoverable
		}
		used[best] = true
		chosen = append(chosen, best)
		totalCost += costs[best]
		for _, e := range sets[best] {
			if e >= 0 && e < universe && !covered[e] {
				covered[e] = true
				remaining--
			}
		}
	}
	return chosen, totalCost, nil
}

// IsSetCover reports whether the chosen sets together cover the whole universe
// {0,...,universe-1}.
func IsSetCover(universe int, sets [][]int, chosen []int) bool {
	covered := make([]bool, universe)
	remaining := universe
	for _, j := range chosen {
		if j < 0 || j >= len(sets) {
			continue
		}
		for _, e := range sets[j] {
			if e >= 0 && e < universe && !covered[e] {
				covered[e] = true
				remaining--
			}
		}
	}
	return remaining == 0
}

// HarmonicNumber returns the n-th harmonic number H(n) = 1 + 1/2 + ... + 1/n,
// with H(0) = 0. It is computed by direct summation.
func HarmonicNumber(n int) float64 {
	s := 0.0
	for k := 1; k <= n; k++ {
		s += 1 / float64(k)
	}
	return s
}

// SetCoverGreedyBound returns the greedy set-cover approximation guarantee
// H(k), the k-th harmonic number, where k is the size of the largest set. The
// greedy solution uses at most H(k) times as many sets (or, weighted, H(k)
// times the cost) as the optimum.
func SetCoverGreedyBound(largestSetSize int) float64 {
	return HarmonicNumber(largestSetSize)
}

// MaxSetSize returns the size of the largest set, used as the parameter of
// [SetCoverGreedyBound].
func MaxSetSize(sets [][]int) int {
	m := 0
	for _, s := range sets {
		if len(s) > m {
			m = len(s)
		}
	}
	return m
}

// ModifiedGreedyKnapsackApproxRatio returns the worst-case approximation ratio
// of the modified greedy 0/1 knapsack heuristic [GreedyKnapsackBest], which is
// 2: it always achieves at least half of the optimal value.
func ModifiedGreedyKnapsackApproxRatio() float64 { return 2.0 }
