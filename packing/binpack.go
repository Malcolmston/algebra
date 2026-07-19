package packing

import (
	"errors"
	"math"
	"sort"
)

// eps is the tolerance used when comparing floating point loads against a bin
// capacity, so that an item that exactly fills the remaining space is accepted
// despite rounding noise.
const eps = 1e-9

// ErrCapacity is returned when a bin capacity is not strictly positive.
var ErrCapacity = errors.New("packing: bin capacity must be positive")

// ErrNegativeSize is returned when an item has negative size.
var ErrNegativeSize = errors.New("packing: item size must be non-negative")

// ErrItemTooLarge is returned when an item is larger than the bin capacity and
// therefore cannot be packed into any bin.
var ErrItemTooLarge = errors.New("packing: item larger than bin capacity")

// Bin is a single bin of a [Packing]. Items holds the indices (into the
// original size slice passed to the packing routine) of the items assigned to
// this bin, in the order they were placed. Load is the sum of their sizes.
type Bin struct {
	Items []int
	Load  float64
}

// Packing is the result of a bin packing routine: an assignment of every item
// to exactly one [Bin] together with the capacity and the original item sizes.
type Packing struct {
	Bins     []Bin
	Capacity float64
	Sizes    []float64
}

// NumBins reports the number of bins used by the packing.
func (p Packing) NumBins() int { return len(p.Bins) }

// NumItems reports the total number of packed items.
func (p Packing) NumItems() int {
	n := 0
	for _, b := range p.Bins {
		n += len(b.Items)
	}
	return n
}

// Loads returns the load (sum of item sizes) of each bin in bin order.
func (p Packing) Loads() []float64 {
	out := make([]float64, len(p.Bins))
	for i, b := range p.Bins {
		out[i] = b.Load
	}
	return out
}

// TotalSize returns the sum of all item sizes.
func (p Packing) TotalSize() float64 {
	s := 0.0
	for _, x := range p.Sizes {
		s += x
	}
	return s
}

// MaxLoad returns the largest bin load, or 0 for an empty packing.
func (p Packing) MaxLoad() float64 {
	m := 0.0
	for _, b := range p.Bins {
		if b.Load > m {
			m = b.Load
		}
	}
	return m
}

// MinLoad returns the smallest bin load, or 0 for an empty packing.
func (p Packing) MinLoad() float64 {
	if len(p.Bins) == 0 {
		return 0
	}
	m := math.Inf(1)
	for _, b := range p.Bins {
		if b.Load < m {
			m = b.Load
		}
	}
	return m
}

// Waste returns the total unused capacity, NumBins*Capacity - TotalSize.
func (p Packing) Waste() float64 {
	return float64(p.NumBins())*p.Capacity - p.TotalSize()
}

// Fullness returns the average fraction of bin capacity that is used, in the
// range [0,1]. An empty packing has fullness 0.
func (p Packing) Fullness() float64 {
	if len(p.Bins) == 0 {
		return 0
	}
	return p.TotalSize() / (float64(p.NumBins()) * p.Capacity)
}

// Valid reports whether the packing is internally consistent: every item index
// appears exactly once, every recorded load matches the sum of its items, and
// no bin load exceeds the capacity (within tolerance).
func (p Packing) Valid() bool {
	seen := make([]bool, len(p.Sizes))
	for _, b := range p.Bins {
		load := 0.0
		for _, idx := range b.Items {
			if idx < 0 || idx >= len(p.Sizes) || seen[idx] {
				return false
			}
			seen[idx] = true
			load += p.Sizes[idx]
		}
		if math.Abs(load-b.Load) > 1e-6 || load > p.Capacity+eps {
			return false
		}
	}
	for _, s := range seen {
		if !s {
			return false
		}
	}
	return true
}

// validateInput checks the shared preconditions of every packing routine.
func validateInput(sizes []float64, capacity float64) error {
	if capacity <= 0 {
		return ErrCapacity
	}
	for _, s := range sizes {
		if s < 0 {
			return ErrNegativeSize
		}
		if s > capacity+eps {
			return ErrItemTooLarge
		}
	}
	return nil
}

// order returns the indices 0..n-1 either in original order or sorted by
// non-increasing size (decreasing == true), ties broken by index.
func order(sizes []float64, decreasing bool) []int {
	idx := make([]int, len(sizes))
	for i := range idx {
		idx[i] = i
	}
	if decreasing {
		sort.SliceStable(idx, func(a, b int) bool {
			return sizes[idx[a]] > sizes[idx[b]]
		})
	}
	return idx
}

// place appends item idx to bin b, updating its load.
func (p *Packing) place(b int, idx int) {
	p.Bins[b].Items = append(p.Bins[b].Items, idx)
	p.Bins[b].Load += p.Sizes[idx]
}

// newPacking builds an empty packing shell sharing the given sizes.
func newPacking(sizes []float64, capacity float64) Packing {
	return Packing{Capacity: capacity, Sizes: sizes}
}

// nextFit packs the items in the given index order using next-fit: keep a
// single open bin and close it (open a new one) as soon as an item does not
// fit.
func nextFit(sizes []float64, capacity float64, idx []int) Packing {
	p := newPacking(sizes, capacity)
	cur := -1
	for _, i := range idx {
		if cur >= 0 && p.Bins[cur].Load+sizes[i] <= capacity+eps {
			p.place(cur, i)
			continue
		}
		p.Bins = append(p.Bins, Bin{})
		cur = len(p.Bins) - 1
		p.place(cur, i)
	}
	return p
}

// firstFit packs the items in the given index order into the lowest-indexed
// bin that has room, opening a new bin only when none fits.
func firstFit(sizes []float64, capacity float64, idx []int) Packing {
	p := newPacking(sizes, capacity)
	for _, i := range idx {
		placed := false
		for b := range p.Bins {
			if p.Bins[b].Load+sizes[i] <= capacity+eps {
				p.place(b, i)
				placed = true
				break
			}
		}
		if !placed {
			p.Bins = append(p.Bins, Bin{})
			p.place(len(p.Bins)-1, i)
		}
	}
	return p
}

// bestFit packs each item into the fullest bin that still has room (leaving the
// least residual capacity), opening a new bin when none fits.
func bestFit(sizes []float64, capacity float64, idx []int) Packing {
	p := newPacking(sizes, capacity)
	for _, i := range idx {
		best := -1
		bestResidual := math.Inf(1)
		for b := range p.Bins {
			res := capacity - p.Bins[b].Load - sizes[i]
			if res >= -eps && res < bestResidual {
				bestResidual = res
				best = b
			}
		}
		if best < 0 {
			p.Bins = append(p.Bins, Bin{})
			best = len(p.Bins) - 1
		}
		p.place(best, i)
	}
	return p
}

// worstFit packs each item into the emptiest bin that still has room (leaving
// the most residual capacity), opening a new bin when none fits.
func worstFit(sizes []float64, capacity float64, idx []int) Packing {
	p := newPacking(sizes, capacity)
	for _, i := range idx {
		worst := -1
		worstResidual := -1.0
		for b := range p.Bins {
			res := capacity - p.Bins[b].Load - sizes[i]
			if res >= -eps && res > worstResidual {
				worstResidual = res
				worst = b
			}
		}
		if worst < 0 {
			p.Bins = append(p.Bins, Bin{})
			worst = len(p.Bins) - 1
		}
		p.place(worst, i)
	}
	return p
}

// almostWorstFit is like worst-fit but places each item into the second
// emptiest feasible bin, which empirically improves on worst-fit.
func almostWorstFit(sizes []float64, capacity float64, idx []int) Packing {
	p := newPacking(sizes, capacity)
	for _, i := range idx {
		// Find the two emptiest feasible bins.
		best1, best2 := -1, -1
		res1, res2 := -1.0, -1.0
		for b := range p.Bins {
			res := capacity - p.Bins[b].Load - sizes[i]
			if res < -eps {
				continue
			}
			if res > res1 {
				best2, res2 = best1, res1
				best1, res1 = b, res
			} else if res > res2 {
				best2, res2 = b, res
			}
		}
		target := best2
		if target < 0 {
			target = best1
		}
		if target < 0 {
			p.Bins = append(p.Bins, Bin{})
			target = len(p.Bins) - 1
		}
		p.place(target, i)
	}
	return p
}

// NextFit packs the items using the next-fit heuristic. It returns an error if
// the capacity is not positive or an item is negative or larger than capacity.
func NextFit(sizes []float64, capacity float64) (Packing, error) {
	if err := validateInput(sizes, capacity); err != nil {
		return Packing{}, err
	}
	return nextFit(sizes, capacity, order(sizes, false)), nil
}

// FirstFit packs the items using the first-fit heuristic.
func FirstFit(sizes []float64, capacity float64) (Packing, error) {
	if err := validateInput(sizes, capacity); err != nil {
		return Packing{}, err
	}
	return firstFit(sizes, capacity, order(sizes, false)), nil
}

// BestFit packs the items using the best-fit heuristic.
func BestFit(sizes []float64, capacity float64) (Packing, error) {
	if err := validateInput(sizes, capacity); err != nil {
		return Packing{}, err
	}
	return bestFit(sizes, capacity, order(sizes, false)), nil
}

// WorstFit packs the items using the worst-fit heuristic.
func WorstFit(sizes []float64, capacity float64) (Packing, error) {
	if err := validateInput(sizes, capacity); err != nil {
		return Packing{}, err
	}
	return worstFit(sizes, capacity, order(sizes, false)), nil
}

// AlmostWorstFit packs the items using the almost-worst-fit heuristic.
func AlmostWorstFit(sizes []float64, capacity float64) (Packing, error) {
	if err := validateInput(sizes, capacity); err != nil {
		return Packing{}, err
	}
	return almostWorstFit(sizes, capacity, order(sizes, false)), nil
}

// NextFitDecreasing sorts the items by non-increasing size and packs them with
// next-fit.
func NextFitDecreasing(sizes []float64, capacity float64) (Packing, error) {
	if err := validateInput(sizes, capacity); err != nil {
		return Packing{}, err
	}
	return nextFit(sizes, capacity, order(sizes, true)), nil
}

// FirstFitDecreasing sorts the items by non-increasing size and packs them with
// first-fit. It is the classical offline heuristic with asymptotic ratio 11/9.
func FirstFitDecreasing(sizes []float64, capacity float64) (Packing, error) {
	if err := validateInput(sizes, capacity); err != nil {
		return Packing{}, err
	}
	return firstFit(sizes, capacity, order(sizes, true)), nil
}

// BestFitDecreasing sorts the items by non-increasing size and packs them with
// best-fit.
func BestFitDecreasing(sizes []float64, capacity float64) (Packing, error) {
	if err := validateInput(sizes, capacity); err != nil {
		return Packing{}, err
	}
	return bestFit(sizes, capacity, order(sizes, true)), nil
}

// WorstFitDecreasing sorts the items by non-increasing size and packs them with
// worst-fit.
func WorstFitDecreasing(sizes []float64, capacity float64) (Packing, error) {
	if err := validateInput(sizes, capacity); err != nil {
		return Packing{}, err
	}
	return worstFit(sizes, capacity, order(sizes, true)), nil
}

// AlmostWorstFitDecreasing sorts the items by non-increasing size and packs
// them with almost-worst-fit.
func AlmostWorstFitDecreasing(sizes []float64, capacity float64) (Packing, error) {
	if err := validateInput(sizes, capacity); err != nil {
		return Packing{}, err
	}
	return almostWorstFit(sizes, capacity, order(sizes, true)), nil
}

// LowerBoundL1 returns the trivial "continuous" lower bound on the number of
// bins: the total size divided by the capacity, rounded up. Every feasible
// packing uses at least this many bins.
func LowerBoundL1(sizes []float64, capacity float64) int {
	if capacity <= 0 {
		return 0
	}
	total := 0.0
	for _, s := range sizes {
		total += s
	}
	return int(math.Ceil(total/capacity - eps))
}

// LowerBoundL2 returns the Martello-Toth L2 lower bound on the number of bins.
// It dominates L1 by counting, for every threshold alpha, the large items that
// must occupy a bin by themselves and the medium items that cannot share a bin
// with one another, then adding the rounded-up residual demand of the small
// items. The maximum over all thresholds is returned.
func LowerBoundL2(sizes []float64, capacity float64) int {
	if capacity <= 0 {
		return 0
	}
	// Candidate thresholds: 0 and every item size <= C/2.
	thresholds := []float64{0}
	for _, s := range sizes {
		if s <= capacity/2+eps {
			thresholds = append(thresholds, s)
		}
	}
	best := LowerBoundL1(sizes, capacity)
	for _, alpha := range thresholds {
		var n1 int             // size > C - alpha  (alone)
		var n2 int             // C/2 < size <= C - alpha
		var sum2, sum3 float64 // total size of the J2 and J3 classes
		for _, s := range sizes {
			switch {
			case s > capacity-alpha+eps:
				n1++
			case s > capacity/2+eps:
				n2++
				sum2 += s
			case s >= alpha-eps:
				sum3 += s
			}
		}
		// Free space in the J2 bins that small items may occupy.
		free := float64(n2)*capacity - sum2
		residual := sum3 - free
		extra := 0
		if residual > eps {
			extra = int(math.Ceil(residual/capacity - eps))
		}
		if v := n1 + n2 + extra; v > best {
			best = v
		}
	}
	return best
}

// NextFitApproxRatio returns the asymptotic worst-case approximation ratio of
// next-fit relative to the optimum, which is 2.
func NextFitApproxRatio() float64 { return 2.0 }

// FirstFitApproxRatio returns the asymptotic worst-case approximation ratio of
// first-fit and best-fit, which is 1.7.
func FirstFitApproxRatio() float64 { return 1.7 }

// FirstFitDecreasingApproxRatio returns the asymptotic worst-case
// approximation ratio of first-fit-decreasing, which is 11/9.
func FirstFitDecreasingApproxRatio() float64 { return 11.0 / 9.0 }

// FirstFitDecreasingAdditiveBound returns the additive constant in the tight
// bound FFD(I) <= (11/9)OPT(I) + 6/9 due to Doesa et al.
func FirstFitDecreasingAdditiveBound() float64 { return 6.0 / 9.0 }

// BestFitApproxRatio returns the asymptotic worst-case approximation ratio of
// best-fit, which is 1.7 (identical to first-fit).
func BestFitApproxRatio() float64 { return 1.7 }

// WorstFitApproxRatio returns the asymptotic worst-case approximation ratio of
// worst-fit, which is 2 (identical to next-fit).
func WorstFitApproxRatio() float64 { return 2.0 }

// NextFitDecreasingApproxRatio returns the asymptotic worst-case approximation
// ratio of next-fit-decreasing, the constant h_infinity ~ 1.6910.
func NextFitDecreasingApproxRatio() float64 { return 1.6910302 }

// SumSizes returns the total of the item sizes.
func SumSizes(sizes []float64) float64 {
	s := 0.0
	for _, x := range sizes {
		s += x
	}
	return s
}

// MaxSize returns the largest item size, or 0 for an empty slice.
func MaxSize(sizes []float64) float64 {
	m := 0.0
	for _, x := range sizes {
		if x > m {
			m = x
		}
	}
	return m
}

// BinPackingLowerBound returns the strongest lower bound implemented here on the
// number of bins, the maximum of [LowerBoundL1] and [LowerBoundL2].
func BinPackingLowerBound(sizes []float64, capacity float64) int {
	l1 := LowerBoundL1(sizes, capacity)
	l2 := LowerBoundL2(sizes, capacity)
	if l2 > l1 {
		return l2
	}
	return l1
}
