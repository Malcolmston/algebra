package tensornetwork

import (
	"errors"
	"fmt"
	"math"
	"math/bits"
	"sort"
)

// Network is a tensor network: a collection of tensors whose axes carry string
// labels. A label shared by two tensors denotes a contracted (summed) index; a
// label appearing on a single tensor is a free output index. Every label must
// appear on at most two tensors and be distinct within a single tensor.
type Network struct {
	tensors []*Tensor
	labels  [][]string
}

// NewNetwork returns an empty tensor network.
func NewNetwork() *Network { return &Network{} }

// Add inserts a tensor into the network with the given axis labels. The number
// of labels must equal the tensor rank and the labels must be distinct. It
// returns an error otherwise.
func (nw *Network) Add(t *Tensor, labels ...string) error {
	if len(labels) != len(t.shape) {
		return fmt.Errorf("tensornetwork: %d labels for rank-%d tensor", len(labels), len(t.shape))
	}
	seen := make(map[string]bool)
	for _, l := range labels {
		if seen[l] {
			return fmt.Errorf("tensornetwork: repeated label %q within a tensor", l)
		}
		seen[l] = true
	}
	nw.tensors = append(nw.tensors, t)
	nw.labels = append(nw.labels, append([]string(nil), labels...))
	return nil
}

// Len returns the number of tensors currently in the network.
func (nw *Network) Len() int { return len(nw.tensors) }

// validate checks that every label appears at most twice and dimensions agree,
// returning the map of label -> dimension.
func (nw *Network) validate() (map[string]int, error) {
	dim := make(map[string]int)
	count := make(map[string]int)
	for ti, labs := range nw.labels {
		for ai, l := range labs {
			d := nw.tensors[ti].shape[ai]
			if prev, ok := dim[l]; ok {
				if prev != d {
					return nil, fmt.Errorf("tensornetwork: label %q dimension mismatch %d != %d", l, prev, d)
				}
			} else {
				dim[l] = d
			}
			count[l]++
			if count[l] > 2 {
				return nil, fmt.Errorf("tensornetwork: label %q appears more than twice", l)
			}
		}
	}
	return dim, nil
}

// FreeLabels returns the sorted list of labels that appear on exactly one
// tensor: the free output indices of the network.
func (nw *Network) FreeLabels() []string {
	count := make(map[string]int)
	for _, labs := range nw.labels {
		for _, l := range labs {
			count[l]++
		}
	}
	var free []string
	for l, c := range count {
		if c == 1 {
			free = append(free, l)
		}
	}
	sort.Strings(free)
	return free
}

// contractPair contracts two labelled tensors over their shared labels and
// returns the resulting tensor and its labels (free labels of a followed by free
// labels of b).
func contractPair(a *Tensor, la []string, b *Tensor, lb []string) (*Tensor, []string, error) {
	posB := make(map[string]int)
	for i, l := range lb {
		posB[l] = i
	}
	var axesA, axesB []int
	commonA := make(map[int]bool)
	for i, l := range la {
		if j, ok := posB[l]; ok {
			axesA = append(axesA, i)
			axesB = append(axesB, j)
			commonA[i] = true
		}
	}
	commonB := make(map[int]bool)
	for _, j := range axesB {
		commonB[j] = true
	}
	res, err := TensorDot(a, b, axesA, axesB)
	if err != nil {
		return nil, nil, err
	}
	var outLabels []string
	for i, l := range la {
		if !commonA[i] {
			outLabels = append(outLabels, l)
		}
	}
	for j, l := range lb {
		if !commonB[j] {
			outLabels = append(outLabels, l)
		}
	}
	return res, outLabels, nil
}

// Contract evaluates the whole network in the greedy contraction order and
// returns the resulting tensor together with its (sorted) free labels. It
// returns an error if the network is empty or invalid.
func (nw *Network) Contract() (*Tensor, []string, error) {
	order, _, err := nw.GreedyOrder()
	if err != nil {
		return nil, nil, err
	}
	return nw.contractWithOrder(order)
}

// ContractStep is one pairwise contraction in an order: the two operand indices
// (referring to the running list of intermediate tensors) that are combined.
type ContractStep struct {
	// A and B are indices into the current working set of tensors.
	A, B int
}

// contractWithOrder performs the pairwise contractions described by order and
// returns the final tensor and its labels.
func (nw *Network) contractWithOrder(order []ContractStep) (*Tensor, []string, error) {
	if len(nw.tensors) == 0 {
		return nil, nil, errors.New("tensornetwork: empty network")
	}
	if _, err := nw.validate(); err != nil {
		return nil, nil, err
	}
	tensors := make([]*Tensor, len(nw.tensors))
	labels := make([][]string, len(nw.labels))
	for i := range nw.tensors {
		tensors[i] = nw.tensors[i]
		labels[i] = append([]string(nil), nw.labels[i]...)
	}
	for _, step := range order {
		a, b := step.A, step.B
		if a < 0 || b < 0 || a >= len(tensors) || b >= len(tensors) || a == b || tensors[a] == nil || tensors[b] == nil {
			return nil, nil, errors.New("tensornetwork: invalid contraction step")
		}
		res, labs, err := contractPair(tensors[a], labels[a], tensors[b], labels[b])
		if err != nil {
			return nil, nil, err
		}
		tensors[a] = res
		labels[a] = labs
		tensors[b] = nil
		labels[b] = nil
	}
	// Combine any remaining (disconnected) pieces by outer product.
	var final *Tensor
	var finalLabels []string
	for i := range tensors {
		if tensors[i] == nil {
			continue
		}
		if final == nil {
			final = tensors[i]
			finalLabels = labels[i]
		} else {
			res, labs, err := contractPair(final, finalLabels, tensors[i], labels[i])
			if err != nil {
				return nil, nil, err
			}
			final = res
			finalLabels = labs
		}
	}
	return final, finalLabels, nil
}

// GreedyOrder returns a contraction order chosen greedily: at each step it
// contracts the pair whose contraction cost (the product of the dimensions of
// the union of their labels) is smallest. It returns the order and its total
// cost, or an error if the network is invalid.
func (nw *Network) GreedyOrder() ([]ContractStep, float64, error) {
	dim, err := nw.validate()
	if err != nil {
		return nil, 0, err
	}
	n := len(nw.tensors)
	if n == 0 {
		return nil, 0, errors.New("tensornetwork: empty network")
	}
	labelSets := make([]map[string]bool, n)
	alive := make([]bool, n)
	for i := 0; i < n; i++ {
		s := make(map[string]bool)
		for _, l := range nw.labels[i] {
			s[l] = true
		}
		labelSets[i] = s
		alive[i] = true
	}
	var order []ContractStep
	total := 0.0
	remaining := n
	for remaining > 1 {
		bestCost := math.Inf(1)
		bestA, bestB := -1, -1
		bestShared := -1
		for a := 0; a < n; a++ {
			if !alive[a] {
				continue
			}
			for b := a + 1; b < n; b++ {
				if !alive[b] {
					continue
				}
				shared := 0
				union := make(map[string]bool)
				for l := range labelSets[a] {
					union[l] = true
				}
				for l := range labelSets[b] {
					if labelSets[a][l] {
						shared++
					}
					union[l] = true
				}
				cost := 1.0
				for l := range union {
					cost *= float64(dim[l])
				}
				// Prefer connected (shared>0) low-cost contractions.
				better := cost < bestCost
				if shared > 0 && bestShared == 0 {
					better = true
				} else if shared == 0 && bestShared > 0 {
					better = false
				}
				if bestA == -1 || better {
					bestCost = cost
					bestA, bestB = a, b
					bestShared = shared
				}
			}
		}
		if bestA == -1 {
			break
		}
		total += bestCost
		// Merge b into a.
		newSet := make(map[string]bool)
		for l := range labelSets[bestA] {
			if !labelSets[bestB][l] {
				newSet[l] = true
			}
		}
		for l := range labelSets[bestB] {
			if !labelSets[bestA][l] {
				newSet[l] = true
			}
		}
		labelSets[bestA] = newSet
		alive[bestB] = false
		order = append(order, ContractStep{A: bestA, B: bestB})
		remaining--
	}
	return order, total, nil
}

// OptimalOrder returns a minimum-cost pairwise contraction order found by
// dynamic programming over subsets of tensors, together with the total cost.
// The cost of a contraction is the product of the dimensions of the union of the
// two operands' open labels. Because it enumerates subsets it is intended for
// networks of up to roughly twenty tensors. It returns an error if the network
// is empty, invalid, or too large.
func (nw *Network) OptimalOrder() ([]ContractStep, float64, error) {
	dim, err := nw.validate()
	if err != nil {
		return nil, 0, err
	}
	n := len(nw.tensors)
	if n == 0 {
		return nil, 0, errors.New("tensornetwork: empty network")
	}
	if n == 1 {
		return nil, 0, nil
	}
	if n > 22 {
		return nil, 0, errors.New("tensornetwork: too many tensors for exact ordering")
	}
	// Index distinct labels.
	labelIdx := make(map[string]int)
	var labelDim []float64
	for l, d := range dim {
		labelIdx[l] = len(labelDim)
		labelDim = append(labelDim, float64(d))
	}
	tensorMask := make([]uint64, n)
	for i := 0; i < n; i++ {
		var m uint64
		for _, l := range nw.labels[i] {
			m |= 1 << uint(labelIdx[l])
		}
		tensorMask[i] = m
	}
	full := 1 << uint(n)
	openMask := make([]uint64, full) // open labels of each tensor subset
	dp := make([]float64, full)      // best cost
	split := make([]int, full)       // chosen left subset
	for s := 1; s < full; s++ {
		// open labels = XOR of member tensor masks.
		var om uint64
		for i := 0; i < n; i++ {
			if s&(1<<uint(i)) != 0 {
				om ^= tensorMask[i]
			}
		}
		openMask[s] = om
		if bits.OnesCount(uint(s)) == 1 {
			dp[s] = 0
			continue
		}
		dp[s] = math.Inf(1)
	}
	prodMask := func(m uint64) float64 {
		p := 1.0
		for m != 0 {
			b := bits.TrailingZeros64(m)
			p *= labelDim[b]
			m &= m - 1
		}
		return p
	}
	for s := 1; s < full; s++ {
		if bits.OnesCount(uint(s)) < 2 {
			continue
		}
		// Enumerate proper non-empty submasks.
		for l := (s - 1) & s; l > 0; l = (l - 1) & s {
			r := s ^ l
			if l > r {
				continue // avoid double counting; ensure l<r canonical
			}
			cost := dp[l] + dp[r] + prodMask(openMask[l]|openMask[r])
			if cost < dp[s] {
				dp[s] = cost
				split[s] = l
			}
		}
	}
	var order []ContractStep
	// Reconstruct order: map each subset to a "result slot" = lowest tensor index.
	var build func(s int) int
	build = func(s int) int {
		if bits.OnesCount(uint(s)) == 1 {
			return bits.TrailingZeros(uint(s))
		}
		l := split[s]
		r := s ^ l
		ia := build(l)
		ib := build(r)
		if ia > ib {
			ia, ib = ib, ia
		}
		order = append(order, ContractStep{A: ia, B: ib})
		return ia
	}
	build(full - 1)
	return order, dp[full-1], nil
}

// ContractOptimal contracts the network following the order from
// [Network.OptimalOrder] and returns the resulting tensor and its labels.
func (nw *Network) ContractOptimal() (*Tensor, []string, error) {
	order, _, err := nw.OptimalOrder()
	if err != nil {
		return nil, nil, err
	}
	return nw.contractWithOrder(order)
}

// ContractionCost returns the total contraction cost of the given order under
// the same cost model as [Network.OptimalOrder]. It returns an error if the
// order is invalid for the network.
func (nw *Network) ContractionCost(order []ContractStep) (float64, error) {
	dim, err := nw.validate()
	if err != nil {
		return 0, err
	}
	n := len(nw.tensors)
	labelSets := make([]map[string]bool, n)
	alive := make([]bool, n)
	for i := 0; i < n; i++ {
		s := make(map[string]bool)
		for _, l := range nw.labels[i] {
			s[l] = true
		}
		labelSets[i] = s
		alive[i] = true
	}
	total := 0.0
	for _, step := range order {
		a, b := step.A, step.B
		if a < 0 || b < 0 || a >= n || b >= n || !alive[a] || !alive[b] || a == b {
			return 0, errors.New("tensornetwork: invalid step in cost computation")
		}
		union := make(map[string]bool)
		for l := range labelSets[a] {
			union[l] = true
		}
		for l := range labelSets[b] {
			union[l] = true
		}
		cost := 1.0
		for l := range union {
			cost *= float64(dim[l])
		}
		total += cost
		newSet := make(map[string]bool)
		for l := range labelSets[a] {
			if !labelSets[b][l] {
				newSet[l] = true
			}
		}
		for l := range labelSets[b] {
			if !labelSets[a][l] {
				newSet[l] = true
			}
		}
		labelSets[a] = newSet
		alive[b] = false
	}
	return total, nil
}
