package auctions

import (
	"errors"
	"sort"
)

// CombinatorialBid is a bid in a combinatorial auction: bidder Bidder offers
// Value for the exact bundle of items Items (item indices). Under XOR semantics
// each bidder wins at most one of their bids; under OR semantics a bidder may
// win several disjoint bids.
type CombinatorialBid struct {
	Bidder int
	Items  []int
	Value  float64
}

// itemsMask returns the bitmask of the bid's items and reports whether every
// item index is within [0, numItems).
func itemsMask(items []int, numItems int) (uint64, bool) {
	var mask uint64
	for _, it := range items {
		if it < 0 || it >= numItems {
			return 0, false
		}
		mask |= uint64(1) << uint(it)
	}
	return mask, true
}

// Allocation is the outcome of a combinatorial winner-determination problem: the
// indices (into the original bid slice) of the accepted bids and the total
// value they realise.
type Allocation struct {
	AcceptedBids []int
	Value        float64
}

// wdpState carries the immutable inputs of the branch-and-bound search.
type wdpState struct {
	masks    []uint64
	values   []float64
	bidders  []int
	numItems int
	xor      bool
}

// WinnerDetermination solves the combinatorial-auction winner-determination
// problem under XOR bidding (each bidder wins at most one bid): it selects a
// set of pairwise item-disjoint bids, no two from the same bidder, maximizing
// total value. It is solved exactly by branch-and-bound on the lowest-indexed
// unallocated item. numItems must be at most 63.
func WinnerDetermination(bids []CombinatorialBid, numItems int) (Allocation, error) {
	return winnerDetermination(bids, numItems, true)
}

// WinnerDeterminationOR solves the winner-determination problem under OR
// bidding, where a bidder may be awarded several of their disjoint bids. It
// selects pairwise item-disjoint bids maximizing total value.
func WinnerDeterminationOR(bids []CombinatorialBid, numItems int) (Allocation, error) {
	return winnerDetermination(bids, numItems, false)
}

func winnerDetermination(bids []CombinatorialBid, numItems int, xor bool) (Allocation, error) {
	if numItems < 0 || numItems > 63 {
		return Allocation{}, errors.New("auctions: numItems must be between 0 and 63")
	}
	st := wdpState{numItems: numItems, xor: xor}
	for _, b := range bids {
		mask, ok := itemsMask(b.Items, numItems)
		if !ok {
			return Allocation{}, errors.New("auctions: bid references an out-of-range item")
		}
		st.masks = append(st.masks, mask)
		st.values = append(st.values, b.Value)
		st.bidders = append(st.bidders, b.Bidder)
	}
	bestVal := 0.0
	var bestSel []int
	var chosen []int
	var used uint64
	usedBidder := map[int]bool{}
	var rec func(itemsUsed uint64, val float64)
	rec = func(itemsUsed uint64, val float64) {
		if val > bestVal+1e-12 {
			bestVal = val
			bestSel = append(bestSel[:0], chosen...)
		}
		// find lowest free item
		free := -1
		for it := 0; it < numItems; it++ {
			if itemsUsed&(uint64(1)<<uint(it)) == 0 {
				free = it
				break
			}
		}
		if free == -1 {
			return
		}
		freeBit := uint64(1) << uint(free)
		// Option: leave item `free` unallocated. To avoid revisiting it, mark it used.
		rec(itemsUsed|freeBit, val)
		// Option: allocate some bid covering `free`.
		for k := range st.masks {
			if st.masks[k]&freeBit == 0 {
				continue // does not cover the free item
			}
			if st.masks[k]&itemsUsed != 0 {
				continue // overlaps already-used items
			}
			if st.values[k] <= 0 {
				continue
			}
			if xor && usedBidder[st.bidders[k]] {
				continue
			}
			chosen = append(chosen, k)
			if xor {
				usedBidder[st.bidders[k]] = true
			}
			rec(itemsUsed|st.masks[k], val+st.values[k])
			if xor {
				usedBidder[st.bidders[k]] = false
			}
			chosen = chosen[:len(chosen)-1]
		}
	}
	_ = used
	rec(0, 0)
	sel := append([]int(nil), bestSel...)
	sort.Ints(sel)
	return Allocation{AcceptedBids: sel, Value: bestVal}, nil
}

// SocialWelfare returns the total value of an allocation over the given bids.
func SocialWelfare(bids []CombinatorialBid, alloc Allocation) float64 {
	var total float64
	for _, i := range alloc.AcceptedBids {
		total += bids[i].Value
	}
	return total
}

// VCGResult is the outcome of running the VCG mechanism on a combinatorial
// auction: the welfare-maximizing allocation, and for every bidder that wins at
// least one bid, their Clarke payment and resulting utility.
type VCGResult struct {
	Allocation Allocation
	// Payment maps a winning bidder to their VCG (Clarke pivot) payment.
	Payment map[int]float64
	// Utility maps a winning bidder to value obtained minus payment.
	Utility map[int]float64
	// Welfare is the total value of the chosen allocation.
	Welfare float64
}

// VCGMechanism runs the Vickrey-Clarke-Groves mechanism under XOR bidding: it
// computes the efficient allocation and charges each winner the externality it
// imposes on the others (the Clarke pivot rule). A bidder's payment equals the
// optimal welfare of the others when the bidder is absent minus the welfare the
// others obtain in the chosen allocation. Truthful bidding is a dominant
// strategy.
func VCGMechanism(bids []CombinatorialBid, numItems int) (VCGResult, error) {
	opt, err := WinnerDetermination(bids, numItems)
	if err != nil {
		return VCGResult{}, err
	}
	// value each bidder obtains in the optimal allocation
	valueOf := map[int]float64{}
	for _, i := range opt.AcceptedBids {
		valueOf[bids[i].Bidder] += bids[i].Value
	}
	res := VCGResult{
		Allocation: opt,
		Payment:    map[int]float64{},
		Utility:    map[int]float64{},
		Welfare:    opt.Value,
	}
	for bidder, gained := range valueOf {
		// welfare of the auction with this bidder's bids removed
		var without []CombinatorialBid
		for _, b := range bids {
			if b.Bidder != bidder {
				without = append(without, b)
			}
		}
		altOpt, err := WinnerDetermination(without, numItems)
		if err != nil {
			return VCGResult{}, err
		}
		welfareOthersInOpt := opt.Value - gained
		payment := altOpt.Value - welfareOthersInOpt
		if payment < 0 {
			payment = 0
		}
		res.Payment[bidder] = payment
		res.Utility[bidder] = gained - payment
	}
	return res, nil
}

// SingleItemToCombinatorial converts single-item sealed bids into combinatorial
// bids for the lone item 0, so that VCGMechanism reproduces the second-price
// (Vickrey) auction.
func SingleItemToCombinatorial(bids []Bid) []CombinatorialBid {
	out := make([]CombinatorialBid, len(bids))
	for i, b := range bids {
		out[i] = CombinatorialBid{Bidder: b.Bidder, Items: []int{0}, Value: b.Value}
	}
	return out
}
