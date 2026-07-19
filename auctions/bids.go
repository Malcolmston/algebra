package auctions

import (
	"errors"
	"sort"
)

// Bid is a single sealed bid: the identity of the bidder and the value they
// bid. For truthful bidders in a Vickrey auction, Value equals the bidder's
// private valuation.
type Bid struct {
	Bidder int
	Value  float64
}

// AuctionOutcome is the result of a single-item sealed-bid auction: the winning
// bidder, the price they pay and the seller's revenue. Winner is -1 when the
// item goes unsold (for example when no bid meets the reserve).
type AuctionOutcome struct {
	Winner  int
	Price   float64
	Revenue float64
}

// ErrNoBids is returned when an auction is run without any bids.
var ErrNoBids = errors.New("auctions: no bids")

// sortedByValue returns the bids sorted by decreasing value, breaking ties by
// increasing bidder index for determinism.
func sortedByValue(bids []Bid) []Bid {
	cp := make([]Bid, len(bids))
	copy(cp, bids)
	sort.SliceStable(cp, func(i, j int) bool {
		if cp[i].Value != cp[j].Value {
			return cp[i].Value > cp[j].Value
		}
		return cp[i].Bidder < cp[j].Bidder
	})
	return cp
}

// HighestBid returns the largest bid value and its bidder. It returns
// ErrNoBids for an empty slice.
func HighestBid(bids []Bid) (Bid, error) {
	if len(bids) == 0 {
		return Bid{}, ErrNoBids
	}
	s := sortedByValue(bids)
	return s[0], nil
}

// SecondHighestBid returns the second-largest bid value and its bidder. If only
// one bid is present it returns that single bid.
func SecondHighestBid(bids []Bid) (Bid, error) {
	if len(bids) == 0 {
		return Bid{}, ErrNoBids
	}
	s := sortedByValue(bids)
	if len(s) == 1 {
		return s[0], nil
	}
	return s[1], nil
}

// FirstPriceAuction runs a first-price sealed-bid auction with an optional
// reserve price: the highest bidder wins and pays their own bid, provided it
// meets the reserve. Ties are broken toward the lower bidder index.
func FirstPriceAuction(bids []Bid, reserve float64) (AuctionOutcome, error) {
	if len(bids) == 0 {
		return AuctionOutcome{Winner: -1}, ErrNoBids
	}
	s := sortedByValue(bids)
	top := s[0]
	if top.Value < reserve {
		return AuctionOutcome{Winner: -1}, nil
	}
	return AuctionOutcome{Winner: top.Bidder, Price: top.Value, Revenue: top.Value}, nil
}

// SecondPriceAuction runs a second-price (Vickrey) sealed-bid auction with an
// optional reserve: the highest bidder wins and pays the larger of the reserve
// and the second-highest bid. In the truthful equilibrium bidding one's value
// is a dominant strategy.
func SecondPriceAuction(bids []Bid, reserve float64) (AuctionOutcome, error) {
	if len(bids) == 0 {
		return AuctionOutcome{Winner: -1}, ErrNoBids
	}
	s := sortedByValue(bids)
	top := s[0]
	if top.Value < reserve {
		return AuctionOutcome{Winner: -1}, nil
	}
	price := reserve
	if len(s) >= 2 && s[1].Value > price {
		price = s[1].Value
	}
	return AuctionOutcome{Winner: top.Bidder, Price: price, Revenue: price}, nil
}

// VickreyAuction is an alias for SecondPriceAuction, the single-item Vickrey
// mechanism.
func VickreyAuction(bids []Bid, reserve float64) (AuctionOutcome, error) {
	return SecondPriceAuction(bids, reserve)
}

// ThirdPriceAuction runs a third-price sealed-bid auction: the highest bidder
// wins and pays the third-highest bid (or the reserve, whichever is larger).
// When fewer than three bids are present the price falls back to the lowest
// available bid or the reserve.
func ThirdPriceAuction(bids []Bid, reserve float64) (AuctionOutcome, error) {
	if len(bids) == 0 {
		return AuctionOutcome{Winner: -1}, ErrNoBids
	}
	s := sortedByValue(bids)
	top := s[0]
	if top.Value < reserve {
		return AuctionOutcome{Winner: -1}, nil
	}
	price := reserve
	idx := 2
	if idx >= len(s) {
		idx = len(s) - 1
	}
	if s[idx].Value > price {
		price = s[idx].Value
	}
	return AuctionOutcome{Winner: top.Bidder, Price: price, Revenue: price}, nil
}

// AllPayAuction runs an all-pay auction: the highest bidder wins the item but
// every bidder pays their own bid, so the seller's revenue is the sum of all
// bids. The reserve determines whether the item is awarded; losing bidders pay
// regardless.
func AllPayAuction(bids []Bid, reserve float64) (AuctionOutcome, error) {
	if len(bids) == 0 {
		return AuctionOutcome{Winner: -1}, ErrNoBids
	}
	s := sortedByValue(bids)
	var revenue float64
	for _, b := range s {
		revenue += b.Value
	}
	top := s[0]
	if top.Value < reserve {
		return AuctionOutcome{Winner: -1, Revenue: revenue}, nil
	}
	return AuctionOutcome{Winner: top.Bidder, Price: top.Value, Revenue: revenue}, nil
}

// BidderUtility returns the payoff of the winner of a single-item auction given
// their private valuation and the outcome: valuation minus price if they won,
// otherwise zero. For an all-pay auction use AllPayUtility instead.
func BidderUtility(valuation float64, out AuctionOutcome, bidder int) float64 {
	if out.Winner == bidder {
		return valuation - out.Price
	}
	return 0
}

// AllPayUtility returns a bidder's payoff in an all-pay auction: their
// valuation minus their own bid if they won, otherwise minus their bid.
func AllPayUtility(valuation, ownBid float64, out AuctionOutcome, bidder int) float64 {
	if out.Winner == bidder {
		return valuation - ownBid
	}
	return -ownBid
}

// IsEfficientAllocation reports whether the auction outcome awards the item to a
// bidder holding the maximum valuation among the supplied valuations, indexed
// by bidder. An unsold item (Winner == -1) is efficient only when no valuation
// exceeds the reserve is not modelled here; the check compares against the
// highest valuation directly.
func IsEfficientAllocation(valuations map[int]float64, out AuctionOutcome) bool {
	if out.Winner < 0 {
		return false
	}
	best := valuations[out.Winner]
	for _, v := range valuations {
		if v > best+1e-12 {
			return false
		}
	}
	return true
}
