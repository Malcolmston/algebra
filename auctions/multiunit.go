package auctions

import "errors"

// MultiUnitOutcome is the result of a multi-unit auction in which every winning
// bidder demands a single identical unit. Winners lists the winning bidder
// indices; Payments is the price each winner pays, aligned with Winners.
type MultiUnitOutcome struct {
	Winners  []int
	Payments []float64
	Revenue  float64
}

// UniformPriceAuction sells units identical goods under single-unit demand: the
// units highest bidders each win one unit and all pay the same clearing price,
// the highest losing bid (the (units+1)-th highest bid), or 0 if there is no
// losing bid. Ties are broken toward lower bidder indices.
func UniformPriceAuction(bids []Bid, units int) (MultiUnitOutcome, error) {
	if units <= 0 {
		return MultiUnitOutcome{}, errors.New("auctions: units must be positive")
	}
	if len(bids) == 0 {
		return MultiUnitOutcome{}, ErrNoBids
	}
	s := sortedByValue(bids)
	k := units
	if k > len(s) {
		k = len(s)
	}
	price := 0.0
	if units < len(s) {
		price = s[units].Value
	}
	out := MultiUnitOutcome{}
	for i := 0; i < k; i++ {
		out.Winners = append(out.Winners, s[i].Bidder)
		out.Payments = append(out.Payments, price)
		out.Revenue += price
	}
	return out, nil
}

// PayAsBidAuction (discriminatory auction) sells units identical goods under
// single-unit demand: the units highest bidders each win one unit and pay their
// own bid.
func PayAsBidAuction(bids []Bid, units int) (MultiUnitOutcome, error) {
	if units <= 0 {
		return MultiUnitOutcome{}, errors.New("auctions: units must be positive")
	}
	if len(bids) == 0 {
		return MultiUnitOutcome{}, ErrNoBids
	}
	s := sortedByValue(bids)
	k := units
	if k > len(s) {
		k = len(s)
	}
	out := MultiUnitOutcome{}
	for i := 0; i < k; i++ {
		out.Winners = append(out.Winners, s[i].Bidder)
		out.Payments = append(out.Payments, s[i].Value)
		out.Revenue += s[i].Value
	}
	return out, nil
}

// MultiUnitVickrey runs the multi-unit Vickrey (VCG) auction under single-unit
// demand: the units highest bidders each win one unit, and each winner pays the
// highest losing bid that would have won in their absence. With one unit this
// reduces to the classic second-price auction.
func MultiUnitVickrey(bids []Bid, units int) (MultiUnitOutcome, error) {
	if units <= 0 {
		return MultiUnitOutcome{}, errors.New("auctions: units must be positive")
	}
	if len(bids) == 0 {
		return MultiUnitOutcome{}, ErrNoBids
	}
	s := sortedByValue(bids)
	k := units
	if k > len(s) {
		k = len(s)
	}
	out := MultiUnitOutcome{}
	// Under single-unit demand every winner's Clarke payment equals the highest
	// losing bid, s[units] (the (units+1)-th highest), or 0 if none exists.
	price := 0.0
	if units < len(s) {
		price = s[units].Value
	}
	for i := 0; i < k; i++ {
		out.Winners = append(out.Winners, s[i].Bidder)
		out.Payments = append(out.Payments, price)
		out.Revenue += price
	}
	return out, nil
}

// TotalSurplus returns the total value realised by an allocation given the
// winners' valuations indexed by bidder.
func TotalSurplus(winners []int, valuations map[int]float64) float64 {
	var total float64
	for _, w := range winners {
		total += valuations[w]
	}
	return total
}
