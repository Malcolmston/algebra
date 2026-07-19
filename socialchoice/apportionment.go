package socialchoice

import (
	"errors"
	"math"
	"sort"
)

// DivisorFunc maps the number of seats a party currently holds to the divisor
// applied to its votes when computing the next-seat priority in a highest-
// averages apportionment. A divisor of zero denotes an infinite priority, giving
// every party with positive votes a seat before any second seat is awarded.
type DivisorFunc func(seatsWon int) float64

// HighestAverages allocates seats among parties by the generic highest-averages
// (divisor) method: it repeatedly awards the next seat to the party with the
// greatest votes/divisor(seatsWon) priority, breaking ties by votes and then by
// the lower party index. votes must be non-negative and seats non-negative.
func HighestAverages(votes []int, seats int, divisor DivisorFunc) ([]int, error) {
	if seats < 0 {
		return nil, errors.New("socialchoice: seats must be non-negative")
	}
	for _, v := range votes {
		if v < 0 {
			return nil, errors.New("socialchoice: votes must be non-negative")
		}
	}
	alloc := make([]int, len(votes))
	for s := 0; s < seats; s++ {
		bestParty := -1
		bestPri := math.Inf(-1)
		bestVotes := -1
		for i, v := range votes {
			d := divisor(alloc[i])
			var pri float64
			if d == 0 {
				if v > 0 {
					pri = math.Inf(1)
				} else {
					pri = math.Inf(-1)
				}
			} else {
				pri = float64(v) / d
			}
			if pri > bestPri || (pri == bestPri && v > bestVotes) {
				bestParty, bestPri, bestVotes = i, pri, v
			}
		}
		if bestParty == -1 {
			break
		}
		alloc[bestParty]++
	}
	return alloc, nil
}

// DHondt allocates seats by the D'Hondt (Jefferson) method, divisor s+1.
func DHondt(votes []int, seats int) ([]int, error) {
	return HighestAverages(votes, seats, func(s int) float64 { return float64(s + 1) })
}

// Jefferson is an alias for DHondt.
func Jefferson(votes []int, seats int) ([]int, error) { return DHondt(votes, seats) }

// SainteLague allocates seats by the Sainte-Laguë (Webster) method, divisor
// 2s+1 (equivalently the odd-number divisors 1, 3, 5, …).
func SainteLague(votes []int, seats int) ([]int, error) {
	return HighestAverages(votes, seats, func(s int) float64 { return float64(2*s + 1) })
}

// Webster is an alias for SainteLague.
func Webster(votes []int, seats int) ([]int, error) { return SainteLague(votes, seats) }

// ModifiedSainteLague allocates seats by the modified Sainte-Laguë method, whose
// first divisor is 1.4 and whose later divisors are 3, 5, 7, ….
func ModifiedSainteLague(votes []int, seats int) ([]int, error) {
	return HighestAverages(votes, seats, func(s int) float64 {
		if s == 0 {
			return 1.4
		}
		return float64(2*s + 1)
	})
}

// Adams allocates seats by the Adams method, divisor s, so every party with
// positive votes receives at least one seat when seats permit.
func Adams(votes []int, seats int) ([]int, error) {
	return HighestAverages(votes, seats, func(s int) float64 { return float64(s) })
}

// Dean allocates seats by the Dean method, whose divisor is the harmonic mean of
// s and s+1.
func Dean(votes []int, seats int) ([]int, error) {
	return HighestAverages(votes, seats, func(s int) float64 {
		if s == 0 {
			return 0
		}
		return 2 * float64(s) * float64(s+1) / float64(2*s+1)
	})
}

// HuntingtonHill allocates seats by the Huntington–Hill (equal proportions)
// method, whose divisor is the geometric mean sqrt(s(s+1)); every party with
// positive votes receives at least one seat when seats permit.
func HuntingtonHill(votes []int, seats int) ([]int, error) {
	return HighestAverages(votes, seats, func(s int) float64 {
		return math.Sqrt(float64(s) * float64(s+1))
	})
}

// LargestRemainder allocates seats by the largest-remainder (quota) method for a
// given quota: each party first receives the floor of its votes/quota, then any
// remaining seats are handed out by descending fractional remainder (surplus
// seats are removed by ascending remainder when the floors overshoot). Ties are
// broken toward the lower party index.
func LargestRemainder(votes []int, seats int, quota float64) ([]int, error) {
	if seats < 0 {
		return nil, errors.New("socialchoice: seats must be non-negative")
	}
	if quota <= 0 {
		return nil, errors.New("socialchoice: quota must be positive")
	}
	for _, v := range votes {
		if v < 0 {
			return nil, errors.New("socialchoice: votes must be non-negative")
		}
	}
	n := len(votes)
	alloc := make([]int, n)
	rem := make([]float64, n)
	assigned := 0
	for i, v := range votes {
		q := float64(v) / quota
		alloc[i] = int(math.Floor(q))
		rem[i] = q - float64(alloc[i])
		assigned += alloc[i]
	}
	idx := make([]int, n)
	for i := range idx {
		idx[i] = i
	}
	if assigned < seats {
		sort.SliceStable(idx, func(a, b int) bool {
			if rem[idx[a]] != rem[idx[b]] {
				return rem[idx[a]] > rem[idx[b]]
			}
			return idx[a] < idx[b]
		})
		for k := 0; assigned < seats && k < n; k++ {
			alloc[idx[k]]++
			assigned++
		}
	} else if assigned > seats {
		sort.SliceStable(idx, func(a, b int) bool {
			if rem[idx[a]] != rem[idx[b]] {
				return rem[idx[a]] < rem[idx[b]]
			}
			return idx[a] < idx[b]
		})
		for k := 0; assigned > seats && k < n; k++ {
			if alloc[idx[k]] > 0 {
				alloc[idx[k]]--
				assigned--
			}
		}
	}
	return alloc, nil
}

// Hamilton allocates seats by Hamilton's method, the largest-remainder method
// with the Hare quota total/seats.
func Hamilton(votes []int, seats int) ([]int, error) {
	total := sumInts(votes)
	return LargestRemainder(votes, seats, HareQuota(total, seats))
}

// DroopLargestRemainder allocates seats by the largest-remainder method with the
// Hagenbach–Bischoff (exact Droop) quota total/(seats+1).
func DroopLargestRemainder(votes []int, seats int) ([]int, error) {
	total := sumInts(votes)
	return LargestRemainder(votes, seats, HagenbachBischoffQuota(total, seats))
}

// sumInts returns the sum of a slice of ints.
func sumInts(xs []int) int {
	s := 0
	for _, x := range xs {
		s += x
	}
	return s
}

// ExactQuotas returns each party's exact seat entitlement votes*seats/total.
func ExactQuotas(votes []int, seats int) []float64 {
	total := sumInts(votes)
	q := make([]float64, len(votes))
	if total == 0 {
		return q
	}
	for i, v := range votes {
		q[i] = float64(v) * float64(seats) / float64(total)
	}
	return q
}

// LowerQuotas returns the floor of each party's exact seat entitlement.
func LowerQuotas(votes []int, seats int) []int {
	q := ExactQuotas(votes, seats)
	out := make([]int, len(q))
	for i := range q {
		out[i] = int(math.Floor(q[i]))
	}
	return out
}

// UpperQuotas returns the ceiling of each party's exact seat entitlement.
func UpperQuotas(votes []int, seats int) []int {
	q := ExactQuotas(votes, seats)
	out := make([]int, len(q))
	for i := range q {
		out[i] = int(math.Ceil(q[i]))
	}
	return out
}

// SatisfiesQuota reports whether every party's allocation lies within its lower
// and upper quota, i.e. the allocation obeys the quota rule.
func SatisfiesQuota(alloc, votes []int, seats int) bool {
	if len(alloc) != len(votes) {
		return false
	}
	lo := LowerQuotas(votes, seats)
	hi := UpperQuotas(votes, seats)
	for i := range alloc {
		if alloc[i] < lo[i] || alloc[i] > hi[i] {
			return false
		}
	}
	return true
}

// ViolatesLowerQuota reports whether any party receives fewer seats than the
// floor of its exact entitlement.
func ViolatesLowerQuota(alloc, votes []int, seats int) bool {
	lo := LowerQuotas(votes, seats)
	for i := range alloc {
		if i < len(lo) && alloc[i] < lo[i] {
			return true
		}
	}
	return false
}

// ViolatesUpperQuota reports whether any party receives more seats than the
// ceiling of its exact entitlement.
func ViolatesUpperQuota(alloc, votes []int, seats int) bool {
	hi := UpperQuotas(votes, seats)
	for i := range alloc {
		if i < len(hi) && alloc[i] > hi[i] {
			return true
		}
	}
	return false
}

// HasAlabamaParadox reports whether increasing the house size from seats to
// seats+1 causes any party to lose a seat under Hamilton's method — the Alabama
// paradox.
func HasAlabamaParadox(votes []int, seats int) bool {
	a, err := Hamilton(votes, seats)
	if err != nil {
		return false
	}
	b, err := Hamilton(votes, seats+1)
	if err != nil {
		return false
	}
	for i := range a {
		if b[i] < a[i] {
			return true
		}
	}
	return false
}

// IsHouseMonotone reports whether the given divisor method is free of the
// Alabama paradox across house sizes 1..maxSeats for these votes; every divisor
// method is house-monotone, so this serves as a verification helper.
func IsHouseMonotone(votes []int, maxSeats int, divisor DivisorFunc) bool {
	prev, err := HighestAverages(votes, 0, divisor)
	if err != nil {
		return false
	}
	for s := 1; s <= maxSeats; s++ {
		cur, err := HighestAverages(votes, s, divisor)
		if err != nil {
			return false
		}
		for i := range cur {
			if cur[i] < prev[i] {
				return false
			}
		}
		prev = cur
	}
	return true
}
