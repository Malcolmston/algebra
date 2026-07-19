package auctions

import "sort"

// Coalition is a subset of players encoded as a bitmask: bit i is set when
// player i belongs to the coalition. Player indices start at 0.
type Coalition uint64

// popcount returns the number of set bits in x.
func popcount(x uint64) int {
	n := 0
	for x != 0 {
		x &= x - 1
		n++
	}
	return n
}

// Contains reports whether player i belongs to the coalition.
func (c Coalition) Contains(i int) bool { return c&(1<<uint(i)) != 0 }

// Add returns the coalition with player i added.
func (c Coalition) Add(i int) Coalition { return c | 1<<uint(i) }

// Remove returns the coalition with player i removed.
func (c Coalition) Remove(i int) Coalition { return c &^ (1 << uint(i)) }

// Union returns the coalition containing every player in c or in other.
func (c Coalition) Union(other Coalition) Coalition { return c | other }

// Intersect returns the coalition of players in both c and other.
func (c Coalition) Intersect(other Coalition) Coalition { return c & other }

// Difference returns the players in c that are not in other.
func (c Coalition) Difference(other Coalition) Coalition { return c &^ other }

// Complement returns the coalition of all players in {0,...,n-1} not in c.
func (c Coalition) Complement(n int) Coalition {
	full := Coalition((uint64(1) << uint(n)) - 1)
	return full &^ c
}

// Size returns the number of players in the coalition.
func (c Coalition) Size() int { return popcount(uint64(c)) }

// IsEmpty reports whether the coalition has no members.
func (c Coalition) IsEmpty() bool { return c == 0 }

// IsSubsetOf reports whether every member of c is also a member of other.
func (c Coalition) IsSubsetOf(other Coalition) bool { return c&other == c }

// Members returns the sorted player indices contained in the coalition, given
// a universe of n players.
func (c Coalition) Members(n int) []int {
	out := make([]int, 0, c.Size())
	for i := 0; i < n; i++ {
		if c.Contains(i) {
			out = append(out, i)
		}
	}
	return out
}

// FullCoalition returns the grand coalition {0,...,n-1}.
func FullCoalition(n int) Coalition { return Coalition((uint64(1) << uint(n)) - 1) }

// EmptyCoalition returns the empty coalition.
func EmptyCoalition() Coalition { return 0 }

// SingletonCoalition returns the coalition {i}.
func SingletonCoalition(i int) Coalition { return Coalition(1) << uint(i) }

// CoalitionFromMembers builds a coalition from an explicit list of player
// indices. Duplicate indices are collapsed.
func CoalitionFromMembers(members []int) Coalition {
	var c Coalition
	for _, i := range members {
		c = c.Add(i)
	}
	return c
}

// AllCoalitions returns every one of the 2^n coalitions over n players, in
// increasing bitmask order (starting with the empty coalition).
func AllCoalitions(n int) []Coalition {
	size := 1 << uint(n)
	out := make([]Coalition, size)
	for m := 0; m < size; m++ {
		out[m] = Coalition(m)
	}
	return out
}

// ProperSubsets returns every non-empty proper subset of the coalition c over a
// universe of n players (excluding the empty set and c itself), sorted by
// bitmask.
func ProperSubsets(c Coalition, n int) []Coalition {
	members := c.Members(n)
	k := len(members)
	var out []Coalition
	for mask := 1; mask < (1<<uint(k))-1; mask++ {
		var s Coalition
		for j := 0; j < k; j++ {
			if mask&(1<<uint(j)) != 0 {
				s = s.Add(members[j])
			}
		}
		out = append(out, s)
	}
	sort.Slice(out, func(a, b int) bool { return out[a] < out[b] })
	return out
}
