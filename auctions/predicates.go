package auctions

// IsSuperadditive reports whether v(S∪T) >= v(S) + v(T) for all disjoint
// coalitions S and T.
func (g CoopGame) IsSuperadditive() bool {
	n := g.Players
	size := 1 << uint(n)
	for s := 0; s < size; s++ {
		// iterate over subsets t disjoint from s
		comp := (size - 1) &^ s
		for t := comp; ; t = (t - 1) & comp {
			if g.Value(Coalition(s|t)) < g.Value(Coalition(s))+g.Value(Coalition(t))-1e-9 {
				return false
			}
			if t == 0 {
				break
			}
		}
	}
	return true
}

// IsMonotone reports whether v(S) <= v(T) whenever S ⊆ T.
func (g CoopGame) IsMonotone() bool {
	n := g.Players
	size := 1 << uint(n)
	for t := 0; t < size; t++ {
		// iterate over subsets s of t
		for s := t; ; s = (s - 1) & t {
			if g.Value(Coalition(s)) > g.Value(Coalition(t))+1e-9 {
				return false
			}
			if s == 0 {
				break
			}
		}
	}
	return true
}

// IsConvex reports whether the game is convex (supermodular):
// v(S∪T) + v(S∩T) >= v(S) + v(T) for all coalitions S and T. Convex games have
// a non-empty core containing the Shapley value.
func (g CoopGame) IsConvex() bool {
	n := g.Players
	size := 1 << uint(n)
	for s := 0; s < size; s++ {
		for t := 0; t < size; t++ {
			lhs := g.Value(Coalition(s|t)) + g.Value(Coalition(s&t))
			rhs := g.Value(Coalition(s)) + g.Value(Coalition(t))
			if lhs < rhs-1e-9 {
				return false
			}
		}
	}
	return true
}

// IsEssential reports whether v(N) > Σ_i v({i}); an essential game admits a
// non-degenerate set of imputations.
func (g CoopGame) IsEssential() bool {
	var sum float64
	for i := 0; i < g.Players; i++ {
		sum += g.Value(SingletonCoalition(i))
	}
	return g.GrandValue() > sum+1e-9
}

// IsConstantSum reports whether v(S) + v(N\S) = v(N) for every coalition S.
func (g CoopGame) IsConstantSum() bool {
	n := g.Players
	size := 1 << uint(n)
	vN := g.GrandValue()
	for s := 0; s < size; s++ {
		comp := (size - 1) &^ s
		if !approxEqual(g.Value(Coalition(s))+g.Value(Coalition(comp)), vN, 1e-9) {
			return false
		}
	}
	return true
}

// IsSimple reports whether the game is a simple game: monotone with worths in
// {0,1} and v(N)=1.
func (g CoopGame) IsSimple() bool {
	if !approxEqual(g.GrandValue(), 1, 1e-9) {
		return false
	}
	n := g.Players
	size := 1 << uint(n)
	for m := 0; m < size; m++ {
		v := g.Value(Coalition(m))
		if !approxEqual(v, 0, 1e-9) && !approxEqual(v, 1, 1e-9) {
			return false
		}
	}
	return g.IsMonotone()
}

// IsZeroNormalized reports whether every singleton coalition is worth zero.
func (g CoopGame) IsZeroNormalized() bool {
	for i := 0; i < g.Players; i++ {
		if !approxEqual(g.Value(SingletonCoalition(i)), 0, 1e-9) {
			return false
		}
	}
	return true
}

// ZeroNormalized returns the strategically equivalent zero-normalized game
// obtained by subtracting each player's singleton worth: w(S) = v(S) - Σ_{i∈S}
// v({i}).
func (g CoopGame) ZeroNormalized() CoopGame {
	n := g.Players
	size := 1 << uint(n)
	vals := make([]float64, size)
	single := make([]float64, n)
	for i := 0; i < n; i++ {
		single[i] = g.Value(SingletonCoalition(i))
	}
	for m := 0; m < size; m++ {
		s := Coalition(m)
		x := g.Value(s)
		for i := 0; i < n; i++ {
			if s.Contains(i) {
				x -= single[i]
			}
		}
		vals[m] = x
	}
	return CoopGame{Players: n, values: vals}
}

// DualGame returns the dual game v*(S) = v(N) - v(N\S). The dual of a simple
// game encodes the "blocking" coalitions of the original.
func (g CoopGame) DualGame() CoopGame {
	n := g.Players
	size := 1 << uint(n)
	vals := make([]float64, size)
	vN := g.GrandValue()
	for m := 0; m < size; m++ {
		comp := (size - 1) &^ m
		vals[m] = vN - g.Value(Coalition(comp))
	}
	return CoopGame{Players: n, values: vals}
}

// IsDummy reports whether player i is a dummy: its marginal contribution equals
// its singleton worth for every coalition, so it adds nothing beyond v({i}).
func (g CoopGame) IsDummy(i int) bool {
	n := g.Players
	size := 1 << uint(n)
	vi := g.Value(SingletonCoalition(i))
	for m := 0; m < size; m++ {
		s := Coalition(m)
		if s.Contains(i) {
			continue
		}
		if !approxEqual(g.Value(s.Add(i))-g.Value(s), vi, 1e-9) {
			return false
		}
	}
	return true
}

// IsNullPlayer reports whether player i is a null player: its marginal
// contribution to every coalition is zero.
func (g CoopGame) IsNullPlayer(i int) bool {
	n := g.Players
	size := 1 << uint(n)
	for m := 0; m < size; m++ {
		s := Coalition(m)
		if s.Contains(i) {
			continue
		}
		if !approxEqual(g.Value(s.Add(i))-g.Value(s), 0, 1e-9) {
			return false
		}
	}
	return true
}

// AreSymmetric reports whether players i and j are interchangeable: swapping
// them leaves the worth of every coalition unchanged.
func (g CoopGame) AreSymmetric(i, j int) bool {
	if i == j {
		return true
	}
	n := g.Players
	size := 1 << uint(n)
	for m := 0; m < size; m++ {
		s := Coalition(m)
		if s.Contains(i) == s.Contains(j) {
			continue
		}
		swapped := s.Remove(i).Remove(j)
		if s.Contains(i) {
			swapped = swapped.Add(j)
		} else {
			swapped = swapped.Add(i)
		}
		if !approxEqual(g.Value(s), g.Value(swapped), 1e-9) {
			return false
		}
	}
	return true
}
