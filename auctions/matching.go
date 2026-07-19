package auctions

import "errors"

// Matching is a two-sided assignment between men and women (or, generally, two
// equal-sized sides). ManToWoman[m] is the woman matched to man m, or -1 if he
// is unmatched; WomanToMan[w] is the symmetric map.
type Matching struct {
	ManToWoman []int
	WomanToMan []int
}

// validatePrefs checks that prefs is an n×n table whose rows are permutations
// of {0,...,n-1}.
func validatePrefs(prefs [][]int, n int) error {
	if len(prefs) != n {
		return errors.New("auctions: preference table has wrong number of rows")
	}
	for _, row := range prefs {
		if len(row) != n {
			return errors.New("auctions: preference row has wrong length")
		}
		seen := make([]bool, n)
		for _, v := range row {
			if v < 0 || v >= n || seen[v] {
				return errors.New("auctions: preference row is not a permutation")
			}
			seen[v] = true
		}
	}
	return nil
}

// rankTable returns ranks[i][j] = the position (0 = most preferred) of item j in
// agent i's preference list.
func rankTable(prefs [][]int, n int) [][]int {
	ranks := make([][]int, n)
	for i := 0; i < n; i++ {
		ranks[i] = make([]int, n)
		for pos, item := range prefs[i] {
			ranks[i][item] = pos
		}
	}
	return ranks
}

// galeShapley runs deferred acceptance with the proposing side's preferences
// proposerPrefs and the receiving side's preferences receiverPrefs. It returns,
// for each proposer, the receiver they are matched to.
func galeShapley(proposerPrefs, receiverPrefs [][]int, n int) []int {
	recRank := rankTable(receiverPrefs, n)
	next := make([]int, n)       // next receiver index each proposer will try
	partner := make([]int, n)    // receiver -> proposer, -1 if free
	proposerOf := make([]int, n) // proposer -> receiver, -1 if free
	for i := 0; i < n; i++ {
		partner[i] = -1
		proposerOf[i] = -1
	}
	free := make([]int, n)
	for i := range free {
		free[i] = i
	}
	for len(free) > 0 {
		m := free[len(free)-1]
		free = free[:len(free)-1]
		if next[m] >= n {
			continue
		}
		w := proposerPrefs[m][next[m]]
		next[m]++
		cur := partner[w]
		if cur == -1 {
			partner[w] = m
			proposerOf[m] = w
		} else if recRank[w][m] < recRank[w][cur] {
			partner[w] = m
			proposerOf[m] = w
			proposerOf[cur] = -1
			free = append(free, cur)
		} else {
			free = append(free, m)
		}
	}
	return proposerOf
}

// GaleShapley computes the man-optimal stable matching by the Gale-Shapley
// deferred-acceptance algorithm with men proposing. menPrefs[m] ranks the
// women and womenPrefs[w] ranks the men; both must be n×n permutation tables.
func GaleShapley(menPrefs, womenPrefs [][]int) (Matching, error) {
	n := len(menPrefs)
	if err := validatePrefs(menPrefs, n); err != nil {
		return Matching{}, err
	}
	if err := validatePrefs(womenPrefs, n); err != nil {
		return Matching{}, err
	}
	manToWoman := galeShapley(menPrefs, womenPrefs, n)
	womanToMan := make([]int, n)
	for i := range womanToMan {
		womanToMan[i] = -1
	}
	for m, w := range manToWoman {
		if w >= 0 {
			womanToMan[w] = m
		}
	}
	return Matching{ManToWoman: manToWoman, WomanToMan: womanToMan}, nil
}

// ManOptimalMatching is Gale-Shapley with men proposing; it returns the stable
// matching most preferred by every man.
func ManOptimalMatching(menPrefs, womenPrefs [][]int) (Matching, error) {
	return GaleShapley(menPrefs, womenPrefs)
}

// WomanOptimalMatching is Gale-Shapley with women proposing; it returns the
// stable matching most preferred by every woman.
func WomanOptimalMatching(menPrefs, womenPrefs [][]int) (Matching, error) {
	n := len(womenPrefs)
	if err := validatePrefs(menPrefs, n); err != nil {
		return Matching{}, err
	}
	if err := validatePrefs(womenPrefs, n); err != nil {
		return Matching{}, err
	}
	womanToMan := galeShapley(womenPrefs, menPrefs, n)
	manToWoman := make([]int, n)
	for i := range manToWoman {
		manToWoman[i] = -1
	}
	for w, m := range womanToMan {
		if m >= 0 {
			manToWoman[m] = w
		}
	}
	return Matching{ManToWoman: manToWoman, WomanToMan: womanToMan}, nil
}

// BlockingPairs returns every (man, woman) pair that blocks the matching: both
// strictly prefer each other to their assigned partner. A matching is stable
// iff this list is empty.
func BlockingPairs(menPrefs, womenPrefs [][]int, m Matching) [][2]int {
	n := len(menPrefs)
	menRank := rankTable(menPrefs, n)
	womenRank := rankTable(womenPrefs, n)
	var out [][2]int
	for man := 0; man < n; man++ {
		curW := m.ManToWoman[man]
		for w := 0; w < n; w++ {
			if w == curW {
				continue
			}
			// man prefers w to his current partner?
			if curW != -1 && menRank[man][w] >= menRank[man][curW] {
				continue
			}
			curM := m.WomanToMan[w]
			// woman prefers man to her current partner?
			if curM != -1 && womenRank[w][man] >= womenRank[w][curM] {
				continue
			}
			out = append(out, [2]int{man, w})
		}
	}
	return out
}

// CountBlockingPairs returns the number of blocking pairs of the matching.
func CountBlockingPairs(menPrefs, womenPrefs [][]int, m Matching) int {
	return len(BlockingPairs(menPrefs, womenPrefs, m))
}

// IsStableMatching reports whether the matching has no blocking pair.
func IsStableMatching(menPrefs, womenPrefs [][]int, m Matching) bool {
	return len(BlockingPairs(menPrefs, womenPrefs, m)) == 0
}

// IsPerfectMatching reports whether every man and woman is matched.
func IsPerfectMatching(m Matching) bool {
	for _, w := range m.ManToWoman {
		if w < 0 {
			return false
		}
	}
	for _, man := range m.WomanToMan {
		if man < 0 {
			return false
		}
	}
	return true
}

// TopTradingCycles runs Gale's top-trading-cycles algorithm for a housing
// market: agent i initially owns house i and prefs[i] ranks all houses. It
// returns the unique core allocation, assignment[i] being the house agent i
// receives. The mechanism is strategy-proof and Pareto-efficient.
func TopTradingCycles(prefs [][]int) ([]int, error) {
	n := len(prefs)
	if err := validatePrefs(prefs, n); err != nil {
		return nil, err
	}
	assignment := make([]int, n)
	for i := range assignment {
		assignment[i] = -1
	}
	houseAvailable := make([]bool, n)
	agentActive := make([]bool, n)
	for i := 0; i < n; i++ {
		houseAvailable[i] = true
		agentActive[i] = true
	}
	// owner[h] = agent currently owning house h (constant here: h).
	owner := make([]int, n)
	for h := range owner {
		owner[h] = h
	}
	remaining := n
	for remaining > 0 {
		// each active agent points to its favourite available house
		points := make([]int, n)
		for a := 0; a < n; a++ {
			if !agentActive[a] {
				points[a] = -1
				continue
			}
			for _, h := range prefs[a] {
				if houseAvailable[h] {
					points[a] = h
					break
				}
			}
		}
		// find a cycle: agent -> house -> owner(house) -> ...
		start := -1
		for a := 0; a < n; a++ {
			if agentActive[a] {
				start = a
				break
			}
		}
		seen := map[int]int{}
		cur := start
		step := 0
		cycleStart := -1
		for {
			if s, ok := seen[cur]; ok {
				cycleStart = s
				break
			}
			seen[cur] = step
			step++
			h := points[cur]
			cur = owner[h]
		}
		// reconstruct agents in the cycle
		order := make([]int, 0, len(seen))
		for a, s := range seen {
			if s >= cycleStart {
				order = append(order, a)
			}
		}
		for _, a := range order {
			h := points[a]
			assignment[a] = h
			houseAvailable[h] = false
			agentActive[a] = false
			remaining--
		}
	}
	return assignment, nil
}

// SerialDictatorship allocates houses by letting agents choose in the given
// priority order, each taking their most-preferred still-available house.
// order is a permutation of agent indices; prefs[i] ranks all houses. It
// returns assignment[i], the house agent i receives.
func SerialDictatorship(prefs [][]int, order []int) ([]int, error) {
	n := len(prefs)
	if err := validatePrefs(prefs, n); err != nil {
		return nil, err
	}
	if len(order) != n {
		return nil, errors.New("auctions: order must list every agent once")
	}
	seen := make([]bool, n)
	for _, a := range order {
		if a < 0 || a >= n || seen[a] {
			return nil, errors.New("auctions: order is not a permutation of agents")
		}
		seen[a] = true
	}
	assignment := make([]int, n)
	for i := range assignment {
		assignment[i] = -1
	}
	available := make([]bool, n)
	for i := range available {
		available[i] = true
	}
	for _, a := range order {
		for _, h := range prefs[a] {
			if available[h] {
				assignment[a] = h
				available[h] = false
				break
			}
		}
	}
	return assignment, nil
}
