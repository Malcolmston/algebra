package gametheory

// ParetoDominates reports whether payoff vector a Pareto-dominates payoff
// vector b: a is at least as large as b in every coordinate and strictly larger
// in at least one. Both slices must have the same length. Larger is better in
// every objective.
func ParetoDominates(a, b []float64) bool {
	if len(a) != len(b) {
		return false
	}
	strict := false
	for i := range a {
		if a[i] < b[i] {
			return false
		}
		if a[i] > b[i] {
			strict = true
		}
	}
	return strict
}

// ParetoOptimalIndices returns the sorted indices of the Pareto-optimal points
// among the given points (each a vector of objectives, all maximized): a point
// is Pareto-optimal when no other point Pareto-dominates it. All points must
// share the same dimension.
func ParetoOptimalIndices(points [][]float64) []int {
	var out []int
	for i := range points {
		dominated := false
		for k := range points {
			if k != i && ParetoDominates(points[k], points[i]) {
				dominated = true
				break
			}
		}
		if !dominated {
			out = append(out, i)
		}
	}
	return out
}

// ParetoFrontier returns the pure-strategy profiles whose payoff pair
// (RowPayoff, ColPayoff) is Pareto-optimal: no other pure profile gives both
// players at least as much with at least one strictly better. The profiles are
// returned in lexicographic (row, column) order.
func (g Game) ParetoFrontier() []PureProfile {
	m, n := g.NumRowStrategies(), g.NumColStrategies()
	var out []PureProfile
	for i := 0; i < m; i++ {
		for j := 0; j < n; j++ {
			if g.IsParetoOptimal(i, j) {
				out = append(out, PureProfile{Row: i, Col: j})
			}
		}
	}
	return out
}

// IsParetoOptimal reports whether the pure-strategy profile (i, j) is
// Pareto-optimal among all pure profiles: no other profile yields both players
// at least as much payoff and at least one player strictly more.
func (g Game) IsParetoOptimal(i, j int) bool {
	ri, ci := g.Row[i][j], g.Col[i][j]
	m, n := g.NumRowStrategies(), g.NumColStrategies()
	for a := 0; a < m; a++ {
		for b := 0; b < n; b++ {
			if a == i && b == j {
				continue
			}
			if g.Row[a][b] >= ri && g.Col[a][b] >= ci &&
				(g.Row[a][b] > ri || g.Col[a][b] > ci) {
				return false
			}
		}
	}
	return true
}

// SocialOptimum returns the pure-strategy profile maximizing the sum of the two
// players' payoffs (the utilitarian social welfare) and that maximal sum. Ties
// are broken by lexicographic (row, column) order.
func (g Game) SocialOptimum() (PureProfile, float64) {
	m, n := g.NumRowStrategies(), g.NumColStrategies()
	best := PureProfile{Row: 0, Col: 0}
	bestVal := g.Row[0][0] + g.Col[0][0]
	for i := 0; i < m; i++ {
		for j := 0; j < n; j++ {
			s := g.Row[i][j] + g.Col[i][j]
			if s > bestVal {
				bestVal = s
				best = PureProfile{Row: i, Col: j}
			}
		}
	}
	return best, bestVal
}
