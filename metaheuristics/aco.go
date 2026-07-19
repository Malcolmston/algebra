package metaheuristics

import "math"

// City is a point in the plane used to build travelling-salesman instances.
type City struct {
	X, Y float64
}

// EuclideanDistance returns the straight-line distance between two cities.
func EuclideanDistance(a, b City) float64 {
	dx := a.X - b.X
	dy := a.Y - b.Y
	return math.Sqrt(dx*dx + dy*dy)
}

// DistanceMatrix builds the symmetric matrix of pairwise Euclidean distances
// for the given cities.
func DistanceMatrix(cities []City) [][]float64 {
	n := len(cities)
	d := make([][]float64, n)
	for i := range d {
		d[i] = make([]float64, n)
	}
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			dist := EuclideanDistance(cities[i], cities[j])
			d[i][j] = dist
			d[j][i] = dist
		}
	}
	return d
}

// TourLength returns the total length of a closed tour visiting the cities in
// the order given by the permutation tour, using the distance matrix dist and
// returning to the start.
func TourLength(tour []int, dist [][]float64) float64 {
	n := len(tour)
	if n < 2 {
		return 0
	}
	total := 0.0
	for i := 0; i < n; i++ {
		total += dist[tour[i]][tour[(i+1)%n]]
	}
	return total
}

// NearestNeighborTour builds a tour greedily starting from the given city,
// always moving to the nearest unvisited city. It returns the tour permutation.
func NearestNeighborTour(dist [][]float64, start int) []int {
	n := len(dist)
	visited := make([]bool, n)
	tour := make([]int, 0, n)
	cur := start
	visited[cur] = true
	tour = append(tour, cur)
	for len(tour) < n {
		next := -1
		best := math.Inf(1)
		for j := 0; j < n; j++ {
			if !visited[j] && dist[cur][j] < best {
				best = dist[cur][j]
				next = j
			}
		}
		if next < 0 {
			break
		}
		visited[next] = true
		tour = append(tour, next)
		cur = next
	}
	return tour
}

// TwoOptSwap returns a new tour with the segment between positions i and j
// (inclusive) reversed. This is the elementary move of 2-opt local search.
func TwoOptSwap(tour []int, i, j int) []int {
	out := append([]int(nil), tour...)
	for i < j {
		out[i], out[j] = out[j], out[i]
		i++
		j--
	}
	return out
}

// TwoOpt improves a tour by repeatedly applying the best 2-opt move until no
// move shortens the tour, returning the locally optimal tour and its length.
func TwoOpt(tour []int, dist [][]float64) ([]int, float64) {
	best := append([]int(nil), tour...)
	bestLen := TourLength(best, dist)
	n := len(best)
	improved := true
	for improved {
		improved = false
		for i := 1; i < n-1; i++ {
			for j := i + 1; j < n; j++ {
				cand := TwoOptSwap(best, i, j)
				l := TourLength(cand, dist)
				if l+1e-12 < bestLen {
					best = cand
					bestLen = l
					improved = true
				}
			}
		}
	}
	return best, bestLen
}

// ACOConfig configures [RunACO], ant colony optimization for the TSP.
type ACOConfig struct {
	// Ants is the number of ants per iteration.
	Ants int
	// Iterations is the number of iterations.
	Iterations int
	// Alpha is the pheromone influence exponent.
	Alpha float64
	// Beta is the heuristic (inverse-distance) influence exponent.
	Beta float64
	// Rho is the pheromone evaporation rate in [0, 1].
	Rho float64
	// Q is the pheromone deposit constant.
	Q float64
	// InitPheromone is the initial pheromone level on every edge.
	InitPheromone float64
	// Elitist, when true, adds an extra deposit along the best-so-far tour each
	// iteration (elitist ant system).
	Elitist bool
	// RecordHistory enables per-iteration best-length recording.
	RecordHistory bool
}

// DefaultACOConfig returns a reasonable configuration.
func DefaultACOConfig() ACOConfig {
	return ACOConfig{
		Ants:          30,
		Iterations:    200,
		Alpha:         1,
		Beta:          3,
		Rho:           0.5,
		Q:             1,
		InitPheromone: 1,
		Elitist:       true,
	}
}

func (c ACOConfig) validate() error {
	if c.Ants <= 0 || c.Iterations <= 0 {
		return ErrInvalidConfig
	}
	if c.Rho < 0 || c.Rho > 1 {
		return ErrInvalidConfig
	}
	return nil
}

// ACOResult holds the outcome of an ant colony run.
type ACOResult struct {
	// Tour is the best tour permutation found.
	Tour []int
	// Length is the length of Tour.
	Length float64
	// History holds the best length per iteration when recording is enabled.
	History []float64
}

// RunACO solves a travelling salesman instance given by its distance matrix
// using ant colony optimization. It is deterministic given rng.
func RunACO(dist [][]float64, cfg ACOConfig, rng *RNG) (ACOResult, error) {
	if err := cfg.validate(); err != nil {
		return ACOResult{}, err
	}
	n := len(dist)
	if n < 2 {
		return ACOResult{Tour: makeIdentity(n)}, nil
	}
	// Pheromone and heuristic matrices.
	pher := make([][]float64, n)
	heur := make([][]float64, n)
	for i := 0; i < n; i++ {
		pher[i] = make([]float64, n)
		heur[i] = make([]float64, n)
		for j := 0; j < n; j++ {
			pher[i][j] = cfg.InitPheromone
			if i != j && dist[i][j] > 0 {
				heur[i][j] = 1 / dist[i][j]
			}
		}
	}

	best := NearestNeighborTour(dist, 0)
	bestLen := TourLength(best, dist)
	var history []float64

	for it := 0; it < cfg.Iterations; it++ {
		iterTours := make([][]int, cfg.Ants)
		iterLens := make([]float64, cfg.Ants)
		for a := 0; a < cfg.Ants; a++ {
			start := rng.Intn(n)
			tour := antWalk(pher, heur, dist, start, cfg, rng)
			iterTours[a] = tour
			iterLens[a] = TourLength(tour, dist)
			if iterLens[a] < bestLen {
				bestLen = iterLens[a]
				best = append([]int(nil), tour...)
			}
		}
		// Evaporate.
		for i := 0; i < n; i++ {
			for j := 0; j < n; j++ {
				pher[i][j] *= (1 - cfg.Rho)
			}
		}
		// Deposit.
		for a := 0; a < cfg.Ants; a++ {
			depositTour(pher, iterTours[a], cfg.Q/iterLens[a])
		}
		if cfg.Elitist {
			depositTour(pher, best, cfg.Q/bestLen)
		}
		if cfg.RecordHistory {
			history = append(history, bestLen)
		}
	}
	return ACOResult{Tour: best, Length: bestLen, History: history}, nil
}

func makeIdentity(n int) []int {
	t := make([]int, n)
	for i := range t {
		t[i] = i
	}
	return t
}

func antWalk(pher, heur, dist [][]float64, start int, cfg ACOConfig, rng *RNG) []int {
	n := len(dist)
	visited := make([]bool, n)
	tour := make([]int, 0, n)
	cur := start
	visited[cur] = true
	tour = append(tour, cur)
	weights := make([]float64, n)
	for len(tour) < n {
		total := 0.0
		for j := 0; j < n; j++ {
			if visited[j] {
				weights[j] = 0
				continue
			}
			w := math.Pow(pher[cur][j], cfg.Alpha) * math.Pow(heur[cur][j], cfg.Beta)
			weights[j] = w
			total += w
		}
		next := -1
		if total <= 0 {
			for j := 0; j < n; j++ {
				if !visited[j] {
					next = j
					break
				}
			}
		} else {
			next = rng.Choice(weights)
			if visited[next] {
				for j := 0; j < n; j++ {
					if !visited[j] {
						next = j
						break
					}
				}
			}
		}
		visited[next] = true
		tour = append(tour, next)
		cur = next
	}
	return tour
}

func depositTour(pher [][]float64, tour []int, amount float64) {
	n := len(tour)
	for i := 0; i < n; i++ {
		a := tour[i]
		b := tour[(i+1)%n]
		pher[a][b] += amount
		pher[b][a] += amount
	}
}
