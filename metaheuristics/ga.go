package metaheuristics

import (
	"math"
	"sort"
)

// Individual is one member of a genetic-algorithm population: a real-valued
// genome together with its cached objective (fitness) value.
type Individual struct {
	Genes   []float64
	Fitness float64
}

// Clone returns a deep copy of the individual.
func (ind Individual) Clone() Individual {
	return Individual{Genes: VecCopy(ind.Genes), Fitness: ind.Fitness}
}

// Population is a slice of individuals.
type Population []Individual

// SortByFitness sorts the population in ascending order of fitness (best, i.e.
// smallest, first) in place.
func (p Population) SortByFitness() {
	sort.Slice(p, func(i, j int) bool { return p[i].Fitness < p[j].Fitness })
}

// Best returns a copy of the fittest (minimum-fitness) individual. It returns a
// zero Individual for an empty population.
func (p Population) Best() Individual {
	if len(p) == 0 {
		return Individual{}
	}
	best := 0
	for i := 1; i < len(p); i++ {
		if p[i].Fitness < p[best].Fitness {
			best = i
		}
	}
	return p[best].Clone()
}

// MeanFitness returns the arithmetic mean of the population's fitness values.
func (p Population) MeanFitness() float64 {
	if len(p) == 0 {
		return 0
	}
	s := 0.0
	for _, ind := range p {
		s += ind.Fitness
	}
	return s / float64(len(p))
}

// TournamentSelect returns the index of the winner (lowest fitness) of a
// tournament among k uniformly-random members of the population.
func TournamentSelect(p Population, k int, rng *RNG) int {
	if k < 1 {
		k = 1
	}
	best := rng.Intn(len(p))
	for i := 1; i < k; i++ {
		c := rng.Intn(len(p))
		if p[c].Fitness < p[best].Fitness {
			best = c
		}
	}
	return best
}

// RouletteSelect returns an index chosen with probability proportional to
// fitness after transforming the minimization objective into non-negative
// weights (worst individual gets weight 0). It falls back to uniform selection
// when the population is uniform.
func RouletteSelect(p Population, rng *RNG) int {
	worst := math.Inf(-1)
	for _, ind := range p {
		if ind.Fitness > worst {
			worst = ind.Fitness
		}
	}
	weights := make([]float64, len(p))
	total := 0.0
	for i, ind := range p {
		w := worst - ind.Fitness
		weights[i] = w
		total += w
	}
	if total <= 0 {
		return rng.Intn(len(p))
	}
	return rng.Choice(weights)
}

// RankSelect returns an index using linear rank selection with selection
// pressure sp in [1, 2]: the best individual receives expected count sp and the
// worst 2-sp. It is robust to fitness scaling.
func RankSelect(p Population, sp float64, rng *RNG) int {
	n := len(p)
	if n == 1 {
		return 0
	}
	if sp < 1 {
		sp = 1
	}
	if sp > 2 {
		sp = 2
	}
	order := make([]int, n)
	for i := range order {
		order[i] = i
	}
	sort.Slice(order, func(a, b int) bool { return p[order[a]].Fitness < p[order[b]].Fitness })
	// rank 0 is best. Weight decreases with rank.
	weights := make([]float64, n)
	for r := 0; r < n; r++ {
		// linear ranking probability
		weights[order[r]] = (sp - (2*sp-2)*float64(r)/float64(n-1))
	}
	return rng.Choice(weights)
}

// StochasticUniversalSampling returns m indices sampled from the population in a
// single sweep with minimum spread, using non-negative weights derived from the
// minimization fitness. It is a lower-variance alternative to repeated roulette.
func StochasticUniversalSampling(p Population, m int, rng *RNG) []int {
	n := len(p)
	worst := math.Inf(-1)
	for _, ind := range p {
		if ind.Fitness > worst {
			worst = ind.Fitness
		}
	}
	weights := make([]float64, n)
	total := 0.0
	for i, ind := range p {
		w := worst - ind.Fitness + 1e-12
		weights[i] = w
		total += w
	}
	out := make([]int, 0, m)
	if total <= 0 {
		for i := 0; i < m; i++ {
			out = append(out, rng.Intn(n))
		}
		return out
	}
	step := total / float64(m)
	start := rng.Float64() * step
	acc := 0.0
	idx := 0
	for i := 0; i < m; i++ {
		pointer := start + float64(i)*step
		for acc+weights[idx] < pointer && idx < n-1 {
			acc += weights[idx]
			idx++
		}
		out = append(out, idx)
	}
	return out
}

// OnePointCrossover recombines two parent genomes at a single random cut point,
// returning two children.
func OnePointCrossover(a, b []float64, rng *RNG) ([]float64, []float64) {
	n := minInt(len(a), len(b))
	c1 := VecCopy(a)
	c2 := VecCopy(b)
	if n < 2 {
		return c1, c2
	}
	cut := 1 + rng.Intn(n-1)
	for i := cut; i < n; i++ {
		c1[i], c2[i] = b[i], a[i]
	}
	return c1, c2
}

// TwoPointCrossover recombines two parents by swapping the segment between two
// random cut points.
func TwoPointCrossover(a, b []float64, rng *RNG) ([]float64, []float64) {
	n := minInt(len(a), len(b))
	c1 := VecCopy(a)
	c2 := VecCopy(b)
	if n < 2 {
		return c1, c2
	}
	i := rng.Intn(n)
	j := rng.Intn(n)
	if i > j {
		i, j = j, i
	}
	for k := i; k <= j; k++ {
		c1[k], c2[k] = b[k], a[k]
	}
	return c1, c2
}

// UniformCrossover recombines two parents by independently swapping each gene
// with probability p.
func UniformCrossover(a, b []float64, p float64, rng *RNG) ([]float64, []float64) {
	n := minInt(len(a), len(b))
	c1 := VecCopy(a)
	c2 := VecCopy(b)
	for i := 0; i < n; i++ {
		if rng.Float64() < p {
			c1[i], c2[i] = b[i], a[i]
		}
	}
	return c1, c2
}

// ArithmeticCrossover returns the two convex combinations
// t*a+(1-t)*b and (1-t)*a+t*b for a random blending factor t in [0,1].
func ArithmeticCrossover(a, b []float64, rng *RNG) ([]float64, []float64) {
	n := minInt(len(a), len(b))
	t := rng.Float64()
	c1 := make([]float64, n)
	c2 := make([]float64, n)
	for i := 0; i < n; i++ {
		c1[i] = t*a[i] + (1-t)*b[i]
		c2[i] = (1-t)*a[i] + t*b[i]
	}
	return c1, c2
}

// BlendCrossover implements BLX-alpha: each child gene is drawn uniformly from
// an interval extended by alpha times the parents' gene spread on both sides.
func BlendCrossover(a, b []float64, alpha float64, rng *RNG) ([]float64, []float64) {
	n := minInt(len(a), len(b))
	c1 := make([]float64, n)
	c2 := make([]float64, n)
	for i := 0; i < n; i++ {
		lo := math.Min(a[i], b[i])
		hi := math.Max(a[i], b[i])
		d := hi - lo
		c1[i] = rng.Float64Range(lo-alpha*d, hi+alpha*d)
		c2[i] = rng.Float64Range(lo-alpha*d, hi+alpha*d)
	}
	return c1, c2
}

// SBXCrossover implements simulated binary crossover with distribution index
// eta (larger eta yields children closer to their parents).
func SBXCrossover(a, b []float64, eta float64, rng *RNG) ([]float64, []float64) {
	n := minInt(len(a), len(b))
	c1 := make([]float64, n)
	c2 := make([]float64, n)
	for i := 0; i < n; i++ {
		u := rng.Float64()
		var beta float64
		if u <= 0.5 {
			beta = math.Pow(2*u, 1/(eta+1))
		} else {
			beta = math.Pow(1/(2*(1-u)), 1/(eta+1))
		}
		c1[i] = 0.5 * ((1+beta)*a[i] + (1-beta)*b[i])
		c2[i] = 0.5 * ((1-beta)*a[i] + (1+beta)*b[i])
	}
	return c1, c2
}

// GaussianMutation adds independent zero-mean Gaussian noise to each gene with
// probability rate; the standard deviation is sigma times the box width. The
// result is clipped to the box.
func GaussianMutation(g []float64, rate, sigma float64, b Bounds, rng *RNG) []float64 {
	width := b.Width()
	out := VecCopy(g)
	for i := range out {
		if rng.Float64() < rate {
			out[i] += rng.NormFloat64() * sigma * width[i]
		}
	}
	return b.ClipInPlace(out)
}

// UniformMutation resets each gene, with probability rate, to a uniform value
// within its box coordinate.
func UniformMutation(g []float64, rate float64, b Bounds, rng *RNG) []float64 {
	out := VecCopy(g)
	for i := range out {
		if rng.Float64() < rate {
			out[i] = rng.Float64Range(b.Lower[i], b.Upper[i])
		}
	}
	return out
}

// PolynomialMutation applies the polynomial mutation operator with distribution
// index eta to each gene with probability rate, respecting the box.
func PolynomialMutation(g []float64, rate, eta float64, b Bounds, rng *RNG) []float64 {
	out := VecCopy(g)
	for i := range out {
		if rng.Float64() >= rate {
			continue
		}
		lo, hi := b.Lower[i], b.Upper[i]
		if hi == lo {
			continue
		}
		x := out[i]
		d1 := (x - lo) / (hi - lo)
		d2 := (hi - x) / (hi - lo)
		u := rng.Float64()
		var dq float64
		mp := 1 / (eta + 1)
		if u < 0.5 {
			dq = math.Pow(2*u+(1-2*u)*math.Pow(1-d1, eta+1), mp) - 1
		} else {
			dq = 1 - math.Pow(2*(1-u)+2*(u-0.5)*math.Pow(1-d2, eta+1), mp)
		}
		out[i] = Clamp(x+dq*(hi-lo), lo, hi)
	}
	return out
}

// BitFlipMutation flips each bit of a binary genome (values 0 or 1) with
// probability rate, returning a new genome.
func BitFlipMutation(g []int, rate float64, rng *RNG) []int {
	out := append([]int(nil), g...)
	for i := range out {
		if rng.Float64() < rate {
			out[i] = 1 - out[i]
		}
	}
	return out
}

// GAConfig configures the real-coded genetic algorithm [RunGA].
type GAConfig struct {
	// Bounds is the search box; genomes have length Bounds.Dim().
	Bounds Bounds
	// PopSize is the number of individuals per generation.
	PopSize int
	// Generations is the number of generations to evolve.
	Generations int
	// CrossoverRate is the probability that a pair of parents is recombined.
	CrossoverRate float64
	// MutationRate is the per-gene mutation probability.
	MutationRate float64
	// MutationSigma is the Gaussian mutation standard deviation as a fraction of
	// each coordinate's box width.
	MutationSigma float64
	// TournamentK is the tournament size for parent selection.
	TournamentK int
	// Elitism is the number of best individuals copied unchanged into the next
	// generation.
	Elitism int
	// SBXEta, when positive, selects simulated binary crossover with that
	// distribution index instead of arithmetic crossover.
	SBXEta float64
	// RecordHistory enables per-generation best-fitness recording.
	RecordHistory bool
}

// DefaultGAConfig returns a reasonable configuration for the given box.
func DefaultGAConfig(b Bounds) GAConfig {
	return GAConfig{
		Bounds:        b,
		PopSize:       50,
		Generations:   200,
		CrossoverRate: 0.9,
		MutationRate:  0.1,
		MutationSigma: 0.1,
		TournamentK:   3,
		Elitism:       2,
		SBXEta:        0,
	}
}

func (c GAConfig) validate() error {
	if !c.Bounds.Valid() {
		return ErrEmptyBounds
	}
	if c.PopSize < 2 || c.Generations <= 0 {
		return ErrInvalidConfig
	}
	if c.Elitism < 0 || c.Elitism >= c.PopSize {
		return ErrInvalidConfig
	}
	return nil
}

// InitPopulation creates a population of n individuals with genomes drawn
// uniformly from the box and their fitnesses evaluated.
func InitPopulation(f ObjectiveFunc, b Bounds, n int, rng *RNG) Population {
	pop := make(Population, n)
	for i := range pop {
		g := rng.UniformVec(b)
		pop[i] = Individual{Genes: g, Fitness: f(g)}
	}
	return pop
}

// RunGA evolves a real-coded genetic algorithm and returns the best individual
// found. It is deterministic given rng.
func RunGA(f ObjectiveFunc, cfg GAConfig, rng *RNG) (Result, error) {
	if err := cfg.validate(); err != nil {
		return Result{}, err
	}
	pop := InitPopulation(f, cfg.Bounds, cfg.PopSize, rng)
	evals := cfg.PopSize
	best := pop.Best()
	res := Result{}
	gen := 0
	for ; gen < cfg.Generations; gen++ {
		pop.SortByFitness()
		if pop[0].Fitness < best.Fitness {
			best = pop[0].Clone()
		}
		next := make(Population, 0, cfg.PopSize)
		for e := 0; e < cfg.Elitism; e++ {
			next = append(next, pop[e].Clone())
		}
		for len(next) < cfg.PopSize {
			pi := TournamentSelect(pop, cfg.TournamentK, rng)
			qi := TournamentSelect(pop, cfg.TournamentK, rng)
			var c1, c2 []float64
			if rng.Float64() < cfg.CrossoverRate {
				if cfg.SBXEta > 0 {
					c1, c2 = SBXCrossover(pop[pi].Genes, pop[qi].Genes, cfg.SBXEta, rng)
				} else {
					c1, c2 = ArithmeticCrossover(pop[pi].Genes, pop[qi].Genes, rng)
				}
			} else {
				c1 = VecCopy(pop[pi].Genes)
				c2 = VecCopy(pop[qi].Genes)
			}
			c1 = GaussianMutation(c1, cfg.MutationRate, cfg.MutationSigma, cfg.Bounds, rng)
			c2 = GaussianMutation(c2, cfg.MutationRate, cfg.MutationSigma, cfg.Bounds, rng)
			cfg.Bounds.ClipInPlace(c1)
			cfg.Bounds.ClipInPlace(c2)
			next = append(next, Individual{Genes: c1, Fitness: f(c1)})
			evals++
			if len(next) < cfg.PopSize {
				next = append(next, Individual{Genes: c2, Fitness: f(c2)})
				evals++
			}
		}
		pop = next
		if cfg.RecordHistory {
			res.History = append(res.History, best.Fitness)
		}
	}
	pop.SortByFitness()
	if pop[0].Fitness < best.Fitness {
		best = pop[0].Clone()
	}
	res.X = best.Genes
	res.F = best.Fitness
	res.Iterations = gen
	res.Evaluations = evals
	return res, nil
}
