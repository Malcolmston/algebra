package metaheuristics

// DEStrategy names a differential-evolution mutation strategy.
type DEStrategy int

const (
	// DERand1 is the DE/rand/1 strategy: base is a random vector.
	DERand1 DEStrategy = iota
	// DEBest1 is the DE/best/1 strategy: base is the current best vector.
	DEBest1
	// DECurrentToBest1 is the DE/current-to-best/1 strategy.
	DECurrentToBest1
	// DERand2 is the DE/rand/2 strategy with two difference vectors.
	DERand2
	// DEBest2 is the DE/best/2 strategy with two difference vectors.
	DEBest2
)

// String returns the conventional name of the strategy, such as "DE/rand/1".
func (s DEStrategy) String() string {
	switch s {
	case DERand1:
		return "DE/rand/1"
	case DEBest1:
		return "DE/best/1"
	case DECurrentToBest1:
		return "DE/current-to-best/1"
	case DERand2:
		return "DE/rand/2"
	case DEBest2:
		return "DE/best/2"
	default:
		return "DE/unknown"
	}
}

// DEConfig configures [RunDE].
type DEConfig struct {
	// Bounds is the search box.
	Bounds Bounds
	// PopSize is the number of population vectors.
	PopSize int
	// Generations is the number of generations.
	Generations int
	// F is the differential weight (scaling factor), typically in [0.4, 1].
	F float64
	// CR is the crossover probability, in [0, 1].
	CR float64
	// Strategy selects the mutation strategy.
	Strategy DEStrategy
	// RecordHistory enables per-generation best-value recording.
	RecordHistory bool
}

// DefaultDEConfig returns a reasonable configuration for the given box.
func DefaultDEConfig(b Bounds) DEConfig {
	return DEConfig{
		Bounds:      b,
		PopSize:     15 * b.Dim(),
		Generations: 300,
		F:           0.8,
		CR:          0.9,
		Strategy:    DERand1,
	}
}

func (c DEConfig) validate() error {
	if !c.Bounds.Valid() {
		return ErrEmptyBounds
	}
	if c.PopSize < 4 || c.Generations <= 0 {
		return ErrInvalidConfig
	}
	if c.CR < 0 || c.CR > 1 {
		return ErrInvalidConfig
	}
	return nil
}

// deMutant builds a donor vector for index i under the chosen strategy.
func deMutant(pop [][]float64, fit []float64, i, bestIdx int, cfg DEConfig, rng *RNG) []float64 {
	n := len(pop)
	dim := len(pop[i])
	pick := func(exclude ...int) int {
		for {
			r := rng.Intn(n)
			ok := true
			for _, e := range exclude {
				if r == e {
					ok = false
					break
				}
			}
			if ok {
				return r
			}
		}
	}
	donor := make([]float64, dim)
	switch cfg.Strategy {
	case DEBest1:
		r1 := pick(i, bestIdx)
		r2 := pick(i, bestIdx, r1)
		for d := 0; d < dim; d++ {
			donor[d] = pop[bestIdx][d] + cfg.F*(pop[r1][d]-pop[r2][d])
		}
	case DECurrentToBest1:
		r1 := pick(i)
		r2 := pick(i, r1)
		for d := 0; d < dim; d++ {
			donor[d] = pop[i][d] + cfg.F*(pop[bestIdx][d]-pop[i][d]) + cfg.F*(pop[r1][d]-pop[r2][d])
		}
	case DERand2:
		r1 := pick(i)
		r2 := pick(i, r1)
		r3 := pick(i, r1, r2)
		r4 := pick(i, r1, r2, r3)
		r5 := pick(i, r1, r2, r3, r4)
		for d := 0; d < dim; d++ {
			donor[d] = pop[r1][d] + cfg.F*(pop[r2][d]-pop[r3][d]) + cfg.F*(pop[r4][d]-pop[r5][d])
		}
	case DEBest2:
		r1 := pick(i, bestIdx)
		r2 := pick(i, bestIdx, r1)
		r3 := pick(i, bestIdx, r1, r2)
		r4 := pick(i, bestIdx, r1, r2, r3)
		for d := 0; d < dim; d++ {
			donor[d] = pop[bestIdx][d] + cfg.F*(pop[r1][d]-pop[r2][d]) + cfg.F*(pop[r3][d]-pop[r4][d])
		}
	default: // DERand1
		r1 := pick(i)
		r2 := pick(i, r1)
		r3 := pick(i, r1, r2)
		for d := 0; d < dim; d++ {
			donor[d] = pop[r1][d] + cfg.F*(pop[r2][d]-pop[r3][d])
		}
	}
	return donor
}

// RunDE minimizes f over the box using differential evolution with binomial
// crossover. It is deterministic given rng.
func RunDE(f ObjectiveFunc, cfg DEConfig, rng *RNG) (Result, error) {
	if err := cfg.validate(); err != nil {
		return Result{}, err
	}
	dim := cfg.Bounds.Dim()
	pop := make([][]float64, cfg.PopSize)
	fit := make([]float64, cfg.PopSize)
	bestIdx := 0
	for i := range pop {
		pop[i] = rng.UniformVec(cfg.Bounds)
		fit[i] = f(pop[i])
		if fit[i] < fit[bestIdx] {
			bestIdx = i
		}
	}
	evals := cfg.PopSize
	res := Result{}
	gen := 0
	for ; gen < cfg.Generations; gen++ {
		for i := 0; i < cfg.PopSize; i++ {
			donor := deMutant(pop, fit, i, bestIdx, cfg, rng)
			trial := make([]float64, dim)
			jrand := rng.Intn(dim)
			for d := 0; d < dim; d++ {
				if rng.Float64() < cfg.CR || d == jrand {
					trial[d] = donor[d]
				} else {
					trial[d] = pop[i][d]
				}
			}
			cfg.Bounds.ClipInPlace(trial)
			tf := f(trial)
			evals++
			if tf <= fit[i] {
				pop[i] = trial
				fit[i] = tf
				if tf < fit[bestIdx] {
					bestIdx = i
				}
			}
		}
		if cfg.RecordHistory {
			res.History = append(res.History, fit[bestIdx])
		}
	}
	res.X = VecCopy(pop[bestIdx])
	res.F = fit[bestIdx]
	res.Iterations = gen
	res.Evaluations = evals
	return res, nil
}
