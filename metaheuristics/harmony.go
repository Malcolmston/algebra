package metaheuristics

// HarmonyConfig configures [RunHarmonySearch].
type HarmonyConfig struct {
	// Bounds is the search box.
	Bounds Bounds
	// MemorySize is the number of solution vectors in harmony memory (HMS).
	MemorySize int
	// Iterations is the number of new harmonies improvised.
	Iterations int
	// MemoryRate is the harmony memory considering rate (HMCR), the probability
	// of drawing a component from memory rather than randomly.
	MemoryRate float64
	// PitchRate is the pitch adjusting rate (PAR), the probability of tweaking a
	// memory-drawn component.
	PitchRate float64
	// Bandwidth is the pitch-adjustment amplitude as a fraction of the box
	// width.
	Bandwidth float64
	// RecordHistory enables per-iteration best-value recording.
	RecordHistory bool
}

// DefaultHarmonyConfig returns a reasonable configuration for the given box.
func DefaultHarmonyConfig(b Bounds) HarmonyConfig {
	return HarmonyConfig{
		Bounds:     b,
		MemorySize: 30,
		Iterations: 2000,
		MemoryRate: 0.95,
		PitchRate:  0.3,
		Bandwidth:  0.02,
	}
}

func (c HarmonyConfig) validate() error {
	if !c.Bounds.Valid() {
		return ErrEmptyBounds
	}
	if c.MemorySize < 1 || c.Iterations <= 0 {
		return ErrInvalidConfig
	}
	if c.MemoryRate < 0 || c.MemoryRate > 1 || c.PitchRate < 0 || c.PitchRate > 1 {
		return ErrInvalidConfig
	}
	return nil
}

// RunHarmonySearch minimizes f over the box using the harmony search
// metaheuristic. It is deterministic given rng.
func RunHarmonySearch(f ObjectiveFunc, cfg HarmonyConfig, rng *RNG) (Result, error) {
	if err := cfg.validate(); err != nil {
		return Result{}, err
	}
	dim := cfg.Bounds.Dim()
	width := cfg.Bounds.Width()

	mem := make([][]float64, cfg.MemorySize)
	memF := make([]float64, cfg.MemorySize)
	worst := 0
	for i := range mem {
		mem[i] = rng.UniformVec(cfg.Bounds)
		memF[i] = f(mem[i])
		if memF[i] > memF[worst] {
			worst = i
		}
	}
	evals := cfg.MemorySize
	best := 0
	for i := range memF {
		if memF[i] < memF[best] {
			best = i
		}
	}
	res := Result{}
	for it := 0; it < cfg.Iterations; it++ {
		newH := make([]float64, dim)
		for d := 0; d < dim; d++ {
			if rng.Float64() < cfg.MemoryRate {
				newH[d] = mem[rng.Intn(cfg.MemorySize)][d]
				if rng.Float64() < cfg.PitchRate {
					newH[d] += (rng.Float64()*2 - 1) * cfg.Bandwidth * width[d]
				}
			} else {
				newH[d] = rng.Float64Range(cfg.Bounds.Lower[d], cfg.Bounds.Upper[d])
			}
		}
		cfg.Bounds.ClipInPlace(newH)
		nf := f(newH)
		evals++
		if nf < memF[worst] {
			mem[worst] = newH
			memF[worst] = nf
			// recompute worst and best
			worst = 0
			for i := range memF {
				if memF[i] > memF[worst] {
					worst = i
				}
			}
			best = 0
			for i := range memF {
				if memF[i] < memF[best] {
					best = i
				}
			}
		}
		if cfg.RecordHistory {
			res.History = append(res.History, memF[best])
		}
	}
	res.X = VecCopy(mem[best])
	res.F = memF[best]
	res.Iterations = cfg.Iterations
	res.Evaluations = evals
	return res, nil
}
