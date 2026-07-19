package metaheuristics

import "math"

// HillClimbConfig configures the hill-climbing optimizers. A candidate step is
// drawn from an isotropic Gaussian of standard deviation StepSize (scaled by
// each coordinate's box width); a step is accepted only if it strictly improves
// the objective (or does not worsen it when acceptance of plateaus is desired).
type HillClimbConfig struct {
	// Bounds is the search box.
	Bounds Bounds
	// MaxIterations is the number of candidate steps per climb.
	MaxIterations int
	// StepSize is the relative standard deviation of the Gaussian step, as a
	// fraction of each coordinate's box width.
	StepSize float64
	// Adaptive, when true, shrinks the step size after consecutive failures
	// and grows it after successes (a simple (1+1) adaptation).
	Adaptive bool
	// Tolerance is the objective improvement below which a step counts as no
	// progress for the stall counter; when StallLimit consecutive steps make no
	// progress the climb stops early.
	Tolerance float64
	// StallLimit is the number of consecutive non-improving steps that ends a
	// climb early. Zero disables the check.
	StallLimit int
	// RecordHistory enables per-iteration best-value recording.
	RecordHistory bool
}

// DefaultHillClimbConfig returns a reasonable configuration for the given
// search box.
func DefaultHillClimbConfig(b Bounds) HillClimbConfig {
	return HillClimbConfig{
		Bounds:        b,
		MaxIterations: 1000,
		StepSize:      0.1,
		Adaptive:      true,
		Tolerance:     1e-12,
		StallLimit:    0,
	}
}

func (c HillClimbConfig) validate() error {
	if !c.Bounds.Valid() {
		return ErrEmptyBounds
	}
	if c.MaxIterations <= 0 || c.StepSize <= 0 {
		return ErrInvalidConfig
	}
	return nil
}

// HillClimb runs stochastic hill climbing from the supplied start point,
// returning the best point found. It uses the given [RNG] for its steps.
func HillClimb(f ObjectiveFunc, cfg HillClimbConfig, start []float64, rng *RNG) (Result, error) {
	if err := cfg.validate(); err != nil {
		return Result{}, err
	}
	if len(start) != cfg.Bounds.Dim() {
		return Result{}, ErrDimMismatch
	}
	width := cfg.Bounds.Width()
	cur := cfg.Bounds.Clip(start)
	curF := f(cur)
	evals := 1
	step := cfg.StepSize
	stall := 0
	res := Result{}
	if cfg.RecordHistory {
		res.History = make([]float64, 0, cfg.MaxIterations)
	}
	iter := 0
	for ; iter < cfg.MaxIterations; iter++ {
		cand := make([]float64, len(cur))
		for i := range cand {
			cand[i] = cur[i] + rng.NormFloat64()*step*width[i]
		}
		cfg.Bounds.ClipInPlace(cand)
		candF := f(cand)
		evals++
		if candF < curF {
			improve := curF - candF
			cur, curF = cand, candF
			if cfg.Adaptive {
				step *= 1.15
				if step > 1 {
					step = 1
				}
			}
			if improve < cfg.Tolerance {
				stall++
			} else {
				stall = 0
			}
		} else {
			stall++
			if cfg.Adaptive {
				step *= 0.9
				if step < 1e-12 {
					step = 1e-12
				}
			}
		}
		if cfg.RecordHistory {
			res.History = append(res.History, curF)
		}
		if cfg.StallLimit > 0 && stall >= cfg.StallLimit {
			iter++
			break
		}
	}
	res.X = cur
	res.F = curF
	res.Iterations = iter
	res.Evaluations = evals
	return res, nil
}

// HillClimbRestarts runs HillClimb from restarts independent random starting
// points drawn uniformly from the box, returning the best result across all
// climbs. This mitigates hill climbing's tendency to stop at local optima.
func HillClimbRestarts(f ObjectiveFunc, cfg HillClimbConfig, restarts int, rng *RNG) (Result, error) {
	if err := cfg.validate(); err != nil {
		return Result{}, err
	}
	if restarts <= 0 {
		return Result{}, ErrInvalidConfig
	}
	best := Result{F: math.Inf(1)}
	totalEval := 0
	totalIter := 0
	for r := 0; r < restarts; r++ {
		start := rng.UniformVec(cfg.Bounds)
		res, err := HillClimb(f, cfg, start, rng)
		if err != nil {
			return Result{}, err
		}
		totalEval += res.Evaluations
		totalIter += res.Iterations
		if res.F < best.F {
			best.X = res.X
			best.F = res.F
		}
	}
	best.Evaluations = totalEval
	best.Iterations = totalIter
	return best, nil
}

// RandomSearch samples MaxIterations points uniformly from the box and returns
// the best. It provides a baseline against which other optimizers are compared.
func RandomSearch(f ObjectiveFunc, b Bounds, iterations int, rng *RNG) (Result, error) {
	if !b.Valid() {
		return Result{}, ErrEmptyBounds
	}
	if iterations <= 0 {
		return Result{}, ErrInvalidConfig
	}
	best := rng.UniformVec(b)
	bestF := f(best)
	for i := 1; i < iterations; i++ {
		cand := rng.UniformVec(b)
		cf := f(cand)
		if cf < bestF {
			best, bestF = cand, cf
		}
	}
	return Result{X: best, F: bestF, Iterations: iterations, Evaluations: iterations}, nil
}

// CoordinateDescent performs axis-parallel line minimization: each iteration it
// probes both directions along every coordinate with the current step and moves
// to any improvement, shrinking the step when a full sweep makes no progress.
// It converges to a coordinate-wise local minimum, which for separable convex
// functions is the global one.
func CoordinateDescent(f ObjectiveFunc, cfg HillClimbConfig, start []float64) (Result, error) {
	if err := cfg.validate(); err != nil {
		return Result{}, err
	}
	if len(start) != cfg.Bounds.Dim() {
		return Result{}, ErrDimMismatch
	}
	width := cfg.Bounds.Width()
	cur := cfg.Bounds.Clip(start)
	curF := f(cur)
	evals := 1
	step := cfg.StepSize
	res := Result{}
	iter := 0
	for ; iter < cfg.MaxIterations; iter++ {
		improved := false
		for i := range cur {
			h := step * width[i]
			for _, dir := range [2]float64{1, -1} {
				cand := VecCopy(cur)
				cand[i] = Clamp(cur[i]+dir*h, cfg.Bounds.Lower[i], cfg.Bounds.Upper[i])
				cf := f(cand)
				evals++
				if cf < curF {
					cur, curF = cand, cf
					improved = true
					break
				}
			}
		}
		if cfg.RecordHistory {
			res.History = append(res.History, curF)
		}
		if !improved {
			step *= 0.5
			if step < 1e-12 {
				iter++
				break
			}
		}
	}
	res.X = cur
	res.F = curF
	res.Iterations = iter
	res.Evaluations = evals
	return res, nil
}
