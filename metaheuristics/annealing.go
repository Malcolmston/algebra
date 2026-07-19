package metaheuristics

import "math"

// CoolingSchedule maps an iteration index k (0-based) to a temperature, given
// the initial temperature T0 and the total number of iterations kMax.
type CoolingSchedule func(T0 float64, k, kMax int) float64

// ExponentialCooling returns a schedule T_k = T0 * alpha^k with 0 < alpha < 1.
// It is the most common geometric cooling law.
func ExponentialCooling(alpha float64) CoolingSchedule {
	return func(T0 float64, k, _ int) float64 {
		return T0 * math.Pow(alpha, float64(k))
	}
}

// GeometricCooling is a synonym for [ExponentialCooling]; each step multiplies
// the temperature by alpha.
func GeometricCooling(alpha float64) CoolingSchedule {
	return ExponentialCooling(alpha)
}

// LinearCooling returns a schedule that decreases the temperature linearly from
// T0 at k=0 to (approximately) 0 at k=kMax.
func LinearCooling() CoolingSchedule {
	return func(T0 float64, k, kMax int) float64 {
		if kMax <= 0 {
			return T0
		}
		t := T0 * (1 - float64(k)/float64(kMax))
		if t < 0 {
			return 0
		}
		return t
	}
}

// LogarithmicCooling returns the classical Boltzmann schedule
// T_k = T0 / ln(k + e), which cools slowly enough to guarantee (in the limit)
// convergence to a global optimum.
func LogarithmicCooling() CoolingSchedule {
	return func(T0 float64, k, _ int) float64 {
		return T0 / math.Log(float64(k)+math.E)
	}
}

// BoltzmannCooling returns the schedule T_k = T0 / (1 + ln(1+k)).
func BoltzmannCooling() CoolingSchedule {
	return func(T0 float64, k, _ int) float64 {
		return T0 / (1 + math.Log(1+float64(k)))
	}
}

// CauchyCooling returns the fast Cauchy schedule T_k = T0 / (1 + k).
func CauchyCooling() CoolingSchedule {
	return func(T0 float64, k, _ int) float64 {
		return T0 / (1 + float64(k))
	}
}

// QuadraticCooling returns a schedule that decreases quadratically from T0 to 0
// over kMax iterations.
func QuadraticCooling() CoolingSchedule {
	return func(T0 float64, k, kMax int) float64 {
		if kMax <= 0 {
			return T0
		}
		r := 1 - float64(k)/float64(kMax)
		if r < 0 {
			r = 0
		}
		return T0 * r * r
	}
}

// AcceptanceProbability returns the Metropolis acceptance probability of moving
// from an energy of curF to a candidate energy candF at temperature T. Downhill
// moves (candF <= curF) are always accepted (probability 1); uphill moves are
// accepted with probability exp(-(candF-curF)/T). A non-positive temperature
// makes uphill moves impossible.
func AcceptanceProbability(curF, candF, T float64) float64 {
	if candF <= curF {
		return 1
	}
	if T <= 0 {
		return 0
	}
	return math.Exp(-(candF - curF) / T)
}

// AnnealConfig configures [SimulatedAnnealing].
type AnnealConfig struct {
	// Bounds is the continuous search box.
	Bounds Bounds
	// T0 is the initial temperature.
	T0 float64
	// TMin is the temperature below which the run stops early. Zero disables it.
	TMin float64
	// MaxIterations is the maximum number of proposed moves.
	MaxIterations int
	// StepSize is the relative standard deviation of the Gaussian proposal, as
	// a fraction of each coordinate's box width.
	StepSize float64
	// Schedule is the cooling schedule; if nil, ExponentialCooling(0.95) is used.
	Schedule CoolingSchedule
	// Reanneal, when true, scales the proposal step with the current
	// temperature, so the search localizes as it cools.
	Reanneal bool
	// RecordHistory enables per-iteration best-value recording.
	RecordHistory bool
}

// DefaultAnnealConfig returns a reasonable configuration for the given box.
func DefaultAnnealConfig(b Bounds) AnnealConfig {
	return AnnealConfig{
		Bounds:        b,
		T0:            10,
		TMin:          1e-8,
		MaxIterations: 5000,
		StepSize:      0.15,
		Schedule:      ExponentialCooling(0.995),
		Reanneal:      true,
	}
}

func (c AnnealConfig) validate() error {
	if !c.Bounds.Valid() {
		return ErrEmptyBounds
	}
	if c.MaxIterations <= 0 || c.StepSize <= 0 || c.T0 <= 0 {
		return ErrInvalidConfig
	}
	return nil
}

// SimulatedAnnealing minimizes f over the continuous box using the Metropolis
// criterion and the configured cooling schedule, starting from the given point
// (or the box center if start is nil). It is deterministic given rng.
func SimulatedAnnealing(f ObjectiveFunc, cfg AnnealConfig, start []float64, rng *RNG) (Result, error) {
	if err := cfg.validate(); err != nil {
		return Result{}, err
	}
	sched := cfg.Schedule
	if sched == nil {
		sched = ExponentialCooling(0.95)
	}
	dim := cfg.Bounds.Dim()
	if start == nil {
		start = cfg.Bounds.Center()
	}
	if len(start) != dim {
		return Result{}, ErrDimMismatch
	}
	width := cfg.Bounds.Width()
	cur := cfg.Bounds.Clip(start)
	curF := f(cur)
	best := VecCopy(cur)
	bestF := curF
	evals := 1
	res := Result{}
	iter := 0
	for ; iter < cfg.MaxIterations; iter++ {
		T := sched(cfg.T0, iter, cfg.MaxIterations)
		if cfg.TMin > 0 && T < cfg.TMin {
			iter++
			break
		}
		scale := cfg.StepSize
		if cfg.Reanneal {
			scale *= math.Sqrt(math.Max(T/cfg.T0, 1e-6))
		}
		cand := make([]float64, dim)
		for i := range cand {
			cand[i] = cur[i] + rng.NormFloat64()*scale*width[i]
		}
		cfg.Bounds.ClipInPlace(cand)
		candF := f(cand)
		evals++
		if AcceptanceProbability(curF, candF, T) >= rng.Float64() {
			cur, curF = cand, candF
			if curF < bestF {
				best = VecCopy(cur)
				bestF = curF
			}
		}
		if cfg.RecordHistory {
			res.History = append(res.History, bestF)
		}
	}
	res.X = best
	res.F = bestF
	res.Iterations = iter
	res.Evaluations = evals
	res.Converged = cfg.TMin > 0 && iter < cfg.MaxIterations
	return res, nil
}

// AnnealDiscrete minimizes an objective over an abstract discrete state space
// using simulated annealing. The caller supplies an initial state, a neighbour
// function that proposes a random successor state, and an energy function. It
// returns the best state found and its energy. States are represented as int
// slices (for example a permutation or a bit vector).
func AnnealDiscrete(
	energy func(state []int) float64,
	neighbor func(state []int, rng *RNG) []int,
	initial []int,
	T0, alpha float64,
	iterations int,
	rng *RNG,
) (best []int, bestEnergy float64, err error) {
	if iterations <= 0 || T0 <= 0 || alpha <= 0 || alpha >= 1 {
		return nil, 0, ErrInvalidConfig
	}
	cur := append([]int(nil), initial...)
	curE := energy(cur)
	best = append([]int(nil), cur...)
	bestEnergy = curE
	T := T0
	for k := 0; k < iterations; k++ {
		cand := neighbor(cur, rng)
		candE := energy(cand)
		if AcceptanceProbability(curE, candE, T) >= rng.Float64() {
			cur = cand
			curE = candE
			if curE < bestEnergy {
				best = append([]int(nil), cur...)
				bestEnergy = curE
			}
		}
		T *= alpha
	}
	return best, bestEnergy, nil
}
