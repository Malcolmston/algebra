package metaheuristics

import "math"

// Particle is one member of a particle swarm: its current position and
// velocity, together with the best position it has personally visited.
type Particle struct {
	Position    []float64
	Velocity    []float64
	Best        []float64
	BestFitness float64
}

// PSOConfig configures [RunPSO].
type PSOConfig struct {
	// Bounds is the search box.
	Bounds Bounds
	// Swarm is the number of particles.
	Swarm int
	// Iterations is the number of update steps.
	Iterations int
	// Inertia is the inertia weight w applied to the previous velocity.
	Inertia float64
	// InertiaDamp multiplies the inertia weight each iteration (1 = constant).
	InertiaDamp float64
	// Cognitive is the personal-best acceleration coefficient c1.
	Cognitive float64
	// Social is the global-best acceleration coefficient c2.
	Social float64
	// VelocityClamp limits |velocity| to VelocityClamp times the box width per
	// coordinate. Zero disables clamping.
	VelocityClamp float64
	// RecordHistory enables per-iteration best-value recording.
	RecordHistory bool
}

// DefaultPSOConfig returns a reasonable configuration for the given box using
// the well-known constriction-inspired coefficients.
func DefaultPSOConfig(b Bounds) PSOConfig {
	return PSOConfig{
		Bounds:        b,
		Swarm:         40,
		Iterations:    300,
		Inertia:       0.729,
		InertiaDamp:   1.0,
		Cognitive:     1.49445,
		Social:        1.49445,
		VelocityClamp: 0.5,
	}
}

func (c PSOConfig) validate() error {
	if !c.Bounds.Valid() {
		return ErrEmptyBounds
	}
	if c.Swarm < 1 || c.Iterations <= 0 {
		return ErrInvalidConfig
	}
	return nil
}

// ConstrictionFactor returns Clerc's constriction coefficient chi for the sum
// of acceleration coefficients phi = c1 + c2 (which must exceed 4). Multiplying
// the velocity update by chi guarantees swarm convergence without velocity
// clamping.
func ConstrictionFactor(c1, c2 float64) float64 {
	phi := c1 + c2
	if phi <= 4 {
		return 1
	}
	return 2 / math.Abs(2-phi-math.Sqrt(phi*phi-4*phi))
}

// InitSwarm creates a swarm of n particles with random positions and small
// random velocities, evaluating each particle's initial fitness.
func InitSwarm(f ObjectiveFunc, b Bounds, n int, rng *RNG) []Particle {
	width := b.Width()
	sw := make([]Particle, n)
	for i := range sw {
		pos := rng.UniformVec(b)
		vel := make([]float64, b.Dim())
		for d := range vel {
			vel[d] = (rng.Float64()*2 - 1) * 0.1 * width[d]
		}
		sw[i] = Particle{
			Position:    pos,
			Velocity:    vel,
			Best:        VecCopy(pos),
			BestFitness: f(pos),
		}
	}
	return sw
}

// RunPSO minimizes f over the box using particle swarm optimization with a
// global (gbest) topology. It is deterministic given rng.
func RunPSO(f ObjectiveFunc, cfg PSOConfig, rng *RNG) (Result, error) {
	if err := cfg.validate(); err != nil {
		return Result{}, err
	}
	dim := cfg.Bounds.Dim()
	width := cfg.Bounds.Width()
	sw := InitSwarm(f, cfg.Bounds, cfg.Swarm, rng)
	evals := cfg.Swarm

	gBest := VecCopy(sw[0].Best)
	gBestF := sw[0].BestFitness
	for _, p := range sw {
		if p.BestFitness < gBestF {
			gBestF = p.BestFitness
			gBest = VecCopy(p.Best)
		}
	}

	w := cfg.Inertia
	res := Result{}
	iter := 0
	for ; iter < cfg.Iterations; iter++ {
		for i := range sw {
			p := &sw[i]
			for d := 0; d < dim; d++ {
				r1 := rng.Float64()
				r2 := rng.Float64()
				p.Velocity[d] = w*p.Velocity[d] +
					cfg.Cognitive*r1*(p.Best[d]-p.Position[d]) +
					cfg.Social*r2*(gBest[d]-p.Position[d])
				if cfg.VelocityClamp > 0 {
					vmax := cfg.VelocityClamp * width[d]
					p.Velocity[d] = Clamp(p.Velocity[d], -vmax, vmax)
				}
				p.Position[d] += p.Velocity[d]
			}
			cfg.Bounds.ClipInPlace(p.Position)
			fit := f(p.Position)
			evals++
			if fit < p.BestFitness {
				p.BestFitness = fit
				copy(p.Best, p.Position)
				if fit < gBestF {
					gBestF = fit
					copy(gBest, p.Position)
				}
			}
		}
		if cfg.InertiaDamp != 0 {
			w *= cfg.InertiaDamp
		}
		if cfg.RecordHistory {
			res.History = append(res.History, gBestF)
		}
	}
	res.X = gBest
	res.F = gBestF
	res.Iterations = iter
	res.Evaluations = evals
	return res, nil
}
