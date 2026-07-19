package stochastic

import (
	"math"
	"math/rand"
)

// Reaction describes a single reaction channel of a stochastic reaction
// network: a stoichiometry vector giving the change applied to the species
// count vector when the reaction fires, and a propensity function giving the
// instantaneous firing rate as a function of the current state.
type Reaction struct {
	Stoichiometry []int
	Propensity    func(state []int) float64
}

// ReactionNetwork is a collection of reaction channels acting on NumSpecies
// molecular species.
type ReactionNetwork struct {
	NumSpecies int
	Reactions  []Reaction
}

// GillespieResult holds the output of a stochastic simulation: the event times
// and the species-count vector after each event (States[0] is the initial
// state at Times[0]).
type GillespieResult struct {
	Times  []float64
	States [][]int
}

// binom returns the number of ways to choose k items from n, as a float, used
// for stochastic mass-action propensities.
func binom(n, k int) float64 {
	if k < 0 || n < 0 || k > n {
		return 0
	}
	r := 1.0
	for i := 0; i < k; i++ {
		r *= float64(n - i)
		r /= float64(i + 1)
	}
	return r
}

// MassAction returns a propensity function for a mass-action reaction with rate
// constant k and the given reactant multiplicities (reactants[i] copies of
// species i are consumed). The propensity is k times the product of the
// combinatorial factors C(state[i], reactants[i]).
func MassAction(k float64, reactants []int) func(state []int) float64 {
	rc := append([]int(nil), reactants...)
	return func(state []int) float64 {
		p := k
		for i, n := range rc {
			if n == 0 {
				continue
			}
			if i >= len(state) {
				return 0
			}
			p *= binom(state[i], n)
		}
		return p
	}
}

// NewReaction builds a Reaction from an explicit stoichiometry vector and
// propensity function.
func NewReaction(stoich []int, propensity func(state []int) float64) Reaction {
	return Reaction{Stoichiometry: append([]int(nil), stoich...), Propensity: propensity}
}

// MassActionReaction builds a mass-action Reaction with the given stoichiometry,
// reactant multiplicities and rate constant.
func MassActionReaction(stoich, reactants []int, k float64) Reaction {
	return Reaction{
		Stoichiometry: append([]int(nil), stoich...),
		Propensity:    MassAction(k, reactants),
	}
}

// TotalPropensity returns the sum of the propensities of all reactions in the
// network at the given state.
func (net ReactionNetwork) TotalPropensity(state []int) float64 {
	s := 0.0
	for _, r := range net.Reactions {
		s += r.Propensity(state)
	}
	return s
}

// Propensities returns the propensity of each reaction at the given state.
func (net ReactionNetwork) Propensities(state []int) []float64 {
	out := make([]float64, len(net.Reactions))
	for i, r := range net.Reactions {
		out[i] = r.Propensity(state)
	}
	return out
}

// Gillespie runs the exact stochastic simulation algorithm (Gillespie's direct
// method) for the network starting from the given initial state until time
// tMax or until no reaction can fire. It returns the full event history.
func Gillespie(rng *rand.Rand, net ReactionNetwork, initial []int, tMax float64) GillespieResult {
	return GillespieMaxSteps(rng, net, initial, tMax, math.MaxInt)
}

// GillespieMaxSteps is like Gillespie but also stops after maxSteps reaction
// events, guarding against unbounded simulations.
func GillespieMaxSteps(rng *rand.Rand, net ReactionNetwork, initial []int, tMax float64, maxSteps int) GillespieResult {
	state := append([]int(nil), initial...)
	res := GillespieResult{
		Times:  []float64{0},
		States: [][]int{append([]int(nil), state...)},
	}
	t := 0.0
	props := make([]float64, len(net.Reactions))
	for step := 0; step < maxSteps; step++ {
		a0 := 0.0
		for i, r := range net.Reactions {
			props[i] = r.Propensity(state)
			a0 += props[i]
		}
		if a0 <= 0 {
			break
		}
		t += ExponentialSample(rng, a0)
		if t > tMax {
			break
		}
		j := DiscreteSample(rng, props)
		if j < 0 {
			break
		}
		for k, dv := range net.Reactions[j].Stoichiometry {
			if k < len(state) {
				state[k] += dv
			}
		}
		res.Times = append(res.Times, t)
		res.States = append(res.States, append([]int(nil), state...))
	}
	return res
}

// TauLeaping runs an explicit tau-leaping approximation of the network, firing
// each reaction a Poisson(propensity*tau) number of times per leap of size tau,
// from time 0 to tMax. Species counts are clamped to be non-negative.
func TauLeaping(rng *rand.Rand, net ReactionNetwork, initial []int, tMax, tau float64) GillespieResult {
	state := append([]int(nil), initial...)
	res := GillespieResult{
		Times:  []float64{0},
		States: [][]int{append([]int(nil), state...)},
	}
	if tau <= 0 {
		return res
	}
	for t := tau; t <= tMax+1e-12; t += tau {
		for _, r := range net.Reactions {
			a := r.Propensity(state)
			if a <= 0 {
				continue
			}
			fires := PoissonSample(rng, a*tau)
			for k, dv := range r.Stoichiometry {
				if k < len(state) {
					state[k] += dv * fires
				}
			}
		}
		for k := range state {
			if state[k] < 0 {
				state[k] = 0
			}
		}
		res.Times = append(res.Times, t)
		res.States = append(res.States, append([]int(nil), state...))
	}
	return res
}

// FinalState returns the last recorded species-count vector of the result.
func (r GillespieResult) FinalState() []int {
	if len(r.States) == 0 {
		return nil
	}
	return r.States[len(r.States)-1]
}

// FinalTime returns the last recorded time of the result.
func (r GillespieResult) FinalTime() float64 {
	if len(r.Times) == 0 {
		return 0
	}
	return r.Times[len(r.Times)-1]
}

// SpeciesTrajectory returns the times and counts of species index i across the
// simulation as a Path (counts cast to float64).
func (r GillespieResult) SpeciesTrajectory(i int) Path {
	n := len(r.Times)
	times := make([]float64, n)
	vals := make([]float64, n)
	for j := 0; j < n; j++ {
		times[j] = r.Times[j]
		if i < len(r.States[j]) {
			vals[j] = float64(r.States[j][i])
		}
	}
	return Path{Times: times, Values: vals}
}

// NumEvents returns the number of reaction events recorded (excluding the
// initial state).
func (r GillespieResult) NumEvents() int {
	if len(r.Times) == 0 {
		return 0
	}
	return len(r.Times) - 1
}

// BirthDeathNetwork returns a one-species linear birth-death network with
// per-capita birth rate birth and per-capita death rate death.
func BirthDeathNetwork(birth, death float64) ReactionNetwork {
	return ReactionNetwork{
		NumSpecies: 1,
		Reactions: []Reaction{
			MassActionReaction([]int{+1}, []int{1}, birth),
			MassActionReaction([]int{-1}, []int{1}, death),
		},
	}
}

// ImmigrationDeathNetwork returns a one-species network with constant
// immigration rate immigration and per-capita death rate death.
func ImmigrationDeathNetwork(immigration, death float64) ReactionNetwork {
	return ReactionNetwork{
		NumSpecies: 1,
		Reactions: []Reaction{
			NewReaction([]int{+1}, func(state []int) float64 { return immigration }),
			MassActionReaction([]int{-1}, []int{1}, death),
		},
	}
}

// SIRNetwork returns the stochastic SIR epidemic network on species
// (S, I, R) with infection rate constant beta and recovery rate gamma.
// Infection S+I -> 2I has propensity beta*S*I; recovery I -> R has propensity
// gamma*I.
func SIRNetwork(beta, gamma float64) ReactionNetwork {
	return ReactionNetwork{
		NumSpecies: 3,
		Reactions: []Reaction{
			MassActionReaction([]int{-1, +1, 0}, []int{1, 1, 0}, beta),
			MassActionReaction([]int{0, -1, +1}, []int{0, 1, 0}, gamma),
		},
	}
}

// LotkaVolterraNetwork returns the stochastic Lotka-Volterra predator-prey
// network on species (prey, predator) with prey birth rate alpha (per prey),
// predation rate beta (per prey-predator pair) and predator death rate gamma
// (per predator).
func LotkaVolterraNetwork(alpha, beta, gamma float64) ReactionNetwork {
	return ReactionNetwork{
		NumSpecies: 2,
		Reactions: []Reaction{
			MassActionReaction([]int{+1, 0}, []int{1, 0}, alpha),
			MassActionReaction([]int{-1, +1}, []int{1, 1}, beta),
			MassActionReaction([]int{0, -1}, []int{0, 1}, gamma),
		},
	}
}
