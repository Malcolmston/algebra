package queueing

import "math"

// MM1K models a single-server queue with Poisson arrivals of rate Lambda,
// exponential service of rate Mu and a finite system capacity K (including the
// customer in service). Arrivals that find K customers present are lost
// (Kendall notation M/M/1/K). Because the state space is finite the queue is
// always stable, even when Lambda >= Mu.
type MM1K struct {
	Lambda float64 // arrival rate
	Mu     float64 // service rate
	K      int     // system capacity (>=1)
}

// NewMM1K constructs an [MM1K] queue, returning an error for non-positive rates
// or a capacity below 1.
func NewMM1K(lambda, mu float64, k int) (MM1K, error) {
	if lambda <= 0 || mu <= 0 {
		return MM1K{}, ErrNonPositiveRate
	}
	if k < 1 {
		return MM1K{}, ErrCapacity
	}
	return MM1K{Lambda: lambda, Mu: mu, K: k}, nil
}

// Rho returns the traffic intensity rho = Lambda/Mu, which may exceed one for a
// finite-capacity queue.
func (q MM1K) Rho() float64 {
	if q.Mu <= 0 {
		return math.NaN()
	}
	return q.Lambda / q.Mu
}

// P0 returns the steady-state probability that the system is empty.
func (q MM1K) P0() float64 {
	rho := q.Rho()
	if math.IsNaN(rho) {
		return math.NaN()
	}
	if rho == 1 {
		return 1 / float64(q.K+1)
	}
	return (1 - rho) / (1 - math.Pow(rho, float64(q.K+1)))
}

// Pn returns the steady-state probability P(N=n) for 0<=n<=K, and 0 outside
// that range.
func (q MM1K) Pn(n int) float64 {
	if n < 0 || n > q.K {
		return 0
	}
	rho := q.Rho()
	if math.IsNaN(rho) {
		return math.NaN()
	}
	if rho == 1 {
		return 1 / float64(q.K+1)
	}
	return q.P0() * math.Pow(rho, float64(n))
}

// BlockingProb returns the probability P(N=K) that an arriving customer is
// blocked and lost. By PASTA this equals the loss probability.
func (q MM1K) BlockingProb() float64 { return q.Pn(q.K) }

// EffectiveArrivalRate returns the rate of accepted customers,
// Lambda*(1-BlockingProb).
func (q MM1K) EffectiveArrivalRate() float64 {
	b := q.BlockingProb()
	if math.IsNaN(b) {
		return math.NaN()
	}
	return q.Lambda * (1 - b)
}

// Throughput returns the departure rate, equal to the effective arrival rate.
func (q MM1K) Throughput() float64 { return q.EffectiveArrivalRate() }

// Utilization returns the fraction of time the server is busy, 1-P0, which also
// equals EffectiveArrivalRate/Mu.
func (q MM1K) Utilization() float64 {
	p0 := q.P0()
	if math.IsNaN(p0) {
		return math.NaN()
	}
	return 1 - p0
}

// L returns the mean number of customers in the system.
func (q MM1K) L() float64 {
	rho := q.Rho()
	if math.IsNaN(rho) {
		return math.NaN()
	}
	if rho == 1 {
		return float64(q.K) / 2
	}
	k := float64(q.K)
	num := rho * (1 - (k+1)*math.Pow(rho, k) + k*math.Pow(rho, k+1))
	den := (1 - rho) * (1 - math.Pow(rho, k+1))
	return num / den
}

// Lq returns the mean number of customers waiting in the queue, L-(1-P0).
func (q MM1K) Lq() float64 {
	l := q.L()
	p0 := q.P0()
	if math.IsNaN(l) || math.IsNaN(p0) {
		return math.NaN()
	}
	return l - (1 - p0)
}

// W returns the mean sojourn time using the effective arrival rate, L/lambdaEff.
func (q MM1K) W() float64 {
	l := q.L()
	le := q.EffectiveArrivalRate()
	if math.IsNaN(l) || le <= 0 {
		return math.NaN()
	}
	return l / le
}

// Wq returns the mean waiting time in the queue, Lq/lambdaEff.
func (q MM1K) Wq() float64 {
	lq := q.Lq()
	le := q.EffectiveArrivalRate()
	if math.IsNaN(lq) || le <= 0 {
		return math.NaN()
	}
	return lq / le
}
