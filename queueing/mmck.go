package queueing

import "math"

// MMcK models a multi-server queue with Poisson arrivals of rate Lambda, C
// exponential servers of rate Mu each and a finite system capacity K>=C.
// Arrivals that find K customers present are lost (Kendall notation M/M/c/K).
// The queue is always stable because the state space is finite.
type MMcK struct {
	Lambda float64 // arrival rate
	Mu     float64 // per-server service rate
	C      int     // number of servers
	K      int     // system capacity (>=C)
}

// NewMMcK constructs an [MMcK] queue, returning an error for non-positive
// rates, a non-positive server count or a capacity smaller than the number of
// servers.
func NewMMcK(lambda, mu float64, c, k int) (MMcK, error) {
	if lambda <= 0 || mu <= 0 {
		return MMcK{}, ErrNonPositiveRate
	}
	if c <= 0 {
		return MMcK{}, ErrServers
	}
	if k < c {
		return MMcK{}, ErrCapacity
	}
	return MMcK{Lambda: lambda, Mu: mu, C: c, K: k}, nil
}

// dist returns the (normalized) stationary distribution over states 0..K.
func (q MMcK) dist() []float64 {
	p := make([]float64, q.K+1)
	p[0] = 1
	total := 1.0
	for n := 1; n <= q.K; n++ {
		servers := n
		if servers > q.C {
			servers = q.C
		}
		p[n] = p[n-1] * q.Lambda / (float64(servers) * q.Mu)
		total += p[n]
	}
	for n := range p {
		p[n] /= total
	}
	return p
}

// OfferedLoad returns the offered load a = Lambda/Mu in erlangs.
func (q MMcK) OfferedLoad() float64 { return OfferedLoad(q.Lambda, q.Mu) }

// Rho returns the nominal per-server traffic intensity Lambda/(C*Mu). For a
// finite-capacity queue it may exceed one.
func (q MMcK) Rho() float64 { return Utilization(q.Lambda, q.Mu, q.C) }

// P0 returns the steady-state probability that the system is empty.
func (q MMcK) P0() float64 { return q.dist()[0] }

// Pn returns the steady-state probability P(N=n) for 0<=n<=K, and 0 otherwise.
func (q MMcK) Pn(n int) float64 {
	if n < 0 || n > q.K {
		return 0
	}
	return q.dist()[n]
}

// BlockingProb returns the loss probability P(N=K).
func (q MMcK) BlockingProb() float64 { return q.dist()[q.K] }

// EffectiveArrivalRate returns the accepted-customer rate
// Lambda*(1-BlockingProb).
func (q MMcK) EffectiveArrivalRate() float64 {
	return q.Lambda * (1 - q.BlockingProb())
}

// Throughput returns the departure rate, equal to the effective arrival rate.
func (q MMcK) Throughput() float64 { return q.EffectiveArrivalRate() }

// L returns the mean number of customers in the system.
func (q MMcK) L() float64 {
	p := q.dist()
	sum := 0.0
	for n, pn := range p {
		sum += float64(n) * pn
	}
	return sum
}

// Lq returns the mean number of customers waiting in the queue.
func (q MMcK) Lq() float64 {
	p := q.dist()
	sum := 0.0
	for n := q.C + 1; n <= q.K; n++ {
		sum += float64(n-q.C) * p[n]
	}
	return sum
}

// MeanBusyServers returns the mean number of busy servers, which equals
// EffectiveArrivalRate/Mu.
func (q MMcK) MeanBusyServers() float64 {
	return q.EffectiveArrivalRate() / q.Mu
}

// Utilization returns the mean fraction of servers busy,
// MeanBusyServers / C.
func (q MMcK) Utilization() float64 {
	return q.MeanBusyServers() / float64(q.C)
}

// W returns the mean sojourn time L/lambdaEff.
func (q MMcK) W() float64 {
	le := q.EffectiveArrivalRate()
	if le <= 0 {
		return math.NaN()
	}
	return q.L() / le
}

// Wq returns the mean waiting time in the queue Lq/lambdaEff.
func (q MMcK) Wq() float64 {
	le := q.EffectiveArrivalRate()
	if le <= 0 {
		return math.NaN()
	}
	return q.Lq() / le
}

// ProbWait returns the probability that an accepted arrival must wait, i.e. it
// finds at least C but fewer than K customers present, normalized by the
// acceptance probability.
func (q MMcK) ProbWait() float64 {
	p := q.dist()
	waiting := 0.0
	for n := q.C; n < q.K; n++ {
		waiting += p[n]
	}
	accept := 1 - p[q.K]
	if accept <= 0 {
		return math.NaN()
	}
	return waiting / accept
}
