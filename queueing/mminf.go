package queueing

import "math"

// MMInf models an infinite-server queue with Poisson arrivals of rate Lambda
// and exponential service of rate Mu (Kendall notation M/M/∞). Every customer
// enters service immediately, so there is never any waiting. In steady state
// the number in system is Poisson with mean a = Lambda/Mu.
type MMInf struct {
	Lambda float64 // arrival rate
	Mu     float64 // per-customer service rate
}

// NewMMInf constructs an [MMInf] queue, returning an error for non-positive
// rates.
func NewMMInf(lambda, mu float64) (MMInf, error) {
	if lambda <= 0 || mu <= 0 {
		return MMInf{}, ErrNonPositiveRate
	}
	return MMInf{Lambda: lambda, Mu: mu}, nil
}

// OfferedLoad returns the offered load a = Lambda/Mu, which equals the mean
// number in system.
func (q MMInf) OfferedLoad() float64 { return OfferedLoad(q.Lambda, q.Mu) }

// Pn returns the steady-state probability P(N=n) = e^{-a} a^n/n!, the Poisson
// mass with mean a = Lambda/Mu. It returns 0 for negative n.
func (q MMInf) Pn(n int) float64 {
	a := q.OfferedLoad()
	if math.IsNaN(a) {
		return math.NaN()
	}
	return PoissonPMF(n, a)
}

// P0 returns the probability that the system is empty, e^{-a}.
func (q MMInf) P0() float64 {
	a := q.OfferedLoad()
	if math.IsNaN(a) {
		return math.NaN()
	}
	return math.Exp(-a)
}

// L returns the mean number of customers in the system, a = Lambda/Mu.
func (q MMInf) L() float64 { return q.OfferedLoad() }

// Lq returns the mean number waiting in the queue, which is always 0.
func (q MMInf) Lq() float64 { return 0 }

// W returns the mean sojourn time, equal to a single service time 1/Mu.
func (q MMInf) W() float64 {
	if q.Mu <= 0 {
		return math.NaN()
	}
	return 1 / q.Mu
}

// Wq returns the mean waiting time in the queue, which is always 0.
func (q MMInf) Wq() float64 { return 0 }

// VarianceN returns the variance of the number in system, equal to the mean a
// for the Poisson distribution.
func (q MMInf) VarianceN() float64 { return q.OfferedLoad() }

// MeanBusyServers returns the mean number of busy servers, equal to a.
func (q MMInf) MeanBusyServers() float64 { return q.OfferedLoad() }
