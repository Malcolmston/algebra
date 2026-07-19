package queueing

import "math"

// MD1 models a single-server queue with Poisson arrivals of rate Lambda and
// deterministic (constant) service time 1/Mu (Kendall notation M/D/1). It is
// the M/G/1 special case with zero service-time variance, so its queue length
// is exactly half that of the corresponding M/M/1 queue.
type MD1 struct {
	Lambda float64 // arrival rate
	Mu     float64 // service rate (service time = 1/Mu)
}

// NewMD1 constructs an [MD1] queue, returning an error for non-positive rates
// or when the stability condition Lambda < Mu is violated.
func NewMD1(lambda, mu float64) (MD1, error) {
	if lambda <= 0 || mu <= 0 {
		return MD1{}, ErrNonPositiveRate
	}
	if lambda >= mu {
		return MD1{}, ErrUnstable
	}
	return MD1{Lambda: lambda, Mu: mu}, nil
}

// Rho returns the utilization rho = Lambda/Mu.
func (q MD1) Rho() float64 {
	if q.Mu <= 0 {
		return math.NaN()
	}
	return q.Lambda / q.Mu
}

// Stable reports whether Lambda < Mu.
func (q MD1) Stable() bool { return q.Lambda > 0 && q.Mu > 0 && q.Lambda < q.Mu }

// Lq returns the mean number waiting in the queue, rho^2/(2(1-rho)).
func (q MD1) Lq() float64 {
	rho := q.Rho()
	if rho >= 1 {
		return math.Inf(1)
	}
	return rho * rho / (2 * (1 - rho))
}

// L returns the mean number in system, Lq + rho.
func (q MD1) L() float64 {
	lq := q.Lq()
	if math.IsInf(lq, 1) {
		return lq
	}
	return lq + q.Rho()
}

// Wq returns the mean waiting time in the queue, rho/(2 Mu (1-rho)).
func (q MD1) Wq() float64 {
	rho := q.Rho()
	if rho >= 1 {
		return math.Inf(1)
	}
	return rho / (2 * q.Mu * (1 - rho))
}

// W returns the mean sojourn time, Wq + 1/Mu.
func (q MD1) W() float64 {
	wq := q.Wq()
	if math.IsInf(wq, 1) {
		return wq
	}
	return wq + 1/q.Mu
}

// P0 returns the probability the system is empty, 1-rho.
func (q MD1) P0() float64 { return 1 - q.Rho() }
