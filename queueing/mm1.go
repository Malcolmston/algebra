package queueing

import "math"

// MM1 models a single-server queue with Poisson arrivals of rate Lambda,
// exponential service of rate Mu, an infinite waiting room and FCFS discipline
// (Kendall notation M/M/1). A steady state exists when Lambda < Mu.
type MM1 struct {
	Lambda float64 // arrival rate
	Mu     float64 // service rate
}

// NewMM1 constructs an [MM1] queue, returning an error for non-positive rates
// or when the stability condition Lambda < Mu is violated.
func NewMM1(lambda, mu float64) (MM1, error) {
	if lambda <= 0 || mu <= 0 {
		return MM1{}, ErrNonPositiveRate
	}
	if lambda >= mu {
		return MM1{}, ErrUnstable
	}
	return MM1{Lambda: lambda, Mu: mu}, nil
}

// Rho returns the server utilization (traffic intensity) rho = Lambda/Mu, which
// also equals the long-run fraction of time the server is busy.
func (q MM1) Rho() float64 {
	if q.Mu <= 0 {
		return math.NaN()
	}
	return q.Lambda / q.Mu
}

// Stable reports whether the queue admits a steady state (Lambda < Mu).
func (q MM1) Stable() bool { return q.Lambda > 0 && q.Mu > 0 && q.Lambda < q.Mu }

// Utilization returns the fraction of time the server is busy, equal to Rho.
func (q MM1) Utilization() float64 { return q.Rho() }

// IdleProb returns the probability the system is empty, P0 = 1 - rho.
func (q MM1) IdleProb() float64 { return q.P0() }

// P0 returns the steady-state probability of an empty system, 1 - rho.
func (q MM1) P0() float64 {
	rho := q.Rho()
	if math.IsNaN(rho) {
		return math.NaN()
	}
	return 1 - rho
}

// Pn returns the steady-state probability P(N=n) = (1-rho) rho^n of finding n
// customers in the system. It returns 0 for negative n.
func (q MM1) Pn(n int) float64 {
	if n < 0 {
		return 0
	}
	rho := q.Rho()
	if math.IsNaN(rho) {
		return math.NaN()
	}
	return (1 - rho) * math.Pow(rho, float64(n))
}

// PAtLeast returns the probability of at least n customers in the system,
// P(N>=n) = rho^n. It returns 1 for n<=0.
func (q MM1) PAtLeast(n int) float64 {
	if n <= 0 {
		return 1
	}
	rho := q.Rho()
	if math.IsNaN(rho) {
		return math.NaN()
	}
	return math.Pow(rho, float64(n))
}

// L returns the mean number of customers in the system, rho/(1-rho).
func (q MM1) L() float64 {
	rho := q.Rho()
	if rho >= 1 {
		return math.Inf(1)
	}
	return rho / (1 - rho)
}

// Lq returns the mean number of customers waiting in the queue,
// rho^2/(1-rho).
func (q MM1) Lq() float64 {
	rho := q.Rho()
	if rho >= 1 {
		return math.Inf(1)
	}
	return rho * rho / (1 - rho)
}

// Ls returns the mean number of customers in service, equal to rho.
func (q MM1) Ls() float64 { return q.Rho() }

// W returns the mean sojourn time (response time) in the system,
// 1/(Mu-Lambda).
func (q MM1) W() float64 {
	if q.Mu <= q.Lambda {
		return math.Inf(1)
	}
	return 1 / (q.Mu - q.Lambda)
}

// Wq returns the mean waiting time in the queue, Lambda/(Mu(Mu-Lambda)).
func (q MM1) Wq() float64 {
	if q.Mu <= q.Lambda {
		return math.Inf(1)
	}
	return q.Lambda / (q.Mu * (q.Mu - q.Lambda))
}

// Throughput returns the departure rate of the system, which in steady state
// equals the arrival rate Lambda.
func (q MM1) Throughput() float64 { return q.Lambda }

// VarianceN returns the variance of the number in system, rho/(1-rho)^2.
func (q MM1) VarianceN() float64 {
	rho := q.Rho()
	if rho >= 1 {
		return math.Inf(1)
	}
	d := 1 - rho
	return rho / (d * d)
}

// WaitTailProb returns P(W>t), the probability the sojourn time exceeds t,
// which is exp(-(Mu-Lambda)t) for the M/M/1 FCFS queue. It returns 1 for t<=0.
func (q MM1) WaitTailProb(t float64) float64 {
	if q.Mu <= q.Lambda {
		return math.NaN()
	}
	if t <= 0 {
		return 1
	}
	return math.Exp(-(q.Mu - q.Lambda) * t)
}

// WaitCDF returns P(W<=t), the sojourn-time distribution 1-exp(-(Mu-Lambda)t).
func (q MM1) WaitCDF(t float64) float64 {
	p := q.WaitTailProb(t)
	if math.IsNaN(p) {
		return math.NaN()
	}
	return 1 - p
}

// WaitqTailProb returns P(Wq>t), the probability the queueing delay exceeds t,
// which is rho*exp(-(Mu-Lambda)t). It returns rho for t<=0 (the probability of
// a positive delay).
func (q MM1) WaitqTailProb(t float64) float64 {
	rho := q.Rho()
	if q.Mu <= q.Lambda {
		return math.NaN()
	}
	if t < 0 {
		t = 0
	}
	return rho * math.Exp(-(q.Mu-q.Lambda)*t)
}

// WaitqCDF returns P(Wq<=t) = 1 - rho*exp(-(Mu-Lambda)t), the queueing-delay
// distribution, which has an atom 1-rho at zero.
func (q MM1) WaitqCDF(t float64) float64 {
	p := q.WaitqTailProb(t)
	if math.IsNaN(p) {
		return math.NaN()
	}
	return 1 - p
}

// WaitPercentile returns the sojourn-time quantile t such that P(W<=t)=p, i.e.
// -ln(1-p)/(Mu-Lambda). It returns NaN for p outside [0,1).
func (q MM1) WaitPercentile(p float64) float64 {
	if p < 0 || p >= 1 || q.Mu <= q.Lambda {
		return math.NaN()
	}
	return -math.Log(1-p) / (q.Mu - q.Lambda)
}

// MeanBusyPeriod returns the mean length of a server busy period, 1/(Mu-Lambda).
func (q MM1) MeanBusyPeriod() float64 {
	if q.Mu <= q.Lambda {
		return math.Inf(1)
	}
	return 1 / (q.Mu - q.Lambda)
}

// MeanBusyCycle returns the mean length of a busy cycle (idle period plus busy
// period), 1/(Lambda(1-rho)).
func (q MM1) MeanBusyCycle() float64 {
	rho := q.Rho()
	if rho >= 1 || q.Lambda <= 0 {
		return math.Inf(1)
	}
	return 1 / (q.Lambda * (1 - rho))
}
