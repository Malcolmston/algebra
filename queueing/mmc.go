package queueing

import "math"

// MMc models a multi-server queue with Poisson arrivals of rate Lambda, C
// identical exponential servers each of rate Mu, an infinite waiting room and
// FCFS discipline (Kendall notation M/M/c). A steady state exists when the
// offered load Lambda/Mu is less than C.
type MMc struct {
	Lambda float64 // arrival rate
	Mu     float64 // per-server service rate
	C      int     // number of servers
}

// NewMMc constructs an [MMc] queue, returning an error for non-positive rates,
// a non-positive server count, or when the stability condition Lambda < C*Mu
// is violated.
func NewMMc(lambda, mu float64, c int) (MMc, error) {
	if lambda <= 0 || mu <= 0 {
		return MMc{}, ErrNonPositiveRate
	}
	if c <= 0 {
		return MMc{}, ErrServers
	}
	if lambda >= float64(c)*mu {
		return MMc{}, ErrUnstable
	}
	return MMc{Lambda: lambda, Mu: mu, C: c}, nil
}

// OfferedLoad returns the offered load a = Lambda/Mu in erlangs, the mean
// number of busy servers in steady state.
func (q MMc) OfferedLoad() float64 { return OfferedLoad(q.Lambda, q.Mu) }

// Rho returns the per-server utilization rho = Lambda/(C*Mu).
func (q MMc) Rho() float64 { return Utilization(q.Lambda, q.Mu, q.C) }

// Stable reports whether the queue admits a steady state (Lambda < C*Mu).
func (q MMc) Stable() bool {
	return q.Lambda > 0 && q.Mu > 0 && q.C > 0 && q.Lambda < float64(q.C)*q.Mu
}

// Utilization returns the per-server utilization, equal to Rho.
func (q MMc) Utilization() float64 { return q.Rho() }

// P0 returns the steady-state probability that the system is empty.
func (q MMc) P0() float64 {
	a := q.OfferedLoad()
	rho := q.Rho()
	if math.IsNaN(a) || rho >= 1 {
		return math.NaN()
	}
	sum := 0.0
	term := 1.0 // a^0/0!
	for n := 0; n < q.C; n++ {
		sum += term
		term *= a / float64(n+1)
	}
	// term now equals a^C / C!
	sum += term / (1 - rho)
	return 1 / sum
}

// Pn returns the steady-state probability P(N=n) of finding n customers in the
// system. It returns 0 for negative n.
func (q MMc) Pn(n int) float64 {
	if n < 0 {
		return 0
	}
	p0 := q.P0()
	a := q.OfferedLoad()
	if math.IsNaN(p0) {
		return math.NaN()
	}
	if n <= q.C {
		return powFactRatio(a, n) * p0
	}
	// a^C/C! * (a/C)^(n-C)
	base := powFactRatio(a, q.C)
	rho := a / float64(q.C)
	return base * math.Pow(rho, float64(n-q.C)) * p0
}

// ErlangC returns the Erlang-C delay probability, the steady-state probability
// that an arriving customer must wait (finds all C servers busy).
func (q MMc) ErlangC() float64 {
	p0 := q.P0()
	a := q.OfferedLoad()
	rho := q.Rho()
	if math.IsNaN(p0) || rho >= 1 {
		return math.NaN()
	}
	base := powFactRatio(a, q.C) // a^C/C!
	return base / (1 - rho) * p0
}

// ProbWait returns the probability that an arriving customer is delayed, equal
// to the Erlang-C value.
func (q MMc) ProbWait() float64 { return q.ErlangC() }

// Lq returns the mean number of customers waiting in the queue,
// ErlangC * rho/(1-rho).
func (q MMc) Lq() float64 {
	c := q.ErlangC()
	rho := q.Rho()
	if math.IsNaN(c) {
		return math.NaN()
	}
	return c * rho / (1 - rho)
}

// L returns the mean number of customers in the system, Lq + a.
func (q MMc) L() float64 {
	lq := q.Lq()
	if math.IsNaN(lq) {
		return math.NaN()
	}
	return lq + q.OfferedLoad()
}

// Wq returns the mean waiting time in the queue, Lq/Lambda.
func (q MMc) Wq() float64 {
	lq := q.Lq()
	if math.IsNaN(lq) || q.Lambda <= 0 {
		return math.NaN()
	}
	return lq / q.Lambda
}

// W returns the mean sojourn time in the system, Wq + 1/Mu.
func (q MMc) W() float64 {
	wq := q.Wq()
	if math.IsNaN(wq) {
		return math.NaN()
	}
	return wq + 1/q.Mu
}

// Throughput returns the departure rate, which equals Lambda in steady state.
func (q MMc) Throughput() float64 { return q.Lambda }

// MeanBusyServers returns the mean number of busy servers, equal to the offered
// load a = Lambda/Mu.
func (q MMc) MeanBusyServers() float64 { return q.OfferedLoad() }

// WaitqTailProb returns P(Wq>t), the probability the queueing delay exceeds t,
// which is ErlangC * exp(-(C*Mu-Lambda) t). It returns the delay probability
// for t<=0.
func (q MMc) WaitqTailProb(t float64) float64 {
	c := q.ErlangC()
	if math.IsNaN(c) {
		return math.NaN()
	}
	if t < 0 {
		t = 0
	}
	rate := float64(q.C)*q.Mu - q.Lambda
	return c * math.Exp(-rate*t)
}

// WaitqCDF returns P(Wq<=t), the queueing-delay distribution.
func (q MMc) WaitqCDF(t float64) float64 {
	p := q.WaitqTailProb(t)
	if math.IsNaN(p) {
		return math.NaN()
	}
	return 1 - p
}

// WaitqPercentile returns the queueing-delay quantile t such that P(Wq<=t)=p.
// It returns 0 when p is at or below the no-wait probability 1-ErlangC, and NaN
// for p outside [0,1).
func (q MMc) WaitqPercentile(p float64) float64 {
	c := q.ErlangC()
	if p < 0 || p >= 1 || math.IsNaN(c) {
		return math.NaN()
	}
	if p <= 1-c {
		return 0
	}
	rate := float64(q.C)*q.Mu - q.Lambda
	return math.Log(c/(1-p)) / rate
}
