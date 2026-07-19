package queueing

import "math"

// MG1 models a single-server queue with Poisson arrivals of rate Lambda and a
// general service-time distribution specified by its mean ServiceMean and
// variance ServiceVar (Kendall notation M/G/1). Mean-value results follow from
// the Pollaczek–Khinchine formula. A steady state exists when the utilization
// rho = Lambda*ServiceMean is below one.
type MG1 struct {
	Lambda      float64 // arrival rate
	ServiceMean float64 // mean service time E[S]
	ServiceVar  float64 // variance of service time Var[S]
}

// NewMG1 constructs an [MG1] queue from an arrival rate and the mean and
// variance of the service time. It returns an error for a non-positive arrival
// rate or mean service time, a negative variance, or when the stability
// condition Lambda*ServiceMean < 1 is violated.
func NewMG1(lambda, serviceMean, serviceVar float64) (MG1, error) {
	if lambda <= 0 || serviceMean <= 0 {
		return MG1{}, ErrNonPositiveRate
	}
	if serviceVar < 0 {
		return MG1{}, ErrNegative
	}
	if lambda*serviceMean >= 1 {
		return MG1{}, ErrUnstable
	}
	return MG1{Lambda: lambda, ServiceMean: serviceMean, ServiceVar: serviceVar}, nil
}

// NewMG1FromSCV constructs an [MG1] queue from the arrival rate, mean service
// time and squared coefficient of variation scv = Var[S]/E[S]^2.
func NewMG1FromSCV(lambda, serviceMean, scv float64) (MG1, error) {
	if scv < 0 {
		return MG1{}, ErrNegative
	}
	return NewMG1(lambda, serviceMean, scv*serviceMean*serviceMean)
}

// Rho returns the utilization rho = Lambda*E[S].
func (q MG1) Rho() float64 { return q.Lambda * q.ServiceMean }

// Stable reports whether the utilization is below one.
func (q MG1) Stable() bool { return q.Rho() < 1 && q.Lambda > 0 && q.ServiceMean > 0 }

// ServiceSCV returns the squared coefficient of variation of the service time,
// Var[S]/E[S]^2.
func (q MG1) ServiceSCV() float64 {
	return SquaredCoefficientOfVariation(q.ServiceMean, q.ServiceVar)
}

// ServiceRate returns the nominal service rate Mu = 1/E[S].
func (q MG1) ServiceRate() float64 { return 1 / q.ServiceMean }

// Lq returns the mean number waiting in the queue from the Pollaczek–Khinchine
// formula, (Lambda^2 Var[S] + rho^2) / (2(1-rho)).
func (q MG1) Lq() float64 {
	rho := q.Rho()
	if rho >= 1 {
		return math.Inf(1)
	}
	num := q.Lambda*q.Lambda*q.ServiceVar + rho*rho
	return num / (2 * (1 - rho))
}

// L returns the mean number in system, Lq + rho.
func (q MG1) L() float64 {
	lq := q.Lq()
	if math.IsInf(lq, 1) {
		return lq
	}
	return lq + q.Rho()
}

// Wq returns the mean waiting time in the queue, Lq/Lambda.
func (q MG1) Wq() float64 {
	lq := q.Lq()
	if math.IsInf(lq, 1) {
		return lq
	}
	return lq / q.Lambda
}

// W returns the mean sojourn time, Wq + E[S].
func (q MG1) W() float64 {
	wq := q.Wq()
	if math.IsInf(wq, 1) {
		return wq
	}
	return wq + q.ServiceMean
}

// P0 returns the probability that the system is empty, 1-rho, which holds for
// any M/G/1 queue.
func (q MG1) P0() float64 { return 1 - q.Rho() }

// MeanResidualService returns the mean residual (forward-recurrence) service
// time seen by an arrival at a busy instant, E[S^2]/(2 E[S]).
func (q MG1) MeanResidualService() float64 {
	es2 := q.ServiceVar + q.ServiceMean*q.ServiceMean
	return es2 / (2 * q.ServiceMean)
}

// PollaczekKhinchineLq returns the Pollaczek–Khinchine mean queue length for an
// M/G/1 queue with utilization rho and service-time squared coefficient of
// variation scv: rho^2(1+scv)/(2(1-rho)). It returns +Inf for rho>=1 and NaN
// for invalid inputs.
func PollaczekKhinchineLq(rho, scv float64) float64 {
	if rho < 0 || scv < 0 {
		return math.NaN()
	}
	if rho >= 1 {
		return math.Inf(1)
	}
	return rho * rho * (1 + scv) / (2 * (1 - rho))
}

// PollaczekKhinchineWq returns the Pollaczek–Khinchine mean waiting time
// rho*(1+scv)/(2(1-rho)) * E[S], where E[S] is the mean service time. It
// returns +Inf for rho>=1 and NaN for invalid inputs.
func PollaczekKhinchineWq(rho, scv, serviceMean float64) float64 {
	if rho < 0 || scv < 0 || serviceMean < 0 {
		return math.NaN()
	}
	if rho >= 1 {
		return math.Inf(1)
	}
	return rho * (1 + scv) / (2 * (1 - rho)) * serviceMean
}
