package queueing

import "math"

// GG1 models a general single-server queue with renewal arrivals and general
// service (Kendall notation G/G/1). Interarrival and service times are
// summarized by their rates and squared coefficients of variation Ca2 and Cs2.
// Exact results are not available, so the type provides the standard heavy-
// traffic diffusion approximations (Kingman and Allen–Cunneen).
type GG1 struct {
	Lambda      float64 // arrival rate
	ServiceMean float64 // mean service time E[S]
	Ca2         float64 // SCV of interarrival time
	Cs2         float64 // SCV of service time
}

// NewGG1 constructs a [GG1] queue from the arrival rate, mean service time and
// the squared coefficients of variation of the interarrival and service times.
// It returns an error for non-positive rates, negative SCVs or when the
// utilization Lambda*ServiceMean is not below one.
func NewGG1(lambda, serviceMean, ca2, cs2 float64) (GG1, error) {
	if lambda <= 0 || serviceMean <= 0 {
		return GG1{}, ErrNonPositiveRate
	}
	if ca2 < 0 || cs2 < 0 {
		return GG1{}, ErrNegative
	}
	if lambda*serviceMean >= 1 {
		return GG1{}, ErrUnstable
	}
	return GG1{Lambda: lambda, ServiceMean: serviceMean, Ca2: ca2, Cs2: cs2}, nil
}

// Rho returns the utilization rho = Lambda*E[S].
func (q GG1) Rho() float64 { return q.Lambda * q.ServiceMean }

// Wq returns the Kingman (VUT) approximation of the mean waiting time,
// (rho/(1-rho)) * ((Ca2+Cs2)/2) * E[S].
func (q GG1) Wq() float64 {
	rho := q.Rho()
	if rho >= 1 {
		return math.Inf(1)
	}
	return (rho / (1 - rho)) * ((q.Ca2 + q.Cs2) / 2) * q.ServiceMean
}

// Lq returns the approximate mean queue length, Lambda*Wq.
func (q GG1) Lq() float64 {
	wq := q.Wq()
	if math.IsInf(wq, 1) {
		return wq
	}
	return q.Lambda * wq
}

// W returns the approximate mean sojourn time, Wq + E[S].
func (q GG1) W() float64 {
	wq := q.Wq()
	if math.IsInf(wq, 1) {
		return wq
	}
	return wq + q.ServiceMean
}

// L returns the approximate mean number in system, Lq + rho.
func (q GG1) L() float64 {
	lq := q.Lq()
	if math.IsInf(lq, 1) {
		return lq
	}
	return lq + q.Rho()
}

// KingmanWq returns the Kingman heavy-traffic approximation of the mean G/G/1
// waiting time given the utilization rho, the interarrival and service SCVs and
// the mean service time: (rho/(1-rho)) * ((ca2+cs2)/2) * serviceMean. It
// returns +Inf for rho>=1 and NaN for invalid inputs.
func KingmanWq(rho, ca2, cs2, serviceMean float64) float64 {
	if rho < 0 || ca2 < 0 || cs2 < 0 || serviceMean < 0 {
		return math.NaN()
	}
	if rho >= 1 {
		return math.Inf(1)
	}
	return (rho / (1 - rho)) * ((ca2 + cs2) / 2) * serviceMean
}

// AllenCunneenWq returns the Allen–Cunneen approximation of the mean waiting
// time in a G/G/c queue with arrival rate lambda, per-server service rate mu, c
// servers and interarrival/service SCVs ca2 and cs2. It combines the exact
// M/M/c waiting time with the variability factor (ca2+cs2)/2. It returns +Inf
// when the queue is unstable and NaN for invalid inputs.
func AllenCunneenWq(lambda, mu float64, c int, ca2, cs2 float64) float64 {
	if lambda <= 0 || mu <= 0 || c <= 0 || ca2 < 0 || cs2 < 0 {
		return math.NaN()
	}
	if lambda >= float64(c)*mu {
		return math.Inf(1)
	}
	mmc := MMc{Lambda: lambda, Mu: mu, C: c}
	return mmc.Wq() * (ca2 + cs2) / 2
}

// AllenCunneenLq returns the Allen–Cunneen approximation of the mean queue
// length in a G/G/c queue, lambda times [AllenCunneenWq].
func AllenCunneenLq(lambda, mu float64, c int, ca2, cs2 float64) float64 {
	wq := AllenCunneenWq(lambda, mu, c, ca2, cs2)
	if math.IsNaN(wq) || math.IsInf(wq, 1) {
		return wq
	}
	return lambda * wq
}
