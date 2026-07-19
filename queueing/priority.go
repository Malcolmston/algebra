package queueing

import "math"

// PriorityMM1 models a single-server queue with Poisson arrivals in several
// priority classes and exponential, class-dependent service. Class 0 has the
// highest priority. Lambda[k] and Mu[k] are the arrival and service rates of
// class k. The type provides the non-preemptive (head-of-line) priority
// results; a steady state exists when the total utilization is below one.
type PriorityMM1 struct {
	Lambda []float64 // per-class arrival rates, index 0 highest priority
	Mu     []float64 // per-class service rates
}

// NewPriorityMM1 constructs a [PriorityMM1] queue from per-class arrival and
// service rates. It returns an error for mismatched or empty slices,
// non-positive rates, or a total utilization at or above one.
func NewPriorityMM1(lambda, mu []float64) (PriorityMM1, error) {
	if len(lambda) != len(mu) || len(lambda) == 0 {
		return PriorityMM1{}, ErrDimension
	}
	total := 0.0
	for i := range lambda {
		if lambda[i] <= 0 || mu[i] <= 0 {
			return PriorityMM1{}, ErrNonPositiveRate
		}
		total += lambda[i] / mu[i]
	}
	if total >= 1 {
		return PriorityMM1{}, ErrUnstable
	}
	l := append([]float64(nil), lambda...)
	m := append([]float64(nil), mu...)
	return PriorityMM1{Lambda: l, Mu: m}, nil
}

// Classes returns the number of priority classes.
func (q PriorityMM1) Classes() int { return len(q.Lambda) }

// ClassRho returns the per-class utilization rho_k = Lambda[k]/Mu[k].
func (q PriorityMM1) ClassRho(k int) float64 {
	if k < 0 || k >= len(q.Lambda) {
		return math.NaN()
	}
	return q.Lambda[k] / q.Mu[k]
}

// Rho returns the total utilization sum_k Lambda[k]/Mu[k].
func (q PriorityMM1) Rho() float64 {
	sum := 0.0
	for i := range q.Lambda {
		sum += q.Lambda[i] / q.Mu[i]
	}
	return sum
}

// meanResidual returns R = sum_i lambda_i E[S_i^2]/2, the aggregate mean
// residual service (delay) with exponential service E[S_i^2]=2/mu_i^2.
func (q PriorityMM1) meanResidual() float64 {
	r := 0.0
	for i := range q.Lambda {
		r += q.Lambda[i] / (q.Mu[i] * q.Mu[i])
	}
	return r
}

// cumRho returns sum_{i=0}^{k} rho_i, the cumulative utilization through class
// k. cumRho(-1) is 0.
func (q PriorityMM1) cumRho(k int) float64 {
	sum := 0.0
	for i := 0; i <= k && i < len(q.Lambda); i++ {
		sum += q.Lambda[i] / q.Mu[i]
	}
	return sum
}

// ClassWq returns the mean waiting time of class k under non-preemptive
// priority, R / ((1-sigma_{k-1})(1-sigma_k)), where sigma_j is the cumulative
// utilization through class j. It returns NaN for an invalid class index.
func (q PriorityMM1) ClassWq(k int) float64 {
	if k < 0 || k >= len(q.Lambda) {
		return math.NaN()
	}
	r := q.meanResidual()
	lower := 1 - q.cumRho(k-1)
	upper := 1 - q.cumRho(k)
	if lower <= 0 || upper <= 0 {
		return math.Inf(1)
	}
	return r / (lower * upper)
}

// ClassW returns the mean sojourn time of class k, ClassWq(k)+1/Mu[k].
func (q PriorityMM1) ClassW(k int) float64 {
	wq := q.ClassWq(k)
	if math.IsNaN(wq) || math.IsInf(wq, 1) {
		return wq
	}
	return wq + 1/q.Mu[k]
}

// ClassLq returns the mean number of class-k customers waiting, Lambda[k]*Wq_k.
func (q PriorityMM1) ClassLq(k int) float64 {
	wq := q.ClassWq(k)
	if math.IsNaN(wq) || math.IsInf(wq, 1) {
		return wq
	}
	return q.Lambda[k] * wq
}

// ClassL returns the mean number of class-k customers in system, Lambda[k]*W_k.
func (q PriorityMM1) ClassL(k int) float64 {
	w := q.ClassW(k)
	if math.IsNaN(w) || math.IsInf(w, 1) {
		return w
	}
	return q.Lambda[k] * w
}

// MeanWait returns the arrival-weighted overall mean waiting time across all
// classes, (sum_k Lambda[k] Wq_k)/(sum_k Lambda[k]).
func (q PriorityMM1) MeanWait() float64 {
	num, den := 0.0, 0.0
	for k := range q.Lambda {
		wq := q.ClassWq(k)
		if math.IsInf(wq, 1) {
			return math.Inf(1)
		}
		num += q.Lambda[k] * wq
		den += q.Lambda[k]
	}
	if den <= 0 {
		return math.NaN()
	}
	return num / den
}
