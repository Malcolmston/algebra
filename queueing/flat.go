package queueing

import "math"

// This file provides a flat, free-function API mirroring the model types for
// callers who prefer to pass rates directly rather than construct a value.

// MM1Rho returns the M/M/1 utilization lambda/mu.
func MM1Rho(lambda, mu float64) float64 { return MM1{Lambda: lambda, Mu: mu}.Rho() }

// MM1L returns the mean number in an M/M/1 system, rho/(1-rho).
func MM1L(lambda, mu float64) float64 { return MM1{Lambda: lambda, Mu: mu}.L() }

// MM1Lq returns the mean M/M/1 queue length rho^2/(1-rho).
func MM1Lq(lambda, mu float64) float64 { return MM1{Lambda: lambda, Mu: mu}.Lq() }

// MM1W returns the mean M/M/1 sojourn time 1/(mu-lambda).
func MM1W(lambda, mu float64) float64 { return MM1{Lambda: lambda, Mu: mu}.W() }

// MM1Wq returns the mean M/M/1 waiting time lambda/(mu(mu-lambda)).
func MM1Wq(lambda, mu float64) float64 { return MM1{Lambda: lambda, Mu: mu}.Wq() }

// MM1P0 returns the M/M/1 empty-system probability 1-rho.
func MM1P0(lambda, mu float64) float64 { return MM1{Lambda: lambda, Mu: mu}.P0() }

// MM1Pn returns the M/M/1 state probability (1-rho)rho^n.
func MM1Pn(lambda, mu float64, n int) float64 { return MM1{Lambda: lambda, Mu: mu}.Pn(n) }

// MM1WaitTailProb returns P(W>t) for an M/M/1 queue.
func MM1WaitTailProb(lambda, mu, t float64) float64 {
	return MM1{Lambda: lambda, Mu: mu}.WaitTailProb(t)
}

// MMcL returns the mean number in an M/M/c system.
func MMcL(lambda, mu float64, c int) float64 { return MMc{Lambda: lambda, Mu: mu, C: c}.L() }

// MMcLq returns the mean M/M/c queue length.
func MMcLq(lambda, mu float64, c int) float64 { return MMc{Lambda: lambda, Mu: mu, C: c}.Lq() }

// MMcW returns the mean M/M/c sojourn time.
func MMcW(lambda, mu float64, c int) float64 { return MMc{Lambda: lambda, Mu: mu, C: c}.W() }

// MMcWq returns the mean M/M/c waiting time.
func MMcWq(lambda, mu float64, c int) float64 { return MMc{Lambda: lambda, Mu: mu, C: c}.Wq() }

// MMcP0 returns the M/M/c empty-system probability.
func MMcP0(lambda, mu float64, c int) float64 { return MMc{Lambda: lambda, Mu: mu, C: c}.P0() }

// MMcPn returns the M/M/c state probability P(N=n).
func MMcPn(lambda, mu float64, c, n int) float64 { return MMc{Lambda: lambda, Mu: mu, C: c}.Pn(n) }

// MMcDelayProb returns the Erlang-C delay probability for an M/M/c queue.
func MMcDelayProb(lambda, mu float64, c int) float64 {
	return MMc{Lambda: lambda, Mu: mu, C: c}.ErlangC()
}

// MMInfL returns the mean number in an M/M/infinity system, lambda/mu.
func MMInfL(lambda, mu float64) float64 { return MMInf{Lambda: lambda, Mu: mu}.L() }

// MMInfPn returns the M/M/infinity state probability (Poisson mass).
func MMInfPn(lambda, mu float64, n int) float64 { return MMInf{Lambda: lambda, Mu: mu}.Pn(n) }

// MM1KBlockingProb returns the M/M/1/K loss probability P(N=K).
func MM1KBlockingProb(lambda, mu float64, k int) float64 {
	return MM1K{Lambda: lambda, Mu: mu, K: k}.BlockingProb()
}

// MM1KL returns the mean number in an M/M/1/K system.
func MM1KL(lambda, mu float64, k int) float64 { return MM1K{Lambda: lambda, Mu: mu, K: k}.L() }

// MMcKBlockingProb returns the M/M/c/K loss probability P(N=K).
func MMcKBlockingProb(lambda, mu float64, c, k int) float64 {
	return MMcK{Lambda: lambda, Mu: mu, C: c, K: k}.BlockingProb()
}

// MMcKL returns the mean number in an M/M/c/K system.
func MMcKL(lambda, mu float64, c, k int) float64 {
	return MMcK{Lambda: lambda, Mu: mu, C: c, K: k}.L()
}

// MG1Lq returns the M/G/1 Pollaczek–Khinchine mean queue length from the
// arrival rate, mean service time and service-time variance.
func MG1Lq(lambda, serviceMean, serviceVar float64) float64 {
	return MG1{Lambda: lambda, ServiceMean: serviceMean, ServiceVar: serviceVar}.Lq()
}

// MG1Wq returns the M/G/1 mean waiting time from the arrival rate, mean service
// time and service-time variance.
func MG1Wq(lambda, serviceMean, serviceVar float64) float64 {
	return MG1{Lambda: lambda, ServiceMean: serviceMean, ServiceVar: serviceVar}.Wq()
}

// MD1Lq returns the M/D/1 mean queue length rho^2/(2(1-rho)).
func MD1Lq(lambda, mu float64) float64 { return MD1{Lambda: lambda, Mu: mu}.Lq() }

// MD1Wq returns the M/D/1 mean waiting time.
func MD1Wq(lambda, mu float64) float64 { return MD1{Lambda: lambda, Mu: mu}.Wq() }

// ExponentialPDF returns the exponential probability density rate*exp(-rate*x)
// for x>=0. It returns NaN for a non-positive rate and 0 for x<0.
func ExponentialPDF(rate, x float64) float64 {
	if rate <= 0 {
		return math.NaN()
	}
	if x < 0 {
		return 0
	}
	return rate * math.Exp(-rate*x)
}

// ExponentialCDF returns the exponential cumulative probability
// 1-exp(-rate*x) for x>=0. It returns NaN for a non-positive rate and 0 for
// x<0.
func ExponentialCDF(rate, x float64) float64 {
	if rate <= 0 {
		return math.NaN()
	}
	if x < 0 {
		return 0
	}
	return 1 - math.Exp(-rate*x)
}

// ExponentialMean returns the mean 1/rate of an exponential distribution.
func ExponentialMean(rate float64) float64 {
	if rate <= 0 {
		return math.NaN()
	}
	return 1 / rate
}

// ErlangPDF returns the density of an Erlang-k distribution (sum of k iid
// exponentials of rate lambda) at x>=0. It returns NaN for invalid parameters.
func ErlangPDF(k int, lambda, x float64) float64 {
	if k <= 0 || lambda <= 0 {
		return math.NaN()
	}
	if x < 0 {
		return 0
	}
	logp := float64(k)*math.Log(lambda) + float64(k-1)*math.Log(x) - lambda*x - LogFactorial(k-1)
	return math.Exp(logp)
}

// ErlangCDF returns the cumulative distribution of an Erlang-k distribution at
// x>=0, 1 - e^{-lambda x} sum_{i=0}^{k-1} (lambda x)^i / i!. It returns NaN for
// invalid parameters.
func ErlangCDF(k int, lambda, x float64) float64 {
	if k <= 0 || lambda <= 0 {
		return math.NaN()
	}
	if x <= 0 {
		return 0
	}
	lx := lambda * x
	sum := 0.0
	term := 1.0
	for i := 0; i < k; i++ {
		sum += term
		term *= lx / float64(i+1)
	}
	return 1 - math.Exp(-lx)*sum
}

// ErlangMean returns the mean k/lambda of an Erlang-k distribution.
func ErlangMean(k int, lambda float64) float64 {
	if k <= 0 || lambda <= 0 {
		return math.NaN()
	}
	return float64(k) / lambda
}

// TrafficFromMetrics returns the arrival rate implied by a mean number in
// system L and a mean sojourn time W via Little's law, an alias for
// [LittleLambda].
func TrafficFromMetrics(l, w float64) float64 { return LittleLambda(l, w) }
