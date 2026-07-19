package queueing

import "math"

// ErlangB returns the Erlang-B loss probability (blocking probability) of an
// M/M/c/c loss system with c servers and offered load a erlangs. It is computed
// with the numerically stable recurrence B(0,a)=1,
// B(n,a) = a*B(n-1,a) / (n + a*B(n-1,a)). It returns NaN for a negative load or
// non-positive server count.
func ErlangB(c int, a float64) float64 {
	if c <= 0 || a < 0 {
		return math.NaN()
	}
	b := 1.0
	for n := 1; n <= c; n++ {
		b = a * b / (float64(n) + a*b)
	}
	return b
}

// ErlangC returns the Erlang-C delay probability, the probability that an
// arriving customer must wait in an M/M/c queue with c servers and offered load
// a erlangs. It is derived from the Erlang-B value via
// C = B / (1 - (a/c)(1-B)). It returns NaN when the offered load is not below
// the number of servers (the queue would be unstable).
func ErlangC(c int, a float64) float64 {
	if c <= 0 || a < 0 {
		return math.NaN()
	}
	if a >= float64(c) {
		return math.NaN()
	}
	b := ErlangB(c, a)
	rho := a / float64(c)
	return b / (1 - rho*(1-b))
}

// ErlangCFromB returns the Erlang-C delay probability computed from a
// precomputed Erlang-B value b for c servers and offered load a. It returns NaN
// for an unstable system.
func ErlangCFromB(c int, a, b float64) float64 {
	if c <= 0 || a < 0 || a >= float64(c) {
		return math.NaN()
	}
	rho := a / float64(c)
	return b / (1 - rho*(1-b))
}

// ErlangBSeries returns the Erlang-B blocking probabilities for all server
// counts from 0 through c at offered load a, produced by a single sweep of the
// stable recurrence. Element i holds B(i,a); B(0,a)=1. It returns nil for a
// negative load or negative c.
func ErlangBSeries(c int, a float64) []float64 {
	if c < 0 || a < 0 {
		return nil
	}
	out := make([]float64, c+1)
	out[0] = 1
	b := 1.0
	for n := 1; n <= c; n++ {
		b = a * b / (float64(n) + a*b)
		out[n] = b
	}
	return out
}

// ErlangBServersFor returns the smallest number of servers whose Erlang-B
// blocking probability at offered load a is at or below the target grade of
// service gos (a probability in (0,1]). It returns 0 for a=0 and -1 if no
// value up to maxServers achieves the target.
func ErlangBServersFor(a, gos float64, maxServers int) int {
	if a < 0 || gos <= 0 || gos > 1 || maxServers < 0 {
		return -1
	}
	if a == 0 {
		return 0
	}
	b := 1.0
	for n := 1; n <= maxServers; n++ {
		b = a * b / (float64(n) + a*b)
		if b <= gos {
			return n
		}
	}
	return -1
}

// ErlangCServersFor returns the smallest number of servers whose Erlang-C delay
// probability at offered load a is at or below the target probability p. The
// search starts at ceil(a)+1 to guarantee stability and stops at maxServers,
// returning -1 on failure.
func ErlangCServersFor(a, p float64, maxServers int) int {
	if a < 0 || p <= 0 || p > 1 || maxServers < 0 {
		return -1
	}
	start := int(math.Floor(a)) + 1
	for n := start; n <= maxServers; n++ {
		c := ErlangC(n, a)
		if !math.IsNaN(c) && c <= p {
			return n
		}
	}
	return -1
}

// CarriedLoad returns the load carried by an M/M/c/c loss system, a(1-B), where
// B is the Erlang-B blocking probability at offered load a.
func CarriedLoad(c int, a float64) float64 {
	b := ErlangB(c, a)
	if math.IsNaN(b) {
		return math.NaN()
	}
	return a * (1 - b)
}

// OfferedLoadFromCarried returns the offered load implied by a measured carried
// load and blocking probability, carried/(1-blocking). It returns NaN for a
// blocking probability outside [0,1).
func OfferedLoadFromCarried(carried, blocking float64) float64 {
	if blocking < 0 || blocking >= 1 || carried < 0 {
		return math.NaN()
	}
	return carried / (1 - blocking)
}

// ErlangCWaitTailProb returns the probability that the queueing delay in an
// M/M/c queue exceeds t, C * exp(-(c*mu - lambda) t), where C is the Erlang-C
// value at offered load lambda/mu. It returns NaN for an unstable system.
func ErlangCWaitTailProb(c int, lambda, mu, t float64) float64 {
	if c <= 0 || lambda <= 0 || mu <= 0 {
		return math.NaN()
	}
	a := lambda / mu
	cc := ErlangC(c, a)
	if math.IsNaN(cc) {
		return math.NaN()
	}
	if t < 0 {
		t = 0
	}
	rate := float64(c)*mu - lambda
	return cc * math.Exp(-rate*t)
}

// ErlangCMeanWait returns the mean queueing delay in an M/M/c queue with c
// servers, arrival rate lambda and per-server service rate mu, computed from
// the Erlang-C formula as C / (c*mu - lambda). It returns NaN for an unstable
// system.
func ErlangCMeanWait(c int, lambda, mu float64) float64 {
	if c <= 0 || lambda <= 0 || mu <= 0 {
		return math.NaN()
	}
	a := lambda / mu
	cc := ErlangC(c, a)
	if math.IsNaN(cc) {
		return math.NaN()
	}
	return cc / (float64(c)*mu - lambda)
}
