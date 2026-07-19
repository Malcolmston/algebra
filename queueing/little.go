package queueing

import "math"

// LittleL returns the mean number in system L = lambda*W from the arrival rate
// and mean sojourn time via Little's law. It returns NaN for negative inputs.
func LittleL(lambda, w float64) float64 {
	if lambda < 0 || w < 0 {
		return math.NaN()
	}
	return lambda * w
}

// LittleW returns the mean sojourn time W = L/lambda from the mean number in
// system and the arrival rate via Little's law. It returns NaN for a
// non-positive arrival rate or negative L.
func LittleW(l, lambda float64) float64 {
	if lambda <= 0 || l < 0 {
		return math.NaN()
	}
	return l / lambda
}

// LittleLambda returns the arrival rate lambda = L/W implied by Little's law
// from the mean number in system and mean sojourn time. It returns NaN for a
// non-positive W or negative L.
func LittleLambda(l, w float64) float64 {
	if w <= 0 || l < 0 {
		return math.NaN()
	}
	return l / w
}

// LittleLq returns the mean queue length Lq = lambda*Wq from the arrival rate
// and mean waiting time via Little's law applied to the queue. It returns NaN
// for negative inputs.
func LittleLq(lambda, wq float64) float64 {
	if lambda < 0 || wq < 0 {
		return math.NaN()
	}
	return lambda * wq
}

// LittleWq returns the mean waiting time in queue Wq = Lq/lambda from the mean
// queue length and arrival rate via Little's law. It returns NaN for a
// non-positive arrival rate or negative Lq.
func LittleWq(lq, lambda float64) float64 {
	if lambda <= 0 || lq < 0 {
		return math.NaN()
	}
	return lq / lambda
}

// SojournFromWait returns the mean sojourn time W = Wq + 1/mu obtained by
// adding a mean service time to the mean waiting time. It returns NaN for a
// non-positive service rate or negative Wq.
func SojournFromWait(wq, mu float64) float64 {
	if mu <= 0 || wq < 0 {
		return math.NaN()
	}
	return wq + 1/mu
}

// WaitFromSojourn returns the mean waiting time Wq = W - 1/mu obtained by
// subtracting the mean service time from the mean sojourn time. It returns NaN
// for a non-positive service rate or negative W.
func WaitFromSojourn(w, mu float64) float64 {
	if mu <= 0 || w < 0 {
		return math.NaN()
	}
	return w - 1/mu
}
