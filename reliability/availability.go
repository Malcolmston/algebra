package reliability

import "math"

// Availability returns the steady-state availability MTBF/(MTBF+MTTR) of a
// repairable item from its mean time between failures and mean time to repair.
func Availability(mtbf, mttr float64) float64 {
	if mtbf < 0 || mttr < 0 || mtbf+mttr == 0 {
		return math.NaN()
	}
	return mtbf / (mtbf + mttr)
}

// InherentAvailability returns MTTF/(MTTF+MTTR), the availability that accounts
// only for corrective maintenance under idealized support conditions.
func InherentAvailability(mttf, mttr float64) float64 {
	if mttf < 0 || mttr < 0 || mttf+mttr == 0 {
		return math.NaN()
	}
	return mttf / (mttf + mttr)
}

// OperationalAvailability returns uptime/(uptime+downtime), the availability
// observed in operation including all sources of downtime.
func OperationalAvailability(uptime, downtime float64) float64 {
	if uptime < 0 || downtime < 0 || uptime+downtime == 0 {
		return math.NaN()
	}
	return uptime / (uptime + downtime)
}

// AchievedAvailability returns MTBM/(MTBM+MMT), the availability based on the
// mean time between maintenance actions (both corrective and preventive) and
// the mean maintenance time.
func AchievedAvailability(mtbm, mmt float64) float64 {
	if mtbm < 0 || mmt < 0 || mtbm+mmt == 0 {
		return math.NaN()
	}
	return mtbm / (mtbm + mmt)
}

// Unavailability returns 1 minus the availability computed from MTBF and MTTR.
func Unavailability(mtbf, mttr float64) float64 {
	a := Availability(mtbf, mttr)
	if math.IsNaN(a) {
		return math.NaN()
	}
	return 1 - a
}

// SteadyStateAvailability returns the limiting availability μ/(λ+μ) of a
// two-state (up/down) Markov model with constant failure rate lambda and
// constant repair rate mu.
func SteadyStateAvailability(lambda, mu float64) float64 {
	if lambda < 0 || mu < 0 || lambda+mu == 0 {
		return math.NaN()
	}
	return mu / (lambda + mu)
}

// InstantaneousAvailability returns the availability at time t of a two-state
// Markov model that starts operational:
// A(t)=μ/(λ+μ)+λ/(λ+μ)·e^{-(λ+μ)t}.
func InstantaneousAvailability(t, lambda, mu float64) float64 {
	if t < 0 || lambda < 0 || mu < 0 || lambda+mu == 0 {
		return math.NaN()
	}
	s := lambda + mu
	return mu/s + (lambda/s)*math.Exp(-s*t)
}

// AvailabilitySeriesSystem returns the availability of a series system, the
// product of its component availabilities.
func AvailabilitySeriesSystem(avails ...float64) float64 {
	prod := 1.0
	for _, a := range avails {
		if a < 0 || a > 1 {
			return math.NaN()
		}
		prod *= a
	}
	return prod
}

// AvailabilityParallelSystem returns the availability of a parallel system,
// 1-∏(1-aᵢ).
func AvailabilityParallelSystem(avails ...float64) float64 {
	q := 1.0
	for _, a := range avails {
		if a < 0 || a > 1 {
			return math.NaN()
		}
		q *= 1 - a
	}
	return 1 - q
}

// AvailabilityKofN returns the availability of a k-out-of-n system of identical
// units each with availability a, using the same binomial sum as
// KofNReliability.
func AvailabilityKofN(k, n int, a float64) float64 {
	return KofNReliability(k, n, a)
}

// ExpectedDowntime returns the expected downtime over a mission of the given
// duration for an item with the supplied availability: (1-availability)·mission.
func ExpectedDowntime(availability, mission float64) float64 {
	if availability < 0 || availability > 1 || mission < 0 {
		return math.NaN()
	}
	return (1 - availability) * mission
}

// ExpectedUptime returns the expected uptime over a mission of the given
// duration for an item with the supplied availability: availability·mission.
func ExpectedUptime(availability, mission float64) float64 {
	if availability < 0 || availability > 1 || mission < 0 {
		return math.NaN()
	}
	return availability * mission
}

// AvailabilityFromRates returns the steady-state availability from a failure
// rate lambda and repair rate mu, an alias of SteadyStateAvailability.
func AvailabilityFromRates(lambda, mu float64) float64 {
	return SteadyStateAvailability(lambda, mu)
}
