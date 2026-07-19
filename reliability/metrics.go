package reliability

import (
	"errors"
	"math"
	"sort"
)

// FailureRate returns the instantaneous hazard rate f/R from a probability
// density f and a reliability (survival) value r. It returns NaN if r is not
// positive.
func FailureRate(f, r float64) float64 {
	if r <= 0 || f < 0 {
		return math.NaN()
	}
	return f / r
}

// HazardFromReliability converts a reliability value r into its cumulative
// hazard H=-ln(r).
func HazardFromReliability(r float64) float64 {
	if r <= 0 || r > 1 {
		return math.NaN()
	}
	return -math.Log(r)
}

// ReliabilityFromHazard converts a cumulative hazard H>=0 into a reliability
// value R=e^{-H}.
func ReliabilityFromHazard(h float64) float64 {
	if h < 0 {
		return math.NaN()
	}
	return math.Exp(-h)
}

// FailuresPerMillionHours converts a failure rate expressed per unit time into
// the FIT-like measure of expected failures per one million hours.
func FailuresPerMillionHours(lambdaPerHour float64) float64 {
	if lambdaPerHour < 0 {
		return math.NaN()
	}
	return lambdaPerHour * 1e6
}

// FIT returns the failure-in-time rate: the number of failures expected per one
// billion (10^9) device-hours for a constant hourly failure rate.
func FIT(lambdaPerHour float64) float64 {
	if lambdaPerHour < 0 {
		return math.NaN()
	}
	return lambdaPerHour * 1e9
}

// FITToRate converts a failure-in-time value (failures per 10^9 hours) back to
// an hourly failure rate.
func FITToRate(fit float64) float64 {
	if fit < 0 {
		return math.NaN()
	}
	return fit / 1e9
}

// MTBFFromFailureRate returns the mean time between failures 1/λ for a constant
// failure rate lambda>0.
func MTBFFromFailureRate(lambda float64) float64 {
	if lambda <= 0 {
		return math.NaN()
	}
	return 1 / lambda
}

// FailureRateFromMTBF returns the constant failure rate 1/MTBF.
func FailureRateFromMTBF(mtbf float64) float64 {
	if mtbf <= 0 {
		return math.NaN()
	}
	return 1 / mtbf
}

// MTBF returns the mean time between failures MTTF+MTTR of a repairable item.
func MTBF(mttf, mttr float64) float64 {
	if mttf < 0 || mttr < 0 {
		return math.NaN()
	}
	return mttf + mttr
}

// MTBFFromData estimates the mean time between failures as the total operating
// time divided by the number of failures observed.
func MTBFFromData(totalOperatingTime float64, failures int) float64 {
	if totalOperatingTime < 0 || failures <= 0 {
		return math.NaN()
	}
	return totalOperatingTime / float64(failures)
}

// MTTFFromData estimates the mean time to failure as the arithmetic mean of a
// set of (complete, uncensored) failure times. It returns an error if the
// slice is empty or contains a negative time.
func MTTFFromData(times []float64) (float64, error) {
	if len(times) == 0 {
		return 0, errors.New("reliability: MTTFFromData requires at least one time")
	}
	sum := 0.0
	for _, t := range times {
		if t < 0 {
			return 0, errors.New("reliability: failure times must be non-negative")
		}
		sum += t
	}
	return sum / float64(len(times)), nil
}

// MTTFFromReliability computes the mean time to failure by numerically
// integrating a reliability function R over [0,∞). The function must be
// non-increasing and eventually decay to zero; scale is a characteristic time
// used to size the integration panels.
func MTTFFromReliability(r func(float64) float64, scale float64) float64 {
	if scale <= 0 {
		scale = 1
	}
	return integrateToInfinity(r, 0, scale/4)
}

// MTTR returns the mean time to repair 1/μ from a constant repair rate mu>0.
func MTTR(mu float64) float64 {
	if mu <= 0 {
		return math.NaN()
	}
	return 1 / mu
}

// RepairRateFromMTTR returns the repair rate μ=1/MTTR.
func RepairRateFromMTTR(mttr float64) float64 {
	if mttr <= 0 {
		return math.NaN()
	}
	return 1 / mttr
}

// MedianRank returns Bernard's approximation to the median rank of the i-th of
// n ordered failures, (i-0.3)/(n+0.4), commonly used to plot empirical
// probabilities on Weibull paper. The rank i is 1-based.
func MedianRank(i, n int) float64 {
	if i < 1 || n < 1 || i > n {
		return math.NaN()
	}
	return (float64(i) - 0.3) / (float64(n) + 0.4)
}

// MeanRank returns the mean plotting position i/(n+1) of the i-th of n ordered
// failures. The rank i is 1-based.
func MeanRank(i, n int) float64 {
	if i < 1 || n < 1 || i > n {
		return math.NaN()
	}
	return float64(i) / (float64(n) + 1)
}

// EmpiricalReliability returns the fraction of a sample of failure times that
// exceed t, an unbiased estimate of R(t) for complete data.
func EmpiricalReliability(times []float64, t float64) float64 {
	n := len(times)
	if n == 0 {
		return math.NaN()
	}
	survivors := 0
	for _, x := range times {
		if x > t {
			survivors++
		}
	}
	return float64(survivors) / float64(n)
}

// EmpiricalCDF returns the fraction of a sample of failure times that are less
// than or equal to t, an estimate of F(t) for complete data.
func EmpiricalCDF(times []float64, t float64) float64 {
	r := EmpiricalReliability(times, t)
	if math.IsNaN(r) {
		return math.NaN()
	}
	return 1 - r
}

// EmpiricalHazardRate estimates the average hazard rate over the interval
// (t, t+dt] from complete data as (survivors past t that fail within dt)
// divided by (dt times survivors past t).
func EmpiricalHazardRate(times []float64, t, dt float64) float64 {
	if dt <= 0 {
		return math.NaN()
	}
	nAtT := 0
	failWithin := 0
	for _, x := range times {
		if x > t {
			nAtT++
			if x <= t+dt {
				failWithin++
			}
		}
	}
	if nAtT == 0 {
		return math.NaN()
	}
	return float64(failWithin) / (dt * float64(nAtT))
}

// SortedFailureTimes returns a sorted copy of a slice of failure times.
func SortedFailureTimes(times []float64) []float64 {
	out := make([]float64, len(times))
	copy(out, times)
	sort.Float64s(out)
	return out
}

// ReliabilityAtTime returns the reliability of a constant-failure-rate item at
// time t: R(t)=e^{-λt}. It is a convenience alias of ExponentialReliability.
func ReliabilityAtTime(t, lambda float64) float64 {
	return ExponentialReliability(t, lambda)
}

// DesignLifeForReliability returns the maximum operating time at which a
// constant-failure-rate item still meets a target reliability rTarget:
// t=-ln(rTarget)/λ.
func DesignLifeForReliability(rTarget, lambda float64) float64 {
	if lambda <= 0 || rTarget <= 0 || rTarget > 1 {
		return math.NaN()
	}
	return -math.Log(rTarget) / lambda
}
