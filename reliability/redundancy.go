package reliability

import "math"

// MinParallelForReliability returns the smallest number of identical parallel
// components, each with reliability r in (0,1), needed to reach a system
// reliability of at least target: n=⌈ln(1-target)/ln(1-r)⌉.
func MinParallelForReliability(target, r float64) int {
	if r <= 0 || r >= 1 || target <= 0 || target >= 1 {
		return -1
	}
	n := math.Log(1-target) / math.Log(1-r)
	return int(math.Ceil(n - 1e-12))
}

// MaxSeriesForReliability returns the largest number of identical series
// components, each with reliability r in (0,1], whose product reliability still
// meets or exceeds target: n=⌊ln(target)/ln(r)⌋.
func MaxSeriesForReliability(target, r float64) int {
	if r <= 0 || r > 1 || target <= 0 || target > 1 {
		return -1
	}
	if r == 1 {
		return math.MaxInt32
	}
	n := math.Log(target) / math.Log(r)
	return int(math.Floor(n + 1e-12))
}

// RedundancyGain returns the improvement in reliability obtained by placing a
// single spare component of reliability r in active parallel with an existing
// component of reliability r: the difference (1-(1-r)²)-r = r(1-r).
func RedundancyGain(r float64) float64 {
	if r < 0 || r > 1 {
		return math.NaN()
	}
	return r * (1 - r)
}

// ReliabilityImportanceBirnbaum returns the Birnbaum reliability importance of
// a component in a system, the partial derivative of system reliability with
// respect to that component's reliability. It is estimated from the system
// reliability with the component forced to work (rWorking) and forced to fail
// (rFailed): I_B = R(working) - R(failed).
func ReliabilityImportanceBirnbaum(rWorking, rFailed float64) float64 {
	if rWorking < 0 || rWorking > 1 || rFailed < 0 || rFailed > 1 {
		return math.NaN()
	}
	return rWorking - rFailed
}

// CriticalityImportance returns the criticality (Lambert) importance of a
// component: its Birnbaum importance weighted by its own unreliability and
// normalized by the overall system unreliability. rComp is the component
// reliability and rSys the system reliability.
func CriticalityImportance(birnbaum, rComp, rSys float64) float64 {
	sysUnrel := 1 - rSys
	if sysUnrel <= 0 {
		return math.NaN()
	}
	return birnbaum * (1 - rComp) / sysUnrel
}

// AllocateReliabilitySeriesEqual applies the equal-apportionment technique to a
// series system: to achieve an overall target reliability across n independent
// subsystems, each subsystem must reach target^{1/n}.
func AllocateReliabilitySeriesEqual(target float64, n int) float64 {
	if target <= 0 || target > 1 || n < 1 {
		return math.NaN()
	}
	return math.Pow(target, 1/float64(n))
}

// AllocateFailureRateSeriesEqual apportions an overall series-system failure
// rate equally across n subsystems, returning the per-subsystem rate
// lambdaTotal/n.
func AllocateFailureRateSeriesEqual(lambdaTotal float64, n int) float64 {
	if lambdaTotal < 0 || n < 1 {
		return math.NaN()
	}
	return lambdaTotal / float64(n)
}
