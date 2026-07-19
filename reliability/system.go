package reliability

import (
	"errors"
	"math"
)

// SeriesReliability returns the reliability of a series system, the product of
// the component reliabilities. A series system works only if every component
// works. With no arguments it returns 1 (the empty product).
func SeriesReliability(rs ...float64) float64 {
	prod := 1.0
	for _, r := range rs {
		if r < 0 || r > 1 {
			return math.NaN()
		}
		prod *= r
	}
	return prod
}

// ParallelReliability returns the reliability of a parallel (active-redundant)
// system, 1-∏(1-rᵢ). A parallel system fails only if every component fails.
// With no arguments it returns 0 (an empty system never works).
func ParallelReliability(rs ...float64) float64 {
	q := 1.0
	for _, r := range rs {
		if r < 0 || r > 1 {
			return math.NaN()
		}
		q *= 1 - r
	}
	return 1 - q
}

// SeriesFailureProbability returns the probability that a series system fails,
// 1 minus its reliability.
func SeriesFailureProbability(rs ...float64) float64 {
	r := SeriesReliability(rs...)
	if math.IsNaN(r) {
		return math.NaN()
	}
	return 1 - r
}

// ParallelFailureProbability returns the probability that a parallel system
// fails, the product ∏(1-rᵢ).
func ParallelFailureProbability(rs ...float64) float64 {
	q := 1.0
	for _, r := range rs {
		if r < 0 || r > 1 {
			return math.NaN()
		}
		q *= 1 - r
	}
	return q
}

// SeriesReliabilityIdentical returns the reliability rⁿ of a series system of n
// identical components each with reliability r.
func SeriesReliabilityIdentical(r float64, n int) float64 {
	if r < 0 || r > 1 || n < 0 {
		return math.NaN()
	}
	return math.Pow(r, float64(n))
}

// ParallelReliabilityIdentical returns the reliability 1-(1-r)ⁿ of a parallel
// system of n identical components each with reliability r.
func ParallelReliabilityIdentical(r float64, n int) float64 {
	if r < 0 || r > 1 || n < 0 {
		return math.NaN()
	}
	return 1 - math.Pow(1-r, float64(n))
}

// SeriesParallelReliability returns the reliability of a system built from
// parallel groups connected in series: each inner slice is a parallel
// sub-system and the sub-systems are combined in series.
func SeriesParallelReliability(groups [][]float64) float64 {
	prod := 1.0
	for _, g := range groups {
		p := ParallelReliability(g...)
		if math.IsNaN(p) {
			return math.NaN()
		}
		prod *= p
	}
	return prod
}

// ParallelSeriesReliability returns the reliability of a system built from
// series strings connected in parallel: each inner slice is a series string
// and the strings are combined in parallel.
func ParallelSeriesReliability(strings [][]float64) float64 {
	q := 1.0
	for _, s := range strings {
		r := SeriesReliability(s...)
		if math.IsNaN(r) {
			return math.NaN()
		}
		q *= 1 - r
	}
	return 1 - q
}

// KofNReliability returns the reliability of a k-out-of-n system of identical
// components, each with reliability r, that works when at least k of its n
// components work: Σ_{i=k}^{n} C(n,i) rⁱ(1-r)^{n-i}.
func KofNReliability(k, n int, r float64) float64 {
	if r < 0 || r > 1 || k < 0 || n < 0 || k > n {
		return math.NaN()
	}
	sum := 0.0
	for i := k; i <= n; i++ {
		sum += binomial(n, i) * math.Pow(r, float64(i)) * math.Pow(1-r, float64(n-i))
	}
	return sum
}

// KofNReliabilityGeneral returns the exact reliability of a k-out-of-n system
// whose n components may have different reliabilities rs. It works when at
// least k components work and is evaluated by a dynamic program over the
// probability of exactly j working components.
func KofNReliabilityGeneral(k int, rs []float64) (float64, error) {
	n := len(rs)
	if k < 0 || k > n {
		return 0, errors.New("reliability: KofNReliabilityGeneral requires 0<=k<=len(rs)")
	}
	for _, r := range rs {
		if r < 0 || r > 1 {
			return 0, errors.New("reliability: component reliabilities must lie in [0,1]")
		}
	}
	// dp[j] = probability that exactly j of the processed components work.
	dp := make([]float64, n+1)
	dp[0] = 1
	processed := 0
	for _, r := range rs {
		for j := processed + 1; j >= 1; j-- {
			dp[j] = dp[j]*(1-r) + dp[j-1]*r
		}
		dp[0] *= 1 - r
		processed++
	}
	sum := 0.0
	for j := k; j <= n; j++ {
		sum += dp[j]
	}
	return sum, nil
}

// BridgeReliability returns the exact reliability of the classic five-component
// bridge network, evaluated by conditioning (pivotal decomposition) on the
// bridging component r5. Components r1 and r2 form the two upper links, r3 and
// r4 the two lower links, and r5 the bridge between the mid-nodes. Sources and
// terminals are assumed perfectly reliable.
func BridgeReliability(r1, r2, r3, r4, r5 float64) float64 {
	for _, r := range []float64{r1, r2, r3, r4, r5} {
		if r < 0 || r > 1 {
			return math.NaN()
		}
	}
	// Condition on the bridge component r5.
	// If r5 works, nodes are merged: (r1 ∥ r3) in series with (r2 ∥ r4).
	up := ParallelReliability(r1, r3)
	down := ParallelReliability(r2, r4)
	withBridge := up * down
	// If r5 fails, two independent series paths in parallel:
	// (r1·r2) ∥ (r3·r4).
	withoutBridge := ParallelReliability(r1*r2, r3*r4)
	return r5*withBridge + (1-r5)*withoutBridge
}

// StandbyRedundancyReliability returns the reliability at time t of a system of
// n identical exponential units (one active, n-1 in cold standby with perfect
// switching) with per-unit failure rate lambda. The number of failures before
// system failure follows a Poisson process, giving
// R(t)=e^{-λt}Σ_{i=0}^{n-1}(λt)^i/i!.
func StandbyRedundancyReliability(t, lambda float64, n int) float64 {
	if t < 0 || lambda <= 0 || n < 1 {
		return math.NaN()
	}
	x := lambda * t
	sum := 0.0
	term := 1.0
	for i := 0; i < n; i++ {
		if i > 0 {
			term *= x / float64(i)
		}
		sum += term
	}
	return math.Exp(-x) * sum
}

// StandbyRedundancyMTTF returns the mean time to failure n/λ of a cold-standby
// system of n identical exponential units with per-unit failure rate lambda
// and perfect switching.
func StandbyRedundancyMTTF(lambda float64, n int) float64 {
	if lambda <= 0 || n < 1 {
		return math.NaN()
	}
	return float64(n) / lambda
}

// SeriesSystemFailureRate returns the total failure rate Σλᵢ of a series system
// of independent constant-rate components.
func SeriesSystemFailureRate(lambdas ...float64) float64 {
	sum := 0.0
	for _, l := range lambdas {
		if l < 0 {
			return math.NaN()
		}
		sum += l
	}
	return sum
}

// SeriesSystemMTTF returns the mean time to failure 1/Σλᵢ of a series system of
// independent constant-rate components.
func SeriesSystemMTTF(lambdas ...float64) float64 {
	total := SeriesSystemFailureRate(lambdas...)
	if math.IsNaN(total) || total <= 0 {
		return math.NaN()
	}
	return 1 / total
}

// ParallelSystemMTTF returns the mean time to failure of a parallel system of n
// identical exponential units with failure rate lambda, which is the harmonic
// sum (1/λ)Σ_{i=1}^{n}1/i.
func ParallelSystemMTTF(lambda float64, n int) float64 {
	if lambda <= 0 || n < 1 {
		return math.NaN()
	}
	h := 0.0
	for i := 1; i <= n; i++ {
		h += 1 / float64(i)
	}
	return h / lambda
}

// binomial returns the binomial coefficient C(n,k) as a float64, computed
// multiplicatively to limit overflow and rounding error.
func binomial(n, k int) float64 {
	if k < 0 || k > n {
		return 0
	}
	if k > n-k {
		k = n - k
	}
	res := 1.0
	for i := 0; i < k; i++ {
		res = res * float64(n-i) / float64(i+1)
	}
	return math.Round(res)
}
