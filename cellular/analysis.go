package cellular

import "math"

// BinaryEntropy returns the Shannon entropy in bits of a Bernoulli(p)
// distribution, -p*log2(p) - (1-p)*log2(1-p). It returns 0 at p = 0 or p = 1.
func BinaryEntropy(p float64) float64 {
	if p <= 0 || p >= 1 {
		return 0
	}
	return -p*math.Log2(p) - (1-p)*math.Log2(1-p)
}

// ShannonEntropyBits returns the Shannon entropy in bits of a discrete
// probability distribution. Negative entries are treated as 0 and the entries
// need not sum to exactly 1 (they are used as given, with 0*log0 taken to be 0).
func ShannonEntropyBits(probs []float64) float64 {
	h := 0.0
	for _, p := range probs {
		if p > 0 {
			h -= p * math.Log2(p)
		}
	}
	return h
}

// SymbolEntropy returns the empirical Shannon entropy in bits of the symbol
// distribution of state over the alphabet 0..k-1.
func SymbolEntropy(state []int, k int) float64 {
	if len(state) == 0 || k < 1 {
		return 0
	}
	counts := make([]float64, k)
	for _, v := range state {
		if v >= 0 && v < k {
			counts[v]++
		}
	}
	n := float64(len(state))
	for i := range counts {
		counts[i] /= n
	}
	return ShannonEntropyBits(counts)
}

// SpatialEntropy returns the binary symbol entropy of state (equivalent to
// SymbolEntropy with k = 2 for a binary configuration).
func SpatialEntropy(state []int) float64 {
	return SymbolEntropy(state, 2)
}

// BlockEntropy returns the Shannon entropy in bits of the distribution of length
// blockLen blocks read with wraparound across state over the alphabet 0..k-1. It
// returns 0 for non-positive blockLen or an empty state.
func BlockEntropy(state []int, blockLen, k int) float64 {
	n := len(state)
	if n == 0 || blockLen <= 0 || k < 2 {
		return 0
	}
	counts := map[int]int{}
	for i := 0; i < n; i++ {
		code := 0
		for j := 0; j < blockLen; j++ {
			code = code*k + state[(i+j)%n]
		}
		counts[code]++
	}
	probs := make([]float64, 0, len(counts))
	total := float64(n)
	for _, c := range counts {
		probs = append(probs, float64(c)/total)
	}
	return ShannonEntropyBits(probs)
}

// SpacetimeDensity returns the mean fraction of non-zero cells across every row
// of a spacetime diagram.
func SpacetimeDensity(rows [][]int) float64 {
	if len(rows) == 0 {
		return 0
	}
	sum := 0.0
	for _, row := range rows {
		sum += Density(row)
	}
	return sum / float64(len(rows))
}

// MeanSpatialEntropy returns the average SpatialEntropy over the rows of a
// spacetime diagram.
func MeanSpatialEntropy(rows [][]int) float64 {
	if len(rows) == 0 {
		return 0
	}
	sum := 0.0
	for _, row := range rows {
		sum += SpatialEntropy(row)
	}
	return sum / float64(len(rows))
}

// DamageSpread evolves two configurations that differ by a single flipped cell
// and returns, for each of the steps time steps, the Hamming distance between
// the two evolutions. It is a proxy for sensitivity to initial conditions
// (Lyapunov-like behaviour). The reference initial condition is used unchanged
// and a copy with the cell at index flipAt toggled (mod k) forms the perturbed
// run.
func DamageSpread(rule Rule1D, initial []int, flipAt, steps int, bc Boundary) []int {
	k := rule.States()
	a := CloneState(initial)
	b := CloneState(initial)
	if flipAt >= 0 && flipAt < len(b) {
		b[flipAt] = (b[flipAt] + 1) % k
	}
	out := make([]int, steps+1)
	d, _ := HammingDistance(a, b)
	out[0] = d
	for t := 1; t <= steps; t++ {
		a = Step1D(rule, a, bc)
		b = Step1D(rule, b, bc)
		d, _ = HammingDistance(a, b)
		out[t] = d
	}
	return out
}

// TemporalPeriod scans a spacetime diagram from the last row backwards and
// returns the smallest positive p (up to maxPeriod) such that the final row
// equals the row p steps earlier, or 0 if no such period is found. A return of 1
// indicates a fixed point.
func TemporalPeriod(rows [][]int, maxPeriod int) int {
	n := len(rows)
	if n == 0 {
		return 0
	}
	last := rows[n-1]
	for p := 1; p <= maxPeriod && p < n; p++ {
		if EqualState(last, rows[n-1-p]) {
			return p
		}
	}
	return 0
}

// Class returns a heuristic Wolfram class (1, 2 or 3) for an elementary rule,
// obtained by evolving a fixed reproducible pseudo-random configuration on a
// prime-width ring and observing the asymptotic behaviour:
//
//   - Class 1: the configuration collapses to a spatially uniform state.
//   - Class 2: it settles into a short-period cycle (period at most 32).
//   - Class 3: it remains aperiodic within the observation window (chaotic).
//
// The genuinely complex rules that Wolfram places in class 4 (such as rule 110)
// are not separated from class 3 by this simple heuristic and are reported as
// class 3; use it as an indicative, not definitive, classifier.
func (r ElementaryRule) Class() int {
	const width = 151 // prime, avoids the special periods of power-of-two rings
	const steps = 600
	init := RandomBinaryState(width, 0xC0FFEE)
	rows := Evolve1D(r, init, steps, Periodic)
	last := rows[steps]
	if Uniform(last) {
		return 1
	}
	if TemporalPeriod(rows, 32) > 0 {
		return 2
	}
	return 3
}

// ClassName returns the textual name of a Wolfram class number.
func ClassName(class int) string {
	switch class {
	case 1:
		return "uniform"
	case 2:
		return "periodic"
	case 3:
		return "chaotic"
	case 4:
		return "complex"
	default:
		return "unknown"
	}
}
