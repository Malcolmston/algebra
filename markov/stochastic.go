package markov

import "math"

// DefaultTol is the tolerance used by predicates that accept none explicitly
// elsewhere and by the convenience constructors.
const DefaultTol = 1e-9

// RowSums returns the sum of each row of a.
func RowSums(a [][]float64) []float64 {
	s := make([]float64, len(a))
	for i := range a {
		for _, x := range a[i] {
			s[i] += x
		}
	}
	return s
}

// ColSums returns the sum of each column of a.
func ColSums(a [][]float64) []float64 {
	if len(a) == 0 {
		return nil
	}
	c := len(a[0])
	s := make([]float64, c)
	for i := range a {
		for j := 0; j < c && j < len(a[i]); j++ {
			s[j] += a[i][j]
		}
	}
	return s
}

// IsNonNegative reports whether every entry of a is >= -tol.
func IsNonNegative(a [][]float64, tol float64) bool {
	for i := range a {
		for _, x := range a[i] {
			if x < -tol {
				return false
			}
		}
	}
	return true
}

// IsStochastic reports whether a is (row) stochastic: square, non-negative, and
// every row sums to 1 within tol.
func IsStochastic(a [][]float64, tol float64) bool {
	if !IsSquare(a) || !IsNonNegative(a, tol) {
		return false
	}
	for _, s := range RowSums(a) {
		if math.Abs(s-1) > tol {
			return false
		}
	}
	return true
}

// IsColumnStochastic reports whether a is square, non-negative, and every
// column sums to 1 within tol.
func IsColumnStochastic(a [][]float64, tol float64) bool {
	if !IsSquare(a) || !IsNonNegative(a, tol) {
		return false
	}
	for _, s := range ColSums(a) {
		if math.Abs(s-1) > tol {
			return false
		}
	}
	return true
}

// IsDoublyStochastic reports whether a is both row- and column-stochastic
// within tol.
func IsDoublyStochastic(a [][]float64, tol float64) bool {
	return IsStochastic(a, tol) && IsColumnStochastic(a, tol)
}

// IsSubStochastic reports whether a is non-negative and every row sums to at
// most 1+tol (rows may sum to less than 1). This is the form of the Q block of
// an absorbing chain.
func IsSubStochastic(a [][]float64, tol float64) bool {
	if !IsNonNegative(a, tol) {
		return false
	}
	for _, s := range RowSums(a) {
		if s > 1+tol {
			return false
		}
	}
	return true
}

// NormalizeRows returns a copy of a with each non-zero row scaled to sum to 1.
// Rows that sum to zero are left unchanged.
func NormalizeRows(a [][]float64) [][]float64 {
	b := CopyMatrix(a)
	for i := range b {
		var s float64
		for _, x := range b[i] {
			s += x
		}
		if s != 0 {
			for j := range b[i] {
				b[i][j] /= s
			}
		}
	}
	return b
}

// Normalize returns a copy of v scaled to sum to 1. If v sums to zero it is
// returned unchanged (copied).
func Normalize(v []float64) []float64 {
	w := CopyVector(v)
	var s float64
	for _, x := range w {
		s += x
	}
	if s != 0 {
		for i := range w {
			w[i] /= s
		}
	}
	return w
}

// UniformDistribution returns the uniform probability vector of length n.
func UniformDistribution(n int) []float64 {
	if n <= 0 {
		return nil
	}
	v := make([]float64, n)
	p := 1.0 / float64(n)
	for i := range v {
		v[i] = p
	}
	return v
}

// PointMass returns the length-n probability vector placing all mass on state i.
func PointMass(n, i int) []float64 {
	if n <= 0 || i < 0 || i >= n {
		return nil
	}
	v := make([]float64, n)
	v[i] = 1
	return v
}

// IsProbabilityVector reports whether v is non-negative and sums to 1 within tol.
func IsProbabilityVector(v []float64, tol float64) bool {
	var s float64
	for _, x := range v {
		if x < -tol {
			return false
		}
		s += x
	}
	return math.Abs(s-1) <= tol
}

// TotalVariationDistance returns the total-variation distance between two
// probability vectors p and q, defined as half their L1 distance. It returns
// NaN if the lengths differ.
func TotalVariationDistance(p, q []float64) float64 {
	if len(p) != len(q) {
		return math.NaN()
	}
	var s float64
	for i := range p {
		s += math.Abs(p[i] - q[i])
	}
	return 0.5 * s
}

// L1Distance returns the L1 (taxicab) distance between p and q. It returns NaN
// if the lengths differ.
func L1Distance(p, q []float64) float64 {
	if len(p) != len(q) {
		return math.NaN()
	}
	var s float64
	for i := range p {
		s += math.Abs(p[i] - q[i])
	}
	return s
}

// EuclideanDistance returns the L2 distance between p and q. It returns NaN if
// the lengths differ.
func EuclideanDistance(p, q []float64) float64 {
	if len(p) != len(q) {
		return math.NaN()
	}
	var s float64
	for i := range p {
		d := p[i] - q[i]
		s += d * d
	}
	return math.Sqrt(s)
}

// HellingerDistance returns the Hellinger distance between probability vectors
// p and q, in [0,1]. It returns NaN if the lengths differ.
func HellingerDistance(p, q []float64) float64 {
	if len(p) != len(q) {
		return math.NaN()
	}
	var s float64
	for i := range p {
		d := math.Sqrt(math.Max(p[i], 0)) - math.Sqrt(math.Max(q[i], 0))
		s += d * d
	}
	return math.Sqrt(s) / math.Sqrt2
}

// KLDivergence returns the Kullback-Leibler divergence D(p‖q) in nats. Terms
// with p[i]=0 contribute 0. If some q[i]=0 while p[i]>0 the divergence is +Inf.
// It returns NaN if the lengths differ.
func KLDivergence(p, q []float64) float64 {
	if len(p) != len(q) {
		return math.NaN()
	}
	var s float64
	for i := range p {
		if p[i] > 0 {
			if q[i] <= 0 {
				return math.Inf(1)
			}
			s += p[i] * math.Log(p[i]/q[i])
		}
	}
	return s
}

// JensenShannonDivergence returns the Jensen-Shannon divergence between p and q
// in nats: the symmetrized, smoothed average of the two KL divergences to their
// mixture. It returns NaN if the lengths differ.
func JensenShannonDivergence(p, q []float64) float64 {
	if len(p) != len(q) {
		return math.NaN()
	}
	m := make([]float64, len(p))
	for i := range p {
		m[i] = 0.5 * (p[i] + q[i])
	}
	return 0.5*KLDivergence(p, m) + 0.5*KLDivergence(q, m)
}

// ShannonEntropy returns the Shannon entropy of the probability vector p in
// nats. Zero-probability outcomes contribute 0.
func ShannonEntropy(p []float64) float64 {
	var s float64
	for _, x := range p {
		if x > 0 {
			s -= x * math.Log(x)
		}
	}
	return s
}
