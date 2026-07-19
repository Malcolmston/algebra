package markov

import "math"

// HadamardProduct returns the elementwise product of a and b. It returns nil if
// the shapes differ.
func HadamardProduct(a, b [][]float64) [][]float64 {
	if len(a) != len(b) {
		return nil
	}
	c := make([][]float64, len(a))
	for i := range a {
		if len(a[i]) != len(b[i]) {
			return nil
		}
		c[i] = make([]float64, len(a[i]))
		for j := range a[i] {
			c[i][j] = a[i][j] * b[i][j]
		}
	}
	return c
}

// OuterProduct returns the outer product u·vᵀ as a len(u)×len(v) matrix.
func OuterProduct(u, v []float64) [][]float64 {
	m := make([][]float64, len(u))
	for i := range u {
		m[i] = make([]float64, len(v))
		for j := range v {
			m[i][j] = u[i] * v[j]
		}
	}
	return m
}

// KroneckerProduct returns the Kronecker product a⊗b. If a is p×q and b is r×s
// the result is (p·r)×(q·s). It is used to build the transition matrix of the
// product of independent chains.
func KroneckerProduct(a, b [][]float64) [][]float64 {
	if len(a) == 0 || len(b) == 0 {
		return nil
	}
	ar, ac := len(a), len(a[0])
	br, bc := len(b), len(b[0])
	out := make([][]float64, ar*br)
	for i := range out {
		out[i] = make([]float64, ac*bc)
	}
	for i := 0; i < ar; i++ {
		for j := 0; j < ac; j++ {
			for k := 0; k < br; k++ {
				for l := 0; l < bc; l++ {
					out[i*br+k][j*bc+l] = a[i][j] * b[k][l]
				}
			}
		}
	}
	return out
}

// IsSymmetricMatrix reports whether the square matrix a equals its transpose
// within tol.
func IsSymmetricMatrix(a [][]float64, tol float64) bool {
	if !IsSquare(a) {
		return false
	}
	for i := range a {
		for j := i + 1; j < len(a); j++ {
			if math.Abs(a[i][j]-a[j][i]) > tol {
				return false
			}
		}
	}
	return true
}

// ArgMax returns the index of the (first) maximal entry of v, or -1 if v is
// empty.
func ArgMax(v []float64) int {
	if len(v) == 0 {
		return -1
	}
	best := 0
	for i := 1; i < len(v); i++ {
		if v[i] > v[best] {
			best = i
		}
	}
	return best
}

// ArgMin returns the index of the (first) minimal entry of v, or -1 if v is
// empty.
func ArgMin(v []float64) int {
	if len(v) == 0 {
		return -1
	}
	best := 0
	for i := 1; i < len(v); i++ {
		if v[i] < v[best] {
			best = i
		}
	}
	return best
}

// MaxValue returns the maximum entry of v, or NaN if v is empty.
func MaxValue(v []float64) float64 {
	if len(v) == 0 {
		return math.NaN()
	}
	m := v[0]
	for _, x := range v[1:] {
		if x > m {
			m = x
		}
	}
	return m
}

// MinValue returns the minimum entry of v, or NaN if v is empty.
func MinValue(v []float64) float64 {
	if len(v) == 0 {
		return math.NaN()
	}
	m := v[0]
	for _, x := range v[1:] {
		if x < m {
			m = x
		}
	}
	return m
}

// LogSumExp returns log(Σ exp(v_i)) computed stably by factoring out the maximum.
func LogSumExp(v []float64) float64 {
	if len(v) == 0 {
		return math.Inf(-1)
	}
	m := MaxValue(v)
	if math.IsInf(m, -1) {
		return math.Inf(-1)
	}
	var s float64
	for _, x := range v {
		s += math.Exp(x - m)
	}
	return m + math.Log(s)
}

// Softmax returns the softmax (normalized exponentials) of v, computed stably.
func Softmax(v []float64) []float64 {
	if len(v) == 0 {
		return nil
	}
	m := MaxValue(v)
	out := make([]float64, len(v))
	var s float64
	for i, x := range v {
		out[i] = math.Exp(x - m)
		s += out[i]
	}
	if s > 0 {
		for i := range out {
			out[i] /= s
		}
	}
	return out
}

// BhattacharyyaCoefficient returns Σ_i sqrt(p_i q_i), a similarity measure in
// [0,1] for probability vectors. It returns NaN if the lengths differ.
func BhattacharyyaCoefficient(p, q []float64) float64 {
	if len(p) != len(q) {
		return math.NaN()
	}
	var s float64
	for i := range p {
		s += math.Sqrt(math.Max(p[i], 0) * math.Max(q[i], 0))
	}
	return s
}

// BhattacharyyaDistance returns -log of the Bhattacharyya coefficient. It is 0
// for identical distributions and grows as they diverge.
func BhattacharyyaDistance(p, q []float64) float64 {
	bc := BhattacharyyaCoefficient(p, q)
	if bc <= 0 {
		return math.Inf(1)
	}
	return -math.Log(bc)
}

// ChiSquareDistance returns the Pearson chi-square distance Σ_i (p_i-q_i)²/q_i
// between probability vectors. Terms with q_i=0 and p_i=0 contribute 0; a term
// with q_i=0 but p_i>0 yields +Inf. It returns NaN if the lengths differ.
func ChiSquareDistance(p, q []float64) float64 {
	if len(p) != len(q) {
		return math.NaN()
	}
	var s float64
	for i := range p {
		d := p[i] - q[i]
		if q[i] == 0 {
			if d != 0 {
				return math.Inf(1)
			}
			continue
		}
		s += d * d / q[i]
	}
	return s
}

// CrossEntropy returns the cross entropy H(p, q) = -Σ_i p_i log q_i in nats.
// Terms with p_i=0 contribute 0; if q_i=0 while p_i>0 the result is +Inf. It
// returns NaN if the lengths differ.
func CrossEntropy(p, q []float64) float64 {
	if len(p) != len(q) {
		return math.NaN()
	}
	var s float64
	for i := range p {
		if p[i] > 0 {
			if q[i] <= 0 {
				return math.Inf(1)
			}
			s -= p[i] * math.Log(q[i])
		}
	}
	return s
}

// MeanVector returns the componentwise mean of a multidimensional chain (a
// slice of equal-length state vectors).
func MeanVector(chain [][]float64) []float64 {
	if len(chain) == 0 {
		return nil
	}
	d := len(chain[0])
	out := make([]float64, d)
	for _, row := range chain {
		for k := 0; k < d && k < len(row); k++ {
			out[k] += row[k]
		}
	}
	inv := 1 / float64(len(chain))
	for k := range out {
		out[k] *= inv
	}
	return out
}

// CovarianceMatrix returns the unbiased sample covariance matrix of a
// multidimensional chain (rows are observations, columns are variables). It
// returns nil if there are fewer than two observations.
func CovarianceMatrix(chain [][]float64) [][]float64 {
	n := len(chain)
	if n < 2 {
		return nil
	}
	d := len(chain[0])
	mean := MeanVector(chain)
	cov := make([][]float64, d)
	for i := range cov {
		cov[i] = make([]float64, d)
	}
	for _, row := range chain {
		for i := 0; i < d; i++ {
			for j := 0; j < d; j++ {
				cov[i][j] += (row[i] - mean[i]) * (row[j] - mean[j])
			}
		}
	}
	inv := 1 / float64(n-1)
	for i := range cov {
		for j := range cov[i] {
			cov[i][j] *= inv
		}
	}
	return cov
}
