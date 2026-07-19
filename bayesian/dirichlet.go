package bayesian

import (
	"math"
)

// Dirichlet is a Dirichlet distribution over the probability simplex with
// concentration vector Alpha, all entries strictly positive. It is the
// conjugate prior for the category probabilities of a Categorical/Multinomial
// likelihood.
type Dirichlet struct {
	Alpha []float64
}

// NewDirichlet constructs a Dirichlet distribution from a concentration vector,
// returning ErrParam if it is empty or has a non-positive entry. The slice is
// copied.
func NewDirichlet(alpha []float64) (Dirichlet, error) {
	if len(alpha) == 0 {
		return Dirichlet{}, ErrParam
	}
	cp := make([]float64, len(alpha))
	for i, a := range alpha {
		if a <= 0 {
			return Dirichlet{}, ErrParam
		}
		cp[i] = a
	}
	return Dirichlet{Alpha: cp}, nil
}

// Dim returns the number of categories.
func (d Dirichlet) Dim() int { return len(d.Alpha) }

// Sum returns the total concentration α₀ = Σαᵢ.
func (d Dirichlet) Sum() float64 {
	var s float64
	for _, a := range d.Alpha {
		s += a
	}
	return s
}

// Mean returns the mean vector αᵢ/α₀.
func (d Dirichlet) Mean() []float64 {
	s := d.Sum()
	out := make([]float64, len(d.Alpha))
	for i, a := range d.Alpha {
		out[i] = a / s
	}
	return out
}

// Variance returns the per-component marginal variances.
func (d Dirichlet) Variance() []float64 {
	s := d.Sum()
	out := make([]float64, len(d.Alpha))
	for i, a := range d.Alpha {
		m := a / s
		out[i] = m * (1 - m) / (s + 1)
	}
	return out
}

// Mode returns the mode vector (αᵢ−1)/(α₀−K); it is defined only when every
// αᵢ > 1, otherwise the returned slice is nil.
func (d Dirichlet) Mode() []float64 {
	k := float64(len(d.Alpha))
	den := d.Sum() - k
	if den <= 0 {
		return nil
	}
	out := make([]float64, len(d.Alpha))
	for i, a := range d.Alpha {
		if a <= 1 {
			return nil
		}
		out[i] = (a - 1) / den
	}
	return out
}

// Covariance returns the K×K covariance matrix of the distribution as a slice
// of rows.
func (d Dirichlet) Covariance() [][]float64 {
	s := d.Sum()
	k := len(d.Alpha)
	cov := make([][]float64, k)
	for i := range cov {
		cov[i] = make([]float64, k)
	}
	for i := 0; i < k; i++ {
		for j := 0; j < k; j++ {
			if i == j {
				ai := d.Alpha[i]
				cov[i][j] = ai * (s - ai) / (s * s * (s + 1))
			} else {
				cov[i][j] = -d.Alpha[i] * d.Alpha[j] / (s * s * (s + 1))
			}
		}
	}
	return cov
}

// MeanLog returns the vector E[ln Xᵢ] = ψ(αᵢ) − ψ(α₀).
func (d Dirichlet) MeanLog() []float64 {
	s0 := Digamma(d.Sum())
	out := make([]float64, len(d.Alpha))
	for i, a := range d.Alpha {
		out[i] = Digamma(a) - s0
	}
	return out
}

// Marginal returns the marginal distribution of component i, which is
// Beta(αᵢ, α₀−αᵢ).
func (d Dirichlet) Marginal(i int) Beta {
	ai := d.Alpha[i]
	return Beta{Alpha: ai, Beta: d.Sum() - ai}
}

// LogPDF returns the natural logarithm of the density at the point x on the
// simplex. It returns math.Inf(-1) when x has the wrong length, a non-positive
// coordinate, or does not sum to one within a small tolerance.
func (d Dirichlet) LogPDF(x []float64) float64 {
	if len(x) != len(d.Alpha) {
		return math.Inf(-1)
	}
	var sum, ll float64
	logNorm := 0.0
	for i, a := range d.Alpha {
		if x[i] <= 0 || x[i] >= 1 {
			return math.Inf(-1)
		}
		sum += x[i]
		ll += (a - 1) * math.Log(x[i])
		logNorm += LogGamma(a)
	}
	if math.Abs(sum-1) > 1e-9 {
		return math.Inf(-1)
	}
	logNorm -= LogGamma(d.Sum())
	return ll - logNorm
}

// PDF returns the density at the point x on the simplex.
func (d Dirichlet) PDF(x []float64) float64 {
	return math.Exp(d.LogPDF(x))
}

// Entropy returns the differential entropy of the distribution in nats.
func (d Dirichlet) Entropy() float64 {
	s := d.Sum()
	k := float64(len(d.Alpha))
	// log B(alpha)
	logB := -LogGamma(s)
	var term float64
	for _, a := range d.Alpha {
		logB += LogGamma(a)
		term += (a - 1) * Digamma(a)
	}
	return logB + (s-k)*Digamma(s) - term
}
