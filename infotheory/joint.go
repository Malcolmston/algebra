package infotheory

import "math"

// MarginalX returns the marginal distribution of the row variable X of the
// joint distribution given as a row-major matrix, where joint[i][j] is the
// probability P(X=i, Y=j). The result has one entry per row, each the sum of
// that row.
func MarginalX(joint [][]float64) []float64 {
	px := make([]float64, len(joint))
	for i, row := range joint {
		var s float64
		for _, v := range row {
			s += v
		}
		px[i] = s
	}
	return px
}

// MarginalY returns the marginal distribution of the column variable Y of the
// joint distribution joint, where joint[i][j] is P(X=i, Y=j). The result has
// one entry per column, each the sum of that column. It returns nil for an
// empty joint distribution.
func MarginalY(joint [][]float64) []float64 {
	if len(joint) == 0 {
		return nil
	}
	py := make([]float64, len(joint[0]))
	for _, row := range joint {
		for j, v := range row {
			py[j] += v
		}
	}
	return py
}

// JointEntropy returns the joint Shannon entropy H(X,Y) = -sum_ij p_ij log2
// p_ij of the joint distribution joint, measured in bits.
func JointEntropy(joint [][]float64) float64 {
	var h float64
	for _, row := range joint {
		for _, v := range row {
			h -= infotheoryXLogX(v)
		}
	}
	return h
}

// ConditionalEntropyYgivenX returns the conditional entropy H(Y|X) =
// H(X,Y) - H(X) of the joint distribution joint, measured in bits. It is the
// average remaining uncertainty about the column variable Y once the row
// variable X is known.
func ConditionalEntropyYgivenX(joint [][]float64) float64 {
	return JointEntropy(joint) - Entropy(MarginalX(joint))
}

// ConditionalEntropyXgivenY returns the conditional entropy H(X|Y) =
// H(X,Y) - H(Y) of the joint distribution joint, measured in bits.
func ConditionalEntropyXgivenY(joint [][]float64) float64 {
	return JointEntropy(joint) - Entropy(MarginalY(joint))
}

// MutualInformation returns the mutual information I(X;Y) =
// H(X) + H(Y) - H(X,Y) of the joint distribution joint, measured in bits. It is
// non-negative, symmetric in X and Y, and zero exactly when X and Y are
// independent.
func MutualInformation(joint [][]float64) float64 {
	i := Entropy(MarginalX(joint)) + Entropy(MarginalY(joint)) - JointEntropy(joint)
	if i < 0 { // clamp tiny negative round-off
		return 0
	}
	return i
}

// NormalizedMutualInformation returns the mutual information of joint scaled to
// [0,1] by the geometric mean of the marginal entropies:
// I(X;Y) / sqrt(H(X) H(Y)). If either marginal entropy is zero the value is
// defined to be zero.
func NormalizedMutualInformation(joint [][]float64) float64 {
	hx := Entropy(MarginalX(joint))
	hy := Entropy(MarginalY(joint))
	denom := math.Sqrt(hx * hy)
	if denom <= 0 {
		return 0
	}
	return MutualInformation(joint) / denom
}

// VariationOfInformation returns the variation of information
// VI(X;Y) = H(X|Y) + H(Y|X) = 2 H(X,Y) - H(X) - H(Y) of the joint distribution
// joint, measured in bits. It is a true metric on the space of partitions and
// is zero exactly when X and Y determine each other.
func VariationOfInformation(joint [][]float64) float64 {
	vi := 2*JointEntropy(joint) - Entropy(MarginalX(joint)) - Entropy(MarginalY(joint))
	if vi < 0 {
		return 0
	}
	return vi
}

// KLDivergence returns the Kullback-Leibler divergence D(p||q) =
// sum_i p_i log2(p_i/q_i) of distribution p from reference distribution q,
// measured in bits. It is non-negative and zero exactly when p equals q. If any
// outcome has p_i > 0 while q_i == 0 the divergence is infinite. p and q must
// have equal length; a mismatch returns NaN.
func KLDivergence(p, q []float64) float64 {
	if len(p) != len(q) {
		return math.NaN()
	}
	var d float64
	for i, pi := range p {
		if pi <= 0 {
			continue
		}
		if q[i] <= 0 {
			return math.Inf(1)
		}
		d += pi * infotheoryLog2(pi/q[i])
	}
	return d
}

// KLDivergenceNat returns the Kullback-Leibler divergence D(p||q) measured in
// nats (natural-logarithm units).
func KLDivergenceNat(p, q []float64) float64 {
	return BitsToNats(KLDivergence(p, q))
}

// CrossEntropy returns the cross entropy H(p,q) = -sum_i p_i log2 q_i between a
// true distribution p and a model distribution q, measured in bits. It equals
// H(p) + D(p||q) and is infinite if some outcome has p_i > 0 while q_i == 0.
// p and q must have equal length; a mismatch returns NaN.
func CrossEntropy(p, q []float64) float64 {
	if len(p) != len(q) {
		return math.NaN()
	}
	var c float64
	for i, pi := range p {
		if pi <= 0 {
			continue
		}
		if q[i] <= 0 {
			return math.Inf(1)
		}
		c -= pi * infotheoryLog2(q[i])
	}
	return c
}

// CrossEntropyNat returns the cross entropy H(p,q) measured in nats.
func CrossEntropyNat(p, q []float64) float64 {
	return BitsToNats(CrossEntropy(p, q))
}

// JensenShannonDivergence returns the Jensen-Shannon divergence between
// distributions p and q, measured in bits: JSD = 0.5 D(p||m) + 0.5 D(q||m)
// with m = (p+q)/2. Unlike the KL divergence it is symmetric, always finite,
// and bounded in [0,1] bits. p and q must have equal length; a mismatch returns
// NaN.
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

// JensenShannonDistance returns the square root of the JensenShannonDivergence
// of p and q. This quantity is a true metric satisfying the triangle
// inequality, measured in sqrt-bits.
func JensenShannonDistance(p, q []float64) float64 {
	return math.Sqrt(JensenShannonDivergence(p, q))
}
