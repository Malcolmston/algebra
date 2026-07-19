package infogeom

import "math"

// ln2 is the natural logarithm of two, used to convert nats to bits.
const ln2 = 0.6931471805599453

// NatsToBits converts an information quantity from nats to bits.
func NatsToBits(nats float64) float64 { return nats / ln2 }

// BitsToNats converts an information quantity from bits to nats.
func BitsToNats(bits float64) float64 { return bits * ln2 }

// IsProbabilityVector reports whether p is a valid probability distribution:
// non-empty, with no negative entry, summing to one within tol.
func IsProbabilityVector(p []float64, tol float64) bool {
	if len(p) == 0 {
		return false
	}
	var s float64
	for _, pi := range p {
		if pi < 0 || math.IsNaN(pi) || math.IsInf(pi, 0) {
			return false
		}
		s += pi
	}
	return math.Abs(s-1) <= tol
}

// checkProbPair validates that p and q are probability vectors of equal length.
func checkProbPair(p, q []float64) error {
	if len(p) != len(q) || len(p) == 0 {
		return ErrDim
	}
	if !IsProbabilityVector(p, probTol) || !IsProbabilityVector(q, probTol) {
		return ErrNotProb
	}
	return nil
}

// Normalize returns a probability vector proportional to the non-negative
// weights w. It returns ErrNotProb when every weight is zero or a weight is
// negative, and ErrDim when w is empty.
func Normalize(w []float64) ([]float64, error) {
	if len(w) == 0 {
		return nil, ErrDim
	}
	var s float64
	for _, wi := range w {
		if wi < 0 {
			return nil, ErrNotProb
		}
		s += wi
	}
	if s == 0 {
		return nil, ErrNotProb
	}
	out := make([]float64, len(w))
	for i, wi := range w {
		out[i] = wi / s
	}
	return out, nil
}

// Entropy returns the Shannon entropy of the distribution p in nats,
// H(p) = -sum p_i ln p_i, using the convention 0 ln 0 = 0. It returns
// ErrNotProb when p is not a probability vector.
func Entropy(p []float64) (float64, error) {
	if !IsProbabilityVector(p, probTol) {
		return 0, ErrNotProb
	}
	var h float64
	for _, pi := range p {
		if pi > 0 {
			h -= pi * math.Log(pi)
		}
	}
	return h, nil
}

// EntropyBits returns the Shannon entropy of p in bits.
func EntropyBits(p []float64) (float64, error) {
	h, err := Entropy(p)
	if err != nil {
		return 0, err
	}
	return NatsToBits(h), nil
}

// CrossEntropy returns the cross entropy H(p,q) = -sum p_i ln q_i in nats. It
// returns math.Inf(1) when p has mass where q is zero, and ErrNotProb / ErrDim
// on malformed input.
func CrossEntropy(p, q []float64) (float64, error) {
	if err := checkProbPair(p, q); err != nil {
		return 0, err
	}
	var h float64
	for i := range p {
		if p[i] > 0 {
			if q[i] <= 0 {
				return math.Inf(1), nil
			}
			h -= p[i] * math.Log(q[i])
		}
	}
	return h, nil
}

// KLDivergence returns the Kullback-Leibler divergence D(p||q) = sum p_i
// ln(p_i/q_i) in nats. It is non-negative and zero iff p == q. It returns
// math.Inf(1) when p is not absolutely continuous with respect to q, and
// ErrNotProb / ErrDim on malformed input.
func KLDivergence(p, q []float64) (float64, error) {
	if err := checkProbPair(p, q); err != nil {
		return 0, err
	}
	var d float64
	for i := range p {
		if p[i] > 0 {
			if q[i] <= 0 {
				return math.Inf(1), nil
			}
			d += p[i] * math.Log(p[i]/q[i])
		}
	}
	return d, nil
}

// KLDivergenceBits returns the Kullback-Leibler divergence D(p||q) in bits.
func KLDivergenceBits(p, q []float64) (float64, error) {
	d, err := KLDivergence(p, q)
	if err != nil {
		return 0, err
	}
	return NatsToBits(d), nil
}

// JeffreysDivergence returns the symmetric Jeffreys divergence
// J(p,q) = D(p||q) + D(q||p) in nats. It returns math.Inf(1) when either
// direction diverges.
func JeffreysDivergence(p, q []float64) (float64, error) {
	a, err := KLDivergence(p, q)
	if err != nil {
		return 0, err
	}
	b, err := KLDivergence(q, p)
	if err != nil {
		return 0, err
	}
	return a + b, nil
}

// mixture returns the equal-weight mixture (p+q)/2.
func mixture(p, q []float64) []float64 {
	m := make([]float64, len(p))
	for i := range p {
		m[i] = 0.5 * (p[i] + q[i])
	}
	return m
}

// JensenShannonDivergence returns the Jensen-Shannon divergence
// JSD(p,q) = 1/2 D(p||m) + 1/2 D(q||m) with m = (p+q)/2, in nats. It is a
// bounded, symmetric smoothing of the KL divergence lying in [0, ln 2].
func JensenShannonDivergence(p, q []float64) (float64, error) {
	if err := checkProbPair(p, q); err != nil {
		return 0, err
	}
	m := mixture(p, q)
	var d float64
	for i := range p {
		if p[i] > 0 {
			d += 0.5 * p[i] * math.Log(p[i]/m[i])
		}
		if q[i] > 0 {
			d += 0.5 * q[i] * math.Log(q[i]/m[i])
		}
	}
	return d, nil
}

// JensenShannonDistance returns the Jensen-Shannon distance, the square root of
// the Jensen-Shannon divergence, which is a true metric on distributions.
func JensenShannonDistance(p, q []float64) (float64, error) {
	d, err := JensenShannonDivergence(p, q)
	if err != nil {
		return 0, err
	}
	if d < 0 {
		d = 0
	}
	return math.Sqrt(d), nil
}

// JensenShannonDivergenceBits returns the Jensen-Shannon divergence in bits;
// it lies in [0, 1].
func JensenShannonDivergenceBits(p, q []float64) (float64, error) {
	d, err := JensenShannonDivergence(p, q)
	if err != nil {
		return 0, err
	}
	return NatsToBits(d), nil
}

// RenyiDivergence returns the Renyi divergence of order alpha,
// D_alpha(p||q) = 1/(alpha-1) ln sum p_i^alpha q_i^(1-alpha), in nats. The
// order alpha must be positive and not equal to one; the limits alpha->1 and
// alpha->0 are handled by KLDivergence and the minus-log Bhattacharyya-type
// expression respectively. It returns ErrDomain for invalid alpha.
func RenyiDivergence(p, q []float64, alpha float64) (float64, error) {
	if err := checkProbPair(p, q); err != nil {
		return 0, err
	}
	if alpha <= 0 || math.IsNaN(alpha) {
		return 0, ErrDomain
	}
	if math.Abs(alpha-1) < 1e-12 {
		return KLDivergence(p, q)
	}
	var s float64
	for i := range p {
		if p[i] > 0 {
			if q[i] <= 0 {
				if alpha > 1 {
					return math.Inf(1), nil
				}
				continue
			}
			s += math.Pow(p[i], alpha) * math.Pow(q[i], 1-alpha)
		}
	}
	if s <= 0 {
		return math.Inf(1), nil
	}
	return math.Log(s) / (alpha - 1), nil
}

// TsallisDivergence returns the Tsallis relative entropy of order q_order,
// D_q(p||r) = (1/(q_order-1)) ( sum p_i^q_order r_i^(1-q_order) - 1 ), a
// deformation of the KL divergence that it recovers as q_order -> 1. It
// returns ErrDomain when q_order is not positive.
func TsallisDivergence(p, r []float64, qOrder float64) (float64, error) {
	if err := checkProbPair(p, r); err != nil {
		return 0, err
	}
	if qOrder <= 0 || math.IsNaN(qOrder) {
		return 0, ErrDomain
	}
	if math.Abs(qOrder-1) < 1e-12 {
		return KLDivergence(p, r)
	}
	var s float64
	for i := range p {
		if p[i] > 0 {
			if r[i] <= 0 {
				if qOrder > 1 {
					return math.Inf(1), nil
				}
				continue
			}
			s += math.Pow(p[i], qOrder) * math.Pow(r[i], 1-qOrder)
		}
	}
	return (s - 1) / (qOrder - 1), nil
}

// ChiSquaredDivergence returns the Pearson chi-squared divergence
// chi^2(p||q) = sum (p_i-q_i)^2 / q_i. It returns math.Inf(1) when q vanishes
// where p does not.
func ChiSquaredDivergence(p, q []float64) (float64, error) {
	if err := checkProbPair(p, q); err != nil {
		return 0, err
	}
	var s float64
	for i := range p {
		if q[i] <= 0 {
			if p[i] > 0 {
				return math.Inf(1), nil
			}
			continue
		}
		d := p[i] - q[i]
		s += d * d / q[i]
	}
	return s, nil
}

// NeymanChiSquaredDivergence returns the Neyman chi-squared divergence
// sum (p_i-q_i)^2 / p_i, the reverse of ChiSquaredDivergence.
func NeymanChiSquaredDivergence(p, q []float64) (float64, error) {
	return ChiSquaredDivergence(q, p)
}

// BhattacharyyaCoefficient returns the Bhattacharyya coefficient
// BC(p,q) = sum sqrt(p_i q_i), which lies in [0,1] and equals one iff p == q.
func BhattacharyyaCoefficient(p, q []float64) (float64, error) {
	if err := checkProbPair(p, q); err != nil {
		return 0, err
	}
	var s float64
	for i := range p {
		s += math.Sqrt(p[i] * q[i])
	}
	return s, nil
}

// BhattacharyyaDistance returns the Bhattacharyya distance
// -ln BC(p,q) in nats, which is non-negative and zero iff p == q.
func BhattacharyyaDistance(p, q []float64) (float64, error) {
	bc, err := BhattacharyyaCoefficient(p, q)
	if err != nil {
		return 0, err
	}
	if bc <= 0 {
		return math.Inf(1), nil
	}
	return -math.Log(bc), nil
}

// HellingerDistance returns the Hellinger distance
// H(p,q) = sqrt( 1/2 sum ( sqrt(p_i) - sqrt(q_i) )^2 ), a metric in [0,1]
// related to the Bhattacharyya coefficient by H = sqrt(1-BC).
func HellingerDistance(p, q []float64) (float64, error) {
	bc, err := BhattacharyyaCoefficient(p, q)
	if err != nil {
		return 0, err
	}
	v := 1 - bc
	if v < 0 {
		v = 0
	}
	return math.Sqrt(v), nil
}

// SquaredHellingerDistance returns the squared Hellinger distance
// 1/2 sum ( sqrt(p_i) - sqrt(q_i) )^2 = 1 - BC(p,q).
func SquaredHellingerDistance(p, q []float64) (float64, error) {
	bc, err := BhattacharyyaCoefficient(p, q)
	if err != nil {
		return 0, err
	}
	v := 1 - bc
	if v < 0 {
		v = 0
	}
	return v, nil
}

// TotalVariationDistance returns the total variation distance
// TV(p,q) = 1/2 sum |p_i - q_i|, which lies in [0,1].
func TotalVariationDistance(p, q []float64) (float64, error) {
	if err := checkProbPair(p, q); err != nil {
		return 0, err
	}
	var s float64
	for i := range p {
		s += math.Abs(p[i] - q[i])
	}
	return 0.5 * s, nil
}

// FDivergence returns the f-divergence D_f(p||q) = sum q_i f(p_i/q_i) for a
// convex generator f with f(1)=0. The generator is evaluated only where q_i>0;
// mass of p where q vanishes contributes q_i f(inf) which the caller's f must
// handle (returning +Inf where appropriate). It returns ErrNotProb / ErrDim on
// malformed input.
func FDivergence(p, q []float64, f func(t float64) float64) (float64, error) {
	if err := checkProbPair(p, q); err != nil {
		return 0, err
	}
	var s float64
	for i := range p {
		if q[i] <= 0 {
			if p[i] > 0 {
				return math.Inf(1), nil
			}
			continue
		}
		s += q[i] * f(p[i]/q[i])
	}
	return s, nil
}

// AlphaDivergence returns Amari's alpha-divergence between the distributions
// p and q,
//
//	D_alpha(p||q) = 4/(1-alpha^2) ( 1 - sum p_i^((1+alpha)/2) q_i^((1-alpha)/2) )
//
// for alpha in (-1,1); the limits alpha->1 and alpha->-1 recover the KL
// divergences D(p||q) and D(q||p) respectively and are handled analytically.
// It returns ErrNotProb / ErrDim on malformed input.
func AlphaDivergence(p, q []float64, alpha float64) (float64, error) {
	if err := checkProbPair(p, q); err != nil {
		return 0, err
	}
	if math.Abs(alpha-1) < 1e-12 {
		return KLDivergence(p, q)
	}
	if math.Abs(alpha+1) < 1e-12 {
		return KLDivergence(q, p)
	}
	var s float64
	for i := range p {
		s += math.Pow(p[i], (1+alpha)/2) * math.Pow(q[i], (1-alpha)/2)
	}
	return 4 / (1 - alpha*alpha) * (1 - s), nil
}

// LpDistance returns the elementwise Lp distance ( sum |p_i-q_i|^r )^(1/r)
// between two equal-length vectors for r >= 1. It does not require the inputs
// to be normalised. It returns ErrDim on a length mismatch and ErrDomain when
// r < 1.
func LpDistance(p, q []float64, r float64) (float64, error) {
	if len(p) != len(q) || len(p) == 0 {
		return 0, ErrDim
	}
	if r < 1 {
		return 0, ErrDomain
	}
	var s float64
	for i := range p {
		s += math.Pow(math.Abs(p[i]-q[i]), r)
	}
	return math.Pow(s, 1/r), nil
}
