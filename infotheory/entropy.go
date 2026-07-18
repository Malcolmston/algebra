package infotheory

import "math"

// infotheoryLog2 returns the base-2 logarithm of x. It is a thin wrapper kept
// so that the entropy routines read in terms of information units.
func infotheoryLog2(x float64) float64 { return math.Log2(x) }

// infotheoryXLogX returns x*log2(x) with the convention 0*log2(0) = 0. Values
// of x that are non-positive contribute zero, matching the limit used
// throughout information theory.
func infotheoryXLogX(x float64) float64 {
	if x <= 0 {
		return 0
	}
	return x * math.Log2(x)
}

// BitsToNats converts an information quantity expressed in bits (base-2 units)
// to nats (natural-logarithm units) by multiplying by ln 2.
func BitsToNats(bits float64) float64 { return bits * math.Ln2 }

// NatsToBits converts an information quantity expressed in nats to bits by
// dividing by ln 2.
func NatsToBits(nats float64) float64 { return nats / math.Ln2 }

// Entropy returns the Shannon entropy H(p) = -sum_i p_i log2 p_i of the
// discrete probability distribution p, measured in bits. Zero-probability
// outcomes contribute nothing (0 log 0 = 0). The distribution is assumed to be
// normalised; entries that are non-positive are ignored.
func Entropy(p []float64) float64 {
	var h float64
	for _, pi := range p {
		h -= infotheoryXLogX(pi)
	}
	return h
}

// EntropyBase returns the Shannon entropy of p measured in units of the given
// logarithm base (for example base 2 for bits, math.E for nats, base 10 for
// bans/hartleys). The base must be greater than one.
func EntropyBase(p []float64, base float64) float64 {
	return Entropy(p) / math.Log2(base)
}

// EntropyNat returns the Shannon entropy of p measured in nats (natural
// logarithm units): H(p) = -sum_i p_i ln p_i.
func EntropyNat(p []float64) float64 {
	return BitsToNats(Entropy(p))
}

// BinaryEntropy returns the binary entropy function
// H_b(p) = -p log2 p - (1-p) log2 (1-p), the entropy in bits of a Bernoulli
// random variable with success probability p. It is defined for p in [0,1],
// with H_b(0) = H_b(1) = 0, and peaks at H_b(0.5) = 1.
func BinaryEntropy(p float64) float64 {
	if p <= 0 || p >= 1 {
		return 0
	}
	return -infotheoryXLogX(p) - infotheoryXLogX(1-p)
}

// Surprisal returns the self-information (surprisal) -log2 p of an outcome with
// probability p, measured in bits. An impossible outcome (p <= 0) has infinite
// surprisal.
func Surprisal(p float64) float64 {
	if p <= 0 {
		return math.Inf(1)
	}
	return -infotheoryLog2(p)
}

// SurprisalBase returns the self-information of an outcome with probability p
// measured in units of the given logarithm base. The base must exceed one.
func SurprisalBase(p, base float64) float64 {
	return Surprisal(p) / math.Log2(base)
}

// SurprisalNat returns the self-information of an outcome with probability p
// measured in nats: -ln p.
func SurprisalNat(p float64) float64 {
	return BitsToNats(Surprisal(p))
}

// Perplexity returns the perplexity 2^H(p) of the distribution p, where H is
// the Shannon entropy in bits. It equals the effective number of equally likely
// outcomes and reaches its maximum, len(p), for the uniform distribution.
func Perplexity(p []float64) float64 {
	return math.Exp2(Entropy(p))
}

// NormalizedEntropy returns the entropy efficiency H(p)/log2(n), where n is the
// number of outcomes len(p). The result lies in [0,1], reaching one for the
// uniform distribution. For n < 2 it returns zero, the maximum entropy being
// itself zero.
func NormalizedEntropy(p []float64) float64 {
	n := len(p)
	if n < 2 {
		return 0
	}
	return Entropy(p) / infotheoryLog2(float64(n))
}

// Redundancy returns 1 - NormalizedEntropy(p), the fraction of the maximum
// entropy that is unused by the distribution p. It lies in [0,1] and is zero
// for the uniform distribution.
func Redundancy(p []float64) float64 {
	return 1 - NormalizedEntropy(p)
}

// GiniImpurity returns the Gini impurity 1 - sum_i p_i^2 of the distribution p,
// the probability that two independent samples fall in different categories. It
// ranges from 0 (a point mass) to 1 - 1/n (the uniform distribution over n
// outcomes).
func GiniImpurity(p []float64) float64 {
	var s float64
	for _, pi := range p {
		s += pi * pi
	}
	return 1 - s
}

// RenyiEntropy returns the Renyi entropy of order alpha of the distribution p,
// measured in bits. For alpha != 1 it is 1/(1-alpha) * log2(sum_i p_i^alpha).
// Special orders are handled by their limits: alpha == 1 yields the Shannon
// Entropy, alpha == 0 yields the HartleyEntropy (log2 of the support size), and
// alpha == +Inf yields the MinEntropy. alpha must be non-negative.
func RenyiEntropy(p []float64, alpha float64) float64 {
	switch {
	case math.IsInf(alpha, 1):
		return MinEntropy(p)
	case alpha == 1:
		return Entropy(p)
	case alpha == 0:
		return HartleyEntropy(p)
	}
	var s float64
	for _, pi := range p {
		if pi > 0 {
			s += math.Pow(pi, alpha)
		}
	}
	return infotheoryLog2(s) / (1 - alpha)
}

// CollisionEntropy returns the Renyi entropy of order two,
// -log2(sum_i p_i^2), measured in bits. It is the negative log of the
// collision (repeat) probability of two independent samples.
func CollisionEntropy(p []float64) float64 {
	var s float64
	for _, pi := range p {
		s += pi * pi
	}
	return -infotheoryLog2(s)
}

// MinEntropy returns the min-entropy -log2(max_i p_i) of the distribution p,
// measured in bits. It is the Renyi entropy in the limit of infinite order and
// gives a worst-case (most conservative) measure of uncertainty.
func MinEntropy(p []float64) float64 {
	var m float64
	for _, pi := range p {
		if pi > m {
			m = pi
		}
	}
	if m <= 0 {
		return 0
	}
	return -infotheoryLog2(m)
}

// HartleyEntropy returns the Hartley (max-) entropy log2(k), where k is the
// number of outcomes with strictly positive probability in p, measured in
// bits. It is the Renyi entropy of order zero and depends only on the support.
func HartleyEntropy(p []float64) float64 {
	var k int
	for _, pi := range p {
		if pi > 0 {
			k++
		}
	}
	if k == 0 {
		return 0
	}
	return infotheoryLog2(float64(k))
}

// TsallisEntropy returns the Tsallis entropy of order q of the distribution p,
// S_q = (1 - sum_i p_i^q) / (q - 1). For q == 1 it reduces, by its limit, to
// the Shannon entropy expressed in nats. The Tsallis entropy is a
// non-logarithmic generalisation of Shannon entropy used in non-extensive
// statistics.
func TsallisEntropy(p []float64, q float64) float64 {
	if q == 1 {
		return EntropyNat(p)
	}
	var s float64
	for _, pi := range p {
		if pi > 0 {
			s += math.Pow(pi, q)
		}
	}
	return (1 - s) / (q - 1)
}
