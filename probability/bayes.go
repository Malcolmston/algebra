package probability

import "math"

// Complement returns the probability of the complementary event, 1 - p. It
// returns NaN if p is outside [0, 1].
func Complement(p float64) float64 {
	if p < 0 || p > 1 || math.IsNaN(p) {
		return math.NaN()
	}
	return 1 - p
}

// ConditionalProbability returns P(A | B) = P(A ∩ B) / P(B). It returns an error
// if either probability is outside [0, 1], if the joint probability exceeds
// P(B), or if P(B) is zero.
func ConditionalProbability(pAandB, pB float64) (float64, error) {
	if pAandB < 0 || pAandB > 1 || pB < 0 || pB > 1 {
		return 0, probabilityErrorf("ConditionalProbability: probability out of range")
	}
	if pB == 0 {
		return 0, probabilityErrorf("ConditionalProbability: conditioning on a zero-probability event")
	}
	if pAandB > pB+probabilityTol {
		return 0, probabilityErrorf("ConditionalProbability: P(A∩B)=%g exceeds P(B)=%g", pAandB, pB)
	}
	return pAandB / pB, nil
}

// JointProbabilityIndependent returns P(A ∩ B) = P(A)·P(B) under the assumption
// that A and B are independent.
func JointProbabilityIndependent(pA, pB float64) float64 {
	return pA * pB
}

// UnionProbability returns P(A ∪ B) = P(A) + P(B) - P(A ∩ B) via the
// inclusion-exclusion principle.
func UnionProbability(pA, pB, pAandB float64) float64 {
	return pA + pB - pAandB
}

// AreIndependent reports whether events A and B are independent, i.e. whether
// P(A ∩ B) equals P(A)·P(B) within [probabilityTol].
func AreIndependent(pA, pB, pAandB float64) bool {
	return probabilityAbs(pAandB-pA*pB) <= probabilityTol
}

// OddsFromProbability converts a probability p in [0, 1) to odds p / (1 - p). It
// returns +Inf as p approaches one and NaN if p is outside [0, 1].
func OddsFromProbability(p float64) float64 {
	if p < 0 || p > 1 || math.IsNaN(p) {
		return math.NaN()
	}
	if p == 1 {
		return math.Inf(1)
	}
	return p / (1 - p)
}

// ProbabilityFromOdds converts non-negative odds o to the probability
// o / (1 + o). It returns NaN for negative odds.
func ProbabilityFromOdds(o float64) float64 {
	if o < 0 || math.IsNaN(o) {
		return math.NaN()
	}
	if math.IsInf(o, 1) {
		return 1
	}
	return o / (1 + o)
}

// Bayes returns the posterior probability P(A | B) via Bayes' theorem,
// P(B | A)·P(A) / P(B). It returns an error if any argument is outside [0, 1] or
// if P(B) is zero.
func Bayes(pBgivenA, pA, pB float64) (float64, error) {
	if pBgivenA < 0 || pBgivenA > 1 || pA < 0 || pA > 1 || pB < 0 || pB > 1 {
		return 0, probabilityErrorf("Bayes: probability out of range")
	}
	if pB == 0 {
		return 0, probabilityErrorf("Bayes: P(B) is zero")
	}
	return pBgivenA * pA / pB, nil
}

// TotalProbability returns P(B) = Σ_i priors[i]·likelihoods[i] via the law of
// total probability, where priors is a partition {A_i} with P(A_i) = priors[i]
// and likelihoods[i] = P(B | A_i). It returns an error if the slices differ in
// length, are empty, or the priors do not sum to one within [probabilityTol].
func TotalProbability(priors, likelihoods []float64) (float64, error) {
	if len(priors) != len(likelihoods) {
		return 0, probabilityErrorf("TotalProbability: length mismatch %d != %d", len(priors), len(likelihoods))
	}
	if len(priors) == 0 {
		return 0, probabilityErrorf("TotalProbability: empty partition")
	}
	if s := probabilitySum(priors); probabilityAbs(s-1) > probabilityTol {
		return 0, probabilityErrorf("TotalProbability: priors sum to %g, not 1", s)
	}
	sum := 0.0
	for i := range priors {
		sum += priors[i] * likelihoods[i]
	}
	return sum, nil
}

// BayesPosterior returns the full posterior distribution over a partition of
// hypotheses {A_i}: posterior[i] = priors[i]·likelihoods[i] / Σ_k
// priors[k]·likelihoods[k], where likelihoods[i] = P(evidence | A_i). It returns
// an error if the slices differ in length, are empty, the priors do not sum to
// one, or the total evidence probability is zero.
func BayesPosterior(priors, likelihoods []float64) ([]float64, error) {
	if len(priors) != len(likelihoods) {
		return nil, probabilityErrorf("BayesPosterior: length mismatch %d != %d", len(priors), len(likelihoods))
	}
	if len(priors) == 0 {
		return nil, probabilityErrorf("BayesPosterior: empty partition")
	}
	if s := probabilitySum(priors); probabilityAbs(s-1) > probabilityTol {
		return nil, probabilityErrorf("BayesPosterior: priors sum to %g, not 1", s)
	}
	joint := make([]float64, len(priors))
	total := 0.0
	for i := range priors {
		joint[i] = priors[i] * likelihoods[i]
		total += joint[i]
	}
	if total <= 0 {
		return nil, probabilityErrorf("BayesPosterior: total evidence probability is zero")
	}
	for i := range joint {
		joint[i] /= total
	}
	return joint, nil
}
