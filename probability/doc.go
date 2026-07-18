// Package probability provides tools for working with discrete probability
// spaces, random variables, and finite-state Markov chains using only the Go
// standard library.
//
// The central type is [Distribution], a finitely-supported discrete probability
// distribution over real-valued outcomes. It offers the full complement of
// summary statistics (expectation, variance, higher moments, skewness, excess
// kurtosis, entropy), the standard generating functions (moment-generating,
// cumulant-generating, probability-generating and the characteristic function),
// cumulative and quantile queries, random-variable transforms (affine maps and
// arbitrary functions), and convolution of independent sums.
//
// [JointDistribution] models a bivariate discrete distribution and exposes its
// marginals, conditionals, covariance and correlation. A small set of
// package-level functions ([Bayes], [BayesPosterior], [TotalProbability],
// [ConditionalProbability], [UnionProbability] and friends) cover the elementary
// laws of probability and Bayesian updating.
//
// [MarkovChain] wraps a row-stochastic transition matrix and computes n-step
// transition matrices, the evolution of a distribution, stationary
// distributions of irreducible chains, and the absorbing-chain quantities
// (fundamental matrix, expected steps to absorption, absorption probabilities).
//
// All routines are deterministic and depend only on the standard library.
// Distributions produced by the constructors and transforms keep their support
// sorted in ascending order with duplicate outcomes merged, so probability-mass
// and cumulative queries are unambiguous.
package probability

import (
	"fmt"
	"math"
	"sort"
)

// probabilityTol is the default absolute tolerance used when validating that a
// probability vector sums to one and when comparing floating-point outcomes.
const probabilityTol = 1e-9

// probabilityErrorf builds a formatted error. It centralizes error creation so
// that every exported routine reports failures in a consistent style.
func probabilityErrorf(format string, args ...any) error {
	return fmt.Errorf(format, args...)
}

// probabilityMerge combines a parallel pair of outcome and probability slices
// into a canonical distribution representation: outcomes are sorted in
// ascending order and probabilities of equal outcomes are summed so that each
// outcome appears exactly once. The input slices must have equal length; the
// originals are not modified.
func probabilityMerge(outcomes, probs []float64) ([]float64, []float64) {
	type pair struct {
		x, p float64
	}
	pairs := make([]pair, len(outcomes))
	for i := range outcomes {
		pairs[i] = pair{outcomes[i], probs[i]}
	}
	sort.Slice(pairs, func(i, j int) bool { return pairs[i].x < pairs[j].x })

	outs := make([]float64, 0, len(pairs))
	ps := make([]float64, 0, len(pairs))
	for i, pr := range pairs {
		if i > 0 && pr.x == outs[len(outs)-1] {
			ps[len(ps)-1] += pr.p
			continue
		}
		outs = append(outs, pr.x)
		ps = append(ps, pr.p)
	}
	return outs, ps
}

// probabilitySum returns the sum of a slice of float64 values.
func probabilitySum(xs []float64) float64 {
	s := 0.0
	for _, x := range xs {
		s += x
	}
	return s
}

// probabilityAbs returns the absolute value of x. It exists so callers avoid a
// dependency edge on math in hot inner loops kept in other files.
func probabilityAbs(x float64) float64 { return math.Abs(x) }
