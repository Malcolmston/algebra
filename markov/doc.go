// Package markov provides tools for finite-state Markov chains, hidden Markov
// models (HMMs), and Markov-chain Monte Carlo (MCMC) using only the Go standard
// library.
//
// The central type for discrete-time chains is [MarkovChain], a wrapper around a
// row-stochastic transition matrix. It computes n-step transition matrices,
// evolves probability distributions, finds stationary distributions (both by
// power iteration and by exact linear solve), classifies states into
// communicating classes, determines periodicity, irreducibility and
// ergodicity, and handles absorbing chains via the fundamental matrix, giving
// expected steps to absorption, absorption probabilities and their variances.
// It also computes mean first-passage and mean recurrence times, hitting
// probabilities to arbitrary target sets, the entropy rate, Kemeny's constant
// and the time-reversed chain.
//
// [HMM] models a discrete hidden Markov model with a transition matrix, an
// emission matrix and an initial distribution. It provides the scaled forward
// and backward algorithms, sequence likelihood and log-likelihood, the Viterbi
// most-likely-path decoder, posterior state marginals, and Baum-Welch
// (expectation-maximization) training on one or many observation sequences.
//
// The MCMC facilities include random-walk and general Metropolis-Hastings
// samplers in one and several dimensions, a Gibbs sampler driven by
// full-conditional samplers, and convergence diagnostics such as the sample
// autocorrelation function, integrated autocorrelation time, effective sample
// size, batch-means variance and the Gelman-Rubin potential scale reduction
// factor.
//
// A small dense linear-algebra layer (LU decomposition, linear solve, inverse,
// determinant, matrix powers and the standard vector/matrix norms) underpins
// the chain computations and is exported for direct use. Probability vectors
// are compared with the usual statistical distances (total variation,
// Hellinger, L1, Euclidean, Kullback-Leibler).
//
// Except for the explicitly random samplers and generators (which take a
// *math/rand.Rand), every routine is deterministic and depends only on the
// standard library.
package markov
