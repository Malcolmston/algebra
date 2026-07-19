// Package stochastic implements stochastic processes and stochastic
// differential equation (SDE) simulation in pure Go.
//
// The package is organized around a small number of building blocks and a
// large collection of concrete processes and estimators built on top of them:
//
//   - Random variate generation: uniform, normal, exponential, gamma, beta,
//     Poisson, binomial, geometric and many other one-dimensional
//     distributions, plus categorical sampling. All draws come from a
//     caller-supplied *math/rand.Rand, so results are fully deterministic given
//     a seed.
//
//   - Point processes: homogeneous and inhomogeneous (thinned) Poisson
//     processes, compound-Poisson processes, superposition and thinning, and
//     the associated counting paths.
//
//   - Random walks: symmetric and biased discrete walks, Gaussian and
//     continuous-time walks, reflecting and absorbing boundaries, gambler's
//     ruin probabilities and durations, and Levy flights.
//
//   - Diffusions: Brownian motion (standard, scaled, with drift), Brownian
//     bridges, geometric Brownian motion, the Ornstein-Uhlenbeck (Vasicek)
//     process, the Cox-Ingersoll-Ross process, fractional Brownian motion and
//     correlated multi-dimensional Brownian motion.
//
//   - SDE integrators: the Euler-Maruyama and Milstein schemes (scalar and
//     diagonal-noise systems), Monte-Carlo expectation with standard error, and
//     ensemble simulation.
//
//   - Martingale and first-passage utilities: stopping times, occupation
//     times, maximum drawdown, the reflection principle and the inverse
//     Gaussian first-passage law of drifted Brownian motion.
//
//   - Gillespie stochastic simulation: exact SSA and tau-leaping for
//     mass-action reaction networks, with ready-made birth-death, SIR and
//     Lotka-Volterra networks.
//
//   - Path statistics and estimation: quadratic and total variation, realized
//     volatility, log-returns, and parameter estimators for GBM and OU models.
//
// Randomness. Every routine that consumes randomness takes an explicit
// *math/rand.Rand (create one with NewRNG). Nothing in the package reads the
// wall clock or the operating-system entropy pool, so a fixed seed always
// reproduces the same output. This makes simulations testable and repeatable.
//
// Conventions. Continuous-time paths are represented by the Path type, a pair
// of equal-length Times and Values slices sampled on a regular grid. Discrete
// processes usually return plain slices. Time always increases along a path.
package stochastic
