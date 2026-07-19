// Package metaheuristics implements a self-contained collection of
// metaheuristic global optimizers over continuous (R^n) and combinatorial
// search spaces, together with a suite of standard benchmark objective
// functions used to validate them.
//
// The package provides the following families of optimizers:
//
//   - Genetic algorithms with configurable selection, crossover, mutation
//     and elitism operators.
//   - Simulated annealing with a family of cooling schedules
//     (exponential, geometric, linear, logarithmic, Boltzmann and Cauchy).
//   - Particle swarm optimization (PSO).
//   - Differential evolution (DE) with the classic mutation strategies.
//   - A compact CMA-ES (covariance matrix adaptation evolution strategy).
//   - Tabu search over neighbourhoods.
//   - Ant colony optimization (ACO) for the travelling salesman problem.
//   - Harmony search.
//   - Hill climbing with random restarts.
//
// For validation the package exposes the classic multimodal and unimodal
// benchmark functions: Sphere, Rastrigin, Ackley, Rosenbrock, Griewank and
// many two-dimensional test functions, each with a documented global
// optimum.
//
// # Determinism
//
// Every stochastic routine draws its randomness from a caller-supplied seed
// via the [RNG] type, which wraps math/rand's deterministic source. No use
// is made of crypto/rand or the wall clock, so a given (seed, configuration)
// pair always produces identical output. This makes the optimizers fully
// reproducible and testable.
//
// # Conventions
//
// All optimizers minimize their objective. To maximize a function f, minimize
// -f (see [Negate]). Continuous search regions are described by a [Bounds]
// value giving per-coordinate lower and upper limits.
package metaheuristics
