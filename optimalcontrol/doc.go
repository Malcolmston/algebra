// Package optimalcontrol implements optimal control and dynamic programming
// algorithms using only the Go standard library.
//
// The package is self-contained: it ships its own small dense linear-algebra
// layer (the Matrix type with LU, Cholesky, least-squares, matrix exponential,
// Faddeev–LeVerrier characteristic polynomials and Durand–Kerner eigenvalues,
// symmetric Jacobi eigensolves, and Lyapunov/Sylvester solvers) and builds the
// control-theoretic algorithms on top of it.
//
// The control content spans the classical pillars of the field:
//
//   - Linear-quadratic regulators in continuous and discrete time, obtained by
//     solving the algebraic Riccati equations. The continuous CARE is solved by
//     the matrix-sign-function method applied to the Hamiltonian matrix and by
//     Kleinman's Newton iteration; the discrete DARE is solved by the Riccati
//     recursion. See SolveCARE, SolveDARE, LQRContinuous and LQRDiscrete.
//
//   - Finite-horizon control via the backward Riccati recursion (discrete) and
//     backward integration of the matrix Riccati differential equation
//     (continuous). See FiniteHorizonLQRDiscrete and FiniteHorizonLQRContinuous.
//
//   - Pontryagin's minimum principle expressed as a Hamiltonian two-point
//     boundary-value problem, solved exactly for linear-quadratic problems and
//     by indirect single shooting for nonlinear ones. See PontryaginLQ and
//     IndirectShooting.
//
//   - Hamilton–Jacobi–Bellman value iteration on state grids via a
//     semi-Lagrangian scheme. See HJBGrid1D.
//
//   - Dynamic programming for finite Markov decision processes: value
//     iteration, Gauss–Seidel value iteration, exact and iterative policy
//     evaluation, Howard's policy iteration and modified policy iteration. See
//     the MDP type.
//
//   - Kalman filtering and the linear-quadratic-Gaussian dual of the LQR,
//     including steady-state continuous and discrete filters and a recursive
//     KalmanFilter. See KalmanContinuous, KalmanDiscrete, LQGContinuous and
//     LQGDiscrete.
//
// Structural analysis (controllability, observability, stabilizability and
// detectability via the Popov–Belevitch–Hautus tests) and Lyapunov Gramians
// round out the toolkit.
//
// All algorithms are deterministic; where randomness is useful the caller
// supplies a seed and uses math/rand. Nothing outside the standard library is
// imported.
package optimalcontrol
