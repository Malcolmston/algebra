package linprog

import "math"

// QP is a convex quadratic program:
//
//	minimize   1/2 x^T Q x + C · x
//	subject to A x = B          (equality constraints)
//	           G x <= H         (inequality constraints)
//
// Q must be symmetric positive semidefinite for the problem to be convex. Any
// of A/B and G/H may be empty. Variables are free (not sign restricted) unless
// bounds are supplied as rows of G.
type QP struct {
	// Q is the symmetric n-by-n Hessian of the quadratic objective.
	Q [][]float64
	// C is the linear objective coefficient vector, length n.
	C []float64
	// A is the equality constraint matrix (may be nil).
	A [][]float64
	// B is the equality right-hand side (may be nil).
	B []float64
	// G is the inequality constraint matrix for G x <= H (may be nil).
	G [][]float64
	// H is the inequality right-hand side (may be nil).
	H []float64
}

// NumVars returns the number of decision variables.
func (qp QP) NumVars() int { return len(qp.C) }

// Objective returns 1/2 x^T Q x + C · x.
func (qp QP) Objective(x []float64) float64 {
	qx := MatVec(qp.Q, x)
	return 0.5*Dot(x, qx) + Dot(qp.C, x)
}

// Gradient returns the objective gradient Q x + C at x.
func (qp QP) Gradient(x []float64) []float64 {
	g := MatVec(qp.Q, x)
	for i := range g {
		g[i] += qp.C[i]
	}
	return g
}

// Feasible reports whether x satisfies A x = B (within tol) and G x <= H
// (within tol).
func (qp QP) Feasible(x []float64, tol float64) bool {
	if len(qp.A) > 0 {
		ax := MatVec(qp.A, x)
		for i := range qp.B {
			if math.Abs(ax[i]-qp.B[i]) > tol {
				return false
			}
		}
	}
	if len(qp.G) > 0 {
		gx := MatVec(qp.G, x)
		for i := range qp.H {
			if gx[i] > qp.H[i]+tol {
				return false
			}
		}
	}
	return true
}

// SolveQPEquality solves the equality-constrained convex QP
//
//	minimize 1/2 x^T Q x + c · x  subject to  A x = b
//
// by forming and solving the single symmetric KKT linear system
//
//	[ Q  A^T ] [ x ]   [ -c ]
//	[ A   0  ] [ λ ] = [  b ]
//
// It returns the optimizer x and the Lagrange multipliers λ (one per equality
// row). The Q block must make the KKT matrix nonsingular (Q positive definite
// on the null space of A); otherwise [ErrSingular] is returned.
func SolveQPEquality(q [][]float64, c []float64, a [][]float64, b []float64) (x, lambda []float64, err error) {
	n := len(c)
	me := len(a)
	dim := n + me
	kkt := make([][]float64, dim)
	for i := range kkt {
		kkt[i] = make([]float64, dim)
	}
	rhs := make([]float64, dim)
	// Q block and -c.
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			kkt[i][j] = q[i][j]
		}
		rhs[i] = -c[i]
	}
	// A and A^T blocks and b.
	for i := 0; i < me; i++ {
		for j := 0; j < n; j++ {
			kkt[n+i][j] = a[i][j] // A
			kkt[j][n+i] = a[i][j] // A^T
		}
		rhs[n+i] = b[i]
	}
	sol, err := SolveLinear(kkt, rhs)
	if err != nil {
		return nil, nil, err
	}
	return sol[:n], sol[n:], nil
}

// QPSolution is the outcome of solving a [QP] with [SolveQP].
type QPSolution struct {
	// Status classifies the outcome.
	Status Status
	// X is the optimizer.
	X []float64
	// LambdaEq holds the multipliers of the equality constraints.
	LambdaEq []float64
	// MuIneq holds the (nonnegative) multipliers of the inequality
	// constraints; inactive constraints have multiplier zero.
	MuIneq []float64
	// Objective is the optimal objective value.
	Objective float64
	// Iterations counts active-set iterations performed.
	Iterations int
}

// SolveQP solves the convex [QP] with inequality (and optional equality)
// constraints by a primal active-set method starting from the feasible point
// x0. At each iteration it solves an equality-constrained QP over the current
// working set, either stepping toward the subproblem solution (adding a
// blocking constraint) or dropping a constraint whose multiplier is negative,
// until the KKT conditions hold.
//
// x0 must be feasible (see [QP.Feasible]); this is the standard requirement of
// a primal active-set method. The returned [QPSolution] reports the optimizer
// and the equality/inequality multipliers.
func SolveQP(qp QP, x0 []float64) QPSolution {
	n := qp.NumVars()
	me := len(qp.A)
	mi := len(qp.G)
	maxIter := 50*(n+mi) + 200

	x := append([]float64(nil), x0...)
	// Working set of inequality constraints: those active at x0.
	working := make([]bool, mi)
	for i := 0; i < mi; i++ {
		if math.Abs(Dot(qp.G[i], x)-qp.H[i]) <= 1e-9 {
			working[i] = true
		}
	}

	iters := 0
	for iters < maxIter {
		iters++
		// Assemble the equality system: all A rows plus active G rows.
		var arows [][]float64
		var actIdx []int // inequality index for each active row (-1 for equality)
		for i := 0; i < me; i++ {
			arows = append(arows, qp.A[i])
			actIdx = append(actIdx, -1)
		}
		for i := 0; i < mi; i++ {
			if working[i] {
				arows = append(arows, qp.G[i])
				actIdx = append(actIdx, i)
			}
		}
		// Subproblem: minimize 1/2 p^T Q p + g^T p s.t. arows p = 0.
		g := qp.Gradient(x)
		p, mult, err := SolveQPEquality(qp.Q, g, arows, make([]float64, len(arows)))
		if err != nil {
			return QPSolution{Status: StatusIterations, X: x, Objective: qp.Objective(x), Iterations: iters}
		}
		if InfNorm(p) <= 1e-9 {
			// p ~ 0: candidate optimum. Check inequality multipliers.
			minMu := 0.0
			minRow := -1
			for r, gi := range actIdx {
				if gi < 0 {
					continue
				}
				if mult[r] < minMu-1e-12 {
					minMu = mult[r]
					minRow = r
				}
			}
			if minRow == -1 {
				// All active inequality multipliers nonnegative: optimal.
				mu := make([]float64, mi)
				lam := make([]float64, me)
				for r, gi := range actIdx {
					if gi < 0 {
						lam[r] = mult[r]
					} else {
						mu[gi] = mult[r]
					}
				}
				return QPSolution{
					Status:     StatusOptimal,
					X:          x,
					LambdaEq:   lam,
					MuIneq:     mu,
					Objective:  qp.Objective(x),
					Iterations: iters,
				}
			}
			// Drop the most-negative-multiplier inequality from the working set.
			working[actIdx[minRow]] = false
			continue
		}
		// Ratio test: shrink the step so no inactive inequality is violated.
		alpha := 1.0
		blocking := -1
		for i := 0; i < mi; i++ {
			if working[i] {
				continue
			}
			gp := Dot(qp.G[i], p)
			if gp > 1e-12 {
				slack := qp.H[i] - Dot(qp.G[i], x)
				a := slack / gp
				if a < alpha {
					alpha = a
					blocking = i
				}
			}
		}
		for k := range x {
			x[k] += alpha * p[k]
		}
		if blocking >= 0 {
			working[blocking] = true
		}
	}
	return QPSolution{Status: StatusIterations, X: x, Objective: qp.Objective(x), Iterations: iters}
}

// KKTResidual holds the first-order (Karush-Kuhn-Tucker) optimality residuals
// of a candidate solution to a [QP]. Every field is a nonnegative violation
// magnitude; all are zero at a KKT point.
type KKTResidual struct {
	// Stationarity is the infinity norm of Q x + C + A^T λ + G^T μ.
	Stationarity float64
	// PrimalEq is the maximum violation of A x = B.
	PrimalEq float64
	// PrimalIneq is the maximum violation of G x <= H.
	PrimalIneq float64
	// DualFeas is the maximum violation of μ >= 0.
	DualFeas float64
	// CompSlack is the maximum |μ_i (G x - H)_i| product.
	CompSlack float64
}

// Max returns the largest of the five residuals.
func (k KKTResidual) Max() float64 {
	m := k.Stationarity
	for _, v := range []float64{k.PrimalEq, k.PrimalIneq, k.DualFeas, k.CompSlack} {
		if v > m {
			m = v
		}
	}
	return m
}

// Satisfied reports whether every residual is within tol.
func (k KKTResidual) Satisfied(tol float64) bool { return k.Max() <= tol }

// ComputeKKTResidual evaluates the [KKTResidual] of the candidate primal point
// x with equality multipliers lambda and inequality multipliers mu for the
// quadratic program qp.
func ComputeKKTResidual(qp QP, x, lambda, mu []float64) KKTResidual {
	var k KKTResidual
	// Stationarity: Q x + C + A^T lambda + G^T mu.
	grad := qp.Gradient(x)
	if len(qp.A) > 0 && len(lambda) == len(qp.A) {
		atl := MatTVec(qp.A, lambda)
		for i := range grad {
			grad[i] += atl[i]
		}
	}
	if len(qp.G) > 0 && len(mu) == len(qp.G) {
		gtm := MatTVec(qp.G, mu)
		for i := range grad {
			grad[i] += gtm[i]
		}
	}
	k.Stationarity = InfNorm(grad)
	// Primal equality feasibility.
	if len(qp.A) > 0 {
		ax := MatVec(qp.A, x)
		for i := range qp.B {
			if v := math.Abs(ax[i] - qp.B[i]); v > k.PrimalEq {
				k.PrimalEq = v
			}
		}
	}
	// Primal inequality feasibility and complementary slackness.
	if len(qp.G) > 0 {
		gx := MatVec(qp.G, x)
		for i := range qp.H {
			if v := gx[i] - qp.H[i]; v > k.PrimalIneq {
				k.PrimalIneq = v
			}
			if i < len(mu) {
				if v := -mu[i]; v > k.DualFeas {
					k.DualFeas = v
				}
				if v := math.Abs(mu[i] * (gx[i] - qp.H[i])); v > k.CompSlack {
					k.CompSlack = v
				}
			}
		}
	}
	return k
}
