package optimalcontrol

// SolveSylvester solves the Sylvester equation A X + X B = C for X, where A is
// m×m, B is n×n and C is m×n. It uses the Kronecker-product formulation
// (I⊗A + Bᵀ⊗I) vec(X) = vec(C) and a dense linear solve, which is robust for
// the small systems arising in control design.
func SolveSylvester(a, b, c *Matrix) (*Matrix, error) {
	m := a.rows
	n := b.rows
	if a.cols != m || b.cols != n || c.rows != m || c.cols != n {
		return nil, ErrDim
	}
	im := Eye(m)
	in := Eye(n)
	lhs := Kron(in, a).Plus(Kron(b.Transpose(), im))
	x, err := Solve(lhs, c.Vec())
	if err != nil {
		return nil, err
	}
	return Unvec(x, m, n), nil
}

// SolveLyapunovContinuous solves the continuous-time Lyapunov equation
// Aᵀ X + X A + Q = 0 for the symmetric matrix X.
func SolveLyapunovContinuous(a, q *Matrix) (*Matrix, error) {
	n := a.rows
	if !a.IsSquare() || q.rows != n || q.cols != n {
		return nil, ErrDim
	}
	in := Eye(n)
	at := a.Transpose()
	lhs := Kron(in, at).Plus(Kron(at, in))
	rhs := q.Neg().Vec()
	x, err := Solve(lhs, rhs)
	if err != nil {
		return nil, err
	}
	return Unvec(x, n, n).Symmetrize(), nil
}

// SolveLyapunovDiscrete solves the discrete-time (Stein) Lyapunov equation
// Aᵀ X A − X + Q = 0 for the symmetric matrix X.
func SolveLyapunovDiscrete(a, q *Matrix) (*Matrix, error) {
	n := a.rows
	if !a.IsSquare() || q.rows != n || q.cols != n {
		return nil, ErrDim
	}
	at := a.Transpose()
	lhs := Eye(n * n).Minus(Kron(at, at))
	x, err := Solve(lhs, q.Vec())
	if err != nil {
		return nil, err
	}
	return Unvec(x, n, n).Symmetrize(), nil
}

// LyapunovContinuousResidual returns Aᵀ X + X A + Q, the residual of the
// continuous Lyapunov equation, useful for verifying a solution.
func LyapunovContinuousResidual(a, q, x *Matrix) *Matrix {
	return a.Transpose().Mul(x).Plus(x.Mul(a)).Plus(q)
}

// LyapunovDiscreteResidual returns Aᵀ X A − X + Q, the residual of the discrete
// Lyapunov equation.
func LyapunovDiscreteResidual(a, q, x *Matrix) *Matrix {
	return a.Transpose().Mul(x).Mul(a).Minus(x).Plus(q)
}

// ControllabilityGramianContinuous returns the infinite-horizon controllability
// Gramian Wc solving A Wc + Wc Aᵀ + B Bᵀ = 0 for a stable A.
func ControllabilityGramianContinuous(a, b *Matrix) (*Matrix, error) {
	bbt := b.Mul(b.Transpose())
	// A Wc + Wc Aᵀ = -B Bᵀ  ==>  (Aᵀ)ᵀ Wc + Wc (Aᵀ) = -BBᵀ, i.e. Lyapunov with Aᵀ.
	return SolveLyapunovContinuous(a.Transpose(), bbt)
}

// ObservabilityGramianContinuous returns the infinite-horizon observability
// Gramian Wo solving Aᵀ Wo + Wo A + Cᵀ C = 0 for a stable A.
func ObservabilityGramianContinuous(a, c *Matrix) (*Matrix, error) {
	ctc := c.Transpose().Mul(c)
	return SolveLyapunovContinuous(a, ctc)
}

// ControllabilityGramianDiscrete returns the discrete controllability Gramian
// solving A Wc Aᵀ − Wc + B Bᵀ = 0 for a Schur-stable A.
func ControllabilityGramianDiscrete(a, b *Matrix) (*Matrix, error) {
	bbt := b.Mul(b.Transpose())
	return SolveLyapunovDiscrete(a.Transpose(), bbt)
}

// ObservabilityGramianDiscrete returns the discrete observability Gramian
// solving Aᵀ Wo A − Wo + Cᵀ C = 0 for a Schur-stable A.
func ObservabilityGramianDiscrete(a, c *Matrix) (*Matrix, error) {
	ctc := c.Transpose().Mul(c)
	return SolveLyapunovDiscrete(a, ctc)
}
