package optimalcontrol

// LinearSystem is a continuous- or discrete-time linear time-invariant model
// x' = A x + B u (continuous) or x_{k+1} = A x_k + B u_k (discrete), with output
// y = C x + D u.
type LinearSystem struct {
	// A is the n×n state matrix.
	A *Matrix
	// B is the n×m input matrix.
	B *Matrix
	// C is the p×n output matrix.
	C *Matrix
	// D is the p×m feedthrough matrix.
	D *Matrix
}

// NewLinearSystem builds a LinearSystem, filling C with the identity and D with
// zeros when they are nil.
func NewLinearSystem(a, b, c, d *Matrix) *LinearSystem {
	n := a.rows
	m := b.cols
	if c == nil {
		c = Eye(n)
	}
	if d == nil {
		d = Zeros(c.rows, m)
	}
	return &LinearSystem{A: a, B: b, C: c, D: d}
}

// Output returns y = C x + D u for the system.
func (s *LinearSystem) Output(x, u []float64) []float64 {
	y := s.C.MulVec(x)
	if s.D != nil && u != nil {
		du := s.D.MulVec(u)
		for i := range y {
			y[i] += du[i]
		}
	}
	return y
}

// StepDiscrete advances the discrete-time state one step: x_{k+1} = A x + B u.
func (s *LinearSystem) StepDiscrete(x, u []float64) []float64 {
	nx := s.A.MulVec(x)
	if u != nil {
		bu := s.B.MulVec(u)
		for i := range nx {
			nx[i] += bu[i]
		}
	}
	return nx
}

// SimulateDiscrete simulates the discrete-time system for the given control
// sequence, returning the state trajectory x[0..len(us)].
func (s *LinearSystem) SimulateDiscrete(x0 []float64, us [][]float64) [][]float64 {
	traj := make([][]float64, len(us)+1)
	x := VecCopy(x0)
	traj[0] = VecCopy(x)
	for k, u := range us {
		x = s.StepDiscrete(x, u)
		traj[k+1] = VecCopy(x)
	}
	return traj
}

// SimulateContinuous integrates the continuous-time system with a fixed control
// u(t) held over each step of size dt using classical RK4, returning the state
// trajectory at the sample times 0, dt, …, steps·dt.
func (s *LinearSystem) SimulateContinuous(x0 []float64, control func(t float64) []float64, dt float64, steps int) [][]float64 {
	traj := make([][]float64, steps+1)
	x := VecCopy(x0)
	traj[0] = VecCopy(x)
	deriv := func(t float64, x []float64) []float64 {
		dx := s.A.MulVec(x)
		if control != nil {
			u := control(t)
			if u != nil {
				bu := s.B.MulVec(u)
				for i := range dx {
					dx[i] += bu[i]
				}
			}
		}
		return dx
	}
	for k := 0; k < steps; k++ {
		t := float64(k) * dt
		k1 := deriv(t, x)
		k2 := deriv(t+dt/2, VecAxpy(dt/2, k1, x))
		k3 := deriv(t+dt/2, VecAxpy(dt/2, k2, x))
		k4 := deriv(t+dt, VecAxpy(dt, k3, x))
		nx := make([]float64, len(x))
		for i := range nx {
			nx[i] = x[i] + dt/6*(k1[i]+2*k2[i]+2*k3[i]+k4[i])
		}
		x = nx
		traj[k+1] = VecCopy(x)
	}
	return traj
}

// DiscretizeZOH returns the exact zero-order-hold discretization (Ad, Bd) of the
// continuous pair (A, B) for sample time dt, computed from the matrix
// exponential of the augmented block matrix [[A, B],[0, 0]].
func DiscretizeZOH(a, b *Matrix, dt float64) (ad, bd *Matrix) {
	n := a.rows
	m := b.cols
	aug := Zeros(n+m, n+m)
	aug.SetBlock(0, 0, a)
	aug.SetBlock(0, n, b)
	e := MatrixExp(aug.Scale(dt))
	ad = e.Submatrix(0, n, 0, n)
	bd = e.Submatrix(0, n, n, n+m)
	return ad, bd
}

// QuadraticCostDiscrete evaluates the discrete quadratic cost
// Σ_k (xₖᵀ Q xₖ + uₖᵀ R uₖ) + x_Nᵀ Qf x_N over a state trajectory and control
// sequence. The control sequence has one fewer entry than the state trajectory.
func QuadraticCostDiscrete(xs [][]float64, us [][]float64, q, r, qf *Matrix) float64 {
	var cost float64
	for k := 0; k < len(us); k++ {
		cost += quadForm(q, xs[k]) + quadForm(r, us[k])
	}
	cost += quadForm(qf, xs[len(xs)-1])
	return cost
}

// ClosedLoopEigenvaluesContinuous returns the eigenvalues of A − B K, the
// continuous closed-loop dynamics under state feedback K.
func ClosedLoopEigenvaluesContinuous(a, b, k *Matrix) []complex128 {
	return Eigenvalues(a.Minus(b.Mul(k)))
}

// ClosedLoopEigenvaluesDiscrete returns the eigenvalues of A − B K, the discrete
// closed-loop dynamics under state feedback K.
func ClosedLoopEigenvaluesDiscrete(a, b, k *Matrix) []complex128 {
	return Eigenvalues(a.Minus(b.Mul(k)))
}
