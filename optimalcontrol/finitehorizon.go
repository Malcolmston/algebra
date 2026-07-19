package optimalcontrol

// FiniteHorizonDiscrete holds the time-varying gains and cost matrices of a
// finite-horizon discrete LQR problem, indexed from stage 0 (initial) to N
// (terminal).
type FiniteHorizonDiscrete struct {
	// P holds the cost-to-go matrices P[0..N]; P[N] is the terminal weight.
	P []*Matrix
	// K holds the feedback gains K[0..N-1] with u_k = −K[k] x_k.
	K []*Matrix
}

// FiniteHorizonLQRDiscrete solves the finite-horizon discrete LQR problem for
// x_{k+1} = A x_k + B u_k minimizing xₙᵀ Qf xₙ + Σ_{k<N} (xₖᵀ Q xₖ + uₖᵀ R uₖ)
// by backward Riccati recursion. It returns the sequence of cost-to-go matrices
// and time-varying gains.
func FiniteHorizonLQRDiscrete(a, b, q, r, qf *Matrix, n int) (*FiniteHorizonDiscrete, error) {
	if n < 1 {
		return nil, ErrDim
	}
	ps := make([]*Matrix, n+1)
	ks := make([]*Matrix, n)
	ps[n] = qf.Clone()
	for k := n - 1; k >= 0; k-- {
		gain, err := DiscreteGain(a, b, r, ps[k+1])
		if err != nil {
			return nil, err
		}
		ks[k] = gain
		next, err := dareStep(a, b, q, r, ps[k+1])
		if err != nil {
			return nil, err
		}
		ps[k] = next.Symmetrize()
	}
	return &FiniteHorizonDiscrete{P: ps, K: ks}, nil
}

// SimulateDiscreteLQR simulates the closed-loop discrete system under the
// time-varying gains produced by FiniteHorizonLQRDiscrete, starting from x0 and
// returning the state trajectory x[0..N].
func SimulateDiscreteLQR(a, b *Matrix, fh *FiniteHorizonDiscrete, x0 []float64) [][]float64 {
	n := len(fh.K)
	traj := make([][]float64, n+1)
	x := append([]float64{}, x0...)
	traj[0] = append([]float64{}, x...)
	for k := 0; k < n; k++ {
		u := fh.K[k].MulVec(x)
		for i := range u {
			u[i] = -u[i]
		}
		ax := a.MulVec(x)
		bu := b.MulVec(u)
		next := make([]float64, len(x))
		for i := range next {
			next[i] = ax[i] + bu[i]
		}
		x = next
		traj[k+1] = append([]float64{}, x...)
	}
	return traj
}

// FiniteHorizonContinuous holds the sampled solution of the continuous
// finite-horizon Riccati differential equation on a time grid.
type FiniteHorizonContinuous struct {
	// Times holds the sample times, ascending from 0 to T.
	Times []float64
	// P holds the Riccati matrix at each sample time.
	P []*Matrix
}

// FiniteHorizonLQRContinuous integrates the matrix Riccati differential equation
//
//	−dP/dt = Aᵀ P + P A − P B R⁻¹ Bᵀ P + Q,   P(T) = Qf,
//
// backward from the terminal time T to 0 using classical fourth-order
// Runge–Kutta with the given number of steps. Samples are returned in ascending
// time order.
func FiniteHorizonLQRContinuous(a, b, q, r, qf *Matrix, tFinal float64, steps int) (*FiniteHorizonContinuous, error) {
	if steps < 1 {
		return nil, ErrDim
	}
	rinv, err := Inverse(r)
	if err != nil {
		return nil, err
	}
	s := b.Mul(rinv).Mul(b.Transpose())
	at := a.Transpose()
	// dP/dt = -(AᵀP + PA - P S P + Q). Integrate backward with negative dt.
	deriv := func(p *Matrix) *Matrix {
		return at.Mul(p).Plus(p.Mul(a)).Minus(p.Mul(s).Mul(p)).Plus(q).Neg()
	}
	dt := tFinal / float64(steps)
	times := make([]float64, steps+1)
	ps := make([]*Matrix, steps+1)
	// Integrate from t=T (index steps) down to t=0 (index 0).
	p := qf.Clone()
	times[steps] = tFinal
	ps[steps] = p.Clone()
	h := -dt
	for i := steps; i > 0; i-- {
		k1 := deriv(p)
		k2 := deriv(p.Plus(k1.Scale(h / 2)))
		k3 := deriv(p.Plus(k2.Scale(h / 2)))
		k4 := deriv(p.Plus(k3.Scale(h)))
		incr := k1.Plus(k2.Scale(2)).Plus(k3.Scale(2)).Plus(k4).Scale(h / 6)
		p = p.Plus(incr).Symmetrize()
		times[i-1] = float64(i-1) * dt
		ps[i-1] = p.Clone()
	}
	return &FiniteHorizonContinuous{Times: times, P: ps}, nil
}

// GainAt returns the continuous LQR feedback gain K(t) = R⁻¹ Bᵀ P(t) at the
// grid index i of a FiniteHorizonContinuous solution.
func (fh *FiniteHorizonContinuous) GainAt(b, r *Matrix, i int) (*Matrix, error) {
	return ContinuousGain(b, r, fh.P[i])
}
