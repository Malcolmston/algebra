package optimalcontrol

import "math"

// HamiltonianMatrix returns the 2n×2n Hamiltonian matrix
//
//	H = [[ A,  -S ],
//	     [ -Q, -Aᵀ ]],  with  S = B R⁻¹ Bᵀ,
//
// whose stable invariant subspace yields the stabilizing solution of the
// continuous-time algebraic Riccati equation.
func HamiltonianMatrix(a, b, q, r *Matrix) (*Matrix, error) {
	rinv, err := Inverse(r)
	if err != nil {
		return nil, err
	}
	s := b.Mul(rinv).Mul(b.Transpose())
	n := a.rows
	h := Zeros(2*n, 2*n)
	h.SetBlock(0, 0, a)
	h.SetBlock(0, n, s.Neg())
	h.SetBlock(n, 0, q.Neg())
	h.SetBlock(n, n, a.Transpose().Neg())
	return h, nil
}

// MatrixSign computes the matrix sign function of a square matrix with no purely
// imaginary eigenvalues, using the Newton iteration with determinantal scaling
//
//	Z_{k+1} = ½ ( c_k Z_k + c_k⁻¹ Z_k⁻¹ ),   c_k = |det Z_k|^{-1/N}.
//
// The returned matrix S satisfies S² = I and shares the eigenvectors of A, with
// eigenvalues ±1 according to the sign of the real part of A's eigenvalues.
func MatrixSign(a *Matrix, maxIter int, tol float64) (*Matrix, error) {
	if !a.IsSquare() {
		return nil, ErrDim
	}
	n := a.rows
	z := a.Clone()
	for iter := 0; iter < maxIter; iter++ {
		f, err := Factor(z)
		if err != nil {
			return nil, err
		}
		zinv, err := f.SolveMatrix(Eye(n))
		if err != nil {
			return nil, err
		}
		det := f.Det()
		c := 1.0
		if det != 0 {
			c = math.Pow(math.Abs(det), -1.0/float64(n))
		}
		next := z.Scale(c).Plus(zinv.Scale(1.0 / c)).Scale(0.5)
		diff := next.Minus(z).MaxAbs()
		z = next
		if diff < tol {
			return z, nil
		}
	}
	return z, ErrNotConverged
}

// SolveCARESign solves the continuous-time algebraic Riccati equation
//
//	Aᵀ X + X A − X B R⁻¹ Bᵀ X + Q = 0
//
// for the symmetric stabilizing solution X using the matrix-sign-function
// method applied to the Hamiltonian matrix. Q must be symmetric positive
// semidefinite and R symmetric positive definite; the pair (A, B) must be
// stabilizable and (A, Q) detectable.
func SolveCARESign(a, b, q, r *Matrix) (*Matrix, error) {
	n := a.rows
	h, err := HamiltonianMatrix(a, b, q, r)
	if err != nil {
		return nil, err
	}
	w, err := MatrixSign(h, 100, 1e-13)
	if err != nil {
		return nil, err
	}
	w11 := w.Submatrix(0, n, 0, n)
	w12 := w.Submatrix(0, n, n, 2*n)
	w21 := w.Submatrix(n, 2*n, 0, n)
	w22 := w.Submatrix(n, 2*n, n, 2*n)
	// Solve X [ (I-W11) | W12 ] = [ -W21 | (W22-I) ] in least squares by
	// transposing: [ (I-W11)ᵀ ; W12ᵀ ] Xᵀ = [ (-W21)ᵀ ; (W22-I)ᵀ ].
	imw11 := Eye(n).Minus(w11)
	w22mi := w22.Minus(Eye(n))
	pt := VStack(imw11.Transpose(), w12.Transpose())
	yt := VStack(w21.Neg().Transpose(), w22mi.Transpose())
	xt, err := LeastSquaresMatrix(pt, yt)
	if err != nil {
		return nil, err
	}
	return xt.Transpose().Symmetrize(), nil
}

// SolveCARE solves the continuous-time algebraic Riccati equation for its
// stabilizing solution. It is an alias for the robust matrix-sign method
// SolveCARESign.
func SolveCARE(a, b, q, r *Matrix) (*Matrix, error) {
	return SolveCARESign(a, b, q, r)
}

// SolveCAREKleinman solves the continuous-time algebraic Riccati equation by
// Kleinman's Newton iteration starting from a stabilizing gain k0 (so that
// A − B k0 is Hurwitz). Each step solves a Lyapunov equation, giving quadratic
// convergence to the stabilizing solution.
func SolveCAREKleinman(a, b, q, r, k0 *Matrix, maxIter int, tol float64) (*Matrix, error) {
	rinv, err := Inverse(r)
	if err != nil {
		return nil, err
	}
	k := k0.Clone()
	var x *Matrix
	for iter := 0; iter < maxIter; iter++ {
		ak := a.Minus(b.Mul(k))
		qk := q.Plus(k.Transpose().Mul(r).Mul(k))
		xNew, err := SolveLyapunovContinuous(ak, qk)
		if err != nil {
			return nil, err
		}
		kNew := rinv.Mul(b.Transpose()).Mul(xNew)
		if x != nil && xNew.Minus(x).MaxAbs() < tol {
			return xNew.Symmetrize(), nil
		}
		x = xNew
		k = kNew
	}
	if x == nil {
		return nil, ErrNotConverged
	}
	return x.Symmetrize(), ErrNotConverged
}

// CAREResidual returns Aᵀ X + X A − X B R⁻¹ Bᵀ X + Q, the residual of the
// continuous algebraic Riccati equation for a candidate solution X.
func CAREResidual(a, b, q, r, x *Matrix) (*Matrix, error) {
	rinv, err := Inverse(r)
	if err != nil {
		return nil, err
	}
	s := b.Mul(rinv).Mul(b.Transpose())
	return a.Transpose().Mul(x).Plus(x.Mul(a)).Minus(x.Mul(s).Mul(x)).Plus(q), nil
}

// SolveDAREIter solves the discrete-time algebraic Riccati equation
//
//	X = Aᵀ X A − Aᵀ X B (R + Bᵀ X B)⁻¹ Bᵀ X A + Q
//
// by fixed-point iteration of the Riccati recursion, starting from X₀ = Q.
func SolveDAREIter(a, b, q, r *Matrix, maxIter int, tol float64) (*Matrix, error) {
	x := q.Clone()
	for iter := 0; iter < maxIter; iter++ {
		xNew, err := dareStep(a, b, q, r, x)
		if err != nil {
			return nil, err
		}
		if xNew.Minus(x).MaxAbs() < tol {
			return xNew.Symmetrize(), nil
		}
		x = xNew
	}
	return x.Symmetrize(), ErrNotConverged
}

// dareStep performs one backward Riccati recursion step for the DARE.
func dareStep(a, b, q, r, x *Matrix) (*Matrix, error) {
	bt := b.Transpose()
	at := a.Transpose()
	rbxb := r.Plus(bt.Mul(x).Mul(b))
	inner, err := Inverse(rbxb)
	if err != nil {
		return nil, err
	}
	term := at.Mul(x).Mul(b).Mul(inner).Mul(bt).Mul(x).Mul(a)
	return at.Mul(x).Mul(a).Minus(term).Plus(q), nil
}

// SolveDARE solves the discrete-time algebraic Riccati equation for its
// stabilizing solution. It uses the iterative recursion SolveDAREIter with a
// generous iteration budget.
func SolveDARE(a, b, q, r *Matrix) (*Matrix, error) {
	return SolveDAREIter(a, b, q, r, 100000, 1e-13)
}

// DAREResidual returns the residual X − (Aᵀ X A − Aᵀ X B (R+BᵀXB)⁻¹ BᵀXA + Q)
// of the discrete algebraic Riccati equation for a candidate solution X.
func DAREResidual(a, b, q, r, x *Matrix) (*Matrix, error) {
	step, err := dareStep(a, b, q, r, x)
	if err != nil {
		return nil, err
	}
	return x.Minus(step), nil
}

// LQRResult bundles the gain matrix, Riccati solution and closed-loop matrix of
// a linear-quadratic regulator design.
type LQRResult struct {
	// K is the optimal feedback gain (u = −K x).
	K *Matrix
	// P is the stabilizing Riccati solution (the cost-to-go Hessian).
	P *Matrix
	// ClosedLoop is A − B K.
	ClosedLoop *Matrix
}

// LQRContinuous designs a continuous-time infinite-horizon LQR for the system
// x' = A x + B u minimizing ∫ (xᵀ Q x + uᵀ R u) dt. The optimal control is
// u = −K x with K = R⁻¹ Bᵀ P and P the stabilizing CARE solution.
func LQRContinuous(a, b, q, r *Matrix) (*LQRResult, error) {
	p, err := SolveCARE(a, b, q, r)
	if err != nil {
		return nil, err
	}
	rinv, err := Inverse(r)
	if err != nil {
		return nil, err
	}
	k := rinv.Mul(b.Transpose()).Mul(p)
	return &LQRResult{K: k, P: p, ClosedLoop: a.Minus(b.Mul(k))}, nil
}

// LQRDiscrete designs a discrete-time infinite-horizon LQR for the system
// x_{k+1} = A x_k + B u_k minimizing Σ (xᵀ Q x + uᵀ R u). The optimal control is
// u = −K x with K = (R + Bᵀ P B)⁻¹ Bᵀ P A and P the stabilizing DARE solution.
func LQRDiscrete(a, b, q, r *Matrix) (*LQRResult, error) {
	p, err := SolveDARE(a, b, q, r)
	if err != nil {
		return nil, err
	}
	k, err := DiscreteGain(a, b, r, p)
	if err != nil {
		return nil, err
	}
	return &LQRResult{K: k, P: p, ClosedLoop: a.Minus(b.Mul(k))}, nil
}

// DiscreteGain returns the discrete LQR feedback gain
// K = (R + Bᵀ P B)⁻¹ Bᵀ P A for a given Riccati solution P.
func DiscreteGain(a, b, r, p *Matrix) (*Matrix, error) {
	bt := b.Transpose()
	inner, err := Inverse(r.Plus(bt.Mul(p).Mul(b)))
	if err != nil {
		return nil, err
	}
	return inner.Mul(bt).Mul(p).Mul(a), nil
}

// ContinuousGain returns the continuous LQR feedback gain K = R⁻¹ Bᵀ P.
func ContinuousGain(b, r, p *Matrix) (*Matrix, error) {
	rinv, err := Inverse(r)
	if err != nil {
		return nil, err
	}
	return rinv.Mul(b.Transpose()).Mul(p), nil
}

// LQRCostToGo returns the optimal quadratic cost xᵀ P x of steering the state x
// to the origin under an infinite-horizon LQR with Riccati solution P.
func LQRCostToGo(p *Matrix, x []float64) float64 {
	return quadForm(p, x)
}

// quadForm evaluates the quadratic form xᵀ M x.
func quadForm(m *Matrix, x []float64) float64 {
	mx := m.MulVec(x)
	var s float64
	for i := range x {
		s += x[i] * mx[i]
	}
	return s
}
