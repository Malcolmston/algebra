package optimalcontrol

// KalmanContinuousResult holds the steady-state error covariance and estimator
// gain of a continuous-time Kalman–Bucy filter.
type KalmanContinuousResult struct {
	// P is the steady-state estimation-error covariance.
	P *Matrix
	// L is the Kalman gain L = P Cᵀ V⁻¹.
	L *Matrix
}

// KalmanContinuous designs the steady-state continuous-time Kalman–Bucy filter
// for x' = A x + w, y = C x + v with process covariance W and measurement
// covariance V. It solves the filter algebraic Riccati equation
// A P + P Aᵀ − P Cᵀ V⁻¹ C P + W = 0 (the dual of the LQR CARE) and returns the
// covariance and gain.
func KalmanContinuous(a, c, w, v *Matrix) (*KalmanContinuousResult, error) {
	p, err := SolveCARE(a.Transpose(), c.Transpose(), w, v)
	if err != nil {
		return nil, err
	}
	vinv, err := Inverse(v)
	if err != nil {
		return nil, err
	}
	l := p.Mul(c.Transpose()).Mul(vinv)
	return &KalmanContinuousResult{P: p, L: l}, nil
}

// KalmanDiscreteResult holds the steady-state a-priori covariance and gain of a
// discrete-time Kalman filter.
type KalmanDiscreteResult struct {
	// P is the steady-state a-priori (predicted) error covariance.
	P *Matrix
	// L is the steady-state Kalman update gain L = P Cᵀ (C P Cᵀ + V)⁻¹.
	L *Matrix
}

// KalmanDiscrete designs the steady-state discrete-time Kalman filter for
// x_{k+1} = A x_k + w_k, y_k = C x_k + v_k with process covariance W and
// measurement covariance V. It solves the dual DARE
// P = A P Aᵀ − A P Cᵀ (C P Cᵀ + V)⁻¹ C P Aᵀ + W and returns the a-priori
// covariance and update gain.
func KalmanDiscrete(a, c, w, v *Matrix) (*KalmanDiscreteResult, error) {
	p, err := SolveDARE(a.Transpose(), c.Transpose(), w, v)
	if err != nil {
		return nil, err
	}
	ct := c.Transpose()
	inner, err := Inverse(c.Mul(p).Mul(ct).Plus(v))
	if err != nil {
		return nil, err
	}
	l := p.Mul(ct).Mul(inner)
	return &KalmanDiscreteResult{P: p, L: l}, nil
}

// LQGResult bundles the regulator and estimator designs of a linear-quadratic
// Gaussian controller.
type LQGResult struct {
	// K is the LQR state-feedback gain (u = −K x̂).
	K *Matrix
	// L is the Kalman estimator gain.
	L *Matrix
	// P is the control Riccati solution.
	P *Matrix
	// Sigma is the estimation-error covariance.
	Sigma *Matrix
}

// LQGContinuous designs a continuous-time LQG controller by combining a
// continuous LQR (weights Q, R on the plant (A, B)) with a continuous Kalman
// filter (process covariance W, measurement covariance V on (A, C)), invoking
// the separation principle.
func LQGContinuous(a, b, c, q, r, w, v *Matrix) (*LQGResult, error) {
	lqr, err := LQRContinuous(a, b, q, r)
	if err != nil {
		return nil, err
	}
	kf, err := KalmanContinuous(a, c, w, v)
	if err != nil {
		return nil, err
	}
	return &LQGResult{K: lqr.K, L: kf.L, P: lqr.P, Sigma: kf.P}, nil
}

// LQGDiscrete designs a discrete-time LQG controller by combining a discrete
// LQR with a discrete Kalman filter via the separation principle.
func LQGDiscrete(a, b, c, q, r, w, v *Matrix) (*LQGResult, error) {
	lqr, err := LQRDiscrete(a, b, q, r)
	if err != nil {
		return nil, err
	}
	kf, err := KalmanDiscrete(a, c, w, v)
	if err != nil {
		return nil, err
	}
	return &LQGResult{K: lqr.K, L: kf.L, P: lqr.P, Sigma: kf.P}, nil
}

// KalmanFilter is a recursive discrete-time Kalman filter that tracks a state
// estimate and its covariance. Construct one with NewKalmanFilter and drive it
// with alternating Predict and Update calls.
type KalmanFilter struct {
	// A is the state-transition matrix.
	A *Matrix
	// B is the (optional) control-input matrix; may be nil for no input.
	B *Matrix
	// C is the measurement matrix.
	C *Matrix
	// Q is the process-noise covariance.
	Q *Matrix
	// R is the measurement-noise covariance.
	R *Matrix
	// X is the current state estimate.
	X []float64
	// P is the current estimate covariance.
	P *Matrix
}

// NewKalmanFilter constructs a Kalman filter with the given model matrices and
// initial estimate x0 with covariance p0. B may be nil when the system has no
// control input.
func NewKalmanFilter(a, b, c, q, r *Matrix, x0 []float64, p0 *Matrix) *KalmanFilter {
	return &KalmanFilter{
		A: a, B: b, C: c, Q: q, R: r,
		X: append([]float64{}, x0...),
		P: p0.Clone(),
	}
}

// Predict advances the estimate through the process model using control input u
// (which may be nil when B is nil): x̂ ← A x̂ + B u, P ← A P Aᵀ + Q.
func (kf *KalmanFilter) Predict(u []float64) {
	x := kf.A.MulVec(kf.X)
	if kf.B != nil && u != nil {
		bu := kf.B.MulVec(u)
		for i := range x {
			x[i] += bu[i]
		}
	}
	kf.X = x
	kf.P = kf.A.Mul(kf.P).Mul(kf.A.Transpose()).Plus(kf.Q).Symmetrize()
}

// Update corrects the prediction with measurement y using the Kalman gain:
// L = P Cᵀ (C P Cᵀ + R)⁻¹, x̂ ← x̂ + L (y − C x̂), P ← (I − L C) P. It returns
// an error if the innovation covariance is singular.
func (kf *KalmanFilter) Update(y []float64) error {
	ct := kf.C.Transpose()
	s := kf.C.Mul(kf.P).Mul(ct).Plus(kf.R)
	sinv, err := Inverse(s)
	if err != nil {
		return err
	}
	l := kf.P.Mul(ct).Mul(sinv)
	// Innovation.
	cx := kf.C.MulVec(kf.X)
	innov := make([]float64, len(y))
	for i := range y {
		innov[i] = y[i] - cx[i]
	}
	corr := l.MulVec(innov)
	for i := range kf.X {
		kf.X[i] += corr[i]
	}
	n := kf.P.rows
	imlc := Eye(n).Minus(l.Mul(kf.C))
	kf.P = imlc.Mul(kf.P).Symmetrize()
	return nil
}
