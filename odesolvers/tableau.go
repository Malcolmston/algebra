package odesolvers

// ButcherTableau describes a Runge-Kutta method by its coefficient arrays.
//
// A stage derivative k_i is defined by
//
//	k_i = f(t + C[i]*h, y + h * sum_j A[i][j] * k_j),
//
// and the step advances by y_next = y + h * sum_i B[i]*k_i. When BStar is
// non-nil the method is an embedded pair: y_star = y + h * sum_i BStar[i]*k_i
// provides a lower-order companion used for error estimation.
type ButcherTableau struct {
	A     [][]float64
	B     []float64
	C     []float64
	BStar []float64 // embedded (lower-order) weights, or nil
	// Order is the order of the primary solution B. OrderStar is the order of
	// the embedded solution BStar (0 when not embedded).
	Order     int
	OrderStar int
	Name      string
}

// Stages returns the number of stages s of the method.
func (bt *ButcherTableau) Stages() int { return len(bt.C) }

// IsEmbedded reports whether the tableau carries an embedded error estimate.
func (bt *ButcherTableau) IsEmbedded() bool { return bt.BStar != nil }

// IsExplicit reports whether the tableau is explicit, i.e. A is strictly lower
// triangular (A[i][j] == 0 for j >= i).
func (bt *ButcherTableau) IsExplicit() bool {
	for i := range bt.A {
		for j := i; j < len(bt.A[i]); j++ {
			if bt.A[i][j] != 0 {
				return false
			}
		}
	}
	return true
}

// stageDerivatives computes the explicit stage derivatives k_i at (t, y) with
// step h. It assumes the tableau is explicit.
func (bt *ButcherTableau) stageDerivatives(f Field, t float64, y []float64, h float64) [][]float64 {
	s := bt.Stages()
	k := make([][]float64, s)
	n := len(y)
	for i := 0; i < s; i++ {
		yi := Clone(y)
		for j := 0; j < i; j++ {
			a := bt.A[i][j]
			if a == 0 {
				continue
			}
			for m := 0; m < n; m++ {
				yi[m] += h * a * k[j][m]
			}
		}
		k[i] = f(t+bt.C[i]*h, yi)
	}
	return k
}

// Step advances one explicit Runge-Kutta step from (t, y) by h and returns the
// new state. The tableau must be explicit.
func (bt *ButcherTableau) Step(f Field, t float64, y []float64, h float64) []float64 {
	k := bt.stageDerivatives(f, t, y, h)
	out := Clone(y)
	for i := 0; i < len(k); i++ {
		if bt.B[i] == 0 {
			continue
		}
		AXPYInPlace(h*bt.B[i], k[i], out)
	}
	return out
}

// StepEmbedded advances one explicit step and, for an embedded tableau, also
// returns the local error estimate err = y_primary - y_star. When the tableau
// is not embedded err is nil.
func (bt *ButcherTableau) StepEmbedded(f Field, t float64, y []float64, h float64) (ynext, err []float64) {
	k := bt.stageDerivatives(f, t, y, h)
	n := len(y)
	ynext = Clone(y)
	for i := 0; i < len(k); i++ {
		if bt.B[i] != 0 {
			AXPYInPlace(h*bt.B[i], k[i], ynext)
		}
	}
	if bt.BStar == nil {
		return ynext, nil
	}
	err = make([]float64, n)
	for i := 0; i < len(k); i++ {
		d := bt.B[i] - bt.BStar[i]
		if d == 0 {
			continue
		}
		for m := 0; m < n; m++ {
			err[m] += h * d * k[i][m]
		}
	}
	return ynext, err
}

// --- Explicit tableau constructors -----------------------------------------

// EulerTableau returns the tableau of the forward (explicit) Euler method,
// order 1.
func EulerTableau() *ButcherTableau {
	return &ButcherTableau{
		A:     [][]float64{{0}},
		B:     []float64{1},
		C:     []float64{0},
		Order: 1,
		Name:  "Euler",
	}
}

// MidpointTableau returns the tableau of the explicit midpoint method, order 2.
func MidpointTableau() *ButcherTableau {
	return &ButcherTableau{
		A:     [][]float64{{0, 0}, {0.5, 0}},
		B:     []float64{0, 1},
		C:     []float64{0, 0.5},
		Order: 2,
		Name:  "Midpoint",
	}
}

// HeunTableau returns the tableau of Heun's method (the explicit trapezoidal
// / improved Euler method), order 2.
func HeunTableau() *ButcherTableau {
	return &ButcherTableau{
		A:     [][]float64{{0, 0}, {1, 0}},
		B:     []float64{0.5, 0.5},
		C:     []float64{0, 1},
		Order: 2,
		Name:  "Heun",
	}
}

// RalstonTableau returns the tableau of Ralston's method, the second-order
// explicit method with minimal truncation error.
func RalstonTableau() *ButcherTableau {
	return &ButcherTableau{
		A:     [][]float64{{0, 0}, {2.0 / 3.0, 0}},
		B:     []float64{0.25, 0.75},
		C:     []float64{0, 2.0 / 3.0},
		Order: 2,
		Name:  "Ralston",
	}
}

// SSPRK3Tableau returns the strong-stability-preserving third-order explicit
// Runge-Kutta tableau (Shu-Osher).
func SSPRK3Tableau() *ButcherTableau {
	return &ButcherTableau{
		A: [][]float64{
			{0, 0, 0},
			{1, 0, 0},
			{0.25, 0.25, 0},
		},
		B:     []float64{1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0},
		C:     []float64{0, 1, 0.5},
		Order: 3,
		Name:  "SSPRK3",
	}
}

// RK4Tableau returns the classical fourth-order Runge-Kutta tableau.
func RK4Tableau() *ButcherTableau {
	return &ButcherTableau{
		A: [][]float64{
			{0, 0, 0, 0},
			{0.5, 0, 0, 0},
			{0, 0.5, 0, 0},
			{0, 0, 1, 0},
		},
		B:     []float64{1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0},
		C:     []float64{0, 0.5, 0.5, 1},
		Order: 4,
		Name:  "RK4",
	}
}

// RK38Tableau returns the "3/8-rule" fourth-order Runge-Kutta tableau, an
// alternative to RK4 with smaller error coefficients.
func RK38Tableau() *ButcherTableau {
	return &ButcherTableau{
		A: [][]float64{
			{0, 0, 0, 0},
			{1.0 / 3.0, 0, 0, 0},
			{-1.0 / 3.0, 1, 0, 0},
			{1, -1, 1, 0},
		},
		B:     []float64{1.0 / 8.0, 3.0 / 8.0, 3.0 / 8.0, 1.0 / 8.0},
		C:     []float64{0, 1.0 / 3.0, 2.0 / 3.0, 1},
		Order: 4,
		Name:  "RK38",
	}
}
