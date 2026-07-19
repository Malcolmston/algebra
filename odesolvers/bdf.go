package odesolvers

// BDFCoefficients returns the coefficients of the fixed-step backward-
// differentiation formula of the given order (1..6).
//
// The formula is
//
//	sum_{j=0}^{k} alpha[j] * y_{n+1-j} = beta * h * f(t_{n+1}, y_{n+1}),
//
// with alpha[0] normalized to 1. The returned alpha slice has length order+1
// (alpha[0] is the coefficient of the new, unknown value). An error is returned
// for orders outside 1..6.
func BDFCoefficients(order int) (alpha []float64, beta float64, err error) {
	switch order {
	case 1:
		return []float64{1, -1}, 1, nil
	case 2:
		return []float64{1, -4.0 / 3.0, 1.0 / 3.0}, 2.0 / 3.0, nil
	case 3:
		return []float64{1, -18.0 / 11.0, 9.0 / 11.0, -2.0 / 11.0}, 6.0 / 11.0, nil
	case 4:
		return []float64{1, -48.0 / 25.0, 36.0 / 25.0, -16.0 / 25.0, 3.0 / 25.0}, 12.0 / 25.0, nil
	case 5:
		return []float64{1, -300.0 / 137.0, 300.0 / 137.0, -200.0 / 137.0, 75.0 / 137.0, -12.0 / 137.0}, 60.0 / 137.0, nil
	case 6:
		return []float64{1, -360.0 / 147.0, 450.0 / 147.0, -400.0 / 147.0, 225.0 / 147.0, -72.0 / 147.0, 10.0 / 147.0}, 60.0 / 147.0, nil
	default:
		return nil, 0, ErrInvalidOrder
	}
}

// BDFStep advances one fixed-step BDF step of the given order. The slice hist
// holds the most recent states in reverse chronological order
// (hist[0] = y_n, hist[1] = y_{n-1}, ...) and must contain at least order
// entries. It solves the implicit relation for y_{n+1} by Newton's method and
// returns the new state at time tNew = t_n + h.
func BDFStep(f Field, order int, tNew float64, hist [][]float64, h float64) ([]float64, error) {
	alpha, beta, err := BDFCoefficients(order)
	if err != nil {
		return nil, err
	}
	if len(hist) < order {
		return nil, ErrInvalidInput
	}
	n := len(hist[0])
	// Constant part: sum_{j=1}^{k} alpha[j] * y_{n+1-j}.
	c := make([]float64, n)
	for j := 1; j <= order; j++ {
		yj := hist[j-1] // y_{n+1-j}
		for m := 0; m < n; m++ {
			c[m] += alpha[j] * yj[m]
		}
	}
	// Residual G(Y) = alpha0*Y + c - beta*h*f(tNew, Y).
	residual := func(Y []float64) []float64 {
		fv := f(tNew, Y)
		out := make([]float64, n)
		for m := 0; m < n; m++ {
			out[m] = alpha[0]*Y[m] + c[m] - beta*h*fv[m]
		}
		return out
	}
	guess := Clone(hist[0])
	return NewtonSolve(residual, guess, 1e-11, 60)
}

// SolveBDF integrates y' = f(t, y) from t0 to tEnd with the fixed-step BDF
// method of the given order (1..6) and step h. The first order-1 steps are
// bootstrapped with RK4 substepping so that the multistep recurrence has enough
// history. A non-nil error aborts integration, returning the partial Solution.
func SolveBDF(f Field, order int, t0 float64, y0 []float64, tEnd, h float64) (*Solution, error) {
	if order < 1 || order > 6 {
		return nil, ErrInvalidOrder
	}
	sol := newSolution("BDF", t0, y0)
	nSteps, step := stepCount(t0, tEnd, h)
	if nSteps == 0 {
		return sol, nil
	}

	// history[0] = latest state, history[k] older.
	history := [][]float64{Clone(y0)}
	t := t0

	for i := 0; i < nSteps; i++ {
		tNew := t0 + float64(i+1)*step
		k := order
		if len(history) < order {
			k = len(history) // ramp up the order as history accumulates
		}
		ynext, err := BDFStep(f, k, tNew, history, step)
		if err != nil {
			return sol, err
		}
		// Prepend the new state.
		history = append([][]float64{Clone(ynext)}, history...)
		if len(history) > order {
			history = history[:order]
		}
		t = tNew
		sol.push(t, ynext)
	}
	return sol, nil
}
