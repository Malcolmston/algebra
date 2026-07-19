package odesolvers

// AdamsBashforthCoefficients returns the coefficients of the explicit
// Adams-Bashforth method of the given order (1..5). The step is
//
//	y_{n+1} = y_n + h * sum_j b[j] * f_{n-j},
//
// so b[0] multiplies the most recent derivative f_n, b[1] multiplies f_{n-1},
// and so on. An error is returned for orders outside 1..5.
func AdamsBashforthCoefficients(order int) ([]float64, error) {
	switch order {
	case 1:
		return []float64{1}, nil
	case 2:
		return []float64{3.0 / 2.0, -1.0 / 2.0}, nil
	case 3:
		return []float64{23.0 / 12.0, -16.0 / 12.0, 5.0 / 12.0}, nil
	case 4:
		return []float64{55.0 / 24.0, -59.0 / 24.0, 37.0 / 24.0, -9.0 / 24.0}, nil
	case 5:
		return []float64{1901.0 / 720.0, -2774.0 / 720.0, 2616.0 / 720.0, -1274.0 / 720.0, 251.0 / 720.0}, nil
	default:
		return nil, ErrInvalidOrder
	}
}

// AdamsMoultonCoefficients returns the coefficients of the implicit
// Adams-Moulton method of the given order (1..5). The step is
//
//	y_{n+1} = y_n + h * ( b[-1]*f_{n+1} + sum_{j>=0} b[j]*f_{n-j} ),
//
// returned as (bNew, b) where bNew multiplies the new derivative f_{n+1} and
// b[0], b[1], ... multiply f_n, f_{n-1}, .... An error is returned for orders
// outside 1..5.
func AdamsMoultonCoefficients(order int) (bNew float64, b []float64, err error) {
	switch order {
	case 1: // backward Euler
		return 1, []float64{}, nil
	case 2: // trapezoidal
		return 1.0 / 2.0, []float64{1.0 / 2.0}, nil
	case 3:
		return 5.0 / 12.0, []float64{8.0 / 12.0, -1.0 / 12.0}, nil
	case 4:
		return 9.0 / 24.0, []float64{19.0 / 24.0, -5.0 / 24.0, 1.0 / 24.0}, nil
	case 5:
		return 251.0 / 720.0, []float64{646.0 / 720.0, -264.0 / 720.0, 106.0 / 720.0, -19.0 / 720.0}, nil
	default:
		return 0, nil, ErrInvalidOrder
	}
}

// AdamsBashforthStep advances one explicit Adams-Bashforth step of the given
// order. fhist holds the most recent derivatives in reverse chronological order
// (fhist[0] = f_n, fhist[1] = f_{n-1}, ...) and must contain at least order
// entries. It returns y_{n+1}.
func AdamsBashforthStep(order int, y []float64, fhist [][]float64, h float64) ([]float64, error) {
	b, err := AdamsBashforthCoefficients(order)
	if err != nil {
		return nil, err
	}
	if len(fhist) < order {
		return nil, ErrInvalidInput
	}
	out := Clone(y)
	for j := 0; j < order; j++ {
		AXPYInPlace(h*b[j], fhist[j], out)
	}
	return out, nil
}

// SolveAdamsBashforth integrates with the explicit Adams-Bashforth method of the
// given order (1..5) and fixed step h, bootstrapping the required derivative
// history with RK4.
func SolveAdamsBashforth(f Field, order int, t0 float64, y0 []float64, tEnd, h float64) (*Solution, error) {
	if order < 1 || order > 5 {
		return nil, ErrInvalidOrder
	}
	sol := newSolution("Adams-Bashforth", t0, y0)
	nSteps, step := stepCount(t0, tEnd, h)
	if nSteps == 0 {
		return sol, nil
	}
	y := Clone(y0)
	t := t0
	fhist := [][]float64{} // fhist[0] newest
	for i := 0; i < nSteps; i++ {
		fhist = append([][]float64{f(t, y)}, fhist...)
		if len(fhist) > order {
			fhist = fhist[:order]
		}
		k := order
		if len(fhist) < order {
			k = len(fhist)
		}
		ynext, err := AdamsBashforthStep(k, y, fhist, step)
		if err != nil {
			return sol, err
		}
		y = ynext
		t = t0 + float64(i+1)*step
		sol.push(t, y)
	}
	return sol, nil
}

// SolveABM integrates with the Adams-Bashforth-Moulton predictor-corrector of
// the given order (1..5) and fixed step h, applying corrections iterations of
// the Adams-Moulton corrector (a value of 1 gives the standard PECE scheme).
// The derivative history is bootstrapped with RK4.
func SolveABM(f Field, order, corrections int, t0 float64, y0 []float64, tEnd, h float64) (*Solution, error) {
	if order < 1 || order > 5 {
		return nil, ErrInvalidOrder
	}
	if corrections < 1 {
		corrections = 1
	}
	abCoef, err := AdamsBashforthCoefficients(order)
	if err != nil {
		return nil, err
	}
	amNew, amOld, err := AdamsMoultonCoefficients(order)
	if err != nil {
		return nil, err
	}

	sol := newSolution("Adams-Bashforth-Moulton", t0, y0)
	nSteps, step := stepCount(t0, tEnd, h)
	if nSteps == 0 {
		return sol, nil
	}
	y := Clone(y0)
	t := t0
	fhist := [][]float64{} // newest first

	for i := 0; i < nSteps; i++ {
		fn := f(t, y)
		fhist = append([][]float64{fn}, fhist...)
		if len(fhist) > order {
			fhist = fhist[:order]
		}
		k := order
		if len(fhist) < order {
			k = len(fhist)
		}
		tNew := t0 + float64(i+1)*step

		var yp []float64
		if k == order {
			// Predict with Adams-Bashforth of full order.
			yp = Clone(y)
			for j := 0; j < order; j++ {
				AXPYInPlace(step*abCoef[j], fhist[j], yp)
			}
			// Correct with Adams-Moulton of full order.
			for c := 0; c < corrections; c++ {
				fNew := f(tNew, yp)
				yc := Clone(y)
				AXPYInPlace(step*amNew, fNew, yc)
				for j := 0; j < len(amOld); j++ {
					AXPYInPlace(step*amOld[j], fhist[j], yc)
				}
				yp = yc
			}
		} else {
			// Ramp-up: use RK4 until enough history exists.
			yp = RK4Step(f, t, y, step)
		}
		y = yp
		t = tNew
		sol.push(t, y)
	}
	return sol, nil
}
