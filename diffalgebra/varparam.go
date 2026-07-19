package diffalgebra

// wronskianColumnReplaced returns the Wronskian matrix of ys with column i
// replaced by the unit vector e_n = (0,...,0,1).
func wronskianColumnReplaced(ys []RatFunc, i int) [][]RatFunc {
	m := WronskianMatrixRatFunc(ys)
	n := len(ys)
	for r := 0; r < n; r++ {
		if r == n-1 {
			m[r][i] = OneRatFunc()
		} else {
			m[r][i] = ZeroRatFunc()
		}
	}
	return m
}

// VariationOfParametersIntegrands returns the derivatives u_i'(x) of the
// variation-of-parameters coefficient functions for the monic linear ODE
// L[y] = g whose fundamental system is ys and whose forcing term is g. The
// particular solution is y_p = sum_i (integral of u_i') * y_i. It returns
// ErrEmpty for no fundamental solutions and ErrSingular when the Wronskian
// vanishes identically.
func VariationOfParametersIntegrands(ys []RatFunc, g RatFunc) ([]RatFunc, error) {
	if len(ys) == 0 {
		return nil, ErrEmpty
	}
	w, err := WronskianRatFunc(ys)
	if err != nil {
		return nil, err
	}
	if w.IsZero() {
		return nil, ErrSingular
	}
	out := make([]RatFunc, len(ys))
	for i := range ys {
		wi, err := DeterminantRatFunc(wronskianColumnReplaced(ys, i))
		if err != nil {
			return nil, err
		}
		ui, err := g.Mul(wi).Div(w)
		if err != nil {
			return nil, err
		}
		out[i] = ui
	}
	return out, nil
}

// VariationOfParameters attempts to build a particular solution of the monic
// linear ODE L[y] = g from the fundamental system ys by integrating the
// variation-of-parameters integrands. It succeeds when every integrand
// integrates to a rational function (no logarithmic part); otherwise it returns
// ErrNoSolution together with the integrands via
// VariationOfParametersIntegrands. The returned RatFunc is the particular
// solution y_p.
func VariationOfParameters(ys []RatFunc, g RatFunc) (RatFunc, error) {
	integrands, err := VariationOfParametersIntegrands(ys, g)
	if err != nil {
		return ZeroRatFunc(), err
	}
	yp := ZeroRatFunc()
	for i, ui := range integrands {
		res, err := IntegrateRational(ui)
		if err != nil {
			return ZeroRatFunc(), err
		}
		if len(res.Logs) != 0 {
			return ZeroRatFunc(), ErrNoSolution
		}
		yp = yp.Add(res.Rational.Mul(ys[i]))
	}
	return yp, nil
}
