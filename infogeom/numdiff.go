package infogeom

// NumericalGradient returns a central finite-difference approximation of the
// gradient of the scalar field f at x using the step h. The returned slice has
// the same length as x.
func NumericalGradient(f func(x []float64) float64, x []float64, h float64) []float64 {
	n := len(x)
	g := make([]float64, n)
	xp := CloneVector(x)
	for i := 0; i < n; i++ {
		orig := xp[i]
		xp[i] = orig + h
		fp := f(xp)
		xp[i] = orig - h
		fm := f(xp)
		xp[i] = orig
		g[i] = (fp - fm) / (2 * h)
	}
	return g
}

// NumericalHessian returns a central finite-difference approximation of the
// Hessian matrix of the scalar field f at x using the step h. The result is
// symmetrised so that H[i][j] == H[j][i].
func NumericalHessian(f func(x []float64) float64, x []float64, h float64) [][]float64 {
	n := len(x)
	xp := CloneVector(x)
	hmat := make([][]float64, n)
	for i := range hmat {
		hmat[i] = make([]float64, n)
	}
	f0 := f(x)
	for i := 0; i < n; i++ {
		for j := i; j < n; j++ {
			var v float64
			if i == j {
				oi := xp[i]
				xp[i] = oi + h
				fp := f(xp)
				xp[i] = oi - h
				fm := f(xp)
				xp[i] = oi
				v = (fp - 2*f0 + fm) / (h * h)
			} else {
				oi, oj := xp[i], xp[j]
				xp[i], xp[j] = oi+h, oj+h
				fpp := f(xp)
				xp[i], xp[j] = oi+h, oj-h
				fpm := f(xp)
				xp[i], xp[j] = oi-h, oj+h
				fmp := f(xp)
				xp[i], xp[j] = oi-h, oj-h
				fmm := f(xp)
				xp[i], xp[j] = oi, oj
				v = (fpp - fpm - fmp + fmm) / (4 * h * h)
			}
			hmat[i][j] = v
			hmat[j][i] = v
		}
	}
	return hmat
}

// NumericalJacobian returns a central finite-difference approximation of the
// Jacobian of the vector field f (from R^n to R^m) at x using the step h. The
// result has m rows and n columns. It returns ErrDim when f returns an empty
// vector.
func NumericalJacobian(f func(x []float64) []float64, x []float64, h float64) ([][]float64, error) {
	n := len(x)
	f0 := f(x)
	m := len(f0)
	if m == 0 {
		return nil, ErrDim
	}
	jac := make([][]float64, m)
	for i := range jac {
		jac[i] = make([]float64, n)
	}
	xp := CloneVector(x)
	for j := 0; j < n; j++ {
		orig := xp[j]
		xp[j] = orig + h
		fp := f(xp)
		xp[j] = orig - h
		fm := f(xp)
		xp[j] = orig
		for i := 0; i < m; i++ {
			jac[i][j] = (fp[i] - fm[i]) / (2 * h)
		}
	}
	return jac, nil
}
