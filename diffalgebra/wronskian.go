package diffalgebra

// DeterminantRatFunc returns the determinant of a square matrix of rational
// functions using fraction-free-free Gaussian elimination over the field Q(x).
// It returns ErrNotSquare for a non-square input and ErrEmpty for an empty one.
func DeterminantRatFunc(m [][]RatFunc) (RatFunc, error) {
	n := len(m)
	if n == 0 {
		return ZeroRatFunc(), ErrEmpty
	}
	for _, row := range m {
		if len(row) != n {
			return ZeroRatFunc(), ErrNotSquare
		}
	}
	// copy
	a := make([][]RatFunc, n)
	for i := range m {
		a[i] = make([]RatFunc, n)
		copy(a[i], m[i])
	}
	det := OneRatFunc()
	for col := 0; col < n; col++ {
		// find pivot
		piv := -1
		for r := col; r < n; r++ {
			if !a[r][col].IsZero() {
				piv = r
				break
			}
		}
		if piv == -1 {
			return ZeroRatFunc(), nil
		}
		if piv != col {
			a[piv], a[col] = a[col], a[piv]
			det = det.Neg()
		}
		det = det.Mul(a[col][col])
		inv, _ := a[col][col].Inv()
		for r := col + 1; r < n; r++ {
			if a[r][col].IsZero() {
				continue
			}
			factor := a[r][col].Mul(inv)
			for c := col; c < n; c++ {
				a[r][c] = a[r][c].Sub(factor.Mul(a[col][c]))
			}
		}
	}
	return det, nil
}

// DeterminantPoly returns the determinant of a square matrix of polynomials by
// evaluating over the field Q(x). It returns a polynomial when the true
// determinant is a polynomial (which is always the case for a polynomial
// matrix).
func DeterminantPoly(m [][]Poly) (Poly, error) {
	n := len(m)
	if n == 0 {
		return ZeroPoly(), ErrEmpty
	}
	rm := make([][]RatFunc, n)
	for i := range m {
		if len(m[i]) != n {
			return ZeroPoly(), ErrNotSquare
		}
		rm[i] = make([]RatFunc, n)
		for j := range m[i] {
			rm[i][j] = RatFuncFromPoly(m[i][j])
		}
	}
	d, err := DeterminantRatFunc(rm)
	if err != nil {
		return ZeroPoly(), err
	}
	return d.Num(), nil
}

// WronskianMatrixRatFunc returns the Wronskian matrix of the functions fs: row
// i holds the (i)-th derivatives of every function.
func WronskianMatrixRatFunc(fs []RatFunc) [][]RatFunc {
	n := len(fs)
	m := make([][]RatFunc, n)
	derivs := make([]RatFunc, n)
	copy(derivs, fs)
	for i := 0; i < n; i++ {
		m[i] = make([]RatFunc, n)
		copy(m[i], derivs)
		for j := range derivs {
			derivs[j] = derivs[j].Derivative()
		}
	}
	return m
}

// WronskianMatrixPoly returns the Wronskian matrix of the polynomials ps.
func WronskianMatrixPoly(ps []Poly) [][]Poly {
	n := len(ps)
	m := make([][]Poly, n)
	derivs := make([]Poly, n)
	copy(derivs, ps)
	for i := 0; i < n; i++ {
		m[i] = make([]Poly, n)
		copy(m[i], derivs)
		for j := range derivs {
			derivs[j] = derivs[j].Derivative()
		}
	}
	return m
}

// WronskianRatFunc returns the Wronskian determinant of the rational functions
// fs. An empty input returns ErrEmpty.
func WronskianRatFunc(fs []RatFunc) (RatFunc, error) {
	if len(fs) == 0 {
		return ZeroRatFunc(), ErrEmpty
	}
	return DeterminantRatFunc(WronskianMatrixRatFunc(fs))
}

// WronskianPoly returns the Wronskian determinant of the polynomials ps.
func WronskianPoly(ps []Poly) (Poly, error) {
	if len(ps) == 0 {
		return ZeroPoly(), ErrEmpty
	}
	return DeterminantPoly(WronskianMatrixPoly(ps))
}

// LinearlyIndependentRatFunc reports whether the rational functions are linearly
// independent over the constants, tested by a non-vanishing Wronskian.
func LinearlyIndependentRatFunc(fs []RatFunc) bool {
	w, err := WronskianRatFunc(fs)
	if err != nil {
		return false
	}
	return !w.IsZero()
}

// LinearlyIndependentPoly reports whether the polynomials are linearly
// independent over Q, tested by a non-vanishing Wronskian.
func LinearlyIndependentPoly(ps []Poly) bool {
	w, err := WronskianPoly(ps)
	if err != nil {
		return false
	}
	return !w.IsZero()
}
