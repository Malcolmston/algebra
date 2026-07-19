package exterioralgebra

// ScalarField is a smooth real-valued function on Rⁿ, used as the coefficient
// of a numerically differentiated differential form.
type ScalarField func(x []float64) float64

// FormField is a differential form on Rⁿ given pointwise: it maps a position x
// to the constant [Form] of coefficients at that point. All returned Forms are
// expected to share the same ambient dimension and grade structure.
type FormField func(x []float64) *Form

// NumGradient returns the exterior derivative of the 0-form f at the point x,
// namely the grade-1 Form Σ_i (∂f/∂x_i) dx^i, using a symmetric finite
// difference of step h on each coordinate. A typical choice is h ≈ 1e-5.
func NumGradient(f ScalarField, x []float64, h float64) *Form {
	n := len(x)
	res := New(n)
	for i := 0; i < n; i++ {
		xp := append([]float64(nil), x...)
		xm := append([]float64(nil), x...)
		xp[i] += h
		xm[i] -= h
		d := (f(xp) - f(xm)) / (2 * h)
		if d != 0 {
			res.terms[uint(1)<<uint(i)] = d
		}
	}
	return res
}

// NumExteriorDerivative returns the exterior derivative of the form field at the
// point x, evaluated with symmetric finite differences of step h. For a field
// Σ_I a_I(x) dx^I it forms Σ_I Σ_j (∂a_I/∂x_j) dx^j∧dx^I by differentiating each
// blade coefficient. The ambient dimension is taken from field(x).
func NumExteriorDerivative(field FormField, x []float64, h float64) *Form {
	base := field(x)
	n := base.n
	res := New(n)
	for j := 0; j < n; j++ {
		xp := append([]float64(nil), x...)
		xm := append([]float64(nil), x...)
		xp[j] += h
		xm[j] -= h
		deriv := field(xp).Sub(field(xm)).Scale(1 / (2 * h))
		bit := uint(1) << uint(j)
		for m, c := range deriv.terms {
			if m&bit != 0 {
				continue
			}
			res.addTerm(bit|m, float64(reorderSign(bit, m))*c)
		}
	}
	return res
}

// NumDirectionalDerivative returns the directional derivative of the scalar
// field f at x along the direction vector v, computed as a symmetric finite
// difference of step h. It equals ⟨∇f, v⟩ in the limit h→0.
func NumDirectionalDerivative(f ScalarField, x, v []float64, h float64) float64 {
	xp := append([]float64(nil), x...)
	xm := append([]float64(nil), x...)
	for i := range x {
		xp[i] += h * v[i]
		xm[i] -= h * v[i]
	}
	return (f(xp) - f(xm)) / (2 * h)
}
