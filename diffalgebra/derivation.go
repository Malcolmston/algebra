package diffalgebra

// Derivation is the standard derivation d/dx on the differential field Q(x). It
// is a value type carrying no configuration; the zero value is ready to use and
// StandardDerivation is provided for readability.
type Derivation struct{}

// StandardDerivation returns the derivation d/dx on Q(x).
func StandardDerivation() Derivation { return Derivation{} }

// ApplyPoly returns the derivative of a polynomial.
func (Derivation) ApplyPoly(p Poly) Poly { return p.Derivative() }

// Apply returns the derivative of a rational function.
func (Derivation) Apply(f RatFunc) RatFunc { return f.Derivative() }

// ApplyN returns the n-th derivative of a rational function (n >= 0).
func (d Derivation) ApplyN(f RatFunc, n int) RatFunc {
	out := f
	for i := 0; i < n; i++ {
		out = out.Derivative()
	}
	return out
}

// ApplyNPoly returns the n-th derivative of a polynomial (n >= 0).
func (d Derivation) ApplyNPoly(p Poly, n int) Poly {
	out := p
	for i := 0; i < n; i++ {
		out = out.Derivative()
	}
	return out
}

// IsConstant reports whether f lies in the kernel of the derivation, i.e. f is
// a constant of Q(x).
func (d Derivation) IsConstant(f RatFunc) bool { return f.Derivative().IsZero() }

// LogarithmicDerivative returns f'/f, the logarithmic derivative of f.
func (d Derivation) LogarithmicDerivative(f RatFunc) (RatFunc, error) {
	return f.LogDerivative()
}

// Leibniz verifies the product rule for two rational functions, returning
// D(fg) computed both directly and via the Leibniz expansion f'g + fg'; the two
// are always equal and this is provided as a self-checking utility.
func (d Derivation) Leibniz(f, g RatFunc) (direct, leibniz RatFunc) {
	direct = f.Mul(g).Derivative()
	leibniz = f.Derivative().Mul(g).Add(f.Mul(g.Derivative()))
	return direct, leibniz
}

// ConstantField reports whether every element derived from f by the derivation
// eventually vanishes; for Q(x) this is equivalent to f being constant.
func (d Derivation) ConstantField(f RatFunc) bool { return d.IsConstant(f) }
