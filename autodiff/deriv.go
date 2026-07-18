package autodiff

// This file exposes the user-facing forward-mode differentiation drivers built
// on [Dual] numbers: single-variable derivatives, gradients, Jacobians,
// directional derivatives and partial derivatives.

// Derivative returns f'(x), the first derivative of the scalar function f at
// the point x, obtained by seeding a unit derivative and reading the output's
// derivative slot.
func Derivative(f func(Dual) Dual, x float64) float64 {
	return f(Variable(x)).Der
}

// ValueAndDerivative returns both f(x) and f'(x) from a single forward-mode
// evaluation of f at x.
func ValueAndDerivative(f func(Dual) Dual, x float64) (value, derivative float64) {
	r := f(Variable(x))
	return r.Val, r.Der
}

// Gradient returns the gradient ∇f of the scalar field f: Rⁿ → R evaluated at
// x. It performs n forward-mode passes, one per coordinate, seeding a unit
// derivative in a single input at a time. The returned slice has the same
// length as x and is freshly allocated.
func Gradient(f func([]Dual) Dual, x []float64) []float64 {
	n := len(x)
	grad := make([]float64, n)
	buf := make([]Dual, n)
	for j := 0; j < n; j++ {
		for i := 0; i < n; i++ {
			buf[i] = Constant(x[i])
		}
		buf[j] = Variable(x[j])
		grad[j] = f(buf).Der
	}
	return grad
}

// ValueAndGradient returns f(x) together with its gradient ∇f(x). The value is
// taken from the first pass, so f is assumed deterministic.
func ValueAndGradient(f func([]Dual) Dual, x []float64) (value float64, grad []float64) {
	n := len(x)
	grad = make([]float64, n)
	buf := make([]Dual, n)
	for j := 0; j < n; j++ {
		for i := 0; i < n; i++ {
			buf[i] = Constant(x[i])
		}
		buf[j] = Variable(x[j])
		r := f(buf)
		grad[j] = r.Der
		if j == 0 {
			value = r.Val
		}
	}
	return value, grad
}

// PartialDerivative returns ∂f/∂xᵢ, the partial derivative of the scalar field
// f: Rⁿ → R with respect to its i-th argument, evaluated at x. It seeds a unit
// derivative in coordinate i only and needs a single evaluation.
func PartialDerivative(f func([]Dual) Dual, x []float64, i int) float64 {
	n := len(x)
	buf := make([]Dual, n)
	for k := 0; k < n; k++ {
		buf[k] = Constant(x[k])
	}
	buf[i] = Variable(x[i])
	return f(buf).Der
}

// DirectionalDerivative returns the derivative of the scalar field f: Rⁿ → R at
// x along the direction dir, that is ∇f(x)·dir. It seeds each input's
// derivative with the corresponding component of dir, so the whole quantity is
// obtained from one evaluation. The direction need not be a unit vector; to get
// the rate of change per unit length, normalize dir first. x and dir must have
// equal length.
func DirectionalDerivative(f func([]Dual) Dual, x, dir []float64) float64 {
	n := len(x)
	buf := make([]Dual, n)
	for k := 0; k < n; k++ {
		buf[k] = NewDual(x[k], dir[k])
	}
	return f(buf).Der
}

// Jacobian returns the Jacobian matrix of the vector field f: Rⁿ → Rᵐ at x. The
// result J is an m×n matrix with J[i][j] = ∂fᵢ/∂xⱼ, where m is the length of the
// slice f returns. It performs n forward-mode passes. The number of outputs is
// inferred from the first evaluation and assumed constant.
func Jacobian(f func([]Dual) []Dual, x []float64) [][]float64 {
	n := len(x)
	buf := make([]Dual, n)
	var jac [][]float64
	for j := 0; j < n; j++ {
		for i := 0; i < n; i++ {
			buf[i] = Constant(x[i])
		}
		buf[j] = Variable(x[j])
		out := f(buf)
		if jac == nil {
			jac = make([][]float64, len(out))
			for i := range jac {
				jac[i] = make([]float64, n)
			}
		}
		for i := range out {
			jac[i][j] = out[i].Der
		}
	}
	return jac
}

// JacobianVectorProduct returns J(x)·v, the product of the Jacobian of
// f: Rⁿ → Rᵐ at x with the vector v, without forming J. It seeds the inputs
// with v and reads the derivative slots of all outputs in a single pass, which
// is the efficient forward-mode primitive for tangent propagation. x and v must
// have equal length; the result has length m.
func JacobianVectorProduct(f func([]Dual) []Dual, x, v []float64) []float64 {
	n := len(x)
	buf := make([]Dual, n)
	for k := 0; k < n; k++ {
		buf[k] = NewDual(x[k], v[k])
	}
	out := f(buf)
	res := make([]float64, len(out))
	for i := range out {
		res[i] = out[i].Der
	}
	return res
}
