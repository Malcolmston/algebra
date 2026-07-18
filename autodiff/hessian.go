package autodiff

// This file exposes second-order differentiation drivers built on [HyperDual]
// numbers: single-variable second derivatives, Hessian matrices, combined
// gradient/Hessian evaluation and Hessian-vector products.

// SecondDerivative returns f”(x), the second derivative of the scalar function
// f at x, from a single hyper-dual evaluation with both first-order slots
// seeded to one.
func SecondDerivative(f func(HyperDual) HyperDual, x float64) float64 {
	return f(HyperVariable(x)).E12
}

// Derivatives2 returns f(x), f'(x) and f”(x) together from one hyper-dual
// evaluation of the scalar function f at x.
func Derivatives2(f func(HyperDual) HyperDual, x float64) (value, first, second float64) {
	r := f(HyperVariable(x))
	return r.Val, r.E1, r.E12
}

// Hessian returns the Hessian matrix of second partial derivatives of the
// scalar field f: Rⁿ → R at x. The result H is n×n with H[i][j] = ∂²f/∂xᵢ∂xⱼ.
// It seeds coordinate i in the ε₁ slot and coordinate j in the ε₂ slot for each
// pair, reading the mixed slot of the output. Symmetry is exploited so only the
// upper triangle is evaluated, taking n(n+1)/2 passes; the matrix returned is
// fully populated and symmetric.
func Hessian(f func([]HyperDual) HyperDual, x []float64) [][]float64 {
	n := len(x)
	h := make([][]float64, n)
	for i := range h {
		h[i] = make([]float64, n)
	}
	buf := make([]HyperDual, n)
	for i := 0; i < n; i++ {
		for j := i; j < n; j++ {
			for k := 0; k < n; k++ {
				buf[k] = HyperConstant(x[k])
			}
			buf[i].E1 = 1
			buf[j].E2 = 1
			v := f(buf).E12
			h[i][j] = v
			h[j][i] = v
		}
	}
	return h
}

// GradientHessian returns both the gradient and the Hessian of the scalar field
// f: Rⁿ → R at x from the same set of hyper-dual passes. The gradient's i-th
// component is read from the ε₁ slot when coordinate i carries the ε₁ seed.
func GradientHessian(f func([]HyperDual) HyperDual, x []float64) (grad []float64, hess [][]float64) {
	n := len(x)
	grad = make([]float64, n)
	hess = make([][]float64, n)
	for i := range hess {
		hess[i] = make([]float64, n)
	}
	buf := make([]HyperDual, n)
	for i := 0; i < n; i++ {
		for j := i; j < n; j++ {
			for k := 0; k < n; k++ {
				buf[k] = HyperConstant(x[k])
			}
			buf[i].E1 = 1
			buf[j].E2 = 1
			r := f(buf)
			hess[i][j] = r.E12
			hess[j][i] = r.E12
			if j == i {
				grad[i] = r.E1
			}
		}
	}
	return grad, hess
}

// HessianVectorProduct returns H(x)·v, the product of the Hessian of the scalar
// field f: Rⁿ → R at x with the vector v, without forming the full Hessian. It
// seeds the direction v into the ε₁ slots and a unit ε₂ into one coordinate at
// a time, so each of the n passes yields one component of the result. x and v
// must have equal length.
func HessianVectorProduct(f func([]HyperDual) HyperDual, x, v []float64) []float64 {
	n := len(x)
	res := make([]float64, n)
	buf := make([]HyperDual, n)
	for j := 0; j < n; j++ {
		for k := 0; k < n; k++ {
			buf[k] = HyperDual{Val: x[k], E1: v[k]}
		}
		buf[j].E2 = 1
		res[j] = f(buf).E12
	}
	return res
}
