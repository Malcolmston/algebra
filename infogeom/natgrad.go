package infogeom

import "math"

// NaturalGradient returns the natural gradient F^{-1} g, the Euclidean gradient
// g preconditioned by the inverse Fisher information matrix F. It is the
// direction of steepest ascent measured in the Fisher-Rao metric. It solves the
// linear system F v = g directly rather than forming the inverse. It returns
// ErrDim on a shape mismatch and ErrSingular when F is singular.
func NaturalGradient(fisher [][]float64, grad []float64) ([]float64, error) {
	return Solve(fisher, grad)
}

// NaturalGradientStep returns the parameter update theta - lr * F^{-1} g of
// natural-gradient descent, moving against the natural gradient with learning
// rate lr. It returns ErrDim on a shape mismatch and ErrSingular when F is
// singular.
func NaturalGradientStep(theta, grad []float64, fisher [][]float64, lr float64) ([]float64, error) {
	ng, err := NaturalGradient(fisher, grad)
	if err != nil {
		return nil, err
	}
	if len(ng) != len(theta) {
		return nil, ErrDim
	}
	out := make([]float64, len(theta))
	for i := range theta {
		out[i] = theta[i] - lr*ng[i]
	}
	return out, nil
}

// NaturalGradientStepInverse returns the natural-gradient update using an
// already inverted Fisher matrix fisherInv, theta - lr * fisherInv * g. It
// returns ErrDim on a shape mismatch.
func NaturalGradientStepInverse(theta, grad []float64, fisherInv [][]float64, lr float64) ([]float64, error) {
	ng, err := MatVec(fisherInv, grad)
	if err != nil {
		return nil, err
	}
	if len(ng) != len(theta) {
		return nil, ErrDim
	}
	out := make([]float64, len(theta))
	for i := range theta {
		out[i] = theta[i] - lr*ng[i]
	}
	return out, nil
}

// DampedNaturalGradient returns the Tikhonov-damped natural gradient
// (F + lambda I)^{-1} g, which regularises an ill-conditioned or singular
// Fisher matrix by adding lambda to its diagonal. It returns ErrDim on a shape
// mismatch and ErrDomain when lambda is negative.
func DampedNaturalGradient(fisher [][]float64, grad []float64, lambda float64) ([]float64, error) {
	n, c, ok := isRectangular(fisher)
	if !ok || n != c || n != len(grad) {
		return nil, ErrDim
	}
	if lambda < 0 {
		return nil, ErrDomain
	}
	damped := CloneMatrix(fisher)
	for i := 0; i < n; i++ {
		damped[i][i] += lambda
	}
	return Solve(damped, grad)
}

// GradientDescentStep returns the ordinary Euclidean gradient-descent update
// theta - lr * g. It returns ErrDim on a length mismatch.
func GradientDescentStep(theta, grad []float64, lr float64) ([]float64, error) {
	if len(theta) != len(grad) {
		return nil, ErrDim
	}
	out := make([]float64, len(theta))
	for i := range theta {
		out[i] = theta[i] - lr*grad[i]
	}
	return out, nil
}

// NaturalGradientNorm returns the intrinsic (dual) norm sqrt( g^T F^{-1} g ) of
// a gradient g under the Fisher metric F, the magnitude of the natural gradient
// measured in the Riemannian metric. It equals sqrt of the directional
// derivative achieved by a unit natural-gradient step and controls the
// second-order change of the objective. It returns ErrSingular when F is
// singular and ErrDomain when the form is negative.
func NaturalGradientNorm(fisher [][]float64, grad []float64) (float64, error) {
	ng, err := NaturalGradient(fisher, grad)
	if err != nil {
		return 0, err
	}
	inner, err := Dot(grad, ng)
	if err != nil {
		return 0, err
	}
	if inner < 0 {
		return 0, ErrDomain
	}
	return math.Sqrt(inner), nil
}

// MirrorDescentStep performs one step of mirror descent with the negative
// entropy mirror map on the probability simplex: it multiplies each coordinate
// of the distribution p by exp(-lr * g_i) and renormalises. This is the
// entropic (exponentiated-gradient) update whose geometry is dual to the
// Fisher metric. It returns ErrDim on a length mismatch and ErrNotProb when p
// is not a probability vector.
func MirrorDescentStep(p, grad []float64, lr float64) ([]float64, error) {
	if len(p) != len(grad) {
		return nil, ErrDim
	}
	if !IsProbabilityVector(p, probTol) {
		return nil, ErrNotProb
	}
	logits := make([]float64, len(p))
	for i := range p {
		if p[i] <= 0 {
			logits[i] = -1e308
			continue
		}
		logits[i] = math.Log(p[i]) - lr*grad[i]
	}
	return Softmax(logits), nil
}
