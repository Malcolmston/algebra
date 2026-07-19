package gaussproc

import (
	"fmt"
	"math"
)

// Kernel is a positive-semidefinite covariance function on pairs of real
// vectors of equal dimension. Implementations must be symmetric, that is
// Eval(x, y) must equal Eval(y, x).
type Kernel interface {
	// Eval returns the covariance k(x, y).
	Eval(x, y []float64) float64
}

// KernelFunc adapts an ordinary function to the [Kernel] interface.
type KernelFunc func(x, y []float64) float64

// Eval reports k(x, y) by calling the underlying function.
func (f KernelFunc) Eval(x, y []float64) float64 { return f(x, y) }

// Constant is the kernel k(x, y) = Value for all inputs. It models an unknown
// but constant offset.
type Constant struct {
	Value float64
}

// NewConstant returns a [Constant] kernel with the given value.
func NewConstant(value float64) Constant { return Constant{Value: value} }

// Eval returns the constant value irrespective of x and y.
func (k Constant) Eval(x, y []float64) float64 { return k.Value }

// String returns a human-readable description of the kernel.
func (k Constant) String() string { return fmt.Sprintf("Constant(%g)", k.Value) }

// WhiteNoise is the kernel k(x, y) = Variance when x and y are identical
// vectors and 0 otherwise. It models independent observation noise.
type WhiteNoise struct {
	Variance float64
}

// NewWhiteNoise returns a [WhiteNoise] kernel with the given variance.
func NewWhiteNoise(variance float64) WhiteNoise { return WhiteNoise{Variance: variance} }

// Eval returns the noise variance when x equals y and 0 otherwise.
func (k WhiteNoise) Eval(x, y []float64) float64 {
	if VectorsEqual(x, y) {
		return k.Variance
	}
	return 0
}

// String returns a human-readable description of the kernel.
func (k WhiteNoise) String() string { return fmt.Sprintf("WhiteNoise(%g)", k.Variance) }

// RBF is the squared-exponential (radial-basis-function) kernel
// k(x, y) = Variance·exp(-‖x-y‖²/(2·LengthScale²)). It is infinitely
// differentiable and hence very smooth.
type RBF struct {
	Variance    float64
	LengthScale float64
}

// NewRBF returns an [RBF] kernel with the given signal variance and length
// scale.
func NewRBF(variance, lengthScale float64) RBF {
	return RBF{Variance: variance, LengthScale: lengthScale}
}

// Eval returns the squared-exponential covariance of x and y.
func (k RBF) Eval(x, y []float64) float64 {
	d2 := SquaredDistance(x, y)
	return k.Variance * math.Exp(-d2/(2*k.LengthScale*k.LengthScale))
}

// String returns a human-readable description of the kernel.
func (k RBF) String() string {
	return fmt.Sprintf("RBF(var=%g, ls=%g)", k.Variance, k.LengthScale)
}

// SquaredExponential is an alias constructor for the [RBF] kernel.
func SquaredExponential(variance, lengthScale float64) RBF {
	return NewRBF(variance, lengthScale)
}

// Matern12 is the Matérn kernel with smoothness ν = 1/2, equal to the
// exponential (Ornstein-Uhlenbeck) kernel
// k(x, y) = Variance·exp(-r/LengthScale) with r = ‖x-y‖.
type Matern12 struct {
	Variance    float64
	LengthScale float64
}

// NewMatern12 returns a [Matern12] kernel with the given variance and length
// scale.
func NewMatern12(variance, lengthScale float64) Matern12 {
	return Matern12{Variance: variance, LengthScale: lengthScale}
}

// Eval returns the Matérn-1/2 covariance of x and y.
func (k Matern12) Eval(x, y []float64) float64 {
	r := Distance(x, y)
	return k.Variance * math.Exp(-r/k.LengthScale)
}

// String returns a human-readable description of the kernel.
func (k Matern12) String() string {
	return fmt.Sprintf("Matern12(var=%g, ls=%g)", k.Variance, k.LengthScale)
}

// Matern32 is the Matérn kernel with smoothness ν = 3/2,
// k(x, y) = Variance·(1+√3·r/ℓ)·exp(-√3·r/ℓ) with r = ‖x-y‖ and ℓ =
// LengthScale. Sample paths are once differentiable.
type Matern32 struct {
	Variance    float64
	LengthScale float64
}

// NewMatern32 returns a [Matern32] kernel with the given variance and length
// scale.
func NewMatern32(variance, lengthScale float64) Matern32 {
	return Matern32{Variance: variance, LengthScale: lengthScale}
}

// Eval returns the Matérn-3/2 covariance of x and y.
func (k Matern32) Eval(x, y []float64) float64 {
	r := Distance(x, y)
	a := math.Sqrt(3) * r / k.LengthScale
	return k.Variance * (1 + a) * math.Exp(-a)
}

// String returns a human-readable description of the kernel.
func (k Matern32) String() string {
	return fmt.Sprintf("Matern32(var=%g, ls=%g)", k.Variance, k.LengthScale)
}

// Matern52 is the Matérn kernel with smoothness ν = 5/2,
// k(x, y) = Variance·(1+√5·r/ℓ+5r²/(3ℓ²))·exp(-√5·r/ℓ) with r = ‖x-y‖ and
// ℓ = LengthScale. Sample paths are twice differentiable.
type Matern52 struct {
	Variance    float64
	LengthScale float64
}

// NewMatern52 returns a [Matern52] kernel with the given variance and length
// scale.
func NewMatern52(variance, lengthScale float64) Matern52 {
	return Matern52{Variance: variance, LengthScale: lengthScale}
}

// Eval returns the Matérn-5/2 covariance of x and y.
func (k Matern52) Eval(x, y []float64) float64 {
	r := Distance(x, y)
	a := math.Sqrt(5) * r / k.LengthScale
	return k.Variance * (1 + a + a*a/3) * math.Exp(-a)
}

// String returns a human-readable description of the kernel.
func (k Matern52) String() string {
	return fmt.Sprintf("Matern52(var=%g, ls=%g)", k.Variance, k.LengthScale)
}

// MaternKernel returns the [Kernel] of the Matérn family for the smoothness
// parameter nu, which must be one of 0.5, 1.5 or 2.5. It panics for any other
// value.
func MaternKernel(nu, variance, lengthScale float64) Kernel {
	switch nu {
	case 0.5:
		return NewMatern12(variance, lengthScale)
	case 1.5:
		return NewMatern32(variance, lengthScale)
	case 2.5:
		return NewMatern52(variance, lengthScale)
	default:
		panic(fmt.Sprintf("gaussproc: unsupported Matern smoothness nu=%g", nu))
	}
}

// RationalQuadratic is the kernel
// k(x, y) = Variance·(1+r²/(2·Alpha·ℓ²))^(-Alpha) with r = ‖x-y‖ and ℓ =
// LengthScale. It is a scale mixture of [RBF] kernels; as Alpha→∞ it tends to
// the [RBF] kernel.
type RationalQuadratic struct {
	Variance    float64
	LengthScale float64
	Alpha       float64
}

// NewRationalQuadratic returns a [RationalQuadratic] kernel with the given
// variance, length scale and shape parameter alpha.
func NewRationalQuadratic(variance, lengthScale, alpha float64) RationalQuadratic {
	return RationalQuadratic{Variance: variance, LengthScale: lengthScale, Alpha: alpha}
}

// Eval returns the rational-quadratic covariance of x and y.
func (k RationalQuadratic) Eval(x, y []float64) float64 {
	d2 := SquaredDistance(x, y)
	base := 1 + d2/(2*k.Alpha*k.LengthScale*k.LengthScale)
	return k.Variance * math.Pow(base, -k.Alpha)
}

// String returns a human-readable description of the kernel.
func (k RationalQuadratic) String() string {
	return fmt.Sprintf("RationalQuadratic(var=%g, ls=%g, alpha=%g)", k.Variance, k.LengthScale, k.Alpha)
}

// Periodic is the exp-sine-squared kernel
// k(x, y) = Variance·exp(-2·sin²(π·r/Period)/ℓ²) with r = ‖x-y‖ and ℓ =
// LengthScale. It models functions that repeat with the given period.
type Periodic struct {
	Variance    float64
	LengthScale float64
	Period      float64
}

// NewPeriodic returns a [Periodic] kernel with the given variance, length scale
// and period.
func NewPeriodic(variance, lengthScale, period float64) Periodic {
	return Periodic{Variance: variance, LengthScale: lengthScale, Period: period}
}

// Eval returns the periodic covariance of x and y.
func (k Periodic) Eval(x, y []float64) float64 {
	r := Distance(x, y)
	s := math.Sin(math.Pi * r / k.Period)
	return k.Variance * math.Exp(-2*s*s/(k.LengthScale*k.LengthScale))
}

// String returns a human-readable description of the kernel.
func (k Periodic) String() string {
	return fmt.Sprintf("Periodic(var=%g, ls=%g, period=%g)", k.Variance, k.LengthScale, k.Period)
}

// Cosine is the kernel k(x, y) = Variance·cos(2π·r/Period) with r = ‖x-y‖. It
// is not strictly positive definite for all data but is useful as a component
// of composite kernels.
type Cosine struct {
	Variance float64
	Period   float64
}

// NewCosine returns a [Cosine] kernel with the given variance and period.
func NewCosine(variance, period float64) Cosine {
	return Cosine{Variance: variance, Period: period}
}

// Eval returns the cosine covariance of x and y.
func (k Cosine) Eval(x, y []float64) float64 {
	r := Distance(x, y)
	return k.Variance * math.Cos(2*math.Pi*r/k.Period)
}

// String returns a human-readable description of the kernel.
func (k Cosine) String() string {
	return fmt.Sprintf("Cosine(var=%g, period=%g)", k.Variance, k.Period)
}

// GammaExponential is the kernel
// k(x, y) = Variance·exp(-(r/ℓ)^Gamma) with r = ‖x-y‖ and ℓ = LengthScale.
// Gamma must lie in (0, 2]; Gamma = 2 recovers the [RBF] kernel and Gamma = 1
// recovers [Matern12].
type GammaExponential struct {
	Variance    float64
	LengthScale float64
	Gamma       float64
}

// NewGammaExponential returns a [GammaExponential] kernel with the given
// variance, length scale and exponent gamma.
func NewGammaExponential(variance, lengthScale, gamma float64) GammaExponential {
	return GammaExponential{Variance: variance, LengthScale: lengthScale, Gamma: gamma}
}

// Eval returns the gamma-exponential covariance of x and y.
func (k GammaExponential) Eval(x, y []float64) float64 {
	r := Distance(x, y)
	return k.Variance * math.Exp(-math.Pow(r/k.LengthScale, k.Gamma))
}

// String returns a human-readable description of the kernel.
func (k GammaExponential) String() string {
	return fmt.Sprintf("GammaExponential(var=%g, ls=%g, gamma=%g)", k.Variance, k.LengthScale, k.Gamma)
}

// Linear is the non-stationary dot-product kernel
// k(x, y) = Bias + Variance·(x-c)·(y-c) with c = Offset broadcast to every
// coordinate. Gaussian processes with this kernel are Bayesian linear
// regression models.
type Linear struct {
	Bias     float64
	Variance float64
	Offset   float64
}

// NewLinear returns a [Linear] kernel with the given bias, slope variance and
// offset.
func NewLinear(bias, variance, offset float64) Linear {
	return Linear{Bias: bias, Variance: variance, Offset: offset}
}

// Eval returns the linear covariance of x and y.
func (k Linear) Eval(x, y []float64) float64 {
	if len(x) != len(y) {
		panic(ErrDimensionMismatch)
	}
	var d float64
	for i := range x {
		d += (x[i] - k.Offset) * (y[i] - k.Offset)
	}
	return k.Bias + k.Variance*d
}

// String returns a human-readable description of the kernel.
func (k Linear) String() string {
	return fmt.Sprintf("Linear(bias=%g, var=%g, offset=%g)", k.Bias, k.Variance, k.Offset)
}

// Polynomial is the kernel k(x, y) = (Alpha·x·y + Coef0)^Degree. Degree must be
// a positive integer.
type Polynomial struct {
	Alpha  float64
	Coef0  float64
	Degree int
}

// NewPolynomial returns a [Polynomial] kernel with the given scale, constant
// offset and degree.
func NewPolynomial(alpha, coef0 float64, degree int) Polynomial {
	return Polynomial{Alpha: alpha, Coef0: coef0, Degree: degree}
}

// Eval returns the polynomial covariance of x and y.
func (k Polynomial) Eval(x, y []float64) float64 {
	base := k.Alpha*Dot(x, y) + k.Coef0
	return math.Pow(base, float64(k.Degree))
}

// String returns a human-readable description of the kernel.
func (k Polynomial) String() string {
	return fmt.Sprintf("Polynomial(alpha=%g, coef0=%g, degree=%d)", k.Alpha, k.Coef0, k.Degree)
}

// Exponential returns the exponential (Ornstein-Uhlenbeck) kernel, which is the
// Matérn kernel with ν = 1/2.
func Exponential(variance, lengthScale float64) Matern12 {
	return NewMatern12(variance, lengthScale)
}
