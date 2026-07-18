package physics

import (
	"math"
	"strconv"
)

// Measurement is a physical quantity paired with its standard (1σ)
// uncertainty. The arithmetic methods implement first-order (linear) error
// propagation: contributions combine in quadrature under the assumption that
// the operands are statistically INDEPENDENT (uncorrelated). Any correlation
// between operands is ignored, so results for correlated inputs are only
// approximate.
type Measurement struct {
	// Value is the central (best-estimate) value of the quantity.
	Value float64
	// Sigma is the standard (1σ) uncertainty. It is conventionally
	// non-negative; NewMeasurement enforces this for constructed values.
	Sigma float64
}

// NewMeasurement returns a Measurement with the given value and 1σ standard
// uncertainty. The uncertainty is stored as its magnitude, so a negative
// argument is treated as |sigma|.
func NewMeasurement(value, sigma float64) Measurement {
	return Measurement{value, math.Abs(sigma)}
}

// Relative returns the fractional (relative) uncertainty Sigma/|Value|. It is
// undefined when Value is 0 and follows IEEE-754 semantics there (+Inf for a
// non-zero Sigma, NaN when Sigma is also 0).
func (m Measurement) Relative() float64 {
	return m.Sigma / math.Abs(m.Value)
}

// Add returns the sum m + o. The values add and the standard uncertainties of
// the (independent) operands combine in quadrature: σ = √(σ₁² + σ₂²).
// Correlation between m and o is ignored.
func (m Measurement) Add(o Measurement) Measurement {
	return Measurement{m.Value + o.Value, physicsUncHypot(m.Sigma, o.Sigma)}
}

// Sub returns the difference m - o. The values subtract and the standard
// uncertainties of the (independent) operands combine in quadrature:
// σ = √(σ₁² + σ₂²). Correlation between m and o is ignored.
func (m Measurement) Sub(o Measurement) Measurement {
	return Measurement{m.Value - o.Value, physicsUncHypot(m.Sigma, o.Sigma)}
}

// Mul returns the product m·o. The values multiply and the RELATIVE
// uncertainties add in quadrature, so the absolute uncertainty is
// |xy|·√((σ₁/x)² + (σ₂/y)²), evaluated in the numerically robust form
// √((y·σ₁)² + (x·σ₂)²). The operands are assumed independent; correlation is
// ignored.
func (m Measurement) Mul(o Measurement) Measurement {
	return Measurement{m.Value * o.Value, physicsUncHypot(o.Value*m.Sigma, m.Value*o.Sigma)}
}

// Div returns the quotient m/o. The values divide and the RELATIVE
// uncertainties add in quadrature, giving an absolute uncertainty of
// |x/y|·√((σ₁/x)² + (σ₂/y)²). o.Value must be non-zero. The operands are
// assumed independent; correlation is ignored.
func (m Measurement) Div(o Measurement) Measurement {
	val := m.Value / o.Value
	sigma := physicsUncHypot(m.Sigma/o.Value, m.Value*o.Sigma/(o.Value*o.Value))
	return Measurement{val, sigma}
}

// Scale returns m multiplied by the exact constant k (a value carrying no
// uncertainty of its own). Both the value and the standard uncertainty scale,
// the latter by |k|.
func (m Measurement) Scale(k float64) Measurement {
	return Measurement{m.Value * k, m.Sigma * math.Abs(k)}
}

// Pow returns m raised to the exact power n. The value is Valueⁿ and the
// relative uncertainty scales by |n|, giving an absolute uncertainty of
// |Valueⁿ|·|n|·(Sigma/|Value|).
func (m Measurement) Pow(n float64) Measurement {
	val := math.Pow(m.Value, n)
	sigma := math.Abs(val) * math.Abs(n) * m.Relative()
	return Measurement{val, sigma}
}

// Apply propagates uncertainty through an arbitrary single-variable function f.
// The output value is f(Value) and the output uncertainty is |f'(Value)|·Sigma,
// where the derivative is estimated by a central finite difference with a
// deterministic, adaptively sized step. It assumes f is smooth near Value and
// that a linear approximation is adequate over the ±Sigma range.
func (m Measurement) Apply(f func(float64) float64) Measurement {
	h := physicsUncStep(m.Value)
	deriv := (f(m.Value+h) - f(m.Value-h)) / (2 * h)
	return Measurement{f(m.Value), math.Abs(deriv) * m.Sigma}
}

// Propagate propagates uncertainty through an arbitrary multi-variable function
// f evaluated at means, with per-input standard uncertainties sigmas. It
// returns a Measurement whose value is f(means) and whose uncertainty is
// √(Σ (∂f/∂xᵢ)² σᵢ²), each partial derivative estimated by a central finite
// difference. The inputs are assumed independent; correlation between them is
// ignored. It returns an error if len(means) != len(sigmas).
//
// A single scratch slice (a copy of means) is allocated once per call and
// reused across every partial-derivative evaluation, so the derivative loop
// performs no additional allocation.
func Propagate(f func([]float64) float64, means, sigmas []float64) (Measurement, error) {
	if len(means) != len(sigmas) {
		return Measurement{}, physicsErrLengthMismatch(len(means), len(sigmas))
	}
	value := f(means)
	scratch := make([]float64, len(means))
	copy(scratch, means)
	var sumSq float64
	for i := range means {
		h := physicsUncStep(means[i])
		scratch[i] = means[i] + h
		fp := f(scratch)
		scratch[i] = means[i] - h
		fm := f(scratch)
		scratch[i] = means[i] // restore before the next iteration
		deriv := (fp - fm) / (2 * h)
		contrib := deriv * sigmas[i]
		sumSq += contrib * contrib
	}
	return Measurement{value, math.Sqrt(sumSq)}, nil
}

// String formats the measurement as "value ± sigma" using the shortest decimal
// representation that round-trips each component.
func (m Measurement) String() string {
	v := strconv.FormatFloat(m.Value, 'g', -1, 64)
	s := strconv.FormatFloat(m.Sigma, 'g', -1, 64)
	return v + " ± " + s
}

// physicsUncHypot returns √(a² + b²) without the overflow guarding of
// math.Hypot, which is unnecessary for the magnitudes typical of measurement
// uncertainties and keeps propagation to pure stack arithmetic.
func physicsUncHypot(a, b float64) float64 { return math.Sqrt(a*a + b*b) }

// physicsUncStep returns a deterministic finite-difference step for evaluating
// a central derivative at x. It scales with the magnitude of x (never below 1)
// and uses the cube root of machine epsilon, which is the order that minimises
// the combined truncation and round-off error of a central difference.
func physicsUncStep(x float64) float64 {
	const cbrtEps = 6.055454452393339e-06 // math.Cbrt(0x1p-52)
	scale := math.Abs(x)
	if scale < 1 {
		scale = 1
	}
	return cbrtEps * scale
}

// physicsErrLengthMismatch builds the error returned by Propagate when the
// means and sigmas slices differ in length.
func physicsErrLengthMismatch(nMeans, nSigmas int) error {
	return &physicsLengthMismatchError{nMeans: nMeans, nSigmas: nSigmas}
}

// physicsLengthMismatchError reports mismatched means/sigmas lengths passed to
// Propagate.
type physicsLengthMismatchError struct {
	nMeans  int
	nSigmas int
}

// Error implements the error interface.
func (e *physicsLengthMismatchError) Error() string {
	return "physics: Propagate means/sigmas length mismatch: " +
		strconv.Itoa(e.nMeans) + " != " + strconv.Itoa(e.nSigmas)
}
