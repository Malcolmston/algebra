package physics

import (
	"math"
	"strconv"
)

// Measurement is a physical quantity paired with its standard (1σ) uncertainty.
// The arithmetic methods propagate uncertainty assuming the operands are
// independent, combining contributions in quadrature per the standard first-
// order (linear) error-propagation formulas.
type Measurement struct {
	// Value is the central (best-estimate) value of the quantity.
	Value float64
	// Uncertainty is the standard uncertainty, always stored non-negative.
	Uncertainty float64
}

// NewMeasurement returns a Measurement with the given value and uncertainty.
// The uncertainty is stored as its absolute value, so a negative argument is
// treated as its magnitude.
func NewMeasurement(value, uncertainty float64) Measurement {
	return Measurement{value, math.Abs(uncertainty)}
}

// RelativeUncertainty returns the fractional (relative) uncertainty
// |Uncertainty/Value| of the measurement. It returns 0 when Value is 0, since
// the ratio is undefined there.
func (m Measurement) RelativeUncertainty() float64 {
	if m.Value == 0 {
		return 0
	}
	return math.Abs(m.Uncertainty / m.Value)
}

// Add returns the sum m + other. The values add and the absolute uncertainties
// combine in quadrature: σ = √(σ₁² + σ₂²).
func (m Measurement) Add(other Measurement) Measurement {
	return Measurement{m.Value + other.Value, physicsHypot(m.Uncertainty, other.Uncertainty)}
}

// Sub returns the difference m - other. The values subtract and the absolute
// uncertainties combine in quadrature: σ = √(σ₁² + σ₂²).
func (m Measurement) Sub(other Measurement) Measurement {
	return Measurement{m.Value - other.Value, physicsHypot(m.Uncertainty, other.Uncertainty)}
}

// Mul returns the product m·other. The values multiply and the relative
// uncertainties combine in quadrature, so the absolute uncertainty is
// |xy|·√((σ₁/x)² + (σ₂/y)²).
func (m Measurement) Mul(other Measurement) Measurement {
	val := m.Value * other.Value
	rel := physicsHypot(m.RelativeUncertainty(), other.RelativeUncertainty())
	return Measurement{val, math.Abs(val) * rel}
}

// Div returns the quotient m/other. The values divide and the relative
// uncertainties combine in quadrature. other.Value must be non-zero.
func (m Measurement) Div(other Measurement) Measurement {
	val := m.Value / other.Value
	rel := physicsHypot(m.RelativeUncertainty(), other.RelativeUncertainty())
	return Measurement{val, math.Abs(val) * rel}
}

// Scale returns the measurement multiplied by an exact constant k (a value with
// no uncertainty of its own). Both the value and the uncertainty scale by |k|.
func (m Measurement) Scale(k float64) Measurement {
	return Measurement{m.Value * k, m.Uncertainty * math.Abs(k)}
}

// Pow returns m raised to the exact power n. The value is Valueⁿ and the
// relative uncertainty scales by |n|, giving an absolute uncertainty of
// |Valueⁿ|·|n|·(σ/|Value|).
func (m Measurement) Pow(n float64) Measurement {
	val := math.Pow(m.Value, n)
	rel := math.Abs(n) * m.RelativeUncertainty()
	return Measurement{val, math.Abs(val) * rel}
}

// String formats the measurement as "value ± uncertainty" using the shortest
// decimal representation that round-trips each component.
func (m Measurement) String() string {
	v := strconv.FormatFloat(m.Value, 'g', -1, 64)
	u := strconv.FormatFloat(m.Uncertainty, 'g', -1, 64)
	return v + " ± " + u
}

// physicsHypot returns √(a² + b²) without the intermediate overflow guarding of
// math.Hypot, which is unnecessary for the typical magnitudes of measurement
// uncertainties and keeps propagation cheap.
func physicsHypot(a, b float64) float64 { return math.Sqrt(a*a + b*b) }
