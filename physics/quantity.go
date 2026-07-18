package physics

import (
	"fmt"
	"math"
)

// Dim is the exponent vector of the seven SI base dimensions. The array indices
// are fixed and given by the exported basis vectors, in the order
// {Length, Mass, Time, Current, Temperature, Amount, Luminosity}. A Dim is a
// value array, so composing or comparing dimensions never allocates.
type Dim [7]int8

// Indices of the seven SI base dimensions within a [Dim] array. They are
// unexported implementation detail; use the exported basis vectors
// ([DimLength] and friends) to build dimensions.
const (
	physicsIdxLength = iota
	physicsIdxMass
	physicsIdxTime
	physicsIdxCurrent
	physicsIdxTemperature
	physicsIdxAmount
	physicsIdxLuminosity
)

// The seven SI base dimensions as unit basis vectors. Each has a single
// exponent of one in its own slot and zero elsewhere.
var (
	// DimLength is the base dimension of length (metre, m).
	DimLength = Dim{physicsIdxLength: 1}
	// DimMass is the base dimension of mass (kilogram, kg).
	DimMass = Dim{physicsIdxMass: 1}
	// DimTime is the base dimension of time (second, s).
	DimTime = Dim{physicsIdxTime: 1}
	// DimCurrent is the base dimension of electric current (ampere, A).
	DimCurrent = Dim{physicsIdxCurrent: 1}
	// DimTemperature is the base dimension of thermodynamic temperature
	// (kelvin, K).
	DimTemperature = Dim{physicsIdxTemperature: 1}
	// DimAmount is the base dimension of amount of substance (mole, mol).
	DimAmount = Dim{physicsIdxAmount: 1}
	// DimLuminosity is the base dimension of luminous intensity (candela, cd).
	DimLuminosity = Dim{physicsIdxLuminosity: 1}
)

// Common composed dimensions, expressed in the SI base dimensions.
var (
	// DimVelocity is length per time, L·T⁻¹ (m·s⁻¹).
	DimVelocity = Dim{physicsIdxLength: 1, physicsIdxTime: -1}
	// DimAccel is length per time squared, L·T⁻² (m·s⁻²).
	DimAccel = Dim{physicsIdxLength: 1, physicsIdxTime: -2}
	// DimForce is mass·length per time squared, M·L·T⁻² (the newton, N).
	DimForce = Dim{physicsIdxLength: 1, physicsIdxMass: 1, physicsIdxTime: -2}
	// DimEnergy is mass·length² per time squared, M·L²·T⁻² (the joule, J).
	DimEnergy = Dim{physicsIdxLength: 2, physicsIdxMass: 1, physicsIdxTime: -2}
	// DimPower is mass·length² per time cubed, M·L²·T⁻³ (the watt, W).
	DimPower = Dim{physicsIdxLength: 2, physicsIdxMass: 1, physicsIdxTime: -3}
	// DimPressure is mass per length per time squared, M·L⁻¹·T⁻² (the pascal, Pa).
	DimPressure = Dim{physicsIdxLength: -1, physicsIdxMass: 1, physicsIdxTime: -2}
	// DimCharge is current·time, I·T (the coulomb, C).
	DimCharge = Dim{physicsIdxTime: 1, physicsIdxCurrent: 1}
	// DimVoltage is mass·length² per time cubed per current, M·L²·T⁻³·I⁻¹
	// (the volt, V).
	DimVoltage = Dim{physicsIdxLength: 2, physicsIdxMass: 1, physicsIdxTime: -3, physicsIdxCurrent: -1}
)

// physicsPrintOrder lists the base-dimension indices and their unit symbols in
// the conventional order used to render SI derived units (kg·m·s⁻²…).
var physicsPrintOrder = [7]struct {
	idx    int
	symbol string
}{
	{physicsIdxMass, "kg"},
	{physicsIdxLength, "m"},
	{physicsIdxTime, "s"},
	{physicsIdxCurrent, "A"},
	{physicsIdxTemperature, "K"},
	{physicsIdxAmount, "mol"},
	{physicsIdxLuminosity, "cd"},
}

// physicsSuperscript renders an integer as Unicode superscript digits (with a
// leading superscript minus for negatives), e.g. -2 becomes "⁻²".
func physicsSuperscript(n int) string {
	if n == 0 {
		return "⁰" // ⁰
	}
	digits := [10]rune{'⁰', '¹', '²', '³', '⁴', '⁵', '⁶', '⁷', '⁸', '⁹'}
	neg := n < 0
	if neg {
		n = -n
	}
	var buf []rune
	for n > 0 {
		buf = append([]rune{digits[n%10]}, buf...)
		n /= 10
	}
	if neg {
		buf = append([]rune{'⁻'}, buf...) // ⁻
	}
	return string(buf)
}

// String renders the dimension as a middle-dot separated product of base-unit
// symbols with Unicode superscript exponents, for example "kg·m·s⁻²" for force.
// The dimensionless dimension renders as "1".
func (d Dim) String() string {
	out := ""
	for _, p := range physicsPrintOrder {
		e := int(d[p.idx])
		if e == 0 {
			continue
		}
		if out != "" {
			out += "·" // ·
		}
		out += p.symbol
		if e != 1 {
			out += physicsSuperscript(e)
		}
	}
	if out == "" {
		return "1"
	}
	return out
}

// Quantity is a physical quantity: a numeric value together with its
// dimension. Value is always expressed in coherent SI base units (metre,
// kilogram, second, ampere, kelvin, mole, candela), so quantities of the same
// [Dim] are directly comparable and additive.
type Quantity struct {
	// Value is the magnitude in coherent SI base units.
	Value float64
	// Dim is the physical dimension of the quantity.
	Dim Dim
}

// NewQuantity returns a Quantity with the given value (in coherent SI base
// units) and dimension.
func NewQuantity(value float64, dim Dim) Quantity {
	return Quantity{Value: value, Dim: dim}
}

// Mul returns the product of q and r: values multiply and dimension exponents
// add. It never allocates.
func (q Quantity) Mul(r Quantity) Quantity {
	var d Dim
	for i := range d {
		d[i] = q.Dim[i] + r.Dim[i]
	}
	return Quantity{Value: q.Value * r.Value, Dim: d}
}

// Div returns the quotient of q and r: values divide and dimension exponents
// subtract. It never allocates.
func (q Quantity) Div(r Quantity) Quantity {
	var d Dim
	for i := range d {
		d[i] = q.Dim[i] - r.Dim[i]
	}
	return Quantity{Value: q.Value / r.Value, Dim: d}
}

// Add returns the sum of q and r. It returns an error if the two quantities
// have different dimensions.
func (q Quantity) Add(r Quantity) (Quantity, error) {
	if q.Dim != r.Dim {
		return Quantity{}, fmt.Errorf("physics: cannot add %v and %v (different dimensions)", q.Dim, r.Dim)
	}
	return Quantity{Value: q.Value + r.Value, Dim: q.Dim}, nil
}

// Sub returns the difference q-r. It returns an error if the two quantities
// have different dimensions.
func (q Quantity) Sub(r Quantity) (Quantity, error) {
	if q.Dim != r.Dim {
		return Quantity{}, fmt.Errorf("physics: cannot subtract %v from %v (different dimensions)", r.Dim, q.Dim)
	}
	return Quantity{Value: q.Value - r.Value, Dim: q.Dim}, nil
}

// Pow returns q raised to the integer power n: the value is raised to n and
// every dimension exponent is multiplied by n. It never allocates.
func (q Quantity) Pow(n int) Quantity {
	var d Dim
	for i := range d {
		d[i] = q.Dim[i] * int8(n)
	}
	return Quantity{Value: math.Pow(q.Value, float64(n)), Dim: d}
}

// IsDimensionless reports whether q has no physical dimension (all base
// exponents zero), such as a pure number or an angle in radians.
func (q Quantity) IsDimensionless() bool {
	return q.Dim == Dim{}
}

// symbolQuantities maps every recognised unit symbol to a Quantity whose Value
// is the size of one such unit in coherent SI base units and whose Dim is the
// unit's dimension. It is built once at package load, so [ParseQuantity] and
// [Quantity.To] cost only a single map lookup plus one multiply, with no
// allocation. Affine units (Celsius, Fahrenheit) are intentionally absent
// because they cannot be represented as a pure scale factor.
var symbolQuantities = map[string]Quantity{
	// Length, base metre (m). Mirrors the length rows of the units map plus
	// the imperial yard.
	"m":  {1, DimLength},
	"km": {1000, DimLength},
	"cm": {0.01, DimLength},
	"mm": {0.001, DimLength},
	"mi": {1609.344, DimLength},
	"ft": {0.3048, DimLength},
	"in": {0.0254, DimLength},
	"yd": {0.9144, DimLength},
	"nm": {1e-9, DimLength},
	"Å":  {1e-10, DimLength},

	// Mass, base kilogram (kg).
	"kg": {1, DimMass},
	"g":  {0.001, DimMass},
	"mg": {1e-6, DimMass},
	"lb": {0.45359237, DimMass},
	"oz": {0.028349523125, DimMass},

	// Time, base second (s).
	"s":   {1, DimTime},
	"min": {60, DimTime},
	"hr":  {3600, DimTime},
	"day": {86400, DimTime},
	"yr":  {31557600, DimTime},

	// Energy, base joule (J).
	"J":     {1, DimEnergy},
	"eV":    {ElectronVolt, DimEnergy},
	"cal":   {4.184, DimEnergy},
	"kWh":   {3.6e6, DimEnergy},
	"BTU":   {1055.05585262, DimEnergy},
	"therm": {1.05505585262e8, DimEnergy},

	// Angle is dimensionless in the SI; the radian is the base unit.
	"rad": {1, Dim{}},
	"deg": {math.Pi / 180, Dim{}},

	// Volume, a derived length³ dimension.
	"gal": {0.003785411784, Dim{physicsIdxLength: 3}},

	// Force, base newton (N).
	"lbf": {4.4482216152605, DimForce},

	// Power, base watt (W). Mechanical horsepower (550 ft·lbf/s).
	"hp": {745.6998715822702, DimPower},

	// Pressure, base pascal (Pa).
	"psi": {6894.757293168361, DimPressure},
	"atm": {101325, DimPressure},
	"bar": {100000, DimPressure},

	// Velocity, base metre per second (m·s⁻¹).
	"mph":  {0.44704, DimVelocity},
	"knot": {0.5144444444444445, DimVelocity},
}

// ParseQuantity returns the Quantity corresponding to value expressed in the
// unit named by unitSymbol. The result's Value is in coherent SI base units.
// Recognised symbols extend those of [Convert] (length, mass, time, energy and
// angle) with imperial and derived units:
//
//	length:   yd
//	volume:   gal
//	force:    lbf
//	power:    hp
//	pressure: psi, atm, bar
//	velocity: mph, knot
//	energy:   BTU, therm
//
// It returns an error for an unknown symbol. Affine temperature units (C, F)
// are not supported here; use [Convert] or the temperature helpers for those.
func ParseQuantity(value float64, unitSymbol string) (Quantity, error) {
	base, ok := symbolQuantities[unitSymbol]
	if !ok {
		return Quantity{}, fmt.Errorf("physics: unknown unit %q", unitSymbol)
	}
	return Quantity{Value: value * base.Value, Dim: base.Dim}, nil
}

// To expresses q in the unit named by unitSymbol and returns the numeric value
// in that unit. It returns an error if the symbol is unknown or if its
// dimension differs from q's dimension. The operation is a single map lookup
// plus one division and allocates nothing.
func (q Quantity) To(unitSymbol string) (float64, error) {
	base, ok := symbolQuantities[unitSymbol]
	if !ok {
		return 0, fmt.Errorf("physics: unknown unit %q", unitSymbol)
	}
	if base.Dim != q.Dim {
		return 0, fmt.Errorf("physics: incompatible dimensions %v and %v (unit %q)", q.Dim, base.Dim, unitSymbol)
	}
	return q.Value / base.Value, nil
}
