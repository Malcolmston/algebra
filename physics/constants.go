package physics

import "math"

// Physical constants in SI units.
//
// Values follow the 2019 SI redefinition and CODATA recommended values.
// Seven constants are exact by definition of the SI base units and are noted
// as "exact (SI 2019)" in their documentation; the remaining values are
// experimentally determined (CODATA) and carry the precision of the recommended
// value.
const (
	// SpeedOfLight is the speed of light in vacuum c, in metres per second
	// (m/s). Exact (SI 2019): it defines the metre.
	SpeedOfLight = 299792458.0

	// PlanckConstant is the Planck constant h, in joule-seconds (J·s).
	// Exact (SI 2019): it defines the kilogram.
	PlanckConstant = 6.62607015e-34

	// BoltzmannConstant is the Boltzmann constant k_B, in joules per kelvin
	// (J/K). Exact (SI 2019): it defines the kelvin.
	BoltzmannConstant = 1.380649e-23

	// AvogadroConstant is the Avogadro constant N_A, in reciprocal moles
	// (1/mol). Exact (SI 2019): it defines the mole.
	AvogadroConstant = 6.02214076e23

	// ElementaryCharge is the elementary charge e, in coulombs (C).
	// Exact (SI 2019): it defines the ampere.
	ElementaryCharge = 1.602176634e-19

	// StandardGravity is the standard acceleration due to gravity g, in metres
	// per second squared (m/s²). Exact: conventional value defined by the CGPM.
	StandardGravity = 9.80665

	// ElectronVolt is one electronvolt eV expressed in joules (J). Exact
	// (SI 2019): equal to the elementary charge times one volt.
	ElectronVolt = 1.602176634e-19

	// GravitationalConstant is the Newtonian constant of gravitation G, in
	// cubic metres per kilogram per second squared (m³·kg⁻¹·s⁻²). CODATA.
	GravitationalConstant = 6.67430e-11

	// ElectronMass is the rest mass of the electron m_e, in kilograms (kg).
	// CODATA.
	ElectronMass = 9.1093837015e-31

	// ProtonMass is the rest mass of the proton m_p, in kilograms (kg). CODATA.
	ProtonMass = 1.67262192369e-27

	// NeutronMass is the rest mass of the neutron m_n, in kilograms (kg).
	// CODATA.
	NeutronMass = 1.67492749804e-27

	// VacuumPermittivity is the electric constant (vacuum permittivity) ε0, in
	// farads per metre (F/m). CODATA.
	VacuumPermittivity = 8.8541878128e-12

	// VacuumPermeability is the magnetic constant (vacuum permeability) μ0, in
	// newtons per ampere squared (N/A², equivalently H/m). CODATA.
	VacuumPermeability = 1.25663706212e-6

	// FineStructureConstant is the dimensionless fine-structure constant α.
	// CODATA.
	FineStructureConstant = 7.2973525693e-3

	// GasConstant is the molar gas constant R, in joules per mole per kelvin
	// (J·mol⁻¹·K⁻¹). Exact (SI 2019): equal to N_A · k_B.
	GasConstant = 8.314462618

	// StefanBoltzmann is the Stefan–Boltzmann constant σ, in watts per square
	// metre per kelvin to the fourth (W·m⁻²·K⁻⁴). Exact (SI 2019): derived from
	// h, c and k_B.
	StefanBoltzmann = 5.670374419e-8

	// FaradayConstant is the Faraday constant F, in coulombs per mole (C/mol).
	// Exact (SI 2019): equal to N_A · e.
	FaradayConstant = 96485.33212

	// RydbergConstant is the Rydberg constant R∞, in reciprocal metres (1/m).
	// CODATA.
	RydbergConstant = 10973731.568160

	// BohrRadius is the Bohr radius a0, in metres (m). CODATA.
	BohrRadius = 5.29177210903e-11

	// AtomicMassUnit is the unified atomic mass unit u (dalton), in kilograms
	// (kg). CODATA.
	AtomicMassUnit = 1.66053906660e-27
)

// ReducedPlanck is the reduced Planck constant ħ = h/(2π), in joule-seconds
// (J·s). Derived from the exact Planck constant.
var ReducedPlanck = PlanckConstant / (2 * math.Pi)

// Constant is a self-describing physical constant record. The set of all
// records is enumerable through [Constants] and searchable through [Lookup],
// which makes the constant table usable for generated documentation and
// tooling.
type Constant struct {
	// Name is the human-readable name, e.g. "Speed of light".
	Name string
	// Symbol is the conventional symbol, e.g. "c".
	Symbol string
	// Value is the numeric value in SI units.
	Value float64
	// Unit is the SI unit string, e.g. "m/s". Dimensionless constants use "".
	Unit string
	// Description is a short note, including whether the value is exact.
	Description string
}

// constants is the canonical table backing Constants and Lookup.
var constants = []Constant{
	{"Speed of light", "c", SpeedOfLight, "m/s", "Speed of light in vacuum; exact (SI 2019)."},
	{"Planck constant", "h", PlanckConstant, "J·s", "Exact (SI 2019); defines the kilogram."},
	{"Reduced Planck constant", "ħ", ReducedPlanck, "J·s", "h/(2π); derived from the exact Planck constant."},
	{"Gravitational constant", "G", GravitationalConstant, "m³·kg⁻¹·s⁻²", "Newtonian constant of gravitation; CODATA."},
	{"Boltzmann constant", "k_B", BoltzmannConstant, "J/K", "Exact (SI 2019); defines the kelvin."},
	{"Avogadro constant", "N_A", AvogadroConstant, "1/mol", "Exact (SI 2019); defines the mole."},
	{"Elementary charge", "e", ElementaryCharge, "C", "Exact (SI 2019); defines the ampere."},
	{"Electron mass", "m_e", ElectronMass, "kg", "Electron rest mass; CODATA."},
	{"Proton mass", "m_p", ProtonMass, "kg", "Proton rest mass; CODATA."},
	{"Neutron mass", "m_n", NeutronMass, "kg", "Neutron rest mass; CODATA."},
	{"Vacuum permittivity", "ε0", VacuumPermittivity, "F/m", "Electric constant; CODATA."},
	{"Vacuum permeability", "μ0", VacuumPermeability, "N/A²", "Magnetic constant; CODATA."},
	{"Fine-structure constant", "α", FineStructureConstant, "", "Dimensionless; CODATA."},
	{"Molar gas constant", "R", GasConstant, "J·mol⁻¹·K⁻¹", "Exact (SI 2019); equals N_A · k_B."},
	{"Stefan–Boltzmann constant", "σ", StefanBoltzmann, "W·m⁻²·K⁻⁴", "Exact (SI 2019); derived from h, c, k_B."},
	{"Faraday constant", "F", FaradayConstant, "C/mol", "Exact (SI 2019); equals N_A · e."},
	{"Rydberg constant", "R∞", RydbergConstant, "1/m", "CODATA."},
	{"Bohr radius", "a0", BohrRadius, "m", "CODATA."},
	{"Standard gravity", "g", StandardGravity, "m/s²", "Conventional value; exact by definition."},
	{"Atomic mass unit", "u", AtomicMassUnit, "kg", "Unified atomic mass unit (dalton); CODATA."},
	{"Electronvolt", "eV", ElectronVolt, "J", "One electronvolt in joules; exact (SI 2019)."},
}

// Constants returns a copy of the table of physical constants known to the
// package, in a stable order. The returned slice is safe to modify; it does
// not alias the internal table.
func Constants() []Constant {
	out := make([]Constant, len(constants))
	copy(out, constants)
	return out
}

// Lookup returns the [Constant] whose Symbol matches symbol exactly, and true
// if found. If no constant has that symbol it returns the zero Constant and
// false.
func Lookup(symbol string) (Constant, bool) {
	for _, c := range constants {
		if c.Symbol == symbol {
			return c, true
		}
	}
	return Constant{}, false
}
