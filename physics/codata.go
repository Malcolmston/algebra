package physics

// Extended CODATA / derived physical constants.
//
// This file augments the core catalogue in constants.go with a second tier of
// physical constants in SI units, together with O(1) lookup helpers. It adds
// only new identifiers and never modifies the base table.
//
// Each value is annotated as either "exact" — fixed by the 2019 SI
// redefinition or derived by closed form from constants that are themselves
// exact — or "measured" — an experimentally determined CODATA recommended
// value carrying the uncertainty of that recommendation.
const (
	// MolarMassConstant is the molar mass constant M_u, in kilograms per mole
	// (kg/mol). Measured (CODATA): since the 2019 SI redefinition it is no
	// longer exactly 1 g/mol.
	MolarMassConstant = 0.99999999965e-3

	// ClassicalElectronRadius is the classical electron radius r_e, in metres
	// (m). Measured (CODATA).
	ClassicalElectronRadius = 2.8179403262e-15

	// ComptonWavelength is the Compton wavelength of the electron λ_C, in
	// metres (m). Measured (CODATA).
	ComptonWavelength = 2.42631023867e-12

	// HartreeEnergy is the Hartree energy E_h, the atomic unit of energy, in
	// joules (J). Measured (CODATA).
	HartreeEnergy = 4.3597447222071e-18

	// MagneticFluxQuantum is the magnetic flux quantum Φ0 = h/(2e), in webers
	// (Wb). Exact: derived from the exact h and e.
	MagneticFluxQuantum = 2.067833848e-15

	// ConductanceQuantum is the conductance quantum G0 = 2e²/h, in siemens (S).
	// Exact: derived from the exact h and e.
	ConductanceQuantum = 7.748091729e-5

	// JosephsonConstant is the Josephson constant K_J = 2e/h, in hertz per volt
	// (Hz/V). Exact: derived from the exact h and e.
	JosephsonConstant = 483597.8484e9

	// VonKlitzingConstant is the von Klitzing constant R_K = h/e², in ohms (Ω).
	// Exact: derived from the exact h and e.
	VonKlitzingConstant = 25812.80745

	// WienDisplacement is the Wien wavelength displacement-law constant b, in
	// metre-kelvin (m·K). Exact: derived from h, c and k_B (transcendental
	// root, quoted truncated).
	WienDisplacement = 2.897771955e-3

	// FirstRadiation is the first radiation constant c_1 = 2πhc², in
	// watt square-metres (W·m²). Exact: derived from the exact h and c.
	FirstRadiation = 3.741771852e-16

	// SecondRadiation is the second radiation constant c_2 = hc/k_B, in
	// metre-kelvin (m·K). Exact: derived from the exact h, c and k_B.
	SecondRadiation = 1.438776877e-2

	// ThomsonCrossSection is the Thomson cross section σ_e, in square metres
	// (m²). Measured (CODATA).
	ThomsonCrossSection = 6.6524587321e-29

	// NuclearMagneton is the nuclear magneton μ_N, in joules per tesla (J/T).
	// Measured (CODATA).
	NuclearMagneton = 5.0507837461e-27

	// BohrMagneton is the Bohr magneton μ_B, in joules per tesla (J/T).
	// Measured (CODATA).
	BohrMagneton = 9.2740100783e-24

	// MuonMass is the rest mass of the muon m_μ, in kilograms (kg). Measured
	// (CODATA).
	MuonMass = 1.883531627e-28

	// TauMass is the rest mass of the tau lepton m_τ, in kilograms (kg).
	// Measured (CODATA).
	TauMass = 3.16754e-27

	// PlanckLength is the Planck length l_P = √(ħG/c³), in metres (m).
	// Measured (CODATA): its uncertainty is dominated by G.
	PlanckLength = 1.616255e-35

	// PlanckTime is the Planck time t_P = √(ħG/c⁵), in seconds (s).
	// Measured (CODATA): its uncertainty is dominated by G.
	PlanckTime = 5.391247e-44

	// PlanckMass is the Planck mass m_P = √(ħc/G), in kilograms (kg).
	// Measured (CODATA): its uncertainty is dominated by G.
	PlanckMass = 2.176434e-8

	// PlanckTemperature is the Planck temperature T_P = √(ħc⁵/G)/k_B, in
	// kelvin (K). Measured (CODATA): its uncertainty is dominated by G.
	PlanckTemperature = 1.416784e32

	// StandardAtmosphere is the standard atmosphere atm, in pascals (Pa).
	// Exact: conventional value defined as 101325 Pa.
	StandardAtmosphere = 101325.0

	// MolarVolumeSTP is the molar volume of an ideal gas V_m at STP
	// (T = 273.15 K, p = 101.325 kPa), in cubic metres per mole (m³/mol).
	// Exact: derived as RT/p from the exact molar gas constant.
	MolarVolumeSTP = 22.41396954e-3
)

// extendedConstants is the canonical table of the additional constants defined
// in this file, in a stable literal order. It backs [ExtendedConstants] and,
// together with the base table exposed by [Constants], the O(1) lookup indexes.
var extendedConstants = []Constant{
	{"Molar mass constant", "M_u", MolarMassConstant, "kg/mol", "Measured (CODATA); no longer exactly 1 g/mol since 2019."},
	{"Classical electron radius", "r_e", ClassicalElectronRadius, "m", "Measured (CODATA)."},
	{"Compton wavelength", "λ_C", ComptonWavelength, "m", "Electron Compton wavelength; measured (CODATA)."},
	{"Hartree energy", "E_h", HartreeEnergy, "J", "Atomic unit of energy; measured (CODATA)."},
	{"Magnetic flux quantum", "Φ0", MagneticFluxQuantum, "Wb", "h/(2e); exact (derived from exact h, e)."},
	{"Conductance quantum", "G0", ConductanceQuantum, "S", "2e²/h; exact (derived from exact h, e)."},
	{"Josephson constant", "K_J", JosephsonConstant, "Hz/V", "2e/h; exact (derived from exact h, e)."},
	{"von Klitzing constant", "R_K", VonKlitzingConstant, "Ω", "h/e²; exact (derived from exact h, e)."},
	{"Wien displacement constant", "b", WienDisplacement, "m·K", "Exact; derived from h, c, k_B."},
	{"First radiation constant", "c_1", FirstRadiation, "W·m²", "2πhc²; exact (derived from exact h, c)."},
	{"Second radiation constant", "c_2", SecondRadiation, "m·K", "hc/k_B; exact (derived from exact h, c, k_B)."},
	{"Thomson cross section", "σ_e", ThomsonCrossSection, "m²", "Measured (CODATA)."},
	{"Nuclear magneton", "μ_N", NuclearMagneton, "J/T", "Measured (CODATA)."},
	{"Bohr magneton", "μ_B", BohrMagneton, "J/T", "Measured (CODATA)."},
	{"Muon mass", "m_μ", MuonMass, "kg", "Muon rest mass; measured (CODATA)."},
	{"Tau mass", "m_τ", TauMass, "kg", "Tau lepton rest mass; measured (CODATA)."},
	{"Planck length", "l_P", PlanckLength, "m", "√(ħG/c³); measured (CODATA), limited by G."},
	{"Planck time", "t_P", PlanckTime, "s", "√(ħG/c⁵); measured (CODATA), limited by G."},
	{"Planck mass", "m_P", PlanckMass, "kg", "√(ħc/G); measured (CODATA), limited by G."},
	{"Planck temperature", "T_P", PlanckTemperature, "K", "√(ħc⁵/G)/k_B; measured (CODATA), limited by G."},
	{"Standard atmosphere", "atm", StandardAtmosphere, "Pa", "Conventional value; exact by definition."},
	{"Molar volume (STP)", "V_m", MolarVolumeSTP, "m³/mol", "Ideal gas at 273.15 K, 101.325 kPa; exact (RT/p)."},
}

// bySymbol indexes every constant — base ([Constants]) and extended
// ([extendedConstants]) — by its Symbol, giving O(1) symbol lookups. It is
// populated once at package initialisation.
var bySymbol map[string]Constant

// byName indexes every constant — base ([Constants]) and extended
// ([extendedConstants]) — by its Name, giving O(1) name lookups. It is
// populated once at package initialisation.
var byName map[string]Constant

// physicsIndexConstants builds the bySymbol and byName indexes once, from the
// base records reachable through the exported [Constants] accessor and from the
// extended records. Base records are indexed first so that the extended table
// can supplement but never silently shadow an existing base entry on the rare
// event of a shared key.
func physicsIndexConstants() {
	base := Constants()
	bySymbol = make(map[string]Constant, len(base)+len(extendedConstants))
	byName = make(map[string]Constant, len(base)+len(extendedConstants))
	for _, c := range base {
		if _, ok := bySymbol[c.Symbol]; !ok {
			bySymbol[c.Symbol] = c
		}
		if _, ok := byName[c.Name]; !ok {
			byName[c.Name] = c
		}
	}
	for _, c := range extendedConstants {
		if _, ok := bySymbol[c.Symbol]; !ok {
			bySymbol[c.Symbol] = c
		}
		if _, ok := byName[c.Name]; !ok {
			byName[c.Name] = c
		}
	}
}

func init() {
	physicsIndexConstants()
}

// ExtendedConstants returns a copy of the table of additional physical
// constants defined in this file, in a stable order. The returned slice is
// safe to modify; it does not alias the internal table.
func ExtendedConstants() []Constant {
	out := make([]Constant, len(extendedConstants))
	copy(out, extendedConstants)
	return out
}

// LookupByName returns the [Constant] whose Name matches name exactly, and true
// if found, searching both the base and the extended tables. If no constant has
// that name it returns the zero Constant and false. The lookup is O(1).
func LookupByName(name string) (Constant, bool) {
	c, ok := byName[name]
	return c, ok
}

// LookupAny returns the [Constant] whose Symbol matches symbol exactly, and
// true if found, searching both the base and the extended tables. Unlike the
// O(n) [Lookup], which scans only the base table, this is an O(1) index lookup
// that also spans the extended set. If no constant has that symbol it returns
// the zero Constant and false.
func LookupAny(symbol string) (Constant, bool) {
	c, ok := bySymbol[symbol]
	return c, ok
}
