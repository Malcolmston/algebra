package physics

import (
	"math"
	"testing"
)

// approx reports whether got and want agree to a relative tolerance of tol,
// falling back to absolute comparison when want is zero.
func approx(got, want, tol float64) bool {
	if want == 0 {
		return math.Abs(got) <= tol
	}
	return math.Abs(got-want)/math.Abs(want) <= tol
}

func assertApprox(t *testing.T, name string, got, want, tol float64) {
	t.Helper()
	if !approx(got, want, tol) {
		t.Errorf("%s = %v, want %v (tol %v)", name, got, want, tol)
	}
}

func TestConstantsValues(t *testing.T) {
	// Exact SI-2019 defining constants must match to the bit.
	if SpeedOfLight != 299792458.0 {
		t.Errorf("SpeedOfLight = %v, want 299792458", SpeedOfLight)
	}
	if PlanckConstant != 6.62607015e-34 {
		t.Errorf("PlanckConstant = %v", PlanckConstant)
	}
	if BoltzmannConstant != 1.380649e-23 {
		t.Errorf("BoltzmannConstant = %v", BoltzmannConstant)
	}
	if ElementaryCharge != 1.602176634e-19 {
		t.Errorf("ElementaryCharge = %v", ElementaryCharge)
	}
	if StandardGravity != 9.80665 {
		t.Errorf("StandardGravity = %v", StandardGravity)
	}

	// Derived exact relations.
	assertApprox(t, "ħ = h/2π", ReducedPlanck, 1.054571817e-34, 1e-9)
	assertApprox(t, "R = N_A·k_B", AvogadroConstant*BoltzmannConstant, GasConstant, 1e-9)
	assertApprox(t, "F = N_A·e", AvogadroConstant*ElementaryCharge, FaradayConstant, 1e-9)

	// A few CODATA values sanity-checked against their known magnitudes.
	assertApprox(t, "G", GravitationalConstant, 6.67430e-11, 1e-6)
	assertApprox(t, "m_e", ElectronMass, 9.1093837015e-31, 1e-9)
	assertApprox(t, "m_p", ProtonMass, 1.67262192369e-27, 1e-9)
	assertApprox(t, "α", FineStructureConstant, 1.0/137.035999, 1e-6)
}

func TestConstantsTable(t *testing.T) {
	cs := Constants()
	if len(cs) < 21 {
		t.Fatalf("Constants() returned %d entries, want at least 21", len(cs))
	}
	// The returned slice must be a copy, not an alias.
	cs[0].Value = -1
	if got := Constants()[0].Value; got == -1 {
		t.Error("Constants() aliases internal table; mutation leaked")
	}

	c, ok := Lookup("c")
	if !ok {
		t.Fatal("Lookup(c) not found")
	}
	if c.Value != SpeedOfLight || c.Unit != "m/s" {
		t.Errorf("Lookup(c) = %+v", c)
	}
	if _, ok := Lookup("nonsense"); ok {
		t.Error("Lookup(nonsense) should not be found")
	}
	// Every table value must equal its corresponding exported constant where
	// symbols are unambiguous.
	if g, _ := Lookup("g"); g.Value != StandardGravity {
		t.Errorf("Lookup(g).Value = %v, want %v", g.Value, StandardGravity)
	}
}

func TestConvertLength(t *testing.T) {
	got, err := Convert(1, "mi", "m")
	if err != nil {
		t.Fatal(err)
	}
	assertApprox(t, "1 mi -> m", got, 1609.344, 1e-12)

	got, _ = Convert(1000, "m", "km")
	assertApprox(t, "1000 m -> km", got, 1, 1e-12)

	got, _ = Convert(1, "Å", "nm")
	assertApprox(t, "1 Å -> nm", got, 0.1, 1e-12)

	got, _ = Convert(12, "in", "ft")
	assertApprox(t, "12 in -> ft", got, 1, 1e-12)
}

func TestConvertMassTimeEnergyAngle(t *testing.T) {
	got, _ := Convert(1, "lb", "kg")
	assertApprox(t, "1 lb -> kg", got, 0.45359237, 1e-12)

	got, _ = Convert(1, "hr", "s")
	assertApprox(t, "1 hr -> s", got, 3600, 1e-12)

	got, _ = Convert(1, "kWh", "J")
	assertApprox(t, "1 kWh -> J", got, 3.6e6, 1e-12)

	got, _ = Convert(1, "eV", "J")
	assertApprox(t, "1 eV -> J", got, ElectronVolt, 1e-12)

	got, _ = Convert(180, "deg", "rad")
	assertApprox(t, "180 deg -> rad", got, math.Pi, 1e-12)
}

func TestConvertTemperature(t *testing.T) {
	got, _ := Convert(0, "C", "K")
	assertApprox(t, "0 C -> K", got, 273.15, 1e-12)

	got, _ = Convert(100, "C", "F")
	assertApprox(t, "100 C -> F", got, 212, 1e-12)

	got, _ = Convert(32, "F", "C")
	assertApprox(t, "32 F -> C", got, 0, 1e-12)

	got, _ = Convert(-40, "C", "F")
	assertApprox(t, "-40 C -> F", got, -40, 1e-12)

	got, _ = Convert(273.15, "K", "C")
	assertApprox(t, "273.15 K -> C", got, 0, 1e-12)
}

func TestConvertErrors(t *testing.T) {
	if _, err := Convert(1, "m", "kg"); err == nil {
		t.Error("expected error for incompatible dimensions m -> kg")
	}
	if _, err := Convert(1, "furlong", "m"); err == nil {
		t.Error("expected error for unknown from-unit")
	}
	if _, err := Convert(1, "m", "furlong"); err == nil {
		t.Error("expected error for unknown to-unit")
	}
}

func TestTypedConversions(t *testing.T) {
	assertApprox(t, "CelsiusToKelvin(0)", CelsiusToKelvin(0), 273.15, 1e-12)
	assertApprox(t, "KelvinToCelsius(273.15)", KelvinToCelsius(273.15), 0, 1e-12)
	assertApprox(t, "CelsiusToFahrenheit(100)", CelsiusToFahrenheit(100), 212, 1e-12)
	assertApprox(t, "FahrenheitToCelsius(212)", FahrenheitToCelsius(212), 100, 1e-12)
	assertApprox(t, "KelvinToFahrenheit(273.15)", KelvinToFahrenheit(273.15), 32, 1e-12)
	assertApprox(t, "FahrenheitToKelvin(32)", FahrenheitToKelvin(32), 273.15, 1e-12)
	assertApprox(t, "DegToRad(180)", DegToRad(180), math.Pi, 1e-12)
	assertApprox(t, "RadToDeg(π)", RadToDeg(math.Pi), 180, 1e-12)
	assertApprox(t, "EVToJoules(1)", EVToJoules(1), ElectronVolt, 1e-12)
	assertApprox(t, "JoulesToEV(e)", JoulesToEV(ElectronVolt), 1, 1e-12)
	assertApprox(t, "MilesToMeters(1)", MilesToMeters(1), 1609.344, 1e-12)
	assertApprox(t, "MetersToMiles(1609.344)", MetersToMiles(1609.344), 1, 1e-12)
	assertApprox(t, "PoundsToKilograms(1)", PoundsToKilograms(1), 0.45359237, 1e-12)
	assertApprox(t, "KilogramsToPounds(0.45359237)", KilogramsToPounds(0.45359237), 1, 1e-12)
}

func TestKinematics(t *testing.T) {
	assertApprox(t, "KineticEnergy(2,3)", KineticEnergy(2, 3), 9, 1e-12) // ½·2·9
	assertApprox(t, "PotentialEnergyGravity(1,10)", PotentialEnergyGravity(1, 10), 98.0665, 1e-12)
	assertApprox(t, "Force(2,3)", Force(2, 3), 6, 1e-12)
	assertApprox(t, "Momentum(2,3)", Momentum(2, 3), 6, 1e-12)
	assertApprox(t, "Work(10,2)", Work(10, 2), 20, 1e-12)
	assertApprox(t, "Power(100,4)", Power(100, 4), 25, 1e-12)
	assertApprox(t, "FreeFallDistance(2)", FreeFallDistance(2), 2*StandardGravity, 1e-12) // ½·g·4 = 2g
	// 45° launch gives maximum range v²/g.
	assertApprox(t, "ProjectileRange(10,45°)", ProjectileRange(10, math.Pi/4), 100/StandardGravity, 1e-12)
	// Straight up gives zero range.
	assertApprox(t, "ProjectileRange(10,90°)", ProjectileRange(10, math.Pi/2), 0, 1e-9)
}

func TestWavesOptics(t *testing.T) {
	assertApprox(t, "WaveSpeed(2,3)", WaveSpeed(2, 3), 6, 1e-12)
	// A wave at wavelength c travelling at 1 Hz moves at c.
	assertApprox(t, "PhotonEnergy(1e15)", PhotonEnergy(1e15), PlanckConstant*1e15, 1e-12)
	// Green light ~5.45e14 Hz ~ 3.6e-19 J ~ 2.25 eV.
	assertApprox(t, "PhotonEnergy green (eV)", JoulesToEV(PhotonEnergy(5.45e14)), 2.254, 1e-2)
	assertApprox(t, "DeBroglieWavelength", DeBroglieWavelength(2, 3), PlanckConstant/6, 1e-12)
}

func TestRelativity(t *testing.T) {
	assertApprox(t, "LorentzFactor(0)", LorentzFactor(0), 1, 1e-15)
	// At v = 0.6c, γ = 1.25.
	assertApprox(t, "LorentzFactor(0.6c)", LorentzFactor(0.6*SpeedOfLight), 1.25, 1e-12)
	// At v = 0.8c, γ = 5/3.
	assertApprox(t, "LorentzFactor(0.8c)", LorentzFactor(0.8*SpeedOfLight), 5.0/3.0, 1e-12)
	assertApprox(t, "MassEnergy(1)", MassEnergy(1), SpeedOfLight*SpeedOfLight, 1e-15)
	assertApprox(t, "MassEnergy(1) value", MassEnergy(1), 8.987551787368176e16, 1e-12)
}

func TestThermoEM(t *testing.T) {
	// 1 mol at 273.15 K in 22.414e-3 m³ ~ 101325 Pa (approx STP).
	p := IdealGasPressure(1, 273.15, 22.414e-3)
	assertApprox(t, "IdealGasPressure STP", p, 101325, 5e-4)

	// Two 1 C charges 1 m apart repel with Coulomb's constant ~8.9876e9 N.
	f := CoulombForce(1, 1, 1)
	assertApprox(t, "CoulombForce(1,1,1)", f, 8.9875e9, 1e-3)
	// Opposite charges attract (negative).
	if CoulombForce(1, -1, 1) >= 0 {
		t.Error("CoulombForce for opposite charges should be negative (attractive)")
	}

	assertApprox(t, "VoltageOhm(2,3)", VoltageOhm(2, 3), 6, 1e-12)
	assertApprox(t, "CurrentOhm(6,3)", CurrentOhm(6, 3), 2, 1e-12)
	assertApprox(t, "ResistanceOhm(6,2)", ResistanceOhm(6, 2), 3, 1e-12)
}
