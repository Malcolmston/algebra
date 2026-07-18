package physics

import "testing"

func TestCoulombConstant(t *testing.T) {
	assertApprox(t, "CoulombConstant", CoulombConstant, 8.9875517873681764e9, 1e-9)
	// Consistent with the existing CoulombForce at unit charges and distance.
	assertApprox(t, "k vs CoulombForce", CoulombConstant, CoulombForce(1, 1, 1), 1e-12)
}

func TestElectricFieldPotential(t *testing.T) {
	// Field of a 1 C charge at 1 m equals the Coulomb constant.
	assertApprox(t, "E field", ElectricFieldPointCharge(1, 1), CoulombConstant, 1e-12)
	assertApprox(t, "V potential", ElectricPotentialPointCharge(1, 1), CoulombConstant, 1e-12)

	// Force on a test charge equals field times charge: F = E·q2.
	e := ElectricFieldPointCharge(2, 3)
	assertApprox(t, "F = E·q", e*5, CoulombForce(2, 5, 3), 1e-9)
}

func TestCircuits(t *testing.T) {
	assertApprox(t, "PowerDissipated", PowerDissipated(2, 3), 12, 1e-12)

	assertApprox(t, "ResistorsSeries", ResistorsSeries(1, 2, 3), 6, 1e-12)
	assertApprox(t, "ResistorsParallel equal", ResistorsParallel(2, 2), 1, 1e-12)
	assertApprox(t, "ResistorsParallel three", ResistorsParallel(1, 1, 1), 1.0/3.0, 1e-12)
	if ResistorsSeries() != 0 || ResistorsParallel() != 0 {
		t.Error("empty resistor networks should be 0")
	}
	if ResistorsParallel(0, 5) != 0 {
		t.Error("a zero resistance should short the parallel network")
	}
}

func TestCapacitors(t *testing.T) {
	assertApprox(t, "CapacitorEnergy", CapacitorEnergy(2, 3), 9, 1e-12)
	assertApprox(t, "CapacitorCharge", CapacitorCharge(2, 3), 6, 1e-12)
	assertApprox(t, "CapacitorsSeries", CapacitorsSeries(2, 2), 1, 1e-12)
	assertApprox(t, "CapacitorsParallel", CapacitorsParallel(1, 2, 3), 6, 1e-12)
	if CapacitorsSeries() != 0 || CapacitorsParallel() != 0 {
		t.Error("empty capacitor networks should be 0")
	}
}

func BenchmarkElectricFieldPointCharge(b *testing.B) {
	var acc float64
	for i := 0; i < b.N; i++ {
		acc += ElectricFieldPointCharge(1.6e-19, 5.29e-11)
	}
	_ = acc
}
