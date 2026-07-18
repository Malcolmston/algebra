package physics

import "testing"

func TestIdealGas(t *testing.T) {
	// Molar volume at STP (1 atm, 0 °C) is about 22.414 L = 0.022414 m³.
	v := IdealGasVolume(1, 273.15, 101325)
	assertApprox(t, "IdealGasVolume STP", v, 0.022414, 5e-4)

	// Round-trips against the existing IdealGasPressure.
	n := IdealGasMoles(101325, v, 273.15)
	assertApprox(t, "IdealGasMoles", n, 1, 1e-9)
	temp := IdealGasTemperature(101325, v, 1)
	assertApprox(t, "IdealGasTemperature", temp, 273.15, 1e-9)
	assertApprox(t, "IdealGasPressure round-trip", IdealGasPressure(1, 273.15, v), 101325, 1e-9)
}

func TestRMSSpeed(t *testing.T) {
	// Nitrogen (M = 0.028 kg/mol) at 300 K: v_rms ≈ 517 m/s.
	assertApprox(t, "RMSSpeed N2", RMSSpeed(300, 0.028), 517, 5e-3)
}

func TestHeatTransfer(t *testing.T) {
	// Raising 1 kg of water (c = 4186 J/kg·K) by 1 K takes 4186 J.
	assertApprox(t, "HeatEnergy", HeatEnergy(1, 4186, 1), 4186, 1e-12)
	// Negative ΔT releases heat.
	assertApprox(t, "HeatEnergy release", HeatEnergy(1, 4186, -1), -4186, 1e-12)

	// Unit slab with unit conductivity and unit gradient conducts 1 W.
	assertApprox(t, "ConductionRate", ConductionRate(1, 1, 1, 1), 1, 1e-12)

	// Blackbody (ε=1) of 1 m² at 100 K radiates σ·10⁸ W.
	assertApprox(t, "RadiatedPower", RadiatedPower(1, 1, 100), StefanBoltzmann*1e8, 1e-12)

	assertApprox(t, "ThermalExpansion", ThermalExpansion(1e-5, 1, 100), 1e-3, 1e-12)
}

func BenchmarkRadiatedPower(b *testing.B) {
	var acc float64
	for i := 0; i < b.N; i++ {
		acc += RadiatedPower(0.9, 2.5, 500)
	}
	_ = acc
}

func BenchmarkRMSSpeed(b *testing.B) {
	var acc float64
	for i := 0; i < b.N; i++ {
		acc += RMSSpeed(300, 0.028)
	}
	_ = acc
}
