package physics

import "testing"

// assertApprox and approx are defined in physics_test.go and shared here.

func TestIdealGas(t *testing.T) {
	// Molar volume at STP (1 atm, 0 °C) is about 22.414 L = 0.022414 m³.
	v := IdealGasVolume(1, 273.15, 101325)
	assertApprox(t, "IdealGasVolume STP", v, 0.022414, 5e-4)

	// Round-trips against the other ideal-gas relations and IdealGasPressure.
	assertApprox(t, "IdealGasMoles", IdealGasMoles(101325, v, 273.15), 1, 1e-9)
	assertApprox(t, "IdealGasTemperature", IdealGasTemperature(101325, v, 1), 273.15, 1e-9)
	assertApprox(t, "IdealGasPressure round-trip", IdealGasPressure(1, 273.15, v), 101325, 1e-9)

	// Number density at STP is the Loschmidt constant, ≈ 2.6868e25 m⁻³.
	assertApprox(t, "NumberDensity STP", NumberDensity(101325, 273.15), 2.686780e25, 1e-3)
}

func TestKineticTheory(t *testing.T) {
	// Nitrogen (M = 0.028 kg/mol) at 300 K.
	assertApprox(t, "RMSSpeed N2", RMSSpeed(300, 0.028), 516.96, 1e-3)
	assertApprox(t, "MostProbableSpeed N2", MostProbableSpeed(300, 0.028), 422.10, 1e-3)

	// Mean translational kinetic energy per particle: 3/2·k_B·T.
	assertApprox(t, "MeanKineticEnergy 300K", MeanKineticEnergy(300), 6.2129205e-21, 1e-9)

	// v_p / v_rms = √(2/3) exactly, independent of gas and temperature.
	ratio := MostProbableSpeed(300, 0.028) / RMSSpeed(300, 0.028)
	assertApprox(t, "vp/vrms", ratio, 0.8164965809, 1e-9)
}

func TestHeat(t *testing.T) {
	// Raising 1 kg of water (c = 4186 J/kg·K) by 1 K takes 4186 J.
	assertApprox(t, "HeatCapacity", HeatCapacity(1, 4186, 1), 4186, 1e-12)
	// Negative ΔT releases heat.
	assertApprox(t, "HeatCapacity release", HeatCapacity(1, 4186, -1), -4186, 1e-12)

	// Melting 10 g of ice: L_fusion = 334 kJ/kg → 3340 J.
	assertApprox(t, "LatentHeat", LatentHeat(0.01, 334000), 3340, 1e-12)

	// Unit slab with unit conductivity and unit gradient conducts 1 W.
	assertApprox(t, "ThermalConduction unit", ThermalConduction(1, 1, 1, 1), 1, 1e-12)
	// Doubling area doubles the rate; doubling thickness halves it.
	assertApprox(t, "ThermalConduction slab", ThermalConduction(2, 4, 10, 0.5), 160, 1e-12)
}

func TestCyclesEfficiency(t *testing.T) {
	// Engine between 300 K and 600 K: Carnot efficiency 1 − 1/2 = 0.5.
	assertApprox(t, "CarnotEfficiency", CarnotEfficiency(300, 600), 0.5, 1e-12)
	assertApprox(t, "COPRefrigerator", COPRefrigerator(300, 600), 1, 1e-12)
	assertApprox(t, "COPHeatPump", COPHeatPump(300, 600), 2, 1e-12)

	// COP_heatpump = COP_refrigerator + 1 always holds.
	assertApprox(t, "COP identity", COPHeatPump(250, 300)-COPRefrigerator(250, 300), 1, 1e-12)
}

func TestRadiation(t *testing.T) {
	// Blackbody (ε=1) of 1 m² at 100 K radiates σ·10⁸ W exactly.
	assertApprox(t, "StefanBoltzmannPower", StefanBoltzmannPower(1, 1, 100), StefanBoltzmann*1e8, 1e-12)

	// Precomputed Wien displacement constant.
	assertApprox(t, "wienB", wienB, 2.897771955e-3, 1e-6)

	// Sun (T ≈ 5778 K) peaks near 501 nm.
	assertApprox(t, "WienPeakWavelength Sun", WienPeakWavelength(5778), 5.0152e-7, 1e-3)

	// Planck spectral radiance at 500 nm, 5778 K.
	assertApprox(t, "PlanckSpectralRadiance", PlanckSpectralRadiance(500e-9, 5778), 2.637e13, 1e-2)
}

func BenchmarkPlanckSpectralRadiance(b *testing.B) {
	var acc float64
	for i := 0; i < b.N; i++ {
		acc += PlanckSpectralRadiance(500e-9, 5778)
	}
	_ = acc
}

func BenchmarkStefanBoltzmannPower(b *testing.B) {
	var acc float64
	for i := 0; i < b.N; i++ {
		acc += StefanBoltzmannPower(0.9, 2.5, 500)
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

func BenchmarkMostProbableSpeed(b *testing.B) {
	var acc float64
	for i := 0; i < b.N; i++ {
		acc += MostProbableSpeed(300, 0.028)
	}
	_ = acc
}
