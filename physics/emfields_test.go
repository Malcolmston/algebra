package physics

import (
	"errors"
	"math"
	"testing"
)

// physicsEMSink prevents the compiler from optimising away benchmark loop bodies.
var physicsEMSink float64

func TestElectricField(t *testing.T) {
	// Field of a 1 C charge at 1 m equals the Coulomb constant k_e.
	assertApprox(t, "E(1,1)=k_e", ElectricField(1, 1), CoulombConstant, 1e-3)

	cases := []struct {
		name string
		q, r float64
		want float64
		tol  float64
	}{
		{"unit", 1, 1, 8.9875517873681764e9, 1e3},
		{"q=2,r=3", 2, 3, 8.9875517873681764e9 * 2 / 9, 1e3},
		{"negative q", -1, 1, -8.9875517873681764e9, 1e3},
	}
	for _, c := range cases {
		assertApprox(t, "ElectricField/"+c.name, ElectricField(c.q, c.r), c.want, c.tol)
	}

	// Invariant: E·r == V for the same charge and distance.
	assertApprox(t, "E*r=V", ElectricField(3, 5)*5, ElectricPotential(3, 5), 1e-6)
}

func TestElectricPotential(t *testing.T) {
	cases := []struct {
		name string
		q, r float64
		want float64
		tol  float64
	}{
		{"unit", 1, 1, 8.9875517873681764e9, 1e3},
		{"q=2,r=4", 2, 4, 8.9875517873681764e9 / 2, 1e3},
		{"negative", -2, 1, -2 * 8.9875517873681764e9, 1e3},
	}
	for _, c := range cases {
		assertApprox(t, "ElectricPotential/"+c.name, ElectricPotential(c.q, c.r), c.want, c.tol)
	}
}

func TestFieldEnergyDensity(t *testing.T) {
	assertApprox(t, "E=2", FieldEnergyDensity(2), 0.5*VacuumPermittivity*4, 1e-20)
	assertApprox(t, "E=1e6", FieldEnergyDensity(1e6), 4.4270939064, 1e-9)
	// Even symmetry: same density for +E and -E.
	assertApprox(t, "symmetry", FieldEnergyDensity(-3), FieldEnergyDensity(3), 1e-20)
}

func TestMagneticFieldStraightWire(t *testing.T) {
	// Expected values are derived from the same μ0 the package exports (CODATA,
	// not the pre-2019 exact 4π×10⁻⁷), so B = μ0·I/(2·π·r) is verified to a
	// tight relative tolerance rather than the legacy 2×10⁻⁷ idealization.
	assertApprox(t, "I=1,r=1", MagneticFieldStraightWire(1, 1), VacuumPermeability/(2*math.Pi), 1e-13)
	assertApprox(t, "I=10,r=2", MagneticFieldStraightWire(10, 2), VacuumPermeability*10/(2*math.Pi*2), 1e-13)
	// Inverse scaling in r and linear in I.
	assertApprox(t, "scaling", MagneticFieldStraightWire(2, 4), MagneticFieldStraightWire(1, 1)*0.5, 1e-16)
}

func TestMagneticForceOnCharge(t *testing.T) {
	cases := []struct {
		name           string
		q, v, B, theta float64
		want           float64
	}{
		{"perp", 2, 3, 4, math.Pi / 2, 24},
		{"parallel", 2, 3, 4, 0, 0},
		{"30deg", 1, 2, 4, math.Pi / 6, 4},
	}
	for _, c := range cases {
		assertApprox(t, "MagneticForceOnCharge/"+c.name,
			MagneticForceOnCharge(c.q, c.v, c.B, c.theta), c.want, 1e-9)
	}
}

func TestLorentzForceMag(t *testing.T) {
	assertApprox(t, "perp", LorentzForceMag(2, 5, 3, 4, math.Pi/2), 34, 1e-9)
	assertApprox(t, "no motion", LorentzForceMag(2, 5, 0, 4, math.Pi/2), 10, 1e-9)
	assertApprox(t, "no field", LorentzForceMag(3, 0, 2, 0, math.Pi/2), 0, 1e-9)
}

func TestParallelPlateCapacitance(t *testing.T) {
	assertApprox(t, "vacuum", ParallelPlateCapacitance(2, 1, 1), 2*VacuumPermittivity, 1e-20)
	assertApprox(t, "sep=2", ParallelPlateCapacitance(1, 2, 1), VacuumPermittivity/2, 1e-20)
	assertApprox(t, "eps_r=2", ParallelPlateCapacitance(1, 1, 2), 2*VacuumPermittivity, 1e-20)
}

func TestInductorEnergy(t *testing.T) {
	assertApprox(t, "L=2,I=3", InductorEnergy(2, 3), 9, 1e-12)
	assertApprox(t, "zero current", InductorEnergy(5, 0), 0, 1e-12)
}

func TestSeriesResistance(t *testing.T) {
	assertApprox(t, "three", SeriesResistance(1, 2, 3), 6, 1e-12)
	assertApprox(t, "empty", SeriesResistance(), 0, 1e-12)
	assertApprox(t, "one", SeriesResistance(5), 5, 1e-12)
}

func TestParallelResistance(t *testing.T) {
	got, err := ParallelResistance(2, 2)
	if err != nil {
		t.Fatalf("ParallelResistance(2,2) unexpected error: %v", err)
	}
	assertApprox(t, "two equal", got, 1, 1e-12)

	got, err = ParallelResistance(4, 4, 4, 4)
	if err != nil {
		t.Fatalf("ParallelResistance four unexpected error: %v", err)
	}
	assertApprox(t, "four equal", got, 1, 1e-12)

	if _, err := ParallelResistance(); !errors.Is(err, errEmptyNetwork) {
		t.Errorf("ParallelResistance() error = %v, want errEmptyNetwork", err)
	}
	if _, err := ParallelResistance(1, 0, 2); !errors.Is(err, errZeroElement) {
		t.Errorf("ParallelResistance(1,0,2) error = %v, want errZeroElement", err)
	}
}

func TestSeriesCapacitance(t *testing.T) {
	got, err := SeriesCapacitance(2, 2)
	if err != nil {
		t.Fatalf("SeriesCapacitance(2,2) unexpected error: %v", err)
	}
	assertApprox(t, "two equal", got, 1, 1e-12)

	got, err = SeriesCapacitance(3, 6)
	if err != nil {
		t.Fatalf("SeriesCapacitance(3,6) unexpected error: %v", err)
	}
	assertApprox(t, "3 and 6", got, 2, 1e-12)

	if _, err := SeriesCapacitance(); !errors.Is(err, errEmptyNetwork) {
		t.Errorf("SeriesCapacitance() error = %v, want errEmptyNetwork", err)
	}
	if _, err := SeriesCapacitance(1, 0); !errors.Is(err, errZeroElement) {
		t.Errorf("SeriesCapacitance(1,0) error = %v, want errZeroElement", err)
	}
}

func TestParallelCapacitance(t *testing.T) {
	assertApprox(t, "three", ParallelCapacitance(1, 2, 3), 6, 1e-12)
	assertApprox(t, "empty", ParallelCapacitance(), 0, 1e-12)
}

func TestRCTimeConstant(t *testing.T) {
	assertApprox(t, "1k,1uF", RCTimeConstant(1000, 1e-6), 1e-3, 1e-15)
}

func TestRCTransients(t *testing.T) {
	const vs, r, c = 10.0, 1000.0, 1e-6 // RC = 1e-3 s

	assertApprox(t, "charge t=0", RCChargeVoltage(vs, r, c, 0), 0, 1e-12)
	assertApprox(t, "charge t=RC", RCChargeVoltage(vs, r, c, 1e-3), 6.321205588285577, 1e-9)

	assertApprox(t, "discharge t=0", RCDischargeVoltage(vs, r, c, 0), 10, 1e-12)
	assertApprox(t, "discharge t=RC", RCDischargeVoltage(vs, r, c, 1e-3), 3.678794411714423, 1e-9)

	// Complementarity: charging + discharging voltages sum to the source.
	assertApprox(t, "sum=Vs",
		RCChargeVoltage(vs, r, c, 1e-3)+RCDischargeVoltage(vs, r, c, 1e-3), vs, 1e-9)
}

func TestResonantFrequency(t *testing.T) {
	// L = 1 H, C = 1/(4π²) F gives f = 1 Hz exactly.
	assertApprox(t, "f=1", ResonantFrequency(1, 1/(4*math.Pi*math.Pi)), 1, 1e-9)
	assertApprox(t, "1mH,1uF", ResonantFrequency(1e-3, 1e-6), 5032.921210448704, 1e-6)
}

// --- Benchmarks for the performance-sensitive functions ---

func BenchmarkElectricField(b *testing.B) {
	var s float64
	for i := 0; i < b.N; i++ {
		s += ElectricField(1.6e-19, 5.29e-11)
	}
	physicsEMSink = s
}

func BenchmarkElectricPotential(b *testing.B) {
	var s float64
	for i := 0; i < b.N; i++ {
		s += ElectricPotential(1.6e-19, 5.29e-11)
	}
	physicsEMSink = s
}

func BenchmarkMagneticFieldStraightWire(b *testing.B) {
	var s float64
	for i := 0; i < b.N; i++ {
		s += MagneticFieldStraightWire(10, 0.05)
	}
	physicsEMSink = s
}

func BenchmarkSeriesResistance(b *testing.B) {
	var s float64
	for i := 0; i < b.N; i++ {
		s += SeriesResistance(1, 2, 3, 4, 5, 6, 7, 8)
	}
	physicsEMSink = s
}

func BenchmarkParallelResistance(b *testing.B) {
	var s float64
	for i := 0; i < b.N; i++ {
		r, _ := ParallelResistance(1, 2, 3, 4, 5, 6, 7, 8)
		s += r
	}
	physicsEMSink = s
}

func BenchmarkSeriesCapacitance(b *testing.B) {
	var s float64
	for i := 0; i < b.N; i++ {
		c, _ := SeriesCapacitance(1, 2, 3, 4, 5, 6, 7, 8)
		s += c
	}
	physicsEMSink = s
}

func BenchmarkParallelCapacitance(b *testing.B) {
	var s float64
	for i := 0; i < b.N; i++ {
		s += ParallelCapacitance(1, 2, 3, 4, 5, 6, 7, 8)
	}
	physicsEMSink = s
}
