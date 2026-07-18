package physics

import "testing"

func TestRelativisticEnergyMomentum(t *testing.T) {
	const v = 0.6 * SpeedOfLight // γ = 1.25

	assertApprox(t, "Beta", Beta(v), 0.6, 1e-12)
	assertApprox(t, "RelMomentum", RelativisticMomentum(1, v), 1.25*1*v, 1e-12)
	assertApprox(t, "RelEnergy", RelativisticEnergy(1, v), 1.25*physicsC2, 1e-12)
	assertApprox(t, "RelKE", RelativisticKineticEnergy(1, v), 0.25*physicsC2, 1e-12)

	// Rest energy: total energy at v=0 equals m·c².
	assertApprox(t, "RelEnergy at rest", RelativisticEnergy(2, 0), MassEnergy(2), 1e-12)
}

func TestRelativisticNewtonianLimit(t *testing.T) {
	// For v ≪ c the relativistic kinetic energy reduces to ½·m·v².
	// The relativistic (γ−1) term loses precision by cancellation at low speed,
	// so this agreement is limited by floating-point rounding, not physics.
	const m, v = 3.0, 1000.0
	assertApprox(t, "KE limit", RelativisticKineticEnergy(m, v), KineticEnergy(m, v), 1e-4)
}

func TestEnergyMomentumRelation(t *testing.T) {
	// At rest, E = m·c².
	assertApprox(t, "EMR rest", EnergyMomentumRelation(1, 0), physicsC2, 1e-12)

	// Consistency: E-p relation reproduces γ·m·c² for a moving mass.
	const m, v = 1.0, 0.8 * SpeedOfLight // γ = 5/3
	p := RelativisticMomentum(m, v)
	assertApprox(t, "EMR moving", EnergyMomentumRelation(m, p), RelativisticEnergy(m, v), 1e-9)
}

func TestVelocityAddition(t *testing.T) {
	// 0.5c ⊕ 0.5c = 0.8c.
	assertApprox(t, "0.5c+0.5c", RelativisticVelocityAddition(0.5*SpeedOfLight, 0.5*SpeedOfLight),
		0.8*SpeedOfLight, 1e-9)
	// Adding c to anything yields c.
	assertApprox(t, "c+0.5c", RelativisticVelocityAddition(SpeedOfLight, 0.5*SpeedOfLight),
		SpeedOfLight, 1e-9)
}

func TestTimeDilationLengthContraction(t *testing.T) {
	const v = 0.6 * SpeedOfLight // γ = 1.25
	assertApprox(t, "TimeDilation", TimeDilation(1, v), 1.25, 1e-12)
	assertApprox(t, "LengthContraction", LengthContraction(1, v), 0.8, 1e-12)

	// Doppler factor is 1 at rest and < 1 when separating.
	assertApprox(t, "Doppler rest", RelativisticDopplerFactor(0), 1, 1e-12)
	if RelativisticDopplerFactor(0.5*SpeedOfLight) >= 1 {
		t.Error("separating source should be redshifted (factor < 1)")
	}
	if RelativisticDopplerFactor(-0.5*SpeedOfLight) <= 1 {
		t.Error("approaching source should be blueshifted (factor > 1)")
	}
}

func BenchmarkRelativisticEnergy(b *testing.B) {
	var acc float64
	for i := 0; i < b.N; i++ {
		acc += RelativisticEnergy(1, 0.6*SpeedOfLight)
	}
	_ = acc
}
