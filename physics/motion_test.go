package physics

import (
	"errors"
	"math"
	"testing"
)

// motionApproxEqual reports whether got and want agree to within a combined
// absolute/relative tolerance, so table tests can use tidy decimal expectations
// without being defeated by floating-point rounding.
func motionApproxEqual(got, want, tol float64) bool {
	diff := math.Abs(got - want)
	if diff <= tol {
		return true
	}
	return diff <= tol*math.Max(math.Abs(got), math.Abs(want))
}

func TestVelocity(t *testing.T) {
	cases := []struct {
		name            string
		v0, a, tm, want float64
	}{
		{"accelerating", 2, 3, 4, 14},
		{"constant speed", 5, 0, 10, 5},
		{"decelerating to rest", 10, -2, 5, 0},
		{"from rest", 0, 9.8, 2, 19.6},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			if got := Velocity(c.v0, c.a, c.tm); !motionApproxEqual(got, c.want, 1e-12) {
				t.Fatalf("Velocity(%v,%v,%v) = %v, want %v", c.v0, c.a, c.tm, got, c.want)
			}
		})
	}
}

func TestVelocityFromDistance(t *testing.T) {
	cases := []struct {
		name           string
		v0, a, d, want float64
	}{
		{"3-4-5 relation", 3, 2, 4, 5},
		{"from rest", 0, 2, 8, 4 * math.Sqrt2}, // v = √(v0²+2ad) = √32 = 4√2
		{"no motion", 7, 0, 100, 7},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			if got := VelocityFromDistance(c.v0, c.a, c.d); !motionApproxEqual(got, c.want, 1e-12) {
				t.Fatalf("VelocityFromDistance(%v,%v,%v) = %v, want %v", c.v0, c.a, c.d, got, c.want)
			}
		})
	}

	// A negative radicand is unreachable and must surface as NaN.
	if got := VelocityFromDistance(1, -1, 10); !math.IsNaN(got) {
		t.Fatalf("VelocityFromDistance(1,-1,10) = %v, want NaN", got)
	}
}

func TestTimeToStop(t *testing.T) {
	cases := []struct {
		name        string
		v0, a, want float64
	}{
		{"deceleration", 10, -2, 5},
		{"reverse braking", -10, 2, 5},
		{"already at rest", 0, -3, 0},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			if got := TimeToStop(c.v0, c.a); !motionApproxEqual(got, c.want, 1e-12) {
				t.Fatalf("TimeToStop(%v,%v) = %v, want %v", c.v0, c.a, got, c.want)
			}
		})
	}
}

func TestImpulseMomentum(t *testing.T) {
	cases := []struct {
		name            string
		force, dt, want float64
	}{
		{"positive impulse", 5, 3, 15},
		{"zero interval", 100, 0, 0},
		{"retarding force", -4, 2.5, -10},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			if got := ImpulseMomentum(c.force, c.dt); !motionApproxEqual(got, c.want, 1e-12) {
				t.Fatalf("ImpulseMomentum(%v,%v) = %v, want %v", c.force, c.dt, got, c.want)
			}
		})
	}
}

func TestCenterOfMass1D(t *testing.T) {
	cases := []struct {
		name      string
		masses    []float64
		positions []float64
		want      float64
	}{
		{"two unequal masses", []float64{1, 3}, []float64{0, 4}, 3},
		{"symmetric equal masses", []float64{2, 2}, []float64{-1, 1}, 0},
		{"single mass", []float64{5}, []float64{7}, 7},
		{"three masses", []float64{1, 1, 2}, []float64{0, 3, 3}, 2.25},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			got, err := CenterOfMass1D(c.masses, c.positions)
			if err != nil {
				t.Fatalf("CenterOfMass1D returned unexpected error: %v", err)
			}
			if !motionApproxEqual(got, c.want, 1e-12) {
				t.Fatalf("CenterOfMass1D(%v,%v) = %v, want %v", c.masses, c.positions, got, c.want)
			}
		})
	}

	errCases := []struct {
		name      string
		masses    []float64
		positions []float64
	}{
		{"length mismatch", []float64{1, 2}, []float64{0}},
		{"empty", []float64{}, []float64{}},
		{"zero total mass", []float64{1, -1}, []float64{0, 5}},
	}
	for _, c := range errCases {
		t.Run("err/"+c.name, func(t *testing.T) {
			if _, err := CenterOfMass1D(c.masses, c.positions); !errors.Is(err, ErrCenterOfMass) {
				t.Fatalf("CenterOfMass1D(%v,%v) error = %v, want ErrCenterOfMass", c.masses, c.positions, err)
			}
		})
	}
}

func TestGravitationalPotentialEnergy(t *testing.T) {
	cases := []struct {
		name            string
		m1, m2, r, want float64
	}{
		{"attractive pair", 1000, 1000, 2, -3.33715e-5},
		{"unit separation", 1000, 1000, 1, -6.67430e-5},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			if got := GravitationalPotentialEnergy(c.m1, c.m2, c.r); !motionApproxEqual(got, c.want, 1e-9) {
				t.Fatalf("GravitationalPotentialEnergy(%v,%v,%v) = %v, want %v", c.m1, c.m2, c.r, got, c.want)
			}
		})
	}
}

func TestVisViva(t *testing.T) {
	// Circular orbit (r == a): vis-viva must reduce to the circular orbital
	// speed √(G·M/r).
	mCentral, r, a := 1e10, 1e5, 1e5
	want := math.Sqrt(GravitationalConstant * mCentral / r)
	if got := VisViva(mCentral, r, a); !motionApproxEqual(got, want, 1e-12) {
		t.Fatalf("VisViva circular = %v, want %v", got, want)
	}

	// Elliptical orbit at periapsis: independent closed-form expectation.
	mCentral, r, a = 1e10, 5e4, 1e5
	want = math.Sqrt(GravitationalConstant * mCentral * (2/r - 1/a))
	if got := VisViva(mCentral, r, a); !motionApproxEqual(got, want, 1e-12) {
		t.Fatalf("VisViva elliptical = %v, want %v", got, want)
	}
}

// --- Benchmarks for the performance-sensitive routines ---

func BenchmarkCenterOfMass1D(b *testing.B) {
	masses := make([]float64, 1024)
	positions := make([]float64, 1024)
	for i := range masses {
		masses[i] = float64(i%7 + 1)
		positions[i] = float64(i) * 0.5
	}
	b.ReportAllocs()
	b.ResetTimer()
	var sink float64
	for i := 0; i < b.N; i++ {
		v, err := CenterOfMass1D(masses, positions)
		if err != nil {
			b.Fatal(err)
		}
		sink += v
	}
	_ = sink
}

func BenchmarkVisViva(b *testing.B) {
	var sink float64
	for i := 0; i < b.N; i++ {
		sink += VisViva(1e24, 7.0e6, 8.0e6)
	}
	_ = sink
}

func BenchmarkVelocityFromDistance(b *testing.B) {
	var sink float64
	for i := 0; i < b.N; i++ {
		sink += VelocityFromDistance(3, 2, 4)
	}
	_ = sink
}
