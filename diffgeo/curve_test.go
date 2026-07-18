package diffgeo

import (
	"math"
	"testing"
)

func TestCurvatureKnown(t *testing.T) {
	tests := []struct {
		name   string
		c      Curve
		t      float64
		kappa  float64
		tau    float64
		tolK   float64
		tolTau float64
	}{
		{"circle R=3", Circle(3), 0.7, 1.0 / 3.0, 0, 1e-5, 1e-5},
		{"circle R=0.5", Circle(0.5), 2.1, 2.0, 0, 1e-4, 1e-5},
		// Helix a=2,b=1: kappa=a/(a^2+b^2)=2/5, tau=b/(a^2+b^2)=1/5.
		{"helix a=2 b=1", Helix(2, 1), 1.3, 2.0 / 5.0, 1.0 / 5.0, 1e-4, 1e-3},
		{"helix a=1 b=2", Helix(1, 2), 0.4, 1.0 / 5.0, 2.0 / 5.0, 1e-4, 1e-3},
		// Straight line: zero curvature and torsion.
		{"line", Line(Vec3{1, 2, 3}, Vec3{2, -1, 4}), 5.0, 0, 0, 1e-6, 1e-6},
	}
	for _, tc := range tests {
		if got := Curvature(tc.c, tc.t); math.Abs(got-tc.kappa) > tc.tolK {
			t.Errorf("%s: Curvature = %v want %v", tc.name, got, tc.kappa)
		}
		if got := Torsion(tc.c, tc.t); math.Abs(got-tc.tau) > tc.tolTau {
			t.Errorf("%s: Torsion = %v want %v", tc.name, got, tc.tau)
		}
	}
}

func TestRadiusOfCurvature(t *testing.T) {
	if got := RadiusOfCurvature(Circle(4), 1.1); math.Abs(got-4) > 1e-4 {
		t.Errorf("radius = %v want 4", got)
	}
	if got := RadiusOfCurvature(Line(Vec3{}, Vec3{1, 1, 1}), 0.5); !math.IsInf(got, 1) {
		t.Errorf("line radius = %v want +Inf", got)
	}
}

func TestFrenetFrameOrthonormal(t *testing.T) {
	c := Helix(2, 1)
	for _, tt := range []float64{0.0, 0.9, 2.7, 4.4} {
		f := FrenetFrame(c, tt)
		for _, u := range []Vec3{f.T, f.N, f.B} {
			if math.Abs(u.Norm()-1) > 1e-6 {
				t.Errorf("t=%v: non-unit basis vector %v", tt, u)
			}
		}
		if math.Abs(f.T.Dot(f.N)) > 1e-6 || math.Abs(f.T.Dot(f.B)) > 1e-6 || math.Abs(f.N.Dot(f.B)) > 1e-6 {
			t.Errorf("t=%v: frame not orthogonal", tt)
		}
		// Right-handed: T x N = B.
		if !f.T.Cross(f.N).Equal(f.B, 1e-6) {
			t.Errorf("t=%v: T x N != B", tt)
		}
	}
}

func TestUnitTangentAndVelocity(t *testing.T) {
	// Circle tangent at t=0 points in +y.
	if got := UnitTangent(Circle(3), 0); !got.Equal(Vec3{0, 1, 0}, 1e-6) {
		t.Errorf("tangent = %v want +y", got)
	}
	// Speed of circle radius R is R (constant).
	if got := Speed(Circle(3), 1.7); math.Abs(got-3) > 1e-5 {
		t.Errorf("speed = %v want 3", got)
	}
	// Principal normal of a circle points toward the center: at t=0 that is -x.
	if got := PrincipalNormal(Circle(3), 0); !got.Equal(Vec3{-1, 0, 0}, 1e-5) {
		t.Errorf("normal = %v want -x", got)
	}
}

func TestArcLengthKnown(t *testing.T) {
	// Circle radius 3 over one full turn has length 2*pi*3.
	if got := ArcLength(Circle(3), 0, 2*math.Pi, 200); math.Abs(got-6*math.Pi) > 1e-6 {
		t.Errorf("circle arc = %v want %v", got, 6*math.Pi)
	}
	// Helix a=2,b=1 over one turn: 2*pi*sqrt(a^2+b^2).
	want := 2 * math.Pi * math.Sqrt(5)
	if got := ArcLength(Helix(2, 1), 0, 2*math.Pi, 200); math.Abs(got-want) > 1e-6 {
		t.Errorf("helix arc = %v want %v", got, want)
	}
	// Straight line length equals |d|*(b-a).
	if got := ArcLength(Line(Vec3{}, Vec3{3, 4, 0}), 0, 2, 10); math.Abs(got-10) > 1e-6 {
		t.Errorf("line arc = %v want 10", got)
	}
}

func TestCurvatureVectorMagnitude(t *testing.T) {
	// |dT/ds| equals the curvature.
	c := Ellipse(3, 2)
	for _, tt := range []float64{0.3, 1.2, 2.5} {
		if math.Abs(CurvatureVector(c, tt).Norm()-Curvature(c, tt)) > 1e-4 {
			t.Errorf("t=%v: |curvature vector| != curvature", tt)
		}
	}
}
