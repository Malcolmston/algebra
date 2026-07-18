package diffgeo

import (
	"math"
	"testing"
)

func TestSphereCurvature(t *testing.T) {
	R := 2.0
	s := Sphere(R)
	pts := [][2]float64{{0, 0}, {0.7, 0.3}, {2.1, -0.6}, {-1.0, 0.9}}
	for _, p := range pts {
		u, v := p[0], p[1]
		if got := GaussianCurvature(s, u, v); math.Abs(got-1/(R*R)) > 1e-3 {
			t.Errorf("(%v,%v): K = %v want %v", u, v, got, 1/(R*R))
		}
		if got := MeanCurvature(s, u, v); math.Abs(math.Abs(got)-1/R) > 1e-3 {
			t.Errorf("(%v,%v): |H| = %v want %v", u, v, math.Abs(got), 1/R)
		}
		k1, k2 := PrincipalCurvatures(s, u, v)
		if math.Abs(math.Abs(k1)-1/R) > 2e-3 || math.Abs(math.Abs(k2)-1/R) > 2e-3 {
			t.Errorf("(%v,%v): principal = %v,%v want ±%v", u, v, k1, k2, 1/R)
		}
		// Normal curvature is 1/R in every tangent direction (umbilic point).
		if got := math.Abs(NormalCurvature(s, u, v, 1, 0.5)); math.Abs(got-1/R) > 1e-3 {
			t.Errorf("(%v,%v): |normal curvature| = %v want %v", u, v, got, 1/R)
		}
	}
	// Outward unit normal at (0,0) is +x.
	if got := SurfaceNormal(s, 0, 0); !got.Equal(Vec3{1, 0, 0}, 1e-4) {
		t.Errorf("normal(0,0) = %v want +x", got)
	}
}

func TestPlaneCurvature(t *testing.T) {
	s := SurfacePlane(Vec3{1, 2, 3}, Vec3{2, 0, 0}, Vec3{0, 3, 1})
	for _, p := range [][2]float64{{0, 0}, {1.5, -2}, {3, 4}} {
		if got := GaussianCurvature(s, p[0], p[1]); math.Abs(got) > 1e-4 {
			t.Errorf("plane K = %v want 0", got)
		}
		if got := MeanCurvature(s, p[0], p[1]); math.Abs(got) > 1e-4 {
			t.Errorf("plane H = %v want 0", got)
		}
	}
}

func TestCylinderCurvature(t *testing.T) {
	R := 1.5
	s := Cylinder(R)
	for _, p := range [][2]float64{{0, 0}, {1.0, 2.0}, {-0.5, 1.0}} {
		if got := GaussianCurvature(s, p[0], p[1]); math.Abs(got) > 1e-3 {
			t.Errorf("cylinder K = %v want 0", got)
		}
		if got := math.Abs(MeanCurvature(s, p[0], p[1])); math.Abs(got-1/(2*R)) > 1e-3 {
			t.Errorf("cylinder |H| = %v want %v", got, 1/(2*R))
		}
	}
}

func TestTorusGaussianCurvature(t *testing.T) {
	Rr, rr := 3.0, 1.0
	s := Torus(Rr, rr)
	// K = cos v / (r (R + r cos v)).
	for _, v := range []float64{0, math.Pi, math.Pi / 2, 0.6, -1.2} {
		want := math.Cos(v) / (rr * (Rr + rr*math.Cos(v)))
		if got := GaussianCurvature(s, 1.1, v); math.Abs(got-want) > 1e-3 {
			t.Errorf("torus K(v=%v) = %v want %v", v, got, want)
		}
	}
}

func TestFundamentalFormsSphere(t *testing.T) {
	R := 2.0
	s := Sphere(R)
	u, v := 0.5, 0.3
	I := FirstFundamental(s, u, v)
	// For this parametrization: E = R^2 cos^2 v, F = 0, G = R^2.
	if math.Abs(I.E-R*R*math.Cos(v)*math.Cos(v)) > 1e-4 {
		t.Errorf("E = %v", I.E)
	}
	if math.Abs(I.F) > 1e-4 {
		t.Errorf("F = %v want 0", I.F)
	}
	if math.Abs(I.G-R*R) > 1e-4 {
		t.Errorf("G = %v want %v", I.G, R*R)
	}
	if math.Abs(I.Determinant()-(R*R*R*R*math.Cos(v)*math.Cos(v))) > 1e-3 {
		t.Errorf("det = %v", I.Determinant())
	}
}

func TestSurfaceAreaSphere(t *testing.T) {
	R := 2.0
	s := Sphere(R)
	// Full sphere area = 4*pi*R^2 (parametrized by u in [0,2pi], v in [-pi/2,pi/2]).
	got := SurfaceArea(s, 0, 2*math.Pi, -math.Pi/2, math.Pi/2, 80)
	want := 4 * math.Pi * R * R
	if math.Abs(got-want) > 1e-2 {
		t.Errorf("sphere area = %v want %v", got, want)
	}
}

func TestAreaElementPlane(t *testing.T) {
	// Unit orthogonal spanning vectors give area element 1 everywhere.
	s := SurfacePlane(Vec3{}, Vec3{1, 0, 0}, Vec3{0, 1, 0})
	if got := AreaElement(s, 2, 3); math.Abs(got-1) > 1e-6 {
		t.Errorf("area element = %v want 1", got)
	}
}

func TestChristoffelPlaneZero(t *testing.T) {
	s := SurfacePlane(Vec3{}, Vec3{2, 0, 0}, Vec3{0, 3, 0})
	c := ChristoffelSymbols(s, 1.0, 2.0)
	for k := 0; k < 2; k++ {
		for i := 0; i < 2; i++ {
			for j := 0; j < 2; j++ {
				if math.Abs(c.At(k, i, j)) > 1e-4 {
					t.Errorf("Gamma[%d][%d][%d] = %v want 0", k, i, j, c.At(k, i, j))
				}
			}
		}
	}
}

func TestChristoffelSymmetry(t *testing.T) {
	// Symmetric in the lower indices for any surface.
	s := Torus(3, 1)
	c := ChristoffelSymbols(s, 0.8, 1.4)
	for k := 0; k < 2; k++ {
		if math.Abs(c.At(k, 0, 1)-c.At(k, 1, 0)) > 1e-6 {
			t.Errorf("Gamma[%d] not symmetric: %v vs %v", k, c.At(k, 0, 1), c.At(k, 1, 0))
		}
	}
}
