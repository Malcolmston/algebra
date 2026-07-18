package complexanalysis

import (
	"math"
	"testing"
)

func TestIntegrateCircle(t *testing.T) {
	// integral of 1/z around unit circle = 2*pi*i.
	got := IntegrateCircle(func(z complex128) complex128 { return 1 / z }, 0, 1, 512)
	if !closeC(got, complex(0, 2*math.Pi), 1e-9) {
		t.Errorf("closed 1/z = %v, want 2*pi*i", got)
	}
	// integral of z (analytic) around any closed contour = 0.
	got = IntegrateCircle(func(z complex128) complex128 { return z }, complex(1, 1), 2, 256)
	if !closeC(got, 0, 1e-9) {
		t.Errorf("closed z = %v, want 0", got)
	}
}

func TestIntegrateSegment(t *testing.T) {
	// integral of z from 0 to 1+i equals (1+i)^2/2 = i.
	got := IntegrateSegment(func(z complex128) complex128 { return z }, 0, complex(1, 1), 100)
	if !closeC(got, complex(0, 1), 1e-9) {
		t.Errorf("segment = %v, want i", got)
	}
}

func TestIntegratePolygonAndContour(t *testing.T) {
	// A closed square around the origin: integral of 1/z = 2*pi*i.
	square := []complex128{complex(1, 1), complex(-1, 1), complex(-1, -1), complex(1, -1)}
	got := IntegratePolygon(func(z complex128) complex128 { return 1 / z }, square, 400)
	if !closeC(got, complex(0, 2*math.Pi), 1e-6) {
		t.Errorf("polygon 1/z = %v, want 2*pi*i", got)
	}
	// ContourIntegral over unit circle parametrization of 1/z.
	gamma := func(tt float64) complex128 { return complex(math.Cos(2*math.Pi*tt), math.Sin(2*math.Pi*tt)) }
	gc := ContourIntegral(func(z complex128) complex128 { return 1 / z }, gamma, 2000)
	if !closeC(gc, complex(0, 2*math.Pi), 1e-4) {
		t.Errorf("ContourIntegral 1/z = %v, want 2*pi*i", gc)
	}
}

func TestResidues(t *testing.T) {
	// Res 1/(z^2+1) at i = -i/2.
	f := func(z complex128) complex128 { return 1 / (z*z + 1) }
	if !closeC(Residue(f, complex(0, 1), 0.3, 400), complex(0, -0.5), 1e-9) {
		t.Error("Residue 1/(z^2+1) at i wrong")
	}
	if !closeC(ResidueSimplePole(f, complex(0, 1), 1e-4), complex(0, -0.5), 1e-6) {
		t.Error("ResidueSimplePole wrong")
	}
	// Res exp(z)/z^2 at 0 (order 2) = 1.
	g := func(z complex128) complex128 { return Exp(z) / (z * z) }
	if !closeC(ResidueOrderM(g, 0, 2, 1, 600), 1, 1e-8) {
		t.Error("ResidueOrderM exp/z^2 wrong")
	}
}

func TestCauchyFormulas(t *testing.T) {
	// f(0) = 1 for exp via Cauchy integral.
	if !closeC(CauchyIntegralValue(Exp, 0, 0, 1, 400), 1, 1e-9) {
		t.Error("CauchyIntegralValue exp(0) wrong")
	}
	// f''(0) = 1 for exp.
	if !closeC(CauchyDerivative(Exp, 0, 0, 2, 1, 400), 1, 1e-8) {
		t.Error("CauchyDerivative exp''(0) wrong")
	}
	// derivative of z^3 at z=1 of order 1 is 3.
	cube := func(z complex128) complex128 { return z * z * z }
	if !closeC(CauchyDerivative(cube, 1, 1, 1, 0.5, 400), 3, 1e-8) {
		t.Error("CauchyDerivative (z^3)'(1) wrong")
	}
}

func TestWindingNumber(t *testing.T) {
	unit := func(tt float64) complex128 { return complex(math.Cos(2*math.Pi*tt), math.Sin(2*math.Pi*tt)) }
	if WindingNumber(unit, 0, 200) != 1 {
		t.Error("winding about 0 should be 1")
	}
	if WindingNumber(unit, 2, 200) != 0 {
		t.Error("winding about exterior point should be 0")
	}
	twice := func(tt float64) complex128 { return complex(math.Cos(4*math.Pi*tt), math.Sin(4*math.Pi*tt)) }
	if WindingNumber(twice, 0, 400) != 2 {
		t.Error("double loop winding should be 2")
	}
	// negative orientation
	rev := func(tt float64) complex128 { return complex(math.Cos(-2*math.Pi*tt), math.Sin(-2*math.Pi*tt)) }
	if WindingNumber(rev, 0, 200) != -1 {
		t.Error("reversed loop winding should be -1")
	}
}

func TestArgumentPrinciple(t *testing.T) {
	// z^2 has a double zero inside the unit disk.
	if CountZeros(func(z complex128) complex128 { return z * z }, 0, 1, 400) != 2 {
		t.Error("CountZeros z^2 should be 2")
	}
	// (z-0.5)(z+0.5) has two simple zeros inside the unit disk.
	f := func(z complex128) complex128 { return (z - 0.5) * (z + 0.5) }
	if CountZeros(f, 0, 1, 400) != 2 {
		t.Error("CountZeros should be 2")
	}
	// A function analytic and non-zero inside: exp has no zeros.
	if CountZeros(Exp, 0, 1, 400) != 0 {
		t.Error("CountZeros exp should be 0")
	}
}
