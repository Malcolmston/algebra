package complexanalysis

import (
	"math"
	"math/cmplx"
	"testing"
)

const testTol = 1e-9

func closeC(a, b complex128, tol float64) bool { return cmplx.Abs(a-b) <= tol }

func TestElementary(t *testing.T) {
	tests := []struct {
		name string
		got  complex128
		want complex128
	}{
		{"Sqrt(-1)", Sqrt(-1), complex(0, 1)},
		{"Sqrt(4)", Sqrt(4), 2},
		{"Cbrt(8)", Cbrt(8), 2},
		{"Cbrt(0)", Cbrt(0), 0},
		{"Exp(i*pi)", Exp(complex(0, math.Pi)), -1},
		{"Log(-1)", Log(-1), complex(0, math.Pi)},
		{"LogBranch(-1,1)", LogBranch(-1, 1), complex(0, 3*math.Pi)},
		{"Pow(i,2)", Pow(complex(0, 1), 2), -1},
		{"Pow(0,0)", Pow(0, 0), 1},
		{"Reciprocal(2i)", Reciprocal(complex(0, 2)), complex(0, -0.5)},
		{"Sign(3)", Sign(3), 1},
		{"Sign(0)", Sign(0), 0},
		{"Sin(pi/2)", Sin(math.Pi / 2), 1},
		{"Cos(0)", Cos(0), 1},
		{"Tan(0)", Tan(0), 0},
		{"Cot(pi/4)", Cot(math.Pi / 4), 1},
		{"Sec(0)", Sec(0), 1},
		{"Csc(pi/2)", Csc(math.Pi / 2), 1},
		{"Sinh(0)", Sinh(0), 0},
		{"Cosh(0)", Cosh(0), 1},
		{"Tanh(0)", Tanh(0), 0},
		{"Asin(1)", Asin(1), math.Pi / 2},
		{"Acos(1)", Acos(1), 0},
		{"Atan(0)", Atan(0), 0},
		{"Asinh(0)", Asinh(0), 0},
		{"Acosh(1)", Acosh(1), 0},
		{"Atanh(0)", Atanh(0), 0},
	}
	for _, tc := range tests {
		if !closeC(tc.got, tc.want, testTol) {
			t.Errorf("%s = %v, want %v", tc.name, tc.got, tc.want)
		}
	}
}

func TestNthRoots(t *testing.T) {
	roots := NthRoots(1, 3)
	if len(roots) != 3 {
		t.Fatalf("expected 3 roots, got %d", len(roots))
	}
	want := []complex128{1, cmplx.Rect(1, 2*math.Pi/3), cmplx.Rect(1, 4*math.Pi/3)}
	for i := range roots {
		if !closeC(roots[i], want[i], testTol) {
			t.Errorf("root %d = %v, want %v", i, roots[i], want[i])
		}
		// Each root cubed must return to 1.
		if !closeC(roots[i]*roots[i]*roots[i], 1, testTol) {
			t.Errorf("root %d cubed = %v, want 1", i, roots[i]*roots[i]*roots[i])
		}
	}
	if NthRoots(1, 0) != nil {
		t.Error("NthRoots with n=0 should be nil")
	}
}

func TestPythagoreanIdentities(t *testing.T) {
	z := complex(0.7, -0.4)
	if !closeC(Sin(z)*Sin(z)+Cos(z)*Cos(z), 1, testTol) {
		t.Error("sin^2+cos^2 != 1")
	}
	if !closeC(Cosh(z)*Cosh(z)-Sinh(z)*Sinh(z), 1, testTol) {
		t.Error("cosh^2-sinh^2 != 1")
	}
}
