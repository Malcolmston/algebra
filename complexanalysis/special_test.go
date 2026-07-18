package complexanalysis

import (
	"math"
	"testing"
)

func TestSpecialFunctions(t *testing.T) {
	tests := []struct {
		name string
		got  complex128
		want complex128
		tol  float64
	}{
		{"Gamma(5)", Gamma(5), 24, 1e-9},
		{"Gamma(0.5)", Gamma(0.5), complex(math.Sqrt(math.Pi), 0), 1e-10},
		{"Gamma(1+i)", Gamma(complex(1, 1)), complex(0.4980156681, -0.1549498283), 1e-9},
		{"LogGamma(10)", LogGamma(10), complex(math.Log(362880), 0), 1e-9},
		{"Digamma(1)", Digamma(1), complex(-0.5772156649, 0), 1e-9},
		{"Digamma(2)", Digamma(2), complex(0.4227843351, 0), 1e-9},
		{"Beta(2,3)", Beta(2, 3), complex(1.0/12, 0), 1e-10},
		{"Erf(1)", Erf(1), complex(0.8427007929, 0), 1e-10},
		{"Erfc(1)", Erfc(1), complex(1-0.8427007929, 0), 1e-10},
		{"Zeta(2)", Zeta(2), complex(math.Pi*math.Pi/6, 0), 1e-9},
		{"Zeta(4)", Zeta(4), complex(math.Pow(math.Pi, 4)/90, 0), 1e-9},
		{"Zeta(0)", Zeta(0), -0.5, 1e-9},
		{"Zeta(-1)", Zeta(-1), complex(-1.0/12, 0), 1e-8},
		{"Factorial(5)", Factorial(5), 120, 1e-8},
		{"RisingFactorial(1,5)", RisingFactorial(1, 5), 120, 1e-9},
		{"RisingFactorial(2,3)", RisingFactorial(2, 3), 24, 1e-9},
		{"Binomial(5,2)", Binomial(5, 2), 10, 1e-9},
		{"Binomial(6,3)", Binomial(6, 3), 20, 1e-9},
	}
	for _, tc := range tests {
		if !closeC(tc.got, tc.want, tc.tol) {
			t.Errorf("%s = %v, want %v", tc.name, tc.got, tc.want)
		}
	}
}

func TestGammaReflection(t *testing.T) {
	// Gamma(z)Gamma(1-z) = pi/sin(pi z).
	z := complex(0.3, 0.2)
	lhs := Gamma(z) * Gamma(1-z)
	rhs := complex(math.Pi, 0) / Sin(math.Pi*z)
	if !closeC(lhs, rhs, 1e-9) {
		t.Errorf("reflection: %v vs %v", lhs, rhs)
	}
}

func TestErfImaginary(t *testing.T) {
	// erf(i) = i * 1.6504257587975428... (a real multiple of i).
	got := Erf(complex(0, 1))
	want := complex(0, 1.6504257587975428)
	if !closeC(got, want, 1e-10) {
		t.Errorf("Erf(i) = %v, want %v", got, want)
	}
}
