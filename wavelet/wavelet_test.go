package wavelet

import (
	"math"
	"testing"
)

// approxEqual reports whether a and b agree within tol.
func approxEqual(a, b, tol float64) bool { return math.Abs(a-b) <= tol }

// sliceApproxEqual reports whether two slices agree elementwise within tol.
func sliceApproxEqual(a, b []float64, tol float64) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if math.Abs(a[i]-b[i]) > tol {
			return false
		}
	}
	return true
}

func TestFilterSumsAndOrthogonality(t *testing.T) {
	for _, w := range []Wavelet{Haar(), DB2(), DB4()} {
		// Low-pass coefficients sum to sqrt(2).
		var sumLo float64
		for _, c := range w.DecLo() {
			sumLo += c
		}
		if !approxEqual(sumLo, math.Sqrt2, 1e-12) {
			t.Errorf("%s: sum of low-pass = %v want sqrt(2)", w.Name(), sumLo)
		}
		// High-pass coefficients sum to zero (>= 1 vanishing moment).
		var sumHi float64
		for _, c := range w.DecHi() {
			sumHi += c
		}
		if !approxEqual(sumHi, 0, 1e-12) {
			t.Errorf("%s: sum of high-pass = %v want 0", w.Name(), sumHi)
		}
		// Unit energy of the scaling filter.
		if !approxEqual(Energy(w.DecLo()), 1, 1e-12) {
			t.Errorf("%s: scaling filter energy = %v want 1", w.Name(), Energy(w.DecLo()))
		}
		if !w.IsOrthogonal(1e-10) {
			t.Errorf("%s: reported non-orthogonal", w.Name())
		}
	}
}

func TestDaubechiesReferenceCoefficients(t *testing.T) {
	// db2 closed-form scaling coefficients (1 +/- sqrt3, 3 +/- sqrt3)/(4 sqrt2).
	s3 := math.Sqrt(3)
	s2 := math.Sqrt2
	wantDB2 := []float64{
		(1 + s3) / (4 * s2),
		(3 + s3) / (4 * s2),
		(3 - s3) / (4 * s2),
		(1 - s3) / (4 * s2),
	}
	if !sliceApproxEqual(DB2().DecLo(), wantDB2, 1e-12) {
		t.Errorf("db2 DecLo = %v want %v", DB2().DecLo(), wantDB2)
	}
	// db4 tabulated reference values.
	wantDB4 := []float64{
		0.23037781330885523, 0.7148465705525415, 0.6308807679295904, -0.02798376941698385,
		-0.18703481171888114, 0.030841381835986965, 0.032883011666982945, -0.010597401784997278,
	}
	if !sliceApproxEqual(DB4().DecLo(), wantDB4, 1e-15) {
		t.Errorf("db4 DecLo mismatch")
	}
}

func TestDaubechiesConstructor(t *testing.T) {
	if _, err := Daubechies(5); err == nil {
		t.Error("Daubechies(5) should error")
	}
	if w, err := Daubechies(1); err != nil || w.Name() != "haar" {
		t.Errorf("Daubechies(1) = %v, %v want haar", w.Name(), err)
	}
	for _, p := range []int{2, 3, 4} {
		w, err := Daubechies(p)
		if err != nil {
			t.Fatalf("Daubechies(%d) error: %v", p, err)
		}
		if w.Taps() != 2*p {
			t.Errorf("Daubechies(%d) taps = %d want %d", p, w.Taps(), 2*p)
		}
		if w.VanishingMoments() != p {
			t.Errorf("Daubechies(%d) vanishing moments = %d want %d", p, w.VanishingMoments(), p)
		}
	}
}

func TestWaveletAccessorsAreCopies(t *testing.T) {
	w := Haar()
	lo := w.DecLo()
	lo[0] = 999
	if approxEqual(w.DecLo()[0], 999, 1e-9) {
		t.Error("DecLo returned a mutable reference to internal state")
	}
}
