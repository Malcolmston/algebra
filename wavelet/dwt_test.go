package wavelet

import (
	"math"
	"testing"
)

func TestHaarKnownAnswer(t *testing.T) {
	// Haar of [1,2,3,4]: approx = (x0+x1)/sqrt2, (x2+x3)/sqrt2;
	// detail = (x0-x1)/sqrt2, (x2-x3)/sqrt2.
	x := []float64{1, 2, 3, 4}
	a, d := DWT(x, Haar())
	wantA := []float64{3 / math.Sqrt2, 7 / math.Sqrt2}
	wantD := []float64{-1 / math.Sqrt2, -1 / math.Sqrt2}
	if !sliceApproxEqual(a, wantA, 1e-12) {
		t.Errorf("approx = %v want %v", a, wantA)
	}
	if !sliceApproxEqual(d, wantD, 1e-12) {
		t.Errorf("detail = %v want %v", d, wantD)
	}
}

func TestPerfectReconstruction(t *testing.T) {
	signals := [][]float64{
		{1, 2, 3, 4, 5, 6, 7, 8},
		{-3, 7, 0.5, -2.25, 11, -8, 4, 4, 1, 0, -6, 9, 2, 2, 2, -1},
	}
	for _, w := range []Wavelet{Haar(), DB2(), DB4()} {
		for _, x := range signals {
			a, d := DWT(x, w)
			r := IDWT(a, d, w)
			if !sliceApproxEqual(r, x, 1e-10) {
				t.Errorf("%s: reconstruction mismatch\n got %v\nwant %v", w.Name(), r, x)
			}
		}
	}
}

func TestParseval(t *testing.T) {
	x := []float64{1, -2, 3, -4, 5, -6, 7, -8}
	for _, w := range []Wavelet{Haar(), DB2(), DB4()} {
		c := Transform(x, w)
		if !approxEqual(c.Energy(), Energy(x), 1e-10) {
			t.Errorf("%s: coefficient energy = %v want %v", w.Name(), c.Energy(), Energy(x))
		}
	}
}

func TestConstantSignalZeroDetail(t *testing.T) {
	// A constant signal is annihilated by every high-pass filter with >= 1
	// vanishing moment, so all detail coefficients must vanish.
	c := []float64{2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5}
	for _, w := range []Wavelet{Haar(), DB2(), DB4()} {
		_, d := DWT(c, w)
		for i, v := range d {
			if math.Abs(v) > 1e-12 {
				t.Errorf("%s: detail[%d] = %v want 0", w.Name(), i, v)
			}
		}
	}
}

func TestLinearSignalVanishingMoments(t *testing.T) {
	// db2 and db4 have >= 2 vanishing moments, so a linear ramp produces zero
	// detail away from the periodic boundary. Check an interior coefficient.
	n := 32
	x := make([]float64, n)
	for i := range x {
		x[i] = 0.5*float64(i) + 3
	}
	for _, w := range []Wavelet{DB2(), DB4()} {
		_, d := DWT(x, w)
		// Interior index safely away from the wrap-around boundary.
		mid := len(d) / 2
		if math.Abs(d[mid]) > 1e-9 {
			t.Errorf("%s: interior detail = %v want ~0 (linear ramp)", w.Name(), d[mid])
		}
	}
}

func TestDWTPanicsOddLength(t *testing.T) {
	defer func() {
		if recover() == nil {
			t.Error("DWT of odd-length signal should panic")
		}
	}()
	DWT([]float64{1, 2, 3}, Haar())
}

func TestCoefficientsRoundTrip(t *testing.T) {
	x := []float64{9, -1, 4, 4, 2, 0, -7, 3}
	c := Transform(x, DB2())
	if c.Length() != len(x) {
		t.Errorf("Length = %d want %d", c.Length(), len(x))
	}
	if !sliceApproxEqual(c.Inverse(DB2()), x, 1e-10) {
		t.Error("Coefficients.Inverse round trip failed")
	}
}
