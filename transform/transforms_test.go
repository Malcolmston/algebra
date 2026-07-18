package transform

import (
	"math"
	"math/cmplx"
	"testing"
)

func TestDCTKnown(t *testing.T) {
	// The DCT-II of a constant vector is a single DC coefficient.
	got := DCT([]float64{1, 1, 1, 1})
	if !floatClose(got, []float64{4, 0, 0, 0}, 1e-9) {
		t.Errorf("DCT const = %v", got)
	}
}

func TestDCTIDCTRoundTrip(t *testing.T) {
	x := []float64{3, 1, 4, 1, 5, 9}
	if !floatClose(IDCT(DCT(x)), x, 1e-9) {
		t.Errorf("IDCT(DCT(x)) != x")
	}
}

func TestDCT4Involution(t *testing.T) {
	// DCT-IV applied twice equals the identity scaled by N/2.
	x := []float64{2, -1, 0.5, 3}
	y := DCT4(DCT4(x))
	scale := 2.0 / float64(len(x))
	for i := range y {
		y[i] *= scale
	}
	if !floatClose(y, x, 1e-9) {
		t.Errorf("DCT4 involution failed: %v", y)
	}
}

func TestDCT1RoundTrip(t *testing.T) {
	x := []float64{1, 2, 3, 4, 5}
	y := DCT1(DCT1(x))
	scale := 2.0 / float64(len(x)-1)
	for i := range y {
		y[i] *= scale
	}
	if !floatClose(y, x, 1e-9) {
		t.Errorf("DCT1 round trip failed: %v", y)
	}
}

func TestDSTRoundTrip(t *testing.T) {
	x := []float64{3, 1, 4, 1, 5}
	if !floatClose(IDST(DST(x)), x, 1e-9) {
		t.Errorf("IDST(DST(x)) != x")
	}
}

func TestConvolveKnown(t *testing.T) {
	got := Convolve([]float64{1, 1, 1}, []float64{1, 1, 1})
	if !floatClose(got, []float64{1, 2, 3, 2, 1}, 1e-12) {
		t.Errorf("Convolve = %v", got)
	}
	// Multiplying the polynomials (1+2x)(1+3x) = 1+5x+6x^2.
	got = Convolve([]float64{1, 2}, []float64{1, 3})
	if !floatClose(got, []float64{1, 5, 6}, 1e-12) {
		t.Errorf("Convolve poly = %v", got)
	}
}

func TestConvolveFFTMatches(t *testing.T) {
	a := []float64{1, 2, 3, 4, 5}
	b := []float64{2, -1, 0.5}
	if !floatClose(ConvolveFFT(a, b), Convolve(a, b), 1e-9) {
		t.Errorf("ConvolveFFT != Convolve")
	}
}

func TestConvolveComplex(t *testing.T) {
	got := ConvolveComplex([]complex128{1, 1i}, []complex128{1, 1i})
	want := []complex128{1, 2i, -1}
	if !cmplxClose(got, want, 1e-12) {
		t.Errorf("ConvolveComplex = %v", got)
	}
}

func TestCircularConvolve(t *testing.T) {
	// Circular convolution with a shift filter [1,1,0,0].
	got := CircularConvolve([]float64{1, 2, 3, 4}, []float64{1, 1, 0, 0})
	want := []float64{5, 3, 5, 7}
	if !floatClose(got, want, 1e-12) {
		t.Errorf("CircularConvolve = %v", got)
	}
}

func TestCorrelate(t *testing.T) {
	got := Correlate([]float64{1, 2, 3}, []float64{0, 1, 0.5})
	if !floatClose(got, CrossCorrelateFFT([]float64{1, 2, 3}, []float64{0, 1, 0.5}), 1e-9) {
		t.Errorf("CrossCorrelateFFT != Correlate")
	}
	// Autocorrelation is symmetric and peaks at the center (signal energy).
	ac := AutoCorrelate([]float64{1, 2, 3})
	center := len(ac) / 2
	if math.Abs(ac[center]-14) > 1e-12 {
		t.Errorf("AutoCorrelate center = %v, want 14", ac[center])
	}
	if math.Abs(ac[0]-ac[len(ac)-1]) > 1e-12 {
		t.Errorf("AutoCorrelate not symmetric")
	}
	_ = got
}

func TestLaplaceForward(t *testing.T) {
	// L{1}(s) = 1/s.
	got := Laplace(func(t float64) float64 { return 1 }, complex(2, 0), 60, 4000)
	if cmplx.Abs(got-0.5) > 1e-5 {
		t.Errorf("Laplace{1}(2) = %v", got)
	}
	// L{e^{-t}}(1) = 1/(1+1) = 0.5.
	got = Laplace(func(t float64) float64 { return math.Exp(-t) }, complex(1, 0), 60, 4000)
	if cmplx.Abs(got-0.5) > 1e-5 {
		t.Errorf("Laplace{e^-t}(1) = %v", got)
	}
}

func TestInverseLaplaceTalbot(t *testing.T) {
	cases := []struct {
		name string
		F    func(complex128) complex128
		t    float64
		want float64
	}{
		{"1/s -> 1", func(s complex128) complex128 { return 1 / s }, 3, 1},
		{"1/(s+1) -> e^-t", func(s complex128) complex128 { return 1 / (s + 1) }, 1.5, math.Exp(-1.5)},
		{"1/(s^2+1) -> sin", func(s complex128) complex128 { return 1 / (s*s + 1) }, 2, math.Sin(2)},
		{"1/s^2 -> t", func(s complex128) complex128 { return 1 / (s * s) }, 2.5, 2.5},
	}
	for _, c := range cases {
		got := InverseLaplaceTalbot(c.F, c.t, 30)
		if math.Abs(got-c.want) > 1e-6 {
			t.Errorf("Talbot %s = %v, want %v", c.name, got, c.want)
		}
	}
}

func TestInverseLaplaceStehfest(t *testing.T) {
	got := InverseLaplaceStehfest(func(s float64) float64 { return 1 / (s + 1) }, 1, 12)
	if math.Abs(got-math.Exp(-1)) > 1e-4 {
		t.Errorf("Stehfest e^-t = %v, want %v", got, math.Exp(-1))
	}
	got = InverseLaplaceStehfest(func(s float64) float64 { return 1 / (s * s) }, 2, 12)
	if math.Abs(got-2) > 1e-4 {
		t.Errorf("Stehfest t = %v, want 2", got)
	}
}

func TestStehfestCoefficients(t *testing.T) {
	// The classical n=10 weights sum to zero (up to rounding).
	v := StehfestCoefficients(10)
	var sum float64
	for _, x := range v {
		sum += x
	}
	if math.Abs(sum) > 1e-6 {
		t.Errorf("Stehfest weights sum = %v, want ~0", sum)
	}
	if math.Abs(v[0]-1.0/12.0) > 1e-9 {
		t.Errorf("V_1 = %v, want 1/12", v[0])
	}
}

func TestInverseLaplaceEuler(t *testing.T) {
	got := InverseLaplaceEuler(func(s complex128) complex128 { return 1 / (s + 1) }, 1, 15)
	if math.Abs(got-math.Exp(-1)) > 1e-6 {
		t.Errorf("Euler e^-t = %v, want %v", got, math.Exp(-1))
	}
	got = InverseLaplaceEuler(func(s complex128) complex128 { return 1 / (s*s + 1) }, 2, 18)
	if math.Abs(got-math.Sin(2)) > 1e-6 {
		t.Errorf("Euler sin = %v, want %v", got, math.Sin(2))
	}
}

func TestZTransform(t *testing.T) {
	// For a finite geometric sequence a^n, X(z) = sum (a/z)^n.
	x := []float64{1, 2, 4, 8}
	z := complex(2, 0)
	got := ZTransform(x, z)
	// 1 + 2/2 + 4/4 + 8/8 = 4.
	if cmplx.Abs(got-4) > 1e-12 {
		t.Errorf("ZTransform = %v, want 4", got)
	}
}

func TestInverseZTransform(t *testing.T) {
	// X(z) = 1/(1 - 0.5 z^{-1}) is the transform of 0.5^n.
	X := func(z complex128) complex128 { return 1 / (1 - 0.5/z) }
	got := InverseZTransform(X, 6, 1.0, 512)
	for n := 0; n < 6; n++ {
		want := math.Pow(0.5, float64(n))
		if math.Abs(real(got[n])-want) > 1e-6 || math.Abs(imag(got[n])) > 1e-6 {
			t.Errorf("InverseZTransform[%d] = %v, want %v", n, got[n], want)
		}
	}
}

func TestChirpZTransformIsDFT(t *testing.T) {
	x := []complex128{1, 2, 3, 4, 5}
	n := len(x)
	w := cmplx.Rect(1, -2*math.Pi/float64(n))
	got := ChirpZTransform(x, n, w, complex(1, 0))
	if !cmplxClose(got, DFT(x), 1e-8) {
		t.Errorf("ChirpZTransform != DFT")
	}
}
