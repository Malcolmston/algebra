package transform

import (
	"math"
	"math/cmplx"
	"testing"
)

const tol = 1e-9

func cmplxClose(a, b []complex128, eps float64) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if cmplx.Abs(a[i]-b[i]) > eps {
			return false
		}
	}
	return true
}

func floatClose(a, b []float64, eps float64) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if math.Abs(a[i]-b[i]) > eps {
			return false
		}
	}
	return true
}

func TestIsAndNextPow2(t *testing.T) {
	for _, n := range []int{1, 2, 4, 8, 1024} {
		if !IsPow2(n) {
			t.Errorf("IsPow2(%d) = false", n)
		}
	}
	for _, n := range []int{0, 3, 6, 100} {
		if IsPow2(n) {
			t.Errorf("IsPow2(%d) = true", n)
		}
	}
	cases := map[int]int{0: 1, 1: 1, 2: 2, 3: 4, 5: 8, 17: 32, 1024: 1024}
	for in, want := range cases {
		if got := NextPow2(in); got != want {
			t.Errorf("NextPow2(%d) = %d, want %d", in, got, want)
		}
	}
}

func TestDFTKnown(t *testing.T) {
	// The DFT of a constant vector is a single non-zero DC bin.
	x := []complex128{1, 1, 1, 1}
	want := []complex128{4, 0, 0, 0}
	if !cmplxClose(DFT(x), want, tol) {
		t.Errorf("DFT const = %v", DFT(x))
	}
	// A pure impulse transforms to the all-ones spectrum.
	imp := []complex128{1, 0, 0, 0, 0, 0, 0, 0}
	got := DFT(imp)
	for k := range got {
		if cmplx.Abs(got[k]-1) > tol {
			t.Errorf("impulse DFT[%d] = %v", k, got[k])
		}
	}
}

func TestFFTMatchesDFT(t *testing.T) {
	x := []complex128{1, 2, 3, 4, 5, 6, 7, 8}
	if !cmplxClose(FFT(x), DFT(x), 1e-9) {
		t.Errorf("FFT != DFT")
	}
}

func TestFFTIFFTRoundTrip(t *testing.T) {
	x := []complex128{3, 1, 4, 1, 5, 9, 2, 6}
	if !cmplxClose(IFFT(FFT(x)), x, 1e-9) {
		t.Errorf("IFFT(FFT(x)) != x")
	}
}

func TestFFTPanicsNonPow2(t *testing.T) {
	defer func() {
		if recover() == nil {
			t.Errorf("FFT did not panic on non-power-of-two length")
		}
	}()
	FFT([]complex128{1, 2, 3})
}

func TestFFTAnyMatchesDFT(t *testing.T) {
	for _, n := range []int{5, 6, 7, 9, 12} {
		x := make([]complex128, n)
		for i := range x {
			x[i] = complex(float64(i*i%7)-3, float64(i%3))
		}
		if !cmplxClose(FFTAny(x), DFT(x), 1e-8) {
			t.Errorf("FFTAny != DFT for n=%d", n)
		}
		if !cmplxClose(IFFTAny(FFTAny(x)), x, 1e-8) {
			t.Errorf("IFFTAny round trip failed for n=%d", n)
		}
	}
}

func TestBluesteinMatchesDFT(t *testing.T) {
	x := []complex128{1, 2, 3, 4, 5}
	if !cmplxClose(Bluestein(x), DFT(x), 1e-8) {
		t.Errorf("Bluestein != DFT")
	}
}

func TestRFFTRoundTrip(t *testing.T) {
	for _, n := range []int{8, 9, 16, 15} {
		x := make([]float64, n)
		for i := range x {
			x[i] = math.Sin(float64(i)) + 0.5*float64(i%4)
		}
		half := RFFT(x)
		if len(half) != n/2+1 {
			t.Errorf("RFFT length = %d, want %d", len(half), n/2+1)
		}
		back := IRFFT(half, n)
		if !floatClose(back, x, 1e-8) {
			t.Errorf("IRFFT(RFFT(x)) != x for n=%d", n)
		}
	}
}

func TestFFTRealMatchesDFT(t *testing.T) {
	x := []float64{1, -2, 3, -4}
	c := []complex128{1, -2, 3, -4}
	if !cmplxClose(FFTReal(x), DFT(c), 1e-9) {
		t.Errorf("FFTReal != DFT")
	}
}

func TestFFTShiftRoundTrip(t *testing.T) {
	for _, n := range []int{4, 5} {
		x := make([]complex128, n)
		for i := range x {
			x[i] = complex(float64(i), 0)
		}
		if !cmplxClose(IFFTShift(FFTShift(x)), x, tol) {
			t.Errorf("shift round trip failed for n=%d", n)
		}
	}
	// Known ordering for an even length.
	got := FFTShift([]complex128{0, 1, 2, 3})
	want := []complex128{2, 3, 0, 1}
	if !cmplxClose(got, want, tol) {
		t.Errorf("FFTShift even = %v", got)
	}
}

func TestFFTFreq(t *testing.T) {
	got := FFTFreq(4, 1)
	want := []float64{0, 0.25, -0.5, -0.25}
	if !floatClose(got, want, tol) {
		t.Errorf("FFTFreq(4,1) = %v", got)
	}
	r := RFFTFreq(4, 1)
	if !floatClose(r, []float64{0, 0.25, 0.5}, tol) {
		t.Errorf("RFFTFreq(4,1) = %v", r)
	}
}

func TestFFT2DRoundTrip(t *testing.T) {
	m := [][]complex128{
		{1, 2, 3},
		{4, 5, 6},
		{7, 8, 9},
	}
	back := IFFT2D(FFT2D(m))
	for i := range m {
		if !cmplxClose(back[i], m[i], 1e-8) {
			t.Errorf("FFT2D round trip failed on row %d: %v", i, back[i])
		}
	}
}

func TestDFTMatrix(t *testing.T) {
	n := 4
	w := DFTMatrix(n)
	x := []complex128{1, 2, 3, 4}
	got := make([]complex128, n)
	for j := 0; j < n; j++ {
		for k := 0; k < n; k++ {
			got[j] += w[j][k] * x[k]
		}
	}
	if !cmplxClose(got, DFT(x), 1e-9) {
		t.Errorf("DFTMatrix product != DFT")
	}
}

func TestSpectralHelpers(t *testing.T) {
	X := []complex128{complex(3, 4), complex(0, -2)}
	if !floatClose(Magnitude(X), []float64{5, 2}, tol) {
		t.Errorf("Magnitude = %v", Magnitude(X))
	}
	if !floatClose(PowerSpectrum(X), []float64{25, 4}, tol) {
		t.Errorf("PowerSpectrum = %v", PowerSpectrum(X))
	}
	if math.Abs(Phase(X)[0]-math.Atan2(4, 3)) > tol {
		t.Errorf("Phase wrong")
	}
}

func TestGoertzelMatchesDFT(t *testing.T) {
	x := []float64{1, 2, 3, 4, 5, 0, -1, 2}
	D := FFTReal(x)
	for k := range x {
		if cmplx.Abs(Goertzel(x, k, len(x))-D[k]) > 1e-9 {
			t.Errorf("Goertzel[%d] = %v, want %v", k, Goertzel(x, k, len(x)), D[k])
		}
	}
}

func TestGoertzelPower(t *testing.T) {
	n := 32
	x := make([]float64, n)
	freq, fs := 4.0, float64(n)
	for i := range x {
		x[i] = math.Cos(2 * math.Pi * freq * float64(i) / fs)
	}
	p := GoertzelPower(x, freq, fs)
	// Compare with the DFT magnitude squared at the same bin.
	k := int(freq)
	ref := cmplx.Abs(FFTReal(x)[k])
	if math.Abs(math.Sqrt(p)-ref) > 1e-6 {
		t.Errorf("GoertzelPower sqrt=%v, DFT mag=%v", math.Sqrt(p), ref)
	}
}

func TestFFTPlan(t *testing.T) {
	p := NewFFTPlan(8)
	if p.Len() != 8 {
		t.Errorf("Len = %d", p.Len())
	}
	x := []complex128{1, 2, 3, 4, 5, 6, 7, 8}
	if !cmplxClose(p.Forward(x), FFT(x), 1e-9) {
		t.Errorf("FFTPlan.Forward != FFT")
	}
	if !cmplxClose(p.Inverse(p.Forward(x)), x, 1e-9) {
		t.Errorf("FFTPlan round trip failed")
	}
}

func BenchmarkFFT(b *testing.B) {
	n := 4096
	x := make([]complex128, n)
	for i := range x {
		x[i] = complex(math.Sin(float64(i)), math.Cos(float64(i)/3))
	}
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = FFT(x)
	}
}
