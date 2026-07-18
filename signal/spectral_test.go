package signal

import (
	"math"
	"testing"
)

func TestDFTKnownValues(t *testing.T) {
	// DFT of a constant is a single DC spike.
	X := DFT([]float64{1, 1, 1, 1})
	if !approx(cmag(X[0]), 4, tol) {
		t.Errorf("DFT DC bin = %v, want 4", cmag(X[0]))
	}
	for k := 1; k < 4; k++ {
		if cmag(X[k]) > tol {
			t.Errorf("DFT bin %d = %v, want 0", k, cmag(X[k]))
		}
	}
	// DFT of a unit impulse is flat.
	Y := DFT([]float64{1, 0, 0, 0})
	for k := 0; k < 4; k++ {
		if !approx(cmag(Y[k]), 1, tol) {
			t.Errorf("impulse DFT bin %d = %v, want 1", k, cmag(Y[k]))
		}
	}
	// Alternating signal concentrates at Nyquist.
	Z := DFT([]float64{1, -1, 1, -1})
	if !approx(cmag(Z[2]), 4, tol) {
		t.Errorf("alternating DFT Nyquist bin = %v, want 4", cmag(Z[2]))
	}
}

func TestIDFTRoundTrip(t *testing.T) {
	x := []float64{3, 1, 4, 1, 5, 9, 2, 6}
	rec := IDFT(DFT(x))
	for i := range x {
		if !approx(real(rec[i]), x[i], 1e-9) || !approx(imag(rec[i]), 0, 1e-9) {
			t.Errorf("IDFT(DFT(x))[%d] = %v, want %v", i, rec[i], x[i])
		}
	}
}

func TestMagnitudePhasePower(t *testing.T) {
	X := []complex128{complex(3, 4), complex(0, -2), complex(-1, 0)}
	mag := Magnitude(X)
	if !approxSlice(mag, []float64{5, 2, 1}, tol) {
		t.Errorf("Magnitude = %v", mag)
	}
	pw := PowerSpectrum(X)
	if !approxSlice(pw, []float64{25, 4, 1}, tol) {
		t.Errorf("PowerSpectrum = %v", pw)
	}
	ph := Phase(X)
	if !approx(ph[0], math.Atan2(4, 3), tol) || !approx(ph[1], -math.Pi/2, tol) || !approx(ph[2], math.Pi, tol) {
		t.Errorf("Phase = %v", ph)
	}
}

func TestFrequencyBins(t *testing.T) {
	got := FrequencyBins(4, 8)
	want := []float64{0, 2, 4, 6}
	if !approxSlice(got, want, tol) {
		t.Errorf("FrequencyBins = %v, want %v", got, want)
	}
}

func TestPeriodogramPeakLocation(t *testing.T) {
	n := 32
	fs := 32.0
	k0 := 5 // integer bin -> exactly on-grid tone
	x := make([]float64, n)
	for i := range x {
		x[i] = math.Cos(2 * math.Pi * float64(k0) * float64(i) / float64(n))
	}
	freqs, psd := Periodogram(x, nil, fs)
	if len(freqs) != n/2+1 {
		t.Fatalf("periodogram length = %d, want %d", len(freqs), n/2+1)
	}
	// Peak must be at bin k0.
	maxi := 0
	for i, p := range psd {
		if p > psd[maxi] {
			maxi = i
		}
	}
	if maxi != k0 {
		t.Errorf("periodogram peak at bin %d, want %d", maxi, k0)
	}
	if !approx(freqs[k0], float64(k0)*fs/float64(n), tol) {
		t.Errorf("peak frequency = %v, want %v", freqs[k0], float64(k0)*fs/float64(n))
	}
}

func TestPeriodogramParseval(t *testing.T) {
	// Integral of the one-sided PSD (rectangular window) times bin spacing
	// equals the mean-square of the signal (Parseval's theorem).
	n := 64
	fs := float64(n)
	x := make([]float64, n)
	for i := range x {
		x[i] = math.Sin(2*math.Pi*3*float64(i)/float64(n)) + 0.5*math.Cos(2*math.Pi*7*float64(i)/float64(n))
	}
	_, psd := Periodogram(x, nil, fs)
	df := fs / float64(n)
	var integral float64
	for _, p := range psd {
		integral += p * df
	}
	var ms float64
	for _, v := range x {
		ms += v * v
	}
	ms /= float64(n)
	if !approx(integral, ms, 1e-6) {
		t.Errorf("PSD integral = %v, mean-square = %v", integral, ms)
	}
}

func TestWelchPSDPeak(t *testing.T) {
	n := 256
	fs := 256.0
	k0 := 20
	x := make([]float64, n)
	for i := range x {
		x[i] = math.Cos(2 * math.Pi * float64(k0) * float64(i) / float64(n))
	}
	segLen := 64
	freqs, psd := WelchPSD(x, segLen, 32, fs)
	if len(freqs) != segLen/2+1 {
		t.Fatalf("welch length = %d", len(freqs))
	}
	// Peak frequency should be near k0 (== bin k0*segLen/n = 5 in segment bins).
	maxi := 0
	for i, p := range psd {
		if p > psd[maxi] {
			maxi = i
		}
	}
	peakFreq := freqs[maxi]
	if math.Abs(peakFreq-float64(k0)) > fs/float64(segLen) {
		t.Errorf("Welch peak at %v Hz, want near %v", peakFreq, float64(k0))
	}
}
