package transform

import (
	"math"
	"math/cmplx"
	"testing"
)

func TestHilbertCosine(t *testing.T) {
	// The Hilbert transform of cos(w n) is sin(w n); the analytic envelope is
	// unity for a pure tone.
	n := 128
	k := 8.0
	sig := make([]float64, n)
	for i := 0; i < n; i++ {
		sig[i] = math.Cos(2 * math.Pi * k * float64(i) / float64(n))
	}
	ht := HilbertTransform(sig)
	env := Envelope(sig)
	// Check interior samples to avoid edge transients.
	for i := 16; i < n-16; i++ {
		want := math.Sin(2 * math.Pi * k * float64(i) / float64(n))
		if math.Abs(ht[i]-want) > 1e-6 {
			t.Fatalf("HilbertTransform[%d] = %v, want %v", i, ht[i], want)
		}
		if math.Abs(env[i]-1) > 1e-6 {
			t.Fatalf("Envelope[%d] = %v, want 1", i, env[i])
		}
	}
}

func TestHilbertRealPart(t *testing.T) {
	// The real part of the analytic signal equals the input.
	sig := []float64{1, -2, 3, 0, -1, 4, 2, -3}
	a := Hilbert(sig)
	for i := range sig {
		if math.Abs(real(a[i])-sig[i]) > 1e-9 {
			t.Errorf("Re(Hilbert)[%d] = %v, want %v", i, real(a[i]), sig[i])
		}
	}
}

func TestInstantaneousFrequency(t *testing.T) {
	// A complex-exponential-like real tone has a constant instantaneous
	// frequency equal to its tone frequency.
	n := 256
	fs := 256.0
	freq := 10.0
	sig := make([]float64, n)
	for i := 0; i < n; i++ {
		sig[i] = math.Cos(2 * math.Pi * freq * float64(i) / fs)
	}
	inst := InstantaneousFrequency(sig, fs)
	// Average over the interior where edge effects are small.
	var sum float64
	cnt := 0
	for i := 40; i < n-40; i++ {
		sum += inst[i]
		cnt++
	}
	avg := sum / float64(cnt)
	if math.Abs(avg-freq) > 0.5 {
		t.Errorf("mean instantaneous frequency = %v, want ~%v", avg, freq)
	}
}

func TestPhaseUnwrap(t *testing.T) {
	wrapped := []float64{0, 3, -3, 3, -3}
	got := PhaseUnwrap(wrapped)
	// Each step of ~ -6 (crossing 2pi) is corrected to a continuous ramp.
	for i := 1; i < len(got); i++ {
		if math.Abs((got[i]-got[i-1])-(3-(3-2*math.Pi))) > 1 {
			// Just check differences are bounded by pi after unwrapping.
		}
		if math.Abs(got[i]-got[i-1]) > math.Pi+1e-9 {
			t.Errorf("PhaseUnwrap step %d too large: %v", i, got[i]-got[i-1])
		}
	}
}

func TestDTFTImpulseAndDelta(t *testing.T) {
	// The DTFT of a unit impulse is 1 at every frequency.
	imp := []float64{1, 0, 0, 0}
	for _, w := range []float64{0, 0.5, 1.3, math.Pi} {
		if cmplx.Abs(DTFT(imp, w)-1) > 1e-12 {
			t.Errorf("DTFT(impulse, %v) = %v", w, DTFT(imp, w))
		}
	}
	// The DTFT at omega=0 equals the sum of the samples.
	x := []float64{1, 2, 3, 4}
	if cmplx.Abs(DTFT(x, 0)-10) > 1e-12 {
		t.Errorf("DTFT(x,0) = %v, want 10", DTFT(x, 0))
	}
}

func TestSampleDTFTMatchesDFT(t *testing.T) {
	x := []float64{1, 2, 3, 4, 5, 6}
	got := SampleDTFT(x, len(x))
	if !cmplxClose(got, FFTReal(x), 1e-9) {
		t.Errorf("SampleDTFT != DFT")
	}
	// DTFTSample at explicit frequencies matches DTFT.
	omegas := []float64{0.1, 0.7, 2.0}
	s := DTFTSample(x, omegas)
	for i, w := range omegas {
		if cmplx.Abs(s[i]-DTFT(x, w)) > 1e-12 {
			t.Errorf("DTFTSample[%d] mismatch", i)
		}
	}
}

func TestDTFTComplex(t *testing.T) {
	x := []complex128{1, 1i, -1, -1i}
	if cmplx.Abs(DTFTComplex(x, 0)-0) > 1e-12 {
		t.Errorf("DTFTComplex sum = %v, want 0", DTFTComplex(x, 0))
	}
}

func TestWindows(t *testing.T) {
	n := 16
	// The Hann and Bartlett windows vanish at both ends.
	for _, w := range [][]float64{Hann(n), Bartlett(n)} {
		if math.Abs(w[0]) > 1e-12 || math.Abs(w[n-1]) > 1e-12 {
			t.Errorf("window endpoints nonzero")
		}
	}
	// All standard windows are symmetric.
	for _, w := range [][]float64{Hann(n), Hamming(n), Blackman(n), Bartlett(n), Welch(n)} {
		for i := 0; i < n/2; i++ {
			if math.Abs(w[i]-w[n-1-i]) > 1e-12 {
				t.Errorf("window not symmetric")
				break
			}
		}
	}
	// Hamming endpoints are 0.08.
	h := Hamming(n)
	if math.Abs(h[0]-0.08) > 1e-12 {
		t.Errorf("Hamming[0] = %v, want 0.08", h[0])
	}
	// ApplyWindow multiplies element-wise.
	got := ApplyWindow([]float64{2, 4}, []float64{0.5, 0.25})
	if !floatClose(got, []float64{1, 1}, 1e-12) {
		t.Errorf("ApplyWindow = %v", got)
	}
}

func TestPeriodogram(t *testing.T) {
	// A cosine at bin k concentrates its power in bins k and N-k.
	n := 16
	k := 3
	x := make([]float64, n)
	for i := range x {
		x[i] = math.Cos(2 * math.Pi * float64(k) * float64(i) / float64(n))
	}
	p := Periodogram(x)
	// Bins k and N-k should each hold N/4 of the power.
	if math.Abs(p[k]-float64(n)/4) > 1e-6 {
		t.Errorf("Periodogram[%d] = %v, want %v", k, p[k], float64(n)/4)
	}
	if math.Abs(p[n-k]-float64(n)/4) > 1e-6 {
		t.Errorf("Periodogram[%d] = %v, want %v", n-k, p[n-k], float64(n)/4)
	}
}
