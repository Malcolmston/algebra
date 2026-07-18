package signal

import (
	"math"
	"testing"
)

func TestSinc(t *testing.T) {
	cases := []struct{ x, want float64 }{
		{0, 1},
		{1, 0},
		{2, 0},
		{0.5, 2 / math.Pi},
		{-0.5, 2 / math.Pi},
	}
	for _, c := range cases {
		if got := Sinc(c.x); !approx(got, c.want, 1e-12) {
			t.Errorf("Sinc(%v) = %v, want %v", c.x, got, c.want)
		}
	}
}

func TestFIRLowpassUnitDCGain(t *testing.T) {
	taps := FIRLowpass(21, 0.3)
	var sum float64
	for _, h := range taps {
		sum += h
	}
	if !approx(sum, 1, 1e-12) {
		t.Errorf("low-pass DC gain = %v, want 1", sum)
	}
}

func TestFIRLinearPhaseSymmetry(t *testing.T) {
	taps := FIRLowpass(15, 0.25)
	for i := 0; i < len(taps)/2; i++ {
		if !approx(taps[i], taps[len(taps)-1-i], 1e-12) {
			t.Errorf("taps not symmetric at %d", i)
		}
	}
}

func TestFIRHighpassGains(t *testing.T) {
	taps := FIRHighpass(21, 0.3)
	// DC gain (sum of taps) must be ~0.
	var dc float64
	for _, h := range taps {
		dc += h
	}
	if !approx(dc, 0, 1e-9) {
		t.Errorf("high-pass DC gain = %v, want 0", dc)
	}
	// Gain at Nyquist is sum h[n]*(-1)^n and must be ~1.
	var nyq float64
	for n, h := range taps {
		if n%2 == 0 {
			nyq += h
		} else {
			nyq -= h
		}
	}
	// Windowed-sinc gain is close to but not exactly 1 at Nyquist.
	if !approx(nyq, 1, 5e-3) {
		t.Errorf("high-pass Nyquist gain = %v, want ~1", nyq)
	}
}

func TestFIRBandpassStopbandNulls(t *testing.T) {
	taps := FIRBandpass(41, 0.3, 0.6)
	resp := FIRFrequencyResponse(taps, []float64{0, 1})
	// Stop-band gain of a windowed-sinc design is small but not exactly zero.
	if m := cmag(resp[0]); m > 5e-3 {
		t.Errorf("band-pass DC gain = %v, want ~0", m)
	}
	if m := cmag(resp[1]); m > 5e-3 {
		t.Errorf("band-pass Nyquist gain = %v, want ~0", m)
	}
	// Mid-band gain should be close to 1.
	mid := FIRFrequencyResponse(taps, []float64{0.45})
	if m := cmag(mid[0]); m < 0.9 {
		t.Errorf("band-pass mid gain = %v, want ~1", m)
	}
}

func TestFIRBandstopGains(t *testing.T) {
	taps := FIRBandstop(41, 0.3, 0.6)
	resp := FIRFrequencyResponse(taps, []float64{0, 1})
	if m := cmag(resp[0]); m < 0.9 {
		t.Errorf("band-stop DC gain = %v, want ~1", m)
	}
	if m := cmag(resp[1]); m < 0.9 {
		t.Errorf("band-stop Nyquist gain = %v, want ~1", m)
	}
}

func TestFIRGroupDelayConstant(t *testing.T) {
	numtaps := 11
	taps := FIRLowpass(numtaps, 0.3)
	freqs := []float64{0.1, 0.2, 0.35}
	gd := FIRGroupDelay(taps, freqs)
	want := float64(numtaps-1) / 2
	for i, g := range gd {
		if !approx(g, want, 1e-6) {
			t.Errorf("group delay at %v = %v, want %v", freqs[i], g, want)
		}
	}
}

func TestFIRFrequencyResponseDC(t *testing.T) {
	taps := FIRLowpass(21, 0.3)
	resp := FIRFrequencyResponse(taps, []float64{0})
	if !approx(cmag(resp[0]), 1, 1e-9) {
		t.Errorf("low-pass DC response = %v", cmag(resp[0]))
	}
}

func TestApplyFIRMatchesConvolve(t *testing.T) {
	taps := []float64{0.5, 0.5}
	x := []float64{1, 2, 3, 4}
	got := ApplyFIR(taps, x)
	want := Convolve(x, taps)
	if !approxSlice(got, want, tol) {
		t.Errorf("ApplyFIR = %v, want %v", got, want)
	}
}

func TestFIRFilterStreamingMatchesBatch(t *testing.T) {
	taps := FIRLowpass(9, 0.4)
	x := make([]float64, 50)
	for i := range x {
		x[i] = math.Sin(0.3 * float64(i))
	}
	// Full batch via convolution, truncated to input length.
	full := Convolve(x, taps)
	f := NewFIRFilter(taps)
	stream := f.Process(x)
	for i := range x {
		if !approx(stream[i], full[i], 1e-9) {
			t.Fatalf("streaming vs batch differ at %d: %v vs %v", i, stream[i], full[i])
		}
	}
	// Splitting the block must give the same result as processing it whole.
	f.Reset()
	a := f.Process(x[:20])
	b := f.Process(x[20:])
	for i := 0; i < 20; i++ {
		if !approx(a[i], stream[i], 1e-12) {
			t.Fatalf("split part A differs at %d", i)
		}
	}
	for i := 0; i < 30; i++ {
		if !approx(b[i], stream[20+i], 1e-12) {
			t.Fatalf("split part B differs at %d", i)
		}
	}
}

func TestFIRHighpassEvenPanics(t *testing.T) {
	defer func() {
		if recover() == nil {
			t.Errorf("FIRHighpass with even numtaps should panic")
		}
	}()
	FIRHighpass(10, 0.3)
}
