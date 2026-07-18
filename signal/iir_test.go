package signal

import (
	"math"
	"testing"
)

// sosMagnitude returns the magnitude of a cascade of sections at a normalized
// frequency f (1.0 == Nyquist), for test assertions.
func sosMagnitude(sos []Biquad, f float64) float64 {
	m := 1.0
	for i := range sos {
		r := sos[i].FrequencyResponse([]float64{f})
		m *= cmag(r[0])
	}
	return m
}

func TestBiquadLowpassGains(t *testing.T) {
	bq := BiquadLowpass(1000, 8000, math.Sqrt2/2)
	// DC gain unity, Nyquist gain ~0.
	dc := bq.FrequencyResponse([]float64{0})
	if !approx(cmag(dc[0]), 1, 1e-9) {
		t.Errorf("low-pass DC gain = %v", cmag(dc[0]))
	}
	nyq := bq.FrequencyResponse([]float64{1})
	if cmag(nyq[0]) > 1e-6 {
		t.Errorf("low-pass Nyquist gain = %v, want ~0", cmag(nyq[0]))
	}
	// At the corner frequency (normalized 0.25 == 1000/4000) magnitude is 1/√2.
	corner := bq.FrequencyResponse([]float64{0.25})
	if !approx(cmag(corner[0]), math.Sqrt2/2, 1e-6) {
		t.Errorf("low-pass corner gain = %v, want %v", cmag(corner[0]), math.Sqrt2/2)
	}
}

func TestBiquadHighpassGains(t *testing.T) {
	bq := BiquadHighpass(1000, 8000, math.Sqrt2/2)
	dc := bq.FrequencyResponse([]float64{0})
	if cmag(dc[0]) > 1e-6 {
		t.Errorf("high-pass DC gain = %v, want ~0", cmag(dc[0]))
	}
	nyq := bq.FrequencyResponse([]float64{1})
	if !approx(cmag(nyq[0]), 1, 1e-9) {
		t.Errorf("high-pass Nyquist gain = %v, want 1", cmag(nyq[0]))
	}
}

func TestBiquadNotch(t *testing.T) {
	bq := BiquadNotch(1000, 8000, 5)
	// Zero at centre frequency (normalized 0.25).
	c := bq.FrequencyResponse([]float64{0.25})
	if cmag(c[0]) > 1e-6 {
		t.Errorf("notch centre gain = %v, want ~0", cmag(c[0]))
	}
	// Unity far away (DC).
	dc := bq.FrequencyResponse([]float64{0})
	if !approx(cmag(dc[0]), 1, 1e-9) {
		t.Errorf("notch DC gain = %v, want 1", cmag(dc[0]))
	}
}

func TestBiquadAllpassFlat(t *testing.T) {
	bq := BiquadAllpass(1200, 8000, 2)
	for _, f := range []float64{0, 0.1, 0.3, 0.5, 0.8, 1} {
		r := bq.FrequencyResponse([]float64{f})
		if !approx(cmag(r[0]), 1, 1e-9) {
			t.Errorf("all-pass magnitude at %v = %v, want 1", f, cmag(r[0]))
		}
	}
}

func TestBiquadBandpassPeak(t *testing.T) {
	bq := BiquadBandpass(1000, 8000, 4)
	// Unity peak at centre (normalized 0.25), zero at DC and Nyquist.
	c := bq.FrequencyResponse([]float64{0.25})
	if !approx(cmag(c[0]), 1, 1e-6) {
		t.Errorf("band-pass peak = %v, want 1", cmag(c[0]))
	}
	dc := bq.FrequencyResponse([]float64{0})
	if cmag(dc[0]) > 1e-9 {
		t.Errorf("band-pass DC gain = %v, want 0", cmag(dc[0]))
	}
}

func TestBiquadPeakingGain(t *testing.T) {
	bq := BiquadPeaking(1000, 8000, 2, 6)
	c := bq.FrequencyResponse([]float64{0.25})
	wantLin := math.Pow(10, 6.0/20)
	if !approx(cmag(c[0]), wantLin, 1e-6) {
		t.Errorf("peaking gain at centre = %v, want %v", cmag(c[0]), wantLin)
	}
	// Far from centre the gain returns to unity.
	dc := bq.FrequencyResponse([]float64{0})
	if !approx(cmag(dc[0]), 1, 1e-6) {
		t.Errorf("peaking DC gain = %v, want 1", cmag(dc[0]))
	}
}

func TestBiquadShelfGains(t *testing.T) {
	ls := BiquadLowShelf(1000, 8000, math.Sqrt2/2, 6)
	wantLin := math.Pow(10, 6.0/20)
	dc := ls.FrequencyResponse([]float64{0})
	if !approx(cmag(dc[0]), wantLin, 1e-6) {
		t.Errorf("low-shelf DC gain = %v, want %v", cmag(dc[0]), wantLin)
	}
	nyq := ls.FrequencyResponse([]float64{1})
	if !approx(cmag(nyq[0]), 1, 1e-6) {
		t.Errorf("low-shelf Nyquist gain = %v, want 1", cmag(nyq[0]))
	}

	hs := BiquadHighShelf(1000, 8000, math.Sqrt2/2, 6)
	hdc := hs.FrequencyResponse([]float64{0})
	if !approx(cmag(hdc[0]), 1, 1e-6) {
		t.Errorf("high-shelf DC gain = %v, want 1", cmag(hdc[0]))
	}
	hnyq := hs.FrequencyResponse([]float64{1})
	if !approx(cmag(hnyq[0]), wantLin, 1e-6) {
		t.Errorf("high-shelf Nyquist gain = %v, want %v", cmag(hnyq[0]), wantLin)
	}
}

func TestButterworthLowpassResponse(t *testing.T) {
	for _, order := range []int{1, 2, 3, 4, 5} {
		sos := ButterworthLowpass(order, 1000, 8000)
		// DC gain unity.
		if m := sosMagnitude(sos, 0); !approx(m, 1, 1e-6) {
			t.Errorf("order %d DC gain = %v, want 1", order, m)
		}
		// −3 dB at the corner (normalized 0.25 == 1000/4000).
		if m := sosMagnitude(sos, 0.25); !approx(m, math.Sqrt2/2, 1e-6) {
			t.Errorf("order %d corner gain = %v, want %v", order, m, math.Sqrt2/2)
		}
		// Monotonic decrease past the corner.
		if m1, m2 := sosMagnitude(sos, 0.3), sosMagnitude(sos, 0.6); m2 >= m1 {
			t.Errorf("order %d not monotonic: %v then %v", order, m1, m2)
		}
	}
}

func TestButterworthHighpassResponse(t *testing.T) {
	for _, order := range []int{2, 3, 4} {
		sos := ButterworthHighpass(order, 1000, 8000)
		if m := sosMagnitude(sos, 0); m > 1e-6 {
			t.Errorf("order %d DC gain = %v, want ~0", order, m)
		}
		if m := sosMagnitude(sos, 0.25); !approx(m, math.Sqrt2/2, 1e-6) {
			t.Errorf("order %d corner gain = %v, want %v", order, m, math.Sqrt2/2)
		}
		if m := sosMagnitude(sos, 1); !approx(m, 1, 1e-6) {
			t.Errorf("order %d Nyquist gain = %v, want 1", order, m)
		}
	}
}

func TestFilterSOSMatchesResponse(t *testing.T) {
	// Feeding a unit impulse through the cascade yields the impulse response;
	// its DFT magnitude must match the analytic frequency response.
	sos := ButterworthLowpass(3, 1000, 8000)
	imp := make([]float64, 128)
	imp[0] = 1
	h := FilterSOS(sos, imp)
	X := DFT(h)
	// Bin k corresponds to normalized frequency 2k/N.
	k := 16 // normalized 2*16/128 = 0.25 == corner
	f := 2 * float64(k) / float64(len(h))
	want := sosMagnitude(sos, f)
	if !approx(cmag(X[k]), want, 1e-6) {
		t.Errorf("impulse-response spectrum = %v, want %v", cmag(X[k]), want)
	}
}

func TestFilterSOSDoesNotMutate(t *testing.T) {
	sos := ButterworthLowpass(2, 1000, 8000)
	before := sos[0]
	FilterSOS(sos, []float64{1, 2, 3, 4})
	if sos[0] != before {
		t.Errorf("FilterSOS mutated caller's sections")
	}
}
