package wavelet

import (
	"math"
	"testing"
)

func TestSoftThreshold(t *testing.T) {
	cases := []struct{ x, lambda, want float64 }{
		{5, 2, 3},
		{-5, 2, -3},
		{1, 2, 0},
		{-1, 2, 0},
		{2, 2, 0},
		{0, 1, 0},
	}
	for _, c := range cases {
		if got := SoftThreshold(c.x, c.lambda); !approxEqual(got, c.want, 1e-12) {
			t.Errorf("SoftThreshold(%v,%v) = %v want %v", c.x, c.lambda, got, c.want)
		}
	}
}

func TestHardThreshold(t *testing.T) {
	cases := []struct{ x, lambda, want float64 }{
		{5, 2, 5},
		{-5, 2, -5},
		{1, 2, 0},
		{2, 2, 0},
		{2.0001, 2, 2.0001},
	}
	for _, c := range cases {
		if got := HardThreshold(c.x, c.lambda); !approxEqual(got, c.want, 1e-12) {
			t.Errorf("HardThreshold(%v,%v) = %v want %v", c.x, c.lambda, got, c.want)
		}
	}
}

func TestThresholdModeDispatch(t *testing.T) {
	if got := Threshold(5, 2, Soft); !approxEqual(got, 3, 1e-12) {
		t.Errorf("Threshold Soft = %v want 3", got)
	}
	if got := Threshold(5, 2, Hard); !approxEqual(got, 5, 1e-12) {
		t.Errorf("Threshold Hard = %v want 5", got)
	}
	out := ThresholdSlice([]float64{5, 1, -4}, 2, Soft)
	want := []float64{3, 0, -2}
	if !sliceApproxEqual(out, want, 1e-12) {
		t.Errorf("ThresholdSlice = %v want %v", out, want)
	}
}

func TestUniversalThreshold(t *testing.T) {
	// sigma * sqrt(2 ln n).
	got := UniversalThreshold(16, 2)
	want := 2 * math.Sqrt(2*math.Log(16))
	if !approxEqual(got, want, 1e-12) {
		t.Errorf("UniversalThreshold = %v want %v", got, want)
	}
	if UniversalThreshold(1, 5) != 0 {
		t.Error("UniversalThreshold(1,..) should be 0")
	}
}

func TestEstimateNoiseSigma(t *testing.T) {
	// median(|d|) / 0.6745. For [1,-1,1,-1] median|.|=1.
	got := EstimateNoiseSigma([]float64{1, -1, 1, -1})
	want := 1.0 / 0.6745
	if !approxEqual(got, want, 1e-12) {
		t.Errorf("EstimateNoiseSigma = %v want %v", got, want)
	}
	if EstimateNoiseSigma(nil) != 0 {
		t.Error("EstimateNoiseSigma(nil) should be 0")
	}
}

func TestDenoiseConstantIsIdentity(t *testing.T) {
	// A noise-free constant signal has zero detail, so VisuShrink leaves it
	// unchanged.
	c := []float64{7, 7, 7, 7, 7, 7, 7, 7}
	for _, w := range []Wavelet{Haar(), DB2(), DB4()} {
		out := Denoise(c, w, 2, Soft)
		if !sliceApproxEqual(out, c, 1e-9) {
			t.Errorf("%s: Denoise(constant) = %v want unchanged", w.Name(), out)
		}
	}
}

func TestDenoiseReducesEnergyOfDetail(t *testing.T) {
	// A signal with a large smooth part and small oscillations: denoising must
	// not increase total energy and must remain finite/deterministic.
	x := make([]float64, 32)
	for i := range x {
		x[i] = 10 + 0.01*math.Sin(3*float64(i))
	}
	out := Denoise(x, DB2(), 3, Soft)
	if Energy(out) > Energy(x)+1e-6 {
		t.Errorf("denoised energy %v exceeds input energy %v", Energy(out), Energy(x))
	}
	if len(out) != len(x) {
		t.Errorf("denoised length = %d want %d", len(out), len(x))
	}
}
