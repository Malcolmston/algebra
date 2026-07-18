package signal

import (
	"math"
	"testing"
)

func TestUpsample(t *testing.T) {
	got := Upsample([]float64{1, 2, 3}, 2)
	want := []float64{1, 0, 2, 0, 3, 0}
	if !approxSlice(got, want, tol) {
		t.Errorf("Upsample = %v, want %v", got, want)
	}
}

func TestDownsample(t *testing.T) {
	got := Downsample([]float64{1, 2, 3, 4, 5}, 2)
	want := []float64{1, 3, 5}
	if !approxSlice(got, want, tol) {
		t.Errorf("Downsample = %v, want %v", got, want)
	}
}

func TestUpDownRoundTrip(t *testing.T) {
	x := []float64{7, -3, 4, 9}
	if !approxSlice(Downsample(Upsample(x, 3), 3), x, tol) {
		t.Errorf("down(up(x)) != x")
	}
}

func TestResampleLinear(t *testing.T) {
	got := ResampleLinear([]float64{0, 10}, 5)
	want := []float64{0, 2.5, 5, 7.5, 10}
	if !approxSlice(got, want, tol) {
		t.Errorf("ResampleLinear = %v, want %v", got, want)
	}
	// Decimating positions of an exact linear ramp are exact.
	ramp := []float64{0, 1, 2, 3, 4}
	if !approxSlice(ResampleLinear(ramp, 3), []float64{0, 2, 4}, tol) {
		t.Errorf("ResampleLinear ramp wrong")
	}
}

func TestResampleLinearEdge(t *testing.T) {
	if got := ResampleLinear([]float64{5}, 4); !approxSlice(got, []float64{5, 5, 5, 5}, tol) {
		t.Errorf("single-sample resample = %v", got)
	}
	if len(ResampleLinear(nil, 3)) != 0 {
		t.Errorf("empty input should give empty output")
	}
}

func TestInterpolateConstant(t *testing.T) {
	// A constant signal must stay (approximately) constant after interpolation.
	x := make([]float64, 30)
	for i := range x {
		x[i] = 4
	}
	y := Interpolate(x, 2)
	if len(y) != len(x)*2 {
		t.Fatalf("interpolate length = %d, want %d", len(y), len(x)*2)
	}
	// Check the settled middle region (away from filter transients).
	for i := 30; i < len(y)-30; i++ {
		if !approx(y[i], 4, 0.05) {
			t.Errorf("interpolated constant off at %d: %v", i, y[i])
		}
	}
}

func TestDecimateConstant(t *testing.T) {
	x := make([]float64, 60)
	for i := range x {
		x[i] = 3
	}
	y := Decimate(x, 3)
	if len(y) != (len(x)+2)/3 {
		t.Fatalf("decimate length = %d", len(y))
	}
	for i := 5; i < len(y)-5; i++ {
		if !approx(y[i], 3, 0.05) {
			t.Errorf("decimated constant off at %d: %v", i, y[i])
		}
	}
}

func TestDecimateRemovesHighFrequency(t *testing.T) {
	// A tone above the post-decimation Nyquist should be strongly attenuated.
	n := 200
	x := make([]float64, n)
	for i := range x {
		// Near Nyquist of the original signal.
		x[i] = math.Cos(0.9 * math.Pi * float64(i))
	}
	y := Decimate(x, 2)
	// Compare energy in the settled region against the input amplitude.
	var e float64
	c := 0
	for i := 20; i < len(y)-20; i++ {
		e += y[i] * y[i]
		c++
	}
	if c > 0 {
		rms := math.Sqrt(e / float64(c))
		if rms > 0.2 {
			t.Errorf("high-frequency tone not attenuated: rms=%v", rms)
		}
	}
}
