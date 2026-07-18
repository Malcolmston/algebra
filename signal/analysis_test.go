package signal

import (
	"math"
	"testing"
)

func TestMovingAverage(t *testing.T) {
	got := MovingAverage([]float64{1, 2, 3, 4}, 2)
	want := []float64{1, 1.5, 2.5, 3.5}
	if !approxSlice(got, want, tol) {
		t.Errorf("MovingAverage = %v, want %v", got, want)
	}
}

func TestMovingAverageCentered(t *testing.T) {
	got := MovingAverageCentered([]float64{1, 2, 3, 4, 5}, 3)
	want := []float64{1.5, 2, 3, 4, 4.5}
	if !approxSlice(got, want, tol) {
		t.Errorf("MovingAverageCentered = %v, want %v", got, want)
	}
}

func TestExponentialMovingAverage(t *testing.T) {
	got := ExponentialMovingAverage([]float64{1, 2, 3}, 0.5)
	want := []float64{1, 1.5, 2.25}
	if !approxSlice(got, want, tol) {
		t.Errorf("ExponentialMovingAverage = %v, want %v", got, want)
	}
}

func TestCumulativeSum(t *testing.T) {
	got := CumulativeSum([]float64{1, 2, 3, 4})
	want := []float64{1, 3, 6, 10}
	if !approxSlice(got, want, tol) {
		t.Errorf("CumulativeSum = %v, want %v", got, want)
	}
}

func TestDiff(t *testing.T) {
	got := Diff([]float64{1, 4, 9, 16})
	want := []float64{3, 5, 7}
	if !approxSlice(got, want, tol) {
		t.Errorf("Diff = %v, want %v", got, want)
	}
	if len(Diff([]float64{5})) != 0 {
		t.Errorf("Diff of single element should be empty")
	}
}

func TestRMSEnergy(t *testing.T) {
	if got := RMS([]float64{3, 4}); !approx(got, math.Sqrt(12.5), tol) {
		t.Errorf("RMS = %v", got)
	}
	if got := Energy([]float64{3, 4}); !approx(got, 25, tol) {
		t.Errorf("Energy = %v", got)
	}
	if RMS(nil) != 0 || Energy(nil) != 0 {
		t.Errorf("RMS/Energy of empty should be 0")
	}
}

func TestZeroPad(t *testing.T) {
	got := ZeroPad([]float64{1, 2}, 4)
	if !approxSlice(got, []float64{1, 2, 0, 0}, tol) {
		t.Errorf("ZeroPad = %v", got)
	}
	// Truncation when n < len.
	if got := ZeroPad([]float64{1, 2, 3}, 2); !approxSlice(got, []float64{1, 2}, tol) {
		t.Errorf("ZeroPad truncate = %v", got)
	}
}

func TestCumSumDiffInverse(t *testing.T) {
	x := []float64{2, -1, 5, 3, 0, 7}
	// Diff of CumulativeSum recovers x[1:].
	d := Diff(CumulativeSum(x))
	if !approxSlice(d, x[1:], tol) {
		t.Errorf("Diff(CumSum(x)) = %v, want %v", d, x[1:])
	}
}

// BenchmarkWelchPSD exercises the heaviest routine in the package: Welch's
// method runs a direct O(N²) DFT over each overlapping segment.
func BenchmarkWelchPSD(b *testing.B) {
	n := 4096
	x := make([]float64, n)
	for i := range x {
		x[i] = math.Sin(0.11*float64(i)) + 0.3*math.Cos(0.37*float64(i))
	}
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		WelchPSD(x, 256, 128, float64(n))
	}
}
