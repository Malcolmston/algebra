package wavelet

import (
	"math"
	"testing"
)

func TestExtensions(t *testing.T) {
	x := []float64{1, 2, 3, 4}
	if got := PeriodicExtend(x, 2); !sliceApproxEqual(got, []float64{3, 4, 1, 2, 3, 4, 1, 2}, 1e-12) {
		t.Errorf("PeriodicExtend = %v", got)
	}
	if got := SymmetricExtend(x, 2); !sliceApproxEqual(got, []float64{2, 1, 1, 2, 3, 4, 4, 3}, 1e-12) {
		t.Errorf("SymmetricExtend = %v", got)
	}
	if got := ZeroExtend(x, 2); !sliceApproxEqual(got, []float64{0, 0, 1, 2, 3, 4, 0, 0}, 1e-12) {
		t.Errorf("ZeroExtend = %v", got)
	}
}

func TestConvolve(t *testing.T) {
	// [1,2,3] * [1,1] = [1,3,5,3].
	got := Convolve([]float64{1, 2, 3}, []float64{1, 1})
	if !sliceApproxEqual(got, []float64{1, 3, 5, 3}, 1e-12) {
		t.Errorf("Convolve = %v", got)
	}
	if len(Convolve(nil, []float64{1})) != 0 {
		t.Error("Convolve with empty input should be empty")
	}
}

func TestDownUpsample(t *testing.T) {
	if got := Downsample([]float64{1, 2, 3, 4, 5}, 2); !sliceApproxEqual(got, []float64{1, 3, 5}, 1e-12) {
		t.Errorf("Downsample = %v", got)
	}
	if got := Upsample([]float64{1, 2, 3}, 2); !sliceApproxEqual(got, []float64{1, 0, 2, 0, 3, 0}, 1e-12) {
		t.Errorf("Upsample = %v", got)
	}
}

func TestNorms(t *testing.T) {
	x := []float64{3, 4}
	if !approxEqual(Energy(x), 25, 1e-12) {
		t.Errorf("Energy = %v want 25", Energy(x))
	}
	if !approxEqual(L2Norm(x), 5, 1e-12) {
		t.Errorf("L2Norm = %v want 5", L2Norm(x))
	}
	if !approxEqual(MaxAbs([]float64{-7, 2, 5}), 7, 1e-12) {
		t.Errorf("MaxAbs = %v want 7", MaxAbs([]float64{-7, 2, 5}))
	}
}

func TestMedianAndMAD(t *testing.T) {
	if !approxEqual(Median([]float64{3, 1, 2}), 2, 1e-12) {
		t.Error("Median odd failed")
	}
	if !approxEqual(Median([]float64{4, 1, 2, 3}), 2.5, 1e-12) {
		t.Error("Median even failed")
	}
	if math.IsNaN(Median([]float64{5})) {
		t.Error("Median single should not be NaN")
	}
	// MAD of [1,2,3,4,5]: median=3, deviations=[2,1,0,1,2], median=1.
	if !approxEqual(MeanAbsoluteDeviation([]float64{1, 2, 3, 4, 5}), 1, 1e-12) {
		t.Errorf("MAD = %v want 1", MeanAbsoluteDeviation([]float64{1, 2, 3, 4, 5}))
	}
}

func TestPowersOfTwo(t *testing.T) {
	if !IsPowerOfTwo(16) || IsPowerOfTwo(12) || IsPowerOfTwo(0) {
		t.Error("IsPowerOfTwo failed")
	}
	if NextPowerOfTwo(5) != 8 || NextPowerOfTwo(8) != 8 || NextPowerOfTwo(1) != 1 {
		t.Error("NextPowerOfTwo failed")
	}
	padded := ZeroPadToPowerOfTwo([]float64{1, 2, 3})
	if len(padded) != 4 || padded[3] != 0 {
		t.Errorf("ZeroPadToPowerOfTwo = %v", padded)
	}
}

func BenchmarkWaveletPacketDB4(b *testing.B) {
	// Heaviest routine: full wavelet packet decomposition and reconstruction.
	n := 1024
	x := make([]float64, n)
	for i := range x {
		x[i] = math.Sin(0.05*float64(i)) + 0.3*math.Cos(0.11*float64(i))
	}
	w := DB4()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		tree, err := WaveletPacket(x, w, 6)
		if err != nil {
			b.Fatal(err)
		}
		_ = tree.Reconstruct()
	}
}
