package signal

import "testing"

func TestConvolve(t *testing.T) {
	a := []float64{1, 2, 3}
	b := []float64{0, 1, 0.5}
	got := Convolve(a, b)
	want := []float64{0, 1, 2.5, 4, 1.5}
	if !approxSlice(got, want, tol) {
		t.Errorf("Convolve = %v, want %v", got, want)
	}
}

func TestConvolveIdentity(t *testing.T) {
	a := []float64{5, -2, 7, 3}
	got := Convolve(a, []float64{1})
	if !approxSlice(got, a, tol) {
		t.Errorf("convolution with unit impulse = %v, want %v", got, a)
	}
}

func TestConvolveSame(t *testing.T) {
	got := ConvolveSame([]float64{1, 2, 3}, []float64{0, 1, 0.5})
	want := []float64{1, 2.5, 4}
	if !approxSlice(got, want, tol) {
		t.Errorf("ConvolveSame = %v, want %v", got, want)
	}
}

func TestConvolveValid(t *testing.T) {
	got := ConvolveValid([]float64{1, 2, 3}, []float64{0, 1, 0.5})
	want := []float64{2.5}
	if !approxSlice(got, want, tol) {
		t.Errorf("ConvolveValid = %v, want %v", got, want)
	}
}

func TestCrossCorrelate(t *testing.T) {
	a := []float64{1, 2, 3}
	b := []float64{0, 1, 0.5}
	got := CrossCorrelate(a, b)
	want := []float64{0.5, 2, 3.5, 3, 0}
	if !approxSlice(got, want, tol) {
		t.Errorf("CrossCorrelate = %v, want %v", got, want)
	}
	// Zero-lag element equals the inner product.
	if !approx(got[len(b)-1], 1*0+2*1+3*0.5, tol) {
		t.Errorf("zero-lag correlation wrong: %v", got[len(b)-1])
	}
}

func TestAutoCorrelate(t *testing.T) {
	a := []float64{1, 2, 3}
	got := AutoCorrelate(a)
	want := []float64{3, 8, 14, 8, 3}
	if !approxSlice(got, want, tol) {
		t.Errorf("AutoCorrelate = %v, want %v", got, want)
	}
	// Centre equals the signal energy and is the maximum.
	if !approx(got[2], Energy(a), tol) {
		t.Errorf("autocorrelation peak != energy")
	}
}

func TestConvolveEmpty(t *testing.T) {
	if len(Convolve(nil, []float64{1, 2})) != 0 {
		t.Errorf("convolution with empty should be empty")
	}
}
