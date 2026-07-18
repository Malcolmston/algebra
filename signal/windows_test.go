package signal

import (
	"math"
	"testing"
)

func TestWindowsKnownValues(t *testing.T) {
	cases := []struct {
		name string
		got  []float64
		want []float64
	}{
		{"Rectangular4", Rectangular(4), []float64{1, 1, 1, 1}},
		{"Hann5", Hann(5), []float64{0, 0.5, 1, 0.5, 0}},
		{"Hamming5", Hamming(5), []float64{0.08, 0.54, 1.0, 0.54, 0.08}},
		{"Blackman5", Blackman(5), []float64{0, 0.34, 1.0, 0.34, 0}},
		{"Bartlett5", Bartlett(5), []float64{0, 0.5, 1, 0.5, 0}},
		{"Welch5", Welch(5), []float64{0, 0.75, 1, 0.75, 0}},
		{"Cosine5", Cosine(5), []float64{0, math.Sqrt2 / 2, 1, math.Sqrt2 / 2, 0}},
		{"Kaiser5beta0", Kaiser(5, 0), []float64{1, 1, 1, 1, 1}},
		{"Tukey9alpha0", Tukey(9, 0), Rectangular(9)},
	}
	for _, c := range cases {
		if !approxSlice(c.got, c.want, 1e-12) {
			t.Errorf("%s = %v, want %v", c.name, c.got, c.want)
		}
	}
}

func TestWindowsSymmetry(t *testing.T) {
	fns := map[string]func(int) []float64{
		"Hann": Hann, "Hamming": Hamming, "Blackman": Blackman,
		"BlackmanHarris": BlackmanHarris, "Nuttall": Nuttall, "FlatTop": FlatTop,
		"Bartlett": Bartlett, "Welch": Welch, "Cosine": Cosine,
	}
	for name, fn := range fns {
		w := fn(17)
		for i := 0; i < len(w)/2; i++ {
			if !approx(w[i], w[len(w)-1-i], 1e-12) {
				t.Errorf("%s not symmetric at %d: %v vs %v", name, i, w[i], w[len(w)-1-i])
			}
		}
	}
}

func TestTukeyEqualsHann(t *testing.T) {
	// Tukey with alpha == 1 is exactly the Hann window.
	if !approxSlice(Tukey(11, 1), Hann(11), 1e-12) {
		t.Errorf("Tukey(11,1) != Hann(11)")
	}
}

func TestBesselI0(t *testing.T) {
	cases := []struct{ x, want float64 }{
		{0, 1},
		{1, 1.2660658777520084},
		{2, 2.2795853023360673},
		{5, 27.239871823604442},
	}
	for _, c := range cases {
		if got := BesselI0(c.x); !approx(got, c.want, 1e-9) {
			t.Errorf("BesselI0(%v) = %v, want %v", c.x, got, c.want)
		}
	}
}

func TestKaiserBeta(t *testing.T) {
	if got := KaiserBeta(60); !approx(got, 0.1102*(60-8.7), 1e-12) {
		t.Errorf("KaiserBeta(60) = %v", got)
	}
	if got := KaiserBeta(10); got != 0 {
		t.Errorf("KaiserBeta(10) = %v, want 0", got)
	}
	// beta must be strictly positive in the 21..50 band.
	if got := KaiserBeta(35); got <= 0 {
		t.Errorf("KaiserBeta(35) = %v, want > 0", got)
	}
}

func TestApplyWindow(t *testing.T) {
	x := []float64{2, 4, 6}
	w := []float64{0.5, 1, 0}
	got := ApplyWindow(x, w)
	if !approxSlice(got, []float64{1, 4, 0}, tol) {
		t.Errorf("ApplyWindow = %v", got)
	}
}

func TestWindowEdgeCases(t *testing.T) {
	if len(Hann(0)) != 0 {
		t.Errorf("Hann(0) should be empty")
	}
	if got := Hann(1); len(got) != 1 || got[0] != 1 {
		t.Errorf("Hann(1) = %v", got)
	}
}
