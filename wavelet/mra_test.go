package wavelet

import (
	"math"
	"testing"
)

func TestMaxLevel(t *testing.T) {
	cases := map[int]int{1: 0, 2: 1, 3: 0, 4: 2, 8: 3, 12: 2, 16: 4, 0: 0}
	for n, want := range cases {
		if got := MaxLevel(n); got != want {
			t.Errorf("MaxLevel(%d) = %d want %d", n, got, want)
		}
	}
}

func TestWaveDecRecRoundTrip(t *testing.T) {
	x := make([]float64, 64)
	for i := range x {
		x[i] = math.Sin(0.3*float64(i)) + 0.1*float64(i%5)
	}
	for _, w := range []Wavelet{Haar(), DB2(), DB4()} {
		for _, lvl := range []int{1, 2, 3} {
			dec, err := WaveDec(x, w, lvl)
			if err != nil {
				t.Fatalf("%s L%d: %v", w.Name(), lvl, err)
			}
			if dec.Levels() != lvl {
				t.Errorf("%s: Levels() = %d want %d", w.Name(), dec.Levels(), lvl)
			}
			r := WaveRec(dec)
			if !sliceApproxEqual(r, x, 1e-9) {
				t.Errorf("%s L%d: WaveRec round trip failed", w.Name(), lvl)
			}
		}
	}
}

func TestWaveDecErrors(t *testing.T) {
	x := make([]float64, 8)
	if _, err := WaveDec(x, Haar(), 0); err == nil {
		t.Error("levels=0 should error")
	}
	if _, err := WaveDec(x, Haar(), 4); err == nil {
		t.Error("levels beyond MaxLevel should error")
	}
}

func TestDecompositionEnergyAndFlatten(t *testing.T) {
	x := []float64{3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5, 8, 9, 7, 9, 3}
	dec, err := WaveDec(x, DB2(), 3)
	if err != nil {
		t.Fatal(err)
	}
	if !approxEqual(dec.Energy(), Energy(x), 1e-9) {
		t.Errorf("Energy = %v want %v", dec.Energy(), Energy(x))
	}
	if len(dec.Flatten()) != len(x) {
		t.Errorf("Flatten length = %d want %d", len(dec.Flatten()), len(x))
	}
}

func TestMRAComponentsSumToSignal(t *testing.T) {
	// The additive multiresolution components must sum back to the signal.
	x := make([]float64, 32)
	for i := range x {
		x[i] = float64((i*i)%13) - 6
	}
	for _, w := range []Wavelet{Haar(), DB2(), DB4()} {
		dec, err := WaveDec(x, w, 3)
		if err != nil {
			t.Fatal(err)
		}
		approx, details := dec.MRAComponents()
		sum := append([]float64(nil), approx...)
		for _, comp := range details {
			if len(comp) != len(x) {
				t.Fatalf("%s: component length = %d want %d", w.Name(), len(comp), len(x))
			}
			for i := range sum {
				sum[i] += comp[i]
			}
		}
		if !sliceApproxEqual(sum, x, 1e-9) {
			t.Errorf("%s: MRA components do not sum to signal", w.Name())
		}
	}
}

func TestDetailAt(t *testing.T) {
	x := make([]float64, 16)
	for i := range x {
		x[i] = float64(i)
	}
	dec, _ := WaveDec(x, Haar(), 2)
	if len(dec.DetailAt(1)) != 8 {
		t.Errorf("DetailAt(1) length = %d want 8", len(dec.DetailAt(1)))
	}
	if len(dec.DetailAt(2)) != 4 {
		t.Errorf("DetailAt(2) length = %d want 4", len(dec.DetailAt(2)))
	}
	defer func() {
		if recover() == nil {
			t.Error("DetailAt out of range should panic")
		}
	}()
	dec.DetailAt(3)
}
