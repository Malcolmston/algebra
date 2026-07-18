package wavelet

import (
	"math"
	"testing"
)

func TestPacketReconstruct(t *testing.T) {
	x := []float64{1, 2, 3, 4, 5, 6, 7, 8, 8, 7, 6, 5, 4, 3, 2, 1}
	for _, w := range []Wavelet{Haar(), DB2(), DB4()} {
		for _, lvl := range []int{1, 2, 3} {
			tree, err := WaveletPacket(x, w, lvl)
			if err != nil {
				t.Fatalf("%s L%d: %v", w.Name(), lvl, err)
			}
			if tree.Levels() != lvl {
				t.Errorf("%s: Levels() = %d want %d", w.Name(), tree.Levels(), lvl)
			}
			if tree.NodeCount(lvl) != 1<<lvl {
				t.Errorf("%s: NodeCount(%d) = %d want %d", w.Name(), lvl, tree.NodeCount(lvl), 1<<lvl)
			}
			r := tree.Reconstruct()
			if !sliceApproxEqual(r, x, 1e-9) {
				t.Errorf("%s L%d: packet reconstruction failed", w.Name(), lvl)
			}
		}
	}
}

func TestPacketEnergyConservation(t *testing.T) {
	x := []float64{2, -1, 3, 0, 4, -2, 1, 5}
	tree, _ := WaveletPacket(x, DB2(), 2)
	var e float64
	for _, leaf := range tree.Leaves() {
		e += Energy(leaf)
	}
	if !approxEqual(e, Energy(x), 1e-9) {
		t.Errorf("packet leaf energy = %v want %v", e, Energy(x))
	}
}

func TestPacketLeafLayout(t *testing.T) {
	x := make([]float64, 8)
	tree, _ := WaveletPacket(x, Haar(), 3)
	leaves := tree.Leaves()
	if len(leaves) != 8 {
		t.Fatalf("leaf count = %d want 8", len(leaves))
	}
	for i, leaf := range leaves {
		if len(leaf) != 1 {
			t.Errorf("leaf %d length = %d want 1", i, len(leaf))
		}
	}
	// Level-0 node equals the input.
	if !sliceApproxEqual(tree.Node(0, 0), x, 1e-12) {
		t.Error("root node should equal the input")
	}
}

func TestShannonEntropyKnownValues(t *testing.T) {
	// Fully concentrated energy -> entropy 0.
	if got := ShannonEntropy([]float64{5, 0, 0, 0}); !approxEqual(got, 0, 1e-12) {
		t.Errorf("ShannonEntropy(concentrated) = %v want 0", got)
	}
	// Uniform energy over n bins -> ln(n).
	if got := ShannonEntropy([]float64{1, 1, 1, 1}); !approxEqual(got, math.Log(4), 1e-12) {
		t.Errorf("ShannonEntropy(uniform4) = %v want ln4", got)
	}
	if got := ShannonEntropy([]float64{2, 2}); !approxEqual(got, math.Log(2), 1e-12) {
		t.Errorf("ShannonEntropy(uniform2) = %v want ln2", got)
	}
	// Zero-energy input.
	if got := ShannonEntropy([]float64{0, 0}); got != 0 {
		t.Errorf("ShannonEntropy(zero) = %v want 0", got)
	}
}

func TestLogEnergy(t *testing.T) {
	// ln(2^2) + ln(3^2) = 2 ln2 + 2 ln3, zeros skipped.
	got := LogEnergy([]float64{2, 0, 3})
	want := 2*math.Log(2) + 2*math.Log(3)
	if !approxEqual(got, want, 1e-12) {
		t.Errorf("LogEnergy = %v want %v", got, want)
	}
}
