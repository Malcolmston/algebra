package infotheory

import (
	"math"
	"testing"
)

const tol = 1e-9

func approx(a, b, eps float64) bool {
	if math.IsInf(a, 0) || math.IsInf(b, 0) {
		return a == b
	}
	return math.Abs(a-b) <= eps
}

func TestEntropy(t *testing.T) {
	tests := []struct {
		name string
		p    []float64
		want float64
	}{
		{"fair coin", []float64{0.5, 0.5}, 1},
		{"uniform 4", []float64{0.25, 0.25, 0.25, 0.25}, 2},
		{"uniform 8", UniformDistribution(8), 3},
		{"deterministic", []float64{1, 0, 0}, 0},
		{"biased", []float64{0.9, 0.1}, 0.4689955935892812},
	}
	for _, tt := range tests {
		if got := Entropy(tt.p); !approx(got, tt.want, 1e-9) {
			t.Errorf("Entropy(%s) = %v, want %v", tt.name, got, tt.want)
		}
	}
}

func TestEntropyNatAndBase(t *testing.T) {
	p := []float64{0.5, 0.5}
	if got := EntropyNat(p); !approx(got, math.Ln2, tol) {
		t.Errorf("EntropyNat fair coin = %v, want ln2 = %v", got, math.Ln2)
	}
	if got := EntropyBase(p, 2); !approx(got, 1, tol) {
		t.Errorf("EntropyBase base 2 = %v, want 1", got)
	}
	if got := EntropyBase(UniformDistribution(10), 10); !approx(got, 1, tol) {
		t.Errorf("EntropyBase base 10 uniform-10 = %v, want 1", got)
	}
}

func TestBinaryEntropy(t *testing.T) {
	tests := []struct {
		p, want float64
	}{
		{0, 0},
		{1, 0},
		{0.5, 1},
		{0.11, 0.499915958164528},
	}
	for _, tt := range tests {
		if got := BinaryEntropy(tt.p); !approx(got, tt.want, 1e-9) {
			t.Errorf("BinaryEntropy(%v) = %v, want %v", tt.p, got, tt.want)
		}
	}
}

func TestSurprisalAndPerplexity(t *testing.T) {
	if got := Surprisal(0.5); !approx(got, 1, tol) {
		t.Errorf("Surprisal(0.5) = %v, want 1", got)
	}
	if got := Surprisal(0.25); !approx(got, 2, tol) {
		t.Errorf("Surprisal(0.25) = %v, want 2", got)
	}
	if got := Surprisal(0); !math.IsInf(got, 1) {
		t.Errorf("Surprisal(0) = %v, want +Inf", got)
	}
	if got := SurprisalNat(math.Exp(-1)); !approx(got, 1, tol) {
		t.Errorf("SurprisalNat(1/e) = %v, want 1", got)
	}
	if got := Perplexity(UniformDistribution(5)); !approx(got, 5, 1e-9) {
		t.Errorf("Perplexity(uniform 5) = %v, want 5", got)
	}
}

func TestNormalizedEntropyAndRedundancy(t *testing.T) {
	p := UniformDistribution(4)
	if got := NormalizedEntropy(p); !approx(got, 1, tol) {
		t.Errorf("NormalizedEntropy(uniform) = %v, want 1", got)
	}
	if got := Redundancy(p); !approx(got, 0, tol) {
		t.Errorf("Redundancy(uniform) = %v, want 0", got)
	}
	if got := NormalizedEntropy([]float64{1, 0}); !approx(got, 0, tol) {
		t.Errorf("NormalizedEntropy(point mass) = %v, want 0", got)
	}
}

func TestGiniImpurity(t *testing.T) {
	if got := GiniImpurity([]float64{0.5, 0.5}); !approx(got, 0.5, tol) {
		t.Errorf("Gini fair coin = %v, want 0.5", got)
	}
	if got := GiniImpurity(UniformDistribution(4)); !approx(got, 0.75, tol) {
		t.Errorf("Gini uniform-4 = %v, want 0.75", got)
	}
	if got := GiniImpurity([]float64{1, 0}); !approx(got, 0, tol) {
		t.Errorf("Gini point mass = %v, want 0", got)
	}
}

func TestRenyiFamily(t *testing.T) {
	u4 := UniformDistribution(4)
	// For the uniform distribution every Renyi order equals log2(n).
	for _, a := range []float64{0, 0.5, 1, 2, math.Inf(1)} {
		if got := RenyiEntropy(u4, a); !approx(got, 2, 1e-9) {
			t.Errorf("RenyiEntropy(uniform4, %v) = %v, want 2", a, got)
		}
	}
	p := []float64{0.5, 0.5}
	if got := CollisionEntropy(p); !approx(got, 1, tol) {
		t.Errorf("CollisionEntropy fair coin = %v, want 1", got)
	}
	if got := MinEntropy([]float64{0.9, 0.1}); !approx(got, -math.Log2(0.9), tol) {
		t.Errorf("MinEntropy = %v, want %v", got, -math.Log2(0.9))
	}
	if got := HartleyEntropy([]float64{0.5, 0.3, 0.2, 0}); !approx(got, math.Log2(3), tol) {
		t.Errorf("HartleyEntropy = %v, want log2(3)", got)
	}
	// Renyi(2) of a non-uniform distribution equals CollisionEntropy.
	q := []float64{0.7, 0.2, 0.1}
	if got, want := RenyiEntropy(q, 2), CollisionEntropy(q); !approx(got, want, tol) {
		t.Errorf("RenyiEntropy(q,2) = %v, want CollisionEntropy = %v", got, want)
	}
}

func TestTsallisEntropy(t *testing.T) {
	// Tsallis at q=1 reduces to Shannon entropy in nats.
	p := []float64{0.5, 0.5}
	if got := TsallisEntropy(p, 1); !approx(got, math.Ln2, tol) {
		t.Errorf("Tsallis(q=1) fair coin = %v, want ln2", got)
	}
	// Tsallis at q=2 equals the Gini impurity.
	q := []float64{0.7, 0.2, 0.1}
	if got, want := TsallisEntropy(q, 2), GiniImpurity(q); !approx(got, want, tol) {
		t.Errorf("Tsallis(q=2) = %v, want Gini = %v", got, want)
	}
}

func TestBitsNatsRoundTrip(t *testing.T) {
	if got := NatsToBits(BitsToNats(3.5)); !approx(got, 3.5, tol) {
		t.Errorf("round trip = %v, want 3.5", got)
	}
}
