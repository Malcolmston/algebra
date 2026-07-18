package infotheory

import (
	"testing"
)

func TestBSCCapacity(t *testing.T) {
	tests := []struct {
		p, want float64
	}{
		{0, 1},
		{0.5, 0},
		{1, 1}, // a channel that always flips is a relabelled noiseless channel
		{0.11, 0.500084041835472},
	}
	for _, tt := range tests {
		if got := BSCCapacity(tt.p); !approx(got, tt.want, 1e-9) {
			t.Errorf("BSCCapacity(%v) = %v, want %v", tt.p, got, tt.want)
		}
	}
}

func TestBECCapacity(t *testing.T) {
	tests := []struct {
		e, want float64
	}{
		{0, 1},
		{0.5, 0.5},
		{1, 0},
	}
	for _, tt := range tests {
		if got := BECCapacity(tt.e); !approx(got, tt.want, 1e-9) {
			t.Errorf("BECCapacity(%v) = %v, want %v", tt.e, got, tt.want)
		}
	}
}

func TestChannelMutualInformation(t *testing.T) {
	// BSC with uniform input: I = capacity = 1 - H_b(p).
	ch := BSC{Crossover: 0.11}.Channel()
	mi, err := ch.MutualInformation(UniformDistribution(2))
	if err != nil {
		t.Fatal(err)
	}
	if !approx(mi, BSCCapacity(0.11), 1e-9) {
		t.Errorf("BSC uniform MI = %v, want %v", mi, BSCCapacity(0.11))
	}
	// Noiseless identity channel with uniform input: I = 1 bit.
	id := Channel{Transition: [][]float64{{1, 0}, {0, 1}}}
	mi, _ = id.MutualInformation(UniformDistribution(2))
	if !approx(mi, 1, 1e-9) {
		t.Errorf("identity channel MI = %v, want 1", mi)
	}
	if _, err := ch.MutualInformation([]float64{1, 0, 0}); err == nil {
		t.Error("expected shape error for mismatched input")
	}
}

func TestBlahutArimotoBSC(t *testing.T) {
	ch := BSC{Crossover: 0.11}.Channel()
	c, px, err := BlahutArimoto(ch, 1e-12, 10000)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(c, BSCCapacity(0.11), 1e-6) {
		t.Errorf("Blahut-Arimoto BSC capacity = %v, want %v", c, BSCCapacity(0.11))
	}
	// Capacity of the BSC is achieved by a uniform input.
	if !approx(px[0], 0.5, 1e-4) || !approx(px[1], 0.5, 1e-4) {
		t.Errorf("optimal input = %v, want ~uniform", px)
	}
}

func TestBlahutArimotoBEC(t *testing.T) {
	ch := BEC{Erasure: 0.3}.Channel()
	c := ChannelCapacity(ch)
	if !approx(c, 0.7, 1e-6) {
		t.Errorf("Blahut-Arimoto BEC capacity = %v, want 0.7", c)
	}
}

func TestBlahutArimotoNoiseless(t *testing.T) {
	// A noiseless channel over 4 symbols has capacity log2(4) = 2 bits.
	id := Channel{Transition: [][]float64{
		{1, 0, 0, 0},
		{0, 1, 0, 0},
		{0, 0, 1, 0},
		{0, 0, 0, 1},
	}}
	if c := ChannelCapacity(id); !approx(c, 2, 1e-6) {
		t.Errorf("noiseless 4-ary capacity = %v, want 2", c)
	}
}

func BenchmarkBlahutArimoto(b *testing.B) {
	// A fixed, moderately sized channel so the benchmark is deterministic.
	const n = 8
	tr := make([][]float64, n)
	for i := range tr {
		tr[i] = make([]float64, n)
		for j := range tr[i] {
			if i == j {
				tr[i][j] = 0.9
			} else {
				tr[i][j] = 0.1 / float64(n-1)
			}
		}
	}
	ch := Channel{Transition: tr}
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_, _, _ = BlahutArimoto(ch, 1e-10, 1000)
	}
}
