package ntheory

import "testing"

func TestDivisorSigma(t *testing.T) {
	tests := []struct {
		k, n, want int64
	}{
		{0, 1, 1},
		{0, 12, 6},  // τ(12)
		{1, 6, 12},  // σ(6) = 1+2+3+6
		{1, 28, 56}, // perfect number
		{2, 6, 50},  // 1+4+9+36
		{3, 4, 73},  // 1+8+64
		{1, 1, 1},
		{2, 1, 1},
		{0, 0, 0},
		{5, 0, 0},
	}
	for _, tt := range tests {
		if got := DivisorSigma(tt.k, tt.n); got != tt.want {
			t.Errorf("DivisorSigma(%d, %d) = %d, want %d", tt.k, tt.n, got, tt.want)
		}
	}
}

func TestDivisorSigmaMatchesExisting(t *testing.T) {
	for n := int64(1); n <= 1000; n++ {
		if got, want := DivisorSigma(0, n), CountDivisors(n); got != want {
			t.Errorf("DivisorSigma(0, %d) = %d, want CountDivisors = %d", n, got, want)
		}
		if got, want := DivisorSigma(1, n), SumDivisors(n); got != want {
			t.Errorf("DivisorSigma(1, %d) = %d, want SumDivisors = %d", n, got, want)
		}
	}
}

func TestTotientSieve(t *testing.T) {
	const n = 2000
	phi := TotientSieve(n)
	if len(phi) != n+1 {
		t.Fatalf("TotientSieve length = %d, want %d", len(phi), n+1)
	}
	if phi[0] != 0 || phi[1] != 1 {
		t.Errorf("phi[0]=%d phi[1]=%d, want 0 and 1", phi[0], phi[1])
	}
	for i := int64(2); i <= n; i++ {
		if phi[i] != EulerPhi(i) {
			t.Errorf("TotientSieve[%d] = %d, want %d", i, phi[i], EulerPhi(i))
		}
	}
}

func TestMobiusSieve(t *testing.T) {
	const n = 2000
	mu := MobiusSieve(n)
	if len(mu) != n+1 {
		t.Fatalf("MobiusSieve length = %d, want %d", len(mu), n+1)
	}
	if mu[0] != 0 || mu[1] != 1 {
		t.Errorf("mu[0]=%d mu[1]=%d, want 0 and 1", mu[0], mu[1])
	}
	for i := int64(2); i <= n; i++ {
		if int(mu[i]) != MobiusMu(i) {
			t.Errorf("MobiusSieve[%d] = %d, want %d", i, mu[i], MobiusMu(i))
		}
	}
}

func TestMertensFunction(t *testing.T) {
	tests := []struct {
		n, want int64
	}{
		{0, 0},
		{1, 1},
		{2, 0},
		{3, -1},
		{4, -1},
		{10, -1},
		{100, 1},
		{1000, 2},
	}
	for _, tt := range tests {
		if got := MertensFunction(tt.n); got != tt.want {
			t.Errorf("MertensFunction(%d) = %d, want %d", tt.n, got, tt.want)
		}
	}
}

func BenchmarkTotientSieve(b *testing.B) {
	for i := 0; i < b.N; i++ {
		TotientSieve(100000)
	}
}

func BenchmarkMobiusSieve(b *testing.B) {
	for i := 0; i < b.N; i++ {
		MobiusSieve(100000)
	}
}
