package infotheory

import "testing"

func TestKraft(t *testing.T) {
	tests := []struct {
		name           string
		lengths        []int
		base           int
		sum            float64
		ineq, equality bool
	}{
		{"complete binary", []int{1, 2, 3, 3}, 2, 1, true, true},
		{"uniform 4x2", []int{2, 2, 2, 2}, 2, 1, true, true},
		{"incomplete", []int{2, 2, 3}, 2, 0.625, true, false},
		{"over-complete", []int{1, 1, 1}, 2, 1.5, false, false},
		{"ternary", []int{1, 1, 1}, 3, 1, true, true},
	}
	for _, tt := range tests {
		if got := KraftSum(tt.lengths, tt.base); !approx(got, tt.sum, 1e-9) {
			t.Errorf("%s: KraftSum = %v, want %v", tt.name, got, tt.sum)
		}
		if got := KraftInequality(tt.lengths, tt.base); got != tt.ineq {
			t.Errorf("%s: KraftInequality = %v, want %v", tt.name, got, tt.ineq)
		}
		if got := KraftEquality(tt.lengths, tt.base); got != tt.equality {
			t.Errorf("%s: KraftEquality = %v, want %v", tt.name, got, tt.equality)
		}
		if got := McMillanInequality(tt.lengths, tt.base); got != tt.ineq {
			t.Errorf("%s: McMillanInequality = %v, want %v", tt.name, got, tt.ineq)
		}
	}
}

func TestIsPrefixFree(t *testing.T) {
	if !IsPrefixFree([]string{"0", "10", "110", "111"}) {
		t.Error("expected prefix-free")
	}
	if IsPrefixFree([]string{"0", "01", "10"}) {
		t.Error("0 is a prefix of 01; expected not prefix-free")
	}
	if IsPrefixFree([]string{"1", "1"}) {
		t.Error("duplicate codewords are not prefix-free")
	}
}

func TestFanoInequalityBound(t *testing.T) {
	// Zero error probability gives a zero bound.
	if got := FanoInequalityBound(0, 4); !approx(got, 0, tol) {
		t.Errorf("Fano bound at Pe=0 = %v, want 0", got)
	}
	// Binary alphabet: bound reduces to H_b(Pe).
	if got := FanoInequalityBound(0.1, 2); !approx(got, BinaryEntropy(0.1), 1e-9) {
		t.Errorf("Fano bound binary = %v, want H_b(0.1)", got)
	}
	// General: H_b(Pe) + Pe log2(|X|-1).
	want := BinaryEntropy(0.25) + 0.25*infotheoryLog2(3)
	if got := FanoInequalityBound(0.25, 4); !approx(got, want, 1e-9) {
		t.Errorf("Fano bound = %v, want %v", got, want)
	}
}

func TestRateHelpers(t *testing.T) {
	if got := CodeRate(4, 7); !approx(got, 4.0/7.0, tol) {
		t.Errorf("CodeRate(4,7) = %v, want 4/7", got)
	}
	if got := CodeRate(1, 0); got != 0 {
		t.Errorf("CodeRate(1,0) = %v, want 0", got)
	}
	if got := EntropyRate([]float64{0.5, 0.5}, 2); !approx(got, 0.5, tol) {
		t.Errorf("EntropyRate = %v, want 0.5", got)
	}
	if got := CodingEfficiency(2.2, 2.2); !approx(got, 1, tol) {
		t.Errorf("CodingEfficiency optimal = %v, want 1", got)
	}
	if got := CodingEfficiency(1, 0); got != 0 {
		t.Errorf("CodingEfficiency(_,0) = %v, want 0", got)
	}
}
