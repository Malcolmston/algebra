package discretemath

import "testing"

// discretemathHasAllWindows checks that every length-n window over the alphabet
// {0..k-1} appears exactly once as a cyclic substring of seq.
func discretemathHasAllWindows(seq []int, k, n int) bool {
	L := len(seq)
	if L != discretemathPow(k, n) {
		return false
	}
	seen := make(map[int]bool, L)
	for i := 0; i < L; i++ {
		code := 0
		for j := 0; j < n; j++ {
			code = code*k + seq[(i+j)%L]
		}
		if seen[code] {
			return false
		}
		seen[code] = true
	}
	return len(seen) == L
}

func TestDeBruijnSequence(t *testing.T) {
	// Known B(2,3): prefer-smallest generation yields 0,0,0,1,0,1,1,1.
	seq, err := DeBruijnSequence(2, 3)
	if err != nil {
		t.Fatal(err)
	}
	want := []int{0, 0, 0, 1, 0, 1, 1, 1}
	if len(seq) != len(want) {
		t.Fatalf("B(2,3) = %v, want %v", seq, want)
	}
	for i := range want {
		if seq[i] != want[i] {
			t.Fatalf("B(2,3) = %v, want %v", seq, want)
		}
	}

	// Structural check across several parameters.
	for _, p := range []struct{ k, n int }{
		{2, 1}, {2, 3}, {2, 5}, {3, 2}, {3, 3}, {4, 2}, {5, 2}, {10, 2},
	} {
		s, err := DeBruijnSequence(p.k, p.n)
		if err != nil {
			t.Fatalf("DeBruijnSequence(%d,%d): %v", p.k, p.n, err)
		}
		if !discretemathHasAllWindows(s, p.k, p.n) {
			t.Errorf("B(%d,%d) is not a valid De Bruijn sequence", p.k, p.n)
		}
	}
}

func TestDeBruijnErrors(t *testing.T) {
	if _, err := DeBruijnSequence(0, 3); err == nil {
		t.Error("expected error for k=0")
	}
	if _, err := DeBruijnSequence(2, 0); err == nil {
		t.Error("expected error for n=0")
	}
	if _, err := DeBruijnString("", 2); err == nil {
		t.Error("expected error for empty alphabet")
	}
}

func TestDeBruijnString(t *testing.T) {
	s, err := DeBruijnString("01", 3)
	if err != nil {
		t.Fatal(err)
	}
	if s != "00010111" {
		t.Errorf("DeBruijnString(01,3) = %q, want %q", s, "00010111")
	}
	// Length must be k^n.
	s2, err := DeBruijnString("abc", 2)
	if err != nil {
		t.Fatal(err)
	}
	if len(s2) != 9 {
		t.Errorf("len = %d, want 9", len(s2))
	}
}
