package seq

import (
	"reflect"
	"testing"
)

func TestCollatz(t *testing.T) {
	// The canonical trajectory of 6.
	want := []uint64{6, 3, 10, 5, 16, 8, 4, 2, 1}
	if got := CollatzSequence(6); !reflect.DeepEqual(got, want) {
		t.Errorf("CollatzSequence(6) = %v, want %v", got, want)
	}
	if got := CollatzSteps(6); got != 8 {
		t.Errorf("CollatzSteps(6) = %d, want 8", got)
	}
	if got := CollatzSteps(27); got != 111 {
		t.Errorf("CollatzSteps(27) = %d, want 111", got)
	}
	if got := CollatzMax(27); got != 9232 {
		t.Errorf("CollatzMax(27) = %d, want 9232", got)
	}
	if CollatzSteps(1) != 0 {
		t.Errorf("CollatzSteps(1) should be 0")
	}
}

func TestHappy(t *testing.T) {
	want := []uint64{1, 7, 10, 13, 19, 23, 28, 31, 32, 44}
	if got := HappyNumbers(10); !reflect.DeepEqual(got, want) {
		t.Errorf("HappyNumbers(10) = %v, want %v", got, want)
	}
	if IsHappy(2) {
		t.Errorf("2 is not happy")
	}
	if !IsHappy(7) {
		t.Errorf("7 is happy")
	}
}

func TestKaprekarNumbers(t *testing.T) {
	want := []uint64{1, 9, 45, 55, 99, 297, 703, 999, 2223, 2728, 4879, 4950, 5050, 5292, 7272, 7777, 9999}
	if got := KaprekarNumbers(9999); !reflect.DeepEqual(got, want) {
		t.Errorf("KaprekarNumbers(9999) = %v, want %v", got, want)
	}
	// Spot checks.
	if !IsKaprekar(45) || !IsKaprekar(297) {
		t.Errorf("45 and 297 are Kaprekar")
	}
	if IsKaprekar(100) || IsKaprekar(10) {
		t.Errorf("10 and 100 are not Kaprekar")
	}
}

func TestKaprekarRoutine(t *testing.T) {
	if got := KaprekarStep(3524); got != 3087 {
		t.Errorf("KaprekarStep(3524) = %d, want 3087", got)
	}
	seq := KaprekarSequence(3524)
	if seq[len(seq)-1] != 6174 {
		t.Errorf("KaprekarSequence(3524) should end at 6174, got %v", seq)
	}
	// Every non-repdigit four-digit number converges to 6174.
	for n := 1000; n <= 9999; n++ {
		d := []int{n / 1000, n / 100 % 10, n / 10 % 10, n % 10}
		if d[0] == d[1] && d[1] == d[2] && d[2] == d[3] {
			continue // repdigit converges to 0
		}
		s := KaprekarSequence(n)
		if s[len(s)-1] != 6174 {
			t.Errorf("KaprekarSequence(%d) did not reach 6174: end %d", n, s[len(s)-1])
		}
	}
}

func TestRecaman(t *testing.T) {
	want := []int64{0, 1, 3, 6, 2, 7, 13, 20, 12, 21, 11, 22, 10, 23, 9, 24, 8, 25, 43, 62}
	if got := RecamanSequence(20); !reflect.DeepEqual(got, want) {
		t.Errorf("RecamanSequence(20) = %v, want %v", got, want)
	}
}

func TestLookAndSay(t *testing.T) {
	want := []string{"1", "11", "21", "1211", "111221", "312211", "13112221"}
	if got := LookAndSaySequence("1", 7); !reflect.DeepEqual(got, want) {
		t.Errorf("LookAndSaySequence(\"1\",7) = %v, want %v", got, want)
	}
	if got := LookAndSayStep("111221"); got != "312211" {
		t.Errorf("LookAndSayStep(111221) = %s, want 312211", got)
	}
}

func BenchmarkLookAndSay(b *testing.B) {
	for i := 0; i < b.N; i++ {
		_ = LookAndSaySequence("1", 30)
	}
}
