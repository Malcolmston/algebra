package discretemath

import (
	"strconv"
	"testing"
)

func TestBaseConversionRoundTrip(t *testing.T) {
	values := []int64{0, 1, -1, 7, -7, 255, 1000, -1000, 1 << 40, -(1 << 40),
		9223372036854775807, -9223372036854775808}
	for _, base := range []int{2, 8, 10, 16, 36} {
		for _, v := range values {
			s, err := ToBase(v, base)
			if err != nil {
				t.Fatalf("ToBase(%d,%d): %v", v, base, err)
			}
			// Cross-check magnitude against strconv.
			if want := strconv.FormatInt(v, base); s != want {
				t.Errorf("ToBase(%d,%d) = %q, want %q", v, base, s, want)
			}
			back, err := FromBase(s, base)
			if err != nil {
				t.Fatalf("FromBase(%q,%d): %v", s, base, err)
			}
			if back != v {
				t.Errorf("round trip base %d: %d -> %q -> %d", base, v, s, back)
			}
		}
	}
}

func TestBaseErrors(t *testing.T) {
	if _, err := ToBase(10, 1); err == nil {
		t.Error("expected base error")
	}
	if _, err := ToBase(10, 37); err == nil {
		t.Error("expected base error")
	}
	if _, err := FromBase("12", 37); err == nil {
		t.Error("expected base error")
	}
	if _, err := FromBase("1g", 16); err == nil {
		t.Error("expected invalid-digit error")
	}
	if _, err := FromBase("", 10); err == nil {
		t.Error("expected empty-string error")
	}
	if _, err := FromBaseUint("18446744073709551616", 10); err == nil {
		t.Error("expected overflow error")
	}
}

func TestBaseUintKnown(t *testing.T) {
	cases := []struct {
		n    uint64
		base int
		want string
	}{
		{255, 16, "ff"},
		{255, 2, "11111111"},
		{35, 36, "z"},
		{36, 36, "10"},
		{1000, 10, "1000"},
	}
	for _, c := range cases {
		got, err := ToBaseUint(c.n, c.base)
		if err != nil {
			t.Fatal(err)
		}
		if got != c.want {
			t.Errorf("ToBaseUint(%d,%d) = %q, want %q", c.n, c.base, got, c.want)
		}
	}
}

func TestRoman(t *testing.T) {
	cases := []struct {
		n int
		r string
	}{
		{1, "I"}, {4, "IV"}, {9, "IX"}, {14, "XIV"}, {40, "XL"},
		{90, "XC"}, {400, "CD"}, {900, "CM"}, {1994, "MCMXCIV"},
		{2023, "MMXXIII"}, {3888, "MMMDCCCLXXXVIII"}, {3999, "MMMCMXCIX"},
		{58, "LVIII"}, {49, "XLIX"},
	}
	for _, c := range cases {
		got, err := IntToRoman(c.n)
		if err != nil {
			t.Fatalf("IntToRoman(%d): %v", c.n, err)
		}
		if got != c.r {
			t.Errorf("IntToRoman(%d) = %q, want %q", c.n, got, c.r)
		}
		back, err := RomanToInt(c.r)
		if err != nil {
			t.Fatalf("RomanToInt(%q): %v", c.r, err)
		}
		if back != c.n {
			t.Errorf("RomanToInt(%q) = %d, want %d", c.r, back, c.n)
		}
	}
	// Full round trip across the whole valid range.
	for n := 1; n <= 3999; n++ {
		r, _ := IntToRoman(n)
		if back, err := RomanToInt(r); err != nil || back != n {
			t.Fatalf("round trip failed at %d (%q): back=%d err=%v", n, r, back, err)
		}
	}
}

func TestRomanErrors(t *testing.T) {
	if _, err := IntToRoman(0); err == nil {
		t.Error("expected range error for 0")
	}
	if _, err := IntToRoman(4000); err == nil {
		t.Error("expected range error for 4000")
	}
	for _, bad := range []string{"", "IIII", "VV", "IC", "ABC", "MMMM", "XM"} {
		if _, err := RomanToInt(bad); err == nil {
			t.Errorf("expected error for non-canonical %q", bad)
		}
		if IsValidRoman(bad) {
			t.Errorf("IsValidRoman(%q) should be false", bad)
		}
	}
	if !IsValidRoman("mcmxciv") { // case-insensitive
		t.Error("lowercase canonical numeral should be valid")
	}
}

func TestIntToWords(t *testing.T) {
	cases := []struct {
		n int64
		w string
	}{
		{0, "zero"},
		{1, "one"},
		{7, "seven"},
		{13, "thirteen"},
		{20, "twenty"},
		{21, "twenty-one"},
		{100, "one hundred"},
		{101, "one hundred one"},
		{342, "three hundred forty-two"},
		{1000, "one thousand"},
		{1234, "one thousand two hundred thirty-four"},
		{1000000, "one million"},
		{-5, "negative five"},
		{2000000, "two million"},
		{1000000000000, "one trillion"},
	}
	for _, c := range cases {
		if got := IntToWords(c.n); got != c.w {
			t.Errorf("IntToWords(%d) = %q, want %q", c.n, got, c.w)
		}
	}
}

func TestIntToWordsOrdinal(t *testing.T) {
	cases := []struct {
		n int64
		w string
	}{
		{1, "first"}, {2, "second"}, {3, "third"}, {4, "fourth"},
		{5, "fifth"}, {8, "eighth"}, {9, "ninth"}, {11, "eleventh"},
		{12, "twelfth"}, {13, "thirteenth"}, {20, "twentieth"},
		{21, "twenty-first"}, {30, "thirtieth"}, {100, "one hundredth"},
		{101, "one hundred first"}, {1000, "one thousandth"},
		{1000000, "one millionth"},
	}
	for _, c := range cases {
		if got := IntToWordsOrdinal(c.n); got != c.w {
			t.Errorf("IntToWordsOrdinal(%d) = %q, want %q", c.n, got, c.w)
		}
	}
}
