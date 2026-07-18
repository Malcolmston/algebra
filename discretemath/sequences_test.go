package discretemath

import "testing"

func TestRunLengthSlice(t *testing.T) {
	in := []int{1, 1, 1, 2, 3, 3}
	runs := RunLengthEncode(in)
	want := []Run[int]{{1, 3}, {2, 1}, {3, 2}}
	if len(runs) != len(want) {
		t.Fatalf("runs = %v, want %v", runs, want)
	}
	for i := range want {
		if runs[i] != want[i] {
			t.Errorf("run %d = %v, want %v", i, runs[i], want[i])
		}
	}
	back := RunLengthDecode(runs)
	if len(back) != len(in) {
		t.Fatalf("decode length %d, want %d", len(back), len(in))
	}
	for i := range in {
		if back[i] != in[i] {
			t.Errorf("decode[%d] = %d, want %d", i, back[i], in[i])
		}
	}
	if len(RunLengthEncode([]int{})) != 0 {
		t.Error("empty input should encode to empty")
	}
}

func TestRunLengthString(t *testing.T) {
	cases := []struct {
		in, enc string
	}{
		{"aaabbc", "3a2b1c"},
		{"", ""},
		{"x", "1x"},
		{"aabbaa", "2a2b2a"},
	}
	for _, c := range cases {
		if got := RunLengthEncodeString(c.in); got != c.enc {
			t.Errorf("RunLengthEncodeString(%q) = %q, want %q", c.in, got, c.enc)
		}
		back, err := RunLengthDecodeString(c.enc)
		if err != nil {
			t.Fatalf("decode %q: %v", c.enc, err)
		}
		if back != c.in {
			t.Errorf("RunLengthDecodeString(%q) = %q, want %q", c.enc, back, c.in)
		}
	}
	// Multi-digit counts and unicode.
	if got := RunLengthEncodeString("aaaaaaaaaaaa"); got != "12a" {
		t.Errorf("multi-digit encode = %q, want 12a", got)
	}
	if _, err := RunLengthDecodeString("a3"); err == nil {
		t.Error("expected error: no leading digit")
	}
	if _, err := RunLengthDecodeString("3"); err == nil {
		t.Error("expected error: count with no char")
	}
}

func TestHammingDistanceString(t *testing.T) {
	got, err := HammingDistanceString("karolin", "kathrin")
	if err != nil {
		t.Fatal(err)
	}
	if got != 3 { // classic textbook example
		t.Errorf("HammingDistanceString = %d, want 3", got)
	}
	if _, err := HammingDistanceString("abc", "ab"); err == nil {
		t.Error("expected length-mismatch error")
	}
}

func TestLevenshtein(t *testing.T) {
	cases := []struct {
		a, b string
		want int
	}{
		{"", "", 0},
		{"abc", "", 3},
		{"", "abc", 3},
		{"kitten", "sitting", 3},
		{"flaw", "lawn", 2},
		{"gumbo", "gambol", 2},
		{"book", "back", 2},
		{"abc", "abc", 0},
		{"sunday", "saturday", 3},
	}
	for _, c := range cases {
		if got := LevenshteinDistance(c.a, c.b); got != c.want {
			t.Errorf("LevenshteinDistance(%q,%q) = %d, want %d", c.a, c.b, got, c.want)
		}
	}
}

func TestOptimalStringAlignment(t *testing.T) {
	cases := []struct {
		a, b string
		want int
	}{
		{"ca", "abc", 3},    // OSA distance (true Damerau is 2)
		{"ab", "ba", 1},     // single transposition
		{"abcd", "acbd", 1}, // adjacent transposition
		{"kitten", "sitting", 3},
		{"", "abc", 3},
	}
	for _, c := range cases {
		if got := OptimalStringAlignmentDistance(c.a, c.b); got != c.want {
			t.Errorf("OptimalStringAlignmentDistance(%q,%q) = %d, want %d", c.a, c.b, got, c.want)
		}
	}
}

func TestLongestCommonSubsequence(t *testing.T) {
	cases := []struct {
		a, b   string
		length int
		lcs    string
	}{
		{"ABCBDAB", "BDCAB", 4, "BCAB"},
		{"AGGTAB", "GXTXAYB", 4, "GTAB"},
		{"", "abc", 0, ""},
		{"abc", "abc", 3, "abc"},
		{"abc", "xyz", 0, ""},
	}
	for _, c := range cases {
		if got := LongestCommonSubsequenceLength(c.a, c.b); got != c.length {
			t.Errorf("LCSLength(%q,%q) = %d, want %d", c.a, c.b, got, c.length)
		}
		got := LongestCommonSubsequence(c.a, c.b)
		if len(got) != c.length {
			t.Errorf("LCS(%q,%q) = %q (len %d), want length %d", c.a, c.b, got, len(got), c.length)
		}
		if !discretemathIsSubsequence(got, c.a) || !discretemathIsSubsequence(got, c.b) {
			t.Errorf("LCS %q is not a common subsequence of %q and %q", got, c.a, c.b)
		}
	}
}

// discretemathIsSubsequence reports whether sub is a subsequence of s, used to
// validate LCS results without pinning a specific tie-break.
func discretemathIsSubsequence(sub, s string) bool {
	rs, rb := []rune(sub), []rune(s)
	i := 0
	for j := 0; i < len(rs) && j < len(rb); j++ {
		if rs[i] == rb[j] {
			i++
		}
	}
	return i == len(rs)
}
