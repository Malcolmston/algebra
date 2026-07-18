package discretemath

import (
	"strconv"
	"strings"
	"unicode/utf8"
)

// Run is a single maximal run of an equal value produced by run-length encoding:
// Value repeated Count times.
type Run[T comparable] struct {
	// Value is the element that repeats.
	Value T
	// Count is the number of consecutive repetitions; always at least one.
	Count int
}

// RunLengthEncode compresses a slice into a sequence of runs, each recording a
// value and the number of times it repeats consecutively. Encoding an empty
// slice yields an empty result.
func RunLengthEncode[T comparable](s []T) []Run[T] {
	var out []Run[T]
	for i := 0; i < len(s); {
		j := i + 1
		for j < len(s) && s[j] == s[i] {
			j++
		}
		out = append(out, Run[T]{Value: s[i], Count: j - i})
		i = j
	}
	return out
}

// RunLengthDecode expands a sequence of runs back into the original slice. Runs
// with a non-positive count contribute no elements.
func RunLengthDecode[T comparable](runs []Run[T]) []T {
	total := 0
	for _, r := range runs {
		if r.Count > 0 {
			total += r.Count
		}
	}
	out := make([]T, 0, total)
	for _, r := range runs {
		for k := 0; k < r.Count; k++ {
			out = append(out, r.Value)
		}
	}
	return out
}

// RunLengthEncodeString compresses a string into the form "<count><rune>" for
// each maximal run (for example "aaabbc" becomes "3a2b1c"). The string is
// processed by rune, so multi-byte UTF-8 characters are handled correctly.
func RunLengthEncodeString(s string) string {
	var sb strings.Builder
	runes := []rune(s)
	for i := 0; i < len(runes); {
		j := i + 1
		for j < len(runes) && runes[j] == runes[i] {
			j++
		}
		sb.WriteString(strconv.Itoa(j - i))
		sb.WriteRune(runes[i])
		i = j
	}
	return sb.String()
}

// RunLengthDecodeString expands a string produced by RunLengthEncodeString back
// into its original form. It returns an error if the input is not a sequence of
// decimal-count / single-rune pairs.
func RunLengthDecodeString(s string) (string, error) {
	var sb strings.Builder
	i := 0
	for i < len(s) {
		start := i
		for i < len(s) && s[i] >= '0' && s[i] <= '9' {
			i++
		}
		if i == start {
			return "", discretemathErrorf("RunLengthDecodeString: expected digit at position %d", i)
		}
		count, err := strconv.Atoi(s[start:i])
		if err != nil {
			return "", discretemathErrorf("RunLengthDecodeString: bad count %q: %v", s[start:i], err)
		}
		if i >= len(s) {
			return "", discretemathErrorf("RunLengthDecodeString: count %d not followed by a character", count)
		}
		r, size := utf8.DecodeRuneInString(s[i:])
		i += size
		for k := 0; k < count; k++ {
			sb.WriteRune(r)
		}
	}
	return sb.String(), nil
}

// HammingDistanceString returns the number of positions at which two equal
// length strings differ, comparing by rune. It returns an error when the strings
// contain a different number of runes.
func HammingDistanceString(a, b string) (int, error) {
	ra, rb := []rune(a), []rune(b)
	if len(ra) != len(rb) {
		return 0, discretemathErrorf("HammingDistanceString: rune length mismatch %d != %d", len(ra), len(rb))
	}
	d := 0
	for i := range ra {
		if ra[i] != rb[i] {
			d++
		}
	}
	return d, nil
}

// LevenshteinDistance returns the Levenshtein edit distance between two strings:
// the minimum number of single-character insertions, deletions and substitutions
// needed to transform a into b. Comparison is by rune.
func LevenshteinDistance(a, b string) int {
	ra, rb := []rune(a), []rune(b)
	if len(ra) == 0 {
		return len(rb)
	}
	if len(rb) == 0 {
		return len(ra)
	}
	prev := make([]int, len(rb)+1)
	curr := make([]int, len(rb)+1)
	for j := range prev {
		prev[j] = j
	}
	for i := 1; i <= len(ra); i++ {
		curr[0] = i
		for j := 1; j <= len(rb); j++ {
			cost := 1
			if ra[i-1] == rb[j-1] {
				cost = 0
			}
			curr[j] = discretemathMin3(
				prev[j]+1,      // deletion
				curr[j-1]+1,    // insertion
				prev[j-1]+cost, // substitution
			)
		}
		prev, curr = curr, prev
	}
	return prev[len(rb)]
}

// OptimalStringAlignmentDistance returns the optimal string alignment distance
// between two strings, the restricted form of the Damerau-Levenshtein distance
// that additionally allows the transposition of two adjacent characters as a
// single edit, subject to the constraint that no substring is edited more than
// once. Comparison is by rune.
func OptimalStringAlignmentDistance(a, b string) int {
	ra, rb := []rune(a), []rune(b)
	n, m := len(ra), len(rb)
	if n == 0 {
		return m
	}
	if m == 0 {
		return n
	}
	d := make([][]int, n+1)
	for i := range d {
		d[i] = make([]int, m+1)
		d[i][0] = i
	}
	for j := 0; j <= m; j++ {
		d[0][j] = j
	}
	for i := 1; i <= n; i++ {
		for j := 1; j <= m; j++ {
			cost := 1
			if ra[i-1] == rb[j-1] {
				cost = 0
			}
			d[i][j] = discretemathMin3(
				d[i-1][j]+1,
				d[i][j-1]+1,
				d[i-1][j-1]+cost,
			)
			if i > 1 && j > 1 && ra[i-1] == rb[j-2] && ra[i-2] == rb[j-1] {
				if t := d[i-2][j-2] + 1; t < d[i][j] {
					d[i][j] = t
				}
			}
		}
	}
	return d[n][m]
}

// LongestCommonSubsequenceLength returns the length of the longest common
// subsequence of two strings, comparing by rune. A subsequence keeps relative
// order but need not be contiguous.
func LongestCommonSubsequenceLength(a, b string) int {
	ra, rb := []rune(a), []rune(b)
	prev := make([]int, len(rb)+1)
	curr := make([]int, len(rb)+1)
	for i := 1; i <= len(ra); i++ {
		for j := 1; j <= len(rb); j++ {
			if ra[i-1] == rb[j-1] {
				curr[j] = prev[j-1] + 1
			} else if prev[j] >= curr[j-1] {
				curr[j] = prev[j]
			} else {
				curr[j] = curr[j-1]
			}
		}
		prev, curr = curr, prev
	}
	return prev[len(rb)]
}

// LongestCommonSubsequence returns one longest common subsequence of two strings
// as a string, comparing by rune. When several subsequences share the maximum
// length, the one preferring earlier positions in a is returned.
func LongestCommonSubsequence(a, b string) string {
	ra, rb := []rune(a), []rune(b)
	n, m := len(ra), len(rb)
	table := make([][]int, n+1)
	for i := range table {
		table[i] = make([]int, m+1)
	}
	for i := 1; i <= n; i++ {
		for j := 1; j <= m; j++ {
			if ra[i-1] == rb[j-1] {
				table[i][j] = table[i-1][j-1] + 1
			} else if table[i-1][j] >= table[i][j-1] {
				table[i][j] = table[i-1][j]
			} else {
				table[i][j] = table[i][j-1]
			}
		}
	}
	out := make([]rune, table[n][m])
	k := len(out)
	i, j := n, m
	for i > 0 && j > 0 {
		if ra[i-1] == rb[j-1] {
			k--
			out[k] = ra[i-1]
			i--
			j--
		} else if table[i-1][j] >= table[i][j-1] {
			i--
		} else {
			j--
		}
	}
	return string(out)
}

// discretemathMin3 returns the smallest of three integers.
func discretemathMin3(a, b, c int) int {
	m := a
	if b < m {
		m = b
	}
	if c < m {
		m = c
	}
	return m
}
