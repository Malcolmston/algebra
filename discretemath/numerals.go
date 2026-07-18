package discretemath

import "strings"

// discretemathRomanValues and discretemathRomanSymbols encode the subtractive
// Roman-numeral pairs in descending order of value, which lets both encoding and
// greedy decoding share a single table.
var (
	discretemathRomanValues  = []int{1000, 900, 500, 400, 100, 90, 50, 40, 10, 9, 5, 4, 1}
	discretemathRomanSymbols = []string{"M", "CM", "D", "CD", "C", "XC", "L", "XL", "X", "IX", "V", "IV", "I"}
)

// IntToRoman converts an integer in the range 1..3999 into its canonical Roman
// numeral using subtractive notation (for example 1994 becomes "MCMXCIV"). It
// returns an error for values outside that range.
func IntToRoman(n int) (string, error) {
	if n < 1 || n > 3999 {
		return "", discretemathErrorf("IntToRoman: value %d out of range [1,3999]", n)
	}
	var sb strings.Builder
	for i, v := range discretemathRomanValues {
		for n >= v {
			sb.WriteString(discretemathRomanSymbols[i])
			n -= v
		}
	}
	return sb.String(), nil
}

// RomanToInt converts a Roman numeral into its integer value. Input is treated
// case-insensitively. The numeral must be in canonical (standard subtractive)
// form; strings such as "IIII" or "VV" are rejected with an error, as are
// strings containing non-Roman characters.
func RomanToInt(s string) (int, error) {
	if s == "" {
		return 0, discretemathErrorf("RomanToInt: empty string")
	}
	upper := strings.ToUpper(s)
	n := 0
	i := 0
	// Greedy longest-prefix match against the descending table.
	for i < len(upper) {
		matched := false
		for k, sym := range discretemathRomanSymbols {
			if strings.HasPrefix(upper[i:], sym) {
				n += discretemathRomanValues[k]
				i += len(sym)
				matched = true
				break
			}
		}
		if !matched {
			return 0, discretemathErrorf("RomanToInt: invalid Roman numeral %q", s)
		}
	}
	// Validate canonical form by re-encoding.
	if n < 1 || n > 3999 {
		return 0, discretemathErrorf("RomanToInt: value %d out of range [1,3999]", n)
	}
	canonical, _ := IntToRoman(n)
	if canonical != upper {
		return 0, discretemathErrorf("RomanToInt: %q is not canonical (expected %q)", s, canonical)
	}
	return n, nil
}

// IsValidRoman reports whether s is a well-formed canonical Roman numeral in the
// range 1..3999. Input is treated case-insensitively.
func IsValidRoman(s string) bool {
	_, err := RomanToInt(s)
	return err == nil
}

// discretemathOnes maps 0..19 to their English words; index 0 ("zero") is used
// only for the standalone-zero case.
var discretemathOnes = []string{
	"zero", "one", "two", "three", "four", "five", "six", "seven", "eight",
	"nine", "ten", "eleven", "twelve", "thirteen", "fourteen", "fifteen",
	"sixteen", "seventeen", "eighteen", "nineteen",
}

// discretemathTens maps a tens digit 2..9 to its English word; indices 0 and 1
// are unused.
var discretemathTens = []string{
	"", "", "twenty", "thirty", "forty", "fifty", "sixty", "seventy",
	"eighty", "ninety",
}

// discretemathScales names the successive thousands groups; index i corresponds
// to 1000**i.
var discretemathScales = []string{
	"", " thousand", " million", " billion", " trillion", " quadrillion",
	" quintillion",
}

// discretemathThreeDigits renders an integer in 1..999 into English words with
// hyphenated compound tens (for example 342 becomes "three hundred forty-two").
func discretemathThreeDigits(n int) string {
	var parts []string
	if n >= 100 {
		parts = append(parts, discretemathOnes[n/100], "hundred")
		n %= 100
	}
	if n >= 20 {
		word := discretemathTens[n/10]
		if n%10 > 0 {
			word += "-" + discretemathOnes[n%10]
		}
		parts = append(parts, word)
	} else if n > 0 {
		parts = append(parts, discretemathOnes[n])
	}
	return strings.Join(parts, " ")
}

// IntToWords converts a signed 64-bit integer into its English cardinal words,
// covering the full int64 range up to quintillions. Zero renders as "zero" and
// negative values are prefixed with "negative".
func IntToWords(n int64) string {
	if n == 0 {
		return "zero"
	}
	neg := n < 0
	var u uint64
	if neg {
		u = uint64(-(n + 1)) + 1
	} else {
		u = uint64(n)
	}
	var groups []int
	for u > 0 {
		groups = append(groups, int(u%1000))
		u /= 1000
	}
	var segs []string
	for i := len(groups) - 1; i >= 0; i-- {
		if groups[i] == 0 {
			continue
		}
		segs = append(segs, discretemathThreeDigits(groups[i])+discretemathScales[i])
	}
	res := strings.Join(segs, " ")
	if neg {
		res = "negative " + res
	}
	return res
}

// discretemathOrdinalSpecials maps cardinal words whose ordinal form is
// irregular to that ordinal form.
var discretemathOrdinalSpecials = map[string]string{
	"one": "first", "two": "second", "three": "third", "five": "fifth",
	"eight": "eighth", "nine": "ninth", "twelve": "twelfth",
	"twenty": "twentieth", "thirty": "thirtieth", "forty": "fortieth",
	"fifty": "fiftieth", "sixty": "sixtieth", "seventy": "seventieth",
	"eighty": "eightieth", "ninety": "ninetieth",
}

// discretemathOrdinalize converts a single cardinal word into its ordinal form,
// applying the irregular table and otherwise appending "th".
func discretemathOrdinalize(word string) string {
	if v, ok := discretemathOrdinalSpecials[word]; ok {
		return v
	}
	return word + "th"
}

// IntToWordsOrdinal converts a signed 64-bit integer into its English ordinal
// words (for example 21 becomes "twenty-first" and 100 becomes "one hundredth").
// Negative values are prefixed with "negative".
func IntToWordsOrdinal(n int64) string {
	cardinal := IntToWords(n)
	head, tail := "", cardinal
	if idx := strings.LastIndex(cardinal, " "); idx >= 0 {
		head, tail = cardinal[:idx+1], cardinal[idx+1:]
	}
	prefix := ""
	if idx := strings.LastIndex(tail, "-"); idx >= 0 {
		prefix, tail = tail[:idx+1], tail[idx+1:]
	}
	return head + prefix + discretemathOrdinalize(tail)
}
