package discretemath

import "strings"

// DeBruijnSequence returns a De Bruijn sequence B(k, n): a sequence of symbols
// drawn from the alphabet {0, 1, ..., k-1} of length k**n in which every possible
// length-n string over the alphabet occurs exactly once as a (cyclic)
// contiguous substring.
//
// It uses the standard prefer-smallest generation via Lyndon words and returns
// an error when k or n is less than one.
func DeBruijnSequence(k, n int) ([]int, error) {
	if k < 1 {
		return nil, discretemathErrorf("DeBruijnSequence: alphabet size k must be >= 1, got %d", k)
	}
	if n < 1 {
		return nil, discretemathErrorf("DeBruijnSequence: word length n must be >= 1, got %d", n)
	}
	a := make([]int, n+1)
	seq := make([]int, 0, discretemathPow(k, n))
	var db func(t, p int)
	db = func(t, p int) {
		if t > n {
			if n%p == 0 {
				for j := 1; j <= p; j++ {
					seq = append(seq, a[j])
				}
			}
			return
		}
		a[t] = a[t-p]
		db(t+1, p)
		for j := a[t-p] + 1; j < k; j++ {
			a[t] = j
			db(t+1, t)
		}
	}
	db(1, 1)
	return seq, nil
}

// DeBruijnString returns a De Bruijn sequence of order n over the given alphabet
// as a string. The alphabet is treated as an ordered list of runes; every
// length-n string over those runes appears exactly once as a cyclic substring.
//
// It returns an error when the alphabet is empty or n is less than one.
func DeBruijnString(alphabet string, n int) (string, error) {
	runes := []rune(alphabet)
	k := len(runes)
	if k < 1 {
		return "", discretemathErrorf("DeBruijnString: alphabet must be non-empty")
	}
	seq, err := DeBruijnSequence(k, n)
	if err != nil {
		return "", err
	}
	var sb strings.Builder
	sb.Grow(len(seq))
	for _, idx := range seq {
		sb.WriteRune(runes[idx])
	}
	return sb.String(), nil
}

// discretemathPow returns base**exp for small non-negative exponents, used to
// size buffers. It assumes the result fits in an int.
func discretemathPow(base, exp int) int {
	result := 1
	for i := 0; i < exp; i++ {
		result *= base
	}
	return result
}
