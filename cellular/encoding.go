package cellular

import (
	"errors"
	"math/big"
)

// EncodeBaseK interprets vals as the digits of a base-k number, most-significant
// digit first, and returns its integer value. Digits outside 0..k-1 yield an
// error.
func EncodeBaseK(vals []int, k int) (int, error) {
	if k < 2 {
		return 0, errors.New("cellular: EncodeBaseK needs k >= 2")
	}
	n := 0
	for i, v := range vals {
		if v < 0 || v >= k {
			return 0, errors.New("cellular: EncodeBaseK digit out of range")
		}
		n = n*k + v
		_ = i
	}
	return n, nil
}

// DecodeBaseK returns the digits digits-long base-k representation of code, most
// significant digit first. It returns an error if code cannot be represented in
// the requested number of digits.
func DecodeBaseK(code, digits, k int) ([]int, error) {
	if k < 2 {
		return nil, errors.New("cellular: DecodeBaseK needs k >= 2")
	}
	if digits < 0 {
		return nil, errors.New("cellular: DecodeBaseK needs digits >= 0")
	}
	if code < 0 {
		return nil, errors.New("cellular: DecodeBaseK needs code >= 0")
	}
	out := make([]int, digits)
	c := code
	for i := digits - 1; i >= 0; i-- {
		out[i] = c % k
		c /= k
	}
	if c != 0 {
		return nil, errors.New("cellular: DecodeBaseK code does not fit in digits")
	}
	return out, nil
}

// PopCountK returns the number of non-zero digits among vals.
func PopCountK(vals []int) int {
	c := 0
	for _, v := range vals {
		if v != 0 {
			c++
		}
	}
	return c
}

// NeighborhoodCount returns the number of distinct neighbourhoods for a k-state,
// radius-r automaton, namely k^(2r+1).
func NeighborhoodCount(k, r int) *big.Int {
	base := big.NewInt(int64(k))
	return new(big.Int).Exp(base, big.NewInt(int64(2*r+1)), nil)
}

// NumRules returns the number of distinct k-state, radius-r rules as an exact
// big integer, k^(k^(2r+1)).
func NumRules(k, r int) *big.Int {
	exp := NeighborhoodCount(k, r)
	return new(big.Int).Exp(big.NewInt(int64(k)), exp, nil)
}

// NumTotalisticRules returns the number of distinct k-state, radius-r totalistic
// rules, k^((2r+1)(k-1)+1).
func NumTotalisticRules(k, r int) *big.Int {
	entries := (2*r+1)*(k-1) + 1
	return new(big.Int).Exp(big.NewInt(int64(k)), big.NewInt(int64(entries)), nil)
}

// NumOuterTotalisticRules returns the number of distinct k-state, radius-r
// outer-totalistic rules, k^(k*(2r(k-1)+1)).
func NumOuterTotalisticRules(k, r int) *big.Int {
	entries := k * (2*r*(k-1) + 1)
	return new(big.Int).Exp(big.NewInt(int64(k)), big.NewInt(int64(entries)), nil)
}
