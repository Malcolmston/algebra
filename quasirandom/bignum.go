package quasirandom

import (
	"math/big"
)

// RadicalInverseBig returns the radical inverse in the given base of a
// non-negative index supplied as a big.Int, as an exact big.Rat in [0,1). It
// lets the van der Corput construction reach indices far beyond the range of
// uint64. It returns an error when base < 2 or n is negative.
func RadicalInverseBig(base int, n *big.Int) (*big.Rat, error) {
	if base < 2 {
		return nil, ErrBadBase
	}
	if n.Sign() < 0 {
		return nil, ErrNonPositive
	}
	b := big.NewInt(int64(base))
	result := new(big.Rat)
	// denom accumulates base^(k+1).
	denom := new(big.Int).Set(b)
	q := new(big.Int).Set(n)
	r := new(big.Int)
	for q.Sign() > 0 {
		q.QuoRem(q, b, r)
		if r.Sign() != 0 {
			term := new(big.Rat).SetFrac(new(big.Int).Set(r), new(big.Int).Set(denom))
			result.Add(result, term)
		}
		denom.Mul(denom, b)
	}
	return result, nil
}

// RadicalInverseBigFloat returns RadicalInverseBig evaluated to a float64,
// convenient when only ordinary precision is needed but the index exceeds
// uint64. It returns an error when base < 2 or n is negative.
func RadicalInverseBigFloat(base int, n *big.Int) (float64, error) {
	r, err := RadicalInverseBig(base, n)
	if err != nil {
		return 0, err
	}
	f, _ := r.Float64()
	return f, nil
}

// HaltonBig returns the n-th Halton point for an index given as a big.Int, using
// the first dim primes as bases and returning each coordinate as a float64. It
// returns an error when dim < 1 or n is negative.
func HaltonBig(dim int, n *big.Int) ([]float64, error) {
	if dim < 1 {
		return nil, ErrDimension
	}
	bases, err := PrimeBases(dim)
	if err != nil {
		return nil, err
	}
	out := make([]float64, dim)
	for i, b := range bases {
		v, err := RadicalInverseBigFloat(b, n)
		if err != nil {
			return nil, err
		}
		out[i] = v
	}
	return out, nil
}

// BigFactorial returns n! as a big.Int, a small helper occasionally useful when
// working with the exact combinatorics behind Faure generator matrices. It
// returns an error when n is negative.
func BigFactorial(n int) (*big.Int, error) {
	if n < 0 {
		return nil, ErrNonPositive
	}
	r := big.NewInt(1)
	for i := int64(2); i <= int64(n); i++ {
		r.Mul(r, big.NewInt(i))
	}
	return r, nil
}

// BigBinomial returns the binomial coefficient C(n,k) as an exact big.Int. It
// returns an error when n or k is negative.
func BigBinomial(n, k int) (*big.Int, error) {
	if n < 0 || k < 0 {
		return nil, ErrNonPositive
	}
	if k > n {
		return big.NewInt(0), nil
	}
	return new(big.Int).Binomial(int64(n), int64(k)), nil
}
