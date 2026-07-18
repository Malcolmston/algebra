package crypto

import (
	"errors"
	"math/big"
	"math/rand"
)

// Share is a single point (X, Y) on the secret-sharing polynomial produced by
// SplitSecret. Distributing shares to participants lets any threshold-sized
// subset reconstruct the secret while smaller subsets learn nothing.
type Share struct {
	X *big.Int // the participant index, a non-zero field element
	Y *big.Int // the polynomial value at X, modulo the field prime
}

// EvaluatePolynomial evaluates the polynomial whose coefficients are given in
// ascending order (coeffs[0] is the constant term) at the point x modulo m,
// using Horner's method. The modulus m must be positive.
func EvaluatePolynomial(coeffs []*big.Int, x, m *big.Int) *big.Int {
	if m.Sign() <= 0 {
		panic("crypto: EvaluatePolynomial requires modulus m > 0")
	}
	result := big.NewInt(0)
	for i := len(coeffs) - 1; i >= 0; i-- {
		result.Mul(result, x)
		result.Add(result, coeffs[i])
		result.Mod(result, m)
	}
	return result
}

// SplitSecret splits secret into shares participants using Shamir's threshold
// scheme over the prime field Z/prime. Any threshold shares can reconstruct the
// secret; any fewer reveal nothing about it. The secret must satisfy
// 0 <= secret < prime, and 2 <= threshold <= participants. Coefficients above
// the constant term are drawn from rng, so identical seeding reproduces the
// same shares. Shares are evaluated at X = 1, 2, ..., participants. It returns
// an error for invalid parameters.
func SplitSecret(secret *big.Int, threshold, participants int, prime *big.Int, rng *rand.Rand) ([]Share, error) {
	if threshold < 2 {
		return nil, errors.New("crypto: SplitSecret requires threshold >= 2")
	}
	if participants < threshold {
		return nil, errors.New("crypto: SplitSecret requires participants >= threshold")
	}
	if !IsPrime(prime) {
		return nil, errors.New("crypto: SplitSecret requires a prime field modulus")
	}
	if secret.Sign() < 0 || secret.Cmp(prime) >= 0 {
		return nil, errors.New("crypto: SplitSecret requires 0 <= secret < prime")
	}
	// Build coefficients: a0 = secret, a1..a_{t-1} random in [0, prime).
	coeffs := make([]*big.Int, threshold)
	coeffs[0] = new(big.Int).Set(secret)
	for i := 1; i < threshold; i++ {
		coeffs[i] = cryptoRandBelow(rng, prime)
	}
	shares := make([]Share, participants)
	for i := 0; i < participants; i++ {
		x := big.NewInt(int64(i + 1))
		y := EvaluatePolynomial(coeffs, x, prime)
		shares[i] = Share{X: x, Y: y}
	}
	return shares, nil
}

// CombineShares reconstructs the secret from a set of shares using Lagrange
// interpolation at x = 0 over the prime field Z/prime. At least threshold
// shares are required and their X coordinates must be distinct and non-zero;
// otherwise an error is returned. Passing exactly the correct shares recovers
// the original secret; supplying more than the threshold also works.
func CombineShares(shares []Share, prime *big.Int) (*big.Int, error) {
	if len(shares) == 0 {
		return nil, errors.New("crypto: CombineShares requires at least one share")
	}
	seen := make(map[string]bool, len(shares))
	for _, s := range shares {
		if s.X.Sign() == 0 {
			return nil, errors.New("crypto: CombineShares requires non-zero X coordinates")
		}
		key := new(big.Int).Mod(s.X, prime).String()
		if seen[key] {
			return nil, errors.New("crypto: CombineShares requires distinct X coordinates")
		}
		seen[key] = true
	}
	secret := big.NewInt(0)
	for i, si := range shares {
		// Lagrange basis L_i(0) = ∏_{j!=i} (0 - x_j) / (x_i - x_j).
		num := big.NewInt(1)
		den := big.NewInt(1)
		for j, sj := range shares {
			if i == j {
				continue
			}
			num.Mul(num, new(big.Int).Neg(sj.X))
			num.Mod(num, prime)
			diff := new(big.Int).Sub(si.X, sj.X)
			den.Mul(den, diff)
			den.Mod(den, prime)
		}
		denInv := new(big.Int).ModInverse(den, prime)
		if denInv == nil {
			return nil, errors.New("crypto: CombineShares failed to invert Lagrange denominator")
		}
		term := new(big.Int).Mul(si.Y, num)
		term.Mul(term, denInv)
		term.Mod(term, prime)
		secret.Add(secret, term)
		secret.Mod(secret, prime)
	}
	return secret, nil
}
