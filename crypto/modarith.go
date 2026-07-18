package crypto

import (
	"errors"
	"math/big"
)

// ErrNotInvertible is returned by ModInverse and related routines when the
// element has no multiplicative inverse modulo the given modulus (that is, when
// the element and modulus are not coprime).
var ErrNotInvertible = errors.New("crypto: element is not invertible modulo m")

// ModExp returns base**exp mod m using right-to-left binary
// (square-and-multiply) exponentiation. The modulus m must be positive.
//
// A negative exponent is supported provided base is invertible modulo m: the
// result is then (base^-1)^(-exp) mod m. base itself may be negative or larger
// than m; it is reduced into the range [0, m) first. ModExp panics if m <= 0
// and if exp < 0 but base has no inverse modulo m.
func ModExp(base, exp, m *big.Int) *big.Int {
	if m.Sign() <= 0 {
		panic("crypto: ModExp requires modulus m > 0")
	}
	if m.Cmp(cryptoOne) == 0 {
		return big.NewInt(0)
	}
	b := new(big.Int).Mod(base, m)
	e := new(big.Int).Set(exp)
	if e.Sign() < 0 {
		inv := new(big.Int).ModInverse(b, m)
		if inv == nil {
			panic("crypto: ModExp with negative exponent requires an invertible base")
		}
		b = inv
		e.Neg(e)
	}
	result := big.NewInt(1)
	for e.Sign() > 0 {
		if e.Bit(0) == 1 {
			result.Mul(result, b)
			result.Mod(result, m)
		}
		e.Rsh(e, 1)
		b.Mul(b, b)
		b.Mod(b, m)
	}
	return result
}

// ModInverse returns the multiplicative inverse of a modulo m, i.e. the unique
// x in [0, m) with a*x ≡ 1 (mod m). The modulus m must be positive. It returns
// ErrNotInvertible when gcd(a, m) != 1.
func ModInverse(a, m *big.Int) (*big.Int, error) {
	if m.Sign() <= 0 {
		panic("crypto: ModInverse requires modulus m > 0")
	}
	inv := new(big.Int).ModInverse(a, m)
	if inv == nil {
		return nil, ErrNotInvertible
	}
	return inv, nil
}

// ModMul returns a*b mod m. The modulus m must be positive; a and b may be
// negative and are reduced into [0, m).
func ModMul(a, b, m *big.Int) *big.Int {
	if m.Sign() <= 0 {
		panic("crypto: ModMul requires modulus m > 0")
	}
	r := new(big.Int).Mul(a, b)
	r.Mod(r, m)
	return r
}

// ModAdd returns a+b mod m. The modulus m must be positive.
func ModAdd(a, b, m *big.Int) *big.Int {
	if m.Sign() <= 0 {
		panic("crypto: ModAdd requires modulus m > 0")
	}
	r := new(big.Int).Add(a, b)
	r.Mod(r, m)
	return r
}

// ModSub returns a-b mod m, normalised into [0, m). The modulus m must be
// positive.
func ModSub(a, b, m *big.Int) *big.Int {
	if m.Sign() <= 0 {
		panic("crypto: ModSub requires modulus m > 0")
	}
	r := new(big.Int).Sub(a, b)
	r.Mod(r, m)
	return r
}

// GCD returns the non-negative greatest common divisor of a and b. GCD(0, 0)
// is 0. The arguments may be negative; the result is always non-negative.
func GCD(a, b *big.Int) *big.Int {
	return new(big.Int).GCD(nil, nil, new(big.Int).Abs(a), new(big.Int).Abs(b))
}

// LCM returns the non-negative least common multiple of a and b. LCM(x, 0) and
// LCM(0, x) are 0 by convention.
func LCM(a, b *big.Int) *big.Int {
	if a.Sign() == 0 || b.Sign() == 0 {
		return big.NewInt(0)
	}
	g := GCD(a, b)
	r := new(big.Int).Div(new(big.Int).Abs(a), g)
	r.Mul(r, new(big.Int).Abs(b))
	return r
}

// ExtendedGCD returns the greatest common divisor g of a and b together with
// Bézout coefficients x and y satisfying a*x + b*y = g. The returned g is
// non-negative.
func ExtendedGCD(a, b *big.Int) (g, x, y *big.Int) {
	g = new(big.Int)
	x = new(big.Int)
	y = new(big.Int)
	g.GCD(x, y, a, b)
	return g, x, y
}

// CRT solves a system of simultaneous congruences x ≡ residues[i] (mod
// moduli[i]) via the Chinese Remainder Theorem and returns the unique solution
// x in [0, M), where M is the product of all moduli. The two slices must be the
// same non-zero length and every modulus must be positive. It returns an error
// if the moduli are not pairwise coprime (in which case no unique solution
// exists over their product).
func CRT(residues, moduli []*big.Int) (*big.Int, error) {
	if len(residues) != len(moduli) {
		return nil, errors.New("crypto: CRT requires len(residues) == len(moduli)")
	}
	if len(moduli) == 0 {
		return nil, errors.New("crypto: CRT requires at least one congruence")
	}
	M := big.NewInt(1)
	for _, m := range moduli {
		if m.Sign() <= 0 {
			return nil, errors.New("crypto: CRT requires positive moduli")
		}
		M.Mul(M, m)
	}
	x := big.NewInt(0)
	for i, m := range moduli {
		Mi := new(big.Int).Div(M, m)
		inv := new(big.Int).ModInverse(Mi, m)
		if inv == nil {
			return nil, errors.New("crypto: CRT requires pairwise coprime moduli")
		}
		term := new(big.Int).Mod(residues[i], m)
		term.Mul(term, Mi)
		term.Mul(term, inv)
		x.Add(x, term)
	}
	x.Mod(x, M)
	return x, nil
}

// CRTPair solves the pair of congruences x ≡ a1 (mod n1) and x ≡ a2 (mod n2)
// and returns the unique solution in [0, n1*n2). Both moduli must be positive
// and coprime; otherwise an error is returned. It is a convenience wrapper
// around CRT for the common two-modulus case.
func CRTPair(a1, n1, a2, n2 *big.Int) (*big.Int, error) {
	return CRT([]*big.Int{a1, a2}, []*big.Int{n1, n2})
}
