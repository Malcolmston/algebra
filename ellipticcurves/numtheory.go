package ellipticcurves

import (
	"errors"
	"math/big"
)

// ErrCRTInconsistent indicates that a system of congruences passed to CRT has
// no simultaneous solution.
var ErrCRTInconsistent = errors.New("ellipticcurves: inconsistent CRT system")

// KroneckerSymbol returns the Kronecker symbol (a/n), the total multiplicative
// extension of the Jacobi and Legendre symbols to all integers n, taking values
// in {-1, 0, 1}.
func KroneckerSymbol(a, n *big.Int) int {
	if n.Sign() == 0 {
		if a.CmpAbs(bigOne) == 0 {
			return 1
		}
		return 0
	}
	result := 1
	nn := new(big.Int).Set(n)
	if nn.Sign() < 0 {
		nn.Neg(nn)
		if a.Sign() < 0 {
			result = -result
		}
	}
	// Factor out powers of two from n.
	twos := 0
	for nn.Bit(0) == 0 {
		nn.Rsh(nn, 1)
		twos++
	}
	if twos > 0 {
		am8 := new(big.Int).Mod(a, big.NewInt(8)).Int64()
		if am8 < 0 {
			am8 += 8
		}
		var k int
		switch am8 {
		case 1, 7:
			k = 1
		case 3, 5:
			k = -1
		default:
			k = 0
		}
		if k == 0 {
			return 0
		}
		if twos%2 == 1 {
			result *= k
		}
	}
	if nn.Cmp(bigOne) == 0 {
		return result
	}
	j := JacobiSymbol(a, nn)
	return result * j
}

// CRT solves the system x = residues[i] (mod moduli[i]) for pairwise coprime
// moduli, returning the least non-negative solution and the product of the
// moduli. It returns ErrCRTInconsistent if the lengths differ or a modulus is
// not positive; for pairwise coprime moduli a solution always exists.
func CRT(residues, moduli []*big.Int) (*big.Int, *big.Int, error) {
	if len(residues) != len(moduli) || len(moduli) == 0 {
		return nil, nil, ErrCRTInconsistent
	}
	x := new(big.Int).Set(Mod(residues[0], moduli[0]))
	m := new(big.Int).Set(moduli[0])
	for i := 1; i < len(moduli); i++ {
		mi := moduli[i]
		if mi.Sign() <= 0 {
			return nil, nil, ErrCRTInconsistent
		}
		// Solve x + m*t = residues[i] (mod mi).
		rhs := new(big.Int).Sub(Mod(residues[i], mi), x)
		inv := new(big.Int).ModInverse(new(big.Int).Mod(m, mi), mi)
		if inv == nil {
			return nil, nil, ErrCRTInconsistent
		}
		t := new(big.Int).Mul(rhs, inv)
		t.Mod(t, mi)
		x.Add(x, new(big.Int).Mul(m, t))
		m.Mul(m, mi)
		x.Mod(x, m)
	}
	return x, m, nil
}

// IsPrimitiveRoot reports whether g is a primitive root modulo the odd prime p,
// i.e. its multiplicative order equals p-1.
func IsPrimitiveRoot(g, p *big.Int) bool {
	return MultiplicativeOrder(g, p).Cmp(new(big.Int).Sub(p, bigOne)) == 0
}

// ModPow returns a^e mod m for a non-negative exponent e, a convenience wrapper
// over big.Int.Exp that reduces a into [0, m) first.
func ModPow(a, e, m *big.Int) *big.Int {
	return new(big.Int).Exp(Mod(a, m), e, m)
}

// OrderModN returns the multiplicative order of a modulo n for a general modulus
// n > 1 with gcd(a, n) = 1, computed by iterated multiplication. It returns 0
// when a is not a unit modulo n.
func OrderModN(a, n *big.Int) *big.Int {
	if Gcd(a, n).Cmp(bigOne) != 0 {
		return big.NewInt(0)
	}
	am := Mod(a, n)
	cur := new(big.Int).Set(am)
	k := big.NewInt(1)
	for cur.Cmp(bigOne) != 0 {
		cur = ModMul(cur, am, n)
		k.Add(k, bigOne)
		if k.Cmp(n) > 0 {
			return big.NewInt(0)
		}
	}
	return k
}
