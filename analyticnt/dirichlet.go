package analyticnt

import (
	"math"
	"math/cmplx"
)

// DirichletCharacter represents a Dirichlet character χ modulo q: a completely
// multiplicative, q-periodic function on the integers, defined by its values on
// a reduced residue system. The zero-valued Values entries mark residues not
// coprime to q (where χ = 0).
type DirichletCharacter struct {
	Q      int64        // the modulus
	Values []complex128 // Values[a] = χ(a) for 0 <= a < q
}

// Eval returns χ(n), reducing n modulo q. For n not coprime to q the result is
// 0.
func (c DirichletCharacter) Eval(n int64) complex128 {
	m := n % c.Q
	if m < 0 {
		m += c.Q
	}
	return c.Values[m]
}

// IsPrincipal reports whether χ is the principal character χ₀ modulo q, which is
// 1 on residues coprime to q and 0 otherwise.
func (c DirichletCharacter) IsPrincipal() bool {
	for a := int64(0); a < c.Q; a++ {
		if gcdInt64(a, c.Q) == 1 {
			if cmplx.Abs(c.Values[a]-1) > 1e-9 {
				return false
			}
		}
	}
	return true
}

// Order returns the multiplicative order of χ in the character group, i.e. the
// least k >= 1 with χ^k = χ₀.
func (c DirichletCharacter) Order() int {
	// χ(a) is a root of unity of order dividing φ(q); find lcm of orders.
	order := 1
	for a := int64(0); a < c.Q; a++ {
		if gcdInt64(a, c.Q) != 1 {
			continue
		}
		v := c.Values[a]
		k := 1
		acc := v
		for cmplx.Abs(acc-1) > 1e-7 && k < 10000 {
			acc *= v
			k++
		}
		order = lcmInt(order, k)
	}
	return order
}

// Conjugate returns the complex-conjugate character χ̄, which is also the
// inverse character χ^{-1}.
func (c DirichletCharacter) Conjugate() DirichletCharacter {
	vals := make([]complex128, len(c.Values))
	for i, v := range c.Values {
		vals[i] = cmplx.Conj(v)
	}
	return DirichletCharacter{Q: c.Q, Values: vals}
}

// IsReal reports whether χ takes only real values (±1, 0), i.e. χ is a
// quadratic or principal character.
func (c DirichletCharacter) IsReal() bool {
	for _, v := range c.Values {
		if math.Abs(imag(v)) > 1e-9 {
			return false
		}
	}
	return true
}

// PrincipalCharacter returns the principal Dirichlet character χ₀ modulo q.
func PrincipalCharacter(q int64) DirichletCharacter {
	if q < 1 {
		panic("analyticnt: PrincipalCharacter requires q >= 1")
	}
	vals := make([]complex128, q)
	for a := int64(0); a < q; a++ {
		if gcdInt64(a, q) == 1 {
			vals[a] = 1
		}
	}
	return DirichletCharacter{Q: q, Values: vals}
}

// LegendreCharacter returns the quadratic character modulo an odd prime p given
// by the Legendre symbol (a|p): +1 for quadratic residues, −1 for
// non-residues, 0 for multiples of p.
func LegendreCharacter(p int64) DirichletCharacter {
	if p < 3 || !IsPrime(p) {
		panic("analyticnt: LegendreCharacter requires an odd prime p")
	}
	vals := make([]complex128, p)
	for a := int64(1); a < p; a++ {
		vals[a] = complex(float64(LegendreSymbol(a, p)), 0)
	}
	return DirichletCharacter{Q: p, Values: vals}
}

// LegendreSymbol returns the Legendre symbol (a|p) for an odd prime p using
// Euler's criterion, evaluated as a^{(p−1)/2} mod p mapped to {−1, 0, 1}.
func LegendreSymbol(a, p int64) int {
	a %= p
	if a < 0 {
		a += p
	}
	if a == 0 {
		return 0
	}
	r := modPow(a, (p-1)/2, p)
	if r == 1 {
		return 1
	}
	return -1
}

// JacobiSymbol returns the Jacobi symbol (a|n) for odd n >= 1, generalizing the
// Legendre symbol multiplicatively over the prime factorization of n.
func JacobiSymbol(a, n int64) int {
	if n <= 0 || n%2 == 0 {
		panic("analyticnt: JacobiSymbol requires odd n >= 1")
	}
	a %= n
	if a < 0 {
		a += n
	}
	result := 1
	for a != 0 {
		for a%2 == 0 {
			a /= 2
			if n%8 == 3 || n%8 == 5 {
				result = -result
			}
		}
		a, n = n, a
		if a%4 == 3 && n%4 == 3 {
			result = -result
		}
		a %= n
	}
	if n == 1 {
		return result
	}
	return 0
}

// CharacterFromPrimitiveRoot builds a Dirichlet character modulo q with cyclic
// group (Z/qZ)* by sending a fixed primitive root g to e^{2πi·j/φ(q)}. The index
// j (0 <= j < φ(q)) selects which character in the cyclic dual group is
// returned. q must have a primitive root (q = 1,2,4,p^k,2p^k).
func CharacterFromPrimitiveRoot(q int64, j int) DirichletCharacter {
	if q < 1 {
		panic("analyticnt: CharacterFromPrimitiveRoot requires q >= 1")
	}
	phi := EulerPhi(q)
	g, ok := primitiveRoot(q)
	if !ok {
		panic("analyticnt: modulus has no primitive root")
	}
	vals := make([]complex128, q)
	// Discrete log table: index of each unit relative to g.
	cur := int64(1)
	for e := int64(0); e < phi; e++ {
		angle := 2 * math.Pi * float64(j) * float64(e) / float64(phi)
		vals[cur] = cmplx.Rect(1, angle)
		cur = (cur * g) % q
	}
	return DirichletCharacter{Q: q, Values: vals}
}

// AllCharactersModulo returns the full group of φ(q) Dirichlet characters modulo
// q when (Z/qZ)* is cyclic. For non-cyclic moduli it returns only the principal
// character together with as many cyclic-root characters as it can build, so
// callers requiring completeness should restrict to cyclic q.
func AllCharactersModulo(q int64) []DirichletCharacter {
	if _, ok := primitiveRoot(q); ok {
		phi := EulerPhi(q)
		out := make([]DirichletCharacter, 0, phi)
		for j := int64(0); j < phi; j++ {
			out = append(out, CharacterFromPrimitiveRoot(q, int(j)))
		}
		return out
	}
	return []DirichletCharacter{PrincipalCharacter(q)}
}

// gcdInt64 returns the greatest common divisor of a and b.
func gcdInt64(a, b int64) int64 {
	if a < 0 {
		a = -a
	}
	if b < 0 {
		b = -b
	}
	for b != 0 {
		a, b = b, a%b
	}
	return a
}

// lcmInt returns the least common multiple of two ints.
func lcmInt(a, b int) int {
	if a == 0 || b == 0 {
		return 0
	}
	g := int(gcdInt64(int64(a), int64(b)))
	return a / g * b
}

// modPow returns base^exp mod m for m > 0 using fast exponentiation.
func modPow(base, exp, m int64) int64 {
	if m == 1 {
		return 0
	}
	base %= m
	if base < 0 {
		base += m
	}
	result := int64(1)
	for exp > 0 {
		if exp&1 == 1 {
			result = mulMod(result, base, m)
		}
		base = mulMod(base, base, m)
		exp >>= 1
	}
	return result
}

// mulMod returns a*b mod m avoiding overflow via 128-bit-style reduction using
// big-free math for moduli up to ~3e9; for larger moduli it falls back on the
// direct product (callers here use small moduli).
func mulMod(a, b, m int64) int64 {
	if a < (1<<31) && b < (1<<31) {
		return (a * b) % m
	}
	// Russian-peasant multiplication modulo m to avoid overflow.
	var result int64
	a %= m
	for b > 0 {
		if b&1 == 1 {
			result = (result + a) % m
		}
		a = (a * 2) % m
		b >>= 1
	}
	return result
}

// primitiveRoot returns a primitive root modulo q and true, or (0, false) if
// none exists.
func primitiveRoot(q int64) (int64, bool) {
	if q == 1 {
		return 0, true
	}
	if q == 2 {
		return 1, true
	}
	if q == 4 {
		return 3, true
	}
	phi := EulerPhi(q)
	// Factor phi.
	factors := distinctPrimeFactors(phi)
	for g := int64(2); g < q; g++ {
		if gcdInt64(g, q) != 1 {
			continue
		}
		ok := true
		for _, p := range factors {
			if modPow(g, phi/p, q) == 1 {
				ok = false
				break
			}
		}
		if ok {
			return g, true
		}
	}
	return 0, false
}

// distinctPrimeFactors returns the distinct prime factors of n.
func distinctPrimeFactors(n int64) []int64 {
	var out []int64
	for p := int64(2); p*p <= n; p++ {
		if n%p == 0 {
			out = append(out, p)
			for n%p == 0 {
				n /= p
			}
		}
	}
	if n > 1 {
		out = append(out, n)
	}
	return out
}
