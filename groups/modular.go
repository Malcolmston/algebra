package groups

// Mod returns the representative of a in Z/nZ lying in [0, n). The modulus n
// must be positive; unlike Go's built-in % operator the result is never
// negative. It panics if n <= 0.
func Mod(a, n int) int {
	if n <= 0 {
		panic("groups: Mod requires n > 0")
	}
	r := a % n
	if r < 0 {
		r += n
	}
	return r
}

// ModAdd returns (a + b) mod n in [0, n). The modulus n must be positive.
func ModAdd(a, b, n int) int {
	return Mod(Mod(a, n)+Mod(b, n), n)
}

// ModSub returns (a - b) mod n in [0, n). The modulus n must be positive.
func ModSub(a, b, n int) int {
	return Mod(Mod(a, n)-Mod(b, n), n)
}

// ModNeg returns (-a) mod n in [0, n). The modulus n must be positive.
func ModNeg(a, n int) int {
	return Mod(-a, n)
}

// ModMul returns (a * b) mod n in [0, n) computed without intermediate
// overflow for operands whose reduced values fit in an int. The modulus n must
// be positive.
func ModMul(a, b, n int) int {
	a = Mod(a, n)
	b = Mod(b, n)
	// Russian-peasant multiplication keeps every intermediate below 2n.
	result := 0
	for b > 0 {
		if b&1 == 1 {
			result = (result + a) % n
		}
		a = (a + a) % n
		b >>= 1
	}
	return result
}

// ModPow returns base^exp mod n in [0, n) via binary exponentiation. The
// exponent exp must be non-negative and the modulus n must be positive. It
// panics if exp < 0.
func ModPow(base, exp, n int) int {
	if exp < 0 {
		panic("groups: ModPow requires exp >= 0")
	}
	if n == 1 {
		return 0
	}
	base = Mod(base, n)
	result := 1
	for exp > 0 {
		if exp&1 == 1 {
			result = ModMul(result, base, n)
		}
		base = ModMul(base, base, n)
		exp >>= 1
	}
	return result
}

// ModInverse returns the multiplicative inverse of a modulo n in [0, n) and
// true when gcd(a, n) == 1. When a is not a unit modulo n it returns (0,
// false). The modulus n must be positive.
func ModInverse(a, n int) (int, bool) {
	if n <= 0 {
		panic("groups: ModInverse requires n > 0")
	}
	g, x, _ := ExtendedGcd(Mod(a, n), n)
	if g != 1 {
		return 0, false
	}
	return Mod(x, n), true
}

// CyclicOrder returns the order (number of elements) of the additive cyclic
// group Z/nZ, which is n. The argument n must be positive.
func CyclicOrder(n int) int {
	if n <= 0 {
		panic("groups: CyclicOrder requires n > 0")
	}
	return n
}

// ElementOrderZn returns the additive order of a in the cyclic group Z/nZ,
// that is the smallest k > 0 with k*a ≡ 0 (mod n). This equals n/gcd(a, n).
// The modulus n must be positive.
func ElementOrderZn(a, n int) int {
	if n <= 0 {
		panic("groups: ElementOrderZn requires n > 0")
	}
	return n / Gcd(Mod(a, n), n)
}

// MultiplicativeOrder returns the multiplicative order of a modulo n, the
// smallest k > 0 with a^k ≡ 1 (mod n), together with true. It returns (0,
// false) when gcd(a, n) != 1, since only units have a finite multiplicative
// order. The modulus n must be positive.
func MultiplicativeOrder(a, n int) (int, bool) {
	if n <= 0 {
		panic("groups: MultiplicativeOrder requires n > 0")
	}
	a = Mod(a, n)
	if Gcd(a, n) != 1 {
		return 0, false
	}
	if n == 1 {
		return 1, true
	}
	k := 1
	cur := a
	for cur != 1 {
		cur = ModMul(cur, a, n)
		k++
	}
	return k, true
}

// EulerTotient returns Euler's totient φ(n): the number of integers in
// [1, n] that are coprime to n, equivalently the order of the unit group
// (Z/nZ)*. The argument n must be positive.
func EulerTotient(n int) int {
	if n <= 0 {
		panic("groups: EulerTotient requires n > 0")
	}
	result := n
	m := n
	for p := 2; p*p <= m; p++ {
		if m%p == 0 {
			for m%p == 0 {
				m /= p
			}
			result -= result / p
		}
	}
	if m > 1 {
		result -= result / m
	}
	return result
}

// UnitsModN returns the elements of the multiplicative group (Z/nZ)* in
// increasing order: every a in [1, n) with gcd(a, n) == 1. For n == 1 it
// returns the single element {0} representing the trivial group. The argument
// n must be positive.
func UnitsModN(n int) []int {
	if n <= 0 {
		panic("groups: UnitsModN requires n > 0")
	}
	if n == 1 {
		return []int{0}
	}
	units := make([]int, 0, EulerTotient(n))
	for a := 1; a < n; a++ {
		if Gcd(a, n) == 1 {
			units = append(units, a)
		}
	}
	return units
}

// IsUnitModN reports whether a is a unit modulo n, i.e. gcd(a mod n, n) == 1.
// The modulus n must be positive.
func IsUnitModN(a, n int) bool {
	if n <= 0 {
		panic("groups: IsUnitModN requires n > 0")
	}
	return Gcd(Mod(a, n), n) == 1
}

// IsPrimitiveRoot reports whether g is a primitive root modulo n, i.e. g
// generates the whole unit group (Z/nZ)*. Equivalently the multiplicative
// order of g equals φ(n). The modulus n must be positive.
func IsPrimitiveRoot(g, n int) bool {
	if n <= 0 {
		panic("groups: IsPrimitiveRoot requires n > 0")
	}
	ord, ok := MultiplicativeOrder(g, n)
	if !ok {
		return false
	}
	return ord == EulerTotient(n)
}

// PrimitiveRoots returns every primitive root modulo n in [1, n) in increasing
// order. The slice is empty when no primitive root exists (which happens
// unless n is 1, 2, 4, p^k, or 2·p^k for an odd prime p). The argument n must
// be positive.
func PrimitiveRoots(n int) []int {
	if n <= 0 {
		panic("groups: PrimitiveRoots requires n > 0")
	}
	phi := EulerTotient(n)
	var roots []int
	for g := 1; g < n; g++ {
		if Gcd(g, n) != 1 {
			continue
		}
		if ord, _ := MultiplicativeOrder(g, n); ord == phi {
			roots = append(roots, g)
		}
	}
	return roots
}
