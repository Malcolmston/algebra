package groups

// groupsPolyModTrim removes trailing zero coefficients (mod p) and reduces every
// remaining coefficient into [0, p).
func groupsPolyModTrim(a []int, p int) []int {
	out := make([]int, len(a))
	for i, v := range a {
		out[i] = Mod(v, p)
	}
	n := len(out)
	for n > 0 && out[n-1] == 0 {
		n--
	}
	return out[:n]
}

// PolyModDegree returns the degree of the polynomial a over GF(p): the index of
// its highest non-zero coefficient after reduction mod p, or -1 for the zero
// polynomial. p must be prime.
func PolyModDegree(a []int, p int) int {
	groupsRequirePrime(p, "PolyModDegree")
	return len(groupsPolyModTrim(a, p)) - 1
}

// PolyModAdd returns the sum a + b of two polynomials over GF(p), with
// coefficients reduced into [0, p). p must be prime.
func PolyModAdd(a, b []int, p int) []int {
	groupsRequirePrime(p, "PolyModAdd")
	n := len(a)
	if len(b) > n {
		n = len(b)
	}
	out := make([]int, n)
	for i := 0; i < n; i++ {
		var av, bv int
		if i < len(a) {
			av = a[i]
		}
		if i < len(b) {
			bv = b[i]
		}
		out[i] = ModAdd(av, bv, p)
	}
	return groupsPolyModTrim(out, p)
}

// PolyModSub returns the difference a - b of two polynomials over GF(p). p must
// be prime.
func PolyModSub(a, b []int, p int) []int {
	groupsRequirePrime(p, "PolyModSub")
	n := len(a)
	if len(b) > n {
		n = len(b)
	}
	out := make([]int, n)
	for i := 0; i < n; i++ {
		var av, bv int
		if i < len(a) {
			av = a[i]
		}
		if i < len(b) {
			bv = b[i]
		}
		out[i] = ModSub(av, bv, p)
	}
	return groupsPolyModTrim(out, p)
}

// PolyModMul returns the product a·b of two polynomials over GF(p). p must be
// prime.
func PolyModMul(a, b []int, p int) []int {
	groupsRequirePrime(p, "PolyModMul")
	a = groupsPolyModTrim(a, p)
	b = groupsPolyModTrim(b, p)
	if len(a) == 0 || len(b) == 0 {
		return []int{}
	}
	out := make([]int, len(a)+len(b)-1)
	for i, av := range a {
		for j, bv := range b {
			out[i+j] = ModAdd(out[i+j], ModMul(av, bv, p), p)
		}
	}
	return groupsPolyModTrim(out, p)
}

// PolyModDivMod returns the quotient q and remainder r of dividing a by b over
// GF(p), satisfying a = q·b + r with deg(r) < deg(b). The divisor b must be
// non-zero mod p; it panics otherwise. p must be prime.
func PolyModDivMod(a, b []int, p int) (q, r []int) {
	groupsRequirePrime(p, "PolyModDivMod")
	b = groupsPolyModTrim(b, p)
	if len(b) == 0 {
		panic("groups: PolyModDivMod division by zero polynomial")
	}
	rem := groupsPolyModTrim(a, p)
	db := len(b) - 1
	invLc := GFInv(b[db], p)
	quo := []int{}
	for len(rem)-1 >= db && len(rem) > 0 {
		dr := len(rem) - 1
		coeff := ModMul(rem[dr], invLc, p)
		shift := dr - db
		for len(quo) <= shift {
			quo = append(quo, 0)
		}
		quo[shift] = coeff
		for i := 0; i <= db; i++ {
			rem[shift+i] = ModSub(rem[shift+i], ModMul(coeff, b[i], p), p)
		}
		rem = groupsPolyModTrim(rem, p)
	}
	return groupsPolyModTrim(quo, p), groupsPolyModTrim(rem, p)
}

// PolyModGCD returns the monic greatest common divisor of a and b over GF(p)
// via the Euclidean algorithm. The GCD of two zero polynomials is the zero
// polynomial. p must be prime.
func PolyModGCD(a, b []int, p int) []int {
	groupsRequirePrime(p, "PolyModGCD")
	a = groupsPolyModTrim(a, p)
	b = groupsPolyModTrim(b, p)
	for len(b) > 0 {
		_, r := PolyModDivMod(a, b, p)
		a, b = b, r
	}
	if len(a) == 0 {
		return []int{}
	}
	return PolyModMonic(a, p)
}

// PolyModMonic returns a scaled so its leading coefficient is 1 over GF(p). The
// zero polynomial is returned unchanged. p must be prime.
func PolyModMonic(a []int, p int) []int {
	groupsRequirePrime(p, "PolyModMonic")
	a = groupsPolyModTrim(a, p)
	if len(a) == 0 {
		return []int{}
	}
	inv := GFInv(a[len(a)-1], p)
	out := make([]int, len(a))
	for i, v := range a {
		out[i] = ModMul(v, inv, p)
	}
	return out
}

// PolyModEval returns the value a(x) mod p evaluated at x over GF(p) using
// Horner's method. p must be prime.
func PolyModEval(a []int, x, p int) int {
	groupsRequirePrime(p, "PolyModEval")
	result := 0
	for i := len(a) - 1; i >= 0; i-- {
		result = ModAdd(ModMul(result, x, p), a[i], p)
	}
	return result
}
