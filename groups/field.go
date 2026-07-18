package groups

// IsPrime reports whether n is a prime number using trial division. Values
// less than 2 are not prime.
func IsPrime(n int) bool {
	if n < 2 {
		return false
	}
	if n%2 == 0 {
		return n == 2
	}
	if n%3 == 0 {
		return n == 3
	}
	for i := 5; i*i <= n; i += 6 {
		if n%i == 0 || n%(i+2) == 0 {
			return false
		}
	}
	return true
}

// groupsRequirePrime panics unless p is prime, naming the calling operation.
func groupsRequirePrime(p int, op string) {
	if !IsPrime(p) {
		panic("groups: " + op + " requires a prime modulus p")
	}
}

// GFAdd returns a + b in the prime field GF(p). p must be prime.
func GFAdd(a, b, p int) int {
	groupsRequirePrime(p, "GFAdd")
	return ModAdd(a, b, p)
}

// GFSub returns a - b in the prime field GF(p). p must be prime.
func GFSub(a, b, p int) int {
	groupsRequirePrime(p, "GFSub")
	return ModSub(a, b, p)
}

// GFNeg returns the additive inverse -a in the prime field GF(p). p must be
// prime.
func GFNeg(a, p int) int {
	groupsRequirePrime(p, "GFNeg")
	return ModNeg(a, p)
}

// GFMul returns a * b in the prime field GF(p). p must be prime.
func GFMul(a, b, p int) int {
	groupsRequirePrime(p, "GFMul")
	return ModMul(a, b, p)
}

// GFInv returns the multiplicative inverse of a in the prime field GF(p). It
// panics if a ≡ 0 (mod p), which has no inverse. p must be prime.
func GFInv(a, p int) int {
	groupsRequirePrime(p, "GFInv")
	if Mod(a, p) == 0 {
		panic("groups: GFInv of zero")
	}
	// Fermat's little theorem: a^(p-2) ≡ a^-1 (mod p).
	return ModPow(a, p-2, p)
}

// GFDiv returns a / b in the prime field GF(p). It panics if b ≡ 0 (mod p). p
// must be prime.
func GFDiv(a, b, p int) int {
	groupsRequirePrime(p, "GFDiv")
	return GFMul(a, GFInv(b, p), p)
}

// GFPow returns a raised to the power e in the prime field GF(p). The exponent
// e must be non-negative and p must be prime.
func GFPow(a, e, p int) int {
	groupsRequirePrime(p, "GFPow")
	return ModPow(a, e, p)
}
