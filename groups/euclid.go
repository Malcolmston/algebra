package groups

// groupsAbs returns the absolute value of an int.
func groupsAbs(a int) int {
	if a < 0 {
		return -a
	}
	return a
}

// Gcd returns the non-negative greatest common divisor of a and b using the
// Euclidean algorithm. Gcd(0, 0) == 0 and the result is invariant under the
// sign of either argument.
func Gcd(a, b int) int {
	a, b = groupsAbs(a), groupsAbs(b)
	for b != 0 {
		a, b = b, a%b
	}
	return a
}

// ExtendedGcd returns g = gcd(a, b) together with Bézout coefficients x and y
// satisfying a*x + b*y == g. The returned g is non-negative.
func ExtendedGcd(a, b int) (g, x, y int) {
	// Track coefficients for the current remainders.
	oldR, r := a, b
	oldS, s := 1, 0
	oldT, t := 0, 1
	for r != 0 {
		q := oldR / r
		oldR, r = r, oldR-q*r
		oldS, s = s, oldS-q*s
		oldT, t = t, oldT-q*t
	}
	if oldR < 0 {
		oldR, oldS, oldT = -oldR, -oldS, -oldT
	}
	return oldR, oldS, oldT
}

// Lcm returns the non-negative least common multiple of a and b. Lcm with a
// zero argument is 0.
func Lcm(a, b int) int {
	if a == 0 || b == 0 {
		return 0
	}
	return groupsAbs(a/Gcd(a, b)) * groupsAbs(b)
}

// GcdMany returns the greatest common divisor of all supplied integers.
// GcdMany with no arguments returns 0.
func GcdMany(nums ...int) int {
	g := 0
	for _, n := range nums {
		g = Gcd(g, n)
	}
	return g
}

// LcmMany returns the least common multiple of all supplied integers.
// LcmMany with no arguments returns 1 (the empty product). If any argument is
// 0 the result is 0.
func LcmMany(nums ...int) int {
	l := 1
	for _, n := range nums {
		l = Lcm(l, n)
	}
	return l
}

// Coprime reports whether a and b are relatively prime, i.e. gcd(a, b) == 1.
func Coprime(a, b int) bool {
	return Gcd(a, b) == 1
}
