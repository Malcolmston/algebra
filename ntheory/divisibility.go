package ntheory

// abs64 returns the absolute value of x. It assumes x != math.MinInt64.
func abs64(x int64) int64 {
	if x < 0 {
		return -x
	}
	return x
}

// GCD returns the greatest common divisor of a and b.
//
// The result is always non-negative. By convention GCD(0, 0) == 0 and
// GCD(n, 0) == GCD(0, n) == |n|. Signs of the inputs are ignored.
func GCD(a, b int64) int64 {
	a, b = abs64(a), abs64(b)
	for b != 0 {
		a, b = b, a%b
	}
	return a
}

// LCM returns the least common multiple of a and b.
//
// The result is non-negative and LCM(a, 0) == LCM(0, b) == 0. The caller is
// responsible for ensuring the true least common multiple fits in an int64.
func LCM(a, b int64) int64 {
	if a == 0 || b == 0 {
		return 0
	}
	// Divide before multiplying to reduce the chance of overflow.
	return abs64(a/GCD(a, b)) * abs64(b)
}

// ExtendedGCD returns g = GCD(a, b) together with Bézout coefficients x and y
// satisfying a*x + b*y == g. The returned g is non-negative.
func ExtendedGCD(a, b int64) (g, x, y int64) {
	oldR, r := a, b
	oldS, s := int64(1), int64(0)
	oldT, t := int64(0), int64(1)
	for r != 0 {
		q := oldR / r
		oldR, r = r, oldR-q*r
		oldS, s = s, oldS-q*s
		oldT, t = t, oldT-q*t
	}
	g, x, y = oldR, oldS, oldT
	if g < 0 {
		g, x, y = -g, -x, -y
	}
	return g, x, y
}

// Divisors returns all positive divisors of n in ascending order.
//
// The sign of n is ignored. Divisors(0) returns nil because zero has
// infinitely many divisors.
func Divisors(n int64) []int64 {
	n = abs64(n)
	if n == 0 {
		return nil
	}
	var small, large []int64
	for d := int64(1); d*d <= n; d++ {
		if n%d == 0 {
			small = append(small, d)
			if q := n / d; q != d {
				large = append(large, q)
			}
		}
	}
	// large currently holds the big divisors in descending order; reverse it.
	for i, j := 0, len(large)-1; i < j; i, j = i+1, j-1 {
		large[i], large[j] = large[j], large[i]
	}
	return append(small, large...)
}

// SumDivisors returns σ(n), the sum of all positive divisors of n
// (including 1 and n itself). The sign of n is ignored. SumDivisors(0) == 0.
func SumDivisors(n int64) int64 {
	n = abs64(n)
	if n == 0 {
		return 0
	}
	// Use the multiplicative formula σ(p^a) = (p^(a+1) - 1)/(p - 1).
	var sum int64 = 1
	factors := Factorize(n)
	for p, a := range factors {
		term := int64(1)
		pk := int64(1)
		for i := 0; i < a; i++ {
			pk *= p
			term += pk
		}
		sum *= term
	}
	return sum
}

// CountDivisors returns τ(n), the number of positive divisors of n.
// The sign of n is ignored. CountDivisors(0) == 0.
func CountDivisors(n int64) int64 {
	n = abs64(n)
	if n == 0 {
		return 0
	}
	var count int64 = 1
	for _, a := range Factorize(n) {
		count *= int64(a + 1)
	}
	return count
}

// IsPerfect reports whether n is a perfect number, i.e. a positive integer
// equal to the sum of its proper divisors (σ(n) == 2n). The smallest perfect
// numbers are 6, 28, 496 and 8128.
func IsPerfect(n int64) bool {
	if n <= 0 {
		return false
	}
	return SumDivisors(n) == 2*n
}
