package contfrac

import "math/big"

// pellTerms builds the list of partial quotients a0, a1, ..., a_{n-1} of the
// continued fraction of sqrt(D) for the first n terms, given the head a0 and the
// repeating period.
func pellTerms(a0 int64, period []int64, n int) []int64 {
	terms := make([]int64, n)
	if n > 0 {
		terms[0] = a0
	}
	r := len(period)
	for i := 1; i < n; i++ {
		terms[i] = period[(i-1)%r]
	}
	return terms
}

// convergentBig returns the numerator and denominator of the last convergent of
// the given partial quotients, computed with big.Int.
func convergentBig(terms []int64) (*big.Int, *big.Int) {
	hPrev, hPrev2 := big.NewInt(1), big.NewInt(0)
	kPrev, kPrev2 := big.NewInt(0), big.NewInt(1)
	for _, a := range terms {
		ab := big.NewInt(a)
		h := new(big.Int).Add(new(big.Int).Mul(ab, hPrev), hPrev2)
		k := new(big.Int).Add(new(big.Int).Mul(ab, kPrev), kPrev2)
		hPrev, hPrev2 = h, hPrev
		kPrev, kPrev2 = k, kPrev
	}
	return hPrev, kPrev
}

// PellFundamental returns the fundamental (smallest positive) solution (x, y) of
// the Pell equation x^2 - D*y^2 = 1 for a non-square D > 1, together with ok. If
// D is a perfect square or D <= 0 there is no non-trivial solution and ok is
// false. The solution is derived from the continued fraction of sqrt(D).
func PellFundamental(D int64) (x, y *big.Int, ok bool) {
	if D <= 0 || IsPerfectSquare(D) {
		return nil, nil, false
	}
	sc := SqrtCF(D)
	a0 := sc.Head[0]
	period := sc.Period
	r := len(period)
	var length int
	if r%2 == 0 {
		length = r
	} else {
		length = 2 * r
	}
	terms := pellTerms(a0, period, length)
	x, y = convergentBig(terms)
	return x, y, true
}

// PellFundamentalInt returns the fundamental solution of x^2 - D*y^2 = 1 as
// int64 values when they fit, together with ok. ok is false when D admits no
// solution or the fundamental solution overflows int64.
func PellFundamentalInt(D int64) (x, y int64, ok bool) {
	bx, by, valid := PellFundamental(D)
	if !valid || !bx.IsInt64() || !by.IsInt64() {
		return 0, 0, false
	}
	return bx.Int64(), by.Int64(), true
}

// PellSolutions returns the first count positive solutions of x^2 - D*y^2 = 1,
// generated from the fundamental solution by the recurrence
// x_{k+1} = x1*x_k + D*y1*y_k, y_{k+1} = x1*y_k + y1*x_k. It returns nil when D
// admits no solution.
func PellSolutions(D int64, count int) [][2]*big.Int {
	x1, y1, ok := PellFundamental(D)
	if !ok || count <= 0 {
		return nil
	}
	Dbig := big.NewInt(D)
	out := make([][2]*big.Int, 0, count)
	xk := new(big.Int).Set(x1)
	yk := new(big.Int).Set(y1)
	for i := 0; i < count; i++ {
		out = append(out, [2]*big.Int{new(big.Int).Set(xk), new(big.Int).Set(yk)})
		// (x_{k+1}, y_{k+1})
		nx := new(big.Int).Add(new(big.Int).Mul(x1, xk), new(big.Int).Mul(Dbig, new(big.Int).Mul(y1, yk)))
		ny := new(big.Int).Add(new(big.Int).Mul(x1, yk), new(big.Int).Mul(y1, xk))
		xk, yk = nx, ny
	}
	return out
}

// PellNthSolution returns the n-th positive solution (n >= 1) of
// x^2 - D*y^2 = 1, with n == 1 the fundamental solution. It returns ok == false
// when D admits no solution or n < 1.
func PellNthSolution(D int64, n int) (x, y *big.Int, ok bool) {
	if n < 1 {
		return nil, nil, false
	}
	sols := PellSolutions(D, n)
	if sols == nil {
		return nil, nil, false
	}
	return sols[n-1][0], sols[n-1][1], true
}

// PellNegative returns the fundamental solution (x, y) of the negative Pell
// equation x^2 - D*y^2 = -1, together with ok. Such a solution exists precisely
// when the period of the continued fraction of sqrt(D) has odd length; ok is
// false otherwise.
func PellNegative(D int64) (x, y *big.Int, ok bool) {
	if D <= 0 || IsPerfectSquare(D) {
		return nil, nil, false
	}
	sc := SqrtCF(D)
	period := sc.Period
	r := len(period)
	if r%2 == 0 {
		return nil, nil, false
	}
	terms := pellTerms(sc.Head[0], period, r) // convergent index r-1
	x, y = convergentBig(terms)
	return x, y, true
}

// PellNegativeSolvable reports whether x^2 - D*y^2 = -1 has a solution, i.e. the
// period of sqrt(D) has odd length.
func PellNegativeSolvable(D int64) bool {
	if D <= 0 || IsPerfectSquare(D) {
		return false
	}
	return SqrtCFPeriodLength(D)%2 == 1
}

// IsPellSolution reports whether (x, y) satisfies x^2 - D*y^2 == n.
func IsPellSolution(D, x, y, n int64) bool {
	return x*x-D*y*y == n
}
