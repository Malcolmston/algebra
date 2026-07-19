package contfrac

import "math/big"

// ECFTerm returns the k-th partial quotient (0-indexed) of the continued
// fraction of Euler's number e = [2; 1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8, ...].
// After the leading 2 the terms follow the pattern 1, 2m, 1 for m = 1, 2, 3, ...
func ECFTerm(k int) int64 {
	if k < 0 {
		panic("contfrac: ECFTerm requires k >= 0")
	}
	if k == 0 {
		return 2
	}
	m := (k + 2) / 3    // block index, 1-based
	if k-(3*m-2) == 1 { // middle element of the block
		return int64(2 * m)
	}
	return 1
}

// ECF returns the first n partial quotients of the continued fraction of e.
func ECF(n int) CF {
	if n < 0 {
		n = 0
	}
	cf := make(CF, n)
	for k := 0; k < n; k++ {
		cf[k] = ECFTerm(k)
	}
	return cf
}

// EConvergents returns the first n convergents of e as reduced fractions.
func EConvergents(n int) []Frac {
	return ECF(n).Convergents()
}

// EApprox returns the n-th convergent of e (using n partial quotients) as a
// single best rational approximation.
func EApprox(n int) Frac {
	if n < 1 {
		n = 1
	}
	return ECF(n).Frac()
}

// GoldenRatioCF returns the first n partial quotients of the continued fraction
// of the golden ratio phi = (1 + sqrt(5))/2 = [1; 1, 1, 1, ...], which is all
// ones.
func GoldenRatioCF(n int) CF {
	if n < 0 {
		n = 0
	}
	cf := make(CF, n)
	for i := range cf {
		cf[i] = 1
	}
	return cf
}

// GoldenRatio returns the golden ratio (1 + sqrt(5))/2 as a float64.
func GoldenRatio() float64 {
	return NewSurd(1, 2, 5).Value()
}

// SqrtTwoCF returns the first n partial quotients of the continued fraction of
// sqrt(2) = [1; 2, 2, 2, ...].
func SqrtTwoCF(n int) CF {
	return SqrtCFExpand(2, n)
}

// piBigFloat computes pi to the given binary precision using the Machin formula
// pi = 16*atan(1/5) - 4*atan(1/239).
func piBigFloat(prec uint) *big.Float {
	a := atanInv(5, prec)
	b := atanInv(239, prec)
	sixteen := new(big.Float).SetPrec(prec).SetInt64(16)
	four := new(big.Float).SetPrec(prec).SetInt64(4)
	term1 := new(big.Float).SetPrec(prec).Mul(sixteen, a)
	term2 := new(big.Float).SetPrec(prec).Mul(four, b)
	return new(big.Float).SetPrec(prec).Sub(term1, term2)
}

// atanInv returns arctan(1/x) to the given precision via its Taylor series.
func atanInv(x int64, prec uint) *big.Float {
	xf := new(big.Float).SetPrec(prec).SetInt64(x)
	x2 := new(big.Float).SetPrec(prec).Mul(xf, xf)
	// term = 1/x, sum = term, then term *= -1/x^2, divide by (2k+1)
	term := new(big.Float).SetPrec(prec).Quo(new(big.Float).SetPrec(prec).SetInt64(1), xf)
	sum := new(big.Float).SetPrec(prec).Set(term)
	sign := int64(-1)
	for k := int64(1); ; k++ {
		term = new(big.Float).SetPrec(prec).Quo(term, x2)
		denom := new(big.Float).SetPrec(prec).SetInt64(2*k + 1)
		t := new(big.Float).SetPrec(prec).Quo(term, denom)
		if sign < 0 {
			sum.Sub(sum, t)
		} else {
			sum.Add(sum, t)
		}
		sign = -sign
		// Stop when the term is below the working precision.
		if t.MantExp(nil) < -int(prec)-4 {
			break
		}
	}
	return sum
}

// PiFloat returns pi computed to roughly the requested number of decimal digits
// as a *big.Float. digits must be positive.
func PiFloat(digits int) *big.Float {
	if digits < 1 {
		digits = 1
	}
	prec := uint(digits*4 + 64)
	return piBigFloat(prec)
}

// PiCF returns the first n partial quotients of the continued fraction of pi,
// which begins [3; 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, ...]. The terms are
// extracted from a high-precision computation of pi, so the result is exact for
// the number of terms requested.
func PiCF(n int) CF {
	if n <= 0 {
		return CF{}
	}
	prec := uint(64*n + 256)
	pi := piBigFloat(prec)
	return FromBigFloat(pi, n)
}

// PiConvergents returns the first n convergents of pi as reduced fractions,
// including the classical 22/7 and 355/113.
func PiConvergents(n int) []Frac {
	return PiCF(n).Convergents()
}

// PiApprox returns the continued-fraction convergent of pi obtained from n
// partial quotients as a single best rational approximation.
func PiApprox(n int) Frac {
	if n < 1 {
		n = 1
	}
	return PiCF(n).Frac()
}
