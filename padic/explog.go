package padic

import "math/big"

// Log returns the p-adic logarithm of x, defined by the series
//
//	log(1 + t) = sum_{n>=1} (-1)^(n+1) t^n / n,   t = x - 1,
//
// which converges when x is a 1-unit, i.e. x == 1 mod p (valuation of x-1 at
// least 1). It returns ErrDomain otherwise. The result is computed to the
// absolute precision of x.
func (x *Padic) Log() (*Padic, error) {
	one := One(x.p, maxInt(x.prec, 1))
	t, err := x.Sub(one)
	if err != nil {
		return nil, err
	}
	if t.IsZero() {
		return newZero(x.p, x.AbsolutePrecision()), nil
	}
	if t.val < 1 {
		return nil, ErrDomain
	}
	return logSeries(t)
}

// Log1p returns log(1 + t) for a p-adic t with valuation at least 1.
func Log1p(t *Padic) (*Padic, error) {
	if t.IsZero() {
		return newZero(t.p, t.AbsolutePrecision()), nil
	}
	if t.val < 1 {
		return nil, ErrDomain
	}
	return logSeries(t)
}

// logSeries sums the log series for t with val(t) >= 1 up to the absolute
// precision of t. The term valuations are not monotonic (they dip whenever n is
// divisible by p), so every term below the target precision is summed rather
// than stopping at the first small one.
func logSeries(t *Padic) (*Padic, error) {
	target := t.AbsolutePrecision()
	nMax := target + 128
	sum := newZero(t.p, target)
	power := t.Copy() // t^n, starting n=1
	var err error
	for n := 1; n <= nMax; n++ {
		// term = (-1)^(n+1) * power / n
		nP := FromInt(t.p, n, target+valGuard(t.p, n)+2)
		termN, e := power.Div(nP)
		if e != nil {
			return nil, e
		}
		if n%2 == 0 {
			termN = termN.Neg()
		}
		if !termN.IsZero() && termN.val < target {
			sum, err = sum.Add(termN)
			if err != nil {
				return nil, err
			}
		}
		power, err = power.Mul(t)
		if err != nil {
			return nil, err
		}
	}
	return sum.ReduceTo(target), nil
}

// valGuard returns the p-adic valuation of n, used to enlarge working precision
// so that division by n does not silently lose digits.
func valGuard(p *big.Int, n int) int {
	if n == 0 {
		return 0
	}
	return ValuationInt(p, big.NewInt(int64(n)))
}

// Exp returns the p-adic exponential of x, defined by sum_{n>=0} x^n / n!,
// which converges when val(x) > 1/(p-1): for odd p this means val(x) >= 1, and
// for p = 2 it means val(x) >= 2. It returns ErrDomain outside that region.
func (x *Padic) Exp() (*Padic, error) {
	minVal := 1
	if x.p.Cmp(bigTwo) == 0 {
		minVal = 2
	}
	if x.IsZero() {
		return One(x.p, maxInt(x.AbsolutePrecision(), 1)), nil
	}
	if x.val < minVal {
		return nil, ErrDomain
	}
	target := x.AbsolutePrecision()
	nMax := 2*target + 16
	sum := One(x.p, target)
	term := One(x.p, target) // running x^n / n!
	var err error
	for n := 1; n <= nMax; n++ {
		// term *= x / n
		term, err = term.Mul(x)
		if err != nil {
			return nil, err
		}
		nP := FromInt(x.p, n, target+valGuard(x.p, n)+2)
		term, err = term.Div(nP)
		if err != nil {
			return nil, err
		}
		if !term.IsZero() && term.val < target {
			sum, err = sum.Add(term)
			if err != nil {
				return nil, err
			}
		}
	}
	return sum.ReduceTo(target), nil
}

// ExpConverges reports whether Exp(x) converges, i.e. val(x) > 1/(p-1).
func (x *Padic) ExpConverges() bool {
	if x.IsZero() {
		return true
	}
	if x.p.Cmp(bigTwo) == 0 {
		return x.val >= 2
	}
	return x.val >= 1
}

// LogConverges reports whether Log(x) converges, i.e. x is a 1-unit.
func (x *Padic) LogConverges() bool {
	one := One(x.p, maxInt(x.prec, 1))
	t, err := x.Sub(one)
	if err != nil {
		return false
	}
	return t.IsZero() || t.val >= 1
}
