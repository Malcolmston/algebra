package padic

import "math/big"

// SqrtUnitMod returns a square root of the unit u modulo p^prec, i.e. an
// integer r in [0, p^prec) with r^2 == u mod p^prec, or ErrNoRoot if u is not
// a square. u must be coprime to p. For odd p the root is found modulo p and
// Hensel lifted; for p = 2 it is built bit by bit and requires u == 1 mod 8.
func SqrtUnitMod(p, u *big.Int, prec int) (*big.Int, error) {
	if prec <= 0 {
		return nil, ErrPrecision
	}
	mod := PPow(p, prec)
	um := new(big.Int).Mod(u, mod)
	if new(big.Int).Mod(um, p).Sign() == 0 {
		return nil, ErrDomain
	}
	if p.Cmp(bigTwo) == 0 {
		return sqrtUnit2(um, prec)
	}
	// odd p: root mod p then Hensel lift of X^2 - u.
	r0, err := SqrtModP(um, p, nil)
	if err != nil {
		return nil, err
	}
	f := []*big.Int{new(big.Int).Neg(um), big.NewInt(0), big.NewInt(1)} // -u + x^2
	return HenselLift(f, r0, p, prec)
}

// sqrtUnit2 returns a square root of the odd unit u modulo 2^prec, requiring
// u == 1 mod 8 (or, for prec < 3, u odd), built one bit at a time.
func sqrtUnit2(u *big.Int, prec int) (*big.Int, error) {
	mod := PPow(bigTwo, prec)
	um := new(big.Int).Mod(u, mod)
	if prec >= 3 && new(big.Int).And(um, big.NewInt(7)).Int64() != 1 {
		return nil, ErrNoRoot
	}
	y := big.NewInt(1)
	for k := 3; k < prec; k++ {
		y2 := new(big.Int).Mul(y, y)
		diff := new(big.Int).Sub(y2, um)
		// c = diff / 2^k
		c := new(big.Int).Rsh(diff, uint(k))
		if c.Bit(0) == 1 {
			y.Add(y, new(big.Int).Lsh(bigOne, uint(k-1)))
		}
	}
	y.Mod(y, mod)
	return y, nil
}

// Sqrt returns a square root of the p-adic number x to relative precision equal
// to that of x, or ErrNoRoot if x is not a square in Q_p. The valuation of x
// must be even. Of the two roots (r and -r) the one with the smaller canonical
// unit is returned; use Neg for the other.
func (x *Padic) Sqrt() (*Padic, error) {
	if x.IsZero() {
		return newZero(x.p, (x.val+1)/2), nil
	}
	if x.val%2 != 0 {
		return nil, ErrNoRoot
	}
	root, err := SqrtUnitMod(x.p, x.unit, x.prec)
	if err != nil {
		return nil, err
	}
	return &Padic{p: new(big.Int).Set(x.p), val: x.val / 2, unit: root, prec: x.prec}, nil
}

// IsSquare reports whether x is a square in Q_p, to the extent its tracked
// precision allows. It requires even valuation and a square unit (for p = 2,
// unit == 1 mod 8, needing at least 3 digits of relative precision).
func (x *Padic) IsSquare() bool {
	if x.IsZero() {
		return true
	}
	if x.val%2 != 0 {
		return false
	}
	if x.p.Cmp(bigTwo) == 0 {
		if x.prec < 3 {
			return x.unit.Bit(0) == 1
		}
		return new(big.Int).And(x.unit, big.NewInt(7)).Int64() == 1
	}
	return LegendreSymbol(x.unit, x.p) == 1
}

// SqrtInt returns a square root of the integer a modulo p^prec (an integer r in
// [0, p^prec) with r^2 == a mod p^prec), or an error when none exists. a may be
// divisible by p, in which case its p-valuation must be even and below prec.
func SqrtInt(p, a *big.Int, prec int) (*big.Int, error) {
	if prec <= 0 {
		return nil, ErrPrecision
	}
	if a.Sign() == 0 {
		return big.NewInt(0), nil
	}
	v := ValuationInt(p, a)
	if v%2 != 0 {
		return nil, ErrNoRoot
	}
	u := UnitPartInt(p, a)
	ru, err := SqrtUnitMod(p, u, prec)
	if err != nil {
		return nil, err
	}
	root := new(big.Int).Mul(ru, PPow(p, v/2))
	root.Mod(root, PPow(p, prec))
	return root, nil
}

// SqrtBothRoots returns both square roots (r and -r) of x, or an error if x is
// not a square. For x = 0 it returns a single zero root.
func (x *Padic) SqrtBothRoots() ([]*Padic, error) {
	r, err := x.Sqrt()
	if err != nil {
		return nil, err
	}
	if x.IsZero() {
		return []*Padic{r}, nil
	}
	return []*Padic{r, r.Neg()}, nil
}
