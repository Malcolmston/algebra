package padic

import "math/big"

// polyEvalMod evaluates the integer polynomial (coeffs low-to-high) at x modulo
// m using Horner's rule.
func polyEvalMod(coeffs []*big.Int, x, m *big.Int) *big.Int {
	acc := new(big.Int)
	for i := len(coeffs) - 1; i >= 0; i-- {
		acc.Mul(acc, x)
		acc.Add(acc, coeffs[i])
		acc.Mod(acc, m)
	}
	return acc
}

// polyDerivative returns the formal derivative of the integer polynomial given
// low-to-high.
func polyDerivative(coeffs []*big.Int) []*big.Int {
	if len(coeffs) <= 1 {
		return []*big.Int{big.NewInt(0)}
	}
	d := make([]*big.Int, len(coeffs)-1)
	for i := 1; i < len(coeffs); i++ {
		d[i-1] = new(big.Int).Mul(coeffs[i], big.NewInt(int64(i)))
	}
	return d
}

// PolyEval evaluates the integer polynomial coeffs (low-to-high) at x exactly.
func PolyEval(coeffs []*big.Int, x *big.Int) *big.Int {
	acc := new(big.Int)
	for i := len(coeffs) - 1; i >= 0; i-- {
		acc.Mul(acc, x)
		acc.Add(acc, coeffs[i])
	}
	return acc
}

// PolyDerivative returns the formal derivative of an integer polynomial given
// low-to-high.
func PolyDerivative(coeffs []*big.Int) []*big.Int {
	return polyDerivative(coeffs)
}

// SimpleRootsModP returns the residues r in [0, p) with f(r) == 0 mod p and
// f'(r) != 0 mod p, i.e. the simple roots of f modulo p that Hensel's lemma can
// lift. f is given by integer coefficients low-to-high.
func SimpleRootsModP(f []*big.Int, p *big.Int) []*big.Int {
	df := polyDerivative(f)
	var roots []*big.Int
	r := new(big.Int)
	for r.Cmp(p) < 0 {
		if polyEvalMod(f, r, p).Sign() == 0 && polyEvalMod(df, r, p).Sign() != 0 {
			roots = append(roots, new(big.Int).Set(r))
		}
		r.Add(r, bigOne)
	}
	return roots
}

// RootsModP returns all residues r in [0, p) with f(r) == 0 mod p, including
// multiple roots.
func RootsModP(f []*big.Int, p *big.Int) []*big.Int {
	var roots []*big.Int
	r := new(big.Int)
	for r.Cmp(p) < 0 {
		if polyEvalMod(f, r, p).Sign() == 0 {
			roots = append(roots, new(big.Int).Set(r))
		}
		r.Add(r, bigOne)
	}
	return roots
}

// HenselLift lifts a simple root r0 of the integer polynomial f modulo p to a
// root modulo p^prec, returning the lifted root in [0, p^prec). It requires
// f(r0) == 0 mod p and f'(r0) invertible mod p; otherwise ErrNoRoot is
// returned. The lift is by linear (one-digit-at-a-time) Newton iteration and
// is unique by Hensel's lemma.
func HenselLift(f []*big.Int, r0, p *big.Int, prec int) (*big.Int, error) {
	if prec <= 0 {
		return nil, ErrPrecision
	}
	df := polyDerivative(f)
	r := new(big.Int).Mod(r0, p)
	if polyEvalMod(f, r, p).Sign() != 0 {
		return nil, ErrNoRoot
	}
	dfr := polyEvalMod(df, r, p)
	if dfr.Sign() == 0 {
		return nil, ErrNoRoot
	}
	dfInv := new(big.Int).ModInverse(dfr, p)
	if dfInv == nil {
		return nil, ErrNoRoot
	}
	for k := 1; k < prec; k++ {
		modK1 := PPow(p, k+1)
		fr := polyEvalMod(f, r, modK1) // divisible by p^k
		pk := PPow(p, k)
		frShift := new(big.Int).Div(fr, pk)
		frShift.Mod(frShift, p)
		// t = -(f(r)/p^k) * f'(r)^{-1} mod p
		t := new(big.Int).Mul(frShift, dfInv)
		t.Neg(t)
		t.Mod(t, p)
		r.Add(r, new(big.Int).Mul(t, pk))
		r.Mod(r, modK1)
	}
	return r, nil
}

// HenselLiftPadic lifts a simple root r0 of f modulo p to a p-adic number to
// relative precision prec.
func HenselLiftPadic(f []*big.Int, r0, p *big.Int, prec int) (*Padic, error) {
	lifted, err := HenselLift(f, r0, p, prec)
	if err != nil {
		return nil, err
	}
	return FromBigInt(p, lifted, prec).ReduceTo(prec), nil
}

// PadicRoots returns all p-adic roots (to relative precision prec) obtainable
// from the simple roots of the integer polynomial f modulo p by Hensel lifting.
// Roots that are not simple modulo p are omitted, since Hensel's lemma does not
// apply directly to them.
func PadicRoots(f []*big.Int, p *big.Int, prec int) ([]*Padic, error) {
	if prec <= 0 {
		return nil, ErrPrecision
	}
	simple := SimpleRootsModP(f, p)
	out := make([]*Padic, 0, len(simple))
	for _, r0 := range simple {
		lifted, err := HenselLift(f, r0, p, prec)
		if err != nil {
			continue
		}
		out = append(out, FromBigInt(p, lifted, prec).ReduceTo(prec))
	}
	return out, nil
}

// NewtonStepPadic performs one full-precision Newton step r' = r - f(r)/f'(r)
// on p-adic numbers, given f and f' evaluated via the integer coefficient
// slices. It returns an error if f'(r) is not a unit.
func NewtonStepPadic(f []*big.Int, r *Padic) (*Padic, error) {
	fr, err := evalPadic(f, r)
	if err != nil {
		return nil, err
	}
	dfr, err := evalPadic(polyDerivative(f), r)
	if err != nil {
		return nil, err
	}
	quot, err := fr.Div(dfr)
	if err != nil {
		return nil, err
	}
	return r.Sub(quot)
}

// evalPadic evaluates an integer-coefficient polynomial at a p-adic argument.
func evalPadic(coeffs []*big.Int, x *Padic) (*Padic, error) {
	prec := maxInt(x.AbsolutePrecision(), 1)
	acc := newZero(x.p, prec)
	for i := len(coeffs) - 1; i >= 0; i-- {
		var err error
		acc, err = acc.Mul(x)
		if err != nil {
			return nil, err
		}
		c := FromBigInt(x.p, coeffs[i], prec)
		acc, err = acc.Add(c)
		if err != nil {
			return nil, err
		}
	}
	return acc, nil
}

// EvalPolyPadic evaluates the integer polynomial f (low-to-high) at the p-adic
// number x.
func EvalPolyPadic(f []*big.Int, x *Padic) (*Padic, error) {
	return evalPadic(f, x)
}
