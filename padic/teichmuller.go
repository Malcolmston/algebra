package padic

import "math/big"

// TeichmullerRep returns the Teichmuller representative of the unit a modulo
// p^prec: the unique (p-1)-th root of unity (for p = 2, the unique element of
// {1, ...} congruent to a) that is congruent to a modulo p. It is computed by
// the contracting iteration x -> x^p, which converges p-adically. a must be
// coprime to p; otherwise ErrDomain is returned.
func TeichmullerRep(p, a *big.Int, prec int) (*big.Int, error) {
	if prec <= 0 {
		return nil, ErrPrecision
	}
	mod := PPow(p, prec)
	am := new(big.Int).Mod(a, mod)
	if am.Sign() == 0 || new(big.Int).Mod(am, p).Sign() == 0 {
		return nil, ErrDomain
	}
	x := new(big.Int).Mod(am, p) // lift of a mod p
	next := new(big.Int)
	// Iterate x <- x^p mod p^prec until the value stabilises. Convergence is
	// guaranteed within prec steps.
	for i := 0; i < prec+2; i++ {
		next.Exp(x, p, mod)
		if next.Cmp(x) == 0 {
			break
		}
		x.Set(next)
	}
	return x, nil
}

// Teichmuller returns the Teichmuller representative of a p-adic unit x as a
// new p-adic number of valuation 0 (a root of unity). x must be a unit.
func (x *Padic) Teichmuller() (*Padic, error) {
	if !x.IsUnit() {
		return nil, ErrDomain
	}
	rep, err := TeichmullerRep(x.p, x.unit, x.prec)
	if err != nil {
		return nil, err
	}
	return &Padic{p: new(big.Int).Set(x.p), val: 0, unit: rep, prec: x.prec}, nil
}

// Teichmuller returns the Teichmuller representative of the integer a as a
// p-adic unit to relative precision prec. It is a package-level convenience
// wrapper over TeichmullerRep.
func Teichmuller(p, a *big.Int, prec int) (*Padic, error) {
	rep, err := TeichmullerRep(p, a, prec)
	if err != nil {
		return nil, err
	}
	return &Padic{p: new(big.Int).Set(p), val: 0, unit: rep, prec: prec}, nil
}

// IsTeichmuller reports whether the integer a, taken modulo p^prec, is its own
// Teichmuller representative, i.e. satisfies a^p == a modulo p^prec.
func IsTeichmuller(p, a *big.Int, prec int) bool {
	if prec <= 0 {
		return false
	}
	mod := PPow(p, prec)
	am := new(big.Int).Mod(a, mod)
	if new(big.Int).Mod(am, p).Sign() == 0 {
		return am.Sign() == 0
	}
	powered := new(big.Int).Exp(am, p, mod)
	return powered.Cmp(am) == 0
}

// TeichmullerDigits returns the Teichmuller representatives of the standard
// p-adic digits of x, i.e. for each digit d in [0, p) it returns
// TeichmullerRep(d) (with 0 mapping to 0). The slice has length equal to the
// relative precision of x. x must be a unit or zero.
func (x *Padic) TeichmullerDigits() ([]*big.Int, error) {
	digits, err := x.Digits()
	if err != nil {
		return nil, err
	}
	out := make([]*big.Int, len(digits))
	for i, d := range digits {
		if d.Sign() == 0 {
			out[i] = big.NewInt(0)
			continue
		}
		rep, err := TeichmullerRep(x.p, d, x.prec-i)
		if err != nil {
			out[i] = new(big.Int).Set(d)
			continue
		}
		out[i] = rep
	}
	return out, nil
}

// RootOfUnity returns a primitive (p-1)-th root of unity in Z_p to relative
// precision prec: the Teichmuller representative of a primitive root modulo p.
// rng supplies randomness used only to locate a primitive root modulo p.
func RootOfUnity(p *big.Int, prec int) (*Padic, error) {
	g := PrimitiveRootModP(p)
	if g == nil {
		return nil, ErrDomain
	}
	return Teichmuller(p, g, prec)
}

// PrimitiveRootModP returns the smallest primitive root modulo the odd prime p,
// or nil if p is 2 (where the group is trivial) or has none.
func PrimitiveRootModP(p *big.Int) *big.Int {
	if p.Cmp(bigTwo) <= 0 {
		return big.NewInt(1)
	}
	phi := new(big.Int).Sub(p, bigOne)
	factors := distinctPrimeFactors(phi)
	g := big.NewInt(2)
	for g.Cmp(p) < 0 {
		ok := true
		for _, q := range factors {
			e := new(big.Int).Div(phi, q)
			if new(big.Int).Exp(g, e, p).Cmp(bigOne) == 0 {
				ok = false
				break
			}
		}
		if ok {
			return new(big.Int).Set(g)
		}
		g.Add(g, bigOne)
	}
	return nil
}

// distinctPrimeFactors returns the distinct prime factors of n by trial
// division. It is used only for small exponents p-1.
func distinctPrimeFactors(n *big.Int) []*big.Int {
	var out []*big.Int
	m := new(big.Int).Set(n)
	d := big.NewInt(2)
	r := new(big.Int)
	for new(big.Int).Mul(d, d).Cmp(m) <= 0 {
		if new(big.Int).Mod(m, d).Sign() == 0 {
			out = append(out, new(big.Int).Set(d))
			for {
				q := new(big.Int)
				q.QuoRem(m, d, r)
				if r.Sign() != 0 {
					break
				}
				m.Set(q)
			}
		}
		d.Add(d, bigOne)
	}
	if m.Cmp(bigOne) > 0 {
		out = append(out, m)
	}
	return out
}

// IsRootOfUnityInt reports whether the integer a modulo p^prec is a root of
// unity in Z_p, i.e. is coprime to p and equal to its Teichmuller
// representative.
func IsRootOfUnityInt(p, a *big.Int, prec int) bool {
	if new(big.Int).Mod(a, p).Sign() == 0 {
		return false
	}
	return IsTeichmuller(p, a, prec)
}

// TeichmullerCharacter returns the value omega(a) mod p^prec of the Teichmuller
// character on the integer a, a synonym for TeichmullerRep that reads naturally
// in character-theoretic contexts. a must be coprime to p.
func TeichmullerCharacter(p, a *big.Int, prec int) (*big.Int, error) {
	return TeichmullerRep(p, a, prec)
}
