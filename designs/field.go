package designs

import (
	"errors"
	"math/big"
)

// polyTrim removes leading (highest-degree) zero coefficients.
func polyTrim(a []int) []int {
	i := len(a) - 1
	for i > 0 && a[i] == 0 {
		i--
	}
	return a[:i+1]
}

// polyDeg returns the degree of a polynomial given as coefficients (index =
// degree). The zero polynomial has degree -1.
func polyDeg(a []int) int {
	a = polyTrim(a)
	if len(a) == 1 && a[0] == 0 {
		return -1
	}
	return len(a) - 1
}

// polyMod reduces polynomial a modulo f over GF(p), returning the remainder.
func polyMod(a, f []int, p int) []int {
	r := make([]int, len(a))
	for i := range a {
		r[i] = ((a[i] % p) + p) % p
	}
	r = polyTrim(r)
	df := polyDeg(f)
	lead := f[df]
	invLead, _ := ModInverse(lead, p)
	for polyDeg(r) >= df && polyDeg(r) >= 0 {
		dr := polyDeg(r)
		c := r[dr] * invLead % p
		shift := dr - df
		for i := 0; i <= df; i++ {
			r[i+shift] = ((r[i+shift]-c*f[i])%p + p) % p
		}
		r = polyTrim(r)
	}
	return polyTrim(r)
}

// polyMulMod multiplies a and b modulo f over GF(p).
func polyMulMod(a, b, f []int, p int) []int {
	prod := make([]int, len(a)+len(b)-1)
	for i := range a {
		if a[i] == 0 {
			continue
		}
		for j := range b {
			prod[i+j] = (prod[i+j] + a[i]*b[j]) % p
		}
	}
	return polyMod(prod, f, p)
}

// polyPowMod computes base**e modulo f over GF(p) for a non-negative big
// exponent e.
func polyPowMod(base []int, e *big.Int, f []int, p int) []int {
	result := []int{1 % p}
	b := polyMod(base, f, p)
	ee := new(big.Int).Set(e)
	zero := big.NewInt(0)
	for ee.Cmp(zero) > 0 {
		if ee.Bit(0) == 1 {
			result = polyMulMod(result, b, f, p)
		}
		ee.Rsh(ee, 1)
		if ee.Cmp(zero) > 0 {
			b = polyMulMod(b, b, f, p)
		}
	}
	return result
}

// polyGCD returns a monic greatest common divisor of a and b over GF(p).
func polyGCD(a, b []int, p int) []int {
	aa := append([]int(nil), a...)
	bb := append([]int(nil), b...)
	for i := range aa {
		aa[i] = ((aa[i] % p) + p) % p
	}
	for i := range bb {
		bb[i] = ((bb[i] % p) + p) % p
	}
	aa = polyTrim(aa)
	bb = polyTrim(bb)
	for polyDeg(bb) >= 0 {
		r := polyMod(aa, bb, p)
		aa, bb = bb, r
	}
	// make monic
	d := polyDeg(aa)
	if d < 0 {
		return []int{0}
	}
	inv, _ := ModInverse(aa[d], p)
	for i := range aa {
		aa[i] = aa[i] * inv % p
	}
	return aa
}

// IsIrreducible reports whether the polynomial f (coefficients indexed by
// degree) is irreducible over the prime field GF(p), using Rabin's test. The
// polynomial must have positive degree.
func IsIrreducible(f []int, p int) bool {
	if !IsPrime(p) {
		return false
	}
	n := polyDeg(f)
	if n < 1 {
		return false
	}
	x := []int{0, 1}
	primes, _ := Factorize(n)
	for _, q := range primes {
		m := n / q
		e := new(big.Int).Exp(big.NewInt(int64(p)), big.NewInt(int64(m)), nil)
		h := polyPowMod(x, e, f, p) // x^(p^m) mod f
		// h - x
		hx := append([]int(nil), h...)
		for len(hx) < 2 {
			hx = append(hx, 0)
		}
		hx[1] = ((hx[1]-1)%p + p) % p
		g := polyGCD(hx, f, p)
		if polyDeg(g) != 0 {
			return false
		}
	}
	e := new(big.Int).Exp(big.NewInt(int64(p)), big.NewInt(int64(n)), nil)
	g := polyPowMod(x, e, f, p) // x^(p^n) mod f, must equal x
	g = polyTrim(g)
	if !(polyDeg(g) == 1 && g[0] == 0 && g[1] == 1%p) {
		return false
	}
	return true
}

// IrreduciblePolynomial returns a monic irreducible polynomial of degree k over
// GF(p), as coefficients indexed by degree (length k+1, leading coefficient 1).
// It reports an error when p is not prime or k<1.
func IrreduciblePolynomial(p, k int) ([]int, error) {
	if !IsPrime(p) {
		return nil, errors.New("designs: p must be prime")
	}
	if k < 1 {
		return nil, errors.New("designs: degree must be positive")
	}
	if k == 1 {
		return []int{0, 1}, nil
	}
	// Enumerate monic polynomials of degree k with nonzero constant term.
	coeffs := make([]int, k+1)
	coeffs[k] = 1
	total := 1
	for i := 0; i < k; i++ {
		total *= p
	}
	for idx := 0; idx < total; idx++ {
		n := idx
		for i := 0; i < k; i++ {
			coeffs[i] = n % p
			n /= p
		}
		if coeffs[0] == 0 {
			continue
		}
		f := append([]int(nil), coeffs...)
		if IsIrreducible(f, p) {
			return f, nil
		}
	}
	return nil, errors.New("designs: no irreducible polynomial found")
}

// GaloisField is a finite field GF(q) with q = p**k elements. Elements are
// represented as integers in [0,q): the base-p digits of an element are the
// coefficients (constant term first) of a polynomial of degree < k in a root of
// the defining irreducible polynomial. For a prime field (k=1) this is ordinary
// arithmetic modulo p.
type GaloisField struct {
	p    int   // characteristic (prime)
	k    int   // extension degree
	q    int   // order p**k
	poly []int // defining monic irreducible polynomial, degree k
}

// NewGaloisField constructs GF(q) for a prime power q = p**k. It reports an
// error when q is not a prime power.
func NewGaloisField(q int) (*GaloisField, error) {
	p, k, err := FactorPrimePower(q)
	if err != nil {
		return nil, err
	}
	f := &GaloisField{p: p, k: k, q: q}
	if k == 1 {
		f.poly = []int{0, 1}
		return f, nil
	}
	poly, err := IrreduciblePolynomial(p, k)
	if err != nil {
		return nil, err
	}
	f.poly = poly
	return f, nil
}

// Order returns the number of elements q = p**k of the field.
func (f *GaloisField) Order() int { return f.q }

// Characteristic returns the prime characteristic p of the field.
func (f *GaloisField) Characteristic() int { return f.p }

// Degree returns the extension degree k, so that q = p**k.
func (f *GaloisField) Degree() int { return f.k }

// Poly returns a copy of the defining monic irreducible polynomial (nil for a
// prime field's trivial polynomial is instead returned as [0,1]).
func (f *GaloisField) Poly() []int { return append([]int(nil), f.poly...) }

// ToCoefficients returns the polynomial coefficients (constant term first,
// length k) of the field element e.
func (f *GaloisField) ToCoefficients(e int) []int {
	c := make([]int, f.k)
	for i := 0; i < f.k; i++ {
		c[i] = e % f.p
		e /= f.p
	}
	return c
}

// FromCoefficients returns the field element whose polynomial coefficients
// (constant term first) are the given values, reduced modulo p.
func (f *GaloisField) FromCoefficients(c []int) int {
	e := 0
	base := 1
	for i := 0; i < f.k; i++ {
		var v int
		if i < len(c) {
			v = ((c[i] % f.p) + f.p) % f.p
		}
		e += v * base
		base *= f.p
	}
	return e
}

// Zero returns the additive identity of the field.
func (f *GaloisField) Zero() int { return 0 }

// One returns the multiplicative identity of the field.
func (f *GaloisField) One() int { return 1 }

// Elements returns all q elements of the field, 0,1,...,q-1.
func (f *GaloisField) Elements() []int {
	out := make([]int, f.q)
	for i := range out {
		out[i] = i
	}
	return out
}

// Add returns the field sum a+b.
func (f *GaloisField) Add(a, b int) int {
	ca := f.ToCoefficients(a)
	cb := f.ToCoefficients(b)
	for i := range ca {
		ca[i] = (ca[i] + cb[i]) % f.p
	}
	return f.FromCoefficients(ca)
}

// Neg returns the field additive inverse -a.
func (f *GaloisField) Neg(a int) int {
	ca := f.ToCoefficients(a)
	for i := range ca {
		ca[i] = (f.p - ca[i]) % f.p
	}
	return f.FromCoefficients(ca)
}

// Sub returns the field difference a-b.
func (f *GaloisField) Sub(a, b int) int { return f.Add(a, f.Neg(b)) }

// Mul returns the field product a*b.
func (f *GaloisField) Mul(a, b int) int {
	if f.k == 1 {
		return a * b % f.p
	}
	ca := f.ToCoefficients(a)
	cb := f.ToCoefficients(b)
	prod := polyMulMod(ca, cb, f.poly, f.p)
	return f.FromCoefficients(prod)
}

// Pow returns a**n for n>=0 in the field.
func (f *GaloisField) Pow(a, n int) int {
	result := 1
	base := a
	for n > 0 {
		if n&1 == 1 {
			result = f.Mul(result, base)
		}
		n >>= 1
		if n > 0 {
			base = f.Mul(base, base)
		}
	}
	return result
}

// Inv returns the multiplicative inverse of a non-zero element a. It reports an
// error when a is zero.
func (f *GaloisField) Inv(a int) (int, error) {
	if a == 0 {
		return 0, errors.New("designs: zero has no inverse")
	}
	// a^(q-2) is the inverse in a field of order q.
	return f.Pow(a, f.q-2), nil
}

// Div returns a/b for non-zero b. It reports an error when b is zero.
func (f *GaloisField) Div(a, b int) (int, error) {
	bi, err := f.Inv(b)
	if err != nil {
		return 0, err
	}
	return f.Mul(a, bi), nil
}

// ElementOrder returns the multiplicative order of a non-zero element a, the
// least e>=1 with a**e == 1. It reports an error when a is zero.
func (f *GaloisField) ElementOrder(a int) (int, error) {
	if a == 0 {
		return 0, errors.New("designs: zero has no multiplicative order")
	}
	e := 1
	cur := a
	for cur != 1 {
		cur = f.Mul(cur, a)
		e++
	}
	return e, nil
}

// IsPrimitiveElement reports whether the non-zero element a generates the
// multiplicative group, i.e. has order q-1.
func (f *GaloisField) IsPrimitiveElement(a int) bool {
	if a == 0 {
		return false
	}
	// order divides q-1; check a^((q-1)/r) != 1 for each prime r | q-1.
	phi := f.q - 1
	primes, _ := Factorize(phi)
	for _, r := range primes {
		if f.Pow(a, phi/r) == 1 {
			return false
		}
	}
	return true
}

// PrimitiveElement returns a generator of the multiplicative group of the
// field. It reports an error when none is found (which cannot happen for a
// valid field).
func (f *GaloisField) PrimitiveElement() (int, error) {
	for a := 1; a < f.q; a++ {
		if f.IsPrimitiveElement(a) {
			return a, nil
		}
	}
	return 0, errors.New("designs: no primitive element found")
}

// Trace returns the field trace of a down to the prime field, defined as
// a + a**p + a**(p^2) + ... + a**(p^(k-1)), a value in {0,...,p-1}.
func (f *GaloisField) Trace(a int) int {
	sum := 0
	cur := a
	for i := 0; i < f.k; i++ {
		sum = f.Add(sum, cur)
		cur = f.Pow(cur, f.p)
	}
	// result lies in the prime field: its constant coefficient is the value.
	return f.ToCoefficients(sum)[0]
}

// Norm returns the field norm of a down to the prime field, defined as the
// product a * a**p * ... * a**(p^(k-1)) = a**((q-1)/(p-1)), a value in
// {0,...,p-1}.
func (f *GaloisField) Norm(a int) int {
	prod := 1
	cur := a
	for i := 0; i < f.k; i++ {
		prod = f.Mul(prod, cur)
		cur = f.Pow(cur, f.p)
	}
	return f.ToCoefficients(prod)[0]
}
