package padic

import (
	"fmt"
	"math"
	"math/big"
)

// Padic is a p-adic number carried to a finite, tracked precision.
//
// A non-zero element equals p^val * unit, where unit is in [0, p^prec) and
// coprime to p. The absolute precision (val + prec) is the exponent modulo
// which the element is known. The zero element carries prec == 0 and a val
// equal to its absolute precision.
type Padic struct {
	p    *big.Int // the prime
	val  int      // valuation (or absolute precision when zero)
	unit *big.Int // unit part in [0, p^prec); zero for the zero element
	prec int      // relative precision (digits of unit); 0 for zero
}

// newZero builds the p-adic zero known to the given absolute precision.
func newZero(p *big.Int, absPrec int) *Padic {
	return &Padic{p: new(big.Int).Set(p), val: absPrec, unit: big.NewInt(0), prec: 0}
}

// makeScaled builds the p-adic number (m * p^base) known to the given absolute
// precision. It extracts the p-power from m and normalises the unit.
func makeScaled(p *big.Int, base int, m *big.Int, absPrec int) *Padic {
	rp := absPrec - base
	if rp <= 0 {
		return newZero(p, absPrec)
	}
	mod := PPow(p, rp)
	mm := new(big.Int).Mod(m, mod)
	if mm.Sign() == 0 {
		return newZero(p, absPrec)
	}
	vs := 0
	q := new(big.Int)
	r := new(big.Int)
	for {
		q.QuoRem(mm, p, r)
		if r.Sign() != 0 {
			break
		}
		mm.Set(q)
		vs++
	}
	newPrec := rp - vs
	unit := new(big.Int).Mod(mm, PPow(p, newPrec))
	return &Padic{p: new(big.Int).Set(p), val: base + vs, unit: unit, prec: newPrec}
}

// New constructs the p-adic number p^val * unit known to relative precision
// prec (so absolute precision val+prec). The unit is reduced and its p-part is
// re-extracted, so callers need not pass a genuine unit. prec must be
// positive; a zero unit yields the zero element at absolute precision val+prec.
func New(p *big.Int, val int, unit *big.Int, prec int) (*Padic, error) {
	if !IsPrime(p) {
		return nil, ErrNotPrime
	}
	if prec <= 0 {
		return nil, ErrPrecision
	}
	if unit.Sign() == 0 {
		return newZero(p, val+prec), nil
	}
	return makeScaled(p, val, unit, val+prec), nil
}

// Zero returns the p-adic zero known to the given absolute precision.
func Zero(p *big.Int, absPrec int) *Padic {
	return newZero(p, absPrec)
}

// One returns the p-adic number 1 to relative precision prec.
func One(p *big.Int, prec int) *Padic {
	if prec <= 0 {
		prec = 1
	}
	return &Padic{p: new(big.Int).Set(p), val: 0, unit: big.NewInt(1), prec: prec}
}

// FromInt64 returns the p-adic value of the machine integer n to relative
// precision prec.
func FromInt64(p *big.Int, n int64, prec int) *Padic {
	return FromBigInt(p, big.NewInt(n), prec)
}

// FromInt is an alias for FromInt64 for callers that pass an int.
func FromInt(p *big.Int, n int, prec int) *Padic {
	return FromBigInt(p, big.NewInt(int64(n)), prec)
}

// FromBigInt returns the p-adic value of the integer n to relative precision
// prec. For n = 0 it returns the zero element at absolute precision prec.
func FromBigInt(p *big.Int, n *big.Int, prec int) *Padic {
	if prec <= 0 {
		prec = 1
	}
	if n.Sign() == 0 {
		return newZero(p, prec)
	}
	return makeScaled(p, 0, n, ValuationInt(p, n)+prec)
}

// FromRational returns the p-adic value of the rational a/b to relative
// precision prec. b must be non-zero.
func FromRational(p *big.Int, a, b *big.Int, prec int) (*Padic, error) {
	if b.Sign() == 0 {
		return nil, ErrZeroDivision
	}
	if prec <= 0 {
		return nil, ErrPrecision
	}
	if a.Sign() == 0 {
		return newZero(p, prec), nil
	}
	val := ValuationRat(p, a, b)
	au := UnitPartInt(p, a)
	bu := UnitPartInt(p, b)
	mod := PPow(p, prec)
	binv := new(big.Int).ModInverse(new(big.Int).Mod(bu, mod), mod)
	if binv == nil {
		return nil, ErrNotInvertible
	}
	unit := new(big.Int).Mul(new(big.Int).Mod(au, mod), binv)
	unit.Mod(unit, mod)
	return &Padic{p: new(big.Int).Set(p), val: val, unit: unit, prec: prec}, nil
}

// FromRat returns the p-adic value of the big.Rat r to relative precision prec.
func FromRat(p *big.Int, r *big.Rat, prec int) (*Padic, error) {
	return FromRational(p, r.Num(), r.Denom(), prec)
}

// Uniformizer returns the p-adic number p (the standard uniformizer) to
// relative precision prec.
func Uniformizer(p *big.Int, prec int) *Padic {
	if prec <= 0 {
		prec = 1
	}
	return &Padic{p: new(big.Int).Set(p), val: 1, unit: big.NewInt(1), prec: prec}
}

// Prime returns a copy of the prime p of x.
func (x *Padic) Prime() *big.Int { return new(big.Int).Set(x.p) }

// Valuation returns the p-adic valuation of x. For the zero element it returns
// the absolute precision, an honest lower bound on the (infinite) valuation;
// use IsZero to distinguish.
func (x *Padic) Valuation() int { return x.val }

// RelativePrecision returns the number of known p-adic digits of the unit part
// of x. It is 0 for the zero element.
func (x *Padic) RelativePrecision() int { return x.prec }

// AbsolutePrecision returns the exponent v such that x is known modulo p^v.
func (x *Padic) AbsolutePrecision() int { return x.val + x.prec }

// Unit returns the unit part of x (an integer coprime to p) reduced modulo
// p^prec, or 0 for the zero element.
func (x *Padic) Unit() *big.Int { return new(big.Int).Set(x.unit) }

// IsZero reports whether x is the p-adic zero to its tracked precision.
func (x *Padic) IsZero() bool { return x.prec == 0 }

// IsUnit reports whether x is a p-adic unit, i.e. has valuation 0.
func (x *Padic) IsUnit() bool { return !x.IsZero() && x.val == 0 }

// IsIntegral reports whether x lies in the ring of p-adic integers, i.e. has
// non-negative valuation. The zero element is integral.
func (x *Padic) IsIntegral() bool { return x.IsZero() || x.val >= 0 }

// IsOne reports whether x equals 1 to its tracked precision.
func (x *Padic) IsOne() bool {
	return x.IsUnit() && x.unit.Cmp(bigOne) == 0
}

// Copy returns an independent deep copy of x.
func (x *Padic) Copy() *Padic {
	return &Padic{p: new(big.Int).Set(x.p), val: x.val, unit: new(big.Int).Set(x.unit), prec: x.prec}
}

// samePrime reports whether x and y share the same prime.
func (x *Padic) samePrime(y *Padic) bool { return x.p.Cmp(y.p) == 0 }

// ReduceTo returns x truncated to the given absolute precision. If absPrec is
// at least the current absolute precision, an exact copy is returned; a p-adic
// number's precision can never be increased without more information.
func (x *Padic) ReduceTo(absPrec int) *Padic {
	if absPrec >= x.AbsolutePrecision() {
		return x.Copy()
	}
	if x.IsZero() {
		return newZero(x.p, minInt(x.val, absPrec))
	}
	return makeScaled(x.p, x.val, x.unit, absPrec)
}

// Add returns x + y. The result is known to the smaller of the two absolute
// precisions. It returns an error only when the primes differ.
func (x *Padic) Add(y *Padic) (*Padic, error) {
	if !x.samePrime(y) {
		return nil, ErrPrimeMismatch
	}
	n := minInt(x.AbsolutePrecision(), y.AbsolutePrecision())
	if x.IsZero() {
		return y.ReduceTo(n), nil
	}
	if y.IsZero() {
		return x.ReduceTo(n), nil
	}
	base := minInt(x.val, y.val)
	xs := new(big.Int).Mul(x.unit, PPow(x.p, x.val-base))
	ys := new(big.Int).Mul(y.unit, PPow(x.p, y.val-base))
	m := xs.Add(xs, ys)
	return makeScaled(x.p, base, m, n), nil
}

// Sub returns x - y to the smaller of the two absolute precisions.
func (x *Padic) Sub(y *Padic) (*Padic, error) {
	return x.Add(y.Neg())
}

// Neg returns -x with the same precision as x.
func (x *Padic) Neg() *Padic {
	if x.IsZero() {
		return x.Copy()
	}
	mod := PPow(x.p, x.prec)
	u := new(big.Int).Sub(mod, x.unit)
	u.Mod(u, mod)
	return &Padic{p: new(big.Int).Set(x.p), val: x.val, unit: u, prec: x.prec}
}

// Mul returns x * y. The relative precision of a product is the smaller of the
// inputs' relative precisions.
func (x *Padic) Mul(y *Padic) (*Padic, error) {
	if !x.samePrime(y) {
		return nil, ErrPrimeMismatch
	}
	if x.IsZero() || y.IsZero() {
		return newZero(x.p, x.val+y.val), nil
	}
	rp := minInt(x.prec, y.prec)
	mod := PPow(x.p, rp)
	u := new(big.Int).Mul(x.unit, y.unit)
	u.Mod(u, mod)
	return &Padic{p: new(big.Int).Set(x.p), val: x.val + y.val, unit: u, prec: rp}, nil
}

// Square returns x*x.
func (x *Padic) Square() *Padic {
	r, _ := x.Mul(x)
	return r
}

// Inv returns the multiplicative inverse of x, or an error if x is zero.
func (x *Padic) Inv() (*Padic, error) {
	if x.IsZero() {
		return nil, ErrZeroDivision
	}
	mod := PPow(x.p, x.prec)
	u := new(big.Int).ModInverse(x.unit, mod)
	if u == nil {
		return nil, ErrNotInvertible
	}
	return &Padic{p: new(big.Int).Set(x.p), val: -x.val, unit: u, prec: x.prec}, nil
}

// Div returns x / y, or an error if y is zero or the primes differ.
func (x *Padic) Div(y *Padic) (*Padic, error) {
	if !x.samePrime(y) {
		return nil, ErrPrimeMismatch
	}
	yi, err := y.Inv()
	if err != nil {
		return nil, err
	}
	return x.Mul(yi)
}

// Pow returns x raised to the integer power n, using binary exponentiation.
// Negative n inverts x first and requires x to be non-zero.
func (x *Padic) Pow(n int) (*Padic, error) {
	if n == 0 {
		return One(x.p, maxInt(x.prec, 1)), nil
	}
	base := x
	if n < 0 {
		inv, err := x.Inv()
		if err != nil {
			return nil, err
		}
		base = inv
		n = -n
	}
	result := One(base.p, maxInt(base.prec, 1))
	for n > 0 {
		if n&1 == 1 {
			var err error
			result, err = result.Mul(base)
			if err != nil {
				return nil, err
			}
		}
		n >>= 1
		if n > 0 {
			base = base.Square()
		}
	}
	return result, nil
}

// Equal reports whether x and y agree as p-adic numbers to the smaller of
// their two absolute precisions.
func (x *Padic) Equal(y *Padic) bool {
	if !x.samePrime(y) {
		return false
	}
	d, err := x.Sub(y)
	if err != nil {
		return false
	}
	return d.IsZero()
}

// IsCloseTo reports whether x and y agree modulo p^n, i.e. the valuation of
// their difference is at least n.
func (x *Padic) IsCloseTo(y *Padic, n int) bool {
	d, err := x.Sub(y)
	if err != nil {
		return false
	}
	return d.IsZero() || d.val >= n
}

// AbsValue returns the p-adic absolute value |x|_p = p^(-val) as an exact
// big.Rat. For the zero element it returns 0.
func (x *Padic) AbsValue() *big.Rat {
	if x.IsZero() {
		return new(big.Rat)
	}
	return absFromVal(x.p, x.val)
}

// NormFloat returns the p-adic absolute value of x as a float64.
func (x *Padic) NormFloat() float64 {
	if x.IsZero() {
		return 0
	}
	pf, _ := new(big.Float).SetInt(x.p).Float64()
	return math.Pow(pf, float64(-x.val))
}

// Rat returns an exact rational representative of x, namely unit * p^val with
// unit in [0, p^prec). For the zero element it returns 0. Distinct p-adic
// numbers within the same precision may share a representative only if equal.
func (x *Padic) Rat() *big.Rat {
	if x.IsZero() {
		return new(big.Rat)
	}
	r := new(big.Rat).SetInt(x.unit)
	if x.val >= 0 {
		r.Mul(r, new(big.Rat).SetInt(PPow(x.p, x.val)))
	} else {
		r.Quo(r, new(big.Rat).SetInt(PPow(x.p, -x.val)))
	}
	return r
}

// BigInt returns the canonical integer representative of x in [0, p^absprec)
// when x is integral (valuation >= 0), or an error otherwise.
func (x *Padic) BigInt() (*big.Int, error) {
	if x.IsZero() {
		return big.NewInt(0), nil
	}
	if x.val < 0 {
		return nil, ErrDomain
	}
	return new(big.Int).Mul(x.unit, PPow(x.p, x.val)), nil
}

// String renders x as "unit*p^val + O(p^absprec)" or "O(p^absprec)" for zero.
func (x *Padic) String() string {
	if x.IsZero() {
		return fmt.Sprintf("O(%s^%d)", x.p.String(), x.val)
	}
	return fmt.Sprintf("%s*%s^%d + O(%s^%d)", x.unit.String(), x.p.String(), x.val, x.p.String(), x.AbsolutePrecision())
}

// Format implements fmt.Formatter for the %v, %s and %d verbs by delegating to
// String, so Padic values print cleanly.
func (x *Padic) Format(f fmt.State, verb rune) {
	fmt.Fprint(f, x.String())
}

// MulPow returns x multiplied by p^k, i.e. shifts the valuation by k while
// preserving the unit and relative precision. Negative k is allowed and
// divides by p^|k|.
func (x *Padic) MulPow(k int) *Padic {
	if x.IsZero() {
		return newZero(x.p, x.val+k)
	}
	return &Padic{p: new(big.Int).Set(x.p), val: x.val + k, unit: new(big.Int).Set(x.unit), prec: x.prec}
}

// ShiftUp returns x * p^k for k >= 0 (a synonym for MulPow with non-negative k).
func (x *Padic) ShiftUp(k int) *Padic { return x.MulPow(k) }

// ShiftDown returns x / p^k for k >= 0 (a synonym for MulPow with negative k).
func (x *Padic) ShiftDown(k int) *Padic { return x.MulPow(-k) }

// AddInt returns x + n for an integer n, tracked at x's precision.
func (x *Padic) AddInt(n int64) *Padic {
	np := FromInt64(x.p, n, maxInt(x.AbsolutePrecision(), 1))
	r, _ := x.Add(np)
	return r
}

// SubInt returns x - n for an integer n.
func (x *Padic) SubInt(n int64) *Padic { return x.AddInt(-n) }

// MulInt returns x * n for an integer n.
func (x *Padic) MulInt(n int64) *Padic {
	np := FromInt64(x.p, n, maxInt(x.AbsolutePrecision(), 1))
	r, _ := x.Mul(np)
	return r
}

// Distance returns the p-adic distance |x - y|_p between x and y as an exact
// big.Rat. Identical values (to precision) have distance 0.
func (x *Padic) Distance(y *Padic) *big.Rat {
	d, err := x.Sub(y)
	if err != nil {
		return nil
	}
	return d.AbsValue()
}

// Digit returns the coefficient of p^i in the standard p-adic expansion of x,
// an integer in [0, p). Positions below the valuation are zero. It returns an
// error if i is at or beyond the absolute precision (the digit is unknown) or x
// has negative valuation below i in a way that is not integral there.
func (x *Padic) Digit(i int) (*big.Int, error) {
	if i < 0 {
		return nil, ErrDomain
	}
	if i >= x.AbsolutePrecision() {
		return nil, ErrPrecision
	}
	if x.IsZero() || i < x.val {
		return big.NewInt(0), nil
	}
	digits := digitsOfUnit(x.p, x.unit, bigOne, x.prec)
	return digits[i-x.val], nil
}

// ConstantResidue returns x reduced modulo p as an integer in [0, p), which
// requires x to be integral (valuation >= 0). Elements of positive valuation
// reduce to 0.
func (x *Padic) ConstantResidue() (*big.Int, error) {
	if x.IsZero() {
		return big.NewInt(0), nil
	}
	if x.val < 0 {
		return nil, ErrDomain
	}
	if x.val > 0 {
		return big.NewInt(0), nil
	}
	return new(big.Int).Mod(x.unit, x.p), nil
}

// IsRootOfUnity reports whether x is a root of unity in Z_p, i.e. a unit equal
// to its own Teichmuller representative (a (p-1)-th root of unity).
func (x *Padic) IsRootOfUnity() bool {
	return x.IsUnit() && IsTeichmuller(x.p, x.unit, x.prec)
}

// Precision returns the absolute precision of x, a synonym for
// AbsolutePrecision that reads naturally at call sites.
func (x *Padic) Precision() int { return x.AbsolutePrecision() }

// FromString parses a decimal integer or rational ("a" or "a/b") and returns
// its p-adic value to relative precision prec. It returns an error for
// malformed input or a zero denominator.
func FromString(p *big.Int, s string, prec int) (*Padic, error) {
	r := new(big.Rat)
	if _, ok := r.SetString(s); !ok {
		return nil, ErrDomain
	}
	return FromRat(p, r, prec)
}

// NewUnit constructs a p-adic unit (valuation 0) from an integer unit value to
// relative precision prec. The value must be coprime to p.
func NewUnit(p, unit *big.Int, prec int) (*Padic, error) {
	if prec <= 0 {
		return nil, ErrPrecision
	}
	if new(big.Int).Mod(unit, p).Sign() == 0 {
		return nil, ErrDomain
	}
	return &Padic{p: new(big.Int).Set(p), val: 0, unit: new(big.Int).Mod(unit, PPow(p, prec)), prec: prec}, nil
}
