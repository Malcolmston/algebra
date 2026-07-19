package contfrac

import (
	"fmt"
	"math"
	"math/big"
	"strconv"
	"strings"
)

// CF is a finite continued fraction [a0; a1, a2, ...] represented by its slice
// of partial quotients. The floor convention is used: a0 may be zero or
// negative, while every later term is a positive integer for a well-formed
// expansion. An empty CF denotes the value 0.
type CF []int64

// FromRational returns the canonical continued-fraction expansion of the
// rational p/q using the Euclidean algorithm. It panics if q == 0. Every
// partial quotient after the first is positive, and the final term is greater
// than one unless the expansion has a single term.
func FromRational(p, q int64) CF {
	if q == 0 {
		panic("contfrac: FromRational requires q != 0")
	}
	if q < 0 {
		p, q = -p, -q
	}
	var cf CF
	for q != 0 {
		a := floorDiv(p, q)
		cf = append(cf, a)
		p, q = q, p-a*q
	}
	return cf
}

// FromFrac returns the continued-fraction expansion of f.
func FromFrac(f Frac) CF {
	return FromRational(f.Num, f.Den)
}

// FromRationalBig returns the continued-fraction expansion of a *big.Rat as a
// CF of int64 terms. The terms of a continued fraction are typically far
// smaller than the numerator and denominator, so int64 terms suffice for
// realistic inputs; callers needing unbounded terms should expand manually.
func FromRationalBig(r *big.Rat) CF {
	p := new(big.Int).Set(r.Num())
	q := new(big.Int).Set(r.Denom()) // always positive for a big.Rat
	var cf CF
	for q.Sign() != 0 {
		a := new(big.Int)
		rem := new(big.Int)
		a.DivMod(p, q, rem) // Euclidean: 0 <= rem < q since q > 0
		cf = append(cf, a.Int64())
		p, q = q, rem
	}
	return cf
}

// FromFloat returns up to maxTerms partial quotients of the continued fraction
// of x, stopping early once the fractional part is negligible. It uses the
// default tolerance 1e-12; see [FromFloatEps] for control over the stopping
// criterion.
func FromFloat(x float64, maxTerms int) CF {
	return FromFloatEps(x, 1e-12, maxTerms)
}

// FromFloatEps returns up to maxTerms partial quotients of the continued
// fraction of x. Expansion stops once the remaining fractional part is at most
// eps (indicating x is, to within tolerance, exactly represented) or once
// maxTerms terms have been produced.
func FromFloatEps(x float64, eps float64, maxTerms int) CF {
	if maxTerms <= 0 {
		maxTerms = 1
	}
	cf := make(CF, 0, maxTerms)
	for i := 0; i < maxTerms; i++ {
		a := int64(math.Floor(x))
		cf = append(cf, a)
		frac := x - float64(a)
		if frac <= eps {
			break
		}
		x = 1 / frac
		if math.IsInf(x, 0) || math.IsNaN(x) {
			break
		}
	}
	return cf
}

// FromBigFloat returns up to n partial quotients of the continued fraction of a
// *big.Float, using the float's own precision for each reciprocal step. It is
// the high-precision counterpart of [FromFloat] and underlies [PiCF].
func FromBigFloat(x *big.Float, n int) CF {
	if n <= 0 {
		return CF{}
	}
	prec := x.Prec()
	if prec == 0 {
		prec = 64
	}
	f := new(big.Float).SetPrec(prec).Set(x)
	one := new(big.Float).SetPrec(prec).SetInt64(1)
	cf := make(CF, 0, n)
	for i := 0; i < n; i++ {
		a := floorBigFloat(f)
		cf = append(cf, a)
		af := new(big.Float).SetPrec(prec).SetInt64(a)
		frac := new(big.Float).SetPrec(prec).Sub(f, af)
		if frac.Sign() <= 0 {
			break
		}
		f = new(big.Float).SetPrec(prec).Quo(one, frac)
	}
	return cf
}

// floorBigFloat returns floor(f) as an int64 for f >= 0 (the only case that
// arises during continued-fraction extraction of a positive irrational).
func floorBigFloat(f *big.Float) int64 {
	i, _ := f.Int(nil)
	return i.Int64()
}

// Len returns the number of partial quotients.
func (c CF) Len() int {
	return len(c)
}

// Terms returns a copy of the underlying slice of partial quotients.
func (c CF) Terms() []int64 {
	out := make([]int64, len(c))
	copy(out, c)
	return out
}

// Clone returns an independent copy of c.
func (c CF) Clone() CF {
	out := make(CF, len(c))
	copy(out, c)
	return out
}

// Value evaluates the continued fraction as a float64 using stable backward
// recurrence, which avoids the overflow that forward convergent recurrence can
// suffer for long expansions.
func (c CF) Value() float64 {
	if len(c) == 0 {
		return 0
	}
	v := float64(c[len(c)-1])
	for i := len(c) - 2; i >= 0; i-- {
		v = float64(c[i]) + 1/v
	}
	return v
}

// Frac evaluates the continued fraction exactly to a reduced [Frac] via the
// convergent recurrence. Overflow is possible for very long expansions; use
// [CF.Rat] for arbitrary precision.
func (c CF) Frac() Frac {
	if len(c) == 0 {
		return Frac{0, 1}
	}
	hPrev, hPrev2 := int64(1), int64(0)
	kPrev, kPrev2 := int64(0), int64(1)
	for _, a := range c {
		h := a*hPrev + hPrev2
		k := a*kPrev + kPrev2
		hPrev, hPrev2 = h, hPrev
		kPrev, kPrev2 = k, kPrev
	}
	return NewFrac(hPrev, kPrev)
}

// Rat evaluates the continued fraction exactly to a *big.Rat, correct for
// expansions of any length.
func (c CF) Rat() *big.Rat {
	hPrev, hPrev2 := big.NewInt(1), big.NewInt(0)
	kPrev, kPrev2 := big.NewInt(0), big.NewInt(1)
	tmp := new(big.Int)
	for _, a := range c {
		ab := big.NewInt(a)
		h := new(big.Int).Add(tmp.Mul(ab, hPrev), hPrev2)
		k := new(big.Int).Add(new(big.Int).Mul(ab, kPrev), kPrev2)
		hPrev, hPrev2 = h, hPrev
		kPrev, kPrev2 = k, kPrev
	}
	if kPrev.Sign() == 0 {
		return new(big.Rat)
	}
	return new(big.Rat).SetFrac(hPrev, kPrev)
}

// Convergent returns the k-th convergent (0-indexed) as a reduced [Frac]. It
// panics if k is out of range.
func (c CF) Convergent(k int) Frac {
	if k < 0 || k >= len(c) {
		panic("contfrac: convergent index out of range")
	}
	hPrev, hPrev2 := int64(1), int64(0)
	kPrev, kPrev2 := int64(0), int64(1)
	for i := 0; i <= k; i++ {
		a := c[i]
		h := a*hPrev + hPrev2
		kk := a*kPrev + kPrev2
		hPrev, hPrev2 = h, hPrev
		kPrev, kPrev2 = kk, kPrev
	}
	return Frac{hPrev, kPrev}
}

// Convergents returns all convergents p0/q0, p1/q1, ... of the continued
// fraction as reduced fractions.
func (c CF) Convergents() []Frac {
	out := make([]Frac, 0, len(c))
	hPrev, hPrev2 := int64(1), int64(0)
	kPrev, kPrev2 := int64(0), int64(1)
	for _, a := range c {
		h := a*hPrev + hPrev2
		k := a*kPrev + kPrev2
		out = append(out, Frac{h, k})
		hPrev, hPrev2 = h, hPrev
		kPrev, kPrev2 = k, kPrev
	}
	return out
}

// ConvergentsBig returns all convergents as *big.Rat, correct for expansions of
// any length.
func (c CF) ConvergentsBig() []*big.Rat {
	out := make([]*big.Rat, 0, len(c))
	hPrev, hPrev2 := big.NewInt(1), big.NewInt(0)
	kPrev, kPrev2 := big.NewInt(0), big.NewInt(1)
	for _, a := range c {
		ab := big.NewInt(a)
		h := new(big.Int).Add(new(big.Int).Mul(ab, hPrev), hPrev2)
		k := new(big.Int).Add(new(big.Int).Mul(ab, kPrev), kPrev2)
		out = append(out, new(big.Rat).SetFrac(new(big.Int).Set(h), new(big.Int).Set(k)))
		hPrev, hPrev2 = h, hPrev
		kPrev, kPrev2 = k, kPrev
	}
	return out
}

// Numerators returns the numerators of the successive convergents.
func (c CF) Numerators() []int64 {
	conv := c.Convergents()
	out := make([]int64, len(conv))
	for i, f := range conv {
		out[i] = f.Num
	}
	return out
}

// Denominators returns the denominators of the successive convergents.
func (c CF) Denominators() []int64 {
	conv := c.Convergents()
	out := make([]int64, len(conv))
	for i, f := range conv {
		out[i] = f.Den
	}
	return out
}

// Semiconvergents returns the sequence of intermediate fractions generated
// along the way to the value: for each partial quotient a_i and each
// m = 1, ..., a_i the fraction (m*h_{i-1}+h_{i-2})/(m*k_{i-1}+k_{i-2}). These
// are the best rational approximations of the first kind and trace the
// Stern-Brocot path to the value.
func (c CF) Semiconvergents() []Frac {
	hPrev, hPrev2 := int64(1), int64(0)
	kPrev, kPrev2 := int64(0), int64(1)
	var out []Frac
	for _, a := range c {
		for m := int64(1); m <= a; m++ {
			out = append(out, Frac{m*hPrev + hPrev2, m*kPrev + kPrev2})
		}
		h := a*hPrev + hPrev2
		k := a*kPrev + kPrev2
		hPrev, hPrev2 = h, hPrev
		kPrev, kPrev2 = k, kPrev
	}
	return out
}

// Truncate returns the continued fraction consisting of the first n partial
// quotients. If n exceeds the length the whole CF is returned.
func (c CF) Truncate(n int) CF {
	if n < 0 {
		n = 0
	}
	if n > len(c) {
		n = len(c)
	}
	return c[:n].Clone()
}

// Reverse returns the continued fraction with its partial quotients in reverse
// order. Reversal leaves the numerator of the final convergent unchanged and is
// a classical symmetry of continuants.
func (c CF) Reverse() CF {
	out := make(CF, len(c))
	for i, v := range c {
		out[len(c)-1-i] = v
	}
	return out
}

// Canonical returns the canonical form of the continued fraction: a trailing
// term equal to one is absorbed into the preceding term ([..., x, 1] becomes
// [..., x+1]), so that finite expansions have a unique representative (except
// for the single-term expansion [1]).
func (c CF) Canonical() CF {
	if len(c) <= 1 {
		return c.Clone()
	}
	out := c.Clone()
	if out[len(out)-1] == 1 {
		out = out[:len(out)-1]
		out[len(out)-1]++
	}
	return out
}

// IsValid reports whether c is a well-formed expansion: it is non-empty and
// every partial quotient after the first is a positive integer.
func (c CF) IsValid() bool {
	if len(c) == 0 {
		return false
	}
	for i := 1; i < len(c); i++ {
		if c[i] <= 0 {
			return false
		}
	}
	return true
}

// Equal reports whether c and d have identical partial quotients.
func (c CF) Equal(d CF) bool {
	if len(c) != len(d) {
		return false
	}
	for i := range c {
		if c[i] != d[i] {
			return false
		}
	}
	return true
}

// EqualValue reports whether c and d denote the same rational number, comparing
// their canonical forms.
func (c CF) EqualValue(d CF) bool {
	return c.Canonical().Equal(d.Canonical())
}

// Compare compares the values of c and d, returning -1, 0 or +1.
func (c CF) Compare(d CF) int {
	return c.Frac().Cmp(d.Frac())
}

// String renders the continued fraction in the usual bracket notation
// "[a0; a1, a2, ...]" (or "[a0]" for a single term).
func (c CF) String() string {
	if len(c) == 0 {
		return "[]"
	}
	var b strings.Builder
	b.WriteByte('[')
	b.WriteString(strconv.FormatInt(c[0], 10))
	if len(c) > 1 {
		b.WriteString("; ")
		for i := 1; i < len(c); i++ {
			if i > 1 {
				b.WriteString(", ")
			}
			b.WriteString(strconv.FormatInt(c[i], 10))
		}
	}
	b.WriteByte(']')
	return b.String()
}

// ParseCF parses the bracket notation produced by [CF.String], for example
// "[3; 7, 15, 1]" or "[3]". Whitespace is ignored and either ';' or ',' may
// separate terms.
func ParseCF(s string) (CF, error) {
	s = strings.TrimSpace(s)
	s = strings.TrimPrefix(s, "[")
	s = strings.TrimSuffix(s, "]")
	s = strings.ReplaceAll(s, ";", ",")
	if strings.TrimSpace(s) == "" {
		return CF{}, nil
	}
	fields := strings.Split(s, ",")
	cf := make(CF, 0, len(fields))
	for _, f := range fields {
		f = strings.TrimSpace(f)
		if f == "" {
			continue
		}
		v, err := strconv.ParseInt(f, 10, 64)
		if err != nil {
			return nil, fmt.Errorf("contfrac: ParseCF: %w", err)
		}
		cf = append(cf, v)
	}
	return cf, nil
}

// Evaluate is a free-function form of [CF.Value].
func Evaluate(cf CF) float64 { return cf.Value() }

// EvaluateRat is a free-function form of [CF.Rat].
func EvaluateRat(cf CF) *big.Rat { return cf.Rat() }

// Convergents is a free-function form of [CF.Convergents].
func Convergents(cf CF) []Frac { return cf.Convergents() }

// Convergent is a free-function form of [CF.Convergent].
func Convergent(cf CF, k int) Frac { return cf.Convergent(k) }

// ConvergentError returns the absolute difference between the value of the full
// continued fraction and its k-th convergent, evaluated in floating point.
func ConvergentError(cf CF, k int) float64 {
	return math.Abs(cf.Value() - cf.Convergent(k).Float())
}

// Continuant returns the continuant K(a0, a1, ..., a_{n-1}), the polynomial
// that equals the numerator of the convergent of the continued fraction with
// those partial quotients. The empty continuant is 1 and K(a0) is a0.
func Continuant(a ...int64) int64 {
	prev, prev2 := int64(1), int64(0)
	for _, x := range a {
		prev, prev2 = x*prev+prev2, prev
	}
	return prev
}

// ContinuantBig is the arbitrary-precision form of [Continuant].
func ContinuantBig(a []int64) *big.Int {
	prev, prev2 := big.NewInt(1), big.NewInt(0)
	for _, x := range a {
		next := new(big.Int).Add(new(big.Int).Mul(big.NewInt(x), prev), prev2)
		prev, prev2 = next, prev
	}
	return prev
}
