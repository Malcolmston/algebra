package contfrac

import (
	"math"
	"strconv"
	"strings"
)

// QuadraticSurd represents the real number (P + sqrt(D))/Q with integer P, Q, D
// and D >= 0, Q != 0. Every quadratic irrational (and every rational) can be
// written in this form, and by Lagrange's theorem its continued fraction is
// eventually periodic; see [QuadraticSurd.CF].
type QuadraticSurd struct {
	P int64 // rational part numerator
	Q int64 // denominator
	D int64 // radicand under the square root (D >= 0)
}

// NewSurd returns the quadratic surd (P + sqrt(D))/Q. It panics if Q == 0 or
// D < 0.
func NewSurd(P, Q, D int64) QuadraticSurd {
	if Q == 0 {
		panic("contfrac: NewSurd requires Q != 0")
	}
	if D < 0 {
		panic("contfrac: NewSurd requires D >= 0")
	}
	return QuadraticSurd{P: P, Q: Q, D: D}
}

// Value returns the surd as a float64.
func (s QuadraticSurd) Value() float64 {
	return (float64(s.P) + math.Sqrt(float64(s.D))) / float64(s.Q)
}

// IsRational reports whether the surd is actually a rational number, i.e. D is a
// perfect square.
func (s QuadraticSurd) IsRational() bool {
	return IsPerfectSquare(s.D)
}

// Conjugate returns the algebraic conjugate (P - sqrt(D))/Q, written in the
// canonical form with a non-negative radical by negating both P and Q.
func (s QuadraticSurd) Conjugate() QuadraticSurd {
	return QuadraticSurd{P: -s.P, Q: -s.Q, D: s.D}
}

// Floor returns the greatest integer not exceeding the surd, computed exactly
// with integer arithmetic (no floating-point rounding error).
func (s QuadraticSurd) Floor() int64 {
	return surdFloor(s.P, s.Q, s.D)
}

// String renders the surd as "(P + sqrt(D))/Q".
func (s QuadraticSurd) String() string {
	return "(" + strconv.FormatInt(s.P, 10) + " + sqrt(" + strconv.FormatInt(s.D, 10) + "))/" + strconv.FormatInt(s.Q, 10)
}

// surdGE reports whether (P + sqrt(D))/Q >= t using exact integer arithmetic
// (D >= 0, Q != 0). It underlies the exact floor of a quadratic surd.
func surdGE(P, Q, D, t int64) bool {
	L := t*Q - P
	if Q > 0 {
		if L <= 0 {
			return true
		}
		return D >= L*L
	}
	// Q < 0: dividing by a negative reverses the inequality.
	if L < 0 {
		return false
	}
	return D <= L*L
}

// surdFloor returns floor((P + sqrt(D))/Q) exactly (D >= 0, Q != 0).
func surdFloor(P, Q, D int64) int64 {
	sq := math.Sqrt(float64(D))
	t := int64(math.Floor((float64(P) + sq) / float64(Q)))
	// Correct any floating-point rounding by exact integer checks.
	for surdGE(P, Q, D, t+1) {
		t++
	}
	for !surdGE(P, Q, D, t) {
		t--
	}
	return t
}

// PeriodicCF is an eventually periodic continued fraction consisting of a
// non-repeating Head followed by an infinitely repeated Period. When Period is
// empty the value is the rational number denoted by Head.
type PeriodicCF struct {
	Head   []int64 // pre-period partial quotients
	Period []int64 // repeating block
}

// PeriodLength returns the length of the repeating block (0 for a rational).
func (pc PeriodicCF) PeriodLength() int {
	return len(pc.Period)
}

// HeadLength returns the length of the non-repeating prefix.
func (pc PeriodicCF) HeadLength() int {
	return len(pc.Head)
}

// IsPurelyPeriodic reports whether the expansion has an empty head and a
// non-empty period, so the whole continued fraction repeats from the start.
func (pc PeriodicCF) IsPurelyPeriodic() bool {
	return len(pc.Head) == 0 && len(pc.Period) > 0
}

// IsRational reports whether the expansion terminates (empty period).
func (pc PeriodicCF) IsRational() bool {
	return len(pc.Period) == 0
}

// Expand returns the first n partial quotients of the (infinite) continued
// fraction: the head followed by as many copies of the period as needed.
func (pc PeriodicCF) Expand(n int) CF {
	if n < 0 {
		n = 0
	}
	out := make(CF, 0, n)
	for _, t := range pc.Head {
		if len(out) >= n {
			return out
		}
		out = append(out, t)
	}
	if len(pc.Period) == 0 {
		return out
	}
	for len(out) < n {
		for _, t := range pc.Period {
			if len(out) >= n {
				break
			}
			out = append(out, t)
		}
	}
	return out
}

// Terms is an alias for Expand returning a plain []int64.
func (pc PeriodicCF) Terms(n int) []int64 {
	return []int64(pc.Expand(n))
}

// Value returns the numeric value of the periodic continued fraction as a
// float64. For a rational (empty period) it is exact up to float rounding; for
// a quadratic irrational it is accurate to full float64 precision.
func (pc PeriodicCF) Value() float64 {
	if len(pc.Period) == 0 {
		return CF(pc.Head).Value()
	}
	return pc.Expand(len(pc.Head) + 64).Value()
}

// Convergents returns the first n convergents of the expansion.
func (pc PeriodicCF) Convergents(n int) []Frac {
	return pc.Expand(n).Convergents()
}

// String renders the expansion as "[h0; h1, ..., (p0, p1, ...)]" with the
// period parenthesised.
func (pc PeriodicCF) String() string {
	var b strings.Builder
	b.WriteByte('[')
	all := append([]int64{}, pc.Head...)
	if len(all) == 0 && len(pc.Period) > 0 {
		// leading term comes from the period
		b.WriteString("(")
		for i, t := range pc.Period {
			if i > 0 {
				b.WriteString(", ")
			}
			b.WriteString(strconv.FormatInt(t, 10))
		}
		b.WriteString(")]")
		return b.String()
	}
	for i, t := range all {
		if i == 0 {
			b.WriteString(strconv.FormatInt(t, 10))
			if len(all) > 1 || len(pc.Period) > 0 {
				b.WriteString("; ")
			}
		} else {
			b.WriteString(strconv.FormatInt(t, 10))
			if i < len(all)-1 || len(pc.Period) > 0 {
				b.WriteString(", ")
			}
		}
	}
	if len(pc.Period) > 0 {
		b.WriteString("(")
		for i, t := range pc.Period {
			if i > 0 {
				b.WriteString(", ")
			}
			b.WriteString(strconv.FormatInt(t, 10))
		}
		b.WriteString(")")
	}
	b.WriteByte(']')
	return b.String()
}

// gsurd is a general quadratic surd (a + b*sqrt(D))/c used for exact arithmetic
// while converting a periodic continued fraction back into a QuadraticSurd.
type gsurd struct {
	a, b, c, D int64
}

func (g gsurd) reduce() gsurd {
	div := GCD(GCD(absInt(g.a), absInt(g.b)), absInt(g.c))
	if div > 1 {
		g.a /= div
		g.b /= div
		g.c /= div
	}
	if g.c < 0 {
		g.a, g.b, g.c = -g.a, -g.b, -g.c
	}
	return g
}

func (g gsurd) addInt(n int64) gsurd {
	g.a += n * g.c
	return g.reduce()
}

func (g gsurd) recip() gsurd {
	// c / (a + b sqrt(D)) = c(a - b sqrt(D)) / (a^2 - b^2 D)
	na := g.c * g.a
	nb := -g.c * g.b
	nc := g.a*g.a - g.b*g.b*g.D
	return gsurd{na, nb, nc, g.D}.reduce()
}

// toSurd converts the general surd to the canonical QuadraticSurd form
// (P + sqrt(Dn))/Q with a non-negative radical.
func (g gsurd) toSurd() QuadraticSurd {
	g = g.reduce()
	Dn := g.b * g.b * g.D
	if g.b >= 0 {
		return QuadraticSurd{P: g.a, Q: g.c, D: Dn}
	}
	return QuadraticSurd{P: -g.a, Q: -g.c, D: Dn}
}

// Surd returns the quadratic surd represented by the periodic continued
// fraction, obtained by solving the quadratic fixed-point equation of the
// period and then folding in the head terms. For a rational expansion the
// result has D == 0. It panics if the expansion is neither rational nor
// eventually periodic (which cannot happen for values built by this package).
func (pc PeriodicCF) Surd() QuadraticSurd {
	if len(pc.Period) == 0 {
		f := CF(pc.Head).Frac()
		return QuadraticSurd{P: f.Num, Q: f.Den, D: 0}
	}
	// Convergent matrix of the period: pn/qn and pn1/qn1.
	hPrev, hPrev2 := int64(1), int64(0)
	kPrev, kPrev2 := int64(0), int64(1)
	for _, a := range pc.Period {
		hPrev, hPrev2 = a*hPrev+hPrev2, hPrev
		kPrev, kPrev2 = a*kPrev+kPrev2, kPrev
	}
	pn, pn1 := hPrev, hPrev2
	qn, qn1 := kPrev, kPrev2
	// Purely periodic value x solves qn x^2 + (qn1 - pn) x - pn1 = 0.
	a := pn - qn1
	disc := a*a + 4*pn1*qn
	x := gsurd{a: a, b: 1, c: 2 * qn, D: disc}.reduce()
	for i := len(pc.Head) - 1; i >= 0; i-- {
		x = x.recip().addInt(pc.Head[i])
	}
	return x.toSurd()
}

// CF returns the eventually periodic continued fraction of the surd. If the
// surd is rational the returned PeriodicCF has an empty period.
func (s QuadraticSurd) CF() PeriodicCF {
	// Rational case.
	if IsPerfectSquare(s.D) {
		r := Isqrt(s.D)
		f := NewFrac(s.P+r, s.Q)
		return PeriodicCF{Head: []int64(FromFrac(f))}
	}
	P, Q, D := s.P, s.Q, s.D
	// Ensure Q | (D - P^2), scaling by |Q| if necessary so the standard
	// recurrence stays in integers. Signs are handled by surdFloor and the
	// recurrence, which work for either sign of Q.
	if mod(D-P*P, Q) != 0 {
		aq := absInt(Q)
		P *= aq
		D *= aq * aq
		Q *= aq
	}
	seen := make(map[[2]int64]int)
	var terms []int64
	for i := 0; ; i++ {
		key := [2]int64{P, Q}
		if idx, ok := seen[key]; ok {
			head := append([]int64{}, terms[:idx]...)
			period := append([]int64{}, terms[idx:]...)
			return PeriodicCF{Head: head, Period: period}
		}
		seen[key] = i
		a := surdFloor(P, Q, D)
		terms = append(terms, a)
		P = a*Q - P
		Q = (D - P*P) / Q
		if i > 1_000_000 { // safety valve
			return PeriodicCF{Head: terms}
		}
	}
}

// mod returns the Euclidean remainder of a modulo b (0 <= result < |b|).
func mod(a, b int64) int64 {
	r := a % b
	if r < 0 {
		if b < 0 {
			r -= b
		} else {
			r += b
		}
	}
	return r
}

// SqrtCF returns the continued fraction of sqrt(n) for n >= 0 as a PeriodicCF.
// For a perfect square the period is empty. Otherwise the head is the single
// term floor(sqrt(n)) and the period is the repeating tail, which is known to
// be a palindrome followed by 2*floor(sqrt(n)). It panics for n < 0.
func SqrtCF(n int64) PeriodicCF {
	if n < 0 {
		panic("contfrac: SqrtCF of negative number")
	}
	a0 := Isqrt(n)
	if a0*a0 == n {
		return PeriodicCF{Head: []int64{a0}}
	}
	var period []int64
	m, d, a := int64(0), int64(1), a0
	for {
		m = d*a - m
		d = (n - m*m) / d
		a = (a0 + m) / d
		period = append(period, a)
		if a == 2*a0 {
			break
		}
	}
	return PeriodicCF{Head: []int64{a0}, Period: period}
}

// SqrtCFPeriod returns just the repeating block of the continued fraction of
// sqrt(n) (nil for a perfect square).
func SqrtCFPeriod(n int64) []int64 {
	return SqrtCF(n).Period
}

// SqrtCFPeriodLength returns the length of the period of the continued fraction
// of sqrt(n) (0 for a perfect square).
func SqrtCFPeriodLength(n int64) int {
	return len(SqrtCFPeriod(n))
}

// SqrtCFExpand returns the first count partial quotients of the continued
// fraction of sqrt(n).
func SqrtCFExpand(n int64, count int) CF {
	return SqrtCF(n).Expand(count)
}

// SqrtConvergents returns the first count convergents of sqrt(n) as reduced
// fractions.
func SqrtConvergents(n int64, count int) []Frac {
	return SqrtCF(n).Expand(count).Convergents()
}

// SqrtConvergent returns the k-th convergent (0-indexed) of sqrt(n).
func SqrtConvergent(n int64, k int) Frac {
	return SqrtCF(n).Expand(k + 1).Convergent(k)
}

// IsPeriodPalindrome reports whether the given period (excluding its final
// element) reads the same forwards and backwards, the classical structure of
// the period of sqrt(n). An empty or single-element period is trivially a
// palindrome.
func IsPeriodPalindrome(period []int64) bool {
	if len(period) <= 1 {
		return true
	}
	body := period[:len(period)-1] // drop the trailing 2*a0
	for i, j := 0, len(body)-1; i < j; i, j = i+1, j-1 {
		if body[i] != body[j] {
			return false
		}
	}
	return true
}
