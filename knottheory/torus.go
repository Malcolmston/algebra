package knottheory

import "fmt"

// This file collects invariants of torus knots and links T(p, q), the curves
// that wind p times around one axis and q times around the other on the surface
// of a torus. When gcd(p, q) = 1 the result is a knot; otherwise it is a link
// with gcd(p, q) components.

// TorusIsKnot reports whether T(p, q) is a knot, that is gcd(|p|,|q|) = 1.
func TorusIsKnot(p, q int) bool { return gcdInt(p, q) == 1 }

// TorusLinkComponents returns the number of components of T(p, q), which is
// gcd(|p|, |q|).
func TorusLinkComponents(p, q int) int {
	g := gcdInt(p, q)
	if g == 0 {
		return 0
	}
	return g
}

// IsTrivialTorusKnot reports whether T(p, q) is the unknot, which happens iff
// |p| <= 1 or |q| <= 1 (and the pair is coprime).
func IsTrivialTorusKnot(p, q int) bool {
	return abs(p) <= 1 || abs(q) <= 1
}

// TorusKnotCrossingNumber returns the crossing number of the torus knot
// T(p, q) with coprime p, q >= 2, namely min(p(q-1), q(p-1)). It returns 0 for
// the unknot.
func TorusKnotCrossingNumber(p, q int) int {
	p, q = abs(p), abs(q)
	if IsTrivialTorusKnot(p, q) {
		return 0
	}
	a := p * (q - 1)
	b := q * (p - 1)
	if a < b {
		return a
	}
	return b
}

// TorusKnotGenus returns the Seifert genus of the torus knot T(p, q) with
// coprime p, q, which is (p-1)(q-1)/2.
func TorusKnotGenus(p, q int) int {
	p, q = abs(p), abs(q)
	return (p - 1) * (q - 1) / 2
}

// TorusKnotSeifertGenus is an alias for TorusKnotGenus.
func TorusKnotSeifertGenus(p, q int) int { return TorusKnotGenus(p, q) }

// TorusKnotUnknottingNumber returns the unknotting number of the torus knot
// T(p, q) with coprime p, q, equal to its genus (p-1)(q-1)/2 by the resolution
// of the Milnor conjecture.
func TorusKnotUnknottingNumber(p, q int) int { return TorusKnotGenus(p, q) }

// TorusKnotBridgeNumber returns the bridge number of T(p, q) with coprime
// p, q >= 2, which is min(|p|, |q|). The unknot has bridge number 1.
func TorusKnotBridgeNumber(p, q int) int {
	p, q = abs(p), abs(q)
	if IsTrivialTorusKnot(p, q) {
		return 1
	}
	if p < q {
		return p
	}
	return q
}

// TorusKnotBraid returns a braid on |p| strands whose closure is T(p, q),
// namely (sigma_1 ... sigma_{p-1})^q.
func TorusKnotBraid(p, q int) (Braid, error) {
	return TorusBraid(abs(p), q)
}

// TorusKnotAlexander returns the Alexander polynomial of the torus knot
// T(p, q) with coprime p, q >= 2, given by the closed form
//
//	Delta(t) = (t^{pq} - 1)(t - 1) / ((t^p - 1)(t^q - 1)),
//
// normalised to the symmetric representative with Delta(1) = 1. It returns an
// error if p and q are not coprime.
func TorusKnotAlexander(p, q int) (Laurent, error) {
	p, q = abs(p), abs(q)
	if gcdInt(p, q) != 1 {
		return Laurent{}, fmt.Errorf("knottheory: TorusKnotAlexander requires coprime p,q")
	}
	if IsTrivialTorusKnot(p, q) {
		return OneLaurent(), nil
	}
	num := tPowMinus1(p * q).Mul(tPowMinus1(1))
	den := tPowMinus1(p).Mul(tPowMinus1(q))
	quo, ok := num.DivExact(den)
	if !ok {
		return Laurent{}, fmt.Errorf("knottheory: torus Alexander division was not exact")
	}
	return normalizeAlexander(quo), nil
}

// TorusKnotJones returns the Jones polynomial of the torus knot T(p, q) with
// coprime p, q >= 2, in the variable t, using the closed form
//
//	V(t) = t^{(p-1)(q-1)/2} (1 - t^{p+1} - t^{q+1} + t^{p+q}) / (1 - t^2).
//
// It returns an error if p and q are not coprime.
func TorusKnotJones(p, q int) (Laurent, error) {
	p, q = abs(p), abs(q)
	if gcdInt(p, q) != 1 {
		return Laurent{}, fmt.Errorf("knottheory: TorusKnotJones requires coprime p,q")
	}
	if IsTrivialTorusKnot(p, q) {
		return OneLaurent(), nil
	}
	// numerator 1 - t^{p+1} - t^{q+1} + t^{p+q}
	num := OneLaurent().
		Sub(Monomial(1, p+1)).
		Sub(Monomial(1, q+1)).
		Add(Monomial(1, p+q))
	den := OneLaurent().Sub(Monomial(1, 2)) // 1 - t^2
	quo, ok := num.DivExact(den)
	if !ok {
		return Laurent{}, fmt.Errorf("knottheory: torus Jones division was not exact")
	}
	shift := (p - 1) * (q - 1) / 2
	return quo.ShiftExp(shift), nil
}

// TorusKnotDeterminant returns the determinant of the torus knot T(p, q), the
// absolute value of its Alexander polynomial at t = -1. For coprime p, q it is
// 1 when both are odd and p (equivalently q) when exactly one is even.
func TorusKnotDeterminant(p, q int) int {
	alex, err := TorusKnotAlexander(p, q)
	if err != nil || alex.IsZero() {
		return 0
	}
	v := alex.EvalUnit(-1)
	if v < 0 {
		return -v
	}
	return v
}

// tPowMinus1 returns the Laurent polynomial t^k - 1.
func tPowMinus1(k int) Laurent {
	return Monomial(1, k).Sub(OneLaurent())
}
