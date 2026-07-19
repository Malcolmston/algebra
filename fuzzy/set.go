package fuzzy

import (
	"errors"
	"fmt"
	"math"
	"sort"
	"strings"
)

// ErrLengthMismatch is returned when a universe and a membership slice have
// different lengths.
var ErrLengthMismatch = errors.New("fuzzy: universe and membership length mismatch")

// ErrEmptySet is returned by operations that require a non-empty universe.
var ErrEmptySet = errors.New("fuzzy: empty universe")

// ErrDomainMismatch is returned when two sets defined on different universes
// are combined by an operation that requires identical domains.
var ErrDomainMismatch = errors.New("fuzzy: universe (domain) mismatch")

// Set is a discrete fuzzy set: a sorted universe of discourse X together with
// the membership grade Mu[i] of X[i]. Grades are kept within [0, 1] by the
// constructors and operations of this package.
type Set struct {
	X  []float64
	Mu []float64
}

// NewSet builds a discrete fuzzy set from parallel universe and membership
// slices. The points are sorted ascending, duplicates keep their maximum grade
// and all grades are clamped to [0, 1]. It returns ErrLengthMismatch when the
// slices differ in length.
func NewSet(x, mu []float64) (Set, error) {
	if len(x) != len(mu) {
		return Set{}, ErrLengthMismatch
	}
	type pair struct{ x, mu float64 }
	ps := make([]pair, len(x))
	for i := range x {
		ps[i] = pair{x[i], clamp01(mu[i])}
	}
	sort.Slice(ps, func(i, j int) bool { return ps[i].x < ps[j].x })
	xs := make([]float64, 0, len(ps))
	ms := make([]float64, 0, len(ps))
	for _, p := range ps {
		if n := len(xs); n > 0 && xs[n-1] == p.x {
			if p.mu > ms[n-1] {
				ms[n-1] = p.mu
			}
			continue
		}
		xs = append(xs, p.x)
		ms = append(ms, p.mu)
	}
	return Set{X: xs, Mu: ms}, nil
}

// FromMF samples the membership function mf over the universe xs and returns
// the resulting discrete fuzzy set.
func FromMF(mf MF, xs []float64) Set {
	s, _ := NewSet(xs, SampleMF(mf, xs))
	return s
}

// Len reports the number of points in the set's universe.
func (s Set) Len() int { return len(s.X) }

// Clone returns a deep copy of the set.
func (s Set) Clone() Set {
	x := make([]float64, len(s.X))
	mu := make([]float64, len(s.Mu))
	copy(x, s.X)
	copy(mu, s.Mu)
	return Set{X: x, Mu: mu}
}

// Membership returns the grade of x. If x is one of the universe points its
// stored grade is returned; otherwise the grade is linearly interpolated
// between neighboring points and is 0 outside the universe.
func (s Set) Membership(x float64) float64 {
	n := len(s.X)
	if n == 0 {
		return 0
	}
	if x <= s.X[0] {
		if x == s.X[0] {
			return s.Mu[0]
		}
		return 0
	}
	if x >= s.X[n-1] {
		if x == s.X[n-1] {
			return s.Mu[n-1]
		}
		return 0
	}
	i := sort.SearchFloat64s(s.X, x)
	if i < n && s.X[i] == x {
		return s.Mu[i]
	}
	// interpolate between i-1 and i
	x0, x1 := s.X[i-1], s.X[i]
	m0, m1 := s.Mu[i-1], s.Mu[i]
	if x1 == x0 {
		return math.Max(m0, m1)
	}
	t := (x - x0) / (x1 - x0)
	return clamp01(m0 + t*(m1-m0))
}

// Height returns the largest membership grade in the set, 0 for an empty set.
func (s Set) Height() float64 {
	h := 0.0
	for _, m := range s.Mu {
		if m > h {
			h = m
		}
	}
	return h
}

// IsNormal reports whether the set attains grade 1 at some point.
func (s Set) IsNormal() bool { return s.Height() >= 1 }

// IsEmpty reports whether every membership grade is 0.
func (s Set) IsEmpty() bool {
	for _, m := range s.Mu {
		if m > 0 {
			return false
		}
	}
	return true
}

// Normalize returns the set rescaled so its height becomes 1, dividing every
// grade by the current height. An empty (all-zero) set is returned unchanged.
func (s Set) Normalize() Set {
	h := s.Height()
	out := s.Clone()
	if h == 0 {
		return out
	}
	for i := range out.Mu {
		out.Mu[i] = clamp01(out.Mu[i] / h)
	}
	return out
}

// Support returns the universe points whose grade is strictly positive, the
// crisp support of the set.
func (s Set) Support() []float64 {
	var out []float64
	for i, m := range s.Mu {
		if m > 0 {
			out = append(out, s.X[i])
		}
	}
	return out
}

// Core returns the universe points whose grade equals 1, the crisp core of the
// set.
func (s Set) Core() []float64 {
	var out []float64
	for i, m := range s.Mu {
		if m >= 1 {
			out = append(out, s.X[i])
		}
	}
	return out
}

// AlphaCut returns the universe points whose grade is at least alpha, the
// (weak) alpha-cut of the set.
func (s Set) AlphaCut(alpha float64) []float64 {
	var out []float64
	for i, m := range s.Mu {
		if m >= alpha {
			out = append(out, s.X[i])
		}
	}
	return out
}

// StrongAlphaCut returns the universe points whose grade is strictly greater
// than alpha, the strong alpha-cut of the set.
func (s Set) StrongAlphaCut(alpha float64) []float64 {
	var out []float64
	for i, m := range s.Mu {
		if m > alpha {
			out = append(out, s.X[i])
		}
	}
	return out
}

// Cardinality returns the scalar (sigma) cardinality of the set, the sum of all
// membership grades.
func (s Set) Cardinality() float64 {
	sum := 0.0
	for _, m := range s.Mu {
		sum += m
	}
	return sum
}

// RelativeCardinality returns the scalar cardinality divided by the number of
// universe points, the mean membership grade. It returns 0 for an empty
// universe.
func (s Set) RelativeCardinality() float64 {
	if len(s.X) == 0 {
		return 0
	}
	return s.Cardinality() / float64(len(s.X))
}

// Complement returns the standard complement 1 - Mu of the set on the same
// universe.
func (s Set) Complement() Set {
	out := s.Clone()
	for i := range out.Mu {
		out.Mu[i] = clamp01(1 - out.Mu[i])
	}
	return out
}

// ComplementWith returns the complement of the set using the fuzzy complement
// operator c.
func (s Set) ComplementWith(c func(float64) float64) Set {
	out := s.Clone()
	for i := range out.Mu {
		out.Mu[i] = clamp01(c(out.Mu[i]))
	}
	return out
}

// sameDomain reports whether two sets share an identical universe.
func sameDomain(a, b Set) bool {
	if len(a.X) != len(b.X) {
		return false
	}
	for i := range a.X {
		if a.X[i] != b.X[i] {
			return false
		}
	}
	return true
}

// Union returns the union (max) of two sets defined on the same universe. It
// returns ErrDomainMismatch if the universes differ.
func (s Set) Union(t Set) (Set, error) { return s.UnionWith(t, TConormMax) }

// Intersection returns the intersection (min) of two sets defined on the same
// universe. It returns ErrDomainMismatch if the universes differ.
func (s Set) Intersection(t Set) (Set, error) { return s.IntersectionWith(t, TNormMin) }

// UnionWith returns the union of two sets using the t-conorm sn. It returns
// ErrDomainMismatch if the universes differ.
func (s Set) UnionWith(t Set, sn TConorm) (Set, error) {
	if !sameDomain(s, t) {
		return Set{}, ErrDomainMismatch
	}
	out := s.Clone()
	for i := range out.Mu {
		out.Mu[i] = clamp01(sn(s.Mu[i], t.Mu[i]))
	}
	return out, nil
}

// IntersectionWith returns the intersection of two sets using the t-norm tn. It
// returns ErrDomainMismatch if the universes differ.
func (s Set) IntersectionWith(t Set, tn TNorm) (Set, error) {
	if !sameDomain(s, t) {
		return Set{}, ErrDomainMismatch
	}
	out := s.Clone()
	for i := range out.Mu {
		out.Mu[i] = clamp01(tn(s.Mu[i], t.Mu[i]))
	}
	return out, nil
}

// Difference returns the bounded fuzzy difference s AND NOT t (min of s and the
// complement of t) on a shared universe. It returns ErrDomainMismatch if the
// universes differ.
func (s Set) Difference(t Set) (Set, error) {
	return s.IntersectionWith(t.Complement(), TNormMin)
}

// Pow returns the set with every grade raised to the power p, the basis of the
// concentration (p>1) and dilation (0<p<1) hedges.
func (s Set) Pow(p float64) Set {
	out := s.Clone()
	for i := range out.Mu {
		out.Mu[i] = clamp01(math.Pow(out.Mu[i], p))
	}
	return out
}

// Concentrate returns the concentration hedge Mu^2 ("very").
func (s Set) Concentrate() Set { return s.Pow(2) }

// Dilate returns the dilation hedge Mu^0.5 ("somewhat").
func (s Set) Dilate() Set { return s.Pow(0.5) }

// Very returns the "very" hedge, Mu^2.
func (s Set) Very() Set { return s.Pow(2) }

// Extremely returns the "extremely" hedge, Mu^3.
func (s Set) Extremely() Set { return s.Pow(3) }

// Somewhat returns the "somewhat" hedge, Mu^0.5.
func (s Set) Somewhat() Set { return s.Pow(0.5) }

// Slightly returns the "slightly" hedge, Mu^(1/3).
func (s Set) Slightly() Set { return s.Pow(1.0 / 3.0) }

// MoreOrLess returns the "more or less" hedge, Mu^0.5.
func (s Set) MoreOrLess() Set { return s.Pow(0.5) }

// Intensify returns the contrast intensification hedge applied to every grade.
func (s Set) Intensify() Set {
	out := s.Clone()
	for i, m := range s.Mu {
		if m <= 0.5 {
			out.Mu[i] = clamp01(2 * m * m)
		} else {
			d := 1 - m
			out.Mu[i] = clamp01(1 - 2*d*d)
		}
	}
	return out
}

// Not returns the standard complement, a synonym for Complement suited to
// linguistic negation.
func (s Set) Not() Set { return s.Complement() }

// IsConvex reports whether the set is fuzzy-convex, that is whether its
// membership sequence over the sorted universe is unimodal (never rises after
// it has fallen). Points must be sorted, which the constructors guarantee.
func (s Set) IsConvex() bool {
	falling := false
	const eps = 1e-12
	for i := 1; i < len(s.Mu); i++ {
		if s.Mu[i] < s.Mu[i-1]-eps {
			falling = true
		} else if s.Mu[i] > s.Mu[i-1]+eps && falling {
			return false
		}
	}
	return true
}

// Equal reports whether two sets share the same universe and have membership
// grades that agree within tol.
func (s Set) Equal(t Set, tol float64) bool {
	if !sameDomain(s, t) {
		return false
	}
	for i := range s.Mu {
		if math.Abs(s.Mu[i]-t.Mu[i]) > tol {
			return false
		}
	}
	return true
}

// IsSubset reports whether s is a fuzzy subset of t, that is Mu_s(x) <= Mu_t(x)
// within tol at every shared universe point. It returns false if the universes
// differ.
func (s Set) IsSubset(t Set, tol float64) bool {
	if !sameDomain(s, t) {
		return false
	}
	for i := range s.Mu {
		if s.Mu[i] > t.Mu[i]+tol {
			return false
		}
	}
	return true
}

// DegreeOfSubsethood returns the Kosko degree to which s is a subset of t,
// cardinality(s AND t) / cardinality(s). It returns 1 for an empty s.
func (s Set) DegreeOfSubsethood(t Set) (float64, error) {
	if !sameDomain(s, t) {
		return 0, ErrDomainMismatch
	}
	num := 0.0
	den := 0.0
	for i := range s.Mu {
		num += math.Min(s.Mu[i], t.Mu[i])
		den += s.Mu[i]
	}
	if den == 0 {
		return 1, nil
	}
	return num / den, nil
}

// String renders the set in Zadeh's notation, mu1/x1 + mu2/x2 + ..., listing
// only points with positive grade.
func (s Set) String() string {
	var b strings.Builder
	first := true
	for i, m := range s.Mu {
		if m <= 0 {
			continue
		}
		if !first {
			b.WriteString(" + ")
		}
		fmt.Fprintf(&b, "%.4g/%.4g", m, s.X[i])
		first = false
	}
	if first {
		return "{}"
	}
	return b.String()
}
