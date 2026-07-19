package fuzzy

import (
	"errors"
	"math"
)

// ErrRelationShape is returned when a matrix of grades does not match the
// declared universes of a relation.
var ErrRelationShape = errors.New("fuzzy: relation matrix shape mismatch")

// ErrComposeShape is returned when two relations cannot be composed because the
// column universe of the first differs in size from the row universe of the
// second.
var ErrComposeShape = errors.New("fuzzy: relation composition shape mismatch")

// Relation is a binary fuzzy relation on the product of universes X and Y. M[i]
// [j] is the grade to which X[i] is related to Y[j]. All grades are within
// [0, 1].
type Relation struct {
	X []float64
	Y []float64
	M [][]float64
}

// NewRelation builds a fuzzy relation from universes x and y and a grade matrix
// m of shape len(x) by len(y). Grades are clamped to [0, 1]. It returns
// ErrRelationShape when the matrix dimensions do not match the universes.
func NewRelation(x, y []float64, m [][]float64) (Relation, error) {
	if len(m) != len(x) {
		return Relation{}, ErrRelationShape
	}
	mm := make([][]float64, len(x))
	for i := range m {
		if len(m[i]) != len(y) {
			return Relation{}, ErrRelationShape
		}
		row := make([]float64, len(y))
		for j := range m[i] {
			row[j] = clamp01(m[i][j])
		}
		mm[i] = row
	}
	xs := make([]float64, len(x))
	ys := make([]float64, len(y))
	copy(xs, x)
	copy(ys, y)
	return Relation{X: xs, Y: ys, M: mm}, nil
}

// RelationFromFunc builds a fuzzy relation whose grade at (X[i], Y[j]) is
// f(X[i], Y[j]), clamped to [0, 1].
func RelationFromFunc(x, y []float64, f func(a, b float64) float64) Relation {
	mm := make([][]float64, len(x))
	for i := range x {
		row := make([]float64, len(y))
		for j := range y {
			row[j] = clamp01(f(x[i], y[j]))
		}
		mm[i] = row
	}
	xs := make([]float64, len(x))
	ys := make([]float64, len(y))
	copy(xs, x)
	copy(ys, y)
	return Relation{X: xs, Y: ys, M: mm}
}

// CartesianProduct returns the fuzzy relation a x b formed by the t-norm tn,
// with grade tn(Mu_a(X[i]), Mu_b(Y[j])). Passing TNormMin gives the classical
// min Cartesian product used to build fuzzy rules.
func CartesianProduct(a, b Set, tn TNorm) Relation {
	mm := make([][]float64, len(a.X))
	for i := range a.X {
		row := make([]float64, len(b.X))
		for j := range b.X {
			row[j] = clamp01(tn(a.Mu[i], b.Mu[j]))
		}
		mm[i] = row
	}
	xs := make([]float64, len(a.X))
	ys := make([]float64, len(b.X))
	copy(xs, a.X)
	copy(ys, b.X)
	return Relation{X: xs, Y: ys, M: mm}
}

// Rows reports the number of rows (size of universe X).
func (r Relation) Rows() int { return len(r.X) }

// Cols reports the number of columns (size of universe Y).
func (r Relation) Cols() int { return len(r.Y) }

// At returns the grade of the (i, j) cell.
func (r Relation) At(i, j int) float64 { return r.M[i][j] }

// Clone returns a deep copy of the relation.
func (r Relation) Clone() Relation {
	x := make([]float64, len(r.X))
	y := make([]float64, len(r.Y))
	copy(x, r.X)
	copy(y, r.Y)
	mm := make([][]float64, len(r.M))
	for i := range r.M {
		row := make([]float64, len(r.M[i]))
		copy(row, r.M[i])
		mm[i] = row
	}
	return Relation{X: x, Y: y, M: mm}
}

// Transpose returns the inverse relation with X and Y swapped.
func (r Relation) Transpose() Relation {
	mm := make([][]float64, len(r.Y))
	for j := range r.Y {
		row := make([]float64, len(r.X))
		for i := range r.X {
			row[i] = r.M[i][j]
		}
		mm[j] = row
	}
	x := make([]float64, len(r.X))
	y := make([]float64, len(r.Y))
	copy(x, r.X)
	copy(y, r.Y)
	return Relation{X: y, Y: x, M: mm}
}

// Complement returns the standard complement 1 - M of the relation.
func (r Relation) Complement() Relation {
	out := r.Clone()
	for i := range out.M {
		for j := range out.M[i] {
			out.M[i][j] = clamp01(1 - out.M[i][j])
		}
	}
	return out
}

// Union returns the elementwise union (max) of two relations of equal shape. It
// returns ErrRelationShape when the shapes differ.
func (r Relation) Union(o Relation) (Relation, error) { return r.combine(o, TConormMax) }

// Intersection returns the elementwise intersection (min) of two relations of
// equal shape. It returns ErrRelationShape when the shapes differ.
func (r Relation) Intersection(o Relation) (Relation, error) { return r.combine(o, TNormMin) }

// combine applies the binary operator op elementwise to two relations.
func (r Relation) combine(o Relation, op func(a, b float64) float64) (Relation, error) {
	if len(r.X) != len(o.X) || len(r.Y) != len(o.Y) {
		return Relation{}, ErrRelationShape
	}
	out := r.Clone()
	for i := range out.M {
		for j := range out.M[i] {
			out.M[i][j] = clamp01(op(r.M[i][j], o.M[i][j]))
		}
	}
	return out, nil
}

// MaxMinComposition returns the max-min composition r o s, whose (i, k) grade is
// max_j min(r[i][j], s[j][k]). It returns ErrComposeShape when the inner
// universes do not agree in size.
func (r Relation) MaxMinComposition(s Relation) (Relation, error) {
	return r.MaxTComposition(s, TNormMin)
}

// MaxProductComposition returns the max-product composition r o s, whose (i, k)
// grade is max_j (r[i][j] * s[j][k]). It returns ErrComposeShape when the inner
// universes do not agree in size.
func (r Relation) MaxProductComposition(s Relation) (Relation, error) {
	return r.MaxTComposition(s, TNormProduct)
}

// MaxTComposition returns the sup-t composition r o s using the t-norm tn,
// whose (i, k) grade is max_j tn(r[i][j], s[j][k]). It returns ErrComposeShape
// when the inner universes do not agree in size.
func (r Relation) MaxTComposition(s Relation, tn TNorm) (Relation, error) {
	if len(r.Y) != len(s.X) {
		return Relation{}, ErrComposeShape
	}
	mm := make([][]float64, len(r.X))
	for i := range r.X {
		row := make([]float64, len(s.Y))
		for k := range s.Y {
			best := 0.0
			for j := range r.Y {
				v := tn(r.M[i][j], s.M[j][k])
				if v > best {
					best = v
				}
			}
			row[k] = clamp01(best)
		}
		mm[i] = row
	}
	x := make([]float64, len(r.X))
	y := make([]float64, len(s.Y))
	copy(x, r.X)
	copy(y, s.Y)
	return Relation{X: x, Y: y, M: mm}, nil
}

// ComposeSet returns the fuzzy set on Y obtained by the max-min composition of
// the fuzzy set a (on X) with the relation r, grade_k = max_i min(Mu_a(X[i]),
// r[i][k]). It returns ErrComposeShape when a's universe size differs from the
// relation's row count.
func (r Relation) ComposeSet(a Set) (Set, error) {
	return r.ComposeSetT(a, TNormMin)
}

// ComposeSetT returns the sup-t composition of the fuzzy set a with the
// relation r using the t-norm tn, grade_k = max_i tn(Mu_a(X[i]), r[i][k]). This
// is the compositional rule of inference. It returns ErrComposeShape when a's
// universe size differs from the relation's row count.
func (r Relation) ComposeSetT(a Set, tn TNorm) (Set, error) {
	if len(a.X) != len(r.X) {
		return Set{}, ErrComposeShape
	}
	mu := make([]float64, len(r.Y))
	for k := range r.Y {
		best := 0.0
		for i := range r.X {
			v := tn(a.Mu[i], r.M[i][k])
			if v > best {
				best = v
			}
		}
		mu[k] = clamp01(best)
	}
	y := make([]float64, len(r.Y))
	copy(y, r.Y)
	return Set{X: y, Mu: mu}, nil
}

// ProjectFirst returns the projection of the relation onto its X universe, the
// fuzzy set whose grade at X[i] is max_j r[i][j].
func (r Relation) ProjectFirst() Set {
	mu := make([]float64, len(r.X))
	for i := range r.X {
		best := 0.0
		for j := range r.Y {
			if r.M[i][j] > best {
				best = r.M[i][j]
			}
		}
		mu[i] = best
	}
	x := make([]float64, len(r.X))
	copy(x, r.X)
	return Set{X: x, Mu: mu}
}

// ProjectSecond returns the projection of the relation onto its Y universe, the
// fuzzy set whose grade at Y[j] is max_i r[i][j].
func (r Relation) ProjectSecond() Set {
	mu := make([]float64, len(r.Y))
	for j := range r.Y {
		best := 0.0
		for i := range r.X {
			if r.M[i][j] > best {
				best = r.M[i][j]
			}
		}
		mu[j] = best
	}
	y := make([]float64, len(r.Y))
	copy(y, r.Y)
	return Set{X: y, Mu: mu}
}

// CylindricalExtension returns the cylindrical extension of the fuzzy set a
// (defined on X) onto the product X x y, with grade r[i][j] = Mu_a(X[i]) for
// every j.
func CylindricalExtension(a Set, y []float64) Relation {
	mm := make([][]float64, len(a.X))
	for i := range a.X {
		row := make([]float64, len(y))
		for j := range y {
			row[j] = a.Mu[i]
		}
		mm[i] = row
	}
	x := make([]float64, len(a.X))
	ys := make([]float64, len(y))
	copy(x, a.X)
	copy(ys, y)
	return Relation{X: x, Y: ys, M: mm}
}

// IsReflexive reports whether the relation is reflexive within tol, that is
// whether it is square (X equals Y) and every diagonal grade is 1.
func (r Relation) IsReflexive(tol float64) bool {
	if len(r.X) != len(r.Y) {
		return false
	}
	for i := range r.X {
		if math.Abs(r.M[i][i]-1) > tol {
			return false
		}
	}
	return true
}

// IsSymmetric reports whether the relation is symmetric within tol. It requires
// a square relation.
func (r Relation) IsSymmetric(tol float64) bool {
	if len(r.X) != len(r.Y) {
		return false
	}
	for i := range r.X {
		for j := range r.Y {
			if math.Abs(r.M[i][j]-r.M[j][i]) > tol {
				return false
			}
		}
	}
	return true
}

// IsMaxMinTransitive reports whether the relation is max-min transitive within
// tol, that is r[i][k] >= max_j min(r[i][j], r[j][k]) for all i, k. It requires
// a square relation.
func (r Relation) IsMaxMinTransitive(tol float64) bool {
	if len(r.X) != len(r.Y) {
		return false
	}
	for i := range r.X {
		for k := range r.Y {
			best := 0.0
			for j := range r.Y {
				v := math.Min(r.M[i][j], r.M[j][k])
				if v > best {
					best = v
				}
			}
			if best > r.M[i][k]+tol {
				return false
			}
		}
	}
	return true
}

// IsFuzzyEquivalence reports whether the relation is a fuzzy equivalence
// (similarity) relation within tol: reflexive, symmetric and max-min
// transitive.
func (r Relation) IsFuzzyEquivalence(tol float64) bool {
	return r.IsReflexive(tol) && r.IsSymmetric(tol) && r.IsMaxMinTransitive(tol)
}

// TransitiveClosure returns the max-min transitive closure of a square relation,
// the union of r, r o r, r o r o r, ... which stabilizes after at most n-1
// compositions. It returns ErrComposeShape when the relation is not square.
func (r Relation) TransitiveClosure() (Relation, error) {
	if len(r.X) != len(r.Y) {
		return Relation{}, ErrComposeShape
	}
	cur := r.Clone()
	n := len(r.X)
	for iter := 0; iter < n; iter++ {
		comp, err := cur.MaxMinComposition(r)
		if err != nil {
			return Relation{}, err
		}
		next, err := cur.Union(comp)
		if err != nil {
			return Relation{}, err
		}
		if next.equalWithin(cur, 1e-12) {
			return next, nil
		}
		cur = next
	}
	return cur, nil
}

// equalWithin reports whether two equally shaped relations agree within tol.
func (r Relation) equalWithin(o Relation, tol float64) bool {
	if len(r.X) != len(o.X) || len(r.Y) != len(o.Y) {
		return false
	}
	for i := range r.M {
		for j := range r.M[i] {
			if math.Abs(r.M[i][j]-o.M[i][j]) > tol {
				return false
			}
		}
	}
	return true
}
