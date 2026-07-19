package grouprep

import (
	"errors"
	"fmt"
)

// Rep is a finite-dimensional complex matrix representation of a [Group]: a map
// assigning to every element i an invertible dim×dim matrix mats[i] such that
// mats[i·j] = mats[i]·mats[j]. Build reps with the constructors below or with
// [NewRep], which validates the homomorphism property.
type Rep struct {
	group *Group
	dim   int
	mats  []Matrix
}

// NewRep builds a representation of g from an explicit matrix per element and
// verifies that it is a genuine homomorphism to within tol. The slice mats must
// have length g.Order() and every matrix must be square of the same size. It
// returns an error otherwise.
func NewRep(g *Group, mats []Matrix, tol float64) (*Rep, error) {
	if len(mats) != g.Order() {
		return nil, errors.New("grouprep: NewRep needs one matrix per element")
	}
	if len(mats) == 0 {
		return nil, errors.New("grouprep: empty representation")
	}
	dim := mats[0].Rows()
	for _, m := range mats {
		if m.Rows() != dim || m.Cols() != dim {
			return nil, errors.New("grouprep: all rep matrices must be square of equal size")
		}
	}
	r := &Rep{group: g, dim: dim, mats: mats}
	if !r.IsHomomorphism(tol) {
		return nil, errors.New("grouprep: matrices do not define a homomorphism")
	}
	if !r.Matrix(g.Identity()).IsIdentity(tol) {
		return nil, errors.New("grouprep: identity does not map to the identity matrix")
	}
	return r, nil
}

// Group returns the group being represented.
func (r *Rep) Group() *Group { return r.group }

// Dim returns the dimension (degree) of the representation.
func (r *Rep) Dim() int { return r.dim }

// Matrix returns the representing matrix of element i.
func (r *Rep) Matrix(i int) Matrix { return r.mats[i] }

// IsHomomorphism reports whether mats[i·j] = mats[i]·mats[j] for all i, j to
// within tol.
func (r *Rep) IsHomomorphism(tol float64) bool {
	n := r.group.Order()
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			prod := mustMul(r.mats[i], r.mats[j])
			if !prod.ApproxEqual(r.mats[r.group.Mul(i, j)], tol) {
				return false
			}
		}
	}
	return true
}

// IsUnitary reports whether every representing matrix is unitary to within tol.
// Every representation of a finite group is equivalent to a unitary one
// (Weyl's unitarian trick); the standard constructions here are unitary.
func (r *Rep) IsUnitary(tol float64) bool {
	for _, m := range r.mats {
		if !m.IsUnitary(tol) {
			return false
		}
	}
	return true
}

// Character returns the character of r, the class function χ(i) = tr mats[i].
func (r *Rep) Character() Character {
	ch := make(Character, r.group.Order())
	for i := range r.mats {
		ch[i] = r.mats[i].Trace()
	}
	return ch
}

// DirectSum returns the direct sum r⊕s, the block-diagonal representation of
// dimension Dim(r)+Dim(s). It panics if the two reps are of different groups.
func (r *Rep) DirectSum(s *Rep) *Rep {
	if r.group != s.group {
		panic("grouprep: DirectSum requires representations of the same group")
	}
	mats := make([]Matrix, r.group.Order())
	for i := range mats {
		mats[i] = r.mats[i].DirectSum(s.mats[i])
	}
	return &Rep{group: r.group, dim: r.dim + s.dim, mats: mats}
}

// TensorProduct returns the tensor (Kronecker) product r⊗s, of dimension
// Dim(r)·Dim(s). Its character is the pointwise product of the two characters.
// It panics if the reps are of different groups.
func (r *Rep) TensorProduct(s *Rep) *Rep {
	if r.group != s.group {
		panic("grouprep: TensorProduct requires representations of the same group")
	}
	mats := make([]Matrix, r.group.Order())
	for i := range mats {
		mats[i] = r.mats[i].Kronecker(s.mats[i])
	}
	return &Rep{group: r.group, dim: r.dim * s.dim, mats: mats}
}

// Dual returns the dual (contragredient) representation, sending g to the
// inverse transpose of mats[g]. Its character is the complex conjugate of the
// original character.
func (r *Rep) Dual() *Rep {
	mats := make([]Matrix, r.group.Order())
	for i := range r.mats {
		inv, err := r.mats[i].Inverse()
		if err != nil {
			panic("grouprep: Dual of a non-invertible representation matrix")
		}
		mats[i] = inv.Transpose()
	}
	return &Rep{group: r.group, dim: r.dim, mats: mats}
}

// DirectSumReps returns the direct sum of several representations of the same
// group, in order. It panics if fewer than one rep is given.
func DirectSumReps(reps ...*Rep) *Rep {
	if len(reps) == 0 {
		panic("grouprep: DirectSumReps needs at least one representation")
	}
	out := reps[0]
	for _, s := range reps[1:] {
		out = out.DirectSum(s)
	}
	return out
}

// TrivialRep returns the one-dimensional representation sending every element
// to [1].
func TrivialRep(g *Group) *Rep {
	mats := make([]Matrix, g.Order())
	for i := range mats {
		mats[i] = Matrix{{1}}
	}
	return &Rep{group: g, dim: 1, mats: mats}
}

// RegularRep returns the left regular representation of g, of dimension |G|:
// element x acts on the basis {e_y} by e_y ↦ e_{x·y}. Its character is |G| at
// the identity and 0 elsewhere.
func RegularRep(g *Group) *Rep {
	n := g.Order()
	mats := make([]Matrix, n)
	for x := 0; x < n; x++ {
		perm := make([]int, n)
		for y := 0; y < n; y++ {
			perm[y] = g.Mul(x, y)
		}
		mats[x] = PermutationMatrix(perm)
	}
	return &Rep{group: g, dim: n, mats: mats}
}

// CyclicRep returns the one-dimensional representation of the cyclic group C_n
// sending the generator r to the k-th root of unity ω^k = exp(2πik/n). For k
// ranging over 0..n-1 these are all the irreducible representations of C_n. It
// panics if n < 1.
func CyclicRep(n, k int) *Rep {
	g := CyclicGroup(n)
	mats := make([]Matrix, n)
	for i := 0; i < n; i++ {
		mats[i] = Matrix{{RootOfUnity(n, (k*i)%n)}}
	}
	return &Rep{group: g, dim: 1, mats: mats}
}

// DihedralRep2D returns the two-dimensional representation of the dihedral
// group D_n in which the rotation r acts by rotation through 2πh/n and the
// reflection s acts by a reflection. For 1 <= h < n/2 these are irreducible;
// its character is 2cos(2πhj/n) on r^j and 0 on every reflection. It panics if
// n < 1.
func DihedralRep2D(n, h int) *Rep {
	g := DihedralGroup(n)
	mats := make([]Matrix, 2*n)
	theta := 2 * 3.141592653589793 * float64(h) / float64(n)
	// Reflection s acts by diag-like reflection across the x-axis.
	refl := Matrix{{1, 0}, {0, -1}}
	for a := 0; a < 2; a++ {
		for b := 0; b < n; b++ {
			rot := RotationMatrix(float64(b) * theta)
			var m Matrix
			if a == 0 {
				m = rot
			} else {
				// ρ(s^a r^b) = F^a R^b, so the reflection acts on the left.
				m = mustMul(refl, rot)
			}
			mats[a*n+b] = m
		}
	}
	return &Rep{group: g, dim: 2, mats: mats}
}

// NaturalRepSymmetric returns the natural n-dimensional permutation
// representation of the symmetric group S_n, in which each permutation acts by
// its permutation matrix. It is reducible for n >= 2, decomposing as the
// trivial representation plus the (n-1)-dimensional standard representation. It
// panics if n < 1.
func NaturalRepSymmetric(n int) *Rep {
	g, perms := symmetricData(n)
	mats := make([]Matrix, len(perms))
	for i, p := range perms {
		mats[i] = p.PermMatrix()
	}
	return &Rep{group: g, dim: n, mats: mats}
}

// SignRepSymmetric returns the one-dimensional sign representation of S_n,
// sending each permutation to its sign (+1 or -1). It panics if n < 1.
func SignRepSymmetric(n int) *Rep {
	g, perms := symmetricData(n)
	mats := make([]Matrix, len(perms))
	for i, p := range perms {
		mats[i] = Matrix{{complex(float64(p.Sign()), 0)}}
	}
	return &Rep{group: g, dim: 1, mats: mats}
}

// QuaternionRep returns the two-dimensional irreducible representation of the
// quaternion group Q8 by Pauli-like matrices, with i, j, k acting as the
// standard unit-quaternion matrices. Its character is (2, -2, 0, 0, 0) on the
// classes {1}, {-1}, {±i}, {±j}, {±k}.
func QuaternionRep() *Rep {
	g := QuaternionGroup()
	I := complex(0, 1)
	one := Matrix{{1, 0}, {0, 1}}
	matI := Matrix{{I, 0}, {0, -I}}
	matJ := Matrix{{0, 1}, {-1, 0}}
	matK := Matrix{{0, I}, {I, 0}}
	neg := func(m Matrix) Matrix { return m.Scale(-1) }
	// index: 0:1,1:-1,2:i,3:-i,4:j,5:-j,6:k,7:-k
	mats := []Matrix{one, neg(one), matI, neg(matI), matJ, neg(matJ), matK, neg(matK)}
	return &Rep{group: g, dim: 2, mats: mats}
}

// RepString renders a short multi-line description of r for debugging.
func (r *Rep) RepString() string {
	return fmt.Sprintf("Rep of %s, dim %d, %d elements", r.group.Name(), r.dim, r.group.Order())
}
