package tensor

import "fmt"

// KroneckerDelta returns the Kronecker delta δ_ij: 1 when i == j and 0
// otherwise.
func KroneckerDelta(i, j int) float64 {
	if i == j {
		return 1
	}
	return 0
}

// Identity returns the n×n identity matrix as a rank-2 tensor (the components of
// the Kronecker delta δ_ij). It panics if n is not positive.
func Identity(n int) *Tensor {
	if n <= 0 {
		panic(fmt.Errorf("%w: Identity needs n > 0, got %d", ErrShape, n))
	}
	t := New(n, n)
	for i := 0; i < n; i++ {
		t.Set(1, i, i)
	}
	return t
}

// KroneckerDeltaTensor is a synonym for [Identity]: the n×n matrix of components
// δ_ij.
func KroneckerDeltaTensor(n int) *Tensor { return Identity(n) }

// LeviCivita returns the value of the Levi-Civita symbol ε for the given
// indices. It is +1 if the indices are an even permutation of 0,1,...,n-1, -1 if
// an odd permutation, and 0 if any index repeats or lies outside [0,n) where n
// is the number of indices. Calling it with no indices returns 1 (the empty
// product convention).
func LeviCivita(indices ...int) int {
	n := len(indices)
	seen := make([]bool, n)
	for _, v := range indices {
		if v < 0 || v >= n || seen[v] {
			return 0
		}
		seen[v] = true
	}
	// Count inversions to determine the sign of the permutation.
	sign := 1
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			if indices[i] > indices[j] {
				sign = -sign
			}
		}
	}
	return sign
}

// LeviCivitaTensor returns the rank-n Levi-Civita (permutation) tensor in n
// dimensions: a tensor of shape [n,n,...,n] (n axes) whose component at a
// multi-index equals [LeviCivita] of that index tuple. It panics if n is not
// positive.
func LeviCivitaTensor(n int) *Tensor {
	if n <= 0 {
		panic(fmt.Errorf("%w: LeviCivitaTensor needs n > 0, got %d", ErrShape, n))
	}
	shape := make([]int, n)
	for i := range shape {
		shape[i] = n
	}
	t := New(shape...)
	for f := 0; f < t.Size(); f++ {
		idx := tensorUnravel(f, shape)
		t.data[f] = float64(LeviCivita(idx...))
	}
	return t
}

// EuclideanMetric returns the n×n Euclidean metric tensor, i.e. the identity
// matrix diag(1,...,1). Raising and lowering indices with it leaves components
// unchanged. It panics if n is not positive.
func EuclideanMetric(n int) *Tensor { return Identity(n) }

// MinkowskiMetric returns the dim×dim Minkowski metric tensor with the
// mostly-plus signature (-,+,+,...): component (0,0) is -1 and the remaining
// diagonal entries are +1. This is the flat-spacetime metric of special
// relativity. It panics if dim is not positive.
func MinkowskiMetric(dim int) *Tensor {
	if dim <= 0 {
		panic(fmt.Errorf("%w: MinkowskiMetric needs dim > 0, got %d", ErrShape, dim))
	}
	g := New(dim, dim)
	g.Set(-1, 0, 0)
	for i := 1; i < dim; i++ {
		g.Set(1, i, i)
	}
	return g
}

// tensorApplyMetric contracts the rank-2 tensor m against axis of t, producing a
// tensor of the same shape as t: out[...,a,...] = sum_b m[a][b] * t[...,b,...].
// It is the common engine behind [LowerIndex] and [RaiseIndex].
func tensorApplyMetric(t, m *Tensor, axis int) (*Tensor, error) {
	if m.Rank() != 2 || m.shape[0] != m.shape[1] {
		return nil, fmt.Errorf("%w: metric must be a square matrix, got %v", ErrShape, m.shape)
	}
	r := t.Rank()
	if axis < 0 || axis >= r {
		return nil, fmt.Errorf("%w: axis %d for rank-%d tensor", ErrAxis, axis, r)
	}
	if t.shape[axis] != m.shape[0] {
		return nil, fmt.Errorf("%w: axis length %d does not match metric dimension %d", ErrShape, t.shape[axis], m.shape[0])
	}
	dim := m.shape[0]
	out := New(t.shape...)
	sidx := make([]int, r)
	for f := 0; f < t.Size(); f++ {
		oidx := tensorUnravel(f, t.shape)
		copy(sidx, oidx)
		sum := 0.0
		for b := 0; b < dim; b++ {
			sidx[axis] = b
			sum += m.At(oidx[axis], b) * t.At(sidx...)
		}
		out.data[f] = sum
	}
	return out, nil
}

// LowerIndex lowers the index on the given axis of t using the metric g,
// returning a tensor of the same shape whose axis component is
// t'_{...a...} = sum_b g_{ab} t^{...b...}. The metric must be a square matrix
// whose dimension matches the length of that axis. It returns [ErrShape] or
// [ErrAxis] on a mismatch.
func LowerIndex(t, g *Tensor, axis int) (*Tensor, error) {
	return tensorApplyMetric(t, g, axis)
}

// RaiseIndex raises the index on the given axis of t using the inverse metric
// gInv, returning a tensor of the same shape whose axis component is
// t'^{...a...} = sum_b g^{ab} t_{...b...}. The inverse metric must be a square
// matrix whose dimension matches the length of that axis. It returns [ErrShape]
// or [ErrAxis] on a mismatch. For the Euclidean metric raising and lowering
// leave the components unchanged.
func RaiseIndex(t, gInv *Tensor, axis int) (*Tensor, error) {
	return tensorApplyMetric(t, gInv, axis)
}
