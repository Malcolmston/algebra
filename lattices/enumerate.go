package lattices

import "math"

// EnumResult holds the outcome of a lattice enumeration search: the lattice
// Vector found, its integer coordinate vector Coeffs in the given basis, and
// the squared Euclidean quantity Norm2 that was minimized (the squared norm for
// a shortest-vector search or the squared distance to the target for a
// closest-vector search).
type EnumResult struct {
	Vector Vec
	Coeffs []int64
	Norm2  float64
}

// Norm returns the Euclidean norm (or distance) sqrt(Norm2) of the result.
func (r EnumResult) Norm() float64 { return math.Sqrt(r.Norm2) }

// enumerate runs a Fincke-Pohst style enumeration over integer coefficient
// vectors. mu and norm2 are the Gram-Schmidt data of the basis, center holds
// the Gram-Schmidt coordinates of the target (all zero for a shortest-vector
// search), r2 is the (squared) search radius, and skipZero forces a nonzero
// coefficient vector. It returns the best coefficient vector, its objective
// value and whether anything was found.
func enumerate(mu [][]float64, norm2, center []float64, r2 float64, skipZero bool) ([]int64, float64, bool) {
	n := len(norm2)
	for i := 0; i < n; i++ {
		if norm2[i] <= 0 {
			return nil, 0, false
		}
	}
	x := make([]int64, n)
	best := make([]int64, n)
	bestVal := math.Inf(1)
	found := false

	var rec func(i int, partial float64)
	rec = func(i int, partial float64) {
		if i < 0 {
			if skipZero {
				allZero := true
				for _, v := range x {
					if v != 0 {
						allZero = false
						break
					}
				}
				if allZero {
					return
				}
			}
			if partial < bestVal {
				bestVal = partial
				copy(best, x)
				found = true
			}
			return
		}
		c := center[i]
		for j := i + 1; j < n; j++ {
			c -= mu[j][i] * float64(x[j])
		}
		rem := r2 - partial
		if rem < 0 {
			return
		}
		width := math.Sqrt(rem / norm2[i])
		lo := int64(math.Ceil(c - width))
		hi := int64(math.Floor(c + width))
		for xi := lo; xi <= hi; xi++ {
			x[i] = xi
			y := float64(xi) - c
			rec(i-1, partial+y*y*norm2[i])
		}
	}
	rec(n-1, 0)
	if !found {
		return nil, 0, false
	}
	return best, bestVal, true
}

// ShortestVector returns the shortest nonzero vector of the lattice by
// Fincke-Pohst enumeration. The basis is LLL reduced first to obtain a good
// starting radius, so the returned Norm2 is the squared first minimum
// lambda_1^2. It returns ErrEmpty for an empty basis and ErrNotFullRank when
// the basis is rank deficient.
func (b Basis) ShortestVector() (EnumResult, error) {
	if len(b) == 0 {
		return EnumResult{}, ErrEmpty
	}
	red := b.LLLDefault()
	gs := red.Orthogonalize()
	for _, n2 := range gs.Norm2 {
		if n2 <= 0 {
			return EnumResult{}, ErrNotFullRank
		}
	}
	// starting radius: shortest reduced basis vector, slightly inflated.
	r2 := math.Inf(1)
	for _, v := range red {
		if v.IsZero() {
			continue
		}
		if s := v.Norm2(); s < r2 {
			r2 = s
		}
	}
	r2 = r2*(1+1e-9) + 1e-12
	center := make([]float64, len(red))
	coeffs, val, ok := enumerate(gs.Mu, gs.Norm2, center, r2, true)
	if !ok {
		return EnumResult{}, ErrNoSolution
	}
	return EnumResult{Vector: red.Point(coeffs), Coeffs: coeffs, Norm2: val}, nil
}

// FirstMinimum returns lambda_1, the length of the shortest nonzero lattice
// vector, computed exactly by enumeration. It returns 0 with an error for an
// empty or rank-deficient basis.
func (b Basis) FirstMinimum() (float64, error) {
	res, err := b.ShortestVector()
	if err != nil {
		return 0, err
	}
	return res.Norm(), nil
}

// ClosestVector returns the lattice vector nearest to target using enumeration
// seeded by Babai's nearest-plane approximation for the search radius. The
// returned Norm2 is the squared distance to the target. It returns ErrEmpty,
// ErrDimMismatch or ErrNotFullRank as appropriate.
func (b Basis) ClosestVector(target Vec) (EnumResult, error) {
	if len(b) == 0 {
		return EnumResult{}, ErrEmpty
	}
	if len(target) != b.Dim() {
		return EnumResult{}, ErrDimMismatch
	}
	gs := b.Orthogonalize()
	for _, n2 := range gs.Norm2 {
		if n2 <= 0 {
			return EnumResult{}, ErrNotFullRank
		}
	}
	// Gram-Schmidt coordinates of the target.
	center := make([]float64, len(b))
	for i := range b {
		center[i] = target.Dot(gs.Star[i]) / gs.Norm2[i]
	}
	// Babai nearest-plane gives an initial upper bound on the distance.
	approx := b.BabaiNearestPlane(target)
	r2 := approx.Dist2(target)*(1+1e-9) + 1e-12
	coeffs, val, ok := enumerate(gs.Mu, gs.Norm2, center, r2, false)
	if !ok {
		return EnumResult{Vector: approx, Norm2: approx.Dist2(target)}, nil
	}
	return EnumResult{Vector: b.Point(coeffs), Coeffs: coeffs, Norm2: val}, nil
}

// BabaiRound returns an approximate closest lattice vector to target using
// Babai's rounding technique: express the target in basis coordinates (in a
// least-squares sense) and round each coordinate to the nearest integer. It is
// fast but only reliable for reasonably orthogonal bases. It panics if the
// dimension of target does not match the ambient dimension.
func (b Basis) BabaiRound(target Vec) Vec {
	if len(target) != b.Dim() {
		panic(ErrDimMismatch)
	}
	coords, err := b.Coordinates(target)
	if err != nil {
		return ZeroVec(b.Dim())
	}
	x := make([]int64, len(coords))
	for i, c := range coords {
		x[i] = int64(math.Round(c))
	}
	return b.Point(x)
}

// BabaiNearestPlane returns an approximate closest lattice vector to target
// using Babai's nearest-plane algorithm, which walks the Gram-Schmidt frame
// from the last vector to the first. It typically produces a closer vector than
// BabaiRound. It panics if the dimension of target does not match the ambient
// dimension.
func (b Basis) BabaiNearestPlane(target Vec) Vec {
	if len(target) != b.Dim() {
		panic(ErrDimMismatch)
	}
	gs := b.Orthogonalize()
	w := target.Clone()
	n := len(b)
	for j := n - 1; j >= 0; j-- {
		if gs.Norm2[j] == 0 {
			continue
		}
		c := math.Round(w.Dot(gs.Star[j]) / gs.Norm2[j])
		w = w.AddScaled(b[j], -c)
	}
	return target.Sub(w)
}

// BabaiCoeffs returns the integer coordinate vector produced by Babai rounding
// of target, so that b.Point(result) equals BabaiRound(target). It panics on a
// dimension mismatch.
func (b Basis) BabaiCoeffs(target Vec) []int64 {
	if len(target) != b.Dim() {
		panic(ErrDimMismatch)
	}
	coords, err := b.Coordinates(target)
	if err != nil {
		return make([]int64, len(b))
	}
	x := make([]int64, len(coords))
	for i, c := range coords {
		x[i] = int64(math.Round(c))
	}
	return x
}

// EnumerateVectors returns all nonzero lattice vectors whose Euclidean norm is
// at most radius, each as an EnumResult, together with a completeness flag that
// is false if the internal limit on the number of returned vectors was reached.
// The search uses Fincke-Pohst enumeration. It returns ErrEmpty or
// ErrNotFullRank as appropriate.
func (b Basis) EnumerateVectors(radius float64, limit int) ([]EnumResult, bool, error) {
	if len(b) == 0 {
		return nil, false, ErrEmpty
	}
	if radius < 0 {
		return nil, false, ErrBadParameter
	}
	gs := b.Orthogonalize()
	for _, n2 := range gs.Norm2 {
		if n2 <= 0 {
			return nil, false, ErrNotFullRank
		}
	}
	if limit <= 0 {
		limit = 1 << 20
	}
	r2 := radius*radius + 1e-12
	n := len(b)
	center := make([]float64, n)
	x := make([]int64, n)
	var out []EnumResult
	complete := true

	var rec func(i int, partial float64)
	rec = func(i int, partial float64) {
		if !complete {
			return
		}
		if i < 0 {
			allZero := true
			for _, v := range x {
				if v != 0 {
					allZero = false
					break
				}
			}
			if allZero {
				return
			}
			coeffs := make([]int64, n)
			copy(coeffs, x)
			out = append(out, EnumResult{Vector: b.Point(coeffs), Coeffs: coeffs, Norm2: partial})
			if len(out) >= limit {
				complete = false
			}
			return
		}
		c := center[i]
		for j := i + 1; j < n; j++ {
			c -= gs.Mu[j][i] * float64(x[j])
		}
		rem := r2 - partial
		if rem < 0 {
			return
		}
		width := math.Sqrt(rem / gs.Norm2[i])
		lo := int64(math.Ceil(c - width))
		hi := int64(math.Floor(c + width))
		for xi := lo; xi <= hi && complete; xi++ {
			x[i] = xi
			y := float64(xi) - c
			rec(i-1, partial+y*y*gs.Norm2[i])
		}
	}
	rec(n-1, 0)
	return out, complete, nil
}

// KissingNumber returns the number of shortest nonzero lattice vectors, that is
// the number of lattice vectors of length exactly lambda_1 (which always comes
// in plus/minus pairs). It returns an error for an empty or rank-deficient
// basis.
func (b Basis) KissingNumber() (int, error) {
	res, err := b.ShortestVector()
	if err != nil {
		return 0, err
	}
	radius := res.Norm() * (1 + 1e-7)
	vecs, complete, err := b.EnumerateVectors(radius, 1<<20)
	if err != nil {
		return 0, err
	}
	if !complete {
		return 0, ErrNoSolution
	}
	tol := res.Norm2 * 1e-6
	count := 0
	for _, v := range vecs {
		if math.Abs(v.Norm2-res.Norm2) <= tol {
			count++
		}
	}
	return count, nil
}
