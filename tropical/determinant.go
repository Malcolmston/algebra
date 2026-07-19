package tropical

import (
	"errors"
	"math"
)

// ErrTooLarge is returned by brute-force routines when the matrix is larger
// than the supported bound.
var ErrTooLarge = errors.New("tropical: matrix too large for brute-force enumeration")

// OptimalAssignment returns the tropical permanent value and an optimal
// permutation for a square matrix. For min-plus the value is the minimum total
// weight of a perfect matching and for max-plus it is the maximum; the returned
// slice maps each row to its assigned column. The value is the tropical zero
// (and the permutation nil) when no finite assignment exists. It uses the
// Hungarian algorithm and returns ErrNotSquare for a non-square matrix.
func (m Matrix) OptimalAssignment() (float64, []int, error) {
	if !m.IsSquare() {
		return 0, nil, ErrNotSquare
	}
	n := m.rows
	if n == 0 {
		return m.sr.One(), []int{}, nil
	}
	// Build a finite minimisation cost matrix. For min-plus minimise A;
	// for max-plus minimise -A. Forbidden entries (tropical zero) get a
	// large penalty so they are avoided when any finite matching exists.
	cost := make([][]float64, n)
	forbidden := make([][]bool, n)
	var mag float64
	for i := 0; i < n; i++ {
		cost[i] = make([]float64, n)
		forbidden[i] = make([]bool, n)
		for j := 0; j < n; j++ {
			a := m.data[i][j]
			if m.sr.IsZero(a) {
				forbidden[i][j] = true
				continue
			}
			c := a
			if m.sr.IsMaxPlus() {
				c = -a
			}
			cost[i][j] = c
			if av := math.Abs(c); av > mag {
				mag = av
			}
		}
	}
	bigM := (mag+1)*float64(n) + 1
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if forbidden[i][j] {
				cost[i][j] = bigM
			}
		}
	}
	_, assign := hungarianMin(cost)
	// Recover the tropical value from the original entries.
	value := m.sr.One()
	for i := 0; i < n; i++ {
		j := assign[i]
		if forbidden[i][j] {
			return m.sr.Zero(), nil, nil
		}
		value = m.sr.Mul(value, m.data[i][j])
	}
	return value, assign, nil
}

// Permanent returns the tropical permanent of a square matrix, the tropical sum
// over all permutations sigma of the tropical product of the entries A[i,
// sigma(i)]. This equals the optimal assignment value: the minimum-weight
// perfect matching for min-plus and the maximum-weight one for max-plus. It
// returns ErrNotSquare for a non-square matrix.
func (m Matrix) Permanent() (float64, error) {
	v, _, err := m.OptimalAssignment()
	return v, err
}

// Determinant returns the tropical determinant of a square matrix. In tropical
// algebra there are no additive inverses, so the determinant carries no signs
// and coincides with the tropical permanent (the optimal assignment value). It
// returns ErrNotSquare for a non-square matrix.
func (m Matrix) Determinant() (float64, error) { return m.Permanent() }

// PermanentBrute returns the tropical permanent computed by explicitly
// enumerating every permutation. It handles arbitrary entries, including the
// tropical zero, and is intended for cross-checking and for matrices whose
// forbidden pattern makes the assignment infeasible. It returns ErrTooLarge for
// matrices larger than 10-by-10 and ErrNotSquare for a non-square matrix.
func (m Matrix) PermanentBrute() (float64, error) {
	if !m.IsSquare() {
		return 0, ErrNotSquare
	}
	n := m.rows
	if n > 10 {
		return 0, ErrTooLarge
	}
	best := m.sr.Zero()
	perm := make([]int, n)
	for i := range perm {
		perm[i] = i
	}
	permute(perm, 0, func(p []int) {
		prod := m.sr.One()
		for i := 0; i < n; i++ {
			prod = m.sr.Mul(prod, m.data[i][p[i]])
		}
		best = m.sr.Add(best, prod)
	})
	return best, nil
}

// IsTropicallySingular reports whether a square matrix is tropically singular:
// its optimal assignment value is the tropical zero (no finite matching) or is
// attained by at least two distinct permutations, agreeing to within tol. It
// enumerates permutations and returns ErrTooLarge above 9-by-9 and
// ErrNotSquare for a non-square matrix.
func (m Matrix) IsTropicallySingular(tol float64) (bool, error) {
	if !m.IsSquare() {
		return false, ErrNotSquare
	}
	n := m.rows
	if n > 9 {
		return false, ErrTooLarge
	}
	if n == 0 {
		return false, nil
	}
	perm := make([]int, n)
	for i := range perm {
		perm[i] = i
	}
	// Pass 1: best (finite) assignment value.
	best := m.sr.Zero()
	permute(perm, 0, func(p []int) {
		prod := m.sr.One()
		for i := 0; i < n; i++ {
			prod = m.sr.Mul(prod, m.data[i][p[i]])
		}
		if !m.sr.IsZero(prod) {
			best = m.sr.Add(best, prod)
		}
	})
	if m.sr.IsZero(best) {
		return true, nil
	}
	// Pass 2: count finite assignments attaining the best value.
	count := 0
	permute(perm, 0, func(p []int) {
		prod := m.sr.One()
		for i := 0; i < n; i++ {
			prod = m.sr.Mul(prod, m.data[i][p[i]])
		}
		if !m.sr.IsZero(prod) && math.Abs(prod-best) <= tol {
			count++
		}
	})
	return count >= 2, nil
}

// hungarianMin solves the square minimum-cost assignment problem with the
// O(n^3) Hungarian algorithm using potentials. It returns the optimal cost and
// the assignment mapping each row to a column. The input must contain finite
// entries only.
func hungarianMin(cost [][]float64) (float64, []int) {
	n := len(cost)
	inf := math.Inf(1)
	u := make([]float64, n+1)
	v := make([]float64, n+1)
	p := make([]int, n+1)
	way := make([]int, n+1)
	for i := 1; i <= n; i++ {
		p[0] = i
		j0 := 0
		minv := make([]float64, n+1)
		used := make([]bool, n+1)
		for j := 0; j <= n; j++ {
			minv[j] = inf
		}
		for {
			used[j0] = true
			i0 := p[j0]
			delta := inf
			j1 := -1
			for j := 1; j <= n; j++ {
				if used[j] {
					continue
				}
				cur := cost[i0-1][j-1] - u[i0] - v[j]
				if cur < minv[j] {
					minv[j] = cur
					way[j] = j0
				}
				if minv[j] < delta {
					delta = minv[j]
					j1 = j
				}
			}
			for j := 0; j <= n; j++ {
				if used[j] {
					u[p[j]] += delta
					v[j] -= delta
				} else {
					minv[j] -= delta
				}
			}
			j0 = j1
			if p[j0] == 0 {
				break
			}
		}
		for j0 != 0 {
			j1 := way[j0]
			p[j0] = p[j1]
			j0 = j1
		}
	}
	assign := make([]int, n)
	for j := 1; j <= n; j++ {
		if p[j] != 0 {
			assign[p[j]-1] = j - 1
		}
	}
	return -v[0], assign
}

// permute invokes fn on every permutation of a, restoring a before returning.
func permute(a []int, k int, fn func([]int)) {
	if k == len(a) {
		fn(a)
		return
	}
	for i := k; i < len(a); i++ {
		a[k], a[i] = a[i], a[k]
		permute(a, k+1, fn)
		a[k], a[i] = a[i], a[k]
	}
}
