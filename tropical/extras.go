package tropical

import "math"

// MinPlusIdentity returns the n-by-n min-plus identity matrix.
func MinPlusIdentity(n int) Matrix { return Identity(MinPlusSemiring(), n) }

// MaxPlusIdentity returns the n-by-n max-plus identity matrix.
func MaxPlusIdentity(n int) Matrix { return Identity(MaxPlusSemiring(), n) }

// MinPlusZeros returns an r-by-c min-plus matrix of tropical zeros (+Inf).
func MinPlusZeros(r, c int) Matrix { return Zeros(MinPlusSemiring(), r, c) }

// MaxPlusZeros returns an r-by-c max-plus matrix of tropical zeros (-Inf).
func MaxPlusZeros(r, c int) Matrix { return Zeros(MaxPlusSemiring(), r, c) }

// EqualScalar reports whether a and b agree to within tol, with infinities
// required to match exactly and in sign.
func (s Semiring) EqualScalar(a, b, tol float64) bool { return closeScalar(a, b, tol) }

// Map returns a new matrix obtained by applying f to every entry.
func (m Matrix) Map(f func(float64) float64) Matrix {
	out := Zeros(m.sr, m.rows, m.cols)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			out.data[i][j] = f(m.data[i][j])
		}
	}
	return out
}

// ElementwiseMul returns the elementwise (Hadamard) tropical product of m and
// n, whose entries are the ordinary sums of the corresponding entries. It
// returns an error if the shapes or semirings differ.
func (m Matrix) ElementwiseMul(n Matrix) (Matrix, error) {
	if m.sr.kind != n.sr.kind || m.rows != n.rows || m.cols != n.cols {
		return Matrix{}, ErrDim
	}
	out := Zeros(m.sr, m.rows, m.cols)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			out.data[i][j] = m.sr.Mul(m.data[i][j], n.data[i][j])
		}
	}
	return out, nil
}

// Dual returns the same matrix reinterpreted in the dual semiring, with every
// entry negated. Negating the entries and switching min-plus with max-plus
// turns any tropical problem into the equivalent dual problem.
func (m Matrix) Dual() Matrix {
	dual := m.sr.Dual()
	out := Zeros(dual, m.rows, m.cols)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			out.data[i][j] = -m.data[i][j]
		}
	}
	return out
}

// MaxEntry returns the largest entry in the usual numeric order, or -Inf for an
// empty matrix.
func (m Matrix) MaxEntry() float64 {
	best := math.Inf(-1)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			if m.data[i][j] > best {
				best = m.data[i][j]
			}
		}
	}
	return best
}

// MinEntry returns the smallest entry in the usual numeric order, or +Inf for
// an empty matrix.
func (m Matrix) MinEntry() float64 {
	best := math.Inf(1)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			if m.data[i][j] < best {
				best = m.data[i][j]
			}
		}
	}
	return best
}

// IsFinite reports whether every entry is finite (no tropical zeros or other
// infinities).
func (m Matrix) IsFinite() bool {
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			if math.IsInf(m.data[i][j], 0) {
				return false
			}
		}
	}
	return true
}

// CountFinite returns the number of finite entries.
func (m Matrix) CountFinite() int {
	c := 0
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			if !math.IsInf(m.data[i][j], 0) {
				c++
			}
		}
	}
	return c
}

// IsSymmetric reports whether the matrix is square and equal to its transpose.
func (m Matrix) IsSymmetric() bool {
	if !m.IsSquare() {
		return false
	}
	for i := 0; i < m.rows; i++ {
		for j := i + 1; j < m.cols; j++ {
			if m.data[i][j] != m.data[j][i] {
				return false
			}
		}
	}
	return true
}

// KroneckerProduct returns the tropical Kronecker product of m and n, the
// (m.rows*n.rows)-by-(m.cols*n.cols) matrix whose block (i,j) is m[i][j]
// tropically scaling n. It returns an error if the semirings differ.
func (m Matrix) KroneckerProduct(n Matrix) (Matrix, error) {
	if m.sr.kind != n.sr.kind {
		return Matrix{}, ErrDim
	}
	out := Zeros(m.sr, m.rows*n.rows, m.cols*n.cols)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			a := m.data[i][j]
			for k := 0; k < n.rows; k++ {
				for l := 0; l < n.cols; l++ {
					out.data[i*n.rows+k][j*n.cols+l] = m.sr.Mul(a, n.data[k][l])
				}
			}
		}
	}
	return out, nil
}

// DirectSum returns the block-diagonal matrix with m in the top-left and n in
// the bottom-right, the off-diagonal blocks filled with tropical zeros. It
// returns an error if the semirings differ.
func (m Matrix) DirectSum(n Matrix) (Matrix, error) {
	if m.sr.kind != n.sr.kind {
		return Matrix{}, ErrDim
	}
	out := Zeros(m.sr, m.rows+n.rows, m.cols+n.cols)
	for i := 0; i < m.rows; i++ {
		copy(out.data[i][:m.cols], m.data[i])
	}
	for i := 0; i < n.rows; i++ {
		copy(out.data[m.rows+i][m.cols:], n.data[i])
	}
	return out, nil
}

// HStack returns the horizontal concatenation [m | n]. It returns an error if
// the row counts or semirings differ.
func (m Matrix) HStack(n Matrix) (Matrix, error) {
	if m.sr.kind != n.sr.kind || m.rows != n.rows {
		return Matrix{}, ErrDim
	}
	out := Zeros(m.sr, m.rows, m.cols+n.cols)
	for i := 0; i < m.rows; i++ {
		copy(out.data[i][:m.cols], m.data[i])
		copy(out.data[i][m.cols:], n.data[i])
	}
	return out, nil
}

// VStack returns the vertical concatenation of m stacked above n. It returns an
// error if the column counts or semirings differ.
func (m Matrix) VStack(n Matrix) (Matrix, error) {
	if m.sr.kind != n.sr.kind || m.cols != n.cols {
		return Matrix{}, ErrDim
	}
	out := Zeros(m.sr, m.rows+n.rows, m.cols)
	for i := 0; i < m.rows; i++ {
		copy(out.data[i], m.data[i])
	}
	for i := 0; i < n.rows; i++ {
		copy(out.data[m.rows+i], n.data[i])
	}
	return out, nil
}

// PathWeight returns the tropical product of the edge weights along the walk
// visiting the given node sequence, that is the sum of m[nodes[k]][nodes[k+1]].
// An empty or single-node walk has weight the tropical one (0). It panics if a
// node index is out of range.
func (m Matrix) PathWeight(nodes []int) float64 {
	w := m.sr.One()
	for k := 0; k+1 < len(nodes); k++ {
		w = m.sr.Mul(w, m.data[nodes[k]][nodes[k+1]])
	}
	return w
}

// AssignmentValue returns the tropical product of the entries m[i][perm[i]] for
// a square matrix, the weight of the permutation as an assignment. It panics if
// perm is not a length-n slice for an n-by-n matrix.
func (m Matrix) AssignmentValue(perm []int) float64 {
	if !m.IsSquare() || len(perm) != m.rows {
		panic("tropical: AssignmentValue requires a permutation of the row indices")
	}
	v := m.sr.One()
	for i := 0; i < m.rows; i++ {
		v = m.sr.Mul(v, m.data[i][perm[i]])
	}
	return v
}

// Cofactor returns the tropical cofactor of entry (i,j): the tropical permanent
// of the minor obtained by deleting row i and column j. It returns
// ErrNotSquare for a non-square matrix.
func (m Matrix) Cofactor(i, j int) (float64, error) {
	if !m.IsSquare() {
		return 0, ErrNotSquare
	}
	return m.Minor(i, j).Permanent()
}

// SingleSourceShortestPaths returns the vector of best walk weights from src to
// every node, computed with the Bellman-Ford relaxation: shortest paths for
// min-plus and longest paths for max-plus. It returns ErrDivergent when a bad
// cycle is reachable and ErrNotSquare for a non-square matrix.
func (m Matrix) SingleSourceShortestPaths(src int) (Vector, error) {
	if !m.IsSquare() {
		return Vector{}, ErrNotSquare
	}
	n := m.rows
	if src < 0 || src >= n {
		return Vector{}, ErrDim
	}
	dist := make([]float64, n)
	z := m.sr.Zero()
	for i := range dist {
		dist[i] = z
	}
	dist[src] = 0
	relax := func() bool {
		changed := false
		for u := 0; u < n; u++ {
			if m.sr.IsZero(dist[u]) {
				continue
			}
			for v := 0; v < n; v++ {
				w := m.data[u][v]
				if m.sr.IsZero(w) {
					continue
				}
				nd := m.sr.Mul(dist[u], w)
				improved := m.sr.Add(dist[v], nd)
				if improved != dist[v] {
					dist[v] = improved
					changed = true
				}
			}
		}
		return changed
	}
	for i := 0; i < n-1; i++ {
		if !relax() {
			break
		}
	}
	if relax() {
		return Vector{}, ErrDivergent
	}
	return Vector{data: dist, sr: m.sr}, nil
}

// Map returns a new vector obtained by applying f to every entry.
func (v Vector) Map(f func(float64) float64) Vector {
	out := make([]float64, len(v.data))
	for i, x := range v.data {
		out[i] = f(x)
	}
	return Vector{data: out, sr: v.sr}
}

// Dual returns the same vector reinterpreted in the dual semiring, with every
// entry negated.
func (v Vector) Dual() Vector {
	out := make([]float64, len(v.data))
	for i, x := range v.data {
		out[i] = -x
	}
	return Vector{data: out, sr: v.sr.Dual()}
}

// Argmin returns the index of the numerically smallest entry, or -1 for an
// empty vector.
func (v Vector) Argmin() int {
	idx := -1
	best := math.Inf(1)
	for i, x := range v.data {
		if idx == -1 || x < best {
			best = x
			idx = i
		}
	}
	return idx
}

// Argmax returns the index of the numerically largest entry, or -1 for an empty
// vector.
func (v Vector) Argmax() int {
	idx := -1
	best := math.Inf(-1)
	for i, x := range v.data {
		if idx == -1 || x > best {
			best = x
			idx = i
		}
	}
	return idx
}

// LeadingCoeff returns the coefficient of the highest power of the polynomial.
func (p Poly) LeadingCoeff() float64 { return p.coeffs[len(p.coeffs)-1] }

// ConstantCoeff returns the coefficient of x^0.
func (p Poly) ConstantCoeff() float64 { return p.coeffs[0] }

// IsMonic reports whether the leading coefficient equals the tropical one (0).
func (p Poly) IsMonic() bool { return p.coeffs[len(p.coeffs)-1] == 0 }

// ShiftUp returns the polynomial multiplied by x^k, prepending k tropical-zero
// coefficients. It panics for negative k.
func (p Poly) ShiftUp(k int) Poly {
	if k < 0 {
		panic("tropical: ShiftUp requires k >= 0")
	}
	out := make([]float64, len(p.coeffs)+k)
	z := p.sr.Zero()
	for i := 0; i < k; i++ {
		out[i] = z
	}
	copy(out[k:], p.coeffs)
	return NewPoly(p.sr, out)
}

// Dual returns the polynomial reinterpreted in the dual semiring with every
// coefficient negated; its roots are the negatives of the original roots.
func (p Poly) Dual() Poly {
	out := make([]float64, len(p.coeffs))
	for i, a := range p.coeffs {
		out[i] = -a
	}
	return NewPoly(p.sr.Dual(), out)
}
