package tropical

import "math"

// MaxCycleMean returns the maximum cycle mean of a square matrix interpreted as
// a weighted digraph (entry (i,j) is the weight of the edge i->j; the tropical
// zero means no edge). It is computed with Karp's algorithm and equals the
// max-plus eigenvalue of an irreducible matrix. The boolean result is false
// when the graph contains no cycle at all. It returns ErrNotSquare for a
// non-square matrix.
func (m Matrix) MaxCycleMean() (float64, bool, error) {
	return m.cycleMean(true)
}

// MinCycleMean returns the minimum cycle mean of a square matrix interpreted as
// a weighted digraph, computed with Karp's algorithm; it equals the min-plus
// eigenvalue of an irreducible matrix. The boolean result is false when the
// graph contains no cycle. It returns ErrNotSquare for a non-square matrix.
func (m Matrix) MinCycleMean() (float64, bool, error) {
	return m.cycleMean(false)
}

// cycleMean implements Karp's algorithm. With wantMax it computes the maximum
// cycle mean; otherwise the minimum cycle mean. The "start anywhere"
// initialisation D[0][v]=0 makes the result range over every cycle of the
// graph rather than only those reachable from a single source.
func (m Matrix) cycleMean(wantMax bool) (float64, bool, error) {
	if !m.IsSquare() {
		return 0, false, ErrNotSquare
	}
	n := m.rows
	if n == 0 {
		return 0, false, nil
	}
	neg := math.Inf(-1)
	pos := math.Inf(1)
	// D[k][v] = best weight of a walk of exactly k edges ending at v.
	d := make([][]float64, n+1)
	for k := range d {
		d[k] = make([]float64, n)
	}
	for v := 0; v < n; v++ {
		d[0][v] = 0
	}
	for k := 1; k <= n; k++ {
		for v := 0; v < n; v++ {
			best := neg
			if !wantMax {
				best = pos
			}
			for u := 0; u < n; u++ {
				w := m.data[u][v] // edge u->v
				if math.IsInf(w, 0) && (m.sr.IsZero(w)) {
					continue
				}
				prev := d[k-1][u]
				if math.IsInf(prev, 0) {
					continue
				}
				cand := prev + w
				if wantMax {
					if cand > best {
						best = cand
					}
				} else {
					if cand < best {
						best = cand
					}
				}
			}
			d[k][v] = best
		}
	}
	found := false
	var lambda float64
	if wantMax {
		lambda = neg
	} else {
		lambda = pos
	}
	for v := 0; v < n; v++ {
		if math.IsInf(d[n][v], 0) {
			continue
		}
		// inner = best over cycle candidates for this v
		var inner float64
		if wantMax {
			inner = pos // we take min over k of the ratios, then max over v
		} else {
			inner = neg
		}
		valid := false
		for k := 0; k < n; k++ {
			if math.IsInf(d[k][v], 0) {
				continue
			}
			ratio := (d[n][v] - d[k][v]) / float64(n-k)
			if wantMax {
				if ratio < inner {
					inner = ratio
				}
			} else {
				if ratio > inner {
					inner = ratio
				}
			}
			valid = true
		}
		if !valid {
			continue
		}
		if wantMax {
			if inner > lambda {
				lambda = inner
			}
		} else {
			if inner < lambda {
				lambda = inner
			}
		}
		found = true
	}
	if !found {
		return 0, false, nil
	}
	return lambda, true, nil
}

// Eigenvalue returns the tropical eigenvalue of an irreducible square matrix:
// the maximum cycle mean for max-plus and the minimum cycle mean for min-plus.
// The boolean result is false when the graph has no cycle. It returns
// ErrNotSquare for a non-square matrix.
func (m Matrix) Eigenvalue() (float64, bool, error) {
	if m.sr.IsMaxPlus() {
		return m.MaxCycleMean()
	}
	return m.MinCycleMean()
}

// Eigenvector returns a tropical eigenvector v of an irreducible square matrix
// together with its eigenvalue lambda, so that A (*) v = lambda (*) v. The
// vector is obtained as a critical column of the normalised matrix closure. The
// boolean result is false when the matrix has no cycle (and hence no finite
// eigenvalue). It returns ErrNotSquare for a non-square matrix.
func (m Matrix) Eigenvector() (Vector, float64, bool, error) {
	if !m.IsSquare() {
		return Vector{}, 0, false, ErrNotSquare
	}
	lambda, ok, err := m.Eigenvalue()
	if err != nil {
		return Vector{}, 0, false, err
	}
	if !ok {
		return Vector{}, 0, false, nil
	}
	n := m.rows
	// Normalise: B = A / lambda so that its critical cycles have mean 0.
	b := m.ScalarMul(-lambda)
	// Bounded closure B* = I (+) B (+) ... (+) B^{n-1}; finite since the
	// normalised matrix has no strictly better-than-zero cycle.
	star, _ := b.StarSeries(n)
	// A critical node c lies on a zero-mean cycle: (B*)_{cc} == 0.
	c := -1
	for i := 0; i < n; i++ {
		if math.Abs(star.data[i][i]) < 1e-9 {
			c = i
			break
		}
	}
	if c == -1 {
		// Fall back to the node whose diagonal closure is best.
		best := m.sr.Zero()
		for i := 0; i < n; i++ {
			if m.sr.AtLeastAsGood(star.data[i][i], best) {
				best = star.data[i][i]
				c = i
			}
		}
	}
	v := star.ColVector(c)
	return v, lambda, true, nil
}

// CriticalNodes returns the indices of the nodes that lie on a cycle attaining
// the tropical eigenvalue (the critical cycles). The boolean result is false
// when the matrix has no cycle. It returns ErrNotSquare for a non-square
// matrix.
func (m Matrix) CriticalNodes() ([]int, bool, error) {
	if !m.IsSquare() {
		return nil, false, ErrNotSquare
	}
	lambda, ok, err := m.Eigenvalue()
	if err != nil || !ok {
		return nil, false, err
	}
	n := m.rows
	b := m.ScalarMul(-lambda)
	star, _ := b.StarSeries(n)
	var nodes []int
	for i := 0; i < n; i++ {
		if math.Abs(star.data[i][i]) < 1e-9 {
			nodes = append(nodes, i)
		}
	}
	return nodes, len(nodes) > 0, nil
}
