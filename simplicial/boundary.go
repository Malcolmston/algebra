package simplicial

// indexOf builds a map from simplex key to position for a sorted slice.
func indexOf(ss []Simplex) map[string]int {
	m := make(map[string]int, len(ss))
	for i, s := range ss {
		m[s.Key()] = i
	}
	return m
}

// BoundaryMatrixGF2 returns the matrix of the boundary map ∂_k : C_k → C_{k-1}
// over GF(2). Columns are indexed by the k-simplices and rows by the
// (k−1)-simplices, both in [CompareSimplices] order. For k ≤ 0 the map is the
// zero map into the trivial group and a matrix with zero rows is returned.
func (c *Complex) BoundaryMatrixGF2(k int) *GF2Matrix {
	cols := c.SimplicesOfDim(k)
	if k <= 0 {
		return NewGF2Matrix(0, len(cols))
	}
	rows := c.SimplicesOfDim(k - 1)
	idx := indexOf(rows)
	m := NewGF2Matrix(len(rows), len(cols))
	for j, s := range cols {
		for _, f := range s.Faces() {
			if i, ok := idx[f.Key()]; ok {
				m.data[i][j] = 1
			}
		}
	}
	return m
}

// BoundaryMatrixZ returns the matrix of the boundary map ∂_k : C_k → C_{k-1}
// over the integers, with the standard orientation signs. Columns are indexed
// by k-simplices and rows by (k−1)-simplices in [CompareSimplices] order. For
// k ≤ 0 a matrix with zero rows is returned.
func (c *Complex) BoundaryMatrixZ(k int) *IntMatrix {
	cols := c.SimplicesOfDim(k)
	if k <= 0 {
		return NewIntMatrix(0, len(cols))
	}
	rows := c.SimplicesOfDim(k - 1)
	idx := indexOf(rows)
	m := NewIntMatrix(len(rows), len(cols))
	for j, s := range cols {
		for _, term := range s.Boundary() {
			if i, ok := idx[term.Face.Key()]; ok {
				m.data[i][j].SetInt64(int64(term.Sign))
			}
		}
	}
	return m
}

// BoundaryMatrixQ returns the boundary map ∂_k over the rationals, identical in
// content to [Complex.BoundaryMatrixZ] but with rational entries.
func (c *Complex) BoundaryMatrixQ(k int) *RatMatrix {
	return c.BoundaryMatrixZ(k).ToRat()
}

// ChainRank returns the dimension of the chain group C_k, i.e. the number of
// k-simplices in the complex.
func (c *Complex) ChainRank(k int) int { return c.NumSimplicesOfDim(k) }

// BoundaryRankGF2 returns the rank over GF(2) of the boundary map ∂_k.
func (c *Complex) BoundaryRankGF2(k int) int {
	if k <= 0 {
		return 0
	}
	return c.BoundaryMatrixGF2(k).Rank()
}

// BoundaryRankQ returns the rank over Q of the boundary map ∂_k.
func (c *Complex) BoundaryRankQ(k int) int {
	if k <= 0 {
		return 0
	}
	return c.BoundaryMatrixZ(k).Rank()
}

// IsCycleGF2 reports whether the given GF(2) chain — a 0/1 coefficient vector
// indexed like [Complex.SimplicesOfDim](k) — is a cycle, i.e. lies in the
// kernel of ∂_k.
func (c *Complex) IsCycleGF2(k int, chain []int) bool {
	if k <= 0 {
		return true
	}
	b := c.BoundaryMatrixGF2(k)
	if len(chain) != b.cols {
		return false
	}
	for i := 0; i < b.rows; i++ {
		s := 0
		for j := 0; j < b.cols; j++ {
			s ^= int(b.data[i][j]) & (chain[j] & 1)
		}
		if s != 0 {
			return false
		}
	}
	return true
}
