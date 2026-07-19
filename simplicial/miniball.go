package simplicial

import "math"

// solveLinear solves the n×n linear system A·x = b by Gaussian elimination with
// partial pivoting. It returns the solution and true, or nil and false if the
// system is singular.
func solveLinear(A [][]float64, b []float64) ([]float64, bool) {
	n := len(b)
	// build an augmented working copy
	m := make([][]float64, n)
	for i := range m {
		m[i] = make([]float64, n+1)
		copy(m[i], A[i])
		m[i][n] = b[i]
	}
	for col := 0; col < n; col++ {
		// partial pivot
		piv := col
		best := math.Abs(m[col][col])
		for r := col + 1; r < n; r++ {
			if v := math.Abs(m[r][col]); v > best {
				best = v
				piv = r
			}
		}
		if best < 1e-12 {
			return nil, false
		}
		m[col], m[piv] = m[piv], m[col]
		for r := 0; r < n; r++ {
			if r == col {
				continue
			}
			f := m[r][col] / m[col][col]
			for c := col; c <= n; c++ {
				m[r][c] -= f * m[col][c]
			}
		}
	}
	x := make([]float64, n)
	for i := 0; i < n; i++ {
		x[i] = m[i][n] / m[i][i]
	}
	return x, true
}

// Circumsphere returns the centre and radius of the unique sphere passing
// through all of the given (affinely independent) points, with the centre lying
// in their affine hull. For a single point the radius is 0; for two points it is
// the sphere with them as a diameter. It returns ok = false when the points are
// affinely dependent so that no such sphere is determined.
func Circumsphere(points [][]float64) (center []float64, radius float64, ok bool) {
	m := len(points)
	if m == 0 {
		return nil, 0, false
	}
	d := len(points[0])
	p0 := points[0]
	if m == 1 {
		return append([]float64(nil), p0...), 0, true
	}
	// v_i = points[i]-p0 for i=1..m-1; solve Gram system G λ = c with
	// G_ij = v_i·v_j and c_i = |v_i|^2/2. Center = p0 + Σ λ_j v_j.
	k := m - 1
	v := make([][]float64, k)
	for i := 0; i < k; i++ {
		vi := make([]float64, d)
		for t := 0; t < d; t++ {
			vi[t] = points[i+1][t] - p0[t]
		}
		v[i] = vi
	}
	G := make([][]float64, k)
	c := make([]float64, k)
	for i := 0; i < k; i++ {
		G[i] = make([]float64, k)
		for j := 0; j < k; j++ {
			G[i][j] = dot(v[i], v[j])
		}
		c[i] = dot(v[i], v[i]) / 2
	}
	lambda, solved := solveLinear(G, c)
	if !solved {
		return nil, 0, false
	}
	center = append([]float64(nil), p0...)
	for j := 0; j < k; j++ {
		for t := 0; t < d; t++ {
			center[t] += lambda[j] * v[j][t]
		}
	}
	radius = EuclideanDistance(center, p0)
	return center, radius, true
}

func dot(a, b []float64) float64 {
	var s float64
	for i := range a {
		s += a[i] * b[i]
	}
	return s
}

const mebEps = 1e-9

// MinimalEnclosingBall returns the centre and radius of the smallest ball
// containing every one of the given points, computed exactly with Welzl's
// algorithm. The empty input yields a zero radius and nil centre.
func MinimalEnclosingBall(points [][]float64) (center []float64, radius float64) {
	if len(points) == 0 {
		return nil, 0
	}
	dim := len(points[0])
	// work on a copy of the index order
	pts := make([][]float64, len(points))
	copy(pts, points)
	c, r := welzl(pts, len(pts), nil, dim)
	return c, r
}

// welzl computes the minimal enclosing ball of the first n points of pts with
// the points in boundary forced onto the sphere.
func welzl(pts [][]float64, n int, boundary [][]float64, dim int) ([]float64, float64) {
	if n == 0 || len(boundary) == dim+1 {
		c, r, ok := Circumsphere(boundary)
		if !ok {
			return trivialBall(boundary)
		}
		return c, r
	}
	// remove the last active point
	p := pts[n-1]
	c, r := welzl(pts, n-1, boundary, dim)
	if c != nil && inBall(p, c, r) {
		return c, r
	}
	nb := append(append([][]float64(nil), boundary...), p)
	return welzl(pts, n-1, nb, dim)
}

// trivialBall returns a small enclosing ball of up to a few boundary points when
// the circumsphere is degenerate, using the farthest pair as a fallback.
func trivialBall(boundary [][]float64) ([]float64, float64) {
	if len(boundary) == 0 {
		return nil, 0
	}
	if len(boundary) == 1 {
		return append([]float64(nil), boundary[0]...), 0
	}
	// center at centroid, radius = max distance
	d := len(boundary[0])
	c := make([]float64, d)
	for _, p := range boundary {
		for t := 0; t < d; t++ {
			c[t] += p[t]
		}
	}
	for t := range c {
		c[t] /= float64(len(boundary))
	}
	var r float64
	for _, p := range boundary {
		if dd := EuclideanDistance(c, p); dd > r {
			r = dd
		}
	}
	return c, r
}

func inBall(p, center []float64, radius float64) bool {
	return EuclideanDistance(p, center) <= radius+mebEps
}

// EnclosingRadius returns the radius of the minimal enclosing ball of the
// points, a convenience wrapper around [MinimalEnclosingBall].
func EnclosingRadius(points [][]float64) float64 {
	_, r := MinimalEnclosingBall(points)
	return r
}
