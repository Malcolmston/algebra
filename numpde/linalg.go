package numpde

// Linspace returns n evenly spaced samples over the closed interval [a, b].
// The first sample equals a and, when n >= 2, the last equals b. It panics if
// n < 1.
func Linspace(a, b float64, n int) []float64 {
	if n < 1 {
		panic("numpde: Linspace requires n >= 1")
	}
	xs := make([]float64, n)
	if n == 1 {
		xs[0] = a
		return xs
	}
	h := (b - a) / float64(n-1)
	for i := range xs {
		xs[i] = a + float64(i)*h
	}
	xs[n-1] = b
	return xs
}

// Clone returns a copy of the slice u. The result never aliases u.
func Clone(u []float64) []float64 {
	c := make([]float64, len(u))
	copy(c, u)
	return c
}

// Clone2D returns a deep copy of the matrix u. Neither the outer slice nor any
// row of the result aliases u.
func Clone2D(u [][]float64) [][]float64 {
	c := make([][]float64, len(u))
	for i := range u {
		c[i] = make([]float64, len(u[i]))
		copy(c[i], u[i])
	}
	return c
}

// Zeros2D allocates an nx by ny matrix with every entry set to zero.
func Zeros2D(nx, ny int) [][]float64 {
	m := make([][]float64, nx)
	for i := range m {
		m[i] = make([]float64, ny)
	}
	return m
}

// MaxAbs returns the largest absolute value among the elements of u, or zero
// when u is empty.
func MaxAbs(u []float64) float64 {
	m := 0.0
	for _, v := range u {
		if a := numpdeAbs(v); a > m {
			m = a
		}
	}
	return m
}

// L2Norm returns the Euclidean (2-norm) length of u, sqrt(sum u_i^2).
func L2Norm(u []float64) float64 {
	s := 0.0
	for _, v := range u {
		s += v * v
	}
	return numpdeSqrt(s)
}

// RMSNorm returns the root-mean-square norm of u, sqrt(mean(u_i^2)). It is the
// discrete L2 norm normalised by the number of samples and is a good measure of
// grid error that is insensitive to mesh refinement. It returns zero for empty
// input.
func RMSNorm(u []float64) float64 {
	if len(u) == 0 {
		return 0
	}
	s := 0.0
	for _, v := range u {
		s += v * v
	}
	return numpdeSqrt(s / float64(len(u)))
}

// LInfNorm returns the maximum-norm (infinity norm) of u, the largest absolute
// element. It is an alias-free synonym for MaxAbs used where norm terminology
// reads more clearly.
func LInfNorm(u []float64) float64 { return MaxAbs(u) }

// MaxAbsDiff returns the infinity norm of the elementwise difference a-b. It
// panics if the slices have different lengths.
func MaxAbsDiff(a, b []float64) float64 {
	if len(a) != len(b) {
		panic("numpde: MaxAbsDiff length mismatch")
	}
	m := 0.0
	for i := range a {
		if d := numpdeAbs(a[i] - b[i]); d > m {
			m = d
		}
	}
	return m
}

// MaxAbsDiff2D returns the infinity norm of the elementwise difference of two
// equally shaped matrices. It panics on a shape mismatch.
func MaxAbsDiff2D(a, b [][]float64) float64 {
	if len(a) != len(b) {
		panic("numpde: MaxAbsDiff2D shape mismatch")
	}
	m := 0.0
	for i := range a {
		if len(a[i]) != len(b[i]) {
			panic("numpde: MaxAbsDiff2D shape mismatch")
		}
		for j := range a[i] {
			if d := numpdeAbs(a[i][j] - b[i][j]); d > m {
				m = d
			}
		}
	}
	return m
}

// ThomasSolve solves the tridiagonal linear system M x = d, where M has
// sub-diagonal a, main diagonal b and super-diagonal c. Element a[0] and
// c[n-1] are ignored. All slices must have equal length n >= 1. The inputs are
// left unmodified; a freshly allocated solution slice is returned. It panics on
// a length mismatch and returns a slice of NaNs is avoided by requiring
// non-zero pivots (a zero pivot triggers a panic, indicating a singular or
// non-diagonally-dominant system).
func ThomasSolve(a, b, c, d []float64) []float64 {
	n := len(b)
	if len(a) != n || len(c) != n || len(d) != n {
		panic("numpde: ThomasSolve length mismatch")
	}
	if n == 0 {
		return []float64{}
	}
	cp := make([]float64, n) // modified super-diagonal
	dp := make([]float64, n) // modified right-hand side
	if b[0] == 0 {
		panic("numpde: ThomasSolve zero pivot")
	}
	cp[0] = c[0] / b[0]
	dp[0] = d[0] / b[0]
	for i := 1; i < n; i++ {
		m := b[i] - a[i]*cp[i-1]
		if m == 0 {
			panic("numpde: ThomasSolve zero pivot")
		}
		cp[i] = c[i] / m
		dp[i] = (d[i] - a[i]*dp[i-1]) / m
	}
	x := make([]float64, n)
	x[n-1] = dp[n-1]
	for i := n - 2; i >= 0; i-- {
		x[i] = dp[i] - cp[i]*x[i+1]
	}
	return x
}

// SecondDerivativeMatrix returns the dense n by n central-difference
// approximation of the second-derivative operator d^2/dx^2 with grid spacing
// dx and homogeneous Dirichlet treatment at the ends (the first and last rows
// use the standard interior stencil truncated at the boundary). Interior rows
// contain [1, -2, 1]/dx^2. It is primarily useful for small problems, teaching,
// and assembling method-of-lines Jacobians.
func SecondDerivativeMatrix(n int, dx float64) [][]float64 {
	m := Zeros2D(n, n)
	inv := 1.0 / (dx * dx)
	for i := 0; i < n; i++ {
		m[i][i] = -2 * inv
		if i-1 >= 0 {
			m[i][i-1] = inv
		}
		if i+1 < n {
			m[i][i+1] = inv
		}
	}
	return m
}

// FirstDerivativeMatrix returns the dense n by n central-difference
// approximation of d/dx with spacing dx. Interior rows contain
// [-1, 0, 1]/(2*dx); the first row uses a forward and the last row a backward
// one-sided difference so the operator remains first-order accurate at the
// boundaries.
func FirstDerivativeMatrix(n int, dx float64) [][]float64 {
	m := Zeros2D(n, n)
	c := 1.0 / (2 * dx)
	for i := 1; i < n-1; i++ {
		m[i][i-1] = -c
		m[i][i+1] = c
	}
	if n >= 2 {
		// forward difference at i=0, backward at i=n-1
		f := 1.0 / dx
		m[0][0] = -f
		m[0][1] = f
		m[n-1][n-2] = -f
		m[n-1][n-1] = f
	}
	return m
}

// MatVec multiplies the dense matrix m by the vector x and returns m*x. It
// panics when the inner dimensions disagree.
func MatVec(m [][]float64, x []float64) []float64 {
	y := make([]float64, len(m))
	for i := range m {
		if len(m[i]) != len(x) {
			panic("numpde: MatVec dimension mismatch")
		}
		s := 0.0
		for j := range m[i] {
			s += m[i][j] * x[j]
		}
		y[i] = s
	}
	return y
}
