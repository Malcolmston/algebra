package optimalcontrol

import "math"

// Rank returns the numerical rank of a matrix using Gaussian elimination with
// partial pivoting and the supplied tolerance for treating pivots as zero.
func Rank(a *Matrix, tol float64) int {
	m := a.Clone()
	rows, cols := m.rows, m.cols
	rank := 0
	for col := 0; col < cols && rank < rows; col++ {
		// Find pivot in column col at or below row rank.
		piv := rank
		max := math.Abs(m.At(rank, col))
		for i := rank + 1; i < rows; i++ {
			if v := math.Abs(m.At(i, col)); v > max {
				max = v
				piv = i
			}
		}
		if max <= tol {
			continue
		}
		// Swap rows.
		for j := 0; j < cols; j++ {
			t := m.At(rank, j)
			m.Set(rank, j, m.At(piv, j))
			m.Set(piv, j, t)
		}
		// Eliminate below.
		for i := 0; i < rows; i++ {
			if i == rank {
				continue
			}
			f := m.At(i, col) / m.At(rank, col)
			for j := col; j < cols; j++ {
				m.Set(i, j, m.At(i, j)-f*m.At(rank, j))
			}
		}
		rank++
	}
	return rank
}

// ControllabilityMatrix returns the controllability matrix
// [B, AB, A²B, …, Aⁿ⁻¹B] for the pair (A, B).
func ControllabilityMatrix(a, b *Matrix) *Matrix {
	n := a.rows
	c := b.Clone()
	term := b.Clone()
	for k := 1; k < n; k++ {
		term = a.Mul(term)
		c = HStack(c, term)
	}
	return c
}

// ObservabilityMatrix returns the observability matrix
// [C; CA; CA²; …; CAⁿ⁻¹] for the pair (A, C).
func ObservabilityMatrix(a, c *Matrix) *Matrix {
	n := a.rows
	o := c.Clone()
	term := c.Clone()
	for k := 1; k < n; k++ {
		term = term.Mul(a)
		o = VStack(o, term)
	}
	return o
}

// IsControllable reports whether the pair (A, B) is controllable, i.e. the
// controllability matrix has full row rank n.
func IsControllable(a, b *Matrix) bool {
	return Rank(ControllabilityMatrix(a, b), 1e-9) == a.rows
}

// IsObservable reports whether the pair (A, C) is observable, i.e. the
// observability matrix has full column rank n.
func IsObservable(a, c *Matrix) bool {
	return Rank(ObservabilityMatrix(a, c), 1e-9) == a.rows
}

// ControllabilityRank returns the rank of the controllability matrix of (A, B).
func ControllabilityRank(a, b *Matrix) int {
	return Rank(ControllabilityMatrix(a, b), 1e-9)
}

// ObservabilityRank returns the rank of the observability matrix of (A, C).
func ObservabilityRank(a, c *Matrix) int {
	return Rank(ObservabilityMatrix(a, c), 1e-9)
}

// IsStabilizableContinuous reports whether the pair (A, B) is stabilizable in
// continuous time: every eigenvalue of A with non-negative real part is
// controllable (Popov–Belevitch–Hautus test).
func IsStabilizableContinuous(a, b *Matrix) bool {
	for _, lam := range Eigenvalues(a) {
		if real(lam) >= 0 {
			if !pbhControllable(a, b, lam) {
				return false
			}
		}
	}
	return true
}

// IsStabilizableDiscrete reports whether the pair (A, B) is stabilizable in
// discrete time: every eigenvalue of A with modulus >= 1 is controllable.
func IsStabilizableDiscrete(a, b *Matrix) bool {
	for _, lam := range Eigenvalues(a) {
		if real(lam)*real(lam)+imag(lam)*imag(lam) >= 1 {
			if !pbhControllable(a, b, lam) {
				return false
			}
		}
	}
	return true
}

// IsDetectableContinuous reports whether the pair (A, C) is detectable in
// continuous time (dual of stabilizability).
func IsDetectableContinuous(a, c *Matrix) bool {
	return IsStabilizableContinuous(a.Transpose(), c.Transpose())
}

// IsDetectableDiscrete reports whether the pair (A, C) is detectable in
// discrete time.
func IsDetectableDiscrete(a, c *Matrix) bool {
	return IsStabilizableDiscrete(a.Transpose(), c.Transpose())
}

// pbhControllable performs the Hautus rank test [A-λI, B] for a single (complex)
// eigenvalue by working with the real and imaginary parts stacked into a real
// matrix.
func pbhControllable(a, b *Matrix, lam complex128) bool {
	n := a.rows
	// Build real 2n × (2n + 2m) matrix representing [A-λI | B] over complexes.
	re := real(lam)
	im := imag(lam)
	m := b.cols
	// Real block form of (A - λI): [[A-reI, imI],[-imI, A-reI]] and B as
	// [[B,0],[0,B]].
	big := Zeros(2*n, 2*n+2*m)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			aij := a.At(i, j)
			big.Set(i, j, aij)
			big.Set(n+i, n+j, aij)
		}
		big.Add(i, i, -re)
		big.Add(n+i, n+i, -re)
		big.Set(i, n+i, im)
		big.Set(n+i, i, -im)
	}
	for i := 0; i < n; i++ {
		for j := 0; j < m; j++ {
			bij := b.At(i, j)
			big.Set(i, 2*n+j, bij)
			big.Set(n+i, 2*n+m+j, bij)
		}
	}
	return Rank(big, 1e-8) == 2*n
}
