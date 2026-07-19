package operatortheory

import (
	"math"
	"math/cmplx"
)

// IsHermitian reports whether the matrix equals its conjugate transpose to
// within tol. A non-positive tol selects a default tolerance.
func (m *Matrix) IsHermitian(tol float64) bool {
	if !m.IsSquare() {
		return false
	}
	tol = orDefault(tol)
	n := m.rows
	for i := 0; i < n; i++ {
		for j := i; j < n; j++ {
			if cmplx.Abs(m.data[i*n+j]-cmplx.Conj(m.data[j*n+i])) > tol {
				return false
			}
		}
	}
	return true
}

// IsSelfAdjoint is a synonym for IsHermitian.
func (m *Matrix) IsSelfAdjoint(tol float64) bool { return m.IsHermitian(tol) }

// IsSkewHermitian reports whether m^H = -m to within tol.
func (m *Matrix) IsSkewHermitian(tol float64) bool {
	if !m.IsSquare() {
		return false
	}
	tol = orDefault(tol)
	n := m.rows
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if cmplx.Abs(m.data[i*n+j]+cmplx.Conj(m.data[j*n+i])) > tol {
				return false
			}
		}
	}
	return true
}

// IsSymmetric reports whether m equals its transpose (no conjugation) to within
// tol.
func (m *Matrix) IsSymmetric(tol float64) bool {
	if !m.IsSquare() {
		return false
	}
	tol = orDefault(tol)
	n := m.rows
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			if cmplx.Abs(m.data[i*n+j]-m.data[j*n+i]) > tol {
				return false
			}
		}
	}
	return true
}

// IsNormal reports whether m commutes with its adjoint, m^H m = m m^H, to
// within tol.
func (m *Matrix) IsNormal(tol float64) bool {
	if !m.IsSquare() {
		return false
	}
	tol = orDefault(tol)
	ah := m.Adjoint()
	c, err := m.Commutator(ah)
	if err != nil {
		return false
	}
	return c.MaxAbs() <= tol*(1+m.FrobeniusNorm())
}

// IsUnitary reports whether m^H m equals the identity to within tol.
func (m *Matrix) IsUnitary(tol float64) bool {
	if !m.IsSquare() {
		return false
	}
	tol = orDefault(tol)
	ah := m.Adjoint()
	p, _ := ah.Mul(m)
	return p.Equal(Identity(m.rows), tol*(1+m.FrobeniusNorm()))
}

// IsOrthogonal reports whether a real matrix satisfies m^T m = I to within tol.
// If the matrix has non-negligible imaginary part it is not orthogonal in this
// sense.
func (m *Matrix) IsOrthogonal(tol float64) bool {
	if !m.IsSquare() {
		return false
	}
	tol = orDefault(tol)
	for _, z := range m.data {
		if math.Abs(imag(z)) > tol {
			return false
		}
	}
	t := m.Transpose()
	p, _ := t.Mul(m)
	return p.Equal(Identity(m.rows), tol*(1+m.FrobeniusNorm()))
}

// IsIsometry reports whether m preserves norms, i.e. m^H m = I. For square
// matrices this coincides with unitarity; for tall matrices it means the
// columns are orthonormal.
func (m *Matrix) IsIsometry(tol float64) bool {
	tol = orDefault(tol)
	ah := m.Adjoint()
	p, err := ah.Mul(m)
	if err != nil {
		return false
	}
	return p.Equal(Identity(m.cols), tol*(1+m.FrobeniusNorm()))
}

// IsPartialIsometry reports whether m^H m is an orthogonal projection, the
// defining property of a partial isometry.
func (m *Matrix) IsPartialIsometry(tol float64) bool {
	tol = orDefault(tol)
	ah := m.Adjoint()
	p, err := ah.Mul(m)
	if err != nil {
		return false
	}
	return p.IsOrthogonalProjection(tol)
}

// IsProjection reports whether m is idempotent, m^2 = m, to within tol.
func (m *Matrix) IsProjection(tol float64) bool {
	if !m.IsSquare() {
		return false
	}
	tol = orDefault(tol)
	sq, _ := m.Mul(m)
	return sq.Equal(m, tol*(1+m.FrobeniusNorm()))
}

// IsIdempotent is a synonym for IsProjection.
func (m *Matrix) IsIdempotent(tol float64) bool { return m.IsProjection(tol) }

// IsOrthogonalProjection reports whether m is both Hermitian and idempotent, so
// that it projects orthogonally onto its range.
func (m *Matrix) IsOrthogonalProjection(tol float64) bool {
	return m.IsHermitian(tol) && m.IsProjection(tol)
}

// IsInvolution reports whether m^2 = I to within tol.
func (m *Matrix) IsInvolution(tol float64) bool {
	if !m.IsSquare() {
		return false
	}
	tol = orDefault(tol)
	sq, _ := m.Mul(m)
	return sq.Equal(Identity(m.rows), tol*(1+m.FrobeniusNorm()))
}

// IsNilpotent reports whether some power m^k (k <= n) is the zero matrix to
// within tol.
func (m *Matrix) IsNilpotent(tol float64) bool {
	if !m.IsSquare() {
		return false
	}
	tol = orDefault(tol)
	p := m.Clone()
	for k := 1; k <= m.rows; k++ {
		if p.MaxAbs() <= tol {
			return true
		}
		p, _ = p.Mul(m)
	}
	return false
}

// IsDiagonal reports whether all off-diagonal entries vanish to within tol.
func (m *Matrix) IsDiagonal(tol float64) bool {
	tol = orDefault(tol)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			if i != j && cmplx.Abs(m.data[i*m.cols+j]) > tol {
				return false
			}
		}
	}
	return true
}

// IsUpperTriangular reports whether all entries below the main diagonal vanish
// to within tol.
func (m *Matrix) IsUpperTriangular(tol float64) bool {
	tol = orDefault(tol)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < i && j < m.cols; j++ {
			if cmplx.Abs(m.data[i*m.cols+j]) > tol {
				return false
			}
		}
	}
	return true
}

// IsLowerTriangular reports whether all entries above the main diagonal vanish
// to within tol.
func (m *Matrix) IsLowerTriangular(tol float64) bool {
	tol = orDefault(tol)
	for i := 0; i < m.rows; i++ {
		for j := i + 1; j < m.cols; j++ {
			if cmplx.Abs(m.data[i*m.cols+j]) > tol {
				return false
			}
		}
	}
	return true
}

// IsPositiveDefinite reports whether m is Hermitian with strictly positive
// eigenvalues, all exceeding tol.
func (m *Matrix) IsPositiveDefinite(tol float64) bool {
	if !m.IsHermitian(orDefault(tol)) {
		return false
	}
	vals, _ := hermitianEigenRaw(m.HermitianPart())
	for _, v := range vals {
		if v <= orDefault(tol) {
			return false
		}
	}
	return true
}

// IsPositiveSemidefinite reports whether m is Hermitian with eigenvalues no
// smaller than -tol.
func (m *Matrix) IsPositiveSemidefinite(tol float64) bool {
	if !m.IsHermitian(orDefault(tol)) {
		return false
	}
	vals, _ := hermitianEigenRaw(m.HermitianPart())
	for _, v := range vals {
		if v < -orDefault(tol) {
			return false
		}
	}
	return true
}

// IsContraction reports whether the operator norm of m is at most 1 + tol.
func (m *Matrix) IsContraction(tol float64) bool {
	return m.OperatorNorm() <= 1+orDefault(tol)
}

// IsStrictContraction reports whether the operator norm of m is strictly less
// than 1 - tol.
func (m *Matrix) IsStrictContraction(tol float64) bool {
	return m.OperatorNorm() < 1-orDefault(tol)
}

// IsStable reports whether every eigenvalue of m has strictly negative real
// part (a continuous-time stable, or Hurwitz, operator).
func (m *Matrix) IsStable(tol float64) bool {
	vals, err := m.Eigenvalues()
	if err != nil {
		return false
	}
	for _, v := range vals {
		if real(v) >= -orDefault(tol) {
			return false
		}
	}
	return true
}

// IsSchurStable reports whether every eigenvalue of m lies strictly inside the
// unit disc (a discrete-time stable operator).
func (m *Matrix) IsSchurStable(tol float64) bool {
	return m.SpectralRadius() < 1-orDefault(tol)
}

// orDefault returns tol when it is positive and the package default otherwise.
func orDefault(tol float64) float64 {
	if tol <= 0 {
		return defaultTol
	}
	return tol
}
