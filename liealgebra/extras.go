package liealgebra

import (
	"fmt"
	"math"
)

// CTraceForm returns the complex trace form tr(AB) of two complex matrices.
func CTraceForm(a, b *CMatrix) (complex128, error) {
	ab, err := a.Mul(b)
	if err != nil {
		return 0, err
	}
	return ab.Trace()
}

// Kronecker returns the Kronecker (tensor) product m ⊗ b of complex matrices.
func (m *CMatrix) Kronecker(b *CMatrix) *CMatrix {
	out := NewCMatrix(m.Rows*b.Rows, m.Cols*b.Cols)
	oc := out.Cols
	for i := 0; i < m.Rows; i++ {
		for j := 0; j < m.Cols; j++ {
			a := m.Data[i*m.Cols+j]
			for p := 0; p < b.Rows; p++ {
				for q := 0; q < b.Cols; q++ {
					out.Data[(i*b.Rows+p)*oc+(j*b.Cols+q)] = a * b.Data[p*b.Cols+q]
				}
			}
		}
	}
	return out
}

// HermitianPart returns (m+mᴴ)/2, the Hermitian part of a square complex matrix.
func (m *CMatrix) HermitianPart() (*CMatrix, error) {
	if !m.IsSquare() {
		return nil, ErrNotSquare
	}
	s, _ := m.Add(m.Dagger())
	return s.Scale(complex(0.5, 0)), nil
}

// AntiHermitianPart returns (m-mᴴ)/2, the anti-Hermitian part.
func (m *CMatrix) AntiHermitianPart() (*CMatrix, error) {
	if !m.IsSquare() {
		return nil, ErrNotSquare
	}
	s, _ := m.Sub(m.Dagger())
	return s.Scale(complex(0.5, 0)), nil
}

// IsAbelian reports whether every pair of basis elements commutes to within
// tolerance tol, i.e. whether the spanned Lie algebra is abelian.
func IsAbelian(basis []*Matrix, tol float64) bool {
	for i := 0; i < len(basis); i++ {
		for j := i + 1; j < len(basis); j++ {
			br, err := Bracket(basis[i], basis[j])
			if err != nil {
				return false
			}
			if br.MaxAbs() > tol {
				return false
			}
		}
	}
	return true
}

// KillingFormRank returns the rank of the Killing form of a basis. It equals the
// dimension of the algebra for a semisimple Lie algebra and is smaller when the
// radical is nonzero.
func KillingFormRank(basis []*Matrix, tol float64) (int, error) {
	k, err := KillingForm(basis)
	if err != nil {
		return 0, err
	}
	return Rank(k, tol), nil
}

// IsNilpotentMatrix reports whether the square matrix m is nilpotent, i.e. some
// power m^k (k up to its dimension) vanishes to within tolerance tol.
func IsNilpotentMatrix(m *Matrix, tol float64) (bool, error) {
	if !m.IsSquare() {
		return false, ErrNotSquare
	}
	p := m.Clone()
	for k := 1; k <= m.Rows; k++ {
		if p.MaxAbs() <= tol {
			return true, nil
		}
		var err error
		p, err = p.Mul(m)
		if err != nil {
			return false, err
		}
	}
	return p.MaxAbs() <= tol, nil
}

// SpectralRadiusBound returns an upper bound on the spectral radius of a real
// matrix via the maximum absolute row sum (an induced ∞-norm).
func SpectralRadiusBound(m *Matrix) float64 {
	max := 0.0
	for i := 0; i < m.Rows; i++ {
		s := 0.0
		for j := 0; j < m.Cols; j++ {
			s += math.Abs(m.Data[i*m.Cols+j])
		}
		if s > max {
			max = s
		}
	}
	return max
}

// CasimirTraceNormalization returns the constant c such that the Killing form
// equals c times the trace form on the defining representation, given a basis;
// it is computed as the ratio of the (0,0) entries of the two Gram matrices and
// returns [ErrRange] when the trace form entry is zero.
func CasimirTraceNormalization(basis []*Matrix) (float64, error) {
	if len(basis) == 0 {
		return 0, ErrRange
	}
	k, err := KillingForm(basis)
	if err != nil {
		return 0, err
	}
	tf, err := TraceForm(basis[0], basis[0])
	if err != nil {
		return 0, err
	}
	if tf == 0 {
		return 0, ErrRange
	}
	return k.Data[0] / tf, nil
}

// LieAlgebraName returns a human-readable classical name for the given Dynkin
// type, such as "sl(4)" for A_3 or "so(7)" for B_3.
func LieAlgebraName(family string, rank int) (string, error) {
	f, err := normFamily(family)
	if err != nil {
		return "", err
	}
	if err := validRank(f, rank); err != nil {
		return "", err
	}
	switch f {
	case "A":
		return fmt.Sprintf("sl(%d)", rank+1), nil
	case "B":
		return fmt.Sprintf("so(%d)", 2*rank+1), nil
	case "C":
		return fmt.Sprintf("sp(%d)", 2*rank), nil
	case "D":
		return fmt.Sprintf("so(%d)", 2*rank), nil
	case "G":
		return "g2", nil
	case "F":
		return "f4", nil
	case "E":
		return fmt.Sprintf("e%d", rank), nil
	}
	return "", ErrType
}
