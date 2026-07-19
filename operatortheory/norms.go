package operatortheory

import (
	"math"
	"math/cmplx"
)

// FrobeniusNorm returns the Frobenius norm, the square root of the sum of the
// squared moduli of the entries.
func (m *Matrix) FrobeniusNorm() float64 {
	var s float64
	for _, z := range m.data {
		s += real(z)*real(z) + imag(z)*imag(z)
	}
	return math.Sqrt(s)
}

// OneNorm returns the induced 1-norm, the maximum absolute column sum.
func (m *Matrix) OneNorm() float64 {
	var mx float64
	for j := 0; j < m.cols; j++ {
		var s float64
		for i := 0; i < m.rows; i++ {
			s += cmplx.Abs(m.data[i*m.cols+j])
		}
		if s > mx {
			mx = s
		}
	}
	return mx
}

// InfNorm returns the induced infinity-norm, the maximum absolute row sum.
func (m *Matrix) InfNorm() float64 {
	var mx float64
	for i := 0; i < m.rows; i++ {
		var s float64
		for j := 0; j < m.cols; j++ {
			s += cmplx.Abs(m.data[i*m.cols+j])
		}
		if s > mx {
			mx = s
		}
	}
	return mx
}

// MaxNorm returns the largest modulus among the entries (the entrywise
// max-norm).
func (m *Matrix) MaxNorm() float64 { return m.MaxAbs() }

// OperatorNorm returns the induced 2-norm (spectral norm), the largest singular
// value. It is computed by power iteration on m^H m for efficiency.
func (m *Matrix) OperatorNorm() float64 {
	n := m.cols
	if n == 0 || m.rows == 0 {
		return 0
	}
	ah := m.Adjoint()
	// Power iteration on A^H A.
	x := make(Vector, n)
	for i := range x {
		x[i] = complex(1+0.1*float64(i), 0)
	}
	x, _ = x.Normalize()
	var lambda float64
	for iter := 0; iter < 1000; iter++ {
		av, _ := m.MulVec(x)
		y, _ := ah.MulVec(av)
		ny := y.Norm()
		if ny == 0 {
			return 0
		}
		newX := y.Scale(complex(1/ny, 0))
		if math.Abs(ny-lambda) <= 1e-14*(1+ny) {
			lambda = ny
			x = newX
			break
		}
		lambda = ny
		x = newX
	}
	return math.Sqrt(lambda)
}

// SpectralNorm is an alias for OperatorNorm.
func (m *Matrix) SpectralNorm() float64 { return m.OperatorNorm() }

// NuclearNorm returns the nuclear (trace) norm, the sum of the singular values.
func (m *Matrix) NuclearNorm() float64 {
	var s float64
	for _, sv := range m.SingularValues() {
		s += sv
	}
	return s
}

// TraceNorm is an alias for NuclearNorm.
func (m *Matrix) TraceNorm() float64 { return m.NuclearNorm() }

// SchattenNorm returns the Schatten p-norm, the l^p norm of the vector of
// singular values, for p >= 1. SchattenNorm(1) is the nuclear norm,
// SchattenNorm(2) the Frobenius norm and the limit p -> infinity the operator
// norm.
func (m *Matrix) SchattenNorm(p float64) float64 {
	if p < 1 {
		p = 1
	}
	var s float64
	for _, sv := range m.SingularValues() {
		s += math.Pow(sv, p)
	}
	return math.Pow(s, 1/p)
}

// KyFanNorm returns the Ky Fan k-norm, the sum of the k largest singular
// values. If k exceeds the number of singular values it uses all of them.
func (m *Matrix) KyFanNorm(k int) float64 {
	s := m.SingularValues()
	if k > len(s) {
		k = len(s)
	}
	var sum float64
	for i := 0; i < k; i++ {
		sum += s[i]
	}
	return sum
}

// DistanceFrobenius returns the Frobenius norm of m - b, a metric on the space
// of matrices of equal shape. It returns ErrDimensionMismatch on a shape
// mismatch.
func (m *Matrix) DistanceFrobenius(b *Matrix) (float64, error) {
	d, err := m.Sub(b)
	if err != nil {
		return 0, err
	}
	return d.FrobeniusNorm(), nil
}
