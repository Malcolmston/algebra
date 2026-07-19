package randommatrix

import (
	"math"
	"math/rand"
)

// Beta values identify the three classical symmetry classes by their Dyson
// index: 1 for the orthogonal (real) class, 2 for the unitary (complex) class
// and 4 for the symplectic (quaternion) class.
const (
	BetaOrthogonal = 1
	BetaUnitary    = 2
	BetaSymplectic = 4
)

// newRNG returns a math/rand generator seeded deterministically by seed.
func newRNG(seed int64) *rand.Rand {
	return rand.New(rand.NewSource(seed))
}

// GaussianVector returns n independent N(mean, std^2) samples.
func GaussianVector(n int, mean, std float64, seed int64) []float64 {
	r := newRNG(seed)
	out := make([]float64, n)
	for i := range out {
		out[i] = mean + std*r.NormFloat64()
	}
	return out
}

// GaussianMatrix returns a rows-by-cols matrix of independent N(mean, std^2)
// entries.
func GaussianMatrix(rows, cols int, mean, std float64, seed int64) *Matrix {
	r := newRNG(seed)
	m := NewMatrix(rows, cols)
	for i := range m.data {
		m.data[i] = mean + std*r.NormFloat64()
	}
	return m
}

// GaussianCMatrix returns a rows-by-cols complex matrix whose real and
// imaginary parts are independent N(0, std^2) samples.
func GaussianCMatrix(rows, cols int, std float64, seed int64) *CMatrix {
	r := newRNG(seed)
	m := NewCMatrix(rows, cols)
	for i := range m.data {
		m.data[i] = complex(std*r.NormFloat64(), std*r.NormFloat64())
	}
	return m
}

// GOE returns a sample from the Gaussian Orthogonal Ensemble: a real symmetric
// n-by-n matrix whose off-diagonal entries are independent N(0,1) and whose
// diagonal entries are independent N(0,2). With this normalisation the
// eigenvalues of GOE(n)/sqrt(n) converge to the semicircle law on [-2,2].
func GOE(n int, seed int64) *Matrix {
	r := newRNG(seed)
	m := NewMatrix(n, n)
	for i := 0; i < n; i++ {
		m.data[i*n+i] = math.Sqrt2 * r.NormFloat64()
		for j := i + 1; j < n; j++ {
			v := r.NormFloat64()
			m.data[i*n+j] = v
			m.data[j*n+i] = v
		}
	}
	return m
}

// GUE returns a sample from the Gaussian Unitary Ensemble: a Hermitian n-by-n
// matrix whose diagonal entries are real N(0,1) and whose off-diagonal entries
// have real and imaginary parts drawn independently from N(0,1/2). The
// eigenvalues of GUE(n)/sqrt(n) converge to the semicircle law on [-2,2].
func GUE(n int, seed int64) *CMatrix {
	r := newRNG(seed)
	m := NewCMatrix(n, n)
	s := 1 / math.Sqrt2
	for i := 0; i < n; i++ {
		m.data[i*n+i] = complex(r.NormFloat64(), 0)
		for j := i + 1; j < n; j++ {
			z := complex(s*r.NormFloat64(), s*r.NormFloat64())
			m.data[i*n+j] = z
			m.data[j*n+i] = complexConj(z)
		}
	}
	return m
}

// complexConj returns the complex conjugate of z.
func complexConj(z complex128) complex128 { return complex(real(z), -imag(z)) }

// GSE returns a sample from the Gaussian Symplectic Ensemble as a 2n-by-2n
// Hermitian, self-dual complex matrix built from an n-by-n array of quaternions.
// Each quaternion is embedded as a 2-by-2 complex block. Diagonal quaternions
// are real N(0,1); the four components of each off-diagonal quaternion are
// independent N(0,1/2). Every eigenvalue is doubly degenerate (Kramers
// degeneracy).
func GSE(n int, seed int64) *CMatrix {
	r := newRNG(seed)
	m := NewCMatrix(2*n, 2*n)
	s := 1 / math.Sqrt2
	set := func(i, j int, q0, q1, q2, q3 float64) {
		// Quaternion q0 + q1 i + q2 j + q3 k as a 2x2 complex block.
		a := complex(q0, q1)
		b := complex(q2, q3)
		bi, bj := 2*i, 2*j
		m.data[bi*2*n+bj] = a
		m.data[bi*2*n+bj+1] = b
		m.data[(bi+1)*2*n+bj] = -complexConj(b)
		m.data[(bi+1)*2*n+bj+1] = complexConj(a)
	}
	for i := 0; i < n; i++ {
		set(i, i, r.NormFloat64(), 0, 0, 0)
		for j := i + 1; j < n; j++ {
			q0 := s * r.NormFloat64()
			q1 := s * r.NormFloat64()
			q2 := s * r.NormFloat64()
			q3 := s * r.NormFloat64()
			set(i, j, q0, q1, q2, q3)
			// Hermitian self-dual conjugate block.
			set(j, i, q0, -q1, -q2, -q3)
		}
	}
	return m
}

// WignerMatrix returns a real symmetric n-by-n Wigner matrix with independent
// off-diagonal entries of standard deviation offStd and diagonal entries of
// standard deviation diagStd.
func WignerMatrix(n int, offStd, diagStd float64, seed int64) *Matrix {
	r := newRNG(seed)
	m := NewMatrix(n, n)
	for i := 0; i < n; i++ {
		m.data[i*n+i] = diagStd * r.NormFloat64()
		for j := i + 1; j < n; j++ {
			v := offStd * r.NormFloat64()
			m.data[i*n+j] = v
			m.data[j*n+i] = v
		}
	}
	return m
}

// GinibreReal returns an n-by-n matrix of independent standard normal real
// entries (the real Ginibre ensemble).
func GinibreReal(n int, seed int64) *Matrix {
	return GaussianMatrix(n, n, 0, 1, seed)
}

// GinibreComplex returns an n-by-n complex matrix whose entries have real and
// imaginary parts drawn independently from N(0,1/2), so that each entry has unit
// expected squared modulus (the complex Ginibre ensemble).
func GinibreComplex(n int, seed int64) *CMatrix {
	return GaussianCMatrix(n, n, 1/math.Sqrt2, seed)
}

// WishartReal returns the real Wishart matrix G Gᵀ, where G is an n-by-p matrix
// of independent standard normal entries. The result is an n-by-n symmetric
// positive semidefinite matrix (the Laguerre orthogonal ensemble).
func WishartReal(n, p int, seed int64) *Matrix {
	g := GaussianMatrix(n, p, 0, 1, seed)
	gt := g.Transpose()
	w, _ := g.Mul(gt)
	return w
}

// LaguerreOrthogonal is a synonym for WishartReal, emphasising its role as the
// beta = 1 Laguerre ensemble.
func LaguerreOrthogonal(n, p int, seed int64) *Matrix { return WishartReal(n, p, seed) }

// WishartComplex returns the complex Wishart matrix G G†, where G is an n-by-p
// complex matrix with real and imaginary parts drawn independently from
// N(0,1/2). The result is an n-by-n Hermitian positive semidefinite matrix (the
// Laguerre unitary ensemble).
func WishartComplex(n, p int, seed int64) *CMatrix {
	g := GaussianCMatrix(n, p, 1/math.Sqrt2, seed)
	gd := g.ConjugateTranspose()
	w, _ := g.Mul(gd)
	return w
}

// LaguerreUnitary is a synonym for WishartComplex.
func LaguerreUnitary(n, p int, seed int64) *CMatrix { return WishartComplex(n, p, seed) }

// SampleCovarianceReal returns (1/p) G Gᵀ where G is n-by-p standard normal.
// With ratio c = n/p its eigenvalues follow the Marchenko-Pastur law with
// variance one and ratio c, supported on [(1-sqrt(c))^2, (1+sqrt(c))^2] when
// c <= 1.
func SampleCovarianceReal(n, p int, seed int64) *Matrix {
	return WishartReal(n, p, seed).Scale(1 / float64(p))
}

// SampleGOE returns the eigenvalues of GOE(n)/scale in ascending order. Passing
// scale = sqrt(n) places the bulk on [-2,2].
func SampleGOE(n int, seed int64, scale float64) ([]float64, error) {
	m := GOE(n, seed)
	if scale != 1 && scale != 0 {
		m = m.Scale(1 / scale)
	}
	return EigenvaluesSymmetric(m)
}

// SampleGUE returns the eigenvalues of GUE(n)/scale in ascending order.
func SampleGUE(n int, seed int64, scale float64) ([]float64, error) {
	m := GUE(n, seed)
	vals, err := EigenvaluesHermitian(m)
	if err != nil {
		return nil, err
	}
	if scale != 1 && scale != 0 {
		for i := range vals {
			vals[i] /= scale
		}
	}
	return vals, nil
}

// SampleWishartReal returns the eigenvalues of the real sample covariance
// matrix (1/p) G Gᵀ in ascending order.
func SampleWishartReal(n, p int, seed int64) ([]float64, error) {
	return EigenvaluesSymmetric(SampleCovarianceReal(n, p, seed))
}
