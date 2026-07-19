package operatortheory

import (
	"math"
	"math/cmplx"
	"sort"
)

// Spectrum returns the eigenvalues of a square matrix sorted lexicographically
// by real part and then imaginary part. It returns ErrNotSquare for a
// non-square matrix.
func (m *Matrix) Spectrum() ([]complex128, error) {
	vals, err := m.Eigenvalues()
	if err != nil {
		return nil, err
	}
	sort.Slice(vals, func(i, j int) bool {
		if real(vals[i]) != real(vals[j]) {
			return real(vals[i]) < real(vals[j])
		}
		return imag(vals[i]) < imag(vals[j])
	})
	return vals, nil
}

// SpectralRadius returns the largest modulus among the eigenvalues. It returns
// 0 for a non-square or empty matrix.
func (m *Matrix) SpectralRadius() float64 {
	if !m.IsSquare() || m.rows == 0 {
		return 0
	}
	vals := eigenvaluesQR(m)
	var r float64
	for _, v := range vals {
		if a := cmplx.Abs(v); a > r {
			r = a
		}
	}
	return r
}

// SpectralAbscissa returns the maximum real part among the eigenvalues, which
// governs the growth rate of the operator semigroup exp(tA).
func (m *Matrix) SpectralAbscissa() float64 {
	vals := eigenvaluesQR(m)
	a := math.Inf(-1)
	for _, v := range vals {
		if real(v) > a {
			a = real(v)
		}
	}
	if math.IsInf(a, -1) {
		return 0
	}
	return a
}

// Determinant is defined in decomp.go; Trace in matrix.go. SpectralGap below
// operates on the Hermitian spectrum.

// SpectralGap returns the difference between the two smallest distinct
// eigenvalues of a Hermitian matrix (the gap above the ground state). If there
// are fewer than two distinct eigenvalues it returns 0.
func (m *Matrix) SpectralGap() float64 {
	vals, _ := hermitianEigenRaw(m.HermitianPart())
	sort.Float64s(vals)
	if len(vals) < 2 {
		return 0
	}
	first := vals[0]
	for _, v := range vals[1:] {
		if v-first > 1e-9*(1+math.Abs(first)) {
			return v - first
		}
	}
	return 0
}

// MaxEigenvalue returns the largest eigenvalue of a Hermitian matrix.
func (m *Matrix) MaxEigenvalue() float64 {
	vals, _ := hermitianEigenRaw(m.HermitianPart())
	mx := math.Inf(-1)
	for _, v := range vals {
		if v > mx {
			mx = v
		}
	}
	return mx
}

// MinEigenvalue returns the smallest eigenvalue of a Hermitian matrix.
func (m *Matrix) MinEigenvalue() float64 {
	vals, _ := hermitianEigenRaw(m.HermitianPart())
	mn := math.Inf(1)
	for _, v := range vals {
		if v < mn {
			mn = v
		}
	}
	return mn
}

// Resolvent returns the resolvent operator (z*I - m)^{-1} at the complex point
// z. It returns ErrSingular when z is (numerically) an eigenvalue.
func (m *Matrix) Resolvent(z complex128) (*Matrix, error) {
	if !m.IsSquare() {
		return nil, ErrNotSquare
	}
	n := m.rows
	a := m.Scale(-1)
	for i := 0; i < n; i++ {
		a.data[i*n+i] += z
	}
	return a.Inverse()
}

// RayleighQuotient returns <v, m v> / <v, v> for a nonzero vector v. For a
// Hermitian operator this is real and lies between the smallest and largest
// eigenvalues.
func (m *Matrix) RayleighQuotient(v Vector) complex128 {
	av, err := m.MulVec(v)
	if err != nil {
		panic("operatortheory: dimension mismatch in RayleighQuotient")
	}
	d := v.Dot(v)
	if d == 0 {
		return 0
	}
	return v.Dot(av) / d
}

// QuadraticForm returns the sesquilinear form <x, m y>.
func (m *Matrix) QuadraticForm(x, y Vector) complex128 {
	my, err := m.MulVec(y)
	if err != nil {
		panic("operatortheory: dimension mismatch in QuadraticForm")
	}
	return x.Dot(my)
}

// NumericalAbscissa returns the maximum real part of the numerical range, equal
// to the largest eigenvalue of the Hermitian part (m + m^H)/2. It bounds the
// initial growth rate of exp(tA).
func (m *Matrix) NumericalAbscissa() float64 {
	return m.HermitianPart().MaxEigenvalue()
}

// NumericalRange returns a convex polygon approximating the boundary of the
// numerical range (field of values) W(m) = { <x, m x> : ||x|| = 1 }, sampled at
// the given number of angles (a minimum of 8 is used). Consecutive points trace
// the boundary counter-clockwise.
func (m *Matrix) NumericalRange(samples int) []complex128 {
	if samples < 8 {
		samples = 8
	}
	pts := make([]complex128, samples)
	for k := 0; k < samples; k++ {
		theta := 2 * math.Pi * float64(k) / float64(samples)
		phase := cmplx.Exp(complex(0, -theta))
		rot := m.Scale(phase)
		h := rot.HermitianPart()
		vals, vecs := hermitianEigenRaw(h)
		// Largest eigenvalue of the Hermitian part.
		best := 0
		for i := 1; i < len(vals); i++ {
			if vals[i] > vals[best] {
				best = i
			}
		}
		x := colOf(vecs, best)
		pts[k] = m.QuadraticForm(x, x)
	}
	return pts
}

// NumericalRadius returns the numerical radius w(m) = max { |z| : z in W(m) },
// estimated by sampling the field of values at the given number of angles.
func (m *Matrix) NumericalRadius(samples int) float64 {
	pts := m.NumericalRange(samples)
	var r float64
	for _, z := range pts {
		if a := cmplx.Abs(z); a > r {
			r = a
		}
	}
	return r
}

// CharacteristicPolynomial returns the coefficients of the characteristic
// polynomial det(x*I - m) in ascending order, so the result has length n+1 with
// a leading coefficient of 1 at index n. It uses the Faddeev-LeVerrier
// algorithm and returns ErrNotSquare for a non-square matrix.
func (m *Matrix) CharacteristicPolynomial() ([]complex128, error) {
	if !m.IsSquare() {
		return nil, ErrNotSquare
	}
	n := m.rows
	coeffs := make([]complex128, n+1)
	coeffs[n] = 1
	if n == 0 {
		return coeffs, nil
	}
	M := Identity(n)
	for k := 1; k <= n; k++ {
		am, _ := m.Mul(M)
		ck := -am.Trace() / complex(float64(k), 0)
		coeffs[n-k] = ck
		if k < n {
			M = am.Clone()
			for i := 0; i < n; i++ {
				M.data[i*n+i] += ck
			}
		}
	}
	return coeffs, nil
}

// EigenvalueMultiplicities groups the eigenvalues of a square matrix into
// clusters that agree to within tol and returns the distinct representatives
// together with their algebraic multiplicities. Representatives are sorted
// lexicographically.
func (m *Matrix) EigenvalueMultiplicities(tol float64) (values []complex128, mult []int, err error) {
	vals, err := m.Spectrum()
	if err != nil {
		return nil, nil, err
	}
	tol = orDefault(tol)
	for _, v := range vals {
		placed := false
		for i := range values {
			if cmplx.Abs(values[i]-v) <= tol {
				mult[i]++
				placed = true
				break
			}
		}
		if !placed {
			values = append(values, v)
			mult = append(mult, 1)
		}
	}
	return values, mult, nil
}

// Inertia returns the numbers of negative, zero and positive eigenvalues of a
// Hermitian matrix, using tol to decide which eigenvalues count as zero.
func (m *Matrix) Inertia(tol float64) (neg, zero, pos int) {
	tol = orDefault(tol)
	vals, _ := hermitianEigenRaw(m.HermitianPart())
	for _, v := range vals {
		switch {
		case v < -tol:
			neg++
		case v > tol:
			pos++
		default:
			zero++
		}
	}
	return neg, zero, pos
}
