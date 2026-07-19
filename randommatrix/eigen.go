package randommatrix

import (
	"math"
	"sort"
)

// EigenSym holds the eigen-decomposition of a real symmetric matrix. Values
// holds the eigenvalues in ascending order and column j of Vectors is the unit
// eigenvector belonging to Values[j].
type EigenSym struct {
	Values  []float64
	Vectors *Matrix
}

// jacobiOffNorm returns the square root of the sum of squared strictly
// upper-triangular entries.
func jacobiOffNorm(a *Matrix, n int) float64 {
	var s float64
	for p := 0; p < n; p++ {
		for q := p + 1; q < n; q++ {
			v := a.data[p*n+q]
			s += v * v
		}
	}
	return math.Sqrt(s)
}

// EigSymmetric computes the full eigen-decomposition of a real symmetric matrix
// using the cyclic Jacobi algorithm. It returns an error if the matrix is not
// square. The matrix is symmetrised internally, so the strictly lower triangle
// is ignored.
func EigSymmetric(m *Matrix) (*EigenSym, error) {
	if !m.IsSquare() {
		return nil, ErrNotSquare
	}
	n := m.rows
	a := m.Symmetrize()
	v := Identity(n)
	if n == 0 {
		return &EigenSym{Values: []float64{}, Vectors: v}, nil
	}
	frob := a.FrobeniusNorm()
	if frob == 0 {
		return &EigenSym{Values: make([]float64, n), Vectors: v}, nil
	}
	eps := 1e-15 * frob
	for sweep := 0; sweep < 100; sweep++ {
		off := jacobiOffNorm(a, n)
		if off <= eps {
			break
		}
		for p := 0; p < n-1; p++ {
			for q := p + 1; q < n; q++ {
				apq := a.data[p*n+q]
				if math.Abs(apq) <= 1e-300 {
					continue
				}
				app := a.data[p*n+p]
				aqq := a.data[q*n+q]
				theta := (aqq - app) / (2 * apq)
				var t float64
				if theta == 0 {
					t = 1
				} else {
					t = math.Copysign(1, theta) / (math.Abs(theta) + math.Sqrt(theta*theta+1))
				}
				c := 1 / math.Sqrt(t*t+1)
				s := t * c
				jacobiRotate(a, v, n, p, q, c, s)
			}
		}
	}
	vals := make([]float64, n)
	for i := 0; i < n; i++ {
		vals[i] = a.data[i*n+i]
	}
	// Sort ascending, permuting eigenvectors to match.
	idx := make([]int, n)
	for i := range idx {
		idx[i] = i
	}
	sort.Slice(idx, func(i, j int) bool { return vals[idx[i]] < vals[idx[j]] })
	sortedVals := make([]float64, n)
	sortedVecs := NewMatrix(n, n)
	for newCol, old := range idx {
		sortedVals[newCol] = vals[old]
		for r := 0; r < n; r++ {
			sortedVecs.data[r*n+newCol] = v.data[r*n+old]
		}
	}
	return &EigenSym{Values: sortedVals, Vectors: sortedVecs}, nil
}

// jacobiRotate applies a Jacobi rotation with cosine c and sine s in the (p,q)
// plane to the symmetric matrix a and accumulates it into the vector matrix v.
func jacobiRotate(a, v *Matrix, n, p, q int, c, s float64) {
	for k := 0; k < n; k++ {
		akp := a.data[k*n+p]
		akq := a.data[k*n+q]
		a.data[k*n+p] = c*akp - s*akq
		a.data[k*n+q] = s*akp + c*akq
	}
	for k := 0; k < n; k++ {
		apk := a.data[p*n+k]
		aqk := a.data[q*n+k]
		a.data[p*n+k] = c*apk - s*aqk
		a.data[q*n+k] = s*apk + c*aqk
	}
	// Force exact symmetry / zero of the eliminated off-diagonal pair.
	a.data[p*n+q] = 0
	a.data[q*n+p] = 0
	for k := 0; k < n; k++ {
		vkp := v.data[k*n+p]
		vkq := v.data[k*n+q]
		v.data[k*n+p] = c*vkp - s*vkq
		v.data[k*n+q] = s*vkp + c*vkq
	}
}

// EigenvaluesSymmetric returns the eigenvalues of a real symmetric matrix in
// ascending order.
func EigenvaluesSymmetric(m *Matrix) ([]float64, error) {
	e, err := EigSymmetric(m)
	if err != nil {
		return nil, err
	}
	return e.Values, nil
}

// EigenvaluesHermitian returns the eigenvalues of a Hermitian matrix in
// ascending order. It uses the real symmetric embedding of the Hermitian
// matrix H = A + iB into [[A, -B], [B, A]], whose 2n eigenvalues are those of H
// each repeated twice; the duplicates are collapsed.
func EigenvaluesHermitian(m *CMatrix) ([]float64, error) {
	if !m.IsSquare() {
		return nil, ErrNotSquare
	}
	n := m.rows
	if n == 0 {
		return []float64{}, nil
	}
	emb := NewMatrix(2*n, 2*n)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			z := m.data[i*n+j]
			a := real(z)
			b := imag(z)
			// top-left A
			emb.data[i*2*n+j] = a
			// bottom-right A
			emb.data[(i+n)*2*n+(j+n)] = a
			// top-right -B
			emb.data[i*2*n+(j+n)] = -b
			// bottom-left B
			emb.data[(i+n)*2*n+j] = b
		}
	}
	all, err := EigenvaluesSymmetric(emb)
	if err != nil {
		return nil, err
	}
	// The eigenvalues come in duplicate pairs; take every second one.
	out := make([]float64, 0, n)
	for i := 0; i < len(all); i += 2 {
		out = append(out, all[i])
	}
	return out, nil
}

// SpectralRadius returns the largest absolute eigenvalue of a real symmetric
// matrix.
func SpectralRadius(m *Matrix) (float64, error) {
	vals, err := EigenvaluesSymmetric(m)
	if err != nil {
		return 0, err
	}
	var r float64
	for _, v := range vals {
		if a := math.Abs(v); a > r {
			r = a
		}
	}
	return r, nil
}

// LargestEigenvalue returns the algebraically largest eigenvalue of a real
// symmetric matrix.
func LargestEigenvalue(m *Matrix) (float64, error) {
	vals, err := EigenvaluesSymmetric(m)
	if err != nil {
		return 0, err
	}
	if len(vals) == 0 {
		return 0, ErrDimensionMismatch
	}
	return vals[len(vals)-1], nil
}

// SmallestEigenvalue returns the algebraically smallest eigenvalue of a real
// symmetric matrix.
func SmallestEigenvalue(m *Matrix) (float64, error) {
	vals, err := EigenvaluesSymmetric(m)
	if err != nil {
		return 0, err
	}
	if len(vals) == 0 {
		return 0, ErrDimensionMismatch
	}
	return vals[0], nil
}

// DeterminantSymmetric returns the determinant of a real symmetric matrix as
// the product of its eigenvalues.
func DeterminantSymmetric(m *Matrix) (float64, error) {
	vals, err := EigenvaluesSymmetric(m)
	if err != nil {
		return 0, err
	}
	det := 1.0
	for _, v := range vals {
		det *= v
	}
	return det, nil
}

// ConditionNumberSymmetric returns the spectral condition number
// |lambda_max| / |lambda_min| of a real symmetric matrix. It returns +Inf when
// the smallest absolute eigenvalue is zero.
func ConditionNumberSymmetric(m *Matrix) (float64, error) {
	vals, err := EigenvaluesSymmetric(m)
	if err != nil {
		return 0, err
	}
	if len(vals) == 0 {
		return 0, ErrDimensionMismatch
	}
	mn, mx := math.Inf(1), 0.0
	for _, v := range vals {
		a := math.Abs(v)
		if a < mn {
			mn = a
		}
		if a > mx {
			mx = a
		}
	}
	if mn == 0 {
		return math.Inf(1), nil
	}
	return mx / mn, nil
}
