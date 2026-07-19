package operatortheory

import (
	"math"
	"math/cmplx"
)

// luResult holds an in-place LU factorisation with partial pivoting.
type luResult struct {
	n        int
	lu       []complex128 // combined L (unit lower) and U (upper)
	piv      []int        // row permutation
	sign     float64      // sign of the permutation
	singular bool
}

// luDecompose computes the LU factorisation of the square matrix m with partial
// pivoting. The input is not modified.
func luDecompose(m *Matrix) *luResult {
	n := m.rows
	lu := make([]complex128, n*n)
	copy(lu, m.data)
	piv := make([]int, n)
	for i := range piv {
		piv[i] = i
	}
	sign := 1.0
	singular := false
	for k := 0; k < n; k++ {
		// Find pivot.
		p := k
		max := cmplx.Abs(lu[k*n+k])
		for i := k + 1; i < n; i++ {
			if a := cmplx.Abs(lu[i*n+k]); a > max {
				max = a
				p = i
			}
		}
		if max == 0 {
			singular = true
			continue
		}
		if p != k {
			for j := 0; j < n; j++ {
				lu[k*n+j], lu[p*n+j] = lu[p*n+j], lu[k*n+j]
			}
			piv[k], piv[p] = piv[p], piv[k]
			sign = -sign
		}
		akk := lu[k*n+k]
		for i := k + 1; i < n; i++ {
			f := lu[i*n+k] / akk
			lu[i*n+k] = f
			for j := k + 1; j < n; j++ {
				lu[i*n+j] -= f * lu[k*n+j]
			}
		}
	}
	return &luResult{n: n, lu: lu, piv: piv, sign: sign, singular: singular}
}

// solve solves A x = b for a single right-hand side using the factorisation.
func (f *luResult) solve(b Vector) Vector {
	n := f.n
	x := make(Vector, n)
	for i := 0; i < n; i++ {
		x[i] = b[f.piv[i]]
	}
	// Forward substitution (unit lower).
	for i := 0; i < n; i++ {
		for j := 0; j < i; j++ {
			x[i] -= f.lu[i*n+j] * x[j]
		}
	}
	// Back substitution (upper).
	for i := n - 1; i >= 0; i-- {
		for j := i + 1; j < n; j++ {
			x[i] -= f.lu[i*n+j] * x[j]
		}
		x[i] /= f.lu[i*n+i]
	}
	return x
}

// det returns the determinant from the factorisation.
func (f *luResult) det() complex128 {
	if f.singular {
		return 0
	}
	d := complex(f.sign, 0)
	for i := 0; i < f.n; i++ {
		d *= f.lu[i*f.n+i]
	}
	return d
}

// qrRaw returns a (thin, reduced) QR factorisation of m using complex
// Householder reflections. Q is rows-by-min(rows,cols) with orthonormal
// columns and R is min-by-cols upper triangular. For square m these are the
// full factors.
func (m *Matrix) qrRaw() (*Matrix, *Matrix) {
	rows, cols := m.rows, m.cols
	k := cols
	if rows < k {
		k = rows
	}
	// Work on a copy in column form.
	r := m.Clone()
	// Accumulate Q as rows-by-rows, then trim.
	q := Identity(rows)
	for j := 0; j < k; j++ {
		// Householder vector for column j below the diagonal.
		var norm float64
		for i := j; i < rows; i++ {
			a := cmplx.Abs(r.data[i*cols+j])
			norm += a * a
		}
		norm = math.Sqrt(norm)
		if norm == 0 {
			continue
		}
		alpha := r.data[j*cols+j]
		var phase complex128 = 1
		if a := cmplx.Abs(alpha); a > 0 {
			phase = alpha / complex(a, 0)
		}
		// v = x + phase*norm*e_j (reflect to -phase*norm to avoid cancellation).
		v := make([]complex128, rows)
		for i := j; i < rows; i++ {
			v[i] = r.data[i*cols+j]
		}
		v[j] += phase * complex(norm, 0)
		var vnorm float64
		for i := j; i < rows; i++ {
			a := cmplx.Abs(v[i])
			vnorm += a * a
		}
		if vnorm == 0 {
			continue
		}
		// Apply H = I - 2 v v^H / (v^H v) to R (from the left).
		for c := 0; c < cols; c++ {
			var dot complex128
			for i := j; i < rows; i++ {
				dot += cmplx.Conj(v[i]) * r.data[i*cols+c]
			}
			f := 2 * dot / complex(vnorm, 0)
			for i := j; i < rows; i++ {
				r.data[i*cols+c] -= f * v[i]
			}
		}
		// Apply H to Q (from the right: Q = Q H, H Hermitian).
		for rr := 0; rr < rows; rr++ {
			var dot complex128
			for i := j; i < rows; i++ {
				dot += q.data[rr*rows+i] * v[i]
			}
			f := 2 * dot / complex(vnorm, 0)
			for i := j; i < rows; i++ {
				q.data[rr*rows+i] -= f * cmplx.Conj(v[i])
			}
		}
	}
	// Trim Q to rows-by-k and R to k-by-cols.
	qt := NewMatrix(rows, k)
	for i := 0; i < rows; i++ {
		for j := 0; j < k; j++ {
			qt.data[i*k+j] = q.data[i*rows+j]
		}
	}
	rt := NewMatrix(k, cols)
	for i := 0; i < k; i++ {
		for j := 0; j < cols; j++ {
			if j >= i {
				rt.data[i*cols+j] = r.data[i*cols+j]
			}
		}
	}
	return qt, rt
}

// hessenbergRaw reduces a square matrix to upper-Hessenberg form H = Q^H A Q by
// Householder reflections and returns H and Q.
func (m *Matrix) hessenbergRaw() (*Matrix, *Matrix) {
	n := m.rows
	h := m.Clone()
	q := Identity(n)
	for k := 0; k < n-2; k++ {
		var norm float64
		for i := k + 1; i < n; i++ {
			a := cmplx.Abs(h.data[i*n+k])
			norm += a * a
		}
		norm = math.Sqrt(norm)
		if norm == 0 {
			continue
		}
		alpha := h.data[(k+1)*n+k]
		var phase complex128 = 1
		if a := cmplx.Abs(alpha); a > 0 {
			phase = alpha / complex(a, 0)
		}
		v := make([]complex128, n)
		for i := k + 1; i < n; i++ {
			v[i] = h.data[i*n+k]
		}
		v[k+1] += phase * complex(norm, 0)
		var vnorm float64
		for i := k + 1; i < n; i++ {
			a := cmplx.Abs(v[i])
			vnorm += a * a
		}
		if vnorm == 0 {
			continue
		}
		// H <- P H (left).
		for c := 0; c < n; c++ {
			var dot complex128
			for i := k + 1; i < n; i++ {
				dot += cmplx.Conj(v[i]) * h.data[i*n+c]
			}
			f := 2 * dot / complex(vnorm, 0)
			for i := k + 1; i < n; i++ {
				h.data[i*n+c] -= f * v[i]
			}
		}
		// H <- H P (right).
		for rr := 0; rr < n; rr++ {
			var dot complex128
			for i := k + 1; i < n; i++ {
				dot += h.data[rr*n+i] * v[i]
			}
			f := 2 * dot / complex(vnorm, 0)
			for i := k + 1; i < n; i++ {
				h.data[rr*n+i] -= f * cmplx.Conj(v[i])
			}
		}
		// Q <- Q P (right).
		for rr := 0; rr < n; rr++ {
			var dot complex128
			for i := k + 1; i < n; i++ {
				dot += q.data[rr*n+i] * v[i]
			}
			f := 2 * dot / complex(vnorm, 0)
			for i := k + 1; i < n; i++ {
				q.data[rr*n+i] -= f * cmplx.Conj(v[i])
			}
		}
	}
	return h, q
}

// eigenvaluesQR computes all eigenvalues of a square matrix via the explicitly
// shifted QR algorithm on the Hessenberg form. It returns the eigenvalues in no
// particular order.
func eigenvaluesQR(m *Matrix) []complex128 {
	n := m.rows
	if n == 0 {
		return nil
	}
	if n == 1 {
		return []complex128{m.data[0]}
	}
	h, _ := m.hessenbergRaw()
	H := h.data
	eig := make([]complex128, n)
	hi := n
	maxIter := 100 * n
	iter := 0
	for hi > 0 {
		if hi == 1 {
			eig[0] = H[0]
			hi = 0
			break
		}
		// Find split point l (largest l with H[l][l-1] negligible).
		l := hi - 1
		for l > 0 {
			s := cmplx.Abs(H[(l-1)*n+(l-1)]) + cmplx.Abs(H[l*n+l])
			if s == 0 {
				s = 1
			}
			if cmplx.Abs(H[l*n+(l-1)]) <= 1e-16*s {
				H[l*n+(l-1)] = 0
				break
			}
			l--
		}
		if l == hi-1 {
			eig[hi-1] = H[(hi-1)*n+(hi-1)]
			hi--
			iter = 0
			continue
		}
		if l == hi-2 {
			a := H[(hi-2)*n+(hi-2)]
			b := H[(hi-2)*n+(hi-1)]
			c := H[(hi-1)*n+(hi-2)]
			d := H[(hi-1)*n+(hi-1)]
			e1, e2 := eig2x2(a, b, c, d)
			eig[hi-2] = e1
			eig[hi-1] = e2
			hi -= 2
			iter = 0
			continue
		}
		iter++
		if iter > maxIter {
			// Fallback: accept the current diagonal of the active block.
			for i := l; i < hi; i++ {
				eig[i] = H[i*n+i]
			}
			hi = l
			iter = 0
			continue
		}
		mu := wilkinsonShift(H, n, hi)
		qrStepHessenberg(H, n, l, hi, mu)
	}
	return eig
}

// eig2x2 returns the two eigenvalues of the 2x2 matrix [[a,b],[c,d]].
func eig2x2(a, b, c, d complex128) (complex128, complex128) {
	tr := a + d
	det := a*d - b*c
	disc := cmplx.Sqrt(tr*tr - 4*det)
	return (tr + disc) / 2, (tr - disc) / 2
}

// wilkinsonShift returns the eigenvalue of the trailing 2x2 block of the active
// window [.., hi) that is closer to H[hi-1][hi-1].
func wilkinsonShift(H []complex128, n, hi int) complex128 {
	a := H[(hi-2)*n+(hi-2)]
	b := H[(hi-2)*n+(hi-1)]
	c := H[(hi-1)*n+(hi-2)]
	d := H[(hi-1)*n+(hi-1)]
	e1, e2 := eig2x2(a, b, c, d)
	if cmplx.Abs(e1-d) <= cmplx.Abs(e2-d) {
		return e1
	}
	return e2
}

// qrStepHessenberg performs one explicitly shifted QR step on the active
// Hessenberg window rows/cols [l, hi) with shift mu.
func qrStepHessenberg(H []complex128, n, l, hi int, mu complex128) {
	for i := l; i < hi; i++ {
		H[i*n+i] -= mu
	}
	cs := make([]float64, hi)
	sn := make([]complex128, hi)
	// QR: zero subdiagonals with Givens rotations.
	for i := l; i < hi-1; i++ {
		a := H[i*n+i]
		b := H[(i+1)*n+i]
		c, s := givens(a, b)
		cs[i] = c
		sn[i] = s
		// Apply to rows i, i+1 across columns [l, hi).
		for j := l; j < hi; j++ {
			h1 := H[i*n+j]
			h2 := H[(i+1)*n+j]
			H[i*n+j] = complex(c, 0)*h1 + s*h2
			H[(i+1)*n+j] = -cmplx.Conj(s)*h1 + complex(c, 0)*h2
		}
	}
	// RQ: apply the transposed rotations from the right.
	for i := l; i < hi-1; i++ {
		c := cs[i]
		s := sn[i]
		for rr := l; rr < hi; rr++ {
			h1 := H[rr*n+i]
			h2 := H[rr*n+i+1]
			H[rr*n+i] = complex(c, 0)*h1 + cmplx.Conj(s)*h2
			H[rr*n+i+1] = -s*h1 + complex(c, 0)*h2
		}
	}
	for i := l; i < hi; i++ {
		H[i*n+i] += mu
	}
}

// givens returns the parameters (c real, s complex) of the 2x2 unitary
// [[c, s], [-conj(s), c]] that maps (a, b)^T to (r, 0)^T.
func givens(a, b complex128) (float64, complex128) {
	if b == 0 {
		return 1, 0
	}
	if a == 0 {
		return 0, cmplx.Conj(b) / complex(cmplx.Abs(b), 0)
	}
	absa := cmplx.Abs(a)
	absb := cmplx.Abs(b)
	denom := math.Hypot(absa, absb)
	c := absa / denom
	s := complex(c, 0) * cmplx.Conj(b) / cmplx.Conj(a)
	return c, s
}

// realSymmetricJacobi diagonalises the n-by-n real symmetric matrix a (given
// row-major) with the cyclic Jacobi method. It returns the eigenvalues and the
// eigenvectors as columns of an orthogonal matrix (row-major, n-by-n).
func realSymmetricJacobi(a []float64, n int) (vals []float64, vecs []float64) {
	A := make([]float64, n*n)
	copy(A, a)
	V := make([]float64, n*n)
	for i := 0; i < n; i++ {
		V[i*n+i] = 1
	}
	for sweep := 0; sweep < 100; sweep++ {
		var off float64
		for p := 0; p < n; p++ {
			for q := p + 1; q < n; q++ {
				off += A[p*n+q] * A[p*n+q]
			}
		}
		if off <= 1e-30 {
			break
		}
		for p := 0; p < n; p++ {
			for q := p + 1; q < n; q++ {
				apq := A[p*n+q]
				if math.Abs(apq) <= 1e-300 {
					continue
				}
				app := A[p*n+p]
				aqq := A[q*n+q]
				phi := 0.5 * math.Atan2(2*apq, aqq-app)
				c := math.Cos(phi)
				s := math.Sin(phi)
				// Rotate rows/cols p and q.
				for k := 0; k < n; k++ {
					akp := A[k*n+p]
					akq := A[k*n+q]
					A[k*n+p] = c*akp - s*akq
					A[k*n+q] = s*akp + c*akq
				}
				for k := 0; k < n; k++ {
					apk := A[p*n+k]
					aqk := A[q*n+k]
					A[p*n+k] = c*apk - s*aqk
					A[q*n+k] = s*apk + c*aqk
				}
				for k := 0; k < n; k++ {
					vkp := V[k*n+p]
					vkq := V[k*n+q]
					V[k*n+p] = c*vkp - s*vkq
					V[k*n+q] = s*vkp + c*vkq
				}
			}
		}
	}
	vals = make([]float64, n)
	for i := 0; i < n; i++ {
		vals[i] = A[i*n+i]
	}
	return vals, V
}

// hermitianEigenRaw computes the eigenvalues (ascending) and orthonormal
// eigenvectors of the Hermitian matrix m via the real symmetric 2n-by-2n
// embedding. Eigenvectors are returned as the columns of the returned matrix.
func hermitianEigenRaw(m *Matrix) ([]float64, *Matrix) {
	n := m.rows
	// Build the real symmetric embedding [[Re, -Im],[Im, Re]].
	N := 2 * n
	B := make([]float64, N*N)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			re := real(m.data[i*n+j])
			im := imag(m.data[i*n+j])
			B[i*N+j] = re
			B[i*N+(j+n)] = -im
			B[(i+n)*N+j] = im
			B[(i+n)*N+(j+n)] = re
		}
	}
	rawVals, rawVecs := realSymmetricJacobi(B, N)
	// Order indices by eigenvalue ascending.
	idx := make([]int, N)
	for i := range idx {
		idx[i] = i
	}
	sortIndicesByValue(idx, rawVals)
	// Build complex candidate eigenvectors x = u + i v (u top half, v bottom).
	type cand struct {
		val float64
		vec Vector
	}
	cands := make([]cand, N)
	for c := 0; c < N; c++ {
		col := idx[c]
		x := make(Vector, n)
		for i := 0; i < n; i++ {
			u := rawVecs[i*N+col]
			w := rawVecs[(i+n)*N+col]
			x[i] = complex(u, w)
		}
		cands[c] = cand{val: rawVals[col], vec: x}
	}
	// Pick n orthonormal eigenvectors via complex Gram-Schmidt in ascending
	// eigenvalue order.
	vals := make([]float64, 0, n)
	vecs := NewMatrix(n, n)
	chosen := make([]Vector, 0, n)
	col := 0
	for _, cd := range cands {
		if len(chosen) == n {
			break
		}
		u := cd.vec.Clone()
		for _, q := range chosen {
			u = u.Sub(q.Scale(q.Dot(u)))
		}
		nrm := u.Norm()
		if nrm <= 1e-9 {
			continue
		}
		u = u.Scale(complex(1/nrm, 0))
		chosen = append(chosen, u)
		vals = append(vals, cd.val)
		for i := 0; i < n; i++ {
			vecs.data[i*n+col] = u[i]
		}
		col++
	}
	return vals, vecs
}

// sortIndicesByValue sorts idx so that vals[idx[0]] <= vals[idx[1]] <= ...
func sortIndicesByValue(idx []int, vals []float64) {
	// Simple insertion sort; n is small in practice.
	for i := 1; i < len(idx); i++ {
		j := i
		for j > 0 && vals[idx[j-1]] > vals[idx[j]] {
			idx[j-1], idx[j] = idx[j], idx[j-1]
			j--
		}
	}
}
