package matrix

import (
	"errors"
	"math"
	"sort"
)

// errQRNoConverge reports that the Francis double-shift QR iteration used
// by [Matrix.EigenvaluesNumeric] failed to isolate an eigenvalue within the
// fixed iteration cap. It should not occur for well-formed real inputs and is
// returned rather than panicking so callers can react.
var errQRNoConverge = errors.New("matrix: QR iteration did not converge")

// Fixed, reproducible tuning constants for the numeric eigen routines. They are
// deliberately hard-coded (rather than derived from the input) so that repeated
// runs on the same matrix produce bit-identical results.
const (
	// matrixJacobiMaxSweeps caps the number of cyclic Jacobi sweeps in
	// [Matrix.EigSym]. Convergence is quadratic and normally completes in well
	// under ten sweeps; the cap only guards against pathological input.
	matrixJacobiMaxSweeps = 100
	// matrixSymRelTol is the scale-relative tolerance used to decide whether a
	// numeric matrix is symmetric. An entry pair (i,j),(j,i) must agree to
	// within matrixSymRelTol times the largest magnitude entry.
	matrixSymRelTol = 1e-9
	// matrixQRMaxIter is the per-eigenvalue iteration cap for the double-shift
	// QR sweep in [Matrix.EigenvaluesNumeric]. Exceptional shifts are injected
	// at fixed iteration counts below this cap.
	matrixQRMaxIter = 30
)

// EigSym computes the eigenvalues and eigenvectors of a real symmetric matrix
// using the cyclic Jacobi eigenvalue algorithm.
//
// The matrix must be square and, after conversion to float64 via
// [Matrix.Floats], symmetric to within a scale-relative tolerance (each pair of
// mirror entries agreeing to within matrixSymRelTol times the largest magnitude
// entry). It returns [ErrNotSquare] for a non-square matrix, and [ErrUnsupported]
// for a matrix with symbolic entries or one that is not numerically symmetric.
//
// The returned values are the eigenvalues in ascending order. The returned
// vectors matrix is square and orthogonal: its i-th column is the orthonormal
// eigenvector belonging to values[i], so that A·vectors[:,i] = values[i]·vectors[:,i].
//
// The Jacobi rotations are applied in place on a single flat symmetric buffer
// and the eigenvector matrix is accumulated in place, with per-sweep scratch
// reused across sweeps rather than reallocated.
func (m *Matrix) EigSym() (values []float64, vectors *Matrix, err error) {
	d, v, n, err := matrixEigSymCore(m)
	if err != nil {
		return nil, nil, err
	}
	rows := make([][]float64, n)
	for i := 0; i < n; i++ {
		rows[i] = make([]float64, n)
		for j := 0; j < n; j++ {
			rows[i][j] = v[i*n+j]
		}
	}
	return d, FromFloats(rows), nil
}

// EigSymValues computes only the eigenvalues of a real symmetric matrix, in
// ascending order, using the cyclic Jacobi eigenvalue algorithm. It shares all
// preconditions and error behavior with [Matrix.EigSym] but avoids materializing
// the eigenvector matrix.
func (m *Matrix) EigSymValues() ([]float64, error) {
	d, _, _, err := matrixEigSymCore(m)
	if err != nil {
		return nil, err
	}
	return d, nil
}

// matrixEigSymCore is the shared implementation behind [Matrix.EigSym] and
// [Matrix.EigSymValues]. It returns the ascending eigenvalues d, the eigenvector
// matrix v as a flat row-major n×n buffer whose i-th column corresponds to d[i],
// and the dimension n.
func matrixEigSymCore(m *Matrix) (d []float64, v []float64, n int, err error) {
	if !m.IsSquare() {
		return nil, nil, 0, ErrNotSquare
	}
	f, ferr := m.Floats()
	if ferr != nil {
		return nil, nil, 0, ErrUnsupported
	}
	n = m.rows
	if n == 0 {
		return nil, nil, 0, nil
	}

	// Largest magnitude entry sets the symmetry tolerance scale.
	var maxAbs float64
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if a := math.Abs(f[i][j]); a > maxAbs {
				maxAbs = a
			}
		}
	}
	tol := matrixSymRelTol * maxAbs

	// Verify symmetry and build the flat symmetric buffer, symmetrizing mirror
	// pairs by averaging so rounding noise below the tolerance is absorbed.
	a := make([]float64, n*n)
	for i := 0; i < n; i++ {
		for j := i; j < n; j++ {
			if math.Abs(f[i][j]-f[j][i]) > tol {
				return nil, nil, 0, ErrUnsupported
			}
			avg := 0.5 * (f[i][j] + f[j][i])
			a[i*n+j] = avg
			a[j*n+i] = avg
		}
	}

	d, vecs := matrixJacobiEigen(a, n)

	// Sort eigenvalues ascending and permute the eigenvector columns to match.
	order := make([]int, n)
	for i := range order {
		order[i] = i
	}
	sort.SliceStable(order, func(i, j int) bool { return d[order[i]] < d[order[j]] })

	sortedD := make([]float64, n)
	sortedV := make([]float64, n*n)
	for newCol, oldCol := range order {
		sortedD[newCol] = d[oldCol]
		for r := 0; r < n; r++ {
			sortedV[r*n+newCol] = vecs[r*n+oldCol]
		}
	}
	return sortedD, sortedV, n, nil
}

// matrixJacobiEigen diagonalizes the symmetric matrix held in the flat row-major
// buffer a (length n*n) using the cyclic Jacobi algorithm. It returns the
// eigenvalues d and an eigenvector matrix v (flat row-major, column j the
// eigenvector for d[j]). Only the upper triangle of a is consulted; a is
// modified in place. The scratch slices b and z are allocated once and reused
// across every sweep rather than inside the rotation loop.
func matrixJacobiEigen(a []float64, n int) (d []float64, v []float64) {
	v = make([]float64, n*n)
	for i := 0; i < n; i++ {
		v[i*n+i] = 1
	}
	d = make([]float64, n)
	b := make([]float64, n)
	z := make([]float64, n)
	for p := 0; p < n; p++ {
		d[p] = a[p*n+p]
		b[p] = d[p]
	}

	for sweep := 0; sweep < matrixJacobiMaxSweeps; sweep++ {
		// Sum of off-diagonal magnitudes; zero means we are diagonal.
		var sm float64
		for p := 0; p < n-1; p++ {
			for q := p + 1; q < n; q++ {
				sm += math.Abs(a[p*n+q])
			}
		}
		if sm == 0 {
			break
		}

		var thresh float64
		if sweep < 3 {
			thresh = 0.2 * sm / float64(n*n)
		}

		for p := 0; p < n-1; p++ {
			for q := p + 1; q < n; q++ {
				apq := a[p*n+q]
				g := 100.0 * math.Abs(apq)
				if sweep > 3 && math.Abs(d[p])+g == math.Abs(d[p]) &&
					math.Abs(d[q])+g == math.Abs(d[q]) {
					a[p*n+q] = 0
					continue
				}
				if math.Abs(apq) <= thresh {
					continue
				}
				h := d[q] - d[p]
				var t float64
				if math.Abs(h)+g == math.Abs(h) {
					t = apq / h
				} else {
					theta := 0.5 * h / apq
					t = 1.0 / (math.Abs(theta) + math.Sqrt(1.0+theta*theta))
					if theta < 0 {
						t = -t
					}
				}
				c := 1.0 / math.Sqrt(1+t*t)
				s := t * c
				tau := s / (1.0 + c)
				h = t * apq
				z[p] -= h
				z[q] += h
				d[p] -= h
				d[q] += h
				a[p*n+q] = 0
				for j := 0; j < p; j++ {
					matrixJacobiRotate(a, n, s, tau, j, p, j, q)
				}
				for j := p + 1; j < q; j++ {
					matrixJacobiRotate(a, n, s, tau, p, j, j, q)
				}
				for j := q + 1; j < n; j++ {
					matrixJacobiRotate(a, n, s, tau, p, j, q, j)
				}
				for j := 0; j < n; j++ {
					matrixJacobiRotate(v, n, s, tau, j, p, j, q)
				}
			}
		}

		for p := 0; p < n; p++ {
			b[p] += z[p]
			d[p] = b[p]
			z[p] = 0
		}
	}
	return d, v
}

// matrixJacobiRotate applies a single Jacobi Givens rotation to the flat
// row-major buffer m (stride n), mixing entries (i,j) and (k,l) using the
// precomputed sine s and tangent-of-half-angle tau. It operates in place.
func matrixJacobiRotate(m []float64, n int, s, tau float64, i, j, k, l int) {
	g := m[i*n+j]
	h := m[k*n+l]
	m[i*n+j] = g - s*(h+g*tau)
	m[k*n+l] = h + s*(g-h*tau)
}

// EigenvaluesNumeric computes all n eigenvalues of a general real square matrix.
//
// The matrix is first converted to float64 via [Matrix.Floats]; symbolic entries
// yield [ErrUnsupported] and a non-square matrix yields [ErrNotSquare]. The
// matrix is reduced to upper-Hessenberg form by Householder reflections and then
// the eigenvalues are extracted by the Francis double-shift QR iteration. Both
// phases run in place on a single row-major buffer, reusing a per-step
// Householder-vector scratch slice rather than allocating inside the iteration.
//
// All n eigenvalues are returned, including complex-conjugate pairs, in a
// deterministic order sorted by real part and then by imaginary part. Because a
// conjugate pair shares a real part, its two members are emitted adjacently. The
// iteration uses fixed exceptional-shift points and a fixed iteration cap for
// reproducibility; exceeding the cap returns a non-nil error.
func (m *Matrix) EigenvaluesNumeric() ([]complex128, error) {
	if !m.IsSquare() {
		return nil, ErrNotSquare
	}
	f, err := m.Floats()
	if err != nil {
		return nil, ErrUnsupported
	}
	n := m.rows
	if n == 0 {
		return nil, nil
	}

	// Work on a 1-based padded row-major buffer (stride n+1, row/col 0 unused)
	// so the classic reduction and QR recurrences index directly, in place.
	stride := n + 1
	a := make([]float64, stride*stride)
	for i := 1; i <= n; i++ {
		for j := 1; j <= n; j++ {
			a[i*stride+j] = f[i-1][j-1]
		}
	}

	ort := make([]float64, n+1) // reused Householder-vector scratch
	matrixHessenberg(a, n, stride, ort)

	wr := make([]float64, n+1)
	wi := make([]float64, n+1)
	if err := matrixHQR(a, n, stride, wr, wi); err != nil {
		return nil, err
	}

	out := make([]complex128, n)
	for i := 1; i <= n; i++ {
		out[i-1] = complex(wr[i], wi[i])
	}
	sort.SliceStable(out, func(i, j int) bool {
		if real(out[i]) != real(out[j]) {
			return real(out[i]) < real(out[j])
		}
		return imag(out[i]) < imag(out[j])
	})
	return out, nil
}

// matrixHessenberg reduces the real matrix stored in the 1-based padded buffer a
// (dimension n, row stride) to upper-Hessenberg form using Householder
// reflections, in place. The scratch slice ort (length >= n+1) holds the current
// Householder vector and is reused for every column rather than reallocated.
// Entries strictly below the first subdiagonal are zeroed on completion so the
// result is a clean Hessenberg matrix.
func matrixHessenberg(a []float64, n, stride int, ort []float64) {
	for mcol := 2; mcol <= n-1; mcol++ {
		var scale float64
		for i := n; i >= mcol; i-- {
			scale += math.Abs(a[i*stride+(mcol-1)])
		}
		if scale == 0 {
			continue
		}
		var h float64
		for i := n; i >= mcol; i-- {
			ort[i] = a[i*stride+(mcol-1)] / scale
			h += ort[i] * ort[i]
		}
		g := -math.Copysign(math.Sqrt(h), ort[mcol])
		h -= ort[mcol] * g
		ort[mcol] -= g
		// Apply (I - u·uᵀ/h) from the left to columns mcol..n.
		for j := mcol; j <= n; j++ {
			var fsum float64
			for i := n; i >= mcol; i-- {
				fsum += ort[i] * a[i*stride+j]
			}
			fsum /= h
			for i := mcol; i <= n; i++ {
				a[i*stride+j] -= fsum * ort[i]
			}
		}
		// Apply (I - u·uᵀ/h) from the right to rows 1..n.
		for i := 1; i <= n; i++ {
			var fsum float64
			for j := n; j >= mcol; j-- {
				fsum += ort[j] * a[i*stride+j]
			}
			fsum /= h
			for j := mcol; j <= n; j++ {
				a[i*stride+j] -= fsum * ort[j]
			}
		}
		a[mcol*stride+(mcol-1)] = scale * g
	}
	// Clean the strictly-below-subdiagonal region.
	for i := 3; i <= n; i++ {
		for j := 1; j <= i-2; j++ {
			a[i*stride+j] = 0
		}
	}
}

// matrixHQR runs the Francis double-shift QR iteration on the upper-Hessenberg
// matrix in the 1-based padded buffer a (dimension n, row stride), extracting the
// real and imaginary parts of every eigenvalue into wr and wi. The buffer is
// consumed in place. Exceptional shifts are injected at fixed iteration counts
// and the per-eigenvalue iteration count is capped at matrixQRMaxIter, returning
// errQRNoConverge if the cap is reached.
func matrixHQR(a []float64, n, stride int, wr, wi []float64) error {
	idx := func(i, j int) int { return i*stride + j }

	var anorm float64
	for i := 1; i <= n; i++ {
		lo := i - 1
		if lo < 1 {
			lo = 1
		}
		for j := lo; j <= n; j++ {
			anorm += math.Abs(a[idx(i, j)])
		}
	}

	nn := n
	var t float64
	for nn >= 1 {
		its := 0
		for {
			var l int
			for l = nn; l >= 2; l-- {
				s := math.Abs(a[idx(l-1, l-1)]) + math.Abs(a[idx(l, l)])
				if s == 0 {
					s = anorm
				}
				if math.Abs(a[idx(l, l-1)])+s == s {
					a[idx(l, l-1)] = 0
					break
				}
			}
			x := a[idx(nn, nn)]
			if l == nn {
				wr[nn] = x + t
				wi[nn] = 0
				nn--
			} else {
				y := a[idx(nn-1, nn-1)]
				w := a[idx(nn, nn-1)] * a[idx(nn-1, nn)]
				if l == nn-1 {
					p := 0.5 * (y - x)
					q := p*p + w
					z := math.Sqrt(math.Abs(q))
					x += t
					if q >= 0 {
						z = p + math.Copysign(z, p)
						wr[nn-1] = x + z
						wr[nn] = x + z
						if z != 0 {
							wr[nn] = x - w/z
						}
						wi[nn-1] = 0
						wi[nn] = 0
					} else {
						wr[nn-1] = x + p
						wr[nn] = x + p
						wi[nn] = z
						wi[nn-1] = -z
					}
					nn -= 2
				} else {
					if its == matrixQRMaxIter {
						return errQRNoConverge
					}
					if its == 10 || its == 20 {
						t += x
						for i := 1; i <= nn; i++ {
							a[idx(i, i)] -= x
						}
						s := math.Abs(a[idx(nn, nn-1)]) + math.Abs(a[idx(nn-1, nn-2)])
						x = 0.75 * s
						y = 0.75 * s
						w = -0.4375 * s * s
					}
					its++

					var m int
					var p, q, r float64
					for m = nn - 2; m >= l; m-- {
						z := a[idx(m, m)]
						rr := x - z
						ss := y - z
						p = (rr*ss-w)/a[idx(m+1, m)] + a[idx(m, m+1)]
						q = a[idx(m+1, m+1)] - z - rr - ss
						r = a[idx(m+2, m+1)]
						sc := math.Abs(p) + math.Abs(q) + math.Abs(r)
						p /= sc
						q /= sc
						r /= sc
						if m == l {
							break
						}
						u := math.Abs(a[idx(m, m-1)]) * (math.Abs(q) + math.Abs(r))
						vv := math.Abs(p) * (math.Abs(a[idx(m-1, m-1)]) +
							math.Abs(z) + math.Abs(a[idx(m+1, m+1)]))
						if u+vv == vv {
							break
						}
					}
					for i := m + 2; i <= nn; i++ {
						a[idx(i, i-2)] = 0
						if i != m+2 {
							a[idx(i, i-3)] = 0
						}
					}
					for k := m; k <= nn-1; k++ {
						if k != m {
							p = a[idx(k, k-1)]
							q = a[idx(k+1, k-1)]
							r = 0
							if k != nn-1 {
								r = a[idx(k+2, k-1)]
							}
							x = math.Abs(p) + math.Abs(q) + math.Abs(r)
							if x != 0 {
								p /= x
								q /= x
								r /= x
							}
						}
						s := math.Copysign(math.Sqrt(p*p+q*q+r*r), p)
						if s == 0 {
							continue
						}
						if k == m {
							if l != m {
								a[idx(k, k-1)] = -a[idx(k, k-1)]
							}
						} else {
							a[idx(k, k-1)] = -s * x
						}
						p += s
						x = p / s
						y = q / s
						z := r / s
						q /= p
						r /= p
						for j := k; j <= nn; j++ {
							pp := a[idx(k, j)] + q*a[idx(k+1, j)]
							if k != nn-1 {
								pp += r * a[idx(k+2, j)]
								a[idx(k+2, j)] -= pp * z
							}
							a[idx(k+1, j)] -= pp * y
							a[idx(k, j)] -= pp * x
						}
						mmin := nn
						if k+3 < nn {
							mmin = k + 3
						}
						for i := l; i <= mmin; i++ {
							pp := x*a[idx(i, k)] + y*a[idx(i, k+1)]
							if k != nn-1 {
								pp += z * a[idx(i, k+2)]
								a[idx(i, k+2)] -= pp * r
							}
							a[idx(i, k+1)] -= pp * q
							a[idx(i, k)] -= pp
						}
					}
				}
			}
			if l >= nn-1 {
				break
			}
		}
	}
	return nil
}
