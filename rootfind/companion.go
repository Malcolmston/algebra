package rootfind

import "math"

// CompanionMatrix returns the companion matrix of the polynomial p as a dense
// row-major n-by-n slice, where n is the degree of p. The characteristic
// polynomial of the returned matrix equals p made monic, so its eigenvalues are
// exactly the roots of p. The matrix is in upper-Hessenberg form with a unit
// subdiagonal. It returns ErrDegreeTooLow for constant polynomials.
func CompanionMatrix(p Poly) ([][]float64, error) {
	n := p.Degree()
	if n < 1 {
		return nil, ErrDegreeTooLow
	}
	lc := p[n]
	a := make([][]float64, n)
	for i := range a {
		a[i] = make([]float64, n)
	}
	// First row: -a_{n-1}/a_n, -a_{n-2}/a_n, ... , -a_0/a_n.
	for j := 0; j < n; j++ {
		a[0][j] = -p[n-1-j] / lc
	}
	// Subdiagonal ones.
	for i := 1; i < n; i++ {
		a[i][i-1] = 1
	}
	return a, nil
}

// CompanionEigenvalues returns all roots of the polynomial p as the eigenvalues
// of its companion matrix, computed with the Francis double-shift QR algorithm
// on the real Hessenberg form (the classic EISPACK hqr routine). This is a
// robust, allocation-light way to obtain every complex root without complex
// arithmetic in the iteration. Roots are returned sorted by real then imaginary
// part.
func CompanionEigenvalues(p Poly) ([]complex128, error) {
	n := p.Degree()
	if n < 1 {
		return nil, ErrDegreeTooLow
	}
	if n == 1 {
		return []complex128{complex(-p[0]/p[1], 0)}, nil
	}
	// Build 1-based Hessenberg companion matrix for hqr.
	a := make([][]float64, n+1)
	for i := range a {
		a[i] = make([]float64, n+1)
	}
	lc := p[n]
	for j := 1; j <= n; j++ {
		a[1][j] = -p[n-j] / lc
	}
	for i := 2; i <= n; i++ {
		a[i][i-1] = 1
	}
	wr := make([]float64, n+1)
	wi := make([]float64, n+1)
	if err := hqr(a, n, wr, wi); err != nil {
		return nil, err
	}
	roots := make([]complex128, n)
	for i := 1; i <= n; i++ {
		roots[i-1] = complex(wr[i], wi[i])
	}
	sortComplex(roots)
	return roots, nil
}

// sign returns |a| with the sign of b, the classic SIGN(a,b) primitive.
func sign(a, b float64) float64 {
	if b >= 0 {
		return math.Abs(a)
	}
	return -math.Abs(a)
}

// hqr finds all eigenvalues of the real upper-Hessenberg matrix a[1..n][1..n]
// using the Francis double-shift QR algorithm, writing the real and imaginary
// parts of the eigenvalues into wr[1..n] and wi[1..n]. The matrix a is
// overwritten. It is a faithful 1-based port of the EISPACK/Numerical-Recipes
// hqr routine and returns ErrNoConvergence if any eigenvalue needs more than 30
// iterations to deflate.
func hqr(a [][]float64, n int, wr, wi []float64) error {
	var nn, m, l, k, j, its, i, mmin int
	var z, y, x, w, v, u, t, s, r, q, p, anorm float64

	anorm = 0.0
	for i = 1; i <= n; i++ {
		lo := i - 1
		if lo < 1 {
			lo = 1
		}
		for j = lo; j <= n; j++ {
			anorm += math.Abs(a[i][j])
		}
	}
	nn = n
	t = 0.0
	for nn >= 1 {
		its = 0
		for {
			for l = nn; l >= 2; l-- {
				s = math.Abs(a[l-1][l-1]) + math.Abs(a[l][l])
				if s == 0.0 {
					s = anorm
				}
				if math.Abs(a[l][l-1])+s == s {
					a[l][l-1] = 0.0
					break
				}
			}
			x = a[nn][nn]
			if l == nn {
				wr[nn] = x + t
				wi[nn] = 0.0
				nn--
			} else {
				y = a[nn-1][nn-1]
				w = a[nn][nn-1] * a[nn-1][nn]
				if l == nn-1 {
					p = 0.5 * (y - x)
					q = p*p + w
					z = math.Sqrt(math.Abs(q))
					x += t
					if q >= 0.0 {
						z = p + sign(z, p)
						wr[nn-1] = x + z
						wr[nn] = x + z
						if z != 0.0 {
							wr[nn] = x - w/z
						}
						wi[nn-1] = 0.0
						wi[nn] = 0.0
					} else {
						wr[nn-1] = x + p
						wr[nn] = x + p
						wi[nn] = z
						wi[nn-1] = -z
					}
					nn -= 2
				} else {
					if its == 30 {
						return ErrNoConvergence
					}
					if its == 10 || its == 20 {
						t += x
						for i = 1; i <= nn; i++ {
							a[i][i] -= x
						}
						s = math.Abs(a[nn][nn-1]) + math.Abs(a[nn-1][nn-2])
						y = 0.75 * s
						x = y
						w = -0.4375 * s * s
					}
					its++
					for m = nn - 2; m >= l; m-- {
						z = a[m][m]
						r = x - z
						s = y - z
						p = (r*s-w)/a[m+1][m] + a[m][m+1]
						q = a[m+1][m+1] - z - r - s
						r = a[m+2][m+1]
						s = math.Abs(p) + math.Abs(q) + math.Abs(r)
						p /= s
						q /= s
						r /= s
						if m == l {
							break
						}
						u = math.Abs(a[m][m-1]) * (math.Abs(q) + math.Abs(r))
						v = math.Abs(p) * (math.Abs(a[m-1][m-1]) + math.Abs(z) + math.Abs(a[m+1][m+1]))
						if u+v == v {
							break
						}
					}
					for i = m + 2; i <= nn; i++ {
						a[i][i-2] = 0.0
						if i != m+2 {
							a[i][i-3] = 0.0
						}
					}
					for k = m; k <= nn-1; k++ {
						if k != m {
							p = a[k][k-1]
							q = a[k+1][k-1]
							r = 0.0
							if k != nn-1 {
								r = a[k+2][k-1]
							}
							x = math.Abs(p) + math.Abs(q) + math.Abs(r)
							if x != 0.0 {
								p /= x
								q /= x
								r /= x
							}
						}
						s = sign(math.Sqrt(p*p+q*q+r*r), p)
						if s != 0.0 {
							if k == m {
								if l != m {
									a[k][k-1] = -a[k][k-1]
								}
							} else {
								a[k][k-1] = -s * x
							}
							p += s
							x = p / s
							y = q / s
							z = r / s
							q /= p
							r /= p
							for j = k; j <= nn; j++ {
								p = a[k][j] + q*a[k+1][j]
								if k != nn-1 {
									p += r * a[k+2][j]
									a[k+2][j] -= p * z
								}
								a[k+1][j] -= p * y
								a[k][j] -= p * x
							}
							mmin = nn
							if k+3 < nn {
								mmin = k + 3
							}
							for i = l; i <= mmin; i++ {
								p = x*a[i][k] + y*a[i][k+1]
								if k != nn-1 {
									p += z * a[i][k+2]
									a[i][k+2] -= p * r
								}
								a[i][k+1] -= p * q
								a[i][k] -= p
							}
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
