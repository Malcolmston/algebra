package optimalcontrol

import (
	"math"
	"math/cmplx"
)

// LU holds the result of an LU factorization with partial pivoting: a combined
// lower/upper factor lu, the row permutation piv, and the sign of the
// permutation (+1 or -1) used for the determinant.
type LU struct {
	lu   *Matrix
	piv  []int
	sign float64
	n    int
}

// Factor computes the LU decomposition of a square matrix with partial
// pivoting. It returns ErrDim for non-square input.
func Factor(a *Matrix) (*LU, error) {
	if !a.IsSquare() {
		return nil, ErrDim
	}
	n := a.rows
	lu := a.Clone()
	piv := make([]int, n)
	for i := range piv {
		piv[i] = i
	}
	sign := 1.0
	for k := 0; k < n; k++ {
		// Pivot selection.
		p := k
		max := math.Abs(lu.At(k, k))
		for i := k + 1; i < n; i++ {
			if v := math.Abs(lu.At(i, k)); v > max {
				max = v
				p = i
			}
		}
		if p != k {
			for j := 0; j < n; j++ {
				tmp := lu.At(k, j)
				lu.Set(k, j, lu.At(p, j))
				lu.Set(p, j, tmp)
			}
			piv[k], piv[p] = piv[p], piv[k]
			sign = -sign
		}
		akk := lu.At(k, k)
		if akk == 0 {
			continue
		}
		for i := k + 1; i < n; i++ {
			f := lu.At(i, k) / akk
			lu.Set(i, k, f)
			for j := k + 1; j < n; j++ {
				lu.Set(i, j, lu.At(i, j)-f*lu.At(k, j))
			}
		}
	}
	return &LU{lu: lu, piv: piv, sign: sign, n: n}, nil
}

// Det returns the determinant computed from the factorization.
func (f *LU) Det() float64 {
	d := f.sign
	for i := 0; i < f.n; i++ {
		d *= f.lu.At(i, i)
	}
	return d
}

// Solve solves A x = b for x using the factorization.
func (f *LU) Solve(b []float64) ([]float64, error) {
	if len(b) != f.n {
		return nil, ErrDim
	}
	for i := 0; i < f.n; i++ {
		if f.lu.At(i, i) == 0 {
			return nil, ErrSingular
		}
	}
	n := f.n
	x := make([]float64, n)
	for i := 0; i < n; i++ {
		x[i] = b[f.piv[i]]
	}
	// Forward substitution (unit lower).
	for i := 0; i < n; i++ {
		for j := 0; j < i; j++ {
			x[i] -= f.lu.At(i, j) * x[j]
		}
	}
	// Back substitution (upper).
	for i := n - 1; i >= 0; i-- {
		for j := i + 1; j < n; j++ {
			x[i] -= f.lu.At(i, j) * x[j]
		}
		x[i] /= f.lu.At(i, i)
	}
	return x, nil
}

// SolveMatrix solves A X = B, treating each column of B independently.
func (f *LU) SolveMatrix(b *Matrix) (*Matrix, error) {
	if b.rows != f.n {
		return nil, ErrDim
	}
	out := Zeros(f.n, b.cols)
	for j := 0; j < b.cols; j++ {
		x, err := f.Solve(b.Col(j))
		if err != nil {
			return nil, err
		}
		out.SetCol(j, x)
	}
	return out, nil
}

// Det returns the determinant of a square matrix.
func Det(a *Matrix) (float64, error) {
	f, err := Factor(a)
	if err != nil {
		return 0, err
	}
	return f.Det(), nil
}

// Solve solves the linear system A x = b.
func Solve(a *Matrix, b []float64) ([]float64, error) {
	f, err := Factor(a)
	if err != nil {
		return nil, err
	}
	return f.Solve(b)
}

// SolveMatrix solves A X = B for the matrix X.
func SolveMatrix(a, b *Matrix) (*Matrix, error) {
	f, err := Factor(a)
	if err != nil {
		return nil, err
	}
	return f.SolveMatrix(b)
}

// Inverse returns the inverse of a square matrix.
func Inverse(a *Matrix) (*Matrix, error) {
	f, err := Factor(a)
	if err != nil {
		return nil, err
	}
	return f.SolveMatrix(Eye(a.rows))
}

// Cholesky computes the lower-triangular Cholesky factor L with A = L Lᵀ for a
// symmetric positive-definite matrix. It returns an error if A is not SPD.
func Cholesky(a *Matrix) (*Matrix, error) {
	if !a.IsSquare() {
		return nil, ErrDim
	}
	n := a.rows
	l := Zeros(n, n)
	for i := 0; i < n; i++ {
		for j := 0; j <= i; j++ {
			s := a.At(i, j)
			for k := 0; k < j; k++ {
				s -= l.At(i, k) * l.At(j, k)
			}
			if i == j {
				if s <= 0 {
					return nil, ErrSingular
				}
				l.Set(i, j, math.Sqrt(s))
			} else {
				l.Set(i, j, s/l.At(j, j))
			}
		}
	}
	return l, nil
}

// IsPositiveDefinite reports whether a symmetric matrix is positive definite.
func IsPositiveDefinite(a *Matrix) bool {
	_, err := Cholesky(a)
	return err == nil
}

// IsPositiveSemiDefinite reports whether a symmetric matrix has all
// eigenvalues >= -tol.
func IsPositiveSemiDefinite(a *Matrix, tol float64) bool {
	ev, err := SymEigenvalues(a)
	if err != nil {
		return false
	}
	for _, x := range ev {
		if x < -tol {
			return false
		}
	}
	return true
}

// LeastSquares solves the overdetermined system A x = b in the least-squares
// sense via the normal equations AᵀA x = Aᵀb. A must have full column rank.
func LeastSquares(a *Matrix, b []float64) ([]float64, error) {
	at := a.Transpose()
	ata := at.Mul(a)
	atb := at.MulVec(b)
	return Solve(ata, atb)
}

// LeastSquaresMatrix solves A X = B in the least-squares sense column by column.
func LeastSquaresMatrix(a, b *Matrix) (*Matrix, error) {
	at := a.Transpose()
	ata := at.Mul(a)
	atb := at.Mul(b)
	return SolveMatrix(ata, atb)
}

// MatrixPow returns A raised to the non-negative integer power p.
func MatrixPow(a *Matrix, p int) *Matrix {
	if !a.IsSquare() {
		panic(ErrDim)
	}
	if p < 0 {
		panic("optimalcontrol: negative power")
	}
	result := Eye(a.rows)
	base := a.Clone()
	for p > 0 {
		if p&1 == 1 {
			result = result.Mul(base)
		}
		base = base.Mul(base)
		p >>= 1
	}
	return result
}

// MatrixExp returns the matrix exponential exp(A) using a scaling-and-squaring
// scheme with a truncated Taylor series. It is accurate for the moderate-sized
// matrices used in control applications.
func MatrixExp(a *Matrix) *Matrix {
	if !a.IsSquare() {
		panic(ErrDim)
	}
	n := a.rows
	// Scaling: choose s so that ||A/2^s|| is small.
	norm := a.InfNorm()
	s := 0
	for norm > 0.5 {
		norm /= 2
		s++
	}
	scaled := a.Scale(1.0 / math.Pow(2, float64(s)))
	// Taylor series.
	result := Eye(n)
	term := Eye(n)
	for k := 1; k <= 20; k++ {
		term = term.Mul(scaled).Scale(1.0 / float64(k))
		result = result.Plus(term)
		if term.MaxAbs() < 1e-18 {
			break
		}
	}
	// Squaring.
	for i := 0; i < s; i++ {
		result = result.Mul(result)
	}
	return result
}

// CharPoly returns the coefficients of the characteristic polynomial
// det(λI − A) of a square matrix, in ascending powers of λ (index i holds the
// coefficient of λ^i). It uses the Faddeev–LeVerrier algorithm.
func CharPoly(a *Matrix) []float64 {
	n := a.rows
	if !a.IsSquare() {
		panic(ErrDim)
	}
	// coeffs holds the char-poly coefficients in ascending powers; the leading
	// coefficient (λ^n) is 1. The Faddeev–LeVerrier recurrence is
	//   M_1 = A,               c_1 = -tr(M_1),
	//   M_k = A (M_{k-1} + c_{k-1} I),  c_k = -tr(M_k)/k.
	coeffs := make([]float64, n+1)
	if n == 0 {
		coeffs[0] = 1
		return coeffs
	}
	coeffs[n] = 1
	M := a.Clone()
	c := -M.Trace()
	coeffs[n-1] = c
	for k := 2; k <= n; k++ {
		tmp := M.Clone()
		for i := 0; i < n; i++ {
			tmp.Add(i, i, c)
		}
		M = a.Mul(tmp)
		c = -M.Trace() / float64(k)
		coeffs[n-k] = c
	}
	return coeffs
}

// Eigenvalues returns all (generally complex) eigenvalues of a square matrix by
// forming its characteristic polynomial and finding its roots with the
// Durand–Kerner method. Results are not returned in any particular order.
func Eigenvalues(a *Matrix) []complex128 {
	coeffs := CharPoly(a)
	return PolyRootsDK(coeffs)
}

// SpectralRadius returns the maximum modulus of the eigenvalues of a matrix.
func SpectralRadius(a *Matrix) float64 {
	ev := Eigenvalues(a)
	var mx float64
	for _, z := range ev {
		if m := cmplx.Abs(z); m > mx {
			mx = m
		}
	}
	return mx
}

// SpectralAbscissa returns the maximum real part of the eigenvalues of a matrix.
func SpectralAbscissa(a *Matrix) float64 {
	ev := Eigenvalues(a)
	mx := math.Inf(-1)
	for _, z := range ev {
		if r := real(z); r > mx {
			mx = r
		}
	}
	return mx
}

// IsStableContinuous reports whether every eigenvalue of A has real part below
// -tol (Hurwitz stability of x' = A x).
func IsStableContinuous(a *Matrix, tol float64) bool {
	for _, z := range Eigenvalues(a) {
		if real(z) >= -tol {
			return false
		}
	}
	return true
}

// IsStableDiscrete reports whether every eigenvalue of A has modulus below
// 1-tol (Schur stability of x_{k+1} = A x_k).
func IsStableDiscrete(a *Matrix, tol float64) bool {
	for _, z := range Eigenvalues(a) {
		if cmplx.Abs(z) >= 1-tol {
			return false
		}
	}
	return true
}

// PolyRootsDK finds all complex roots of a real polynomial using the
// Durand–Kerner (Weierstrass) iteration. Coefficients are given in ascending
// power order: coeffs[i] multiplies x^i. The leading coefficient must be
// nonzero.
func PolyRootsDK(coeffs []float64) []complex128 {
	// Trim trailing (highest-order) zeros.
	deg := len(coeffs) - 1
	for deg > 0 && coeffs[deg] == 0 {
		deg--
	}
	if deg <= 0 {
		return nil
	}
	// Normalize to monic ascending coefficients of length deg+1.
	lead := coeffs[deg]
	a := make([]complex128, deg+1)
	for i := 0; i <= deg; i++ {
		a[i] = complex(coeffs[i]/lead, 0)
	}
	eval := func(x complex128) complex128 {
		// Horner in descending order.
		res := a[deg]
		for i := deg - 1; i >= 0; i-- {
			res = res*x + a[i]
		}
		return res
	}
	// Initial guesses on a spiral.
	roots := make([]complex128, deg)
	seed := complex(0.4, 0.9)
	roots[0] = complex(1, 0)
	for i := 1; i < deg; i++ {
		roots[i] = roots[i-1] * seed
	}
	for iter := 0; iter < 500; iter++ {
		var maxDelta float64
		for i := 0; i < deg; i++ {
			num := eval(roots[i])
			den := complex(1, 0)
			for j := 0; j < deg; j++ {
				if j != i {
					den *= roots[i] - roots[j]
				}
			}
			if den == 0 {
				continue
			}
			delta := num / den
			roots[i] -= delta
			if d := cmplx.Abs(delta); d > maxDelta {
				maxDelta = d
			}
		}
		if maxDelta < 1e-14 {
			break
		}
	}
	// Clean tiny imaginary parts.
	for i := range roots {
		if math.Abs(imag(roots[i])) < 1e-10 {
			roots[i] = complex(real(roots[i]), 0)
		}
	}
	return roots
}

// jacobiRotate performs one symmetric Jacobi eigenvalue sweep helper: not
// exported.
func maxOffDiag(a *Matrix) (int, int, float64) {
	n := a.rows
	p, q := 0, 1
	mx := 0.0
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			if v := math.Abs(a.At(i, j)); v > mx {
				mx = v
				p, q = i, j
			}
		}
	}
	return p, q, mx
}

// SymEigen computes eigenvalues and orthonormal eigenvectors of a symmetric
// matrix using the cyclic Jacobi method. The eigenvalues are returned as a
// slice and the eigenvectors as the columns of V, so that A V = V diag(w).
func SymEigen(a *Matrix) (w []float64, v *Matrix, err error) {
	if !a.IsSquare() {
		return nil, nil, ErrDim
	}
	n := a.rows
	d := a.Symmetrize()
	V := Eye(n)
	for sweep := 0; sweep < 100; sweep++ {
		p, q, off := maxOffDiag(d)
		if off < 1e-300 || off < 1e-15*d.FrobeniusNorm() {
			break
		}
		app := d.At(p, p)
		aqq := d.At(q, q)
		apq := d.At(p, q)
		phi := 0.5 * math.Atan2(2*apq, aqq-app)
		c := math.Cos(phi)
		s := math.Sin(phi)
		for i := 0; i < n; i++ {
			dip := d.At(i, p)
			diq := d.At(i, q)
			d.Set(i, p, c*dip-s*diq)
			d.Set(i, q, s*dip+c*diq)
		}
		for i := 0; i < n; i++ {
			dpi := d.At(p, i)
			dqi := d.At(q, i)
			d.Set(p, i, c*dpi-s*dqi)
			d.Set(q, i, s*dpi+c*dqi)
		}
		for i := 0; i < n; i++ {
			vip := V.At(i, p)
			viq := V.At(i, q)
			V.Set(i, p, c*vip-s*viq)
			V.Set(i, q, s*vip+c*viq)
		}
	}
	w = make([]float64, n)
	for i := 0; i < n; i++ {
		w[i] = d.At(i, i)
	}
	return w, V, nil
}

// SymEigenvalues returns the eigenvalues of a symmetric matrix in ascending
// order.
func SymEigenvalues(a *Matrix) ([]float64, error) {
	w, _, err := SymEigen(a)
	if err != nil {
		return nil, err
	}
	sortFloats(w)
	return w, nil
}

// sortFloats sorts a slice ascending in place (insertion sort, small n).
func sortFloats(w []float64) {
	for i := 1; i < len(w); i++ {
		v := w[i]
		j := i - 1
		for j >= 0 && w[j] > v {
			w[j+1] = w[j]
			j--
		}
		w[j+1] = v
	}
}
