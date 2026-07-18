package controltheory

import "math"

// StateSpace is a SISO continuous-time state-space realization
//
//	x' = A x + B u
//	y  = C x + D u
//
// where A is n×n, B is n×1, C is 1×n, and D is a scalar. Matrices are stored
// as row-major slices of slices.
type StateSpace struct {
	// A is the n×n system (state) matrix.
	A [][]float64
	// B is the n×1 input matrix stored as a column vector of length n.
	B []float64
	// C is the 1×n output matrix stored as a row vector of length n.
	C []float64
	// D is the scalar feedthrough term.
	D float64
}

// NewStateSpace constructs a StateSpace from the given matrices. The slices are
// copied. It panics if the dimensions are inconsistent (A must be square and
// B, C must match its size).
func NewStateSpace(a [][]float64, b, c []float64, d float64) StateSpace {
	n := len(a)
	for _, row := range a {
		if len(row) != n {
			panic("controltheory: A must be square")
		}
	}
	if len(b) != n || len(c) != n {
		panic("controltheory: B and C must have length equal to the order of A")
	}
	ss := StateSpace{
		A: controltheoryCopyMat(a),
		B: append([]float64{}, b...),
		C: append([]float64{}, c...),
		D: d,
	}
	return ss
}

// Order returns the number of states n, the dimension of the A matrix.
func (s StateSpace) Order() int {
	return len(s.A)
}

// controltheoryCopyMat returns a deep copy of a matrix.
func controltheoryCopyMat(m [][]float64) [][]float64 {
	out := make([][]float64, len(m))
	for i, row := range m {
		out[i] = append([]float64{}, row...)
	}
	return out
}

// controltheoryMatVec multiplies matrix m (r×c) by vector v (length c).
func controltheoryMatVec(m [][]float64, v []float64) []float64 {
	out := make([]float64, len(m))
	for i, row := range m {
		var acc float64
		for j, a := range row {
			acc += a * v[j]
		}
		out[i] = acc
	}
	return out
}

// controltheoryMatMul multiplies matrix a (r×k) by matrix b (k×c).
func controltheoryMatMul(a, b [][]float64) [][]float64 {
	r := len(a)
	k := len(b)
	c := 0
	if k > 0 {
		c = len(b[0])
	}
	out := make([][]float64, r)
	for i := 0; i < r; i++ {
		out[i] = make([]float64, c)
		for t := 0; t < k; t++ {
			aij := a[i][t]
			if aij == 0 {
				continue
			}
			for j := 0; j < c; j++ {
				out[i][j] += aij * b[t][j]
			}
		}
	}
	return out
}

// controltheoryIdentity returns the n×n identity matrix.
func controltheoryIdentity(n int) [][]float64 {
	m := make([][]float64, n)
	for i := 0; i < n; i++ {
		m[i] = make([]float64, n)
		m[i][i] = 1
	}
	return m
}

// controltheoryRank returns the numerical rank of matrix m using Gaussian
// elimination with partial pivoting and the given tolerance.
func controltheoryRank(m [][]float64, tol float64) int {
	a := controltheoryCopyMat(m)
	rows := len(a)
	if rows == 0 {
		return 0
	}
	cols := len(a[0])
	rank := 0
	pivotRow := 0
	for col := 0; col < cols && pivotRow < rows; col++ {
		// Find pivot.
		best := pivotRow
		bestAbs := math.Abs(a[pivotRow][col])
		for r := pivotRow + 1; r < rows; r++ {
			if v := math.Abs(a[r][col]); v > bestAbs {
				bestAbs = v
				best = r
			}
		}
		if bestAbs <= tol {
			continue
		}
		a[pivotRow], a[best] = a[best], a[pivotRow]
		pv := a[pivotRow][col]
		for r := 0; r < rows; r++ {
			if r == pivotRow {
				continue
			}
			f := a[r][col] / pv
			if f == 0 {
				continue
			}
			for c := col; c < cols; c++ {
				a[r][c] -= f * a[pivotRow][c]
			}
		}
		rank++
		pivotRow++
	}
	return rank
}

// ControllabilityMatrix returns the controllability matrix
// [B, AB, A^2 B, ..., A^(n-1) B], an n×n matrix for a SISO system.
func (s StateSpace) ControllabilityMatrix() [][]float64 {
	n := len(s.A)
	cols := make([][]float64, n) // cols[k] is the k-th column vector
	v := append([]float64{}, s.B...)
	for k := 0; k < n; k++ {
		cols[k] = append([]float64{}, v...)
		v = controltheoryMatVec(s.A, v)
	}
	// Assemble into row-major n×n matrix.
	out := make([][]float64, n)
	for i := 0; i < n; i++ {
		out[i] = make([]float64, n)
		for k := 0; k < n; k++ {
			out[i][k] = cols[k][i]
		}
	}
	return out
}

// ObservabilityMatrix returns the observability matrix
// [C; CA; CA^2; ...; CA^(n-1)], an n×n matrix for a SISO system.
func (s StateSpace) ObservabilityMatrix() [][]float64 {
	n := len(s.A)
	out := make([][]float64, n)
	row := append([]float64{}, s.C...)
	for k := 0; k < n; k++ {
		out[k] = append([]float64{}, row...)
		// next row = row * A  (1×n times n×n).
		next := make([]float64, n)
		for j := 0; j < n; j++ {
			var acc float64
			for i := 0; i < n; i++ {
				acc += row[i] * s.A[i][j]
			}
			next[j] = acc
		}
		row = next
	}
	return out
}

// ControllabilityRank returns the rank of the controllability matrix.
func (s StateSpace) ControllabilityRank() int {
	return controltheoryRank(s.ControllabilityMatrix(), 1e-9)
}

// ObservabilityRank returns the rank of the observability matrix.
func (s StateSpace) ObservabilityRank() int {
	return controltheoryRank(s.ObservabilityMatrix(), 1e-9)
}

// IsControllable reports whether the system is completely state controllable,
// i.e. the controllability matrix has full rank n.
func (s StateSpace) IsControllable() bool {
	return s.ControllabilityRank() == len(s.A)
}

// IsObservable reports whether the system is completely observable, i.e. the
// observability matrix has full rank n.
func (s StateSpace) IsObservable() bool {
	return s.ObservabilityRank() == len(s.A)
}

// CharacteristicPolynomial returns the characteristic polynomial det(sI - A) of
// the state matrix as an ascending-power monic [Poly], computed with the
// Faddeev-LeVerrier algorithm.
func (s StateSpace) CharacteristicPolynomial() Poly {
	den, _ := controltheoryFaddeevLeverrier(s.A)
	return den
}

// Poles returns the poles of the realization, i.e. the eigenvalues of A, found
// as the roots of the characteristic polynomial.
func (s StateSpace) Poles() []complex128 {
	return s.CharacteristicPolynomial().Roots()
}

// controltheoryFaddeevLeverrier runs the Faddeev-LeVerrier algorithm on the
// n×n matrix a. It returns the characteristic polynomial det(sI-A) as an
// ascending-power monic Poly of degree n, and the sequence of adjugate
// coefficient matrices R_0..R_{n-1} for which
// adj(sI-A) = sum_{k=0}^{n-1} R_k s^{n-1-k}.
func controltheoryFaddeevLeverrier(a [][]float64) (Poly, [][][]float64) {
	n := len(a)
	// Descending coefficients c[0..n] with c[0] = 1 (monic).
	c := make([]float64, n+1)
	c[0] = 1
	rMats := make([][][]float64, n) // R_0 .. R_{n-1}
	R := controltheoryIdentity(n)   // R_0 = I
	for k := 1; k <= n; k++ {
		rMats[k-1] = controltheoryCopyMat(R)
		AR := controltheoryMatMul(a, R)
		// trace(AR)
		var tr float64
		for i := 0; i < n; i++ {
			tr += AR[i][i]
		}
		ck := -tr / float64(k)
		c[k] = ck
		if k < n {
			// R = A*R + ck*I  (next adjugate coefficient matrix).
			R = AR
			for i := 0; i < n; i++ {
				R[i][i] += ck
			}
		}
	}
	// Convert descending c[0..n] (coeff of s^{n-k}) to ascending Poly.
	p := make(Poly, n+1)
	for k := 0; k <= n; k++ {
		p[n-k] = c[k]
	}
	return p, rMats
}

// TransferFunctionToStateSpace converts a proper transfer function into a
// state-space realization in controllable canonical form. Any direct
// feedthrough (when the numerator and denominator have equal degree) is placed
// in the D term. It panics if the transfer function is not proper or the
// denominator is the zero polynomial.
func TransferFunctionToStateSpace(g TransferFunction) StateSpace {
	den := controltheoryTrim(g.Den)
	if len(den) == 0 {
		panic("controltheory: denominator is the zero polynomial")
	}
	if !g.IsProper() {
		panic("controltheory: transfer function is not proper")
	}
	n := len(den) - 1
	// Normalize so denominator is monic.
	lead := den[n]
	a := make(Poly, n+1)
	for i := range a {
		a[i] = den[i] / lead
	}
	num := make(Poly, n+1)
	src := controltheoryTrim(g.Num)
	for i := 0; i < len(src) && i <= n; i++ {
		num[i] = src[i] / lead
	}
	// Extract feedthrough: num = d*a + strictlyProperNum.
	d := num[n] // coefficient of s^n in numerator (0 if strictly proper)
	b := make([]float64, n)
	for i := 0; i < n; i++ {
		b[i] = num[i] - d*a[i]
	}
	if n == 0 {
		return StateSpace{A: [][]float64{}, B: []float64{}, C: []float64{}, D: d}
	}
	// Controllable canonical form.
	A := make([][]float64, n)
	for i := 0; i < n; i++ {
		A[i] = make([]float64, n)
	}
	for i := 0; i < n-1; i++ {
		A[i][i+1] = 1
	}
	for j := 0; j < n; j++ {
		A[n-1][j] = -a[j]
	}
	B := make([]float64, n)
	B[n-1] = 1
	C := make([]float64, n)
	copy(C, b[:n])
	return StateSpace{A: A, B: B, C: C, D: d}
}

// TransferFunction converts the state-space realization to an equivalent SISO
// transfer function G(s) = C(sI-A)^{-1}B + D using the Faddeev-LeVerrier
// algorithm.
func (s StateSpace) TransferFunction() TransferFunction {
	n := len(s.A)
	if n == 0 {
		return TransferFunction{Num: Poly{s.D}, Den: Poly{1}}
	}
	den, rMats := controltheoryFaddeevLeverrier(s.A)
	// den is ascending, degree n, monic. Recover descending c[0..n].
	c := make([]float64, n+1)
	for k := 0; k <= n; k++ {
		c[k] = den[n-k]
	}
	// Numerator descending: coeff of s^{n-1-k} is C R_k B, for k=0..n-1,
	// plus D * den.
	numDesc := make([]float64, n+1) // index j is coeff of s^{n-j}
	for k := 0; k < n; k++ {
		// C * R_k * B  (1×n * n×n * n×1).
		rb := controltheoryMatVec(rMats[k], s.B)
		var val float64
		for i := 0; i < n; i++ {
			val += s.C[i] * rb[i]
		}
		numDesc[k+1] += val // coeff of s^{n-1-k} = s^{n-(k+1)}
	}
	for j := 0; j <= n; j++ {
		numDesc[j] += s.D * c[j]
	}
	// Convert to ascending.
	num := make(Poly, n+1)
	for j := 0; j <= n; j++ {
		num[n-j] = numDesc[j]
	}
	return TransferFunction{Num: controltheoryTrim(num), Den: controltheoryTrim(den)}
}
