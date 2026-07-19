package liealgebra

import "math"

// Bracket returns the Lie bracket (commutator) [A,B] = AB - BA of two square
// real matrices. It returns [ErrDim] if the shapes are incompatible.
func Bracket(a, b *Matrix) (*Matrix, error) {
	ab, err := a.Mul(b)
	if err != nil {
		return nil, err
	}
	ba, err := b.Mul(a)
	if err != nil {
		return nil, err
	}
	return ab.Sub(ba)
}

// Anticommutator returns {A,B} = AB + BA of two square real matrices.
func Anticommutator(a, b *Matrix) (*Matrix, error) {
	ab, err := a.Mul(b)
	if err != nil {
		return nil, err
	}
	ba, err := b.Mul(a)
	if err != nil {
		return nil, err
	}
	return ab.Add(ba)
}

// CBracket returns the Lie bracket [A,B] = AB - BA of two complex matrices.
func CBracket(a, b *CMatrix) (*CMatrix, error) {
	ab, err := a.Mul(b)
	if err != nil {
		return nil, err
	}
	ba, err := b.Mul(a)
	if err != nil {
		return nil, err
	}
	return ab.Sub(ba)
}

// CAnticommutator returns {A,B} = AB + BA of two complex matrices.
func CAnticommutator(a, b *CMatrix) (*CMatrix, error) {
	ab, err := a.Mul(b)
	if err != nil {
		return nil, err
	}
	ba, err := b.Mul(a)
	if err != nil {
		return nil, err
	}
	return ab.Add(ba)
}

// AdjointAction returns ad_X(Y) = [X,Y], the action of X in the adjoint
// representation on the matrix Y.
func AdjointAction(x, y *Matrix) (*Matrix, error) { return Bracket(x, y) }

// NestedBracket returns the left-nested bracket [x0,[x1,[...,xk]...]] for a list
// of matrices. An empty list returns [ErrDim]; a single element returns it.
func NestedBracket(xs ...*Matrix) (*Matrix, error) {
	if len(xs) == 0 {
		return nil, ErrDim
	}
	acc := xs[len(xs)-1].Clone()
	for i := len(xs) - 2; i >= 0; i-- {
		var err error
		acc, err = Bracket(xs[i], acc)
		if err != nil {
			return nil, err
		}
	}
	return acc, nil
}

// JacobiResidual returns [[A,B],C] + [[B,C],A] + [[C,A],B], which is identically
// the zero matrix for any associative product and hence for matrix Lie
// algebras. It is a convenient numerical check of the Jacobi identity.
func JacobiResidual(a, b, c *Matrix) (*Matrix, error) {
	ab, err := Bracket(a, b)
	if err != nil {
		return nil, err
	}
	t1, err := Bracket(ab, c)
	if err != nil {
		return nil, err
	}
	bc, err := Bracket(b, c)
	if err != nil {
		return nil, err
	}
	t2, err := Bracket(bc, a)
	if err != nil {
		return nil, err
	}
	ca, err := Bracket(c, a)
	if err != nil {
		return nil, err
	}
	t3, err := Bracket(ca, b)
	if err != nil {
		return nil, err
	}
	s, err := t1.Add(t2)
	if err != nil {
		return nil, err
	}
	return s.Add(t3)
}

// SatisfiesJacobi reports whether the Jacobi identity holds for A, B, C to
// within tolerance tol (measured by the max absolute entry of the residual).
func SatisfiesJacobi(a, b, c *Matrix, tol float64) bool {
	r, err := JacobiResidual(a, b, c)
	if err != nil {
		return false
	}
	return r.MaxAbs() <= tol
}

// TraceForm returns the trace form tr(AB) of two square matrices. For many
// matrix Lie algebras this is proportional to the Killing form.
func TraceForm(a, b *Matrix) (float64, error) {
	ab, err := a.Mul(b)
	if err != nil {
		return 0, err
	}
	return ab.Trace()
}

// sameShape reports whether all matrices in a basis share one square shape.
func sameShape(basis []*Matrix) (int, bool) {
	if len(basis) == 0 {
		return 0, false
	}
	n := basis[0].Rows
	for _, m := range basis {
		if m.Rows != n || m.Cols != n {
			return 0, false
		}
	}
	return n, true
}

// StructureConstants computes the structure constants c[k][i][j] of a matrix
// Lie algebra with the given basis, defined by [e_i, e_j] = Σ_k c[k][i][j] e_k.
// The basis matrices must be square, of equal size and linearly independent;
// otherwise [ErrRank] or [ErrDim] is returned. The bracket of any two basis
// elements must lie in the span of the basis (a genuine subalgebra).
func StructureConstants(basis []*Matrix) ([][][]float64, error) {
	n, ok := sameShape(basis)
	if !ok {
		return nil, ErrDim
	}
	d := len(basis)
	// Build the n*n by d matrix whose columns are vec(e_k).
	m := NewMatrix(n*n, d)
	for k := 0; k < d; k++ {
		for r := 0; r < n*n; r++ {
			m.Data[r*d+k] = basis[k].Data[r]
		}
	}
	if Rank(m, 1e-9) != d {
		return nil, ErrRank
	}
	c := make([][][]float64, d)
	for k := range c {
		c[k] = make([][]float64, d)
		for i := range c[k] {
			c[k][i] = make([]float64, d)
		}
	}
	for i := 0; i < d; i++ {
		for j := 0; j < d; j++ {
			br, err := Bracket(basis[i], basis[j])
			if err != nil {
				return nil, err
			}
			coeff, err := SolveLeastSquares(m, br.Data)
			if err != nil {
				return nil, err
			}
			for k := 0; k < d; k++ {
				c[k][i][j] = coeff[k]
			}
		}
	}
	return c, nil
}

// IsClosedUnderBracket reports whether the bracket of every pair of basis
// elements lies in the span of the basis to within tolerance tol, i.e. whether
// the basis spans a Lie subalgebra.
func IsClosedUnderBracket(basis []*Matrix, tol float64) bool {
	n, ok := sameShape(basis)
	if !ok {
		return false
	}
	d := len(basis)
	m := NewMatrix(n*n, d)
	for k := 0; k < d; k++ {
		for r := 0; r < n*n; r++ {
			m.Data[r*d+k] = basis[k].Data[r]
		}
	}
	for i := 0; i < d; i++ {
		for j := i + 1; j < d; j++ {
			br, err := Bracket(basis[i], basis[j])
			if err != nil {
				return false
			}
			coeff, err := SolveLeastSquares(m, br.Data)
			if err != nil {
				return false
			}
			recon, _ := m.MatVec(coeff)
			maxd := 0.0
			for r := range recon {
				if v := math.Abs(recon[r] - br.Data[r]); v > maxd {
					maxd = v
				}
			}
			if maxd > tol {
				return false
			}
		}
	}
	return true
}

// AdjointMatrix returns the matrix of ad_{e_index} in the given basis: the
// d-by-d matrix M with M[k][j] = c[k][index][j] where c are the structure
// constants. Column j is the coordinate vector of [e_index, e_j].
func AdjointMatrix(basis []*Matrix, index int) (*Matrix, error) {
	if index < 0 || index >= len(basis) {
		return nil, ErrRange
	}
	c, err := StructureConstants(basis)
	if err != nil {
		return nil, err
	}
	d := len(basis)
	m := NewMatrix(d, d)
	for k := 0; k < d; k++ {
		for j := 0; j < d; j++ {
			m.Data[k*d+j] = c[k][index][j]
		}
	}
	return m, nil
}

// KillingForm returns the Killing form matrix K of a Lie algebra given by a
// basis, with entries K[i][j] = tr(ad_{e_i} ad_{e_j}). The result is symmetric.
func KillingForm(basis []*Matrix) (*Matrix, error) {
	c, err := StructureConstants(basis)
	if err != nil {
		return nil, err
	}
	d := len(basis)
	// Precompute adjoint matrices ad_i[k][j] = c[k][i][j].
	ad := make([]*Matrix, d)
	for i := 0; i < d; i++ {
		m := NewMatrix(d, d)
		for k := 0; k < d; k++ {
			for j := 0; j < d; j++ {
				m.Data[k*d+j] = c[k][i][j]
			}
		}
		ad[i] = m
	}
	k := NewMatrix(d, d)
	for i := 0; i < d; i++ {
		for j := i; j < d; j++ {
			prod, _ := ad[i].Mul(ad[j])
			tr, _ := prod.Trace()
			k.Data[i*d+j] = tr
			k.Data[j*d+i] = tr
		}
	}
	return k, nil
}

// KillingFormValue returns the Killing form value K(X,Y) = tr(ad_X ad_Y) for two
// elements X and Y given by their coordinate vectors in the supplied basis.
func KillingFormValue(basis []*Matrix, x, y []float64) (float64, error) {
	k, err := KillingForm(basis)
	if err != nil {
		return 0, err
	}
	d := len(basis)
	if len(x) != d || len(y) != d {
		return 0, ErrDim
	}
	s := 0.0
	for i := 0; i < d; i++ {
		for j := 0; j < d; j++ {
			s += x[i] * k.Data[i*d+j] * y[j]
		}
	}
	return s, nil
}

// IsSemisimple reports whether the Killing form of the given basis is
// nondegenerate (Cartan's criterion), tested by |det K| > tol.
func IsSemisimple(basis []*Matrix, tol float64) (bool, error) {
	k, err := KillingForm(basis)
	if err != nil {
		return false, err
	}
	det, err := Det(k)
	if err != nil {
		return false, err
	}
	return math.Abs(det) > tol, nil
}

// JacobiResidualConstants returns the maximum absolute value of the Jacobi
// identity applied to structure constants:
// Σ_m ( c[m][i][j] c[n][m][k] + c[m][j][k] c[n][m][i] + c[m][k][i] c[n][m][j] ).
// It is zero for a genuine Lie algebra and provides a coordinate-space Jacobi
// check independent of the matrix realisation.
func JacobiResidualConstants(c [][][]float64) float64 {
	d := len(c)
	max := 0.0
	for i := 0; i < d; i++ {
		for j := 0; j < d; j++ {
			for k := 0; k < d; k++ {
				for n := 0; n < d; n++ {
					s := 0.0
					for m := 0; m < d; m++ {
						s += c[m][i][j]*c[n][m][k] +
							c[m][j][k]*c[n][m][i] +
							c[m][k][i]*c[n][m][j]
					}
					if a := math.Abs(s); a > max {
						max = a
					}
				}
			}
		}
	}
	return max
}
