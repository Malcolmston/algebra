package liealgebra

import (
	"math"
	"math/cmplx"
)

// expTaylorOrder is the truncation order of the Taylor series used after
// scaling. With the argument scaled to spectral-ish norm below 1/2 this is far
// beyond double-precision accuracy.
const expTaylorOrder = 20

// oneNorm returns the maximum absolute column sum of a real matrix.
func oneNorm(m *Matrix) float64 {
	max := 0.0
	for j := 0; j < m.Cols; j++ {
		s := 0.0
		for i := 0; i < m.Rows; i++ {
			s += math.Abs(m.Data[i*m.Cols+j])
		}
		if s > max {
			max = s
		}
	}
	return max
}

// MatExp returns the matrix exponential exp(A) of a square real matrix, computed
// by scaling and squaring with a truncated Taylor series. It returns
// [ErrNotSquare] for non-square input.
func MatExp(a *Matrix) (*Matrix, error) {
	if !a.IsSquare() {
		return nil, ErrNotSquare
	}
	n := a.Rows
	norm := oneNorm(a)
	s := 0
	for norm/math.Pow(2, float64(s)) > 0.5 {
		s++
	}
	scale := 1.0 / math.Pow(2, float64(s))
	b := a.Scale(scale)
	result := IdentityMatrix(n)
	term := IdentityMatrix(n)
	for k := 1; k <= expTaylorOrder; k++ {
		t, _ := term.Mul(b)
		term = t.Scale(1.0 / float64(k))
		result, _ = result.Add(term)
	}
	for ; s > 0; s-- {
		result, _ = result.Mul(result)
	}
	return result, nil
}

// ExpMap is an alias for [MatExp]; for a Lie-algebra element X it returns the
// corresponding group element exp(X).
func ExpMap(a *Matrix) (*Matrix, error) { return MatExp(a) }

// MatPow returns A raised to a nonnegative integer power using exponentiation by
// squaring; A^0 is the identity. It returns [ErrNotSquare] for non-square input
// and [ErrRange] for a negative power.
func MatPow(a *Matrix, p int) (*Matrix, error) {
	if !a.IsSquare() {
		return nil, ErrNotSquare
	}
	if p < 0 {
		return nil, ErrRange
	}
	result := IdentityMatrix(a.Rows)
	base := a.Clone()
	for p > 0 {
		if p&1 == 1 {
			result, _ = result.Mul(base)
		}
		p >>= 1
		if p > 0 {
			base, _ = base.Mul(base)
		}
	}
	return result, nil
}

// oneNormC returns the maximum absolute column sum of a complex matrix.
func oneNormC(m *CMatrix) float64 {
	max := 0.0
	for j := 0; j < m.Cols; j++ {
		s := 0.0
		for i := 0; i < m.Rows; i++ {
			s += cmplx.Abs(m.Data[i*m.Cols+j])
		}
		if s > max {
			max = s
		}
	}
	return max
}

// CMatExp returns the matrix exponential exp(A) of a square complex matrix by
// scaling and squaring with a truncated Taylor series.
func CMatExp(a *CMatrix) (*CMatrix, error) {
	if !a.IsSquare() {
		return nil, ErrNotSquare
	}
	n := a.Rows
	norm := oneNormC(a)
	s := 0
	for norm/math.Pow(2, float64(s)) > 0.5 {
		s++
	}
	scale := complex(1.0/math.Pow(2, float64(s)), 0)
	b := a.Scale(scale)
	result := IdentityCMatrix(n)
	term := IdentityCMatrix(n)
	for k := 1; k <= expTaylorOrder; k++ {
		t, _ := term.Mul(b)
		term = t.Scale(complex(1.0/float64(k), 0))
		result, _ = result.Add(term)
	}
	for ; s > 0; s-- {
		result, _ = result.Mul(result)
	}
	return result, nil
}

// CMatPow returns a complex matrix raised to a nonnegative integer power.
func CMatPow(a *CMatrix, p int) (*CMatrix, error) {
	if !a.IsSquare() {
		return nil, ErrNotSquare
	}
	if p < 0 {
		return nil, ErrRange
	}
	result := IdentityCMatrix(a.Rows)
	base := a.Clone()
	for p > 0 {
		if p&1 == 1 {
			result, _ = result.Mul(base)
		}
		p >>= 1
		if p > 0 {
			base, _ = base.Mul(base)
		}
	}
	return result, nil
}

// BCHApprox returns the second-order Baker-Campbell-Hausdorff approximation to
// log(exp(X)exp(Y)):
//
//	X + Y + ½[X,Y] + (1/12)([X,[X,Y]] + [Y,[Y,X]]).
//
// It is exact through third order in X and Y and is the standard truncation for
// small elements of a matrix Lie algebra.
func BCHApprox(x, y *Matrix) (*Matrix, error) {
	sum, err := x.Add(y)
	if err != nil {
		return nil, err
	}
	xy, err := Bracket(x, y)
	if err != nil {
		return nil, err
	}
	sum, err = sum.Add(xy.Scale(0.5))
	if err != nil {
		return nil, err
	}
	xxy, err := NestedBracket(x, x, y) // [x,[x,y]]
	if err != nil {
		return nil, err
	}
	yyx, err := NestedBracket(y, y, x) // [y,[y,x]]
	if err != nil {
		return nil, err
	}
	third, err := xxy.Add(yyx)
	if err != nil {
		return nil, err
	}
	return sum.Add(third.Scale(1.0 / 12.0))
}

// ExpBracketSeries returns the truncated series for Ad_{exp(X)}(Y) =
// exp(ad_X)(Y) = Σ_{k≥0} (1/k!) ad_X^k(Y), summed to the given number of terms
// (terms>=1). This is the finite version of the identity exp(X)Y exp(-X).
func ExpBracketSeries(x, y *Matrix, terms int) (*Matrix, error) {
	if terms < 1 {
		return nil, ErrRange
	}
	result := y.Clone()
	cur := y.Clone()
	fact := 1.0
	for k := 1; k < terms; k++ {
		var err error
		cur, err = Bracket(x, cur)
		if err != nil {
			return nil, err
		}
		fact *= float64(k)
		result, _ = result.Add(cur.Scale(1.0 / fact))
	}
	return result, nil
}
