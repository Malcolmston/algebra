package diffalgebra

import "math/big"

// solveRatSystem solves the linear system A x = b exactly over Q. A is given as
// rows of rational coefficients (each row length == number of unknowns). It
// returns a particular solution (free variables set to zero) and true when the
// system is consistent, or nil and false when it is inconsistent.
func solveRatSystem(A [][]*big.Rat, b []*big.Rat) ([]*big.Rat, bool) {
	rows := len(A)
	if rows == 0 {
		return nil, false
	}
	cols := len(A[0])
	// augmented copy
	m := make([][]*big.Rat, rows)
	for i := 0; i < rows; i++ {
		m[i] = make([]*big.Rat, cols+1)
		for j := 0; j < cols; j++ {
			m[i][j] = cloneRat(A[i][j])
		}
		m[i][cols] = cloneRat(b[i])
	}
	pivotCol := make([]int, 0, cols)
	r := 0
	for c := 0; c < cols && r < rows; c++ {
		// find pivot in column c at row >= r
		piv := -1
		for i := r; i < rows; i++ {
			if !ratZero(m[i][c]) {
				piv = i
				break
			}
		}
		if piv == -1 {
			continue
		}
		m[r], m[piv] = m[piv], m[r]
		// normalise pivot row
		inv := ratInv(m[r][c])
		for j := c; j <= cols; j++ {
			m[r][j] = ratMul(m[r][j], inv)
		}
		// eliminate other rows
		for i := 0; i < rows; i++ {
			if i == r || ratZero(m[i][c]) {
				continue
			}
			f := cloneRat(m[i][c])
			for j := c; j <= cols; j++ {
				m[i][j] = ratSub(m[i][j], ratMul(f, m[r][j]))
			}
		}
		pivotCol = append(pivotCol, c)
		r++
	}
	// consistency: any all-zero coefficient row with nonzero rhs
	for i := 0; i < rows; i++ {
		allZero := true
		for j := 0; j < cols; j++ {
			if !ratZero(m[i][j]) {
				allZero = false
				break
			}
		}
		if allZero && !ratZero(m[i][cols]) {
			return nil, false
		}
	}
	x := make([]*big.Rat, cols)
	for j := range x {
		x[j] = ratInt(0)
	}
	for i, c := range pivotCol {
		x[c] = cloneRat(m[i][cols])
	}
	return x, true
}

// ratSqrt returns the exact rational square root of r when r is a non-negative
// perfect square of a rational, reporting success.
func ratSqrt(r *big.Rat) (*big.Rat, bool) {
	if r.Sign() < 0 {
		return nil, false
	}
	if r.Sign() == 0 {
		return ratInt(0), true
	}
	num := r.Num()
	den := r.Denom()
	ns := new(big.Int).Sqrt(num)
	ds := new(big.Int).Sqrt(den)
	if new(big.Int).Mul(ns, ns).Cmp(num) != 0 {
		return nil, false
	}
	if new(big.Int).Mul(ds, ds).Cmp(den) != 0 {
		return nil, false
	}
	return new(big.Rat).SetFrac(ns, ds), true
}

// ratIsInteger reports whether r is a (possibly negative) integer and returns
// that integer value.
func ratIsInteger(r *big.Rat) (int, bool) {
	if r.Denom().Cmp(big.NewInt(1)) != 0 {
		return 0, false
	}
	if !r.Num().IsInt64() {
		return 0, false
	}
	return int(r.Num().Int64()), true
}
