package auctions

import "math"

// lpEps is the pivoting tolerance used by the internal simplex solver.
const lpEps = 1e-9

// lpRel encodes the relation of a linear constraint.
type lpRel int

const (
	relLE lpRel = -1 // less-than-or-equal
	relEQ lpRel = 0  // equal
	relGE lpRel = 1  // greater-than-or-equal
)

// lpConstraint is a single linear constraint coef·x (rel) rhs.
type lpConstraint struct {
	coef []float64
	rel  lpRel
	rhs  float64
}

// lpResult is the outcome of solving a linear program.
type lpResult struct {
	x        []float64
	obj      float64
	feasible bool
	bounded  bool
}

// lpMinimize minimizes c·x subject to the given constraints and x >= 0 using a
// two-phase primal simplex with Bland's anti-cycling rule. Free or ranged
// variables must be modelled by the caller (e.g. by splitting x = p - q).
func lpMinimize(nVars int, c []float64, cons []lpConstraint) lpResult {
	m := len(cons)
	rel := make([]lpRel, m)
	rhs := make([]float64, m)
	rows := make([][]float64, m)
	for i, con := range cons {
		row := make([]float64, nVars)
		copy(row, con.coef)
		r := con.rel
		b := con.rhs
		if b < 0 {
			for j := range row {
				row[j] = -row[j]
			}
			b = -b
			switch r {
			case relLE:
				r = relGE
			case relGE:
				r = relLE
			}
		}
		rows[i] = row
		rel[i] = r
		rhs[i] = b
	}

	slackCol := make([]int, m)
	surplus := make([]bool, m)
	artCol := make([]int, m)
	for i := range slackCol {
		slackCol[i] = -1
		artCol[i] = -1
	}
	nCols := nVars
	for i := 0; i < m; i++ {
		switch rel[i] {
		case relLE:
			slackCol[i] = nCols
			nCols++
		case relGE:
			slackCol[i] = nCols
			surplus[i] = true
			nCols++
		}
	}
	artStart := nCols
	for i := 0; i < m; i++ {
		if rel[i] == relGE || rel[i] == relEQ {
			artCol[i] = nCols
			nCols++
		}
	}

	A := make([][]float64, m)
	basis := make([]int, m)
	bb := make([]float64, m)
	copy(bb, rhs)
	for i := 0; i < m; i++ {
		A[i] = make([]float64, nCols)
		copy(A[i], rows[i])
		if slackCol[i] >= 0 {
			if surplus[i] {
				A[i][slackCol[i]] = -1
			} else {
				A[i][slackCol[i]] = 1
			}
		}
		if artCol[i] >= 0 {
			A[i][artCol[i]] = 1
			basis[i] = artCol[i]
		} else {
			basis[i] = slackCol[i]
		}
	}

	isArt := make([]bool, nCols)
	for i := 0; i < m; i++ {
		if artCol[i] >= 0 {
			isArt[artCol[i]] = true
		}
	}

	if artStart < nCols {
		phase1 := make([]float64, nCols)
		for j := artStart; j < nCols; j++ {
			phase1[j] = 1
		}
		lpSimplex(A, bb, basis, phase1, isArt, false)
		obj := 0.0
		for i := 0; i < m; i++ {
			obj += phase1[basis[i]] * bb[i]
		}
		if obj > 1e-7 {
			return lpResult{feasible: false, bounded: true}
		}
		for i := 0; i < m; i++ {
			if isArt[basis[i]] {
				found := -1
				for j := 0; j < nCols; j++ {
					if isArt[j] {
						continue
					}
					if math.Abs(A[i][j]) > lpEps {
						found = j
						break
					}
				}
				if found >= 0 {
					lpPivot(A, bb, basis, i, found)
				}
			}
		}
	}

	cost := make([]float64, nCols)
	copy(cost, c)
	bounded := lpSimplex(A, bb, basis, cost, isArt, true)
	if !bounded {
		return lpResult{feasible: true, bounded: false}
	}
	x := make([]float64, nVars)
	for i := 0; i < m; i++ {
		if basis[i] < nVars {
			x[basis[i]] = bb[i]
		}
	}
	obj := 0.0
	for j := 0; j < nVars; j++ {
		obj += c[j] * x[j]
	}
	return lpResult{x: x, obj: obj, feasible: true, bounded: true}
}

// lpSimplex runs primal simplex minimization in canonical form with Bland's
// rule. forbid marks columns barred from entering the basis (used to keep
// artificial variables out during phase two). It returns false iff the program
// is unbounded.
func lpSimplex(A [][]float64, b []float64, basis []int, cost []float64, forbid []bool, useForbid bool) bool {
	m := len(A)
	if m == 0 {
		return true
	}
	nCols := len(A[0])
	inBasis := make([]bool, nCols)
	for _, bcol := range basis {
		inBasis[bcol] = true
	}
	for iter := 0; iter < 200000; iter++ {
		entering := -1
		for j := 0; j < nCols; j++ {
			if inBasis[j] {
				continue
			}
			if useForbid && forbid[j] {
				continue
			}
			r := cost[j]
			for i := 0; i < m; i++ {
				r -= cost[basis[i]] * A[i][j]
			}
			if r < -1e-9 {
				entering = j
				break
			}
		}
		if entering == -1 {
			return true
		}
		leaving := -1
		best := math.Inf(1)
		for i := 0; i < m; i++ {
			if A[i][entering] > lpEps {
				ratio := b[i] / A[i][entering]
				switch {
				case leaving == -1 || ratio < best-1e-12:
					best = ratio
					leaving = i
				case math.Abs(ratio-best) <= 1e-12:
					if basis[i] < basis[leaving] {
						leaving = i
					}
				}
			}
		}
		if leaving == -1 {
			return false
		}
		inBasis[basis[leaving]] = false
		inBasis[entering] = true
		lpPivot(A, b, basis, leaving, entering)
	}
	return true
}

// lpPivot performs a Gauss-Jordan pivot on the tableau at (row, col), keeping
// the basis columns as an identity and updating the basis assignment.
func lpPivot(A [][]float64, b []float64, basis []int, row, col int) {
	m := len(A)
	nCols := len(A[0])
	p := A[row][col]
	for j := 0; j < nCols; j++ {
		A[row][j] /= p
	}
	b[row] /= p
	for i := 0; i < m; i++ {
		if i == row {
			continue
		}
		f := A[i][col]
		if f == 0 {
			continue
		}
		for j := 0; j < nCols; j++ {
			A[i][j] -= f * A[row][j]
		}
		b[i] -= f * b[row]
	}
	basis[row] = col
}
