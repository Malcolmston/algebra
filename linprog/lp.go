package linprog

import "math"

// Sense selects whether a linear or quadratic objective is minimized or
// maximized.
type Sense int

const (
	// Minimize requests minimization of the objective.
	Minimize Sense = iota
	// Maximize requests maximization of the objective.
	Maximize
)

// String returns "Minimize" or "Maximize".
func (s Sense) String() string {
	if s == Maximize {
		return "Maximize"
	}
	return "Minimize"
}

// Relation is the comparison operator of a linear constraint row.
type Relation int

const (
	// LessEqual denotes a "<=" constraint.
	LessEqual Relation = iota
	// Equal denotes an "=" constraint.
	Equal
	// GreaterEqual denotes a ">=" constraint.
	GreaterEqual
)

// String returns "<=", "=" or ">=".
func (r Relation) String() string {
	switch r {
	case LessEqual:
		return "<="
	case GreaterEqual:
		return ">="
	default:
		return "="
	}
}

// Status classifies the outcome of solving an optimization problem.
type Status int

const (
	// StatusOptimal indicates a finite optimal solution was found.
	StatusOptimal Status = iota
	// StatusInfeasible indicates the feasible region is empty.
	StatusInfeasible
	// StatusUnbounded indicates the objective is unbounded on the feasible
	// region.
	StatusUnbounded
	// StatusIterations indicates the solver hit its iteration limit without
	// converging.
	StatusIterations
)

// String returns a human-readable name for the status.
func (s Status) String() string {
	switch s {
	case StatusOptimal:
		return "Optimal"
	case StatusInfeasible:
		return "Infeasible"
	case StatusUnbounded:
		return "Unbounded"
	case StatusIterations:
		return "IterationLimit"
	default:
		return "Unknown"
	}
}

// LP is a general linear program over nonnegative variables:
//
//	optimize   C · x            (Minimize or Maximize per Sense)
//	subject to A[i] · x  Rel[i]  B[i]   for each constraint i
//	           x >= 0
//
// A is a row-major matrix with one row per constraint and one column per
// variable. C has one entry per variable. Rel and B have one entry per
// constraint.
type LP struct {
	// Sense selects minimization or maximization.
	Sense Sense
	// C is the objective coefficient vector, length NumVars.
	C []float64
	// A is the constraint matrix, NumConstraints rows by NumVars columns.
	A [][]float64
	// Rel holds the relation of each constraint row.
	Rel []Relation
	// B is the right-hand-side vector, length NumConstraints.
	B []float64
}

// NewLP constructs an [LP] from its parts. The slices are copied so the
// returned program is independent of the caller's data. It panics if the
// dimensions are inconsistent.
func NewLP(sense Sense, c []float64, a [][]float64, rel []Relation, b []float64) LP {
	if len(a) != len(rel) || len(a) != len(b) {
		panic(ErrDimension)
	}
	for _, row := range a {
		if len(row) != len(c) {
			panic(ErrDimension)
		}
	}
	return LP{
		Sense: sense,
		C:     append([]float64(nil), c...),
		A:     CloneMatrix(a),
		Rel:   append([]Relation(nil), rel...),
		B:     append([]float64(nil), b...),
	}
}

// NumVars returns the number of decision variables.
func (lp LP) NumVars() int { return len(lp.C) }

// NumConstraints returns the number of constraint rows.
func (lp LP) NumConstraints() int { return len(lp.B) }

// Objective returns the raw objective value C · x, ignoring [Sense].
func (lp LP) Objective(x []float64) float64 { return Dot(lp.C, x) }

// Clone returns an independent deep copy of the program.
func (lp LP) Clone() LP {
	return LP{
		Sense: lp.Sense,
		C:     append([]float64(nil), lp.C...),
		A:     CloneMatrix(lp.A),
		Rel:   append([]Relation(nil), lp.Rel...),
		B:     append([]float64(nil), lp.B...),
	}
}

// Feasible reports whether x satisfies every constraint (within tol) and is
// nonnegative (within tol).
func (lp LP) Feasible(x []float64, tol float64) bool {
	if len(x) != lp.NumVars() {
		return false
	}
	for _, v := range x {
		if v < -tol {
			return false
		}
	}
	for i := range lp.A {
		lhs := Dot(lp.A[i], x)
		switch lp.Rel[i] {
		case LessEqual:
			if lhs > lp.B[i]+tol {
				return false
			}
		case GreaterEqual:
			if lhs < lp.B[i]-tol {
				return false
			}
		default:
			if math.Abs(lhs-lp.B[i]) > tol {
				return false
			}
		}
	}
	return true
}

// AddConstraint returns a copy of the program with one additional constraint
// row appended. The inputs are copied. It panics if row's length does not
// equal NumVars.
func (lp LP) AddConstraint(row []float64, rel Relation, rhs float64) LP {
	if len(row) != lp.NumVars() {
		panic(ErrDimension)
	}
	out := lp.Clone()
	out.A = append(out.A, append([]float64(nil), row...))
	out.Rel = append(out.Rel, rel)
	out.B = append(out.B, rhs)
	return out
}

// StandardLP is a linear program in standard equality form:
//
//	minimize   C · x
//	subject to A x = B
//	           x >= 0
//
// It is the form consumed directly by [Simplex].
type StandardLP struct {
	// C is the objective coefficient vector, length NumVars.
	C []float64
	// A is the equality constraint matrix.
	A [][]float64
	// B is the right-hand-side vector.
	B []float64
}

// NumVars returns the number of variables in the standard program.
func (s StandardLP) NumVars() int { return len(s.C) }

// NumConstraints returns the number of equality rows.
func (s StandardLP) NumConstraints() int { return len(s.B) }

// Objective returns C · x.
func (s StandardLP) Objective(x []float64) float64 { return Dot(s.C, x) }

// Feasible reports whether x satisfies A x = B (within tol) and x >= 0.
func (s StandardLP) Feasible(x []float64, tol float64) bool {
	if len(x) != s.NumVars() {
		return false
	}
	for _, v := range x {
		if v < -tol {
			return false
		}
	}
	for i := range s.A {
		if math.Abs(Dot(s.A[i], x)-s.B[i]) > tol {
			return false
		}
	}
	return true
}

// Standard reduces the general program to [StandardLP] equality form. Each
// "<=" constraint gains a nonnegative slack variable, each ">=" constraint a
// nonnegative surplus variable, and each "=" constraint is copied unchanged.
// The objective is negated when [Sense] is [Maximize] so that the returned
// minimization has the same optimizer.
//
// The first NumVars columns of the returned program are the original
// (structural) variables in the same order, so a solution to the standard
// program restricts to a solution of the original one by truncation.
func (lp LP) Standard() StandardLP {
	n := lp.NumVars()
	m := lp.NumConstraints()
	extra := 0
	for _, r := range lp.Rel {
		if r != Equal {
			extra++
		}
	}
	N := n + extra
	c := make([]float64, N)
	for j := 0; j < n; j++ {
		if lp.Sense == Maximize {
			c[j] = -lp.C[j]
		} else {
			c[j] = lp.C[j]
		}
	}
	a := make([][]float64, m)
	b := append([]float64(nil), lp.B...)
	slack := n
	for i := 0; i < m; i++ {
		row := make([]float64, N)
		copy(row, lp.A[i])
		switch lp.Rel[i] {
		case LessEqual:
			row[slack] = 1
			slack++
		case GreaterEqual:
			row[slack] = -1
			slack++
		}
		a[i] = row
	}
	return StandardLP{C: c, A: a, B: b}
}
