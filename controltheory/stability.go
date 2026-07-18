package controltheory

// RouthResult holds the outcome of a Routh-Hurwitz stability analysis of a
// characteristic polynomial.
type RouthResult struct {
	// Table is the completed Routh array. Row 0 corresponds to the highest
	// power of s. Each row is padded with trailing zeros to a common width.
	Table [][]float64
	// FirstColumn is the first column of the Routh array, whose sign changes
	// count the roots in the right half-plane.
	FirstColumn []float64
	// SignChanges is the number of sign changes in the first column, equal to
	// the number of characteristic roots with positive real part.
	SignChanges int
	// Stable reports whether the polynomial is Hurwitz stable, i.e. every root
	// has a strictly negative real part (no sign changes in the first column).
	Stable bool
}

// controltheoryRouthEpsilon replaces an exact zero in the first column so the
// tabulation can continue (the epsilon method for premature zeros).
const controltheoryRouthEpsilon = 1e-12

// RouthHurwitz performs the Routh-Hurwitz stability test on the given
// characteristic polynomial (ascending-power [Poly]). It builds the Routh
// array, counts sign changes in the first column, and reports stability. A
// zero appearing in the first column is replaced by a small positive epsilon so
// the tabulation can proceed. The polynomial must have degree at least 1.
func RouthHurwitz(p Poly) RouthResult {
	t := controltheoryTrim(p)
	n := t.Degree()
	if n < 1 {
		return RouthResult{Table: nil, FirstColumn: nil, SignChanges: 0, Stable: false}
	}
	// Descending coefficients d[0..n]: d[0] is coeff of s^n.
	d := make([]float64, n+1)
	for i := 0; i <= n; i++ {
		d[i] = t[n-i]
	}
	width := n/2 + 1
	rows := n + 1
	table := make([][]float64, rows)
	for i := range table {
		table[i] = make([]float64, width)
	}
	// First two rows from alternating coefficients.
	for j := 0; 2*j <= n; j++ {
		table[0][j] = d[2*j]
	}
	for j := 0; 2*j+1 <= n; j++ {
		table[1][j] = d[2*j+1]
	}
	for i := 2; i < rows; i++ {
		lead := table[i-1][0]
		if lead == 0 {
			lead = controltheoryRouthEpsilon
			table[i-1][0] = lead
		}
		for j := 0; j < width-1; j++ {
			a := table[i-2][0]
			b := table[i-2][j+1]
			c := table[i-1][0]
			e := table[i-1][j+1]
			table[i][j] = (c*b - a*e) / c
		}
	}
	first := make([]float64, rows)
	for i := 0; i < rows; i++ {
		first[i] = table[i][0]
	}
	changes := 0
	prevSign := controltheorySign(first[0])
	for i := 1; i < rows; i++ {
		s := controltheorySign(first[i])
		if s == 0 {
			s = controltheorySign(controltheoryRouthEpsilon)
		}
		if prevSign != 0 && s != 0 && s != prevSign {
			changes++
		}
		if s != 0 {
			prevSign = s
		}
	}
	return RouthResult{
		Table:       table,
		FirstColumn: first,
		SignChanges: changes,
		Stable:      changes == 0,
	}
}

// controltheorySign returns -1, 0, or +1 according to the sign of x.
func controltheorySign(x float64) int {
	switch {
	case x > 0:
		return 1
	case x < 0:
		return -1
	default:
		return 0
	}
}

// IsHurwitzStable reports whether the polynomial is Hurwitz stable, i.e. all of
// its roots lie strictly in the left half-plane, as determined by the
// Routh-Hurwitz test.
func IsHurwitzStable(p Poly) bool {
	return RouthHurwitz(p).Stable
}

// NumRightHalfPlaneRoots returns the number of roots of the polynomial with
// strictly positive real part, equal to the number of sign changes in the
// first column of the Routh array.
func NumRightHalfPlaneRoots(p Poly) int {
	return RouthHurwitz(p).SignChanges
}
