package designs

import "errors"

// LatinSquare is an n-by-n array of symbols drawn from {0,...,n-1} in which
// every symbol occurs exactly once in each row and once in each column.
type LatinSquare [][]int

// NewLatinSquare validates that the supplied array is a Latin square and
// returns it. It reports an error when the array is not square or the Latin
// property fails.
func NewLatinSquare(rows [][]int) (LatinSquare, error) {
	n := len(rows)
	if n == 0 {
		return nil, errors.New("designs: empty square")
	}
	ls := make(LatinSquare, n)
	for i, r := range rows {
		if len(r) != n {
			return nil, errors.New("designs: array is not square")
		}
		ls[i] = append([]int(nil), r...)
	}
	if !ls.IsLatin() {
		return nil, errors.New("designs: array is not a Latin square")
	}
	return ls, nil
}

// Order returns the side length n of the square.
func (l LatinSquare) Order() int { return len(l) }

// At returns the symbol in row i, column j.
func (l LatinSquare) At(i, j int) int { return l[i][j] }

// IsLatin reports whether the array is a valid Latin square: square shape,
// entries in range, and each symbol appearing once per row and once per column.
func (l LatinSquare) IsLatin() bool {
	n := len(l)
	if n == 0 {
		return false
	}
	for i := 0; i < n; i++ {
		if len(l[i]) != n {
			return false
		}
		rowSeen := make([]bool, n)
		colSeen := make([]bool, n)
		for j := 0; j < n; j++ {
			a := l[i][j]
			if a < 0 || a >= n {
				return false
			}
			if rowSeen[a] {
				return false
			}
			rowSeen[a] = true
			b := l[j][i]
			if b < 0 || b >= n {
				return false
			}
			if colSeen[b] {
				return false
			}
			colSeen[b] = true
		}
	}
	return true
}

// Transpose returns the transpose of the square, which is again a Latin square.
func (l LatinSquare) Transpose() LatinSquare {
	n := len(l)
	t := make(LatinSquare, n)
	for i := 0; i < n; i++ {
		t[i] = make([]int, n)
		for j := 0; j < n; j++ {
			t[i][j] = l[j][i]
		}
	}
	return t
}

// IsSymmetric reports whether the square equals its transpose (a commutative
// quasigroup).
func (l LatinSquare) IsSymmetric() bool {
	n := len(l)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if l[i][j] != l[j][i] {
				return false
			}
		}
	}
	return true
}

// IsIdempotent reports whether l[i][i]==i for every i, i.e. the diagonal is the
// identity.
func (l LatinSquare) IsIdempotent() bool {
	for i := range l {
		if l[i][i] != i {
			return false
		}
	}
	return true
}

// IsReduced reports whether the square is in reduced (standard) form: its first
// row and first column are both 0,1,...,n-1 in order.
func (l LatinSquare) IsReduced() bool {
	n := len(l)
	for i := 0; i < n; i++ {
		if l[0][i] != i || l[i][0] != i {
			return false
		}
	}
	return true
}

// Row returns a copy of row i.
func (l LatinSquare) Row(i int) []int { return append([]int(nil), l[i]...) }

// Col returns a copy of column j.
func (l LatinSquare) Col(j int) []int {
	c := make([]int, len(l))
	for i := range l {
		c[i] = l[i][j]
	}
	return c
}

// CyclicLatinSquare returns the Latin square of order n given by
// l[i][j] = (i+j) mod n, the Cayley table of the cyclic group Z_n.
func CyclicLatinSquare(n int) (LatinSquare, error) {
	if n < 1 {
		return nil, errors.New("designs: order must be positive")
	}
	l := make(LatinSquare, n)
	for i := 0; i < n; i++ {
		l[i] = make([]int, n)
		for j := 0; j < n; j++ {
			l[i][j] = (i + j) % n
		}
	}
	return l, nil
}

// LatinSquareFromField returns the Latin square l[x][y] = a*x + y over GF(q) for
// a fixed non-zero field element a, with rows and columns indexed by field
// elements 0,...,q-1. For a=1 this is the additive Cayley table. It reports an
// error when a is zero or q is not a prime power.
func LatinSquareFromField(q, a int) (LatinSquare, error) {
	f, err := NewGaloisField(q)
	if err != nil {
		return nil, err
	}
	if a == 0 {
		return nil, errors.New("designs: multiplier must be non-zero")
	}
	l := make(LatinSquare, q)
	for x := 0; x < q; x++ {
		l[x] = make([]int, q)
		for y := 0; y < q; y++ {
			l[x][y] = f.Add(f.Mul(a, x), y)
		}
	}
	return l, nil
}

// AreOrthogonal reports whether two Latin squares of the same order are
// orthogonal: superimposing them yields each of the n^2 ordered symbol pairs
// exactly once. It reports an error when the orders differ.
func AreOrthogonal(a, b LatinSquare) (bool, error) {
	n := len(a)
	if len(b) != n {
		return false, errors.New("designs: squares have different orders")
	}
	seen := make([]bool, n*n)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			p := a[i][j]*n + b[i][j]
			if seen[p] {
				return false, nil
			}
			seen[p] = true
		}
	}
	return true, nil
}

// GraecoLatinSquare superimposes two orthogonal Latin squares a and b into a
// single array of ordered pairs encoded as a[i][j]*n + b[i][j]. It reports an
// error when the squares are not orthogonal or differ in order.
func GraecoLatinSquare(a, b LatinSquare) ([][]int, error) {
	ok, err := AreOrthogonal(a, b)
	if err != nil {
		return nil, err
	}
	if !ok {
		return nil, errors.New("designs: squares are not orthogonal")
	}
	n := len(a)
	out := make([][]int, n)
	for i := 0; i < n; i++ {
		out[i] = make([]int, n)
		for j := 0; j < n; j++ {
			out[i][j] = a[i][j]*n + b[i][j]
		}
	}
	return out, nil
}

// MOLS returns a complete set of q-1 mutually orthogonal Latin squares of order
// q, constructed from the finite field GF(q) as L_a[x][y] = a*x + y for each
// non-zero field element a. It reports an error when q is not a prime power.
func MOLS(q int) ([]LatinSquare, error) {
	f, err := NewGaloisField(q)
	if err != nil {
		return nil, err
	}
	out := make([]LatinSquare, 0, q-1)
	for a := 1; a < q; a++ {
		l := make(LatinSquare, q)
		for x := 0; x < q; x++ {
			l[x] = make([]int, q)
			for y := 0; y < q; y++ {
				l[x][y] = f.Add(f.Mul(a, x), y)
			}
		}
		out = append(out, l)
	}
	return out, nil
}

// IsMOLS reports whether the given squares are pairwise mutually orthogonal
// Latin squares of a common order.
func IsMOLS(squares []LatinSquare) bool {
	for _, s := range squares {
		if !s.IsLatin() {
			return false
		}
	}
	for i := 0; i < len(squares); i++ {
		for j := i + 1; j < len(squares); j++ {
			ok, err := AreOrthogonal(squares[i], squares[j])
			if err != nil || !ok {
				return false
			}
		}
	}
	return true
}

// MaxMOLS returns the theoretical upper bound n-1 on the number of mutually
// orthogonal Latin squares of order n, for n>=2.
func MaxMOLS(n int) int {
	if n < 2 {
		return 0
	}
	return n - 1
}

// OrthogonalArrayFromMOLS converts a set of m mutually orthogonal Latin squares
// of order n into the rows of the corresponding orthogonal array OA(n^2, m+2,
// n, 2): each of the n^2 cells (i,j) produces the row [i, j, L_0[i][j], ...,
// L_{m-1}[i][j]].
func OrthogonalArrayFromMOLS(squares []LatinSquare) ([][]int, error) {
	if len(squares) == 0 {
		return nil, errors.New("designs: need at least one square")
	}
	n := len(squares[0])
	for _, s := range squares {
		if len(s) != n {
			return nil, errors.New("designs: squares have different orders")
		}
	}
	rows := make([][]int, 0, n*n)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			row := []int{i, j}
			for _, s := range squares {
				row = append(row, s[i][j])
			}
			rows = append(rows, row)
		}
	}
	return rows, nil
}
