package quasirandom

// Fibonacci returns the n-th Fibonacci number with Fibonacci(0)==0,
// Fibonacci(1)==1. It returns zero for negative n.
func Fibonacci(n int) int {
	if n < 0 {
		return 0
	}
	a, b := 0, 1
	for i := 0; i < n; i++ {
		a, b = b, a+b
	}
	return a
}

// RankOneLatticePoint returns the i-th point (zero-based) of the rank-1 lattice
// rule with integer generating vector gen and total nodes, namely the
// fractional part of i*gen_k/total in each coordinate. It returns an error when
// gen is empty or total < 1.
func RankOneLatticePoint(gen []int, total, i int) ([]float64, error) {
	if len(gen) == 0 {
		return nil, ErrDimension
	}
	if total < 1 {
		return nil, ErrNonPositive
	}
	out := make([]float64, len(gen))
	for k, g := range gen {
		r := ((i*g)%total + total) % total
		out[k] = float64(r) / float64(total)
	}
	return out, nil
}

// RankOneLattice returns all total points of the rank-1 lattice rule with the
// given generating vector. It returns an error when gen is empty or total < 1.
func RankOneLattice(gen []int, total int) ([][]float64, error) {
	if len(gen) == 0 {
		return nil, ErrDimension
	}
	if total < 1 {
		return nil, ErrNonPositive
	}
	out := make([][]float64, total)
	for i := 0; i < total; i++ {
		p, err := RankOneLatticePoint(gen, total, i)
		if err != nil {
			return nil, err
		}
		out[i] = p
	}
	return out, nil
}

// KorobovGenerator returns the Korobov generating vector (1, a, a^2, ..., a^{dim-1})
// reduced modulo total. It returns an error when dim < 1 or total < 1.
func KorobovGenerator(dim, a, total int) ([]int, error) {
	if dim < 1 {
		return nil, ErrDimension
	}
	if total < 1 {
		return nil, ErrNonPositive
	}
	gen := make([]int, dim)
	cur := 1 % total
	for k := 0; k < dim; k++ {
		gen[k] = cur
		cur = (cur * a) % total
	}
	return gen, nil
}

// KorobovLattice returns all total points of the Korobov lattice rule of the
// given dimension with parameter a. It returns an error when dim < 1 or
// total < 1.
func KorobovLattice(dim, a, total int) ([][]float64, error) {
	gen, err := KorobovGenerator(dim, a, total)
	if err != nil {
		return nil, err
	}
	return RankOneLattice(gen, total)
}

// FibonacciLattice returns the two-dimensional Fibonacci lattice with
// Fibonacci(m) points, the rank-1 lattice with generating vector
// (1, Fibonacci(m-1)); it is the classic optimal two-dimensional lattice rule.
// It returns an error when m < 2.
func FibonacciLattice(m int) ([][]float64, error) {
	if m < 2 {
		return nil, ErrDimension
	}
	total := Fibonacci(m)
	gen := []int{1, Fibonacci(m - 1)}
	return RankOneLattice(gen, total)
}

// ShiftedRankOneLatticePoint returns the i-th point of the rank-1 lattice rule
// shifted by the vector shift (each coordinate taken modulo one), a cheap
// randomization-free way to move the lattice off the origin. It returns an
// error when the shapes disagree or total < 1.
func ShiftedRankOneLatticePoint(gen []int, shift []float64, total, i int) ([]float64, error) {
	if len(gen) == 0 || len(shift) != len(gen) {
		return nil, ErrDimension
	}
	if total < 1 {
		return nil, ErrNonPositive
	}
	out := make([]float64, len(gen))
	for k, g := range gen {
		r := ((i*g)%total + total) % total
		out[k] = Frac(float64(r)/float64(total) + shift[k])
	}
	return out, nil
}
