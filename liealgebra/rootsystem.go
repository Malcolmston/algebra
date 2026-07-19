package liealgebra

import (
	"fmt"
	"math/big"
	"strings"
)

// normFamily returns the canonical upper-case single-letter family label.
func normFamily(family string) (string, error) {
	f := strings.ToUpper(strings.TrimSpace(family))
	switch f {
	case "A", "B", "C", "D", "E", "F", "G":
		return f, nil
	}
	return "", ErrType
}

// validRank reports whether rank is admissible for the given family.
func validRank(family string, rank int) error {
	switch family {
	case "A":
		if rank >= 1 {
			return nil
		}
	case "B", "C":
		if rank >= 1 {
			return nil
		}
	case "D":
		if rank >= 2 {
			return nil
		}
	case "G":
		if rank == 2 {
			return nil
		}
	case "F":
		if rank == 4 {
			return nil
		}
	case "E":
		if rank == 6 || rank == 7 || rank == 8 {
			return nil
		}
	}
	return ErrRange
}

// CartanMatrix returns the Cartan matrix of the Lie algebra of the given Dynkin
// type. The family is one of "A","B","C","D" (classical) or "E","F","G"
// (exceptional); rank must be admissible. Entry A[i][j] equals
// 2(α_i,α_j)/(α_j,α_j) in the Bourbaki convention.
func CartanMatrix(family string, rank int) (*Matrix, error) {
	f, err := normFamily(family)
	if err != nil {
		return nil, err
	}
	if err := validRank(f, rank); err != nil {
		return nil, err
	}
	switch f {
	case "A":
		return cartanA(rank), nil
	case "B":
		return cartanB(rank), nil
	case "C":
		return cartanC(rank), nil
	case "D":
		return cartanD(rank), nil
	case "G":
		return NewMatrixFromRows([][]float64{{2, -1}, {-3, 2}})
	case "F":
		return NewMatrixFromRows([][]float64{
			{2, -1, 0, 0},
			{-1, 2, -2, 0},
			{0, -1, 2, -1},
			{0, 0, -1, 2},
		})
	case "E":
		return cartanE(rank), nil
	}
	return nil, ErrType
}

// CartanMatrixA returns the Cartan matrix of type A_n (the Lie algebra sl(n+1)).
func CartanMatrixA(n int) (*Matrix, error) { return CartanMatrix("A", n) }

// CartanMatrixB returns the Cartan matrix of type B_n (the Lie algebra so(2n+1)).
func CartanMatrixB(n int) (*Matrix, error) { return CartanMatrix("B", n) }

// CartanMatrixC returns the Cartan matrix of type C_n (the Lie algebra sp(2n)).
func CartanMatrixC(n int) (*Matrix, error) { return CartanMatrix("C", n) }

// CartanMatrixD returns the Cartan matrix of type D_n (the Lie algebra so(2n)).
func CartanMatrixD(n int) (*Matrix, error) { return CartanMatrix("D", n) }

func cartanA(n int) *Matrix {
	m := NewMatrix(n, n)
	for i := 0; i < n; i++ {
		m.Data[i*n+i] = 2
		if i+1 < n {
			m.Data[i*n+i+1] = -1
			m.Data[(i+1)*n+i] = -1
		}
	}
	return m
}

func cartanB(n int) *Matrix {
	m := cartanA(n)
	if n >= 2 {
		// α_n is the short root: A[n-2][n-1] = 2(α_{n-1},α_n)/(α_n,α_n) = -2.
		m.Data[(n-2)*n+(n-1)] = -2
	}
	return m
}

func cartanC(n int) *Matrix {
	m := cartanA(n)
	if n >= 2 {
		// α_n is the long root: A[n-1][n-2] = 2(α_n,α_{n-1})/(α_{n-1},α_{n-1}) = -2.
		m.Data[(n-1)*n+(n-2)] = -2
	}
	return m
}

func cartanD(n int) *Matrix {
	m := NewMatrix(n, n)
	for i := 0; i < n; i++ {
		m.Data[i*n+i] = 2
	}
	// Chain among first n-1 nodes: 0-1-...-(n-2).
	for i := 0; i+1 <= n-2; i++ {
		m.Data[i*n+(i+1)] = -1
		m.Data[(i+1)*n+i] = -1
	}
	// Fork: node n-1 attaches to node n-3.
	if n >= 3 {
		m.Data[(n-3)*n+(n-1)] = -1
		m.Data[(n-1)*n+(n-3)] = -1
	} else if n == 2 {
		// D2 = A1 x A1: no edges.
	}
	return m
}

func cartanE(n int) *Matrix {
	m := NewMatrix(n, n)
	for i := 0; i < n; i++ {
		m.Data[i*n+i] = 2
	}
	set := func(i, j int) {
		m.Data[i*n+j] = -1
		m.Data[j*n+i] = -1
	}
	// Bourbaki E_n: chain 1-3-4-5-...-n with node 2 attached to node 4.
	// 0-indexed nodes; node 1 is the branch attached to node 3.
	set(0, 2)
	set(2, 3)
	set(1, 3)
	for i := 3; i+1 < n; i++ {
		set(i, i+1)
	}
	return m
}

// CartanMatrixDeterminant returns the determinant of the Cartan matrix, an
// important invariant equal to the order of the center of the simply connected
// group (the index of the root lattice in the weight lattice).
func CartanMatrixDeterminant(family string, rank int) (float64, error) {
	c, err := CartanMatrix(family, rank)
	if err != nil {
		return 0, err
	}
	return Det(c)
}

// InverseCartanMatrix returns the inverse of the Cartan matrix. Its rows give
// the fundamental weights expressed in the basis of simple roots.
func InverseCartanMatrix(family string, rank int) (*Matrix, error) {
	c, err := CartanMatrix(family, rank)
	if err != nil {
		return nil, err
	}
	return Inverse(c)
}

// IsSimplyLaced reports whether the Dynkin type is simply laced, i.e. all
// bonds are single (types A, D, E).
func IsSimplyLaced(family string, rank int) (bool, error) {
	c, err := CartanMatrix(family, rank)
	if err != nil {
		return false, err
	}
	n := c.Rows
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if i != j && c.Data[i*n+j] != 0 && c.Data[i*n+j] != -1 {
				return false, nil
			}
		}
	}
	return true, nil
}

// DynkinDiagramAdjacency returns the 0/1 adjacency matrix of the Dynkin diagram:
// nodes i and j are adjacent when the Cartan entry A[i][j] is nonzero.
func DynkinDiagramAdjacency(family string, rank int) (*Matrix, error) {
	c, err := CartanMatrix(family, rank)
	if err != nil {
		return nil, err
	}
	n := c.Rows
	adj := NewMatrix(n, n)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if i != j && c.Data[i*n+j] != 0 {
				adj.Data[i*n+j] = 1
			}
		}
	}
	return adj, nil
}

// DynkinBondMatrix returns the matrix B[i][j] = A[i][j]*A[j][i], whose value
// (0,1,2,3) counts the number of bonds between nodes i and j in the Dynkin
// diagram.
func DynkinBondMatrix(family string, rank int) (*Matrix, error) {
	c, err := CartanMatrix(family, rank)
	if err != nil {
		return nil, err
	}
	n := c.Rows
	out := NewMatrix(n, n)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if i != j {
				out.Data[i*n+j] = c.Data[i*n+j] * c.Data[j*n+i]
			}
		}
	}
	return out, nil
}

// intKey builds a stable map key for an integer coefficient vector.
func intKey(v []int) string {
	var b strings.Builder
	for i, x := range v {
		if i > 0 {
			b.WriteByte(',')
		}
		fmt.Fprintf(&b, "%d", x)
	}
	return b.String()
}

func cloneInts(v []int) []int {
	out := make([]int, len(v))
	copy(out, v)
	return out
}

// PositiveRootLabels returns every positive root of the given type expressed as
// an integer coefficient vector in the basis of simple roots. The algorithm is
// the standard height-by-height root-string construction driven by the Cartan
// matrix and works uniformly for classical and exceptional types.
func PositiveRootLabels(family string, rank int) ([][]int, error) {
	c, err := CartanMatrix(family, rank)
	if err != nil {
		return nil, err
	}
	return positiveRootLabels(c), nil
}

func positiveRootLabels(c *Matrix) [][]int {
	n := c.Rows
	set := map[string]bool{}
	var roots [][]int
	add := func(v []int) bool {
		k := intKey(v)
		if set[k] {
			return false
		}
		set[k] = true
		roots = append(roots, cloneInts(v))
		return true
	}
	for i := 0; i < n; i++ {
		v := make([]int, n)
		v[i] = 1
		add(v)
	}
	changed := true
	for changed {
		changed = false
		snapshot := make([][]int, len(roots))
		copy(snapshot, roots)
		for _, beta := range snapshot {
			for i := 0; i < n; i++ {
				// p = largest k>=0 with beta - k*alpha_i still a root.
				p := 0
				test := cloneInts(beta)
				for {
					test[i]--
					if set[intKey(test)] {
						p++
					} else {
						break
					}
				}
				// Cartan integer <beta, alpha_i^vee> = sum_j beta_j A[j][i].
				cv := 0
				for j := 0; j < n; j++ {
					cv += beta[j] * int(c.Data[j*n+i])
				}
				q := p - cv
				if q > 0 {
					nb := cloneInts(beta)
					nb[i]++
					if add(nb) {
						changed = true
					}
				}
			}
		}
	}
	return roots
}

// NumPositiveRoots returns the number of positive roots of the given type.
func NumPositiveRoots(family string, rank int) (int, error) {
	labels, err := PositiveRootLabels(family, rank)
	if err != nil {
		return 0, err
	}
	return len(labels), nil
}

// NumRoots returns the total number of roots (positive and negative).
func NumRoots(family string, rank int) (int, error) {
	p, err := NumPositiveRoots(family, rank)
	if err != nil {
		return 0, err
	}
	return 2 * p, nil
}

// LieAlgebraDimension returns the dimension of the semisimple Lie algebra of the
// given type, equal to rank + number of roots.
func LieAlgebraDimension(family string, rank int) (int, error) {
	nr, err := NumRoots(family, rank)
	if err != nil {
		return 0, err
	}
	return rank + nr, nil
}

// HighestRootLabel returns the coefficient vector (in the simple-root basis) of
// the highest root, the unique positive root of maximal height.
func HighestRootLabel(family string, rank int) ([]int, error) {
	labels, err := PositiveRootLabels(family, rank)
	if err != nil {
		return nil, err
	}
	best := labels[0]
	bestH := heightOf(best)
	for _, l := range labels[1:] {
		if h := heightOf(l); h > bestH {
			bestH = h
			best = l
		}
	}
	return cloneInts(best), nil
}

func heightOf(v []int) int {
	s := 0
	for _, x := range v {
		s += x
	}
	return s
}

// RootBilinearForm returns the symmetric matrix B with B[i][j] = (α_i,α_j), the
// inner products of the simple roots normalised so that the longest roots have
// squared length 2. It is derived from the Cartan matrix and underlies the
// Weyl dimension and Casimir formulas.
func RootBilinearForm(family string, rank int) (*Matrix, error) {
	c, err := CartanMatrix(family, rank)
	if err != nil {
		return nil, err
	}
	return rootBilinear(c), nil
}

func rootBilinear(c *Matrix) *Matrix {
	n := c.Rows
	s := make([]float64, n)
	visited := make([]bool, n)
	// Process each connected component of the Dynkin diagram.
	for start := 0; start < n; start++ {
		if visited[start] {
			continue
		}
		s[start] = 1
		visited[start] = true
		queue := []int{start}
		for len(queue) > 0 {
			i := queue[0]
			queue = queue[1:]
			for j := 0; j < n; j++ {
				if j == i || c.Data[i*n+j] == 0 {
					continue
				}
				if !visited[j] {
					// s_i/s_j = A_ij/A_ji => s_j = s_i * A_ji/A_ij.
					s[j] = s[i] * c.Data[j*n+i] / c.Data[i*n+j]
					visited[j] = true
					queue = append(queue, j)
				}
			}
		}
	}
	max := 0.0
	for _, v := range s {
		if v > max {
			max = v
		}
	}
	if max == 0 {
		max = 1
	}
	scale := 2.0 / max
	for i := range s {
		s[i] *= scale
	}
	b := NewMatrix(n, n)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			b.Data[i*n+j] = c.Data[i*n+j] * s[j] / 2
		}
	}
	return b
}

// RootInnerProduct returns the inner product of two vectors expressed in the
// simple-root basis using the bilinear form B (see [RootBilinearForm]).
func RootInnerProduct(b *Matrix, x, y []float64) (float64, error) {
	n := b.Rows
	if len(x) != n || len(y) != n {
		return 0, ErrDim
	}
	s := 0.0
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			s += x[i] * b.Data[i*n+j] * y[j]
		}
	}
	return s, nil
}

// FundamentalWeightsRootBasis returns the fundamental weights expressed in the
// basis of simple roots; row i is ω_i. They are the rows of the inverse Cartan
// matrix.
func FundamentalWeightsRootBasis(family string, rank int) ([][]float64, error) {
	inv, err := InverseCartanMatrix(family, rank)
	if err != nil {
		return nil, err
	}
	n := inv.Rows
	out := make([][]float64, n)
	for i := 0; i < n; i++ {
		out[i] = inv.Row(i)
	}
	return out, nil
}

// WeylVectorRootBasis returns the Weyl vector ρ (the sum of the fundamental
// weights, equivalently half the sum of the positive roots) expressed in the
// simple-root basis.
func WeylVectorRootBasis(family string, rank int) ([]float64, error) {
	inv, err := InverseCartanMatrix(family, rank)
	if err != nil {
		return nil, err
	}
	n := inv.Rows
	rho := make([]float64, n)
	for k := 0; k < n; k++ {
		s := 0.0
		for i := 0; i < n; i++ {
			s += inv.Data[i*n+k]
		}
		rho[k] = s
	}
	return rho, nil
}

// WeylGroupOrder returns the order of the Weyl group of the given Dynkin type as
// an arbitrary-precision integer.
func WeylGroupOrder(family string, rank int) (*big.Int, error) {
	f, err := normFamily(family)
	if err != nil {
		return nil, err
	}
	if err := validRank(f, rank); err != nil {
		return nil, err
	}
	n := int64(rank)
	fact := func(k int64) *big.Int {
		r := big.NewInt(1)
		for i := int64(2); i <= k; i++ {
			r.Mul(r, big.NewInt(i))
		}
		return r
	}
	pow2 := func(k int64) *big.Int {
		return new(big.Int).Lsh(big.NewInt(1), uint(k))
	}
	switch f {
	case "A":
		return fact(n + 1), nil
	case "B", "C":
		return new(big.Int).Mul(pow2(n), fact(n)), nil
	case "D":
		return new(big.Int).Mul(pow2(n-1), fact(n)), nil
	case "G":
		return big.NewInt(12), nil
	case "F":
		return big.NewInt(1152), nil
	case "E":
		switch rank {
		case 6:
			return big.NewInt(51840), nil
		case 7:
			return big.NewInt(2903040), nil
		case 8:
			return big.NewInt(696729600), nil
		}
	}
	return nil, ErrType
}

// CoxeterNumber returns the Coxeter number h of the given Dynkin type, equal to
// (number of roots)/rank and to one plus the height of the highest root.
func CoxeterNumber(family string, rank int) (int, error) {
	f, err := normFamily(family)
	if err != nil {
		return 0, err
	}
	if err := validRank(f, rank); err != nil {
		return 0, err
	}
	switch f {
	case "A":
		return rank + 1, nil
	case "B", "C":
		return 2 * rank, nil
	case "D":
		return 2*rank - 2, nil
	case "G":
		return 6, nil
	case "F":
		return 12, nil
	case "E":
		switch rank {
		case 6:
			return 12, nil
		case 7:
			return 18, nil
		case 8:
			return 30, nil
		}
	}
	return 0, ErrType
}

// DualCoxeterNumber returns the dual Coxeter number h∨ of the given Dynkin type.
func DualCoxeterNumber(family string, rank int) (int, error) {
	f, err := normFamily(family)
	if err != nil {
		return 0, err
	}
	if err := validRank(f, rank); err != nil {
		return 0, err
	}
	switch f {
	case "A":
		return rank + 1, nil
	case "B":
		return 2*rank - 1, nil
	case "C":
		return rank + 1, nil
	case "D":
		return 2*rank - 2, nil
	case "G":
		return 4, nil
	case "F":
		return 9, nil
	case "E":
		switch rank {
		case 6:
			return 12, nil
		case 7:
			return 18, nil
		case 8:
			return 30, nil
		}
	}
	return 0, ErrType
}

// WeylReflection returns the reflection of the vector v across the hyperplane
// orthogonal to root, s_root(v) = v - 2(v,root)/(root,root) root, in Euclidean
// coordinates. It returns [ErrDim] on a length mismatch or a zero root.
func WeylReflection(root, v []float64) ([]float64, error) {
	if len(root) != len(v) {
		return nil, ErrDim
	}
	rr := VecNormSquared(root)
	if rr == 0 {
		return nil, ErrDim
	}
	vr, _ := VecDot(v, root)
	coeff := 2 * vr / rr
	out := make([]float64, len(v))
	for i := range v {
		out[i] = v[i] - coeff*root[i]
	}
	return out, nil
}

// WeylReflectionMatrix returns the orthogonal matrix of the reflection across
// the hyperplane orthogonal to the given root in Euclidean coordinates.
func WeylReflectionMatrix(root []float64) (*Matrix, error) {
	rr := VecNormSquared(root)
	if rr == 0 {
		return nil, ErrDim
	}
	n := len(root)
	m := IdentityMatrix(n)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			m.Data[i*n+j] -= 2 * root[i] * root[j] / rr
		}
	}
	return m, nil
}

// Coroot returns the coroot α∨ = 2α/(α,α) of a Euclidean root vector.
func Coroot(root []float64) ([]float64, error) {
	rr := VecNormSquared(root)
	if rr == 0 {
		return nil, ErrDim
	}
	return VecScale(root, 2/rr), nil
}

// CartanInteger returns the Cartan integer ⟨α,β∨⟩ = 2(α,β)/(β,β) of two
// Euclidean root vectors.
func CartanInteger(alpha, beta []float64) (float64, error) {
	bb := VecNormSquared(beta)
	if bb == 0 {
		return 0, ErrDim
	}
	ab, err := VecDot(alpha, beta)
	if err != nil {
		return 0, err
	}
	return 2 * ab / bb, nil
}
