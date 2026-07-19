package liealgebra

import "math"

// SL2Generators returns the standard basis (E, F, H) of the Lie algebra sl(2)
// of traceless 2x2 real matrices, satisfying [H,E]=2E, [H,F]=-2F and [E,F]=H.
//
//	E = [[0,1],[0,0]], F = [[0,0],[1,0]], H = [[1,0],[0,-1]].
func SL2Generators() (e, f, h *Matrix) {
	e, _ = NewMatrixFromRows([][]float64{{0, 1}, {0, 0}})
	f, _ = NewMatrixFromRows([][]float64{{0, 0}, {1, 0}})
	h, _ = NewMatrixFromRows([][]float64{{1, 0}, {0, -1}})
	return e, f, h
}

// SL2RaisingE returns the raising generator E of sl(2).
func SL2RaisingE() *Matrix { e, _, _ := SL2Generators(); return e }

// SL2LoweringF returns the lowering generator F of sl(2).
func SL2LoweringF() *Matrix { _, f, _ := SL2Generators(); return f }

// SL2CartanH returns the Cartan generator H of sl(2).
func SL2CartanH() *Matrix { _, _, h := SL2Generators(); return h }

// SO3Generators returns the three real antisymmetric generators (Lx, Ly, Lz) of
// the rotation algebra so(3), satisfying [Lx,Ly]=Lz and cyclic permutations.
// Each Li is (L_i)_{jk} = -ε_{ijk}.
func SO3Generators() (lx, ly, lz *Matrix) {
	lx, _ = NewMatrixFromRows([][]float64{{0, 0, 0}, {0, 0, -1}, {0, 1, 0}})
	ly, _ = NewMatrixFromRows([][]float64{{0, 0, 1}, {0, 0, 0}, {-1, 0, 0}})
	lz, _ = NewMatrixFromRows([][]float64{{0, -1, 0}, {1, 0, 0}, {0, 0, 0}})
	return lx, ly, lz
}

// PauliX returns the Pauli matrix σ_x = [[0,1],[1,0]].
func PauliX() *CMatrix {
	m, _ := NewCMatrixFromRows([][]complex128{{0, 1}, {1, 0}})
	return m
}

// PauliY returns the Pauli matrix σ_y = [[0,-i],[i,0]].
func PauliY() *CMatrix {
	m, _ := NewCMatrixFromRows([][]complex128{{0, -1i}, {1i, 0}})
	return m
}

// PauliZ returns the Pauli matrix σ_z = [[1,0],[0,-1]].
func PauliZ() *CMatrix {
	m, _ := NewCMatrixFromRows([][]complex128{{1, 0}, {0, -1}})
	return m
}

// PauliMatrices returns the three Pauli matrices (σ_x, σ_y, σ_z). They are
// Hermitian, traceless, unitary and satisfy σ_iσ_j = δ_ij I + i ε_ijk σ_k.
func PauliMatrices() (sx, sy, sz *CMatrix) {
	return PauliX(), PauliY(), PauliZ()
}

// SU2Generators returns the anti-Hermitian generators of su(2) given by
// T_i = -i σ_i / 2, satisfying [T_i, T_j] = ε_ijk T_k. These are the generators
// for which the exponential map lands in SU(2).
func SU2Generators() (t1, t2, t3 *CMatrix) {
	sx, sy, sz := PauliMatrices()
	half := complex(0, -0.5)
	return sx.Scale(half), sy.Scale(half), sz.Scale(half)
}

// SU2SpinMatrices returns the Hermitian spin-1/2 operators (S_x, S_y, S_z) =
// (σ_x, σ_y, σ_z)/2 satisfying [S_i,S_j] = i ε_ijk S_k.
func SU2SpinMatrices() (sx, sy, sz *CMatrix) {
	x, y, z := PauliMatrices()
	return x.Scale(0.5), y.Scale(0.5), z.Scale(0.5)
}

// GellMannMatrix returns the i-th Gell-Mann matrix, i in 1..8. These are the
// eight Hermitian traceless generators of su(3). It panics if i is out of range.
func GellMannMatrix(i int) *CMatrix {
	var rows [][]complex128
	switch i {
	case 1:
		rows = [][]complex128{{0, 1, 0}, {1, 0, 0}, {0, 0, 0}}
	case 2:
		rows = [][]complex128{{0, -1i, 0}, {1i, 0, 0}, {0, 0, 0}}
	case 3:
		rows = [][]complex128{{1, 0, 0}, {0, -1, 0}, {0, 0, 0}}
	case 4:
		rows = [][]complex128{{0, 0, 1}, {0, 0, 0}, {1, 0, 0}}
	case 5:
		rows = [][]complex128{{0, 0, -1i}, {0, 0, 0}, {1i, 0, 0}}
	case 6:
		rows = [][]complex128{{0, 0, 0}, {0, 0, 1}, {0, 1, 0}}
	case 7:
		rows = [][]complex128{{0, 0, 0}, {0, 0, -1i}, {0, 1i, 0}}
	case 8:
		s := complex(1/math.Sqrt(3), 0)
		rows = [][]complex128{{s, 0, 0}, {0, s, 0}, {0, 0, -2 * s}}
	default:
		panic("liealgebra: GellMannMatrix index must be 1..8")
	}
	m, _ := NewCMatrixFromRows(rows)
	return m
}

// GellMannMatrices returns all eight Gell-Mann matrices λ_1..λ_8.
func GellMannMatrices() []*CMatrix {
	out := make([]*CMatrix, 8)
	for i := 1; i <= 8; i++ {
		out[i-1] = GellMannMatrix(i)
	}
	return out
}

// SU3Generators returns the anti-Hermitian su(3) generators T_a = -i λ_a / 2
// for a in 1..8.
func SU3Generators() []*CMatrix {
	gm := GellMannMatrices()
	out := make([]*CMatrix, 8)
	for a := range gm {
		out[a] = gm[a].Scale(complex(0, -0.5))
	}
	return out
}

// SpinMatrices returns the (2j+1)-dimensional Hermitian angular-momentum
// operators (Jx, Jy, Jz) for spin/quantum number j (a nonnegative half-integer
// such as 0.5, 1, 1.5, ...). They satisfy [Ji,Jj] = i ε_ijk Jk. It returns
// [ErrRange] if 2j is not a nonnegative integer.
func SpinMatrices(j float64) (jx, jy, jz *CMatrix, err error) {
	twoj := int(math.Round(2 * j))
	if twoj < 0 || math.Abs(float64(twoj)-2*j) > 1e-9 {
		return nil, nil, nil, ErrRange
	}
	dim := twoj + 1
	jz = NewCMatrix(dim, dim)
	jp := NewCMatrix(dim, dim) // raising
	jm := NewCMatrix(dim, dim) // lowering
	// Basis ordered by m = j, j-1, ..., -j at indices 0..dim-1.
	for a := 0; a < dim; a++ {
		m := j - float64(a)
		jz.Data[a*dim+a] = complex(m, 0)
	}
	// J+ |j,m> = sqrt(j(j+1)-m(m+1)) |j,m+1>. Index of m is (j-m).
	for a := 1; a < dim; a++ {
		m := j - float64(a) // lower state
		c := math.Sqrt(j*(j+1) - m*(m+1))
		// raises m -> m+1: row (a-1), col a
		jp.Data[(a-1)*dim+a] = complex(c, 0)
	}
	// J- is the conjugate transpose of J+.
	jm = jp.ConjugateTranspose()
	// Jx = (J+ + J-)/2, Jy = (J+ - J-)/(2i).
	sum, _ := jp.Add(jm)
	jx = sum.Scale(0.5)
	diff, _ := jp.Sub(jm)
	jy = diff.Scale(complex(0, -0.5))
	return jx, jy, jz, nil
}

// SpinJz returns the diagonal Jz operator for spin j.
func SpinJz(j float64) (*CMatrix, error) {
	_, _, jz, err := SpinMatrices(j)
	return jz, err
}

// SpinRaising returns the raising operator J+ for spin j.
func SpinRaising(j float64) (*CMatrix, error) {
	twoj := int(math.Round(2 * j))
	if twoj < 0 || math.Abs(float64(twoj)-2*j) > 1e-9 {
		return nil, ErrRange
	}
	dim := twoj + 1
	jp := NewCMatrix(dim, dim)
	for a := 1; a < dim; a++ {
		m := j - float64(a)
		c := math.Sqrt(j*(j+1) - m*(m+1))
		jp.Data[(a-1)*dim+a] = complex(c, 0)
	}
	return jp, nil
}

// SpinLowering returns the lowering operator J- for spin j.
func SpinLowering(j float64) (*CMatrix, error) {
	jp, err := SpinRaising(j)
	if err != nil {
		return nil, err
	}
	return jp.ConjugateTranspose(), nil
}

// HeisenbergGenerators returns the three matrices (X, Y, Z) generating the
// 3-dimensional Heisenberg Lie algebra of strictly upper-triangular 3x3 real
// matrices, satisfying [X,Y]=Z and [X,Z]=[Y,Z]=0.
func HeisenbergGenerators() (x, y, z *Matrix) {
	x, _ = NewMatrixFromRows([][]float64{{0, 1, 0}, {0, 0, 0}, {0, 0, 0}})
	y, _ = NewMatrixFromRows([][]float64{{0, 0, 0}, {0, 0, 1}, {0, 0, 0}})
	z, _ = NewMatrixFromRows([][]float64{{0, 0, 1}, {0, 0, 0}, {0, 0, 0}})
	return x, y, z
}

// SO3ToVector maps a 3x3 antisymmetric matrix to its axial vector (w) such that
// the matrix acts as w×. It returns [ErrDim] if the matrix is not 3x3.
func SO3ToVector(m *Matrix) ([]float64, error) {
	if m.Rows != 3 || m.Cols != 3 {
		return nil, ErrDim
	}
	return []float64{m.At(2, 1), m.At(0, 2), m.At(1, 0)}, nil
}

// SO3FromVector returns the 3x3 antisymmetric matrix [w]_× representing the
// cross product w×(·). It returns [ErrDim] unless w has length 3.
func SO3FromVector(w []float64) (*Matrix, error) {
	if len(w) != 3 {
		return nil, ErrDim
	}
	return NewMatrixFromRows([][]float64{
		{0, -w[2], w[1]},
		{w[2], 0, -w[0]},
		{-w[1], w[0], 0},
	})
}
