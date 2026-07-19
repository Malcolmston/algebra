package tensornetwork

import (
	"errors"
	"fmt"
	"math"
	"math/rand"
)

// CPDecomposition is a rank-R CANDECOMP/PARAFAC factorization of a tensor. The
// approximation is the sum over r of Weights[r] times the outer product of the
// r-th columns of the factor matrices, one factor per mode. Each factor has
// dimensions (Iₙ, R) with unit-norm columns.
type CPDecomposition struct {
	// Factors holds one (Iₙ, R) matrix per tensor mode.
	Factors []*Matrix
	// Weights holds the R nonnegative weights λ that scale the rank-one terms.
	Weights []float64
	// Rank is the number of rank-one components R.
	Rank int
	// Shape is the shape of the reconstructed tensor.
	Shape []int
}

// CPOptions configures [CPALS].
type CPOptions struct {
	// MaxIter is the maximum number of ALS sweeps.
	MaxIter int
	// Tol is the convergence tolerance on the relative change in fit.
	Tol float64
	// Seed seeds the deterministic random initialization of the factors.
	Seed int64
}

// DefaultCPOptions returns reasonable defaults for [CPALS]: 200 iterations, a
// tolerance of 1e-10 and seed 1.
func DefaultCPOptions() CPOptions {
	return CPOptions{MaxIter: 200, Tol: 1e-10, Seed: 1}
}

// CPALS computes a rank-R CP decomposition of t by alternating least squares.
// The factor matrices are initialized from a deterministic random source seeded
// by opts.Seed and refined until the relative fit changes by less than opts.Tol
// or opts.MaxIter sweeps have run. It returns an error if rank < 1 or t has rank
// below 2.
func CPALS(t *Tensor, rank int, opts CPOptions) (*CPDecomposition, error) {
	if rank < 1 {
		return nil, errors.New("tensornetwork: CP rank must be >= 1")
	}
	nd := len(t.shape)
	if nd < 2 {
		return nil, errors.New("tensornetwork: CPALS requires rank >= 2 tensor")
	}
	if opts.MaxIter <= 0 {
		opts.MaxIter = 200
	}
	rng := rand.New(rand.NewSource(opts.Seed))
	factors := make([]*Matrix, nd)
	for n := 0; n < nd; n++ {
		f := NewMatrix(t.shape[n], rank)
		for i := range f.data {
			f.data[i] = rng.NormFloat64()
		}
		factors[n] = f
	}
	weights := make([]float64, rank)
	for r := range weights {
		weights[r] = 1
	}
	normT := t.Norm()
	prevFit := -1.0
	unfoldings := make([]*Matrix, nd)
	for n := 0; n < nd; n++ {
		unfoldings[n], _ = t.Unfold(n)
	}
	for iter := 0; iter < opts.MaxIter; iter++ {
		for n := 0; n < nd; n++ {
			// V = Hadamard over k!=n of (A_k^T A_k).
			v := onesMatrix(rank, rank)
			var krFactors []*Matrix
			for k := 0; k < nd; k++ {
				if k == n {
					continue
				}
				g := factors[k].GramColumns()
				v, _ = v.Hadamard(g)
				krFactors = append(krFactors, factors[k])
			}
			kr, err := KhatriRao(krFactors...)
			if err != nil {
				return nil, err
			}
			// A_n = X_(n) * KR * pinv(V).
			mkr, err := unfoldings[n].Mul(kr)
			if err != nil {
				return nil, err
			}
			an, err := mkr.Mul(v.Pinv(0))
			if err != nil {
				return nil, err
			}
			// Normalize columns; store weights.
			for r := 0; r < rank; r++ {
				var nrm float64
				for i := 0; i < an.rows; i++ {
					nrm += an.data[i*rank+r] * an.data[i*rank+r]
				}
				nrm = math.Sqrt(nrm)
				weights[r] = nrm
				if nrm > 1e-300 {
					for i := 0; i < an.rows; i++ {
						an.data[i*rank+r] /= nrm
					}
				}
			}
			factors[n] = an
		}
		cp := &CPDecomposition{Factors: factors, Weights: weights, Rank: rank, Shape: append([]int(nil), t.shape...)}
		rec := cp.Reconstruct()
		diff, _ := t.Sub(rec)
		fit := 1.0
		if normT > 0 {
			fit = 1 - diff.Norm()/normT
		}
		if prevFit >= 0 && math.Abs(fit-prevFit) < opts.Tol {
			prevFit = fit
			break
		}
		prevFit = fit
	}
	return &CPDecomposition{Factors: factors, Weights: weights, Rank: rank, Shape: append([]int(nil), t.shape...)}, nil
}

// onesMatrix returns an r x c matrix of all ones.
func onesMatrix(r, c int) *Matrix {
	m := NewMatrix(r, c)
	for i := range m.data {
		m.data[i] = 1
	}
	return m
}

// Reconstruct returns the dense tensor represented by the CP decomposition.
func (cp *CPDecomposition) Reconstruct() *Tensor {
	nd := len(cp.Factors)
	out := New(cp.Shape...)
	idx := make([]int, nd)
	for flat := 0; flat < len(out.data); flat++ {
		rem := flat
		for a := 0; a < nd; a++ {
			idx[a] = (rem / out.stride[a]) % cp.Shape[a]
		}
		var s float64
		for r := 0; r < cp.Rank; r++ {
			p := cp.Weights[r]
			for n := 0; n < nd; n++ {
				p *= cp.Factors[n].data[idx[n]*cp.Rank+r]
			}
			s += p
		}
		out.data[flat] = s
	}
	return out
}

// RelError returns the relative Frobenius reconstruction error of the CP
// decomposition against t.
func (cp *CPDecomposition) RelError(t *Tensor) (float64, error) {
	return t.RelError(cp.Reconstruct())
}

// TuckerDecomposition is a Tucker factorization: a core tensor together with one
// orthonormal factor matrix per mode. The approximation is the core multiplied
// along each mode by the corresponding factor.
type TuckerDecomposition struct {
	// Core is the (R₀,…,R_{N-1}) core tensor.
	Core *Tensor
	// Factors holds one (Iₙ, Rₙ) orthonormal factor per mode.
	Factors []*Matrix
}

// HOSVD computes a truncated Tucker decomposition of t by the higher-order
// singular value decomposition. ranks gives the target multilinear rank per
// mode; a value <= 0 or larger than the mode size keeps the full mode. It
// returns an error if len(ranks) != t.Rank().
func HOSVD(t *Tensor, ranks []int) (*TuckerDecomposition, error) {
	nd := len(t.shape)
	if len(ranks) != nd {
		return nil, fmt.Errorf("tensornetwork: ranks length %d != tensor rank %d", len(ranks), nd)
	}
	factors := make([]*Matrix, nd)
	for n := 0; n < nd; n++ {
		unf, err := t.Unfold(n)
		if err != nil {
			return nil, err
		}
		u, _, _, err := unf.SVD()
		if err != nil {
			return nil, err
		}
		r := ranks[n]
		if r <= 0 || r > u.cols {
			r = u.cols
		}
		factors[n] = u.SubMatrix(0, u.rows, 0, r)
	}
	core, err := tuckerCore(t, factors)
	if err != nil {
		return nil, err
	}
	return &TuckerDecomposition{Core: core, Factors: factors}, nil
}

// tuckerCore computes the core tensor G = T ×₀ A₀ᵀ ×₁ A₁ᵀ … for orthonormal
// factors.
func tuckerCore(t *Tensor, factors []*Matrix) (*Tensor, error) {
	cur := t
	var err error
	for n := 0; n < len(factors); n++ {
		cur, err = ModeProduct(cur, factors[n].Transpose(), n)
		if err != nil {
			return nil, err
		}
	}
	return cur, nil
}

// HOOI refines a Tucker decomposition of t by higher-order orthogonal iteration,
// starting from the HOSVD and alternately re-estimating each factor. ranks is
// the target multilinear rank; maxIter bounds the number of sweeps. It returns
// an error if the ranks are invalid.
func HOOI(t *Tensor, ranks []int, maxIter int) (*TuckerDecomposition, error) {
	tuck, err := HOSVD(t, ranks)
	if err != nil {
		return nil, err
	}
	nd := len(t.shape)
	if maxIter <= 0 {
		maxIter = 50
	}
	factors := tuck.Factors
	for iter := 0; iter < maxIter; iter++ {
		for n := 0; n < nd; n++ {
			cur := t
			for k := 0; k < nd; k++ {
				if k == n {
					continue
				}
				cur, err = ModeProduct(cur, factors[k].Transpose(), k)
				if err != nil {
					return nil, err
				}
			}
			unf, err := cur.Unfold(n)
			if err != nil {
				return nil, err
			}
			u, _, _, err := unf.SVD()
			if err != nil {
				return nil, err
			}
			r := ranks[n]
			if r <= 0 || r > u.cols {
				r = u.cols
			}
			factors[n] = u.SubMatrix(0, u.rows, 0, r)
		}
	}
	core, err := tuckerCore(t, factors)
	if err != nil {
		return nil, err
	}
	return &TuckerDecomposition{Core: core, Factors: factors}, nil
}

// Reconstruct returns the dense tensor represented by the Tucker decomposition.
func (tu *TuckerDecomposition) Reconstruct() *Tensor {
	cur := tu.Core
	for n := 0; n < len(tu.Factors); n++ {
		cur, _ = ModeProduct(cur, tu.Factors[n], n)
	}
	return cur
}

// RelError returns the relative Frobenius reconstruction error of the Tucker
// decomposition against t.
func (tu *TuckerDecomposition) RelError(t *Tensor) (float64, error) {
	return t.RelError(tu.Reconstruct())
}

// MultilinearRank returns the numerical multilinear rank of t: the rank of each
// mode-n unfolding, using the given relative singular-value tolerance.
func MultilinearRank(t *Tensor, tol float64) []int {
	nd := len(t.shape)
	out := make([]int, nd)
	for n := 0; n < nd; n++ {
		unf, _ := t.Unfold(n)
		out[n] = unf.Rank(tol)
	}
	return out
}

// TensorTrain is a tensor-train (matrix-product-state) factorization. Core k has
// shape (Rₖ, nₖ, Rₖ₊₁) with R₀ = R_d = 1, and the represented tensor is the
// sequential contraction of the cores.
type TensorTrain struct {
	// Cores holds the d three-index cores of the train.
	Cores []*Tensor
	// Ranks holds the d+1 TT ranks, beginning and ending with 1.
	Ranks []int
	// Shape holds the shape of the represented tensor.
	Shape []int
}

// TTSVD computes a tensor-train decomposition of t using successive truncated
// SVDs. Singular values are discarded so that the total truncation error stays
// within eps (relative to ‖t‖); eps <= 0 gives an essentially exact
// decomposition. It returns an error if t has rank below 1.
func TTSVD(t *Tensor, eps float64) (*TensorTrain, error) {
	return ttsvd(t, eps, nil)
}

// TTSVDRank computes a tensor-train decomposition of t whose internal ranks are
// capped by maxRanks (length t.Rank()-1). A cap <= 0 leaves that rank
// uncapped. It returns an error if maxRanks has the wrong length.
func TTSVDRank(t *Tensor, maxRanks []int) (*TensorTrain, error) {
	if len(maxRanks) != len(t.shape)-1 {
		return nil, fmt.Errorf("tensornetwork: maxRanks length %d != %d", len(maxRanks), len(t.shape)-1)
	}
	return ttsvd(t, 0, maxRanks)
}

// ttsvd is the shared implementation of the tensor-train SVD.
func ttsvd(t *Tensor, eps float64, maxRanks []int) (*TensorTrain, error) {
	nd := len(t.shape)
	if nd < 1 {
		return nil, errors.New("tensornetwork: TTSVD requires rank >= 1 tensor")
	}
	normT := t.Norm()
	delta := 0.0
	if eps > 0 && nd > 1 {
		delta = eps * normT / math.Sqrt(float64(nd-1))
	}
	cores := make([]*Tensor, nd)
	ranks := make([]int, nd+1)
	ranks[0] = 1
	ranks[nd] = 1
	c := t.Clone()
	rprev := 1
	total := len(t.data)
	for k := 0; k < nd-1; k++ {
		nk := t.shape[k]
		rows := rprev * nk
		cols := total / rows
		mat, err := NewMatrixData(rows, cols, c.data)
		if err != nil {
			return nil, err
		}
		u, s, v, err := mat.SVD()
		if err != nil {
			return nil, err
		}
		// Choose truncation rank.
		r := chooseRank(s, delta, rank2Cap(maxRanks, k))
		if r < 1 {
			r = 1
		}
		// Core k = reshape U[:, :r] into (rprev, nk, r).
		coreData := make([]float64, rows*r)
		for i := 0; i < rows; i++ {
			for j := 0; j < r; j++ {
				coreData[i*r+j] = u.data[i*u.cols+j]
			}
		}
		core, err := NewWithData(coreData, rprev, nk, r)
		if err != nil {
			return nil, err
		}
		cores[k] = core
		ranks[k+1] = r
		// Remaining = diag(s[:r]) * V[:, :r]^T  (r x cols).
		rem := make([]float64, r*cols)
		for j := 0; j < r; j++ {
			for col := 0; col < cols; col++ {
				rem[j*cols+col] = s[j] * v.data[col*v.cols+j]
			}
		}
		c = &Tensor{shape: []int{r, cols}, stride: rowMajorStrides([]int{r, cols}), data: rem}
		rprev = r
		total = r * cols
	}
	// Last core: (rprev, n_{d-1}, 1).
	lastData := make([]float64, len(c.data))
	copy(lastData, c.data)
	lastCore, err := NewWithData(lastData, rprev, t.shape[nd-1], 1)
	if err != nil {
		return nil, err
	}
	cores[nd-1] = lastCore
	return &TensorTrain{Cores: cores, Ranks: ranks, Shape: append([]int(nil), t.shape...)}, nil
}

// rank2Cap returns the cap for internal rank index k, or 0 if uncapped.
func rank2Cap(maxRanks []int, k int) int {
	if maxRanks == nil || k >= len(maxRanks) {
		return 0
	}
	return maxRanks[k]
}

// chooseRank picks a truncation rank from descending singular values so that the
// discarded tail has 2-norm at most delta, further capped by cap (if positive).
func chooseRank(s []float64, delta float64, cap int) int {
	n := len(s)
	// Largest rank satisfying the error bound.
	r := n
	if delta > 0 {
		var tail float64
		r = 0
		// Accumulate from the smallest singular values upward.
		for i := n - 1; i >= 0; i-- {
			if math.Sqrt(tail+s[i]*s[i]) <= delta {
				tail += s[i] * s[i]
			} else {
				r = i + 1
				break
			}
		}
	}
	if r < 1 {
		r = 1
	}
	// Drop exact zeros.
	for r > 1 && s[r-1] <= 1e-14*s[0] {
		r--
	}
	if cap > 0 && r > cap {
		r = cap
	}
	return r
}

// Reconstruct returns the dense tensor represented by the tensor train.
func (tt *TensorTrain) Reconstruct() *Tensor {
	nd := len(tt.Cores)
	// Merge cores left to right: acc has shape (n0,...,n_{k}, r_{k+1}).
	acc := tt.Cores[0].Clone() // (1, n0, r1)
	accShape := []int{tt.Shape[0], tt.Ranks[1]}
	accData := make([]float64, len(acc.data))
	copy(accData, acc.data)
	acc = &Tensor{shape: accShape, stride: rowMajorStrides(accShape), data: accData}
	for k := 1; k < nd; k++ {
		core := tt.Cores[k] // (r_k, n_k, r_{k+1})
		rk := tt.Ranks[k]
		nk := tt.Shape[k]
		rk1 := tt.Ranks[k+1]
		leftCols := len(acc.data) / rk // product of previous dims
		// Contract acc last axis (r_k) with core first axis.
		coreMat, _ := NewMatrixData(rk, nk*rk1, core.data)
		accMat, _ := NewMatrixData(leftCols, rk, acc.data)
		prod, _ := accMat.Mul(coreMat) // (leftCols, nk*rk1)
		newShape := append(append([]int(nil), accShape[:len(accShape)-1]...), nk, rk1)
		acc = &Tensor{shape: newShape, stride: rowMajorStrides(newShape), data: prod.data}
		accShape = newShape
	}
	// Final acc shape is (n0,...,n_{d-1}, 1); drop trailing 1.
	final := &Tensor{shape: append([]int(nil), tt.Shape...), stride: rowMajorStrides(tt.Shape), data: acc.data}
	return final
}

// RelError returns the relative Frobenius reconstruction error of the tensor
// train against t.
func (tt *TensorTrain) RelError(t *Tensor) (float64, error) {
	return t.RelError(tt.Reconstruct())
}

// MaxRank returns the largest internal rank of the tensor train.
func (tt *TensorTrain) MaxRank() int {
	m := 0
	for _, r := range tt.Ranks {
		if r > m {
			m = r
		}
	}
	return m
}

// Compression returns the ratio of the number of dense entries to the number of
// parameters stored by the tensor train, a measure of how much the train
// compresses t.
func (tt *TensorTrain) Compression() float64 {
	dense := sizeOf(tt.Shape)
	params := 0
	for _, c := range tt.Cores {
		params += len(c.data)
	}
	if params == 0 {
		return 0
	}
	return float64(dense) / float64(params)
}
