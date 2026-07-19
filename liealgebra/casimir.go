package liealgebra

import "math"

// highestWeightRootBasis converts Dynkin labels (a_1,...,a_n), meaning the
// highest weight λ = Σ a_i ω_i, into coordinates in the simple-root basis using
// the inverse Cartan matrix.
func highestWeightRootBasis(inv *Matrix, dynkin []int) []float64 {
	n := inv.Rows
	out := make([]float64, n)
	for k := 0; k < n; k++ {
		s := 0.0
		for i := 0; i < n; i++ {
			s += float64(dynkin[i]) * inv.Data[i*n+k]
		}
		out[k] = s
	}
	return out
}

// WeylDimension returns the dimension of the irreducible representation of the
// given Dynkin type with highest weight specified by nonnegative Dynkin labels
// (a_1,...,a_rank), via the Weyl dimension formula
//
//	dim = Π_{α>0} (λ+ρ, α) / (ρ, α).
//
// The result is rounded to the nearest integer. It returns [ErrDim] if the
// number of labels does not match the rank.
func WeylDimension(family string, rank int, dynkin []int) (int, error) {
	if len(dynkin) != rank {
		return 0, ErrDim
	}
	inv, err := InverseCartanMatrix(family, rank)
	if err != nil {
		return 0, err
	}
	b, err := RootBilinearForm(family, rank)
	if err != nil {
		return 0, err
	}
	rho, err := WeylVectorRootBasis(family, rank)
	if err != nil {
		return 0, err
	}
	lambda := highestWeightRootBasis(inv, dynkin)
	lamRho, _ := VecAdd(lambda, rho)
	posLabels, err := PositiveRootLabels(family, rank)
	if err != nil {
		return 0, err
	}
	prod := 1.0
	for _, lab := range posLabels {
		alpha := make([]float64, len(lab))
		for i, x := range lab {
			alpha[i] = float64(x)
		}
		num, _ := RootInnerProduct(b, lamRho, alpha)
		den, _ := RootInnerProduct(b, rho, alpha)
		prod *= num / den
	}
	return int(math.Round(prod)), nil
}

// CasimirEigenvalue returns the eigenvalue of the quadratic Casimir operator on
// the irreducible representation with the given highest weight (Dynkin labels),
// computed as (λ, λ+2ρ) in the normalisation where the long roots have squared
// length 2. It returns [ErrDim] if the label count does not match the rank.
func CasimirEigenvalue(family string, rank int, dynkin []int) (float64, error) {
	if len(dynkin) != rank {
		return 0, ErrDim
	}
	inv, err := InverseCartanMatrix(family, rank)
	if err != nil {
		return 0, err
	}
	b, err := RootBilinearForm(family, rank)
	if err != nil {
		return 0, err
	}
	rho, err := WeylVectorRootBasis(family, rank)
	if err != nil {
		return 0, err
	}
	lambda := highestWeightRootBasis(inv, dynkin)
	twoRho := VecScale(rho, 2)
	lamTwoRho, _ := VecAdd(lambda, twoRho)
	return RootInnerProduct(b, lambda, lamTwoRho)
}

// CasimirSU2 returns the physics-convention quadratic Casimir eigenvalue j(j+1)
// of the spin-j irreducible representation of su(2).
func CasimirSU2(j float64) float64 { return j * (j + 1) }

// DimensionSU2 returns the dimension 2j+1 of the spin-j representation of su(2).
// It returns [ErrRange] if 2j is not a nonnegative integer.
func DimensionSU2(j float64) (int, error) {
	twoj := int(math.Round(2 * j))
	if twoj < 0 || math.Abs(float64(twoj)-2*j) > 1e-9 {
		return 0, ErrRange
	}
	return twoj + 1, nil
}

// DimensionSLn returns the dimension n^2-1 of the Lie algebra sl(n) = A_{n-1}.
func DimensionSLn(n int) (int, error) {
	if n < 1 {
		return 0, ErrRange
	}
	return n*n - 1, nil
}

// DimensionSOn returns the dimension n(n-1)/2 of the Lie algebra so(n).
func DimensionSOn(n int) (int, error) {
	if n < 1 {
		return 0, ErrRange
	}
	return n * (n - 1) / 2, nil
}

// DimensionSPn returns the dimension n(2n+1) of the Lie algebra sp(2n) = C_n.
func DimensionSPn(n int) (int, error) {
	if n < 1 {
		return 0, ErrRange
	}
	return n * (2*n + 1), nil
}

// DimensionSUn returns the dimension n^2-1 of the Lie algebra su(n).
func DimensionSUn(n int) (int, error) {
	if n < 1 {
		return 0, ErrRange
	}
	return n*n - 1, nil
}

// AdjointDimension returns the dimension of the adjoint representation, which
// equals the dimension of the Lie algebra itself.
func AdjointDimension(family string, rank int) (int, error) {
	return LieAlgebraDimension(family, rank)
}

// DualCoxeterFromRoots returns the dual Coxeter number computed as
// 1 + (ρ, θ∨) where θ is the highest root; provided as an independent
// cross-check of [DualCoxeterNumber] for classical and exceptional types.
func DualCoxeterFromRoots(family string, rank int) (int, error) {
	b, err := RootBilinearForm(family, rank)
	if err != nil {
		return 0, err
	}
	rho, err := WeylVectorRootBasis(family, rank)
	if err != nil {
		return 0, err
	}
	theta, err := HighestRootLabel(family, rank)
	if err != nil {
		return 0, err
	}
	thetaF := make([]float64, len(theta))
	for i, x := range theta {
		thetaF[i] = float64(x)
	}
	// θ∨ = 2θ/(θ,θ).
	tt, _ := RootInnerProduct(b, thetaF, thetaF)
	coTheta := VecScale(thetaF, 2/tt)
	rt, _ := RootInnerProduct(b, rho, coTheta)
	return int(math.Round(rt)) + 1, nil
}
