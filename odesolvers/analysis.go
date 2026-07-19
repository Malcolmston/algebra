package odesolvers

import "math"

// RowSums returns the sum of each row of the tableau's A matrix. For a
// consistent Runge-Kutta method these equal the nodes C.
func (bt *ButcherTableau) RowSums() []float64 {
	out := make([]float64, len(bt.A))
	for i := range bt.A {
		var s float64
		for _, v := range bt.A[i] {
			s += v
		}
		out[i] = s
	}
	return out
}

// SumB returns the sum of the primary weights B, which must equal 1 for a
// consistent method.
func (bt *ButcherTableau) SumB() float64 {
	var s float64
	for _, v := range bt.B {
		s += v
	}
	return s
}

// IsConsistent reports whether the tableau satisfies the first-order
// consistency conditions to within tol: the weights B sum to 1 and every row of
// A sums to the corresponding node C.
func (bt *ButcherTableau) IsConsistent(tol float64) bool {
	if math.Abs(bt.SumB()-1) > tol {
		return false
	}
	rs := bt.RowSums()
	for i := range rs {
		if math.Abs(rs[i]-bt.C[i]) > tol {
			return false
		}
	}
	return true
}

// RichardsonExtrapolate combines two approximations coarse and fine, computed
// with step sizes differing by the factor ratio (>1) for a method of the given
// order p, into a higher-order estimate. With T(h) = A + C h^p + ..., the
// extrapolated value is (ratio^p * fine - coarse) / (ratio^p - 1).
func RichardsonExtrapolate(coarse, fine []float64, ratio float64, order int) []float64 {
	f := math.Pow(ratio, float64(order))
	out := make([]float64, len(fine))
	for i := range fine {
		out[i] = (f*fine[i] - coarse[i]) / (f - 1)
	}
	return out
}

// EstimateOrder estimates the empirical convergence order of a fixed-step
// integrator from three errors measured at step sizes h, h/2 and h/4:
// order ≈ log2(e(h)/e(h/2)) using the two coarser errors. The finer error is
// used to average two consecutive ratios for robustness. All errors must be
// positive.
func EstimateOrder(errH, errH2, errH4 float64) float64 {
	r1 := math.Log2(errH / errH2)
	r2 := math.Log2(errH2 / errH4)
	return 0.5 * (r1 + r2)
}

// GlobalError returns the infinity-norm difference between the final state of a
// solution and a known exact reference state.
func GlobalError(sol *Solution, exact []float64) float64 {
	return NormInf(Sub(sol.Final(), exact))
}

// KuttaThirdOrderTableau returns Kutta's classical third-order explicit
// Runge-Kutta tableau.
func KuttaThirdOrderTableau() *ButcherTableau {
	return &ButcherTableau{
		A: [][]float64{
			{0, 0, 0},
			{0.5, 0, 0},
			{-1, 2, 0},
		},
		B:     []float64{1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0},
		C:     []float64{0, 0.5, 1},
		Order: 3,
		Name:  "Kutta 3",
	}
}

// HeunThirdOrderTableau returns Heun's third-order explicit Runge-Kutta
// tableau.
func HeunThirdOrderTableau() *ButcherTableau {
	return &ButcherTableau{
		A: [][]float64{
			{0, 0, 0},
			{1.0 / 3.0, 0, 0},
			{0, 2.0 / 3.0, 0},
		},
		B:     []float64{1.0 / 4.0, 0, 3.0 / 4.0},
		C:     []float64{0, 1.0 / 3.0, 2.0 / 3.0},
		Order: 3,
		Name:  "Heun 3",
	}
}

// SolveKutta3 integrates with Kutta's third-order method and a fixed step h.
func SolveKutta3(f Field, t0 float64, y0 []float64, tEnd, h float64) *Solution {
	return IntegrateFixed(f, KuttaThirdOrderTableau(), t0, y0, tEnd, h)
}

// GaussLegendre4Step advances one two-stage Gauss-Legendre (order 4) step.
func GaussLegendre4Step(f Field, t float64, y []float64, h float64) ([]float64, error) {
	return ImplicitRKStep(f, GaussLegendre4Tableau(), t, y, h)
}

// RadauIIA3Step advances one order-3 Radau IIA step and returns the new state.
func RadauIIA3Step(f Field, t float64, y []float64, h float64) ([]float64, error) {
	return ImplicitRKStep(f, RadauIIA3Tableau(), t, y, h)
}
