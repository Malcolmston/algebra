package gaussproc

import "math"

// RKHSFunction represents an element of the reproducing-kernel Hilbert space
// (RKHS) of a positive-definite kernel as a finite expansion
// f(·) = Σ Weights[i]·Kernel(Points[i], ·). Such expansions are dense in the
// RKHS and are exactly the form of the solutions returned by the representer
// theorem.
type RKHSFunction struct {
	Kernel  Kernel
	Points  [][]float64
	Weights []float64
}

// NewRKHSFunction returns an [RKHSFunction] with the given kernel, centres and
// weights. It panics if the number of centres and weights differ.
func NewRKHSFunction(kernel Kernel, points [][]float64, weights []float64) RKHSFunction {
	if len(points) != len(weights) {
		panic(ErrDimensionMismatch)
	}
	return RKHSFunction{Kernel: kernel, Points: points, Weights: weights}
}

// Eval returns the value f(x) = Σ Weights[i]·Kernel(Points[i], x).
func (f RKHSFunction) Eval(x []float64) float64 {
	var s float64
	for i, p := range f.Points {
		s += f.Weights[i] * f.Kernel.Eval(p, x)
	}
	return s
}

// EvalBatch returns f evaluated at each input in xs.
func (f RKHSFunction) EvalBatch(xs [][]float64) []float64 {
	out := make([]float64, len(xs))
	for i, x := range xs {
		out[i] = f.Eval(x)
	}
	return out
}

// InnerProduct returns the RKHS inner product ⟨f, g⟩_H, which for finite
// expansions equals Σ_i Σ_j f.Weights[i]·g.Weights[j]·Kernel(f.Points[i],
// g.Points[j]). Both functions must share the same kernel; the kernel of the
// receiver is used.
func (f RKHSFunction) InnerProduct(g RKHSFunction) float64 {
	var s float64
	for i, pi := range f.Points {
		wi := f.Weights[i]
		if wi == 0 {
			continue
		}
		for j, pj := range g.Points {
			s += wi * g.Weights[j] * f.Kernel.Eval(pi, pj)
		}
	}
	return s
}

// SquaredNorm returns the squared RKHS norm ⟨f, f⟩_H.
func (f RKHSFunction) SquaredNorm() float64 {
	return f.InnerProduct(f)
}

// Norm returns the RKHS norm ‖f‖_H = √⟨f, f⟩_H.
func (f RKHSFunction) Norm() float64 {
	return math.Sqrt(f.SquaredNorm())
}

// Scale returns the RKHS function c·f, sharing the same centres.
func (f RKHSFunction) Scale(c float64) RKHSFunction {
	w := make([]float64, len(f.Weights))
	for i := range f.Weights {
		w[i] = c * f.Weights[i]
	}
	return RKHSFunction{Kernel: f.Kernel, Points: f.Points, Weights: w}
}

// Add returns the RKHS function f+g by concatenating their expansions. Both
// functions must share the same kernel; the kernel of the receiver is used.
func (f RKHSFunction) Add(g RKHSFunction) RKHSFunction {
	points := make([][]float64, 0, len(f.Points)+len(g.Points))
	points = append(points, f.Points...)
	points = append(points, g.Points...)
	weights := make([]float64, 0, len(f.Weights)+len(g.Weights))
	weights = append(weights, f.Weights...)
	weights = append(weights, g.Weights...)
	return RKHSFunction{Kernel: f.Kernel, Points: points, Weights: weights}
}

// RKHSInnerProduct returns the RKHS inner product of two finite kernel
// expansions with weights a at centres xa and weights b at centres xb, all
// under the given kernel. It equals Σ_i Σ_j a[i]·b[j]·kernel(xa[i], xb[j]).
func RKHSInnerProduct(kernel Kernel, xa [][]float64, a []float64, xb [][]float64, b []float64) float64 {
	if len(xa) != len(a) || len(xb) != len(b) {
		panic(ErrDimensionMismatch)
	}
	var s float64
	for i, pi := range xa {
		if a[i] == 0 {
			continue
		}
		for j, pj := range xb {
			s += a[i] * b[j] * kernel.Eval(pi, pj)
		}
	}
	return s
}

// RKHSNorm returns the RKHS norm of the finite kernel expansion with weights a
// at centres x under the given kernel.
func RKHSNorm(kernel Kernel, x [][]float64, a []float64) float64 {
	return math.Sqrt(RKHSInnerProduct(kernel, x, a, x, a))
}

// RKHSDistance returns the RKHS norm of the difference f-g, computed as
// √(⟨f,f⟩ - 2⟨f,g⟩ + ⟨g,g⟩). Both functions must share the same kernel.
func RKHSDistance(f, g RKHSFunction) float64 {
	d := f.SquaredNorm() - 2*f.InnerProduct(g) + g.SquaredNorm()
	if d < 0 {
		d = 0
	}
	return math.Sqrt(d)
}

// KernelRidgeRegression fits a kernel ridge regression model to training data
// (x, y) with regularisation strength lambda ≥ 0, returning the solution as an
// [RKHSFunction]. The weights solve (K + lambda·I)·w = y, where K is the Gram
// matrix; by the representer theorem this minimises the regularised squared
// loss Σ(f(xᵢ)-yᵢ)² + lambda·‖f‖²_H.
func KernelRidgeRegression(kernel Kernel, x [][]float64, y []float64, lambda float64) (RKHSFunction, error) {
	if len(x) == 0 {
		return RKHSFunction{}, ErrEmpty
	}
	if len(x) != len(y) {
		return RKHSFunction{}, ErrDimensionMismatch
	}
	k := AddToDiagonal(GramMatrix(kernel, x), lambda)
	w, err := SolveSPD(k, y)
	if err != nil {
		return RKHSFunction{}, err
	}
	return RKHSFunction{Kernel: kernel, Points: x, Weights: w}, nil
}

// RepresenterCoefficients returns the coefficient vector w solving
// (K + lambda·I)·w = y for the kernel ridge problem, without wrapping it in an
// [RKHSFunction].
func RepresenterCoefficients(kernel Kernel, x [][]float64, y []float64, lambda float64) ([]float64, error) {
	if len(x) != len(y) {
		return nil, ErrDimensionMismatch
	}
	k := AddToDiagonal(GramMatrix(kernel, x), lambda)
	return SolveSPD(k, y)
}

// KernelMatrixRank estimates the numerical rank of the Gram matrix of kernel
// over x by counting Cholesky pivots that exceed tol, adding progressively
// larger jitter is avoided; instead a pivoted diagonal check is used.
func KernelMatrixRank(kernel Kernel, x [][]float64, tol float64) int {
	k := GramMatrix(kernel, x)
	n := k.Rows()
	// Gaussian elimination with partial diagonal pivoting to count nonzero
	// pivots of the symmetric matrix.
	a := k.Clone()
	rank := 0
	used := make([]bool, n)
	for step := 0; step < n; step++ {
		// find pivot: largest remaining diagonal entry.
		piv := -1
		best := tol
		for i := 0; i < n; i++ {
			if used[i] {
				continue
			}
			if math.Abs(a[i][i]) > best {
				best = math.Abs(a[i][i])
				piv = i
			}
		}
		if piv < 0 {
			break
		}
		used[piv] = true
		rank++
		d := a[piv][piv]
		for i := 0; i < n; i++ {
			if used[i] {
				continue
			}
			f := a[i][piv] / d
			for j := 0; j < n; j++ {
				a[i][j] -= f * a[piv][j]
			}
		}
	}
	return rank
}
