package gaussproc

import (
	"fmt"
	"strings"
)

// SumKernel is the pointwise sum of its component kernels:
// k(x, y) = Σ Kernels[i].Eval(x, y). A sum of positive-definite kernels is
// positive definite.
type SumKernel struct {
	Kernels []Kernel
}

// NewSumKernel returns a [SumKernel] over the given components.
func NewSumKernel(kernels ...Kernel) SumKernel {
	return SumKernel{Kernels: kernels}
}

// Eval returns the sum of the component covariances of x and y.
func (k SumKernel) Eval(x, y []float64) float64 {
	var s float64
	for _, c := range k.Kernels {
		s += c.Eval(x, y)
	}
	return s
}

// String returns a human-readable description of the kernel.
func (k SumKernel) String() string {
	parts := make([]string, len(k.Kernels))
	for i, c := range k.Kernels {
		parts[i] = fmt.Sprint(c)
	}
	return "(" + strings.Join(parts, " + ") + ")"
}

// ProductKernel is the pointwise product of its component kernels:
// k(x, y) = Π Kernels[i].Eval(x, y). A product of positive-definite kernels is
// positive definite.
type ProductKernel struct {
	Kernels []Kernel
}

// NewProductKernel returns a [ProductKernel] over the given components.
func NewProductKernel(kernels ...Kernel) ProductKernel {
	return ProductKernel{Kernels: kernels}
}

// Eval returns the product of the component covariances of x and y.
func (k ProductKernel) Eval(x, y []float64) float64 {
	p := 1.0
	for _, c := range k.Kernels {
		p *= c.Eval(x, y)
	}
	return p
}

// String returns a human-readable description of the kernel.
func (k ProductKernel) String() string {
	parts := make([]string, len(k.Kernels))
	for i, c := range k.Kernels {
		parts[i] = fmt.Sprint(c)
	}
	return "(" + strings.Join(parts, " * ") + ")"
}

// ScaledKernel multiplies a base kernel by a non-negative scalar:
// k(x, y) = Factor·Kernel.Eval(x, y). For Factor ≥ 0 it preserves positive
// definiteness.
type ScaledKernel struct {
	Factor float64
	Kernel Kernel
}

// NewScaledKernel returns a [ScaledKernel] scaling base by factor.
func NewScaledKernel(factor float64, base Kernel) ScaledKernel {
	return ScaledKernel{Factor: factor, Kernel: base}
}

// Eval returns the scaled covariance of x and y.
func (k ScaledKernel) Eval(x, y []float64) float64 {
	return k.Factor * k.Kernel.Eval(x, y)
}

// String returns a human-readable description of the kernel.
func (k ScaledKernel) String() string {
	return fmt.Sprintf("%g*%v", k.Factor, k.Kernel)
}

// Sum returns a [Kernel] equal to the pointwise sum a+b.
func Sum(a, b Kernel) Kernel { return NewSumKernel(a, b) }

// Product returns a [Kernel] equal to the pointwise product a·b.
func Product(a, b Kernel) Kernel { return NewProductKernel(a, b) }

// Scale returns a [Kernel] equal to the base kernel multiplied by the scalar
// factor.
func Scale(factor float64, base Kernel) Kernel { return NewScaledKernel(factor, base) }

// AddKernels returns the pointwise sum of any number of kernels.
func AddKernels(kernels ...Kernel) Kernel { return NewSumKernel(kernels...) }

// MultiplyKernels returns the pointwise product of any number of kernels.
func MultiplyKernels(kernels ...Kernel) Kernel { return NewProductKernel(kernels...) }

// Exponentiate returns the kernel raised to a positive integer power p,
// implemented as a [ProductKernel] with p copies of base. It panics for p < 1.
func Exponentiate(base Kernel, p int) Kernel {
	if p < 1 {
		panic("gaussproc: Exponentiate requires p >= 1")
	}
	ks := make([]Kernel, p)
	for i := range ks {
		ks[i] = base
	}
	return NewProductKernel(ks...)
}
