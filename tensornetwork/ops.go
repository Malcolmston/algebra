package tensornetwork

import (
	"errors"
	"math"
)

// binaryElementwise applies f entrywise to two identically shaped tensors.
func binaryElementwise(a, b *Tensor, f func(x, y float64) float64) (*Tensor, error) {
	if !shapeEqual(a.shape, b.shape) {
		return nil, errors.New("tensornetwork: shape mismatch in elementwise op")
	}
	out := New(a.shape...)
	for i := range a.data {
		out.data[i] = f(a.data[i], b.data[i])
	}
	return out, nil
}

// Add returns the entrywise sum t+o. It returns an error on a shape mismatch.
func (t *Tensor) Add(o *Tensor) (*Tensor, error) {
	return binaryElementwise(t, o, func(x, y float64) float64 { return x + y })
}

// Sub returns the entrywise difference t-o. It returns an error on a shape
// mismatch.
func (t *Tensor) Sub(o *Tensor) (*Tensor, error) {
	return binaryElementwise(t, o, func(x, y float64) float64 { return x - y })
}

// Mul returns the entrywise (Hadamard) product t∘o. It returns an error on a
// shape mismatch.
func (t *Tensor) Mul(o *Tensor) (*Tensor, error) {
	return binaryElementwise(t, o, func(x, y float64) float64 { return x * y })
}

// Div returns the entrywise quotient t/o. It returns an error on a shape
// mismatch.
func (t *Tensor) Div(o *Tensor) (*Tensor, error) {
	return binaryElementwise(t, o, func(x, y float64) float64 { return x / y })
}

// Maximum returns the entrywise maximum of t and o.
func (t *Tensor) Maximum(o *Tensor) (*Tensor, error) {
	return binaryElementwise(t, o, math.Max)
}

// Minimum returns the entrywise minimum of t and o.
func (t *Tensor) Minimum(o *Tensor) (*Tensor, error) {
	return binaryElementwise(t, o, math.Min)
}

// Scale returns t with every entry multiplied by s.
func (t *Tensor) Scale(s float64) *Tensor {
	out := t.Clone()
	for i := range out.data {
		out.data[i] *= s
	}
	return out
}

// AddScalar returns t with s added to every entry.
func (t *Tensor) AddScalar(s float64) *Tensor {
	out := t.Clone()
	for i := range out.data {
		out.data[i] += s
	}
	return out
}

// Neg returns the entrywise negation of t.
func (t *Tensor) Neg() *Tensor { return t.Scale(-1) }

// Apply returns a new tensor with f applied to every entry of t.
func (t *Tensor) Apply(f func(float64) float64) *Tensor {
	out := t.Clone()
	for i := range out.data {
		out.data[i] = f(out.data[i])
	}
	return out
}

// Abs returns the entrywise absolute value of t.
func (t *Tensor) Abs() *Tensor { return t.Apply(math.Abs) }

// Sqrt returns the entrywise square root of t.
func (t *Tensor) Sqrt() *Tensor { return t.Apply(math.Sqrt) }

// Exp returns the entrywise exponential of t.
func (t *Tensor) Exp() *Tensor { return t.Apply(math.Exp) }

// Log returns the entrywise natural logarithm of t.
func (t *Tensor) Log() *Tensor { return t.Apply(math.Log) }

// Sum returns the sum of all entries of t.
func (t *Tensor) Sum() float64 {
	s := 0.0
	for _, v := range t.data {
		s += v
	}
	return s
}

// Product returns the product of all entries of t.
func (t *Tensor) Product() float64 {
	p := 1.0
	for _, v := range t.data {
		p *= v
	}
	return p
}

// Mean returns the arithmetic mean of all entries of t.
func (t *Tensor) Mean() float64 {
	if len(t.data) == 0 {
		return 0
	}
	return t.Sum() / float64(len(t.data))
}

// Max returns the largest entry of t. It panics on an empty tensor.
func (t *Tensor) Max() float64 {
	if len(t.data) == 0 {
		panic("tensornetwork: Max of empty tensor")
	}
	m := t.data[0]
	for _, v := range t.data[1:] {
		if v > m {
			m = v
		}
	}
	return m
}

// Min returns the smallest entry of t. It panics on an empty tensor.
func (t *Tensor) Min() float64 {
	if len(t.data) == 0 {
		panic("tensornetwork: Min of empty tensor")
	}
	m := t.data[0]
	for _, v := range t.data[1:] {
		if v < m {
			m = v
		}
	}
	return m
}

// ArgMax returns the flat index of the largest entry of t.
func (t *Tensor) ArgMax() int {
	best := 0
	for i, v := range t.data {
		if v > t.data[best] {
			best = i
		}
	}
	return best
}

// ArgMin returns the flat index of the smallest entry of t.
func (t *Tensor) ArgMin() int {
	best := 0
	for i, v := range t.data {
		if v < t.data[best] {
			best = i
		}
	}
	return best
}

// Norm returns the Frobenius (Euclidean) norm of t: the square root of the sum
// of squares of all entries.
func (t *Tensor) Norm() float64 {
	s := 0.0
	for _, v := range t.data {
		s += v * v
	}
	return math.Sqrt(s)
}

// NormP returns the entrywise p-norm of t for p >= 1. NormP(1) is the sum of
// absolute values and NormP(2) equals [Tensor.Norm].
func (t *Tensor) NormP(p float64) float64 {
	if math.IsInf(p, 1) {
		m := 0.0
		for _, v := range t.data {
			if a := math.Abs(v); a > m {
				m = a
			}
		}
		return m
	}
	s := 0.0
	for _, v := range t.data {
		s += math.Pow(math.Abs(v), p)
	}
	return math.Pow(s, 1/p)
}

// reduceAxis reduces t along axis using the accumulator init and combine, then
// removing the axis.
func (t *Tensor) reduceAxis(axis int, init float64, combine func(acc, x float64) float64) (*Tensor, error) {
	n := len(t.shape)
	if axis < 0 || axis >= n {
		return nil, errors.New("tensornetwork: axis out of range")
	}
	var sh []int
	for i, d := range t.shape {
		if i != axis {
			sh = append(sh, d)
		}
	}
	out := Full(init, sh...)
	idx := make([]int, n)
	redStride := t.stride[axis]
	for flat := 0; flat < len(t.data); flat++ {
		rem := flat
		for a := 0; a < n; a++ {
			idx[a] = (rem / t.stride[a]) % t.shape[a]
		}
		_ = redStride
		var outIdx []int
		for a := 0; a < n; a++ {
			if a != axis {
				outIdx = append(outIdx, idx[a])
			}
		}
		off, _ := out.flatIndex(outIdx)
		out.data[off] = combine(out.data[off], t.data[flat])
	}
	return out, nil
}

// SumAxis returns t summed along the given axis, which is removed from the
// result.
func (t *Tensor) SumAxis(axis int) (*Tensor, error) {
	return t.reduceAxis(axis, 0, func(acc, x float64) float64 { return acc + x })
}

// ProdAxis returns t multiplied along the given axis, which is removed from the
// result.
func (t *Tensor) ProdAxis(axis int) (*Tensor, error) {
	return t.reduceAxis(axis, 1, func(acc, x float64) float64 { return acc * x })
}

// MaxAxis returns the maximum of t along the given axis, which is removed from
// the result.
func (t *Tensor) MaxAxis(axis int) (*Tensor, error) {
	return t.reduceAxis(axis, math.Inf(-1), math.Max)
}

// MinAxis returns the minimum of t along the given axis, which is removed from
// the result.
func (t *Tensor) MinAxis(axis int) (*Tensor, error) {
	return t.reduceAxis(axis, math.Inf(1), math.Min)
}

// MeanAxis returns the mean of t along the given axis, which is removed from the
// result.
func (t *Tensor) MeanAxis(axis int) (*Tensor, error) {
	s, err := t.SumAxis(axis)
	if err != nil {
		return nil, err
	}
	return s.Scale(1 / float64(t.shape[axis])), nil
}

// Normalize returns t scaled so that its Frobenius norm is 1, together with the
// original norm. If the norm is zero t is returned unchanged with norm 0.
func (t *Tensor) Normalize() (*Tensor, float64) {
	nrm := t.Norm()
	if nrm == 0 {
		return t.Clone(), 0
	}
	return t.Scale(1 / nrm), nrm
}

// Dot returns the full inner product of t and o, the sum over all positions of
// the entrywise products. It returns an error on a shape mismatch.
func (t *Tensor) Dot(o *Tensor) (float64, error) {
	if !shapeEqual(t.shape, o.shape) {
		return 0, errors.New("tensornetwork: shape mismatch in Dot")
	}
	s := 0.0
	for i := range t.data {
		s += t.data[i] * o.data[i]
	}
	return s, nil
}

// RelError returns the relative Frobenius error ‖t-o‖/‖t‖ between t and o. If t
// has zero norm the absolute error ‖t-o‖ is returned. It returns an error on a
// shape mismatch.
func (t *Tensor) RelError(o *Tensor) (float64, error) {
	d, err := t.Sub(o)
	if err != nil {
		return 0, err
	}
	nt := t.Norm()
	if nt == 0 {
		return d.Norm(), nil
	}
	return d.Norm() / nt, nil
}
