package quadrature

import (
	"errors"
	"math"
	"sort"
)

// Func is a real-valued function of one real variable.
type Func func(x float64) float64

// Func2 is a real-valued function of two real variables.
type Func2 func(x, y float64) float64

// Func3 is a real-valued function of three real variables.
type Func3 func(x, y, z float64) float64

// FuncN is a real-valued function of a point in n-dimensional space.
type FuncN func(x []float64) float64

// Common sentinel errors returned by constructors that validate their input.
var (
	// ErrNonPositiveN is returned when a rule is requested with fewer than one
	// node.
	ErrNonPositiveN = errors.New("quadrature: number of nodes must be positive")
	// ErrBadParameter is returned when a Jacobi/Laguerre parameter is out of
	// the admissible range.
	ErrBadParameter = errors.New("quadrature: parameter out of range")
	// ErrDimMismatch is returned when slice arguments have inconsistent
	// lengths.
	ErrDimMismatch = errors.New("quadrature: dimension mismatch")
)

// Rule is a quadrature rule: a set of abscissae (Nodes) with matching Weights
// such that the integral of f is approximated by the weighted sum of f at the
// nodes. The two slices always have equal length.
type Rule struct {
	Nodes   []float64
	Weights []float64
}

// NewRule builds a Rule from parallel node and weight slices. It panics if the
// slices have different lengths, since that indicates a programming error.
func NewRule(nodes, weights []float64) Rule {
	if len(nodes) != len(weights) {
		panic("quadrature: NewRule length mismatch")
	}
	n := make([]float64, len(nodes))
	w := make([]float64, len(weights))
	copy(n, nodes)
	copy(w, weights)
	return Rule{Nodes: n, Weights: w}
}

// Len reports the number of nodes in the rule.
func (r Rule) Len() int { return len(r.Nodes) }

// Integrate applies the rule to f, returning the weighted sum of f at the
// nodes. For a rule generated on its canonical interval this is the estimate
// of the corresponding weighted integral.
func (r Rule) Integrate(f Func) float64 {
	var s float64
	for i, x := range r.Nodes {
		s += r.Weights[i] * f(x)
	}
	return s
}

// WeightSum returns the sum of the rule's weights. For a rule that integrates
// the constant weight over its interval this equals the measure of the
// interval.
func (r Rule) WeightSum() float64 {
	var s float64
	for _, w := range r.Weights {
		s += w
	}
	return s
}

// Clone returns a deep copy of the rule.
func (r Rule) Clone() Rule { return NewRule(r.Nodes, r.Weights) }

// Reversed returns a copy of the rule with nodes and weights in reverse order.
func (r Rule) Reversed() Rule {
	n := len(r.Nodes)
	nodes := make([]float64, n)
	weights := make([]float64, n)
	for i := 0; i < n; i++ {
		nodes[i] = r.Nodes[n-1-i]
		weights[i] = r.Weights[n-1-i]
	}
	return Rule{Nodes: nodes, Weights: weights}
}

// Scale maps a rule defined on the canonical interval [-1, 1] to the interval
// [a, b] by the affine change of variables x = (a+b)/2 + (b-a)/2 * t, scaling
// the weights by the Jacobian (b-a)/2. It is the standard way to move a
// Gauss-Legendre or Clenshaw-Curtis rule onto an arbitrary interval.
func (r Rule) Scale(a, b float64) Rule {
	half := 0.5 * (b - a)
	mid := 0.5 * (a + b)
	n := len(r.Nodes)
	nodes := make([]float64, n)
	weights := make([]float64, n)
	for i := 0; i < n; i++ {
		nodes[i] = mid + half*r.Nodes[i]
		weights[i] = half * r.Weights[i]
	}
	return Rule{Nodes: nodes, Weights: weights}
}

// IntegrateOn scales the canonical rule to [a, b] and applies it to f in one
// step. It is equivalent to r.Scale(a, b).Integrate(f) but avoids allocating
// the intermediate rule.
func (r Rule) IntegrateOn(f Func, a, b float64) float64 {
	half := 0.5 * (b - a)
	mid := 0.5 * (a + b)
	var s float64
	for i, t := range r.Nodes {
		s += r.Weights[i] * f(mid+half*t)
	}
	return half * s
}

// Rescale maps parallel node/weight slices from the canonical interval [-1, 1]
// to [a, b], returning fresh slices. The weights are multiplied by the
// Jacobian (b-a)/2.
func Rescale(nodes, weights []float64, a, b float64) (rn, rw []float64) {
	half := 0.5 * (b - a)
	mid := 0.5 * (a + b)
	rn = make([]float64, len(nodes))
	rw = make([]float64, len(weights))
	for i := range nodes {
		rn[i] = mid + half*nodes[i]
		rw[i] = half * weights[i]
	}
	return rn, rw
}

// symTriEig computes the eigenvalues and the first components of the
// (orthonormal) eigenvectors of the symmetric tridiagonal matrix whose main
// diagonal is diag (length n) and whose sub/super-diagonal is sub (length
// n-1). It uses the implicitly shifted QL algorithm, accumulating only the
// first row of the eigenvector matrix, which is all the Golub-Welsch formula
// requires. The returned eigenvalues are not sorted.
func symTriEig(diag, sub []float64) (eig, firstComp []float64) {
	n := len(diag)
	d := make([]float64, n)
	copy(d, diag)
	e := make([]float64, n) // e[n-1] stays zero
	for i := 0; i < n-1; i++ {
		e[i] = sub[i]
	}
	z := make([]float64, n)
	z[0] = 1
	if n == 1 {
		return d, z
	}
	const eps = 2.220446049250313e-16
	for l := 0; l < n; l++ {
		iter := 0
		for {
			var m int
			for m = l; m < n-1; m++ {
				dd := math.Abs(d[m]) + math.Abs(d[m+1])
				if math.Abs(e[m]) <= eps*dd {
					break
				}
			}
			if m == l {
				break
			}
			iter++
			if iter > 60 {
				break
			}
			g := (d[l+1] - d[l]) / (2 * e[l])
			r := math.Hypot(g, 1)
			g = d[m] - d[l] + e[l]/(g+math.Copysign(r, g))
			s, c := 1.0, 1.0
			p := 0.0
			var i int
			broke := false
			for i = m - 1; i >= l; i-- {
				f := s * e[i]
				b := c * e[i]
				r = math.Hypot(f, g)
				e[i+1] = r
				if r == 0 {
					d[i+1] -= p
					e[m] = 0
					broke = true
					break
				}
				s = f / r
				c = g / r
				g = d[i+1] - p
				r = (d[i]-g)*s + 2*c*b
				p = s * r
				d[i+1] = g + p
				g = c*r - b
				f = z[i+1]
				z[i+1] = s*z[i] + c*f
				z[i] = c*z[i] - s*f
			}
			if broke {
				continue
			}
			d[l] -= p
			e[l] = g
			e[m] = 0
		}
	}
	return d, z
}

// GolubWelsch computes Gauss quadrature nodes and weights from the monic
// three-term recurrence coefficients of the associated orthogonal polynomials.
// alpha holds the diagonal recurrence coefficients (length n); beta holds the
// off-diagonal coefficients with beta[0] equal to the zeroth moment mu0 (the
// total mass of the weight) and beta[k] equal to the recurrence coefficient
// beta_k for k >= 1. The returned nodes are sorted ascending and the weights
// are w_k = mu0 * v_{0k}^2 where v_{0k} is the first component of the k-th
// normalized eigenvector of the Jacobi matrix.
func GolubWelsch(alpha, beta []float64) (nodes, weights []float64) {
	n := len(alpha)
	sub := make([]float64, 0)
	if n > 1 {
		sub = make([]float64, n-1)
		for i := 1; i < n; i++ {
			sub[i-1] = math.Sqrt(beta[i])
		}
	}
	eig, z := symTriEig(alpha, sub)
	mu0 := beta[0]
	idx := make([]int, n)
	for i := range idx {
		idx[i] = i
	}
	sort.Slice(idx, func(i, j int) bool { return eig[idx[i]] < eig[idx[j]] })
	nodes = make([]float64, n)
	weights = make([]float64, n)
	for k, j := range idx {
		nodes[k] = eig[j]
		weights[k] = mu0 * z[j] * z[j]
	}
	return nodes, weights
}

// monicEvalPair evaluates the monic orthogonal polynomials defined by the
// recurrence coefficients alpha and beta at x, returning p_{n-1}(x) and
// p_{n-2}(x), where n = len(alpha). beta[0] is unused (it is the zeroth
// moment); beta[k] for k >= 1 is the recurrence coefficient.
func monicEvalPair(alpha, beta []float64, x float64) (pnm1, pnm2 float64) {
	n := len(alpha)
	pkm1 := 0.0 // p_{-1}
	pk := 1.0   // p_0
	for k := 0; k < n-1; k++ {
		pkp1 := (x-alpha[k])*pk - beta[k]*pkm1
		pkm1 = pk
		pk = pkp1
	}
	return pk, pkm1
}
