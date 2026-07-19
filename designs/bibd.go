package designs

import "errors"

// BIBDParams holds the five parameters of a balanced incomplete block design
// (BIBD, or 2-design): v points, b blocks, each point in r blocks, each block
// of size k, and every pair of points together in lambda blocks.
type BIBDParams struct {
	V, B, R, K, Lambda int
}

// NewBIBDParams constructs a parameter set from v, k and lambda, deriving the
// replication number r = lambda*(v-1)/(k-1) and block count b = v*r/k. It
// reports an error when 2<=k<v does not hold or the divisibility conditions
// required for these to be integers fail.
func NewBIBDParams(v, k, lambda int) (BIBDParams, error) {
	if v < 2 || k < 2 || k >= v || lambda < 1 {
		return BIBDParams{}, errors.New("designs: require lambda>=1 and 2<=k<v")
	}
	num := lambda * (v - 1)
	if num%(k-1) != 0 {
		return BIBDParams{}, errors.New("designs: lambda*(v-1) not divisible by k-1")
	}
	r := num / (k - 1)
	if (v*r)%k != 0 {
		return BIBDParams{}, errors.New("designs: v*r not divisible by k")
	}
	b := v * r / k
	return BIBDParams{V: v, B: b, R: r, K: k, Lambda: lambda}, nil
}

// Valid reports whether the parameters satisfy the two fundamental counting
// identities r*(k-1) = lambda*(v-1) and b*k = v*r, along with the basic range
// constraints 2<=k<v, lambda>=1.
func (p BIBDParams) Valid() bool {
	if p.V < 2 || p.K < 2 || p.K >= p.V || p.Lambda < 1 || p.R < 1 || p.B < 1 {
		return false
	}
	if p.R*(p.K-1) != p.Lambda*(p.V-1) {
		return false
	}
	if p.B*p.K != p.V*p.R {
		return false
	}
	return true
}

// SatisfiesNecessary reports whether the parameters meet the necessary
// existence conditions: the divisibility identities (Valid) together with
// Fisher's inequality b>=v.
func (p BIBDParams) SatisfiesNecessary() bool {
	return p.Valid() && p.FisherInequalityHolds()
}

// FisherInequalityHolds reports whether Fisher's inequality b>=v holds for the
// parameters. Every BIBD with more than one block per point satisfies it.
func (p BIBDParams) FisherInequalityHolds() bool { return p.B >= p.V }

// IsSymmetric reports whether the design is symmetric, meaning b==v
// (equivalently r==k).
func (p BIBDParams) IsSymmetric() bool { return p.B == p.V }

// Complement returns the parameters of the complementary design, obtained by
// replacing each block with its complement. It is a
// (v, b, b-r, v-k, b-2r+lambda)-design and reports an error when the resulting
// block size is not at least 2.
func (p BIBDParams) Complement() (BIBDParams, error) {
	cv := p.V
	cb := p.B
	cr := p.B - p.R
	ck := p.V - p.K
	cl := p.B - 2*p.R + p.Lambda
	if ck < 2 || cl < 1 {
		return BIBDParams{}, errors.New("designs: complement is not a 2-design")
	}
	return BIBDParams{V: cv, B: cb, R: cr, K: ck, Lambda: cl}, nil
}

// Derived returns the parameters of the derived design of a symmetric design
// with respect to a block: it is a (k, k-1... )-design on the k points of a
// block. Concretely a symmetric 2-(v,k,lambda) design yields a derived
// 2-(k, lambda, lambda-1) design. It reports an error when the design is not
// symmetric or the result is degenerate.
func (p BIBDParams) Derived() (BIBDParams, error) {
	if !p.IsSymmetric() {
		return BIBDParams{}, errors.New("designs: derived design requires a symmetric design")
	}
	v := p.K
	k := p.Lambda
	lambda := p.Lambda - 1
	if k < 2 || k >= v || lambda < 1 {
		return BIBDParams{}, errors.New("designs: derived design is degenerate")
	}
	return NewBIBDParams(v, k, lambda)
}

// Residual returns the parameters of the residual design of a symmetric
// 2-(v,k,lambda) design: a 2-(v-k, k-lambda, lambda) design. It reports an
// error when the design is not symmetric or the result is degenerate.
func (p BIBDParams) Residual() (BIBDParams, error) {
	if !p.IsSymmetric() {
		return BIBDParams{}, errors.New("designs: residual design requires a symmetric design")
	}
	v := p.V - p.K
	k := p.K - p.Lambda
	lambda := p.Lambda
	if k < 2 || k >= v || lambda < 1 {
		return BIBDParams{}, errors.New("designs: residual design is degenerate")
	}
	return NewBIBDParams(v, k, lambda)
}

// BlockIntersection returns the number of points common to any two distinct
// blocks of a symmetric design, which is the constant lambda. It reports an
// error when the design is not symmetric.
func (p BIBDParams) BlockIntersection() (int, error) {
	if !p.IsSymmetric() {
		return 0, errors.New("designs: block intersection constant requires a symmetric design")
	}
	return p.Lambda, nil
}

// Order returns the order n = r - lambda of the design, a quantity that governs
// much of the theory of symmetric designs.
func (p BIBDParams) Order() int { return p.R - p.Lambda }

// FisherInequality reports whether the design d satisfies Fisher's inequality
// b>=v. It reports an error when d is not a 2-design.
func FisherInequality(d *Design) (bool, error) {
	p, err := d.Parameters()
	if err != nil {
		return false, err
	}
	return p.FisherInequalityHolds(), nil
}

// ProjectivePlaneParams returns the BIBD parameters 2-(n^2+n+1, n+1, 1) of a
// projective plane of order n, for n>=2.
func ProjectivePlaneParams(n int) (BIBDParams, error) {
	if n < 2 {
		return BIBDParams{}, errors.New("designs: plane order must be at least 2")
	}
	return NewBIBDParams(n*n+n+1, n+1, 1)
}

// AffinePlaneParams returns the BIBD parameters 2-(n^2, n, 1) of an affine
// plane of order n, for n>=2.
func AffinePlaneParams(n int) (BIBDParams, error) {
	if n < 2 {
		return BIBDParams{}, errors.New("designs: plane order must be at least 2")
	}
	return NewBIBDParams(n*n, n, 1)
}

// SteinerTripleParams returns the BIBD parameters 2-(v,3,1) of a Steiner triple
// system on v points. It reports an error unless v is congruent to 1 or 3
// modulo 6 (the admissible orders).
func SteinerTripleParams(v int) (BIBDParams, error) {
	if v%6 != 1 && v%6 != 3 {
		return BIBDParams{}, errors.New("designs: STS requires v = 1 or 3 (mod 6)")
	}
	return NewBIBDParams(v, 3, 1)
}

// HadamardDesignParams returns the parameters 2-(4t-1, 2t-1, t-1) of the
// symmetric design associated with a Hadamard matrix of order 4t, for t>=2.
func HadamardDesignParams(t int) (BIBDParams, error) {
	if t < 2 {
		return BIBDParams{}, errors.New("designs: require t>=2")
	}
	return NewBIBDParams(4*t-1, 2*t-1, t-1)
}
