package probability

import "math"

// JointDistribution is a bivariate discrete distribution over a grid of
// outcomes. X holds the distinct outcome values of the first variable (rows) and
// Y holds those of the second (columns); P[i][j] is the joint probability
// P(X = X[i], Y = Y[j]). The probabilities are non-negative and sum to one.
type JointDistribution struct {
	// X holds the distinct outcome values of the first variable.
	X []float64
	// Y holds the distinct outcome values of the second variable.
	Y []float64
	// P is the joint probability grid indexed as P[i][j] = P(X[i], Y[j]).
	P [][]float64
}

// NewJointDistribution builds a JointDistribution and validates it: X and Y must
// be non-empty, P must be a len(X)-by-len(Y) grid of non-negative finite
// probabilities summing to one within [probabilityTol]. The inputs are copied.
func NewJointDistribution(x, y []float64, p [][]float64) (JointDistribution, error) {
	if len(x) == 0 || len(y) == 0 {
		return JointDistribution{}, probabilityErrorf("NewJointDistribution: empty support")
	}
	if len(p) != len(x) {
		return JointDistribution{}, probabilityErrorf("NewJointDistribution: P has %d rows, want %d", len(p), len(x))
	}
	total := 0.0
	xc := make([]float64, len(x))
	yc := make([]float64, len(y))
	copy(xc, x)
	copy(yc, y)
	pc := make([][]float64, len(x))
	for i := range p {
		if len(p[i]) != len(y) {
			return JointDistribution{}, probabilityErrorf("NewJointDistribution: row %d has %d columns, want %d", i, len(p[i]), len(y))
		}
		pc[i] = make([]float64, len(y))
		for j, v := range p[i] {
			if v < 0 || math.IsNaN(v) || math.IsInf(v, 0) {
				return JointDistribution{}, probabilityErrorf("NewJointDistribution: invalid probability %g at (%d,%d)", v, i, j)
			}
			pc[i][j] = v
			total += v
		}
	}
	if probabilityAbs(total-1) > probabilityTol {
		return JointDistribution{}, probabilityErrorf("NewJointDistribution: probabilities sum to %g, not 1", total)
	}
	return JointDistribution{X: xc, Y: yc, P: pc}, nil
}

// IndependentJoint builds the joint distribution of two independent random
// variables from their marginals, setting P[i][j] = a.Probs[i]·b.Probs[j]. The
// resulting joint is independent by construction.
func IndependentJoint(a, b Distribution) JointDistribution {
	p := make([][]float64, len(a.Outcomes))
	for i := range a.Outcomes {
		p[i] = make([]float64, len(b.Outcomes))
		for j := range b.Outcomes {
			p[i][j] = a.Probs[i] * b.Probs[j]
		}
	}
	x := make([]float64, len(a.Outcomes))
	y := make([]float64, len(b.Outcomes))
	copy(x, a.Outcomes)
	copy(y, b.Outcomes)
	return JointDistribution{X: x, Y: y, P: p}
}

// MarginalX returns the marginal distribution of X obtained by summing the joint
// probabilities across the Y axis for each row.
func (j JointDistribution) MarginalX() Distribution {
	probs := make([]float64, len(j.X))
	for i := range j.X {
		s := 0.0
		for k := range j.Y {
			s += j.P[i][k]
		}
		probs[i] = s
	}
	outs := make([]float64, len(j.X))
	copy(outs, j.X)
	mo, mp := probabilityMerge(outs, probs)
	return Distribution{Outcomes: mo, Probs: mp}
}

// MarginalY returns the marginal distribution of Y obtained by summing the joint
// probabilities across the X axis for each column.
func (j JointDistribution) MarginalY() Distribution {
	probs := make([]float64, len(j.Y))
	for k := range j.Y {
		s := 0.0
		for i := range j.X {
			s += j.P[i][k]
		}
		probs[k] = s
	}
	outs := make([]float64, len(j.Y))
	copy(outs, j.Y)
	mo, mp := probabilityMerge(outs, probs)
	return Distribution{Outcomes: mo, Probs: mp}
}

// ConditionalXGivenY returns the conditional distribution of X given Y = Y[k],
// i.e. P(X = X[i] | Y = Y[k]) = P[i][k] / P(Y = Y[k]). It returns an error if k
// is out of range or the conditioning event has zero probability.
func (j JointDistribution) ConditionalXGivenY(k int) (Distribution, error) {
	if k < 0 || k >= len(j.Y) {
		return Distribution{}, probabilityErrorf("ConditionalXGivenY: column %d out of range", k)
	}
	denom := 0.0
	for i := range j.X {
		denom += j.P[i][k]
	}
	if denom <= 0 {
		return Distribution{}, probabilityErrorf("ConditionalXGivenY: P(Y=%g) is zero", j.Y[k])
	}
	outs := make([]float64, len(j.X))
	probs := make([]float64, len(j.X))
	for i := range j.X {
		outs[i] = j.X[i]
		probs[i] = j.P[i][k] / denom
	}
	mo, mp := probabilityMerge(outs, probs)
	return Distribution{Outcomes: mo, Probs: mp}, nil
}

// ConditionalYGivenX returns the conditional distribution of Y given X = X[i],
// i.e. P(Y = Y[k] | X = X[i]) = P[i][k] / P(X = X[i]). It returns an error if i
// is out of range or the conditioning event has zero probability.
func (j JointDistribution) ConditionalYGivenX(i int) (Distribution, error) {
	if i < 0 || i >= len(j.X) {
		return Distribution{}, probabilityErrorf("ConditionalYGivenX: row %d out of range", i)
	}
	denom := 0.0
	for k := range j.Y {
		denom += j.P[i][k]
	}
	if denom <= 0 {
		return Distribution{}, probabilityErrorf("ConditionalYGivenX: P(X=%g) is zero", j.X[i])
	}
	outs := make([]float64, len(j.Y))
	probs := make([]float64, len(j.Y))
	for k := range j.Y {
		outs[k] = j.Y[k]
		probs[k] = j.P[i][k] / denom
	}
	mo, mp := probabilityMerge(outs, probs)
	return Distribution{Outcomes: mo, Probs: mp}, nil
}

// ExpectationXY returns E[g(X, Y)] = Σ_{i,k} g(X[i], Y[k]) · P[i][k] for an
// arbitrary function g of the two variables.
func (j JointDistribution) ExpectationXY(g func(x, y float64) float64) float64 {
	sum := 0.0
	for i := range j.X {
		for k := range j.Y {
			sum += g(j.X[i], j.Y[k]) * j.P[i][k]
		}
	}
	return sum
}

// Covariance returns Cov(X, Y) = E[XY] - E[X]·E[Y], the covariance of the two
// variables under the joint distribution.
func (j JointDistribution) Covariance() float64 {
	exy := j.ExpectationXY(func(x, y float64) float64 { return x * y })
	return exy - j.MarginalX().Mean()*j.MarginalY().Mean()
}

// Correlation returns the Pearson correlation coefficient
// Cov(X, Y) / (σ_X · σ_Y), a value in [-1, 1]. It returns NaN when either
// marginal has zero standard deviation.
func (j JointDistribution) Correlation() float64 {
	sx := j.MarginalX().StdDev()
	sy := j.MarginalY().StdDev()
	if sx == 0 || sy == 0 {
		return math.NaN()
	}
	return j.Covariance() / (sx * sy)
}

// Independent reports whether X and Y are independent, i.e. whether
// P[i][k] equals P(X = X[i])·P(Y = Y[k]) for every cell within [probabilityTol].
func (j JointDistribution) Independent() bool {
	mx := j.MarginalX()
	my := j.MarginalY()
	for i := range j.X {
		for k := range j.Y {
			if probabilityAbs(j.P[i][k]-mx.Probs[i]*my.Probs[k]) > probabilityTol {
				return false
			}
		}
	}
	return true
}
