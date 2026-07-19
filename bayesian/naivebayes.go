package bayesian

import (
	"errors"
	"math"
	"sort"
)

// ErrNotFitted is returned by classifier prediction methods called before Fit.
var ErrNotFitted = errors.New("bayesian: classifier has not been fitted")

// ErrShape is returned when input matrices and label vectors have inconsistent
// or empty shapes.
var ErrShape = errors.New("bayesian: inconsistent input shape")

// uniqueSortedClasses returns the sorted distinct labels and a map from label
// to its dense index.
func uniqueSortedClasses(y []int) ([]int, map[int]int) {
	seen := map[int]bool{}
	for _, v := range y {
		seen[v] = true
	}
	classes := make([]int, 0, len(seen))
	for v := range seen {
		classes = append(classes, v)
	}
	sort.Ints(classes)
	idx := make(map[int]int, len(classes))
	for i, c := range classes {
		idx[c] = i
	}
	return classes, idx
}

// ------------------------------------------------------------------
// Gaussian naive Bayes
// ------------------------------------------------------------------

// GaussianNB is a Gaussian naive-Bayes classifier: each continuous feature is
// modeled as a class-conditional normal distribution with features assumed
// conditionally independent given the class.
type GaussianNB struct {
	classes  []int
	priors   []float64
	means    [][]float64
	vars     [][]float64
	nFeat    int
	varSmall float64
	fitted   bool
}

// NewGaussianNB returns an unfitted Gaussian naive-Bayes classifier.
func NewGaussianNB() *GaussianNB { return &GaussianNB{} }

// Fit estimates class priors and per-feature class-conditional means and
// variances from the training matrix X (rows are samples) and labels y. A small
// variance floor proportional to the largest feature variance is added for
// numerical stability.
func (g *GaussianNB) Fit(X [][]float64, y []int) error {
	if len(X) == 0 || len(X) != len(y) {
		return ErrShape
	}
	nFeat := len(X[0])
	if nFeat == 0 {
		return ErrShape
	}
	classes, idx := uniqueSortedClasses(y)
	k := len(classes)
	counts := make([]int, k)
	means := make([][]float64, k)
	vars := make([][]float64, k)
	for c := 0; c < k; c++ {
		means[c] = make([]float64, nFeat)
		vars[c] = make([]float64, nFeat)
	}
	for i, row := range X {
		if len(row) != nFeat {
			return ErrShape
		}
		c := idx[y[i]]
		counts[c]++
		for j, v := range row {
			means[c][j] += v
		}
	}
	for c := 0; c < k; c++ {
		if counts[c] > 0 {
			for j := 0; j < nFeat; j++ {
				means[c][j] /= float64(counts[c])
			}
		}
	}
	for i, row := range X {
		c := idx[y[i]]
		for j, v := range row {
			d := v - means[c][j]
			vars[c][j] += d * d
		}
	}
	var maxVar float64
	for c := 0; c < k; c++ {
		if counts[c] > 0 {
			for j := 0; j < nFeat; j++ {
				vars[c][j] /= float64(counts[c])
				if vars[c][j] > maxVar {
					maxVar = vars[c][j]
				}
			}
		}
	}
	smooth := 1e-9 * maxVar
	if smooth == 0 {
		smooth = 1e-9
	}
	priors := make([]float64, k)
	for c := 0; c < k; c++ {
		priors[c] = float64(counts[c]) / float64(len(y))
		for j := 0; j < nFeat; j++ {
			vars[c][j] += smooth
		}
	}
	g.classes = classes
	g.priors = priors
	g.means = means
	g.vars = vars
	g.nFeat = nFeat
	g.varSmall = smooth
	g.fitted = true
	return nil
}

// Classes returns the sorted class labels learned by Fit.
func (g *GaussianNB) Classes() []int { return g.classes }

// Priors returns the estimated class prior probabilities aligned with Classes.
func (g *GaussianNB) Priors() []float64 { return g.priors }

// JointLogLikelihood returns, for the sample x, the unnormalized log joint
// probabilities log P(class) + log P(x | class) for each class.
func (g *GaussianNB) JointLogLikelihood(x []float64) ([]float64, error) {
	if !g.fitted {
		return nil, ErrNotFitted
	}
	if len(x) != g.nFeat {
		return nil, ErrShape
	}
	out := make([]float64, len(g.classes))
	for c := range g.classes {
		ll := math.Log(g.priors[c])
		for j := 0; j < g.nFeat; j++ {
			v := g.vars[c][j]
			d := x[j] - g.means[c][j]
			ll += -0.5*math.Log(2*math.Pi*v) - d*d/(2*v)
		}
		out[c] = ll
	}
	return out, nil
}

// PredictProba returns the posterior class probabilities for the sample x.
func (g *GaussianNB) PredictProba(x []float64) ([]float64, error) {
	jll, err := g.JointLogLikelihood(x)
	if err != nil {
		return nil, err
	}
	norm := LogSumExp(jll)
	for i := range jll {
		jll[i] = math.Exp(jll[i] - norm)
	}
	return jll, nil
}

// Predict returns the most probable class label for the sample x.
func (g *GaussianNB) Predict(x []float64) (int, error) {
	jll, err := g.JointLogLikelihood(x)
	if err != nil {
		return 0, err
	}
	best := 0
	for i := 1; i < len(jll); i++ {
		if jll[i] > jll[best] {
			best = i
		}
	}
	return g.classes[best], nil
}

// PredictBatch returns predictions for every row of X.
func (g *GaussianNB) PredictBatch(X [][]float64) ([]int, error) {
	out := make([]int, len(X))
	for i, x := range X {
		p, err := g.Predict(x)
		if err != nil {
			return nil, err
		}
		out[i] = p
	}
	return out, nil
}

// ------------------------------------------------------------------
// Multinomial naive Bayes
// ------------------------------------------------------------------

// MultinomialNB is a multinomial naive-Bayes classifier suitable for count
// features such as word frequencies, using additive (Laplace/Lidstone)
// smoothing controlled by Alpha.
type MultinomialNB struct {
	classes    []int
	priors     []float64
	featLogPrb [][]float64 // log P(feature j | class c)
	nFeat      int
	alpha      float64
	fitted     bool
}

// NewMultinomialNB returns an unfitted multinomial naive-Bayes classifier with
// the given additive smoothing parameter alpha (alpha = 1 is Laplace
// smoothing). A non-positive alpha is replaced by 1.
func NewMultinomialNB(alpha float64) *MultinomialNB {
	if alpha <= 0 {
		alpha = 1
	}
	return &MultinomialNB{alpha: alpha}
}

// Fit estimates class priors and smoothed per-feature log probabilities from
// the count matrix X and labels y.
func (m *MultinomialNB) Fit(X [][]float64, y []int) error {
	if len(X) == 0 || len(X) != len(y) {
		return ErrShape
	}
	nFeat := len(X[0])
	if nFeat == 0 {
		return ErrShape
	}
	classes, idx := uniqueSortedClasses(y)
	k := len(classes)
	counts := make([]int, k)
	featCount := make([][]float64, k)
	for c := 0; c < k; c++ {
		featCount[c] = make([]float64, nFeat)
	}
	for i, row := range X {
		if len(row) != nFeat {
			return ErrShape
		}
		c := idx[y[i]]
		counts[c]++
		for j, v := range row {
			featCount[c][j] += v
		}
	}
	logPrb := make([][]float64, k)
	priors := make([]float64, k)
	for c := 0; c < k; c++ {
		priors[c] = float64(counts[c]) / float64(len(y))
		var total float64
		for j := 0; j < nFeat; j++ {
			total += featCount[c][j]
		}
		denom := total + m.alpha*float64(nFeat)
		logPrb[c] = make([]float64, nFeat)
		for j := 0; j < nFeat; j++ {
			logPrb[c][j] = math.Log((featCount[c][j] + m.alpha) / denom)
		}
	}
	m.classes = classes
	m.priors = priors
	m.featLogPrb = logPrb
	m.nFeat = nFeat
	m.fitted = true
	return nil
}

// Classes returns the sorted class labels learned by Fit.
func (m *MultinomialNB) Classes() []int { return m.classes }

// Priors returns the estimated class prior probabilities aligned with Classes.
func (m *MultinomialNB) Priors() []float64 { return m.priors }

// JointLogLikelihood returns the unnormalized log joint probabilities for the
// count vector x under each class.
func (m *MultinomialNB) JointLogLikelihood(x []float64) ([]float64, error) {
	if !m.fitted {
		return nil, ErrNotFitted
	}
	if len(x) != m.nFeat {
		return nil, ErrShape
	}
	out := make([]float64, len(m.classes))
	for c := range m.classes {
		ll := math.Log(m.priors[c])
		for j := 0; j < m.nFeat; j++ {
			ll += x[j] * m.featLogPrb[c][j]
		}
		out[c] = ll
	}
	return out, nil
}

// PredictProba returns the posterior class probabilities for the count vector x.
func (m *MultinomialNB) PredictProba(x []float64) ([]float64, error) {
	jll, err := m.JointLogLikelihood(x)
	if err != nil {
		return nil, err
	}
	norm := LogSumExp(jll)
	for i := range jll {
		jll[i] = math.Exp(jll[i] - norm)
	}
	return jll, nil
}

// Predict returns the most probable class label for the count vector x.
func (m *MultinomialNB) Predict(x []float64) (int, error) {
	jll, err := m.JointLogLikelihood(x)
	if err != nil {
		return 0, err
	}
	best := 0
	for i := 1; i < len(jll); i++ {
		if jll[i] > jll[best] {
			best = i
		}
	}
	return m.classes[best], nil
}

// ------------------------------------------------------------------
// Bernoulli naive Bayes
// ------------------------------------------------------------------

// BernoulliNB is a Bernoulli naive-Bayes classifier for binary feature vectors;
// each feature contributes both when present and, via its complement, when
// absent. Additive smoothing is controlled by Alpha.
type BernoulliNB struct {
	classes   []int
	priors    []float64
	featProb  [][]float64 // P(feature j = 1 | class c)
	nFeat     int
	alpha     float64
	fitted    bool
	threshold float64
}

// NewBernoulliNB returns an unfitted Bernoulli naive-Bayes classifier. Features
// with value greater than threshold are treated as present (1); alpha is the
// additive smoothing parameter (replaced by 1 if non-positive).
func NewBernoulliNB(alpha, threshold float64) *BernoulliNB {
	if alpha <= 0 {
		alpha = 1
	}
	return &BernoulliNB{alpha: alpha, threshold: threshold}
}

// Fit estimates class priors and smoothed feature-presence probabilities.
func (b *BernoulliNB) Fit(X [][]float64, y []int) error {
	if len(X) == 0 || len(X) != len(y) {
		return ErrShape
	}
	nFeat := len(X[0])
	if nFeat == 0 {
		return ErrShape
	}
	classes, idx := uniqueSortedClasses(y)
	k := len(classes)
	counts := make([]int, k)
	present := make([][]float64, k)
	for c := 0; c < k; c++ {
		present[c] = make([]float64, nFeat)
	}
	for i, row := range X {
		if len(row) != nFeat {
			return ErrShape
		}
		c := idx[y[i]]
		counts[c]++
		for j, v := range row {
			if v > b.threshold {
				present[c][j]++
			}
		}
	}
	prob := make([][]float64, k)
	priors := make([]float64, k)
	for c := 0; c < k; c++ {
		priors[c] = float64(counts[c]) / float64(len(y))
		prob[c] = make([]float64, nFeat)
		denom := float64(counts[c]) + 2*b.alpha
		for j := 0; j < nFeat; j++ {
			prob[c][j] = (present[c][j] + b.alpha) / denom
		}
	}
	b.classes = classes
	b.priors = priors
	b.featProb = prob
	b.nFeat = nFeat
	b.fitted = true
	return nil
}

// Classes returns the sorted class labels learned by Fit.
func (b *BernoulliNB) Classes() []int { return b.classes }

// JointLogLikelihood returns the unnormalized log joint probabilities for the
// binary feature vector x under each class.
func (b *BernoulliNB) JointLogLikelihood(x []float64) ([]float64, error) {
	if !b.fitted {
		return nil, ErrNotFitted
	}
	if len(x) != b.nFeat {
		return nil, ErrShape
	}
	out := make([]float64, len(b.classes))
	for c := range b.classes {
		ll := math.Log(b.priors[c])
		for j := 0; j < b.nFeat; j++ {
			p := b.featProb[c][j]
			if x[j] > b.threshold {
				ll += math.Log(p)
			} else {
				ll += math.Log(1 - p)
			}
		}
		out[c] = ll
	}
	return out, nil
}

// PredictProba returns the posterior class probabilities for x.
func (b *BernoulliNB) PredictProba(x []float64) ([]float64, error) {
	jll, err := b.JointLogLikelihood(x)
	if err != nil {
		return nil, err
	}
	norm := LogSumExp(jll)
	for i := range jll {
		jll[i] = math.Exp(jll[i] - norm)
	}
	return jll, nil
}

// Predict returns the most probable class label for x.
func (b *BernoulliNB) Predict(x []float64) (int, error) {
	jll, err := b.JointLogLikelihood(x)
	if err != nil {
		return 0, err
	}
	best := 0
	for i := 1; i < len(jll); i++ {
		if jll[i] > jll[best] {
			best = i
		}
	}
	return b.classes[best], nil
}

// Accuracy returns the fraction of samples in X whose predicted label from the
// given classifier matches the true label in y. The classifier is any type with
// a Predict([]float64) (int, error) method.
func Accuracy(clf interface {
	Predict([]float64) (int, error)
}, X [][]float64, y []int) (float64, error) {
	if len(X) != len(y) || len(X) == 0 {
		return 0, ErrShape
	}
	correct := 0
	for i, x := range X {
		p, err := clf.Predict(x)
		if err != nil {
			return 0, err
		}
		if p == y[i] {
			correct++
		}
	}
	return float64(correct) / float64(len(y)), nil
}
