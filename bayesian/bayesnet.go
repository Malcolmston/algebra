package bayesian

import (
	"errors"
	"sort"
)

// ErrNetwork is returned for malformed Bayesian-network definitions and
// queries.
var ErrNetwork = errors.New("bayesian: invalid network operation")

// BayesianNetwork is a discrete directed graphical model. Each variable has a
// finite cardinality and a conditional probability table (CPT) given its
// parents. The joint distribution factorizes as the product of the CPTs.
type BayesianNetwork struct {
	names   []string
	card    map[string]int
	parents map[string][]string
	cpts    map[string]Factor
}

// NewBayesianNetwork returns an empty Bayesian network.
func NewBayesianNetwork() *BayesianNetwork {
	return &BayesianNetwork{
		card:    map[string]int{},
		parents: map[string][]string{},
		cpts:    map[string]Factor{},
	}
}

// AddVariable declares a categorical variable with the given cardinality. It
// returns ErrNetwork if the variable already exists or the cardinality is
// non-positive.
func (bn *BayesianNetwork) AddVariable(name string, card int) error {
	if card <= 0 {
		return ErrNetwork
	}
	if _, ok := bn.card[name]; ok {
		return ErrNetwork
	}
	bn.names = append(bn.names, name)
	bn.card[name] = card
	return nil
}

// Variables returns the network's variable names in declaration order.
func (bn *BayesianNetwork) Variables() []string {
	return append([]string(nil), bn.names...)
}

// Cardinality returns the cardinality of the named variable, or 0 if unknown.
func (bn *BayesianNetwork) Cardinality(name string) int { return bn.card[name] }

// SetCPT installs the conditional probability table P(name | parents). The
// table is laid out in row-major order with the parents most significant (in
// the order given) and the child variable least significant, so its length must
// equal card(name)·Πcard(parents). Each parent configuration must define a
// proper distribution over the child. It returns ErrNetwork on any shape or
// declaration error.
func (bn *BayesianNetwork) SetCPT(name string, parents []string, table []float64) error {
	c, ok := bn.card[name]
	if !ok {
		return ErrNetwork
	}
	vars := make([]string, 0, len(parents)+1)
	card := make([]int, 0, len(parents)+1)
	for _, p := range parents {
		pc, ok := bn.card[p]
		if !ok {
			return ErrNetwork
		}
		vars = append(vars, p)
		card = append(card, pc)
	}
	vars = append(vars, name)
	card = append(card, c)
	f, err := NewFactor(vars, card, table)
	if err != nil {
		return err
	}
	bn.parents[name] = append([]string(nil), parents...)
	bn.cpts[name] = f
	return nil
}

// Parents returns the parent variables of name in the order they were declared.
func (bn *BayesianNetwork) Parents(name string) []string {
	return append([]string(nil), bn.parents[name]...)
}

// CPT returns the conditional probability factor for name and reports whether
// it has been set.
func (bn *BayesianNetwork) CPT(name string) (Factor, bool) {
	f, ok := bn.cpts[name]
	return f, ok
}

// JointProbability returns P(assignment) for a full assignment of every
// variable, the product of each CPT evaluated at the assignment. It returns 0
// if any variable is unassigned or out of range.
func (bn *BayesianNetwork) JointProbability(assign map[string]int) float64 {
	p := 1.0
	for _, name := range bn.names {
		f, ok := bn.cpts[name]
		if !ok {
			return 0
		}
		p *= f.Get(assign)
	}
	return p
}

// eliminate performs variable elimination and returns the (unnormalized) factor
// over the query variables given the evidence.
func (bn *BayesianNetwork) eliminate(query []string, evidence map[string]int) (Factor, error) {
	for _, q := range query {
		if _, ok := bn.card[q]; !ok {
			return Factor{}, ErrNetwork
		}
	}
	// Gather CPT factors reduced by evidence.
	factors := make([]Factor, 0, len(bn.names))
	for _, name := range bn.names {
		f, ok := bn.cpts[name]
		if !ok {
			return Factor{}, ErrNetwork
		}
		factors = append(factors, f.Reduce(evidence))
	}
	// Determine hidden variables to eliminate (not queried, not evidence).
	inQuery := map[string]bool{}
	for _, q := range query {
		inQuery[q] = true
	}
	hidden := make([]string, 0)
	for _, name := range bn.names {
		if inQuery[name] {
			continue
		}
		if _, ok := evidence[name]; ok {
			continue
		}
		hidden = append(hidden, name)
	}
	sort.Strings(hidden)
	// Eliminate each hidden variable.
	for _, h := range hidden {
		var involved []Factor
		var rest []Factor
		for _, f := range factors {
			if f.varIndex(h) != -1 {
				involved = append(involved, f)
			} else {
				rest = append(rest, f)
			}
		}
		if len(involved) == 0 {
			factors = rest
			continue
		}
		prod := involved[0]
		for i := 1; i < len(involved); i++ {
			prod = prod.Multiply(involved[i])
		}
		prod = prod.Marginalize(h)
		factors = append(rest, prod)
	}
	// Multiply remaining factors together.
	if len(factors) == 0 {
		return Factor{}, ErrNetwork
	}
	result := factors[0]
	for i := 1; i < len(factors); i++ {
		result = result.Multiply(factors[i])
	}
	return result, nil
}

// Marginal returns the marginal distribution over the query variables as a
// normalized Factor.
func (bn *BayesianNetwork) Marginal(query ...string) (Factor, error) {
	f, err := bn.eliminate(query, map[string]int{})
	if err != nil {
		return Factor{}, err
	}
	return f.Normalize(), nil
}

// Conditional returns the posterior distribution over the query variables given
// the evidence, as a normalized Factor.
func (bn *BayesianNetwork) Conditional(query []string, evidence map[string]int) (Factor, error) {
	f, err := bn.eliminate(query, evidence)
	if err != nil {
		return Factor{}, err
	}
	return f.Normalize(), nil
}

// MarginalProbability returns P(variable = value), a scalar marginal for a
// single variable and category.
func (bn *BayesianNetwork) MarginalProbability(variable string, value int) (float64, error) {
	f, err := bn.Marginal(variable)
	if err != nil {
		return 0, err
	}
	return f.Get(map[string]int{variable: value}), nil
}

// ConditionalProbability returns P(variable = value | evidence), a scalar
// posterior for a single variable and category.
func (bn *BayesianNetwork) ConditionalProbability(variable string, value int, evidence map[string]int) (float64, error) {
	f, err := bn.Conditional([]string{variable}, evidence)
	if err != nil {
		return 0, err
	}
	return f.Get(map[string]int{variable: value}), nil
}

// EvidenceProbability returns P(evidence), the probability of the observed
// assignment obtained by summing the joint distribution consistent with it.
func (bn *BayesianNetwork) EvidenceProbability(evidence map[string]int) (float64, error) {
	factors := make([]Factor, 0, len(bn.names))
	for _, name := range bn.names {
		f, ok := bn.cpts[name]
		if !ok {
			return 0, ErrNetwork
		}
		factors = append(factors, f.Reduce(evidence))
	}
	// Multiply all reduced factors and sum out everything.
	result := factors[0]
	for i := 1; i < len(factors); i++ {
		result = result.Multiply(factors[i])
	}
	// Sum out remaining variables.
	for len(result.Vars) > 0 {
		result = result.Marginalize(result.Vars[0])
	}
	if len(result.Table) == 0 {
		return 0, nil
	}
	return result.Table[0], nil
}
