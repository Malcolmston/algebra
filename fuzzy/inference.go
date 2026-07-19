package fuzzy

import (
	"errors"
	"math"
)

// ErrUnknownVariable is returned when a rule references an input or output
// variable that the system does not define.
var ErrUnknownVariable = errors.New("fuzzy: unknown variable")

// ErrUnknownTerm is returned when a rule references a linguistic term that a
// variable does not define.
var ErrUnknownTerm = errors.New("fuzzy: unknown linguistic term")

// ErrNoInput is returned when a crisp value for a referenced input variable is
// missing.
var ErrNoInput = errors.New("fuzzy: missing crisp input for variable")

// ErrNoRules is returned when an inference system is asked to reason with no
// rules defined.
var ErrNoRules = errors.New("fuzzy: no rules defined")

// Variable is a linguistic variable: a named universe [Min, Max] together with
// a collection of named linguistic terms, each described by a membership
// function.
type Variable struct {
	Name  string
	Min   float64
	Max   float64
	Terms map[string]MF
}

// NewVariable creates a linguistic variable named name over the universe
// [min, max] with no terms.
func NewVariable(name string, min, max float64) *Variable {
	return &Variable{Name: name, Min: min, Max: max, Terms: map[string]MF{}}
}

// AddTerm registers the linguistic term named term with membership function mf
// and returns the variable for chaining.
func (v *Variable) AddTerm(term string, mf MF) *Variable {
	v.Terms[term] = mf
	return v
}

// Fuzzify returns the membership grade of the crisp value x in the term named
// term. It returns ErrUnknownTerm when the term is not defined.
func (v *Variable) Fuzzify(term string, x float64) (float64, error) {
	mf, ok := v.Terms[term]
	if !ok {
		return 0, ErrUnknownTerm
	}
	return clamp01(mf(x)), nil
}

// RuleOp selects how the clauses of a rule antecedent are combined.
type RuleOp int

const (
	// OpAnd combines antecedent clauses with the system t-norm (conjunction).
	OpAnd RuleOp = iota
	// OpOr combines antecedent clauses with the system t-conorm (disjunction).
	OpOr
)

// Clause is one atomic antecedent proposition "Var is [not] Term".
type Clause struct {
	Var  string
	Term string
	Not  bool
}

// If is a convenience constructor for a positive antecedent clause "v is t".
func If(v, t string) Clause { return Clause{Var: v, Term: t} }

// IfNot is a convenience constructor for a negated antecedent clause
// "v is not t".
func IfNot(v, t string) Clause { return Clause{Var: v, Term: t, Not: true} }

// Implication selects how a rule's firing strength shapes its consequent set.
type Implication int

const (
	// ImplicationMin clips the consequent membership at the firing strength
	// (Mamdani min implication).
	ImplicationMin Implication = iota
	// ImplicationProduct scales the consequent membership by the firing
	// strength (Larsen product implication).
	ImplicationProduct
)

// MamdaniRule is one rule of a Mamdani system: an antecedent (clauses combined
// by Op), a consequent term OutTerm of the output variable and a rule weight.
type MamdaniRule struct {
	Antecedents []Clause
	Op          RuleOp
	OutTerm     string
	Weight      float64
}

// MamdaniSystem is a Mamdani (max-min) fuzzy inference system with several
// input variables and a single output variable.
type MamdaniSystem struct {
	Inputs     map[string]*Variable
	Output     *Variable
	Rules      []MamdaniRule
	And        TNorm
	Or         TConorm
	Not        func(float64) float64
	Impl       Implication
	Aggregate  TConorm
	Method     Defuzz
	Resolution int
}

// NewMamdani creates a Mamdani system with the given output variable and the
// standard operator set (min t-norm, max t-conorm, standard complement, min
// implication, max aggregation, centroid defuzzification and a 101 point output
// sampling grid).
func NewMamdani(output *Variable) *MamdaniSystem {
	return &MamdaniSystem{
		Inputs:     map[string]*Variable{},
		Output:     output,
		And:        TNormMin,
		Or:         TConormMax,
		Not:        ComplementStandard,
		Impl:       ImplicationMin,
		Aggregate:  TConormMax,
		Method:     DefuzzCentroid,
		Resolution: 101,
	}
}

// AddInput registers an input variable and returns the system for chaining.
func (m *MamdaniSystem) AddInput(v *Variable) *MamdaniSystem {
	m.Inputs[v.Name] = v
	return m
}

// AddRule appends a Mamdani rule and returns the system for chaining. A weight
// of 0 is treated as the default weight 1.
func (m *MamdaniSystem) AddRule(op RuleOp, outTerm string, weight float64, clauses ...Clause) *MamdaniSystem {
	if weight == 0 {
		weight = 1
	}
	m.Rules = append(m.Rules, MamdaniRule{Antecedents: clauses, Op: op, OutTerm: outTerm, Weight: weight})
	return m
}

// clauseGrade evaluates a single antecedent clause against the crisp inputs.
func (m *MamdaniSystem) clauseGrade(c Clause, inputs map[string]float64) (float64, error) {
	v, ok := m.Inputs[c.Var]
	if !ok {
		return 0, ErrUnknownVariable
	}
	x, ok := inputs[c.Var]
	if !ok {
		return 0, ErrNoInput
	}
	g, err := v.Fuzzify(c.Term, x)
	if err != nil {
		return 0, err
	}
	if c.Not {
		g = clamp01(m.Not(g))
	}
	return g, nil
}

// FiringStrength returns the degree to which the rule's antecedent is satisfied
// by the crisp inputs, combining its clauses with the system's And or Or
// operator and multiplying by the rule weight.
func (m *MamdaniSystem) FiringStrength(rule MamdaniRule, inputs map[string]float64) (float64, error) {
	if len(rule.Antecedents) == 0 {
		return 0, nil
	}
	var acc float64
	for i, c := range rule.Antecedents {
		g, err := m.clauseGrade(c, inputs)
		if err != nil {
			return 0, err
		}
		if i == 0 {
			acc = g
			continue
		}
		if rule.Op == OpAnd {
			acc = m.And(acc, g)
		} else {
			acc = m.Or(acc, g)
		}
	}
	return clamp01(acc * rule.Weight), nil
}

// AggregateOutput returns the aggregated output fuzzy set produced by firing all
// rules against the crisp inputs, sampled on the output universe at Resolution
// points. Each rule clips or scales its consequent term by its firing strength
// and the results are combined with the aggregation t-conorm.
func (m *MamdaniSystem) AggregateOutput(inputs map[string]float64) (Set, error) {
	if len(m.Rules) == 0 {
		return Set{}, ErrNoRules
	}
	xs := Linspace(m.Output.Min, m.Output.Max, m.Resolution)
	agg := make([]float64, len(xs))
	for _, rule := range m.Rules {
		w, err := m.FiringStrength(rule, inputs)
		if err != nil {
			return Set{}, err
		}
		if w == 0 {
			continue
		}
		mf, ok := m.Output.Terms[rule.OutTerm]
		if !ok {
			return Set{}, ErrUnknownTerm
		}
		for i, x := range xs {
			cons := mf(x)
			var shaped float64
			if m.Impl == ImplicationProduct {
				shaped = cons * w
			} else {
				shaped = math.Min(cons, w)
			}
			agg[i] = clamp01(m.Aggregate(agg[i], shaped))
		}
	}
	return NewSet(xs, agg)
}

// Infer runs the full Mamdani pipeline (fuzzify, apply rules, aggregate,
// defuzzify) and returns the crisp output. It returns ErrNoArea when no rule
// fires.
func (m *MamdaniSystem) Infer(inputs map[string]float64) (float64, error) {
	agg, err := m.AggregateOutput(inputs)
	if err != nil {
		return 0, err
	}
	return agg.Defuzzify(m.Method)
}

// SugenoConsequent is a first-order Takagi-Sugeno-Kang consequent
// z = Const + sum(Coeffs[var] * input[var]). A constant (zeroth order)
// consequent uses an empty Coeffs map.
type SugenoConsequent struct {
	Const  float64
	Coeffs map[string]float64
}

// Eval computes the crisp consequent value for the given crisp inputs. Missing
// inputs referenced by a coefficient are treated as ErrNoInput.
func (c SugenoConsequent) Eval(inputs map[string]float64) (float64, error) {
	z := c.Const
	for name, k := range c.Coeffs {
		x, ok := inputs[name]
		if !ok {
			return 0, ErrNoInput
		}
		z += k * x
	}
	return z, nil
}

// SugenoRule is one rule of a Sugeno system: an antecedent (clauses combined by
// Op), a crisp/linear consequent and a rule weight.
type SugenoRule struct {
	Antecedents []Clause
	Op          RuleOp
	Consequent  SugenoConsequent
	Weight      float64
}

// SugenoSystem is a Takagi-Sugeno-Kang fuzzy inference system whose crisp
// output is the firing-strength weighted average of the rule consequents.
type SugenoSystem struct {
	Inputs map[string]*Variable
	Rules  []SugenoRule
	And    TNorm
	Or     TConorm
	Not    func(float64) float64
}

// NewSugeno creates a Sugeno system with the standard operator set (min
// t-norm, max t-conorm and standard complement).
func NewSugeno() *SugenoSystem {
	return &SugenoSystem{
		Inputs: map[string]*Variable{},
		And:    TNormMin,
		Or:     TConormMax,
		Not:    ComplementStandard,
	}
}

// AddInput registers an input variable and returns the system for chaining.
func (s *SugenoSystem) AddInput(v *Variable) *SugenoSystem {
	s.Inputs[v.Name] = v
	return s
}

// AddConstRule appends a zeroth-order Sugeno rule whose consequent is the
// constant z and returns the system for chaining.
func (s *SugenoSystem) AddConstRule(op RuleOp, z, weight float64, clauses ...Clause) *SugenoSystem {
	if weight == 0 {
		weight = 1
	}
	s.Rules = append(s.Rules, SugenoRule{
		Antecedents: clauses,
		Op:          op,
		Consequent:  SugenoConsequent{Const: z, Coeffs: map[string]float64{}},
		Weight:      weight,
	})
	return s
}

// AddRule appends a first-order Sugeno rule with an explicit consequent and
// returns the system for chaining.
func (s *SugenoSystem) AddRule(op RuleOp, cons SugenoConsequent, weight float64, clauses ...Clause) *SugenoSystem {
	if weight == 0 {
		weight = 1
	}
	s.Rules = append(s.Rules, SugenoRule{Antecedents: clauses, Op: op, Consequent: cons, Weight: weight})
	return s
}

// clauseGrade evaluates a single antecedent clause against the crisp inputs.
func (s *SugenoSystem) clauseGrade(c Clause, inputs map[string]float64) (float64, error) {
	v, ok := s.Inputs[c.Var]
	if !ok {
		return 0, ErrUnknownVariable
	}
	x, ok := inputs[c.Var]
	if !ok {
		return 0, ErrNoInput
	}
	g, err := v.Fuzzify(c.Term, x)
	if err != nil {
		return 0, err
	}
	if c.Not {
		g = clamp01(s.Not(g))
	}
	return g, nil
}

// FiringStrength returns the degree to which the rule's antecedent is satisfied
// by the crisp inputs, combining its clauses with the system's And or Or
// operator and multiplying by the rule weight.
func (s *SugenoSystem) FiringStrength(rule SugenoRule, inputs map[string]float64) (float64, error) {
	if len(rule.Antecedents) == 0 {
		return 0, nil
	}
	var acc float64
	for i, c := range rule.Antecedents {
		g, err := s.clauseGrade(c, inputs)
		if err != nil {
			return 0, err
		}
		if i == 0 {
			acc = g
			continue
		}
		if rule.Op == OpAnd {
			acc = s.And(acc, g)
		} else {
			acc = s.Or(acc, g)
		}
	}
	return clamp01(acc * rule.Weight), nil
}

// Infer runs the Sugeno pipeline and returns the crisp output, the
// firing-strength weighted average of the rule consequents,
// sum(w_i * z_i) / sum(w_i). It returns ErrNoArea when no rule fires.
func (s *SugenoSystem) Infer(inputs map[string]float64) (float64, error) {
	if len(s.Rules) == 0 {
		return 0, ErrNoRules
	}
	num := 0.0
	den := 0.0
	for _, rule := range s.Rules {
		w, err := s.FiringStrength(rule, inputs)
		if err != nil {
			return 0, err
		}
		if w == 0 {
			continue
		}
		z, err := rule.Consequent.Eval(inputs)
		if err != nil {
			return 0, err
		}
		num += w * z
		den += w
	}
	if den == 0 {
		return 0, ErrNoArea
	}
	return num / den, nil
}
