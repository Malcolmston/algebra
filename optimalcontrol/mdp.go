package optimalcontrol

import "math"

// MDP is a finite, discounted Markov decision process with the standard
// reward-maximization convention. Transitions for each action are stored as
// row-stochastic matrices and rewards as expected immediate rewards per
// state-action pair.
type MDP struct {
	// States is the number of states.
	States int
	// Actions is the number of actions.
	Actions int
	// Trans[a] is the States×States transition matrix under action a; each row
	// must sum to one.
	Trans []*Matrix
	// Reward is a States×Actions matrix of expected immediate rewards.
	Reward *Matrix
	// Gamma is the discount factor in [0, 1).
	Gamma float64
}

// NewMDP constructs an MDP from per-action transition matrices, a reward matrix
// and a discount factor, validating the dimensions.
func NewMDP(trans []*Matrix, reward *Matrix, gamma float64) (*MDP, error) {
	nA := len(trans)
	if nA == 0 {
		return nil, ErrDim
	}
	nS := trans[0].rows
	for _, t := range trans {
		if t.rows != nS || t.cols != nS {
			return nil, ErrDim
		}
	}
	if reward.rows != nS || reward.cols != nA {
		return nil, ErrDim
	}
	return &MDP{States: nS, Actions: nA, Trans: trans, Reward: reward, Gamma: gamma}, nil
}

// QValues returns the States×Actions matrix of action values
// Q(s, a) = R(s, a) + γ Σ_{s'} P_a(s, s') V(s') for a given value function.
func (m *MDP) QValues(v []float64) *Matrix {
	q := Zeros(m.States, m.Actions)
	for a := 0; a < m.Actions; a++ {
		pv := m.Trans[a].MulVec(v)
		for s := 0; s < m.States; s++ {
			q.Set(s, a, m.Reward.At(s, a)+m.Gamma*pv[s])
		}
	}
	return q
}

// BellmanBackup applies one Bellman optimality update, returning the improved
// value function and the greedy policy that attains it.
func (m *MDP) BellmanBackup(v []float64) (newV []float64, policy []int) {
	q := m.QValues(v)
	newV = make([]float64, m.States)
	policy = make([]int, m.States)
	for s := 0; s < m.States; s++ {
		best := math.Inf(-1)
		bestA := 0
		for a := 0; a < m.Actions; a++ {
			if val := q.At(s, a); val > best {
				best = val
				bestA = a
			}
		}
		newV[s] = best
		policy[s] = bestA
	}
	return newV, policy
}

// ValueIterationResult holds the outcome of value iteration.
type ValueIterationResult struct {
	// Value is the (approximately) optimal value function.
	Value []float64
	// Policy is the greedy policy with respect to Value.
	Policy []int
	// Iterations is the number of sweeps performed.
	Iterations int
	// Converged reports whether the max-norm residual fell below the tolerance.
	Converged bool
}

// ValueIteration runs value iteration until the max-norm change between sweeps
// falls below tol or maxIter sweeps have elapsed.
func (m *MDP) ValueIteration(tol float64, maxIter int) *ValueIterationResult {
	v := make([]float64, m.States)
	res := &ValueIterationResult{}
	for iter := 1; iter <= maxIter; iter++ {
		newV, policy := m.BellmanBackup(v)
		var delta float64
		for s := range newV {
			if d := math.Abs(newV[s] - v[s]); d > delta {
				delta = d
			}
		}
		v = newV
		res.Value = v
		res.Policy = policy
		res.Iterations = iter
		if delta < tol {
			res.Converged = true
			break
		}
	}
	return res
}

// GaussSeidelValueIteration runs value iteration using in-place (Gauss–Seidel)
// updates, which typically converges in fewer sweeps than the Jacobi form.
func (m *MDP) GaussSeidelValueIteration(tol float64, maxIter int) *ValueIterationResult {
	v := make([]float64, m.States)
	res := &ValueIterationResult{}
	for iter := 1; iter <= maxIter; iter++ {
		var delta float64
		for s := 0; s < m.States; s++ {
			best := math.Inf(-1)
			for a := 0; a < m.Actions; a++ {
				val := m.Reward.At(s, a)
				row := m.Trans[a]
				var acc float64
				for sp := 0; sp < m.States; sp++ {
					acc += row.At(s, sp) * v[sp]
				}
				val += m.Gamma * acc
				if val > best {
					best = val
				}
			}
			if d := math.Abs(best - v[s]); d > delta {
				delta = d
			}
			v[s] = best
		}
		res.Iterations = iter
		if delta < tol {
			res.Converged = true
			break
		}
	}
	_, policy := m.BellmanBackup(v)
	res.Value = v
	res.Policy = policy
	return res
}

// GreedyPolicy returns the greedy policy with respect to a value function.
func (m *MDP) GreedyPolicy(v []float64) []int {
	_, policy := m.BellmanBackup(v)
	return policy
}

// PolicyMatrices returns the transition matrix P_π and reward vector r_π induced
// by a deterministic policy.
func (m *MDP) PolicyMatrices(policy []int) (*Matrix, []float64) {
	p := Zeros(m.States, m.States)
	r := make([]float64, m.States)
	for s := 0; s < m.States; s++ {
		a := policy[s]
		p.SetRow(s, m.Trans[a].Row(s))
		r[s] = m.Reward.At(s, a)
	}
	return p, r
}

// PolicyEvaluationExact computes the exact value function of a deterministic
// policy by solving the linear system (I − γ P_π) V = r_π.
func (m *MDP) PolicyEvaluationExact(policy []int) ([]float64, error) {
	p, r := m.PolicyMatrices(policy)
	lhs := Eye(m.States).Minus(p.Scale(m.Gamma))
	return Solve(lhs, r)
}

// PolicyEvaluationIterative computes the value function of a policy by iterative
// application of the Bellman expectation operator.
func (m *MDP) PolicyEvaluationIterative(policy []int, tol float64, maxIter int) []float64 {
	p, r := m.PolicyMatrices(policy)
	v := make([]float64, m.States)
	for iter := 0; iter < maxIter; iter++ {
		pv := p.MulVec(v)
		var delta float64
		for s := 0; s < m.States; s++ {
			nv := r[s] + m.Gamma*pv[s]
			if d := math.Abs(nv - v[s]); d > delta {
				delta = d
			}
			v[s] = nv
		}
		if delta < tol {
			break
		}
	}
	return v
}

// PolicyIterationResult holds the outcome of policy iteration.
type PolicyIterationResult struct {
	// Policy is the optimal deterministic policy.
	Policy []int
	// Value is the value function of the optimal policy.
	Value []float64
	// Iterations is the number of policy-improvement steps performed.
	Iterations int
}

// PolicyIteration runs Howard's policy iteration with exact policy evaluation,
// returning the optimal policy and its value function.
func (m *MDP) PolicyIteration(maxIter int) (*PolicyIterationResult, error) {
	policy := make([]int, m.States)
	res := &PolicyIterationResult{}
	for iter := 1; iter <= maxIter; iter++ {
		v, err := m.PolicyEvaluationExact(policy)
		if err != nil {
			return nil, err
		}
		newPolicy := m.GreedyPolicy(v)
		stable := true
		for s := range policy {
			if newPolicy[s] != policy[s] {
				stable = false
				break
			}
		}
		policy = newPolicy
		res.Policy = policy
		res.Value = v
		res.Iterations = iter
		if stable {
			break
		}
	}
	return res, nil
}

// ModifiedPolicyIteration runs modified policy iteration, evaluating each policy
// with a fixed number k of Bellman expectation sweeps rather than an exact
// solve.
func (m *MDP) ModifiedPolicyIteration(k, maxIter int, tol float64) (*PolicyIterationResult, error) {
	policy := make([]int, m.States)
	v := make([]float64, m.States)
	res := &PolicyIterationResult{}
	for iter := 1; iter <= maxIter; iter++ {
		p, r := m.PolicyMatrices(policy)
		for j := 0; j < k; j++ {
			pv := p.MulVec(v)
			for s := 0; s < m.States; s++ {
				v[s] = r[s] + m.Gamma*pv[s]
			}
		}
		newV, newPolicy := m.BellmanBackup(v)
		var delta float64
		for s := range newV {
			if d := math.Abs(newV[s] - v[s]); d > delta {
				delta = d
			}
		}
		copy(v, newV)
		policy = newPolicy
		res.Policy = policy
		res.Value = v
		res.Iterations = iter
		if delta < tol {
			break
		}
	}
	return res, nil
}
