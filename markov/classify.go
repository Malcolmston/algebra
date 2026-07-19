package markov

import "sort"

// AdjacencyMatrix returns the boolean adjacency matrix of the chain's transition
// graph: entry [i][j] is true iff P[i][j] > tol (state j is reachable from i in
// exactly one step).
func (c *MarkovChain) AdjacencyMatrix(tol float64) [][]bool {
	adj := make([][]bool, c.n)
	for i := 0; i < c.n; i++ {
		adj[i] = make([]bool, c.n)
		for j := 0; j < c.n; j++ {
			adj[i][j] = c.p[i][j] > tol
		}
	}
	return adj
}

// Accessible reports whether state j is accessible from state i, i.e. there is
// a path of non-zero-probability transitions from i to j (including i=j at
// distance zero).
func (c *MarkovChain) Accessible(i, j int) bool {
	if i < 0 || i >= c.n || j < 0 || j >= c.n {
		return false
	}
	seen := c.reachSet(i)
	return seen[j]
}

// reachSet returns the set of states reachable from i (including i) via a BFS
// over positive-probability edges.
func (c *MarkovChain) reachSet(i int) []bool {
	seen := make([]bool, c.n)
	seen[i] = true
	stack := []int{i}
	for len(stack) > 0 {
		u := stack[len(stack)-1]
		stack = stack[:len(stack)-1]
		for v := 0; v < c.n; v++ {
			if c.p[u][v] > 0 && !seen[v] {
				seen[v] = true
				stack = append(stack, v)
			}
		}
	}
	return seen
}

// ReachableFrom returns the sorted list of states reachable from state i
// (including i itself).
func (c *MarkovChain) ReachableFrom(i int) []int {
	if i < 0 || i >= c.n {
		return nil
	}
	seen := c.reachSet(i)
	var out []int
	for j, ok := range seen {
		if ok {
			out = append(out, j)
		}
	}
	return out
}

// Communicates reports whether states i and j communicate, i.e. each is
// accessible from the other. A state always communicates with itself.
func (c *MarkovChain) Communicates(i, j int) bool {
	return c.Accessible(i, j) && c.Accessible(j, i)
}

// ReachabilityMatrix returns the boolean matrix whose [i][j] entry is true iff
// state j is accessible from state i (transitive closure of the one-step graph,
// including the diagonal).
func (c *MarkovChain) ReachabilityMatrix() [][]bool {
	r := make([][]bool, c.n)
	for i := 0; i < c.n; i++ {
		r[i] = c.reachSet(i)
	}
	return r
}

// CommunicatingClasses returns the communicating classes (strongly connected
// components of the transition graph) as a slice of sorted state slices. The
// classes themselves are ordered by their smallest member.
func (c *MarkovChain) CommunicatingClasses() [][]int {
	comps := c.tarjanSCC()
	for _, comp := range comps {
		sort.Ints(comp)
	}
	sort.Slice(comps, func(a, b int) bool { return comps[a][0] < comps[b][0] })
	return comps
}

// tarjanSCC computes the strongly connected components of the positive-edge
// transition graph using Tarjan's algorithm.
func (c *MarkovChain) tarjanSCC() [][]int {
	n := c.n
	index := make([]int, n)
	low := make([]int, n)
	onStack := make([]bool, n)
	for i := range index {
		index[i] = -1
	}
	var stack []int
	var comps [][]int
	counter := 0

	var strongconnect func(v int)
	strongconnect = func(v int) {
		index[v] = counter
		low[v] = counter
		counter++
		stack = append(stack, v)
		onStack[v] = true
		for w := 0; w < n; w++ {
			if c.p[v][w] <= 0 {
				continue
			}
			if index[w] == -1 {
				strongconnect(w)
				if low[w] < low[v] {
					low[v] = low[w]
				}
			} else if onStack[w] {
				if index[w] < low[v] {
					low[v] = index[w]
				}
			}
		}
		if low[v] == index[v] {
			var comp []int
			for {
				w := stack[len(stack)-1]
				stack = stack[:len(stack)-1]
				onStack[w] = false
				comp = append(comp, w)
				if w == v {
					break
				}
			}
			comps = append(comps, comp)
		}
	}

	for v := 0; v < n; v++ {
		if index[v] == -1 {
			strongconnect(v)
		}
	}
	return comps
}

// classMembership returns, for each state, the index of its communicating class
// in the slice returned by CommunicatingClasses.
func (c *MarkovChain) classMembership(classes [][]int) []int {
	member := make([]int, c.n)
	for ci, comp := range classes {
		for _, s := range comp {
			member[s] = ci
		}
	}
	return member
}

// IsClosedClass reports whether the given set of states is closed: no
// transition leaves the set. A closed communicating class is recurrent.
func (c *MarkovChain) IsClosedClass(states []int) bool {
	in := make([]bool, c.n)
	for _, s := range states {
		if s < 0 || s >= c.n {
			return false
		}
		in[s] = true
	}
	for _, s := range states {
		for j := 0; j < c.n; j++ {
			if c.p[s][j] > 0 && !in[j] {
				return false
			}
		}
	}
	return true
}

// RecurrentStates returns the sorted list of recurrent states: those belonging
// to a closed communicating class.
func (c *MarkovChain) RecurrentStates() []int {
	classes := c.CommunicatingClasses()
	var out []int
	for _, comp := range classes {
		if c.IsClosedClass(comp) {
			out = append(out, comp...)
		}
	}
	sort.Ints(out)
	return out
}

// TransientStates returns the sorted list of transient states: those not in any
// closed communicating class.
func (c *MarkovChain) TransientStates() []int {
	rec := make(map[int]bool)
	for _, s := range c.RecurrentStates() {
		rec[s] = true
	}
	var out []int
	for s := 0; s < c.n; s++ {
		if !rec[s] {
			out = append(out, s)
		}
	}
	return out
}

// AbsorbingStates returns the sorted list of absorbing states (states i with
// P[i][i] = 1).
func (c *MarkovChain) AbsorbingStates() []int {
	var out []int
	for i := 0; i < c.n; i++ {
		if c.p[i][i] == 1 {
			out = append(out, i)
		}
	}
	return out
}

// IsAbsorbingState reports whether state i is absorbing (P[i][i] = 1).
func (c *MarkovChain) IsAbsorbingState(i int) bool {
	return i >= 0 && i < c.n && c.p[i][i] == 1
}

// IsAbsorbing reports whether the chain is an absorbing chain: it has at least
// one absorbing state and from every state some absorbing state is accessible.
func (c *MarkovChain) IsAbsorbing() bool {
	abs := c.AbsorbingStates()
	if len(abs) == 0 {
		return false
	}
	absSet := make(map[int]bool)
	for _, a := range abs {
		absSet[a] = true
	}
	for i := 0; i < c.n; i++ {
		reach := c.reachSet(i)
		found := false
		for a := range absSet {
			if reach[a] {
				found = true
				break
			}
		}
		if !found {
			return false
		}
	}
	return true
}

// IsIrreducible reports whether the chain is irreducible: all states form a
// single communicating class.
func (c *MarkovChain) IsIrreducible() bool {
	if c.n <= 1 {
		return c.n == 1
	}
	// State 0 must reach every state and every state must reach state 0.
	fromZero := c.reachSet(0)
	for j := 0; j < c.n; j++ {
		if !fromZero[j] {
			return false
		}
	}
	// Build reverse reachability: does every state reach 0?
	rt := c.Reverse01()
	toZero := rt.reachSet(0)
	for j := 0; j < c.n; j++ {
		if !toZero[j] {
			return false
		}
	}
	return true
}

// Reverse01 returns a chain whose transition graph is the edge-reversal of this
// chain's positive-probability graph, with rows normalized to be stochastic.
// It is used internally for reverse-reachability queries; the numeric
// probabilities are not meaningful beyond their support.
func (c *MarkovChain) Reverse01() *MarkovChain {
	q := make([][]float64, c.n)
	for i := range q {
		q[i] = make([]float64, c.n)
	}
	for i := 0; i < c.n; i++ {
		for j := 0; j < c.n; j++ {
			if c.p[i][j] > 0 {
				q[j][i] = 1
			}
		}
	}
	for i := range q {
		var s float64
		for _, x := range q[i] {
			s += x
		}
		if s == 0 {
			q[i][i] = 1
		} else {
			for j := range q[i] {
				q[i][j] /= s
			}
		}
	}
	return &MarkovChain{p: q, n: c.n}
}

// Period returns the period of state i: the greatest common divisor of the
// lengths of all closed paths (returns) from i to itself. An absorbing or
// self-looping state has period 1. A state that never returns to itself is
// assigned period 0 by convention. The period is a class property for
// communicating states.
func (c *MarkovChain) Period(i int) int {
	if i < 0 || i >= c.n {
		return 0
	}
	// BFS assigning levels; gcd of (level[u]+1 - level[v]) over edges u->v that
	// stay within the set reachable-and-co-reachable-from-i (its class). For
	// simplicity we compute over the whole reachable subgraph restricted to
	// states that can return to i.
	reach := c.reachSet(i)
	rt := c.Reverse01()
	back := rt.reachSet(i)
	inClass := make([]bool, c.n)
	for s := 0; s < c.n; s++ {
		inClass[s] = reach[s] && back[s]
	}
	level := make([]int, c.n)
	for k := range level {
		level[k] = -1
	}
	level[i] = 0
	queue := []int{i}
	g := 0
	for len(queue) > 0 {
		u := queue[0]
		queue = queue[1:]
		for v := 0; v < c.n; v++ {
			if c.p[u][v] <= 0 || !inClass[v] {
				continue
			}
			if level[v] == -1 {
				level[v] = level[u] + 1
				queue = append(queue, v)
			} else {
				diff := level[u] + 1 - level[v]
				g = gcdInt(g, absInt(diff))
			}
		}
	}
	if g == 0 {
		// No cycle found through i.
		return 0
	}
	return g
}

// IsAperiodic reports whether the chain is aperiodic, meaning every state has
// period 1. States with no return path (period 0) are ignored.
func (c *MarkovChain) IsAperiodic() bool {
	for i := 0; i < c.n; i++ {
		if p := c.Period(i); p > 1 {
			return false
		}
	}
	return true
}

// IsRegular reports whether the chain is regular: some power P^k has all
// strictly positive entries. Equivalently the chain is irreducible and
// aperiodic. This is checked directly by examining powers up to n^2+1.
func (c *MarkovChain) IsRegular() bool {
	if c.n == 0 {
		return false
	}
	limit := c.n*c.n + 1
	m := CopyMatrix(c.p)
	for k := 1; k <= limit; k++ {
		allPos := true
		for i := 0; i < c.n && allPos; i++ {
			for j := 0; j < c.n; j++ {
				if m[i][j] <= 0 {
					allPos = false
					break
				}
			}
		}
		if allPos {
			return true
		}
		m = MatMul(m, c.p)
	}
	return false
}

// IsErgodic reports whether the chain is ergodic (irreducible and aperiodic),
// which for a finite chain is equivalent to being regular.
func (c *MarkovChain) IsErgodic() bool {
	return c.IsIrreducible() && c.IsAperiodic()
}

// ClassifyState returns a short label for state i: "absorbing", "recurrent", or
// "transient".
func (c *MarkovChain) ClassifyState(i int) string {
	if c.IsAbsorbingState(i) {
		return "absorbing"
	}
	for _, r := range c.RecurrentStates() {
		if r == i {
			return "recurrent"
		}
	}
	return "transient"
}

func gcdInt(a, b int) int {
	for b != 0 {
		a, b = b, a%b
	}
	if a < 0 {
		return -a
	}
	return a
}

func absInt(a int) int {
	if a < 0 {
		return -a
	}
	return a
}
