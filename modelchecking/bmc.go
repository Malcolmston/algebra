package modelchecking

// Unrolling records the breadth-first unrolling of a Kripke structure's
// transition relation from its initial states up to a fixed bound. Layer i holds
// the states first reached at BFS depth i.
type Unrolling struct {
	Bound  int
	Layers []StateSet
}

// Unroll computes the bounded breadth-first unrolling of k from its initial
// states to the given depth bound. Layer 0 is the set of initial states; layer i
// contains the states whose shortest distance from an initial state is exactly
// i, for i up to bound.
func Unroll(k *Kripke, bound int) Unrolling {
	if bound < 0 {
		bound = 0
	}
	seen := NewStateSet(k.n)
	layers := make([]StateSet, 0, bound+1)
	frontier := NewStateSet(k.n)
	for _, s := range k.InitialStates() {
		frontier.Add(s)
		seen.Add(s)
	}
	layers = append(layers, frontier)
	for d := 1; d <= bound; d++ {
		next := NewStateSet(k.n)
		for _, u := range frontier.Elements() {
			for _, v := range k.succ[u] {
				if !seen.Contains(v) {
					seen.Add(v)
					next.Add(v)
				}
			}
		}
		layers = append(layers, next)
		frontier = next
		if next.IsEmpty() {
			break
		}
	}
	return Unrolling{Bound: bound, Layers: layers}
}

// ReachableWithin returns the union of all layers of the unrolling: the states
// reachable within Bound steps of an initial state.
func (u Unrolling) ReachableWithin() StateSet {
	if len(u.Layers) == 0 {
		return NewStateSet(0)
	}
	acc := u.Layers[0].Clone()
	for i := 1; i < len(u.Layers); i++ {
		acc = acc.Union(u.Layers[i])
	}
	return acc
}

// Depth returns the number of populated layers minus one, i.e. the greatest BFS
// depth at which a new state was discovered.
func (u Unrolling) Depth() int {
	d := 0
	for i, l := range u.Layers {
		if !l.IsEmpty() {
			d = i
		}
	}
	return d
}

// bfsKripkePath returns a shortest path (inclusive of endpoints) from any of the
// sources to some state of target within at most bound steps, or nil if none
// exists. When target contains a source the singleton path is returned.
func bfsKripkePath(k *Kripke, sources []int, target StateSet, bound int) []int {
	parent := make([]int, k.n)
	depth := make([]int, k.n)
	for i := range parent {
		parent[i] = -2
	}
	var queue []int
	for _, s := range sources {
		if s < 0 || s >= k.n || parent[s] != -2 {
			continue
		}
		parent[s] = -1
		depth[s] = 0
		if target.Contains(s) {
			return []int{s}
		}
		queue = append(queue, s)
	}
	for len(queue) > 0 {
		u := queue[0]
		queue = queue[1:]
		if depth[u] >= bound {
			continue
		}
		for _, v := range k.succ[u] {
			if parent[v] == -2 {
				parent[v] = u
				depth[v] = depth[u] + 1
				if target.Contains(v) {
					return kReconstruct(parent, v)
				}
				queue = append(queue, v)
			}
		}
	}
	return nil
}

func kReconstruct(parent []int, target int) []int {
	var rev []int
	for x := target; x != -1; x = parent[x] {
		rev = append(rev, x)
	}
	for i, j := 0, len(rev)-1; i < j; i, j = i+1, j-1 {
		rev[i], rev[j] = rev[j], rev[i]
	}
	return rev
}

// BoundedReachable searches for a path of length at most bound from an initial
// state to some state in target. It returns the path and true when one is found.
func BoundedReachable(k *Kripke, target StateSet, bound int) ([]int, bool) {
	p := bfsKripkePath(k, k.InitialStates(), target, bound)
	return p, p != nil
}

// BMCInvariant performs bounded model checking of the invariant "always inv": it
// searches for a path of length at most bound from an initial state to a state
// violating inv. When such a path exists it returns a counterexample and true;
// otherwise it returns (nil, false), meaning no violation exists within the
// bound (which does not by itself prove the invariant on the unbounded system).
func BMCInvariant(k *Kripke, inv StateSet, bound int) (*Counterexample, bool) {
	bad := inv.Complement()
	path := bfsKripkePath(k, k.InitialStates(), bad, bound)
	if path == nil {
		return nil, false
	}
	return &Counterexample{
		Kripke:  k,
		Lasso:   Lasso{Prefix: path, Loop: nil},
		Formula: "AG(inv)",
	}, true
}

// BMCReachInvariant is the bounded analogue of proving AG p: it reports whether
// no violation of p is reachable within bound steps.
func BMCReachInvariant(k *Kripke, p StateSet, bound int) bool {
	_, found := BMCInvariant(k, p, bound)
	return !found
}

// BMCFindLasso searches for an ultimately periodic path (a lasso) reachable from
// an initial state whose total length (prefix plus loop) is at most bound. It
// returns the lasso and true when one is found. Such a lasso is a bounded
// witness to the existence of an infinite run.
func BMCFindLasso(k *Kripke, bound int) (Lasso, bool) {
	return bmcLassoIn(k, FullStateSet(k.n), bound)
}

// BMCExistsGlobally searches within the bound for a lasso all of whose states
// satisfy p, i.e. a bounded witness to EG p from an initial state.
func BMCExistsGlobally(k *Kripke, p StateSet, bound int) (Lasso, bool) {
	return bmcLassoIn(k, p, bound)
}

// bmcLassoIn does a depth-bounded search from the initial states for a path that
// revisits a state, restricted to the set allowed.
func bmcLassoIn(k *Kripke, allowed StateSet, bound int) (Lasso, bool) {
	var path []int
	onPath := map[int]int{} // state -> position in path
	var dfs func(u, depth int) (Lasso, bool)
	dfs = func(u, depth int) (Lasso, bool) {
		if pos, ok := onPath[u]; ok {
			// found a loop back to position pos
			prefix := append([]int(nil), path[:pos]...)
			loop := append([]int(nil), path[pos:]...)
			return Lasso{Prefix: prefix, Loop: loop}, true
		}
		if depth >= bound {
			return Lasso{}, false
		}
		onPath[u] = len(path)
		path = append(path, u)
		for _, v := range k.succ[u] {
			if !allowed.Contains(v) {
				continue
			}
			if l, ok := dfs(v, depth+1); ok {
				return l, true
			}
		}
		path = path[:len(path)-1]
		delete(onPath, u)
		return Lasso{}, false
	}
	for _, s := range k.InitialStates() {
		if !allowed.Contains(s) {
			continue
		}
		path = path[:0]
		onPath = map[int]int{}
		if l, ok := dfs(s, 0); ok {
			return l, true
		}
	}
	return Lasso{}, false
}

// BMCUnrollReachSet returns the exact set of states reachable within bound
// steps, computed from the unrolling. It is the bounded frontier used by safety
// checking.
func BMCUnrollReachSet(k *Kripke, bound int) StateSet {
	return Unroll(k, bound).ReachableWithin()
}
