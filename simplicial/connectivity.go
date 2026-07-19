package simplicial

import "sort"

// unionFind is a disjoint-set forest over an arbitrary set of integer labels.
type unionFind struct {
	parent map[int]int
	rank   map[int]int
}

func newUnionFind() *unionFind {
	return &unionFind{parent: map[int]int{}, rank: map[int]int{}}
}

func (u *unionFind) add(x int) {
	if _, ok := u.parent[x]; !ok {
		u.parent[x] = x
		u.rank[x] = 0
	}
}

func (u *unionFind) find(x int) int {
	for u.parent[x] != x {
		u.parent[x] = u.parent[u.parent[x]]
		x = u.parent[x]
	}
	return x
}

func (u *unionFind) union(a, b int) {
	ra, rb := u.find(a), u.find(b)
	if ra == rb {
		return
	}
	if u.rank[ra] < u.rank[rb] {
		ra, rb = rb, ra
	}
	u.parent[rb] = ra
	if u.rank[ra] == u.rank[rb] {
		u.rank[ra]++
	}
}

// NumConnectedComponents returns the number of connected components of the
// 1-skeleton of the complex — equivalently the rank b_0 of degree-0 homology.
// The empty complex has zero components.
func (c *Complex) NumConnectedComponents() int {
	uf := newUnionFind()
	for _, s := range c.simplices {
		if s.Dim() == 0 {
			uf.add(s.verts[0])
		}
	}
	for _, s := range c.simplices {
		if s.Dim() == 1 {
			uf.union(s.verts[0], s.verts[1])
		}
	}
	roots := map[int]struct{}{}
	for v := range uf.parent {
		roots[uf.find(v)] = struct{}{}
	}
	return len(roots)
}

// IsConnected reports whether the complex has at most one connected component.
// The empty complex is considered connected.
func (c *Complex) IsConnected() bool { return c.NumConnectedComponents() <= 1 }

// ConnectedComponents partitions the vertices into connected components and
// returns, for each component, its sorted vertex list. Components are ordered
// by their smallest vertex.
func (c *Complex) ConnectedComponents() [][]int {
	uf := newUnionFind()
	for _, s := range c.simplices {
		if s.Dim() == 0 {
			uf.add(s.verts[0])
		}
	}
	for _, s := range c.simplices {
		if s.Dim() == 1 {
			uf.union(s.verts[0], s.verts[1])
		}
	}
	groups := map[int][]int{}
	for v := range uf.parent {
		r := uf.find(v)
		groups[r] = append(groups[r], v)
	}
	out := make([][]int, 0, len(groups))
	for _, g := range groups {
		sort.Ints(g)
		out = append(out, g)
	}
	sort.Slice(out, func(i, j int) bool { return out[i][0] < out[j][0] })
	return out
}

// ComponentOf returns the sorted vertex list of the connected component
// containing vertex v, or nil if v is not a vertex of the complex.
func (c *Complex) ComponentOf(v int) []int {
	if !c.HasVertex(v) {
		return nil
	}
	for _, comp := range c.ConnectedComponents() {
		for _, w := range comp {
			if w == v {
				return comp
			}
		}
	}
	return nil
}
