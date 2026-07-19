package meshgen

// Quad is an axis-aligned rectangular cell defined by its minimum and maximum
// corners.
type Quad struct {
	Min, Max Vec2
}

// NewQuad returns the quad with the given corners.
func NewQuad(min, max Vec2) Quad { return Quad{Min: min, Max: max} }

// Center returns the centre point of the quad.
func (q Quad) Center() Vec2 { return q.Min.Add(q.Max).Div(2) }

// Width returns the extent of the quad along X.
func (q Quad) Width() float64 { return q.Max.X - q.Min.X }

// Height returns the extent of the quad along Y.
func (q Quad) Height() float64 { return q.Max.Y - q.Min.Y }

// Area returns the area of the quad.
func (q Quad) Area() float64 { return q.Width() * q.Height() }

// Contains reports whether p lies within the quad, including its boundary.
func (q Quad) Contains(p Vec2) bool {
	return p.X >= q.Min.X && p.X <= q.Max.X && p.Y >= q.Min.Y && p.Y <= q.Max.Y
}

// Corners returns the four corners of the quad in counterclockwise order
// starting from the minimum corner.
func (q Quad) Corners() [4]Vec2 {
	return [4]Vec2{
		{q.Min.X, q.Min.Y},
		{q.Max.X, q.Min.Y},
		{q.Max.X, q.Max.Y},
		{q.Min.X, q.Max.Y},
	}
}

// Subdivide splits the quad into its four child quadrants in the order
// SW, SE, NE, NW.
func (q Quad) Subdivide() [4]Quad {
	c := q.Center()
	return [4]Quad{
		{Vec2{q.Min.X, q.Min.Y}, Vec2{c.X, c.Y}},
		{Vec2{c.X, q.Min.Y}, Vec2{q.Max.X, c.Y}},
		{Vec2{c.X, c.Y}, Vec2{q.Max.X, q.Max.Y}},
		{Vec2{q.Min.X, c.Y}, Vec2{c.X, q.Max.Y}},
	}
}

// QuadNode is a node of a Quadtree. A leaf has nil Children.
type QuadNode struct {
	Bounds   Quad
	Depth    int
	Children [4]*QuadNode
}

// IsLeaf reports whether the node has no children.
func (n *QuadNode) IsLeaf() bool { return n.Children[0] == nil }

// Quadtree is an adaptively refined spatial subdivision of a rectangular
// region.
type Quadtree struct {
	Root *QuadNode
}

// BuildQuadtree builds a quadtree over bounds, subdividing a cell whenever
// refine reports true for it, down to at most maxDepth levels. The refine
// predicate receives the candidate cell and its depth.
func BuildQuadtree(bounds Quad, maxDepth int, refine func(cell Quad, depth int) bool) *Quadtree {
	var build func(q Quad, depth int) *QuadNode
	build = func(q Quad, depth int) *QuadNode {
		node := &QuadNode{Bounds: q, Depth: depth}
		if depth < maxDepth && refine(q, depth) {
			kids := q.Subdivide()
			for i := 0; i < 4; i++ {
				node.Children[i] = build(kids[i], depth+1)
			}
		}
		return node
	}
	return &Quadtree{Root: build(bounds, 0)}
}

// Leaves returns the leaf cells of the quadtree in depth-first order.
func (t *Quadtree) Leaves() []Quad {
	var out []Quad
	var walk func(n *QuadNode)
	walk = func(n *QuadNode) {
		if n == nil {
			return
		}
		if n.IsLeaf() {
			out = append(out, n.Bounds)
			return
		}
		for _, c := range n.Children {
			walk(c)
		}
	}
	walk(t.Root)
	return out
}

// LeafCenters returns the centre point of every leaf cell.
func (t *Quadtree) LeafCenters() []Vec2 {
	leaves := t.Leaves()
	out := make([]Vec2, len(leaves))
	for i, q := range leaves {
		out[i] = q.Center()
	}
	return out
}

// LeafCount returns the number of leaf cells.
func (t *Quadtree) LeafCount() int { return len(t.Leaves()) }

// Depth returns the maximum depth reached by any node.
func (t *Quadtree) Depth() int {
	var maxD int
	var walk func(n *QuadNode)
	walk = func(n *QuadNode) {
		if n == nil {
			return
		}
		if n.Depth > maxD {
			maxD = n.Depth
		}
		for _, c := range n.Children {
			walk(c)
		}
	}
	walk(t.Root)
	return maxD
}

// LeafCorners returns the deduplicated set of all leaf-cell corner points of the
// quadtree, suitable as vertices for a background mesh.
func (t *Quadtree) LeafCorners() []Vec2 {
	seen := make(map[[2]float64]struct{})
	var out []Vec2
	for _, q := range t.Leaves() {
		for _, c := range q.Corners() {
			key := [2]float64{c.X, c.Y}
			if _, ok := seen[key]; !ok {
				seen[key] = struct{}{}
				out = append(out, c)
			}
		}
	}
	return out
}
