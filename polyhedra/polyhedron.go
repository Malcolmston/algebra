package polyhedra

import (
	"errors"
	"fmt"
	"math"
	"sort"
)

// Edge is an undirected edge between two vertex indices. Its endpoints are
// stored in ascending order (A <= B) so that edges compare and map canonically.
type Edge struct {
	A, B int
}

// MakeEdge returns the canonical Edge connecting vertex indices i and j, with
// endpoints sorted so that A <= B.
func MakeEdge(i, j int) Edge {
	if i <= j {
		return Edge{i, j}
	}
	return Edge{j, i}
}

// Polyhedron is a boundary representation of a polyhedron: a list of vertex
// positions together with faces, each face being an ordered loop of vertex
// indices. For a closed, consistently oriented mesh every face loop is listed
// counter-clockwise when viewed from outside, so that face normals point
// outward.
type Polyhedron struct {
	Verts []Vec3
	Faces [][]int
}

// NewPolyhedron returns a Polyhedron with the given vertices and faces. The
// slices are retained (not copied); callers that intend to keep mutating their
// inputs should pass copies.
func NewPolyhedron(verts []Vec3, faces [][]int) *Polyhedron {
	return &Polyhedron{Verts: verts, Faces: faces}
}

// Clone returns a deep copy of the polyhedron.
func (p *Polyhedron) Clone() *Polyhedron {
	vs := make([]Vec3, len(p.Verts))
	copy(vs, p.Verts)
	fs := make([][]int, len(p.Faces))
	for i, f := range p.Faces {
		g := make([]int, len(f))
		copy(g, f)
		fs[i] = g
	}
	return &Polyhedron{Verts: vs, Faces: fs}
}

// NumVertices returns the number of vertices.
func (p *Polyhedron) NumVertices() int { return len(p.Verts) }

// NumFaces returns the number of faces.
func (p *Polyhedron) NumFaces() int { return len(p.Faces) }

// NumEdges returns the number of distinct undirected edges across all faces.
func (p *Polyhedron) NumEdges() int { return len(p.EdgeList()) }

// Validate reports the first structural problem with the polyhedron, or nil when
// it is well formed: every face must have at least three vertices and every
// index must be within range.
func (p *Polyhedron) Validate() error {
	n := len(p.Verts)
	for fi, f := range p.Faces {
		if len(f) < 3 {
			return fmt.Errorf("polyhedra: face %d has %d vertices, need >= 3", fi, len(f))
		}
		for _, idx := range f {
			if idx < 0 || idx >= n {
				return fmt.Errorf("polyhedra: face %d references vertex %d out of range [0,%d)", fi, idx, n)
			}
		}
	}
	return nil
}

// FaceLoop returns the vertex positions of face i as an ordered slice.
func (p *Polyhedron) FaceLoop(i int) []Vec3 {
	f := p.Faces[i]
	loop := make([]Vec3, len(f))
	for k, idx := range f {
		loop[k] = p.Verts[idx]
	}
	return loop
}

// EdgeList returns all distinct undirected edges of the polyhedron, sorted
// lexicographically by endpoint indices.
func (p *Polyhedron) EdgeList() []Edge {
	set := make(map[Edge]struct{})
	for _, f := range p.Faces {
		m := len(f)
		for i := 0; i < m; i++ {
			set[MakeEdge(f[i], f[(i+1)%m])] = struct{}{}
		}
	}
	out := make([]Edge, 0, len(set))
	for e := range set {
		out = append(out, e)
	}
	sort.Slice(out, func(i, j int) bool {
		if out[i].A != out[j].A {
			return out[i].A < out[j].A
		}
		return out[i].B < out[j].B
	})
	return out
}

// EdgeLengths returns the length of every edge in EdgeList order.
func (p *Polyhedron) EdgeLengths() []float64 {
	es := p.EdgeList()
	out := make([]float64, len(es))
	for i, e := range es {
		out[i] = p.Verts[e.A].Distance(p.Verts[e.B])
	}
	return out
}

// EulerCharacteristic returns V - E + F. For a simply connected closed
// polyhedron (topologically a sphere) it equals 2.
func (p *Polyhedron) EulerCharacteristic() int {
	return p.NumVertices() - p.NumEdges() + p.NumFaces()
}

// IsClosed reports whether the mesh is closed, i.e. every undirected edge is
// shared by exactly two faces.
func (p *Polyhedron) IsClosed() bool {
	count := make(map[Edge]int)
	for _, f := range p.Faces {
		m := len(f)
		for i := 0; i < m; i++ {
			count[MakeEdge(f[i], f[(i+1)%m])]++
		}
	}
	for _, c := range count {
		if c != 2 {
			return false
		}
	}
	return len(count) > 0
}

// FaceNormal returns the outward unit normal of face i (via Newell's method) and
// true, or the zero vector and false when the face is degenerate.
func (p *Polyhedron) FaceNormal(i int) (Vec3, bool) {
	return PolygonUnitNormal(p.FaceLoop(i))
}

// FaceNormals returns the unit normal of every face. Degenerate faces yield the
// zero vector.
func (p *Polyhedron) FaceNormals() []Vec3 {
	out := make([]Vec3, len(p.Faces))
	for i := range p.Faces {
		n, _ := p.FaceNormal(i)
		out[i] = n
	}
	return out
}

// FaceArea returns the area of face i.
func (p *Polyhedron) FaceArea(i int) float64 {
	return PolygonArea(p.FaceLoop(i))
}

// FaceAreas returns the area of every face.
func (p *Polyhedron) FaceAreas() []float64 {
	out := make([]float64, len(p.Faces))
	for i := range p.Faces {
		out[i] = p.FaceArea(i)
	}
	return out
}

// FaceCentroid returns the area-weighted centroid of face i and true, or the
// zero vector and false when face i is empty.
func (p *Polyhedron) FaceCentroid(i int) (Vec3, bool) {
	return PolygonCentroid(p.FaceLoop(i))
}

// FaceCentroids returns the centroid of every face (the zero vector for any
// empty face).
func (p *Polyhedron) FaceCentroids() []Vec3 {
	out := make([]Vec3, len(p.Faces))
	for i := range p.Faces {
		c, _ := p.FaceCentroid(i)
		out[i] = c
	}
	return out
}

// FaceSides returns the number of sides (vertices) of every face.
func (p *Polyhedron) FaceSides() []int {
	out := make([]int, len(p.Faces))
	for i, f := range p.Faces {
		out[i] = len(f)
	}
	return out
}

// SurfaceArea returns the total surface area of the polyhedron, the sum of its
// face areas.
func (p *Polyhedron) SurfaceArea() float64 {
	var s float64
	for i := range p.Faces {
		s += p.FaceArea(i)
	}
	return s
}

// SignedVolume returns the signed enclosed volume computed from the divergence
// theorem by fanning each face into triangles from its first vertex and summing
// signed tetrahedron volumes to the origin. For a closed mesh whose faces are
// oriented counter-clockwise as seen from outside, the result is positive and
// equals the enclosed volume.
func (p *Polyhedron) SignedVolume() float64 {
	var vol float64
	for _, f := range p.Faces {
		m := len(f)
		if m < 3 {
			continue
		}
		a := p.Verts[f[0]]
		for i := 1; i+1 < m; i++ {
			b := p.Verts[f[i]]
			c := p.Verts[f[i+1]]
			vol += a.Triple(b, c) / 6
		}
	}
	return vol
}

// Volume returns the (non-negative) enclosed volume of the polyhedron.
func (p *Polyhedron) Volume() float64 {
	return math.Abs(p.SignedVolume())
}

// VertexCentroid returns the arithmetic mean of the vertex positions and true,
// or the zero vector and false when the polyhedron has no vertices.
func (p *Polyhedron) VertexCentroid() (Vec3, bool) {
	return Centroid(p.Verts)
}

// Centroid returns the centroid of the solid region enclosed by the mesh
// (the center of mass of a uniform-density body) and true, or the vertex
// centroid and true as a fallback for a mesh of negligible volume, or false when
// the mesh has no faces. It uses the divergence-theorem decomposition into
// signed tetrahedra.
func (p *Polyhedron) Centroid() (Vec3, bool) {
	if len(p.Faces) == 0 {
		return Vec3{}, false
	}
	var (
		vol float64
		acc Vec3
	)
	for _, f := range p.Faces {
		m := len(f)
		if m < 3 {
			continue
		}
		a := p.Verts[f[0]]
		for i := 1; i+1 < m; i++ {
			b := p.Verts[f[i]]
			c := p.Verts[f[i+1]]
			w := a.Triple(b, c) / 6 // signed volume of tetra (origin,a,b,c)
			vol += w
			// centroid of that tetra is (origin+a+b+c)/4 = (a+b+c)/4
			acc = acc.Add(a.Add(b).Add(c).Scale(w / 4))
		}
	}
	if math.Abs(vol) < Eps*Eps {
		return p.VertexCentroid()
	}
	return acc.Div(vol), true
}

// BoundingBox returns the axis-aligned minimum and maximum corners of the vertex
// set and true, or two zero vectors and false when there are no vertices.
func (p *Polyhedron) BoundingBox() (min, max Vec3, ok bool) {
	if len(p.Verts) == 0 {
		return Vec3{}, Vec3{}, false
	}
	min, max = p.Verts[0], p.Verts[0]
	for _, v := range p.Verts[1:] {
		min = min.MinComp(v)
		max = max.MaxComp(v)
	}
	return min, max, true
}

// Translate returns a copy of the polyhedron with every vertex shifted by d.
func (p *Polyhedron) Translate(d Vec3) *Polyhedron {
	q := p.Clone()
	for i := range q.Verts {
		q.Verts[i] = q.Verts[i].Add(d)
	}
	return q
}

// ScaledBy returns a copy of the polyhedron with every vertex scaled about the
// origin by the factor s.
func (p *Polyhedron) ScaledBy(s float64) *Polyhedron {
	q := p.Clone()
	for i := range q.Verts {
		q.Verts[i] = q.Verts[i].Scale(s)
	}
	return q
}

// Centered returns a copy of the polyhedron translated so that its solid
// centroid lies at the origin.
func (p *Polyhedron) Centered() *Polyhedron {
	c, ok := p.Centroid()
	if !ok {
		return p.Clone()
	}
	return p.Translate(c.Neg())
}

// Circumradius returns the maximum distance from the centroid c to any vertex.
// It measures how far the outermost vertices reach.
func (p *Polyhedron) Circumradius() float64 {
	c, ok := p.Centroid()
	if !ok {
		return 0
	}
	var r float64
	for _, v := range p.Verts {
		if d := v.Distance(c); d > r {
			r = d
		}
	}
	return r
}

// MinVertexRadius returns the minimum distance from the centroid to any vertex.
func (p *Polyhedron) MinVertexRadius() float64 {
	c, ok := p.Centroid()
	if !ok {
		return 0
	}
	r := math.Inf(1)
	for _, v := range p.Verts {
		if d := v.Distance(c); d < r {
			r = d
		}
	}
	if math.IsInf(r, 1) {
		return 0
	}
	return r
}

// Inradius returns the minimum distance from the centroid to any face plane,
// i.e. the radius of the largest sphere centered at the centroid that fits
// inside the mesh (exact for a convex polyhedron).
func (p *Polyhedron) Inradius() float64 {
	c, ok := p.Centroid()
	if !ok {
		return 0
	}
	r := math.Inf(1)
	for i := range p.Faces {
		n, ok := p.FaceNormal(i)
		if !ok {
			continue
		}
		fc, _ := p.FaceCentroid(i)
		d := math.Abs(fc.Sub(c).Dot(n))
		if d < r {
			r = d
		}
	}
	if math.IsInf(r, 1) {
		return 0
	}
	return r
}

// Midradius returns the mean distance from the centroid to the edge midpoints
// (the midradius of a canonical solid, where all such distances coincide).
func (p *Polyhedron) Midradius() float64 {
	c, ok := p.Centroid()
	if !ok {
		return 0
	}
	es := p.EdgeList()
	if len(es) == 0 {
		return 0
	}
	var s float64
	for _, e := range es {
		mid := p.Verts[e.A].Midpoint(p.Verts[e.B])
		s += mid.Distance(c)
	}
	return s / float64(len(es))
}

// Triangulate returns a new polyhedron whose faces are all triangles, obtained
// by fanning every face from its first vertex. Vertices are shared with the
// original (the slice is copied).
func (p *Polyhedron) Triangulate() *Polyhedron {
	vs := make([]Vec3, len(p.Verts))
	copy(vs, p.Verts)
	var tris [][]int
	for _, f := range p.Faces {
		for i := 1; i+1 < len(f); i++ {
			tris = append(tris, []int{f[0], f[i], f[i+1]})
		}
	}
	return &Polyhedron{Verts: vs, Faces: tris}
}

// VertexDegrees returns, for each vertex index, the number of edges incident to
// it.
func (p *Polyhedron) VertexDegrees() []int {
	deg := make([]int, len(p.Verts))
	for _, e := range p.EdgeList() {
		deg[e.A]++
		deg[e.B]++
	}
	return deg
}

// FacesAroundVertex returns the indices of the faces incident to vertex v, in
// ascending face-index order.
func (p *Polyhedron) FacesAroundVertex(v int) []int {
	var out []int
	for fi, f := range p.Faces {
		for _, idx := range f {
			if idx == v {
				out = append(out, fi)
				break
			}
		}
	}
	return out
}

// OrientOutward returns a copy of the mesh in which every face is oriented so
// that its normal points away from the mesh centroid. This is correct for convex
// polyhedra and any mesh that is star-shaped from its centroid.
func (p *Polyhedron) OrientOutward() *Polyhedron {
	c, ok := p.Centroid()
	if !ok {
		c, _ = p.VertexCentroid()
	}
	q := p.Clone()
	for fi := range q.Faces {
		n, ok := q.FaceNormal(fi)
		if !ok {
			continue
		}
		fc, _ := q.FaceCentroid(fi)
		if fc.Sub(c).Dot(n) < 0 {
			reverse(q.Faces[fi])
		}
	}
	return q
}

// reverse reverses an int slice in place.
func reverse(s []int) {
	for i, j := 0, len(s)-1; i < j; i, j = i+1, j-1 {
		s[i], s[j] = s[j], s[i]
	}
}

// ErrEmptyMesh is returned when an operation requires a non-empty mesh.
var ErrEmptyMesh = errors.New("polyhedra: empty mesh")

// Dual returns the polar/canonical dual polyhedron. Each face of p becomes a
// vertex placed at that face's centroid, and each vertex of p becomes a face
// whose vertices are the duals of the faces surrounding it, ordered
// consistently around the original vertex. The mesh must be closed and convex
// (as produced by the solid constructors and ConvexHull) for the result to be a
// valid polyhedron; otherwise ErrEmptyMesh or ErrDegenerate is returned.
func (p *Polyhedron) Dual() (*Polyhedron, error) {
	if len(p.Faces) == 0 || len(p.Verts) == 0 {
		return nil, ErrEmptyMesh
	}
	center, _ := p.Centroid()

	// New vertex per original face: its centroid.
	dualVerts := make([]Vec3, len(p.Faces))
	faceNormals := make([]Vec3, len(p.Faces))
	for fi := range p.Faces {
		fc, _ := p.FaceCentroid(fi)
		dualVerts[fi] = fc
		n, _ := p.FaceNormal(fi)
		faceNormals[fi] = n
	}

	var dualFaces [][]int
	for v := range p.Verts {
		around := p.FacesAroundVertex(v)
		if len(around) < 3 {
			continue
		}
		// Order the surrounding faces angularly around the direction from the
		// solid center to the vertex.
		axis := p.Verts[v].Sub(center)
		u, ok := axis.Normalize()
		if !ok {
			u = ZAxis()
		}
		// Build an orthonormal basis (e1, e2) spanning the plane orthogonal to u.
		ref := XAxis()
		if math.Abs(u.Dot(ref)) > 0.9 {
			ref = YAxis()
		}
		e1, _ := ref.Sub(u.Scale(ref.Dot(u))).Normalize()
		e2 := u.Cross(e1)
		type fa struct {
			idx int
			ang float64
		}
		fas := make([]fa, len(around))
		for i, fi := range around {
			d := dualVerts[fi].Sub(p.Verts[v])
			fas[i] = fa{fi, math.Atan2(d.Dot(e2), d.Dot(e1))}
		}
		sort.Slice(fas, func(i, j int) bool { return fas[i].ang < fas[j].ang })
		loop := make([]int, len(fas))
		for i, f := range fas {
			loop[i] = f.idx
		}
		dualFaces = append(dualFaces, loop)
	}
	dual := &Polyhedron{Verts: dualVerts, Faces: dualFaces}
	return dual.OrientOutward(), nil
}
